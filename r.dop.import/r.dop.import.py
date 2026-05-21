#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Anika Weinmann, Julia Haas
# PURPOSE:     Downloads Digital Orthophotos (DOPs) wihtin a specified
#              federal state and area of interest
# SPDX-FileCopyrightText: (c) 2024-2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
############################################################################

# %Module
# % description: Downloads and imports DOPs for specified federal state and AOI
# % keyword: raster
# % keyword: import
# % keyword: imagery
# % keyword: download
# % keyword: DOP
# % keyword: open-geodata-germany
# %end

# %option G_OPT_V_INPUT
# % key: aoi
# % description: Vector map to restrict DOPs import to
# % required: no
# %end

# %option
# % key: federal_state
# % multiple: yes
# % required: no
# % description: Federal state(s) related to the area of interest, e.g.: "Nordrhein-Westfalen"
# % options: Baden-Württemberg,BW,Bayern,BY,Berlin,BE,Brandenburg,BB,Hamburg,HH,Nordrhein-Westfalen,NW,Sachsen,SN,Thüringen,TH,Rheinland-Pfalz,RP,Niedersachsen,NI
# %end

# %option G_OPT_F_INPUT
# % key: federal_state_file
# % description: Path to text file containing the federal state(s) related to the area of interest
# % required: no
# %end

# %option G_OPT_M_DIR
# % key: local_data_dir
# % required: no
# % description: Directory with raster map of DOPs to import (e.g. VRT)
# %end

# %option
# % key: download_dir
# % label: Path of output folder
# % description: Path of download folder
# % required: no
# % multiple: no
# %end

# %option
# % key: hist_year
# % description: Download historic data for given year (Note: Currently only supported for NW)
# % required: no
# %end

# %option G_OPT_R_OUTPUT
# % description: Name for output raster map
# %end

# %option
# % key: nprocs
# % type: integer
# % required: no
# % multiple: no
# % label: Number of parallel processes
# % description: Number of cores for multiprocessing, -2 is the number of available cores - 1
# % answer: -2
# %end

# %option G_OPT_MEMORYMB
# %end

# %option
# % key: metadata
# % type: string
# % required: no
# % multiple: no
# % description: Path to metadata output file (markdown format)
# % answer:
# %end

# %flag
# % key: k
# % description: Keep downloaded data in the download directory
# %end

# %flag
# % key: r
# % description: Use native data resolution
# %end

# %rules
# % required: federal_state, federal_state_file
# % excludes: federal_state_file, federal_state
# % requires_all: -k, download_dir
# %end

import atexit
import os
import sys
import re
import pathlib
import traceback
from datetime import datetime
from html.parser import HTMLParser

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
)
from grass_gis_helpers.open_geodata_germany.federal_state import (
    get_federal_states,
    FS_ABBREVIATION,
)
from grass_gis_helpers.raster import create_vrt
from grass_gis_helpers.data_import import import_local_raster_data

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the DOP library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import OPEN_DATA_AVAILABILITY
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

# set global varibales
ID = grass.tempname(12)
ORIG_REGION = f"original_region_{ID}"
rm_rasters = []
SUPPORTED = OPEN_DATA_AVAILABILITY["SUPPORTED"]
NO_OPEN_DATA = OPEN_DATA_AVAILABILITY["NO_OPEN_DATA"]
NOT_YET_SUPPORTED = OPEN_DATA_AVAILABILITY["NOT_YET_SUPPORTED"]


def cleanup():
    """Remove all not needed files at the end"""
    general_cleanup(
        orig_region=ORIG_REGION,
        rm_rasters=rm_rasters,
    )


def import_local_data(aoi, out, local_data_dir, fs, all_dops, native_res_flag):
    """Import local DOP data

    Args:
        aoi (str): Vector map with area of interest
        out (str): Base output name
        local_data_dir (str): Path to local data directory with federal state
                              subfolders
        fs (str): the abbrivation of the federal state
        all_dops (list): empty list where the imported DOP rasters
                         will be appended
        native_res_flag (bool): True if native data resolution should be used

    """
    imported_local_data = import_local_raster_data(
        aoi,
        f"{out}_{fs}",
        os.path.join(local_data_dir, fs),
        native_res_flag,
        all_dops,
        rm_rasters,
        band_dict={1: "red", 2: "green", 3: "blue", 4: "nir"},
    )

    if not imported_local_data and fs in ["BW"]:
        grass.fatal(_("Local data does not overlap with AOI."))
    elif not imported_local_data:
        grass.message(
            _(
                "Local data does not overlap with AOI. Data will be downloaded"
                " from Open Data portal.",
            ),
        )
    return imported_local_data


def get_dop_urls_from_tindex(fs):
    """Extract DOP download URLs from tile index vector

    Args:
        fs (str): Federal state abbreviation

    Returns:
        list: List of download URLs

    """

    dop_urls = []

    # Search for TINDEX in all mapsets
    tindex_found = None
    mapsets = grass.read_command("g.mapsets", flags="p").strip().split()

    for mapset in mapsets:
        vectors = grass.list_grouped("vector").get(mapset, [])
        matching = [v for v in vectors if "tindex" in v.lower()]
        if matching:
            tindex_found = f"{matching[0]}@{mapset}"
            grass.debug(f"Found TINDEX: {tindex_found}")
            break

    if not tindex_found:
        grass.debug(f"No TINDEX found for {fs}")
        return dop_urls

    # Read attributes of TINDEX
    try:
        columns_info = grass.vector_colums(tindex_found)

        if "location" not in columns_info:
            grass.debug(
                f"'location' column not found in TINDEX {tindex_found}",
            )
            grass.debug(f"Available columns: {list(columns_info.keys())}")
            return dop_urls

        grass.debug("Using 'location' column from TINDEX")

        # Read URLs from table
        result = grass.read_command(
            "v.db.select",
            map=tindex_found,
            columns="location",
            flags="c",
        ).strip()

        if result:
            for line in result.split("\n"):
                line = line.strip()
                if not line:
                    continue

                # Parse URLs from location attribute
                urls_in_line = []

                if "," in line:
                    parts = [p.strip() for p in line.split(",")]
                else:
                    parts = [line]

                for part in parts:
                    https_urls = re.findall(r"https://[^\s,]+", part)
                    urls_in_line.extend(https_urls)

            dop_urls.extend(urls_in_line)

        seen = set()
        dop_urls = [
            url for url in dop_urls if not (url in seen or seen.add(url))
        ]

        grass.debug(f"Extracted {len(dop_urls)} URLs from TINDEX")

    except Exception as e:
        grass.debug(f"Error reading TINDEX {tindex_found}: {e}")

        grass.debug(traceback.format_exc())

    return dop_urls


def extract_dop_info_from_url(url):
    """Extract DOP filename from download URL

    Args:
        url (str): Download URL

    Returns:
        str: DOP filename

    """

    url_clean = url.split("?")[0]

    filename = os.path.basename(url_clean)

    if ".zip" in url_clean:
        match = re.search(
            r"([^/]+\.(?:tif|tiff|jp2))(?:/|\.zip|$)",
            url_clean,
            re.IGNORECASE,
        )
        if match:
            filename = match.group(1)

    return filename.replace("/vsizip/", "").replace("vsicurl/", "")


def collect_metadata(
    fs,
    dop_list,
    license_info=None,
    base_url=None,
    dop_names=None,
    dop_urls=None,
):
    """Collect metadata for imported DOPs

    Args:
        fs (str): Federal state abbreviation
        dop_list (list): List of imported DOP raster names
        license_info (str): License information (optional)
        base_url (str): Base URL of the data source (optional)
        dop_names (list): List of original downloaded DOP names/files (optional)
        dop_urls (list): List of download URLs (optional)

    Returns:
        dict: Metadata dictionary for this federal state

    """

    if dop_names and len(dop_names) > 0:
        original_files = dop_names
    elif dop_urls and len(dop_urls) > 0:
        original_files = [extract_dop_info_from_url(url) for url in dop_urls]
    else:
        unique_dops = set()
        for dop in dop_list:
            base_name = re.sub(r"_(red|green|blue|nir)$", "", dop)
            unique_dops.add(base_name)
        original_files = sorted(list(unique_dops))

    return {
        "federal_state": fs,
        "download_date": datetime.now().strftime("%d.%m.%Y"),
        "license": license_info,
        "base_url": base_url,
        "dop_rasters": sorted(list(set(original_files))),
        "dop_urls": dop_urls or [],
        "count": len(original_files),
    }


def get_license_and_url_from_addon(fs):
    """Extract license information and base URL from federal state specific addon

    Args:
        fs (str): Federal state abbreviation

    Returns:
        tuple: (license_info, base_url)

    """

    try:
        if fs in ["BB", "BE"]:
            addon_name = "r.dop.import.bb.be"
        else:
            addon_name = f"r.dop.import.{fs.lower()}"

        # Path to HTML file
        html_file = os.path.join(
            pathlib.Path("~").expanduser(),
            ".grass8",
            "addons",
            "docs",
            "html",
            f"{addon_name}.html",
        )

        if not pathlib.Path(html_file).exists():
            grass.debug(f"HTML file not found: {html_file}")
            return None, None

        # Read HTML file
        html_content = pathlib.Path(html_file).read_text(encoding="utf-8")

        if not html_content:
            grass.debug(f"HTML file is empty: {html_file}")
            return None, None

        license_info = None
        base_url = None

        # Search for license information
        id_match = re.search(
            r"<br>\s*id:\s*([^,\n]+)",
            html_content,
            re.IGNORECASE,
        )
        name_match = re.search(
            r"<br>\s*name:\s*([^,\n]+)",
            html_content,
            re.IGNORECASE,
        )
        url_match = re.search(
            r"<br>\s*url:\s*(https?://[^,\s\n]+)",
            html_content,
            re.IGNORECASE,
        )
        source_match = re.search(
            r"<br>\s*source:\s*(.+?)(?=\n|<h\d>)",
            html_content,
            re.DOTALL | re.IGNORECASE,
        )

        if id_match and name_match and url_match and source_match:
            license_id = id_match.group(1).strip()
            license_name = name_match.group(1).strip()
            license_url = url_match.group(1).strip()
            source_html = source_match.group(1).strip()

            # Remove html tags from source info
            class MLStripper(HTMLParser):

                def __init__(self) -> None:
                    super().__init__()
                    self.strict = False
                    self.convert_charrefs = True
                    self.text = []

                def handle_data(self, data) -> None:
                    self.text.append(data)

                def get_data(self) -> None:
                    return "".join(self.text)

            s = MLStripper()
            s.feed(source_html)
            source_clean = s.get_data().strip()

            # Create formatted license information
            license_info = (
                f"{license_name} ({license_id}), "
                f"{license_url}, "
                f"Quelle: {source_clean}"
            )

            # Get base URL from source
            source_link_match = re.search(
                r'href=["\']([^"\']+)["\']',
                source_html,
            )
            if source_link_match:
                base_url = source_link_match.group(1).strip()

        return license_info, base_url

    except Exception as e:
        grass.warning(f"Could not extract license/URL for {fs}: {e}")

    return None, None


def get_federal_state_name(fs_abbr):
    """Get full name of federal state from abbreviation

    Args:
        fs_abbr (str): Federal state abbreviation

    Returns:
        str: Full federal state name

    """

    for name, abbr in FS_ABBREVIATION.items():
        if abbr == fs_abbr and name != abbr:
            return name
    return fs_abbr


def write_metadata_markdown(metadata_list, metadata_path=None):
    """Write metadata to Markdown file

    Args:
        metadata_list (list): List of metadata dictionaries
        metadata_path (str): Full path to metadata file (optional)

    """

    # Don't write metadata file if path not set
    if not metadata_path or metadata_path == "":
        grass.message(
            _("No metadata path specified. Skipping metadata file creation."),
        )
        return

    # Make sure path ends with ".md"
    if not metadata_path.endswith(".md"):
        metadata_path = f"{metadata_path}.md"

    # Create directory if it does not exist yet
    metadata_dir = os.path.dirname(metadata_path)
    if metadata_dir and not pathlib.Path(metadata_dir).exists():
        try:
            pathlib.Path(metadata_dir).mkdir(parents=True)
            grass.message(_(f"Created directory: {metadata_dir}"))
        except Exception as e:
            grass.warning(_(f"Could not create directory {metadata_dir}: {e}"))
            return

    mapset = grass.gisenv()["MAPSET"]

    # Create markdown file
    try:
        with pathlib.Path(metadata_path).open("w", encoding="utf-8") as f:
            # Header
            f.write(f"# Metadaten der DOPs im Mapset {mapset}\n\n")

            if metadata_list:
                f.write(
                    f"**Downloaddatum:** {metadata_list[0]['download_date']}\n\n",
                )

            # Licenses
            f.write("## Lizenzen\n\n")
            unique_licenses = {}
            for fs_meta in metadata_list:
                fs_name = get_federal_state_name(fs_meta["federal_state"])
                if fs_meta["license"] not in unique_licenses.values():
                    unique_licenses[fs_name] = fs_meta["license"]

            for fs_name, license_text in unique_licenses.items():
                f.write(f"**{fs_name}:** {license_text}\n\n")

            f.write("## Heruntergeladene DOPs\n\n")

            for fs_meta in metadata_list:
                fs_name = get_federal_state_name(fs_meta["federal_state"])
                base_url = fs_meta["base_url"]

                f.write(
                    f"### Folgende DOPs wurden aus {fs_name} ({base_url}) "
                    "bezogen:\n\n",
                )

                # Show file name with link if URLs are found
                if fs_meta.get("dop_urls") and len(fs_meta["dop_urls"]) > 0:
                    for url in fs_meta["dop_urls"]:
                        filename = extract_dop_info_from_url(url)
                        f.write(f"- [{filename}]({url})\n")
                    f.write(f"\n**Anzahl:** {fs_meta['count']}\n\n")

                else:
                    for dop in fs_meta["dop_rasters"]:
                        if "DOP-Kacheln" in dop or "via WMS" in dop:
                            f.write(f"{dop}\n\n")
                        else:
                            f.write(f"- `{dop}`\n")

                    if not any(
                        "DOP-Kacheln" in dop for dop in fs_meta["dop_rasters"]
                    ):
                        f.write(f"\n**Anzahl:** {fs_meta['count']}\n\n")
                    else:
                        f.write("\n")

                f.writelines(f"- `{dop}`\n" for dop in fs_meta["dop_rasters"])

                f.write(f"\n**Anzahl:** {fs_meta['count']}\n\n")

            # Additional info
            f.write("---\n\n")
            f.write(f"*Erstellt am {datetime.now().strftime('%d.%m.%Y')}\n")

        grass.message(_(f"Metadata file created: {metadata_path}"))

    except Exception as e:
        grass.warning(_(f"Could not write metadata file: {e}"))


def main():
    """Main function of r.dop.import"""
    aoi = options["aoi"]
    federal_states = get_federal_states(
        options["federal_state"],
        options["federal_state_file"],
    )
    local_data_dir = options["local_data_dir"]
    download_dir = check_download_dir(options["download_dir"])
    hist_year = options["hist_year"]
    output = options["output"]
    nprocs = options["nprocs"]
    memory = options["memory"]
    metadata_path = options.get("metadata", "")
    keep_data = flags["k"]
    native_res = flags["r"]

    # save original region
    grass.run_command("g.region", save=ORIG_REGION, quiet=True)

    # local DOP files
    local_fs_list = []
    if local_data_dir and local_data_dir != "":
        local_fs_list = os.listdir(local_data_dir)

    # loop over federal states and import data
    all_dops = {"red": [], "green": [], "blue": [], "nir": []}
    metadata_list = []
    for fs in set(federal_states):
        grass.message(_(f"Importing DOPs for {fs}..."))
        fs_dop_list = []
        dop_names = []
        dop_urls = []

        # check if local data for federal state given
        imported_local_data = False
        if fs in local_fs_list:
            imported_local_data = import_local_data(
                aoi,
                output,
                local_data_dir,
                fs,
                all_dops,
                native_res,
            )
            if imported_local_data:
                fs_dop_list = [
                    f"{output}_{fs}_red",
                    f"{output}_{fs}_green",
                    f"{output}_{fs}_blue",
                    f"{output}_{fs}_nir",
                ]
                local_fs_dir = os.path.join(local_data_dir, fs)
                if pathlib.Path(local_fs_dir).exists():
                    for files in os.walk(local_fs_dir):
                        dop_names.extend(
                            file
                            for file in files
                            if file.lower().endswith(
                                ".tif", ".tiff", ".jp2", ".jpeg",
                            )
                        )

        elif fs in NO_OPEN_DATA:
            grass.fatal(
                _(
                    f"No local data for {fs} available. For the federal state "
                    "there are no open data available. Is the path correct?",
                ),
            )

        # import data when local import was not used
        if not imported_local_data:
            # implement data download and import from open data
            out_fs = f"dop_{fs}_{ID}"

            # not yet supported
            if fs in NOT_YET_SUPPORTED:
                grass.fatal(
                    _(
                        "The import of the open data is not yet supported for "
                        f"{fs}.",
                    ),
                )
            # no open data available
            elif fs in NO_OPEN_DATA:
                grass.fatal(
                    _(
                        f"For the federal state {fs} there are no open data "
                        "available. Please use local data <local_data_dir>.",
                    ),
                )
            else:
                r_dop_import_fs_flags = ""
                if keep_data:
                    r_dop_import_fs_flags += "k"
                if native_res:
                    r_dop_import_fs_flags += "r"
                # change fs for BB and BE (within one addon)
                if fs in ["BB", "BE"]:
                    addon = "r.dop.import.bb.be"
                    worker_addon = "r.dop.import.worker.bb.be"
                else:
                    addon = f"r.dop.import.{fs.lower()}"
                    worker_addon = f"r.dop.import.worker.{fs.lower()}"
                params = {
                    "aoi": aoi,
                    "download_dir": download_dir,
                    "output": out_fs,
                    "memory": memory,
                    "flags": r_dop_import_fs_flags,
                    "overwrite": True,
                }
                if fs == "NW":
                    params["hist_year"] = hist_year
                if grass.find_program(worker_addon, "--help"):
                    params["nprocs"] = nprocs

                # Run addon
                grass.run_command(addon, **params)

                dop_urls = get_dop_urls_from_tindex(fs)

                if not dop_urls and (
                    keep_data
                    and download_dir
                    and pathlib.Path(download_dir).exists()
                ):
                    for files in os.walk(download_dir):
                        dop_names.extend(
                            file
                            for file in files
                            if file.lower().endswith(
                                ".tif", ".tiff", ".jp2", ".jpeg",
                            )
                        )

                if not dop_urls and not dop_names:
                    imported_rasters = grass.list_grouped("raster").get(
                        grass.gisenv()["MAPSET"],
                        [],
                    )
                    matching_rasters = [
                        r for r in imported_rasters if r.startswith(out_fs)
                    ]
                    num_dop_tiles = (
                        len(matching_rasters) // 4
                        if len(matching_rasters) >= 4
                        else 0
                    )

                    if num_dop_tiles > 0:
                        if fs in ["BW"]:
                            dop_names = [
                                (
                                    f"{num_dop_tiles} DOP-Kacheln "
                                    "(via WMS heruntergeladen)"
                                ),
                            ]
                        else:
                            dop_names = [f"{num_dop_tiles} DOP-Kacheln"]

            all_dops["red"].append(f"{out_fs}_red")
            all_dops["green"].append(f"{out_fs}_green")
            all_dops["blue"].append(f"{out_fs}_blue")
            all_dops["nir"].append(f"{out_fs}_nir")

            fs_dop_list = [
                f"{out_fs}_red",
                f"{out_fs}_green",
                f"{out_fs}_blue",
                f"{out_fs}_nir",
            ]

        # Collect metadata for this federal state
        license_info, base_url = get_license_and_url_from_addon(fs)
        fs_metadata = collect_metadata(
            fs,
            fs_dop_list,
            license_info,
            base_url,
            dop_names,
            dop_urls,
        )
        metadata_list.append(fs_metadata)

    create_vrt(all_dops["red"], f"{output}_red")
    create_vrt(all_dops["green"], f"{output}_green")
    create_vrt(all_dops["blue"], f"{output}_blue")
    create_vrt(all_dops["nir"], f"{output}_nir")
    grass.run_command(
        "i.group",
        group=output,
        input=[
            f"{output}_red",
            f"{output}_green",
            f"{output}_blue",
            f"{output}_nir",
        ],
        overwrite=True,
    )

    # Write metadata file
    write_metadata_markdown(metadata_list, metadata_path)

    grass.message(_(f"DOP group <{output}> is created."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
