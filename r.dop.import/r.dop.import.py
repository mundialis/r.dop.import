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
from datetime import datetime

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
)
from grass_gis_helpers.open_geodata_germany.federal_state import (
    get_federal_states,
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


def collect_metadata(fs, dop_list, license_info=None, base_url=None):
    """Collect metadata for imported DOPs

    Args:
        fs (str): Federal state abbreviation
        dop_list (list): List of imported DOP raster names
        license_info (str): License information (optional)
        base_url (str): Base URL of the data source (optional)

    Returns:
        dict: Metadata dictionary for this federal state

    """
    return {
        "federal_state": fs,
        "download_date": datetime.now().strftime("%d.%m.%Y"),
        "license": license_info,
        "base_url": base_url,
        "dop_rasters": dop_list,
        "count": len(dop_list),
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

        import subprocess

        result = subprocess.run(
            ["g.manual", "-m", addon_name], capture_output=True, text=True
        )

        license_info = None
        base_url = None

        if result.returncode == 0:
            html_path = result.stdout.strip()
            if os.path.exists(html_path):
                with open(html_path, "r", encoding="utf-8") as f:
                    content = f.read()

                    import re

                    license_patterns = [
                        r"<h[23]>.*?[Ll]icen[sc]e.*?</h[23]>(.*?)<h[23]>",
                        r"<h[23]>.*?[Ll]izenz.*?</h[23]>(.*?)<h[23]>",
                        r"[Ll]icen[sc]e:</?\w*>(.*?)</\w+>",
                        r"[Ll]izenz:</?\w*>(.*?)</\w+>",
                    ]

                    for pattern in license_patterns:
                        license_match = re.search(
                            pattern, content, re.DOTALL | re.IGNORECASE
                        )
                        if license_match:
                            from html.parser import HTMLParser

                            class MLStripper(HTMLParser):

                                def __init__(self):
                                    super().__init__()
                                    self.strict = False
                                    self.convert_charrefs = True
                                    self.text = []

                                def handle_data(self, d):
                                    self.text.append(d)

                                def get_data(self):
                                    return "".join(self.text)

                            s = MLStripper()
                            s.feed(license_match.group(1))
                            license_info = s.get_data().strip
                            if license_info:
                                break

                    url_patterns = [
                        r'[Bb]ase\s*URL[:\s]+<?([https]+://[^\s<>"]+)>?',
                        r'[Dd]ata\s*[Ss]ource[:\s]+<?([https]+://[^\s<>"]+)>?',
                        r'[Dd]ownload\s*URL[:\s]+<?([https]+://[^\s<>"]+)>?',
                        r'href=["\'](https?://[^"\']+geoportal[^"\']+)["\']',
                        r'href=["\'](https?://[^"\']+open.*data[^"\']+)["\']',
                    ]

                    for pattern in url_patterns:
                        url_match = re.search(pattern, content, re.IGNORECASE)
                        if url_match:
                            base_url = url_match.group(1)
                            break

        return license_info, base_url

    except Exception as e:
        grass.warning(f"Could not extract license/URL for {fs}: {e}")

    return None, None


def get_federal_state_name(fs):
    """Get full name of federal state from abbreviation

    Args:
        fs (str): Federal state abbreviation

    Returns:
        str: Full federal state name

    """

    federal_state_names = {
        "BW": "Baden-Württemberg",
        "BY": "Bayern",
        "BE": "Berlin",
        "BB": "Brandenburg",
        "HB": "Bremen",
        "HH": "Hamburg",
        "HE": "Hessen",
        "MV": "Mecklenburg-Vorpommern",
        "NI": "Niedersachsen",
        "NW": "Nordrhein-Westfalen",
        "RP": "Rheinland-Pfalz",
        "SL": "Saarland",
        "SN": "Sachsen",
        "ST": "Sachsen-Anhalt",
        "SH": "Schleswig-Holstein",
        "TH": "Thüringen",
    }
    return federal_state_names.get(fs, fs)


def write_metadata_markdown(metadata_list, output_name, download_dir):
    """Write metadata to Markdown file

    Args:
        metadata_list (list): List of metadata dictionaries
        output_name (str): Base output name
        download_dir (str): Directory where metadata file should be saved

    """

    # Get name of mapset
    mapset = grass.gisenv()["MAPSET"]

    # Create markdown file
    md_file = os.path.join(download_dir, f"{output_name}_metadata.md")

    with open(md_file, "w", encoding="utf-8") as f:
        # Header
        f.write(f"# Metadaten der DOPs im Mapset {mapset}\n\n")

        # Date of Download
        if metadata_list:
            f.write(
                f"**Downloaddatum:** {metadata_list[0]['download_date']}\n\n"
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
                "bezogen:\n\n"
            )

            for dop in fs_meta["dop_rasters"]:
                f.write(f"- `{dop}`\n")

            f.write(f"\n**Anzahl:** {fs_meta['count']}\n\n")

        # Additional info
        f.write("---\n\n")
        f.write(f"*Erstellt am {datetime.now().strftime('%d.%m.%Y')}\n")

    grass.message(_(f"Metadata file created: {md_file}"))


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
                grass.run_command(addon, **params)
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
        fs_metadata = collect_metadata(fs, fs_dop_list, license_info, base_url)
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
    write_metadata_markdown(metadata_list, output, download_dir)

    grass.message(_(f"DOP group <{output}> is created."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
