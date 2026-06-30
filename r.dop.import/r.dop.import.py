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
import pathlib

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
)
from grass_gis_helpers.open_geodata_germany.federal_state import (
    get_federal_states,
)
from grass_gis_helpers.open_geodata_germany.metadata import (
    collect_metadata,
    get_license_and_url_from_addon,
    write_metadata_markdown,
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
    from r_dop_import_metadata_lib import get_download_urls_and_names
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

# set global varibales
ID = grass.tempname(12)
ORIG_REGION = f"original_region_{ID}"
rm_rasters = []
SUPPORTED = OPEN_DATA_AVAILABILITY["SUPPORTED"]
NO_OPEN_DATA = OPEN_DATA_AVAILABILITY["NO_OPEN_DATA"]
NOT_YET_SUPPORTED = OPEN_DATA_AVAILABILITY["NOT_YET_SUPPORTED"]

# Band suffixes used by DOP imports
DOP_BAND_SUFFIXES = ("_red", "_green", "_blue", "_nir")


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


def get_addon_name(fs):
    """Function to get the addon name for the function to get license info"""
    if fs in ["BB", "BE"]:
        return "r.dop.import.bb.be"
    return f"r.dop.import.{fs.lower()}"


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
                    f"{output}_{fs}{band}" for band in DOP_BAND_SUFFIXES
                ]
                local_fs_dir = os.path.join(local_data_dir, fs)
                if pathlib.Path(local_fs_dir).exists():
                    for _root, _dirs, files in os.walk(local_fs_dir):
                        dop_names.extend(
                            file
                            for file in files
                            if file.lower().endswith(
                                (".tif", ".tiff", ".jp2", ".jpeg"),
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

                # Only create a tempfile for URL/metadata exchange with the
                # state-specific addon if a metadata file was actually
                # requested by the user
                metadata_tmpfile = None
                if metadata_path:
                    metadata_tmpfile = grass.tempfile()
                    params["metadata_file"] = metadata_tmpfile

                # Run addon
                grass.run_command(addon, **params)

                if metadata_tmpfile:
                    # Reads URLs from the tempfile written by the addon, with
                    # fallbacks to the download directory and raster count
                    # if no URLs could be determined (see
                    # r_dop_import_metadata_lib.py)
                    dop_urls, dop_names = get_download_urls_and_names(
                        metadata_tmpfile=metadata_tmpfile,
                        keep_data=keep_data,
                        download_dir=download_dir,
                        out_fs=out_fs,
                        fs=fs,
                    )

            for band in ("red", "green", "blue", "nir"):
                all_dops[band].append(f"{out_fs}_{band}")
            fs_dop_list = [
                f"{out_fs}_{band}" for band in ("red", "green", "blue", "nir")
            ]

        # Collect metadata for this federal state (license/source info comes
        # from the addon's HTML documentation, file/URL info from above)
        addon_name = get_addon_name(fs)
        license_info, base_url = get_license_and_url_from_addon(addon_name)
        fs_metadata = collect_metadata(
            fs=fs,
            raster_list=fs_dop_list,
            license_info=license_info,
            base_url=base_url,
            original_names=dop_names,
            download_urls=dop_urls,
            band_suffixes=DOP_BAND_SUFFIXES,
        )
        metadata_list.append(fs_metadata)

    # Build one VRT per band across all federal states, then group them into
    # the final multi-band output
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

    # Write metadata file if metadata_path was set
    write_metadata_markdown(
        metadata_list=metadata_list,
        metadata_path=metadata_path,
        data_label="DOP",
    )

    grass.message(_(f"DOP group <{output}> is created."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
