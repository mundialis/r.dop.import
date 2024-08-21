#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Anika Weinmann, Julia Haas
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) wihtin a specified
#              federal state and area of interest
# COPYRIGHT:   (C) 2024 by mundialis GmbH & Co. KG and the GRASS
#              Development Team
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
############################################################################

# %Module
# % description: Downloads and imports DOPs for specified federal state and aoi
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
# % options: Berlin,BE,Brandenburg,BB,Nordrhein-Westfalen,NW,Sachsen,SN,Thüringen,TH
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
        band_dict=None,
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


def main():
    """Main function of r.dop.import"""
    aoi = options["aoi"]
    federal_states = get_federal_states(
        options["federal_state"],
        options["federal_state_file"],
    )
    local_data_dir = options["local_data_dir"]
    download_dir = check_download_dir(options["download_dir"])
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
    for fs in set(federal_states):
        grass.message(_(f"Importing DOPs for {fs}..."))
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
                if grass.find_program(worker_addon, "--help"):
                    params["nprocs"] = nprocs
                grass.run_command(addon, **params)
            all_dops["red"].append(f"{out_fs}_red")
            all_dops["green"].append(f"{out_fs}_green")
            all_dops["blue"].append(f"{out_fs}_blue")
            all_dops["nir"].append(f"{out_fs}_nir")

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
    grass.message(_(f"DOP group <{output}> is created."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
