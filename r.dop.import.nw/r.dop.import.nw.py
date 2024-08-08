#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.nw
# AUTHOR(S):   Johannes Halbauer, Anika Weinmann
#
# PURPOSE:     Downloads DOPs for Nordrhein-Westfalen and AOI
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

# %module
# % description: Downloads DOPs for Nordrhein-Westfalen and AOI.
# % keyword: raster
# % keyword: import
# % keyword: DOP
# % keyword: open-geodata-germany
# %end

# %option G_OPT_V_INPUT
# % key: aoi
# % description: Polygon of the area of interest to set region
# % required: no
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
# % label: Keep downloaded data in the download directory
# %end

# %flag
# % key: r
# % label: Use native data resolution
# %end

# %rules
# % requires_all: -k,download_dir
# %end

import atexit
import os
import sys

import grass.script as grass
from grass.pygrass.modules import Module, ParallelModuleQueue
from grass.pygrass.utils import get_lib_path
from osgeo import gdal

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.data_import import (
    download_and_import_tindex,
    get_list_of_tindex_locations,
)
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
)
from grass_gis_helpers.raster import create_vrt

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the dop library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import setup_parallel_processing
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

# set global variables
TINDEX = (
    "https://github.com/mundialis/tile-indices/raw/main/DOP/NW/"
    "openNRW_DOP10_tileindex.gpkg.gz"
)

ID = grass.tempname(12)
orig_region = None
rm_rasters = []
rm_vectors = []
download_dir = None
rm_dirs = []


def cleanup():
    general_cleanup(
        orig_region=orig_region,
        rm_rasters=rm_rasters,
        rm_vectors=rm_vectors,
        rm_dirs=rm_dirs,
    )


def main():
    global ID, orig_region, rm_rasters, rm_vectors
    global temp_loc_path, download_dir, new_mapset

    aoi = options["aoi"]
    download_dir = check_download_dir(options["download_dir"])
    nprocs = int(options["nprocs"])
    nprocs = setup_parallel_processing(nprocs)
    output = options["output"]
    fs = "NW"

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # create list for each raster band for building entire raster
    all_raster = {
        "red": [],
        "green": [],
        "blue": [],
        "nir": [],
    }

    # save original region
    orig_region = f"original_region_{ID}"
    grass.run_command("g.region", save=orig_region, quiet=True)

    # get region resolution and check if resolution consistent
    reg = grass.region()
    if reg["nsres"] == reg["ewres"]:
        ns_res = reg["nsres"]
    else:
        grass.fatal("N/S resolution is not the same as E/W resolution!")

    # set region if aoi is given
    if aoi:
        grass.run_command("g.region", vector=aoi, flags="a")
    # if no aoi save region as aoi
    else:
        aoi = f"region_aoi_{ID}"
        grass.run_command(
            "v.in.region",
            output=aoi,
            quiet=True,
        )

    # get tile index
    tindex_vect = f"dop_tindex_{ID}"
    rm_vectors.append(tindex_vect)
    download_and_import_tindex(TINDEX, tindex_vect, download_dir)

    # get download urls which overlap with AOI
    # or current region if no AOI is given
    url_tiles = get_list_of_tindex_locations(tindex_vect, aoi)

    url_tiles_count = 0
    for i in range(len(url_tiles)):
        url_tiles_count += 1
        url_tiles[i] = (url_tiles_count, [url_tiles[i]])
    number_tiles = len(url_tiles)

    # set number of parallel processes to number of tiles
    if number_tiles < nprocs:
        nprocs = number_tiles
    queue = ParallelModuleQueue(nprocs=nprocs)

    # get GISDBASE and Location
    gisenv = grass.gisenv()
    gisdbase = gisenv["GISDBASE"]
    location = gisenv["LOCATION_NAME"]

    # set queue and variables for worker adddon
    try:
        grass.message(
            _(f"Importing {number_tiles} DOPs for NW in parallel...")
        )
        for tile in url_tiles:
            key = tile[0]
            new_mapset = f"tmp_mapset_rdop_import_tile_{key}_{os.getpid()}"
            rm_dirs.append(os.path.join(gisdbase, location, new_mapset))
            b_name = os.path.basename(tile[1][0])
            raster_name = (
                f"{b_name.split('.')[0].replace('-', '_')}" f"_{os.getpid()}"
            )
            for key_rast in all_raster:
                all_raster[key_rast].append(
                    f"{fs}_{raster_name}_{key_rast}@{new_mapset}"
                )
            param = {
                "tile_key": key,
                "tile_url": tile[1][0],
                "raster_name": raster_name,
                "orig_region": orig_region,
                "memory": 1000,
                "new_mapset": new_mapset,
                "flags": "",
            }
            grass.message(_(f"raster name: {raster_name}"))

            # modify params
            if aoi:
                param["aoi"] = aoi
            if options["download_dir"]:
                param["download_dir"] = download_dir
            if flags["k"]:
                param["flags"] += "k"
            if flags["r"]:
                dop_src = gdal.Open(param["tile_url"])
                param["resolution_to_import"] = abs(
                    dop_src.GetGeoTransform()[1]
                )
            else:
                param["resolution_to_import"] = ns_res

            # append raster bands to download to remove list
            rm_red = f"{raster_name}_red"
            rm_green = f"{raster_name}_green"
            rm_blue = f"{raster_name}_blue"
            rm_nir = f"{raster_name}_nir"
            rm_rasters.append(rm_red)
            rm_rasters.append(rm_green)
            rm_rasters.append(rm_blue)
            rm_rasters.append(rm_nir)

            # run worker addon in parallel
            r_dop_import_worker_NW = Module(
                "r.dop.import.worker.nw",
                **param,
                run_=False,
            )
            # catch all GRASS output to stdout and stderr
            r_dop_import_worker_NW.stdout = grass.PIPE
            r_dop_import_worker_NW.stderr = grass.PIPE
            queue.put(r_dop_import_worker_NW)
        queue.wait()
    except Exception:
        for proc_num in range(queue.get_num_run_procs()):
            proc = queue.get(proc_num)
            if proc.returncode != 0:
                # save all stderr to a variable and pass it to a GRASS
                # exception
                errmsg = proc.outputs["stderr"].value.strip()
                grass.fatal(
                    _(f"\nERROR by processing <{proc.get_bash()}>: {errmsg}")
                )

    # create one vrt per band of all imported DOPs
    raster_out = []
    for band, b_list in all_raster.items():
        out = f"{output}_{band}"
        create_vrt(b_list, out)
        raster_out.append(out)

    grass.message(_(f"Generated following raster maps: {raster_out}"))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
