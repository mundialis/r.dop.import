#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.he
# AUTHOR(S):   Johannes Halbauer, Anika Weinmann
#
# PURPOSE:     Downloads DOPs for Hessen and AOI
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
# % description: Downloads DOPs for Hessen and AOI.
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

from grass_gis_helpers.cleanup import general_cleanup
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

ID = grass.tempname(12)
orig_region = None
rm_rasters = []
rm_vectors = []
download_dir = None
rm_dirs = []

WMS = "https://www.gds-srv.hessen.de/cgi-bin/lika-services/ogc-free-images.ows?language=ger&"
NATIVE_DOP_RES = 0.2


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
    fs = "HE"

    # if -k flag is set print warning that it will be ignored because
    # the data will be directly imported into GRASS from WMS
    if flags["k"]:
        grass.warning(
            _(
                "-k flag will be ignored, beacuse TH DOPs will be imported directly from WMS into GRASS. "
                "Use r.out.gdal module to export DOPs into download directory!"
            )
        )

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
        ew_res = reg["ewres"]
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

    # create grid for downloading
    grass.message(_("Creating DOP tiles for HE..."))

    # set tile size in map units (meter)
    tile_size = 1000

    # set grid name
    grid = f"tmp_grid_HE_{ID}"

    # check if aoi is smaller than tile size
    if ns_res <= float(tile_size) and ew_res <= float(tile_size):
        grass.run_command("v.in.region", output=grid, quiet=True)
        rm_vectors.append(grid)
        grass.run_command(
            "v.db.addtable", map=grid, columns="cat int", quiet=True
        )
    else:
        grass.run_command("g.region", res=tile_size, flags="a", quiet=True)

        # create grid
        grass.run_command(
            "v.mkgrid",
            map=grid,
            box=f"{tile_size},{tile_size}",
            quiet=True,
        )
        # reset region
        grass.run_command("g.region", vector=aoi)

    # set grid name
    grid_name = f"tmp_grid_area_{ID}"

    # choose tiles overlapping with aoi
    grass.run_command(
        "v.select",
        ainput=grid,
        binput=aoi,
        output=grid_name,
        operator="overlap",
        quiet=True,
    )
    rm_vectors.append(grid_name)

    # create list of tiles
    tiles_num_list = list(
        grass.parse_command(
            "v.db.select", map=grid_name, columns="cat", flags="c", quiet=True
        ).keys()
    )
    number_tiles = len(tiles_num_list)

    grass.message(_(f"Number of tiles: {number_tiles}"))
    tiles_list = []
    for tile in tiles_num_list:
        tile_area = f"HE_DOP_{tile}"
        grass.run_command(
            "v.extract",
            input=grid_name,
            where=f"cat == {tile}",
            output=tile_area,
            quiet=True,
        )
        tiles_list.append(tile_area)
        rm_vectors.append(tile_area)

    # set number of parallel processes to number of tiles
    if number_tiles < nprocs:
        nprocs = number_tiles
    queue = ParallelModuleQueue(nprocs=nprocs)

    # get GISDBASE and Location
    gisenv = grass.gisenv()
    gisdbase = gisenv["GISDBASE"]
    location = gisenv["LOCATION_NAME"]

    # set queue and variables for worker addon
    try:
        grass.message(
            _(f"Importing {number_tiles} DOPs for HE in parallel...")
        )
        for tile in tiles_list:
            key = tile
            new_mapset = f"tmp_mapset_rdop_import_tile_{key}_{os.getpid()}"
            rm_dirs.append(os.path.join(gisdbase, location, new_mapset))
            raster_name = tile
            for key_rast in all_raster:
                all_raster[key_rast].append(
                    f"{fs}_{raster_name}_{key_rast}@{new_mapset}"
                )
            param = {
                "flags": "",
                "tile_key": key,
                "tile_url": WMS,
                "raster_name": raster_name,
                "orig_region": orig_region,
                "memory": 1000,
                "new_mapset": new_mapset,
            }
            grass.message(f"raster_name: {raster_name}")
            if aoi:
                param["aoi"] = aoi
            if options["download_dir"]:
                param["download_dir"] = download_dir
            if flags["k"]:
                param["flags"] += "k"
            if flags["r"]:
                param["resolution_to_import"] = NATIVE_DOP_RES
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
            r_dop_import_worker_HE = Module(
                "r.dop.import.worker.he",
                **param,
                run_=False,
            )
            # catch all GRASS output to stdout and stderr
            r_dop_import_worker_HE.stdout = grass.PIPE
            r_dop_import_worker_HE.stderr = grass.PIPE
            queue.put(r_dop_import_worker_HE)
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
