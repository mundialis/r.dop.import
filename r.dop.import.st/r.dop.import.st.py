#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.st
# AUTHOR(S):   Johannes Halbauer, Anika Weinmann, Leon Louwarts
# PURPOSE:     Downloads DOPs for Sachsen-Anhalt and AOI
# SPDX-FileCopyrightText: (c) 2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
############################################################################

# %module
# % description: Downloads DOPs for Sachsen-Anhalt and AOI.
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
# % label: Path to output folder
# % description: Path to download folder
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

# %option
# % key: metadata_file
# % type: string
# % required: no
# % description: Temporary file for metadata URLs
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
import pathlib

import grass.script as grass
from grass.pygrass.modules import Module, ParallelModuleQueue
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.data_import import (
    download_and_import_tindex,
    get_list_of_tindex_locations,
)
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
)
from grass_gis_helpers.raster import (
    adjust_raster_resolution,
    create_vrt,
    vrt_to_raster,
)

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
    "https://github.com/mundialis/tile-indices/raw/main/DOP/ST/"
    "st_dop_tindex_proj.gpkg.gz"
)
NATIVE_DOP_RES = 0.2

ID = grass.tempname(12)
ORIG_REGION = f"original_region_{ID}"
rm_rasters = []
rm_vectors = []
download_dir = None
rm_dirs = []


def cleanup():
    """Remove all not needed files at the end"""
    general_cleanup(
        orig_region=ORIG_REGION,
        rm_rasters=rm_rasters,
        rm_vectors=rm_vectors,
        rm_dirs=rm_dirs,
    )


def main():
    """Main function of r.dop.import.st"""
    aoi = options["aoi"]
    download_dir = check_download_dir(options["download_dir"])
    nprocs = int(options["nprocs"])
    nprocs = setup_parallel_processing(nprocs)
    metadata_file = options["metadata_file"]
    output = options["output"]
    fs = "ST"

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
    grass.run_command("g.region", save=ORIG_REGION, quiet=True)

    # get region resolution and check if resolution consistent
    reg = grass.region()
    if reg["nsres"] == reg["ewres"]:
        ns_res = reg["nsres"]
    else:
        grass.fatal("N/S resolution is not the same as E/W resolution!")

    # set region if aoi is given
    if aoi:
        grass.run_command("g.region", vector=aoi, res=ns_res, flags="a")
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

    # get item_ids which overlap with AOI (or region if no AOI given)
    id_tiles = get_list_of_tindex_locations(
        tindex_vect,
        aoi,
        column="item_id",
    )
    id_tiles = [(i + 1, item_id) for i, item_id in enumerate(id_tiles)]
    number_tiles = len(id_tiles)

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
            _(f"Importing {number_tiles} DOPs for ST in parallel..."),
        )
        for tile_key, item_id in id_tiles:
            new_mapset = (
                f"tmp_mapset_rdop_import_tile_{tile_key}_{os.getpid()}"
            )
            rm_dirs.append(os.path.join(gisdbase, location, new_mapset))
            # b_name = parse_qs(urlparse(tile[1][0]).query)["file"][0]
            raster_name = f"dop20_{item_id}_{os.getpid()}"
            for item in all_raster.items():
                item[1].append(f"{fs}_{raster_name}_{item[0]}@{new_mapset}")
            param = {
                "tile_key": tile_key,
                "item_id": item_id,
                "raster_name": raster_name,
                "orig_region": ORIG_REGION,
                "memory": 1000,
                "new_mapset": new_mapset,
                "resolution_to_import": NATIVE_DOP_RES,
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

            rm_rasters.extend(
                f"{fs}_{raster_name}{band}"
                for band in ("red", "green", "blue", "nir")
            )
            # for band in ("red", "green", "blue", "nir"):
            #     rm_rasters.append(f"{fs}_{raster_name}{band}")

            # run worker addon in parallel
            r_dop_import_worker_st = Module(
                "r.dop.import.worker.st",
                **param,
                run_=False,
            )
            # catch all GRASS output to stdout and stderr
            r_dop_import_worker_st.stdout_ = grass.PIPE
            r_dop_import_worker_st.stderr_ = grass.PIPE
            queue.put(r_dop_import_worker_st)
        queue.wait()
    except Exception:
        for proc_num in range(queue.get_num_run_procs()):
            proc = queue.get(proc_num)
            if proc.returncode != 0:
                # save all stderr to a variable and pass it to a GRASS
                # exception
                errmsg = proc.outputs["stderr"].value.strip()
                grass.fatal(
                    _(f"\nERROR by processing <{proc.get_bash()}>: {errmsg}"),
                )

    if metadata_file:
        try:
            with pathlib.Path(metadata_file).open("w", encoding="utf-8") as f:
                f.writelines(f"{item_id}\n" for _, item_id in id_tiles)

            grass.debug(f"Wrote {len(id_tiles)} URLs to tempfile")

        except Exception as e:
            grass.warning(f"Could not write tempfile metadata: {e}")

    # create one vrt per band of all imported DOPs
    raster_out = []
    for band, b_list in all_raster.items():
        vrt = f"vrt_{output}_{band}_{ID}"
        rm_rasters.append(vrt)
        rm_rasters.extend([r.split("@")[0] for r in b_list])
        create_vrt(b_list, vrt)

        out_band = f"{output}_{band}"
        if flags["r"]:
            # Note: Want real raster/no VRT as output
            vrt_to_raster(vrt, out_band)
        else:
            grass.message(_(f"Resampling / interpolating {band} band..."))
            grass.run_command("g.region", raster=vrt)
            grass.run_command("g.region", res=ns_res, flags="a")
            adjust_raster_resolution(vrt, out_band, ns_res)
        raster_out.append(out_band)

    if not raster_out:
        grass.fatal("No output rasters could be created")

    grass.message(_(f"Generated following raster maps: {raster_out}"))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
