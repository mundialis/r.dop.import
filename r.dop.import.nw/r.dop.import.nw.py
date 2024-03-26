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
import grass.script as grass

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.open_geodata_germany.download_data import (
    check_download_dir,
    download_data_using_threadpool,
)
from grass_gis_helpers.raster import adjust_raster_resolution, create_vrt
from grass_gis_helpers.data_import import (
    download_and_import_tindex,
    get_list_of_tindex_locations,
)

# set global varibales
TINDEX = (
    "https://github.com/mundialis/tile-indices/raw/main/DOP/NW/"
    "openNRW_DOP10_tileindex.gpkg.gz"
)
DATA_BASE_URL = (
    "https://www.opengeodata.nrw.de/produkte/geobasis/lusat/akt/dop/"
    "dop_jp2_f10/"
)

ID = grass.tempname(12)
orig_region = None
keep_data = False
rm_rasters = []
rm_vectors = []
download_dir = None


def cleanup():
    rm_dirs = []
    if not keep_data:
        if download_dir:
            rm_dirs.append(download_dir)
    general_cleanup(
        orig_region=orig_region,
        rm_rasters=rm_rasters,
        rm_vectors=rm_vectors,
        rm_dirs=rm_dirs,
    )


def main():
    global ID, orig_region, rm_rasters, rm_vectors, keep_data, download_dir

    aoi = options["aoi"]
    download_dir = check_download_dir(options["download_dir"])
    output = options["output"]
    keep_data = flags["k"]
    native_res = flags["r"]

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # save original region
    orig_region = f"original_region_{ID}"
    grass.run_command("g.region", save=orig_region, quiet=True)
    ns_res = grass.region()["nsres"]

    # set region if aoi is given
    if aoi:
        grass.run_command("g.region", vector=aoi, flags="a")

    # get tile index
    tindex_vect = f"dop_tindex_{ID}"
    rm_vectors.append(tindex_vect)
    download_and_import_tindex(TINDEX, tindex_vect, download_dir)

    # get download urls which overlap with AOI
    # or current region if no AOI is given
    url_tiles = get_list_of_tindex_locations(tindex_vect, aoi)

    # if k flag is set download DOPs to download_dir
    if keep_data:
        download_list = [url.replace("/vsicurl/", "") for url in url_tiles]

        # download with using nprocs=3
        download_data_using_threadpool(download_list, download_dir, 3)

    # import DOPs directly
    grass.message(_("Importing DOPs..."))

    # create list of lists for all DOPs caused by GRASS import structure
    # (one raster map per band)
    all_dops = [[], [], [], []]
    for url in url_tiles:
        # define name for imported DOP
        dop_name = os.path.splitext(os.path.basename(url))[0].replace("-", "")

        # define input parameter for import
        if keep_data:
            input = os.path.join(download_dir, os.path.basename(url))
        else:
            if "/vsicurl/" not in url:
                input = f"/vsicurl/{url}"
            else:
                input = url

        # import DOPs
        grass.run_command(
            "r.import",
            input=input,
            output=dop_name,
            extent="region",
            overwrite=True,
            quiet=True,
        )
        # append DOP name with band suffix
        all_dops[0].append(dop_name + ".1")
        all_dops[1].append(dop_name + ".2")
        all_dops[2].append(dop_name + ".3")
        all_dops[3].append(dop_name + ".4")

    # create VRT
    vrt_outputs = []
    for dops, band in zip(all_dops, ["red", "green", "blue", "nir"]):
        vrt_output = f"{output}_{band}"
        create_vrt(dops, vrt_output)
        vrt_outputs.append(vrt_output)

    # resample / interpolate whole VRT (because interpolating single files leads
    # to empty rows and columns)
    # check resolution and resample / interpolate data if needed
    if not native_res:
        grass.message(_("Resampling / interpolating data..."))
        grass.run_command("g.region", raster=vrt_output, res=ns_res, flags="a")
        for vrt_out in vrt_outputs:
            adjust_raster_resolution(vrt_out, vrt_out, ns_res)

    for result in vrt_outputs:
        grass.message(_(f"DOP raster map <{result}> is created."))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
