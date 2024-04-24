#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.bb
# AUTHOR(S):   Johannes Halbauer & Lina Krisztian
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area in Brandenburg
# COPYRIGHT:   (C) 2024 by mundialis GmbH & Co. KG and the GRASS Development Team
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
#############################################################################

# %Module
# % description: Downloads and imports single Digital Orthophotos (DOPs) in Brandenburg
# % keyword: imagery
# % keyword: download
# % keyword: DOP
# %end

# %option G_OPT_V_INPUT
# % key: aoi
# % required: no
# % description: Vector map to restrict DOP import to
# %end

# %option
# % key: download_dir
# % label: Path to output folder
# % description: Path to download folder
# % required: no
# % multiple: no
# %end

# %option
# % key: tile_key
# % required: yes
# % description: Key of tile-DOP to import
# %end

# %option
# % key: tile_url
# % required: yes
# % description: URL of tile-DOP to import
# %end

# %option
# % key: new_mapset
# % type: string
# % required: yes
# % multiple: no
# % key_desc: name
# % description: Name for new mapset
# %end

# %option
# % key: orig_region
# % required: yes
# % description: Original region
# %end

# %option
# % key: resolution_to_import
# % required: no
# % description: Resolution of region, for which DOP will be imported (only if flag r not set)
# %end

# %option G_OPT_R_OUTPUT
# % key: raster_name
# % description: Name of raster output
# %end

# %option G_OPT_MEMORYMB
# % description: Memory which is used by all processes (it is divided by nprocs for each single parallel process)
# %end

# %flag
# % key: r
# % description: Use native DOP resolution
# %end

# %flag
# % key: k
# % label: Keep downloaded data in the download directory
# %end

# %rules
# % requires_all: -k,download_dir
# %end


import atexit
import os
import sys
from time import sleep

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup, cleaning_tmp_location
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.location import (
    get_current_location,
    create_tmp_location,
    switch_back_original_location,
)
from grass_gis_helpers.mapset import switch_to_new_mapset
from grass_gis_helpers.open_geodata_germany.download_data import (
    download_data_using_threadpool,
    extract_compressed_files,
)
from grass_gis_helpers.raster import adjust_raster_resolution

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the dop library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import rescale_to_1_256
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

rm_rast = []
rm_group = []

gisdbase = None
tmp_loc = None
tmp_gisrc = None

RETRIES = 30
WAITING_TIME = 10


def cleanup():
    cleaning_tmp_location(
        None, tmp_loc=tmp_loc, tmp_gisrc=tmp_gisrc, gisdbase=gisdbase
    )
    general_cleanup(
        rm_rasters=rm_rast,
        rm_groups=rm_group,
    )


def import_and_reproject(
    url,
    raster_name,
    resolution_to_import,
    aoi_map=None,
    download_dir=None,
    epsg=None,
):
    """Import DOPs and reproject them if needed. This is needed for the DOPs
    of Brandenburg and Berlin because GDAL (at least smaller 3.6.3) does not
    support the coordinate reference system in the data.

    Args:
        url (str): The URL of the DOP to import
        raster_name (str): The prefix name for the output rasters
        aoi_map (str): Name of AOI vector map
        download_dir (str): Path to local directory to downlaod DOPs to
        epsg (int): EPSG code which has to be set if the reproduction should be
                    done manually and not by r.import
    """
    global gisdbase, tmp_loc, tmp_gisrc
    aoi_map_to_set_region1 = aoi_map

    # get actual location, mapset, ...
    loc, mapset, gisdbase, gisrc = get_current_location()
    if not aoi_map:
        aoi_map = f"region_aoi_{grass.tempname(8)}"
        aoi_map_to_set_region1 = aoi_map
        grass.run_command("v.in.region", output=aoi_map, quiet=True)
    else:
        aoi_map_mapset = aoi_map.split("@")
        aoi_map_to_set_region1 = aoi_map_mapset[0]
        if len(aoi_map_mapset) == 2:
            mapset = aoi_map_mapset[1]

    # create temporary location with EPSG:25833
    tmp_loc, tmp_gisrc = create_tmp_location(epsg)

    # reproject aoi
    if aoi_map:
        grass.run_command(
            "v.proj",
            location=loc,
            mapset=mapset,
            input=aoi_map_to_set_region1,
            output=aoi_map_to_set_region1,
            quiet=True,
        )
        grass.run_command(
            "g.region",
            vector=aoi_map_to_set_region1,
            res=resolution_to_import,
            flags="a",
        )

    # import data
    # set memory manually to 1000
    # Process stuck, when memory is too large (100000)
    # GDAL_CACHEMAX wird nur als MB interpretiert
    kwargs = {
        "input": url,
        "output": raster_name,
        "memory": 1000,
        "quiet": True,
        "flags": "o",
        "extent": "region",
    }
    # download and keep data to download dir if -k flag ist set
    # and change input parameter in kwargs to local path
    if flags["k"]:
        url = url.replace("/vsizip/vsicurl/", "")
        basename = os.path.basename(url)
        url = url.replace(basename, "")[:-1]
        download_data_using_threadpool([url], download_dir, 1)
        extract_compressed_files(
            [basename.replace(".tif", ".zip")], download_dir
        )
        kwargs["input"] = os.path.join(download_dir, basename)

    import_sucess = False
    tries = 0
    while not import_sucess:
        tries += 1
        if tries > RETRIES:
            grass.fatal(
                _(
                    f"Importing {kwargs['input']} failed after {RETRIES} "
                    "retries."
                )
            )
        try:
            grass.run_command("r.import", **kwargs)
            import_sucess = True
        except Exception:
            sleep(WAITING_TIME)
    if not aoi_map:
        grass.run_command("g.region", raster=f"{raster_name}.1")

    # reproject data
    res = float(
        grass.parse_command("r.info", flags="g", map=f"{raster_name}.1")[
            "nsres"
        ]
    )
    # switch location
    os.environ["GISRC"] = str(gisrc)
    if aoi_map:
        grass.run_command("g.region", vector=aoi_map, res=res, flags="a")
    else:
        grass.run_command("g.region", res=res, flags="a")
    for i in range(1, 5):
        name = f"{raster_name}.{i}"
        # set memory manually to 1000
        # Process stuck, when memory is too large (100000)
        # GDAL_CACHEMAX is only interpreted as MB, if value is <100000
        grass.run_command(
            "r.proj",
            location=tmp_loc,
            mapset="PERMANENT",
            input=name,
            output=name,
            resolution=res,
            flags="n",
            quiet=True,
            memory=1000,
        )


def main():
    # parser options
    tile_key = options["tile_key"]
    tile_url = options["tile_url"]
    raster_name = options["raster_name"]
    resolution_to_import = float(options["resolution_to_import"])
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]
    download_dir = options["download_dir"]

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command("g.region", region=f"{orig_region}@{old_mapset}")
    if options["aoi"]:
        aoi_map = f"{options['aoi']}@{old_mapset}"
    else:
        aoi_map = None

    # import DOP tile with original resolution
    grass.message(
        _(f"Started DOP import for key: {tile_key} and URL: {tile_url}")
    )

    # import and reproject DOP tiles based on tileindex
    import_and_reproject(
        tile_url,
        raster_name,
        resolution_to_import,
        aoi_map,
        download_dir,
        epsg=25833,
    )

    # adjust resolution if required
    if not flags["r"]:
        for band in [1, 2, 3, 4]:
            raster_name_band = f"{raster_name}.{band}"
            adjust_raster_resolution(
                raster_name_band, raster_name_band, resolution_to_import
            )
    rm_group.append(raster_name)
    grass.message(_(f"Finishing raster import for {raster_name}..."))

    # rescale imported DOPs
    new_rm_rast = rescale_to_1_256("BB_BE", raster_name)
    rm_rast.extend(new_rm_rast)

    # switch back to original location
    switch_back_original_location(gisrc)
    grass.utils.try_remove(newgisrc)
    grass.message(
        _(f"DOP import for key: {tile_key} and URL: {tile_url} done!")
    )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
