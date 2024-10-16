#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.rp
# AUTHOR(S):   Victoria-Leandra Brunn
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              in Rheinland-Pfalz
# COPYRIGHT:   (C) 2024 by mundialis GmbH & Co. KG and the GRASS Development
#              Team
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
# % description: Downloads and imports single Digital Orthophotos (DOPs) in Rhineland-Palatine
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
# % required: yes
# % description: Resolution of region, for which DOP will be imported
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
import sys
import ssl
import os
import urllib.request

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup, cleaning_tmp_location
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.location import switch_back_original_location
from grass_gis_helpers.mapset import switch_to_new_mapset
from grass_gis_helpers.raster import adjust_raster_resolution

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the dop library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import rescale_to_1_256, import_and_reproject
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

rm_rast = []
rm_group = []

gisdbase = None
TMP_LOC = None
TMP_GISRC = None


def cleanup():
    """Remove all not needed files at the end"""
    cleaning_tmp_location(
        None,
        tmp_loc=TMP_LOC,
        tmp_gisrc=TMP_GISRC,
        gisdbase=gisdbase,
    )
    general_cleanup(
        rm_rasters=rm_rast,
        rm_groups=rm_group,
    )


def main():
    """Main function of r.dop.import.worker.rp"""
    global gisdbase, TMP_LOC, TMP_GISRC
    # parser options
    tile_key = options["tile_key"]
    tile_url = options["tile_url"]
    raster_name = options["raster_name"]
    resolution_to_import = float(options["resolution_to_import"])
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]
    download_dir = options["download_dir"]
    keep_data = flags["k"]

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command("g.region", region=f"{orig_region}@{old_mapset}")
    aoi_map = f"{options['aoi']}@{old_mapset}" if options["aoi"] else None

    # import DOP tile with original resolution
    grass.message(
        _(f"Started DOP import for key: {tile_key} and URL: {tile_url}"),
    )

    # sth has to happen here as import from website is not possible with r.import used in import_and_reproject
    dir_tmp = grass.tempdir()
    # get around security certificate of geoportal-rlp

    ssl._create_default_https_context = ssl._create_unverified_context
    tile_path = os.path.join(dir_tmp, raster_name + ".jp2")
    urllib.request.urlretrieve(tile_url, tile_path)
    grass.message(
        (f"DOP download to DIR: {dir_tmp} and Raster: {raster_name}"),
    )

    # import and reproject DOP tiles based on tileindex
    gisdbase, TMP_LOC, TMP_GISRC = import_and_reproject(
        tile_path,  # tile_url,
        raster_name,
        resolution_to_import,
        "RP",
        aoi_map,
        download_dir,
        epsg=25832,
        keep_data=keep_data,
    )

    # adjust resolution if required
    for band in [1, 2, 3, 4]:
        raster_name_band = f"{raster_name}.{band}"
        grass.run_command(
            "g.region",
            raster=raster_name_band,
            res=resolution_to_import,
            flags="a",
        )
        adjust_raster_resolution(
            raster_name_band,
            raster_name_band,
            resolution_to_import,
        )
    rm_group.append(raster_name)
    grass.message(_(f"Finishing raster import for {raster_name}..."))

    # rescale imported DOPs
    new_rm_rast = rescale_to_1_256("RP", raster_name)
    rm_rast.extend(new_rm_rast)

    # switch back to original location
    switch_back_original_location(gisrc)
    grass.utils.try_remove(newgisrc)
    grass.message(
        _(f"DOP import for key: {tile_key} and URL: {tile_url} done!"),
    )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
