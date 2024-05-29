#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.th
# AUTHOR(S):   Johannes Halbauer & Lina Krisztian
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area in Thüringen
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
# % description: Downloads and imports single Digital Orthophotos (DOPs) in Thüringen
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

# %flag
# % key: r
# % description: Use native DOP resolution
# %end


import atexit
import sys
from time import sleep

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup, cleaning_tmp_location
from grass_gis_helpers.location import switch_back_original_location
from grass_gis_helpers.mapset import switch_to_new_mapset
from grass_gis_helpers.raster import adjust_raster_resolution, rename_raster

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


def import_dop_from_th_wms(
    tile_key, rastername, tile_url, resolution_to_import
):
    """Import DOP from TH WMS

    Args:
        tile_key (str): Key of current tile
        rastername (str): Name of resulting raster
        tile_url (str): WMS URL to get DOPs from
        resolution_to_import (float): Resolution to resample imported raster to
    """
    # set region and create variable names
    grass.run_command("g.region", vector=tile_key)
    tile_key = tile_key.split("@")[0]
    for name in ["20cir", ""]:
        if name == "20cir":
            out_tmp = f"{name}_{tile_key}_tmp"
            bands = ["red"]
        else:
            out_tmp = f"{tile_key}_tmp"
            bands = ["red", "green", "blue"]
        rm_group.append(out_tmp)
        for band in ["red", "green", "blue"]:
            rm_rast.append(f"{out_tmp}.{band}")

        # import wms data and retry download if wms fails 15 times
        trydownload = True
        count = 0
        while trydownload:
            try:
                count += 1
                grass.run_command(
                    "r.in.wms",
                    url=tile_url,
                    layer=f"th_dop{name}",
                    output=out_tmp,
                    format="tiff",
                    flags="b",
                    overwrite=True,
                )
                trydownload = False
            except Exception:
                # remove maps where wms download failed
                grass.run_command(
                    "g.remove",
                    type="raster",
                    pattern=f"{out_tmp}.*",
                    flags="f",
                )
                grass.message(_("Retry download..."))
                if count > 15:
                    grass.fatal(f"Download of {tile_url} not working.")
                sleep(10)

        # change band name to band number
        for band in bands:
            if name == "20cir":
                oband = "4"
            else:
                oband = band
                if band == "red":
                    oband = "1"
                elif band == "green":
                    oband = "2"
                elif band == "blue":
                    oband = "3"

            # create old and new name for adjusting resolution/renaming
            old_name = f"{out_tmp}.{band}"
            new_name = f"{rastername}.{oband}"

            # adjust resolution/rename raster maps
            if not flags["r"]:
                adjust_raster_resolution(
                    old_name, new_name, resolution_to_import
                )
            else:
                rename_raster(old_name, new_name)
            rm_rast.append(old_name)

    # drop rest of CIR DOP
    grass.run_command(
        "g.remove", type="raster", pattern="cir_*_tmp*", flags="f"
    )


def main():
    # parser options
    tile_key = options["tile_key"]
    tile_url = options["tile_url"]
    raster_name = options["raster_name"]
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]
    if not flags["r"]:
        resolution_to_import = float(options["resolution_to_import"])
    else:
        resolution_to_import = None

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command(
        "g.region",
        region=f"{orig_region}@{old_mapset}",
        res=resolution_to_import,
    )

    # import DOP tile with original resolution
    grass.message(
        _(f"Started DOP import for key: {tile_key} and URL: {tile_url}")
    )
    # import DOPs from WMS
    import_dop_from_th_wms(
        f"{tile_key }@{old_mapset}",
        raster_name,
        tile_url,
        resolution_to_import,
    )

    grass.message(_(f"Finishing raster import for {raster_name}..."))

    # rescale imported DOPs
    new_rm_rast = rescale_to_1_256("TH", raster_name, extension="num")
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
