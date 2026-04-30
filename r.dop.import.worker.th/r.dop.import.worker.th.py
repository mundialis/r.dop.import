#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.th
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Leon Louwarts
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              in Thüringen
# SPDX-FileCopyrightText: (c) 2024-2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
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
# % key: layer_name_cir
# % required: yes
# % description: Layer name of CIR tile-DOP to import
# %end

# %option
# % key: layer_name_rgb
# % required: yes
# % description: Layer name of RGB tile-DOP to import
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

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.location import switch_back_original_location
from grass_gis_helpers.mapset import switch_to_new_mapset
from grass_gis_helpers.raster import adjust_raster_resolution

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the dop library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import rescale_to_1_255, import_dop_from_wms
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

rm_rast = []
rm_group = []

original_nprocs = None

RETRIES = 30
WAITING_TIME = 10


def cleanup():
    """Remove all not needed files at the end"""
    general_cleanup(
        rm_rasters=rm_rast,
        rm_groups=rm_group,
    )
    """Reset nprocs"""
    if original_nprocs:
        grass.run_command("g.gisenv", set=f"NPROCS={original_nprocs}")


def main():
    """Main function of r.dop.import.worker.th"""
    # parser options
    tile_key = options["tile_key"]
    tile_url = options["tile_url"]
    layer_name_cir = options["layer_name_cir"]
    layer_name_rgb = options["layer_name_rgb"]
    raster_name = options["raster_name"]
    resolution_to_import = None
    if options["resolution_to_import"]:
        resolution_to_import = float(options["resolution_to_import"])
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]

    # check number of nprocs used and set to 1, write original value in variable
    try:
        original_nprocs = grass.read_command("g.gisenv", get="NPROCS").strip()
        if int(original_nprocs) > 1:
            grass.run_command("g.gisenv", set="NPROCS=1")
    except (ValueError, AttributeError):
        original_nprocs = None
        grass.run_command("g.gisenv", set="NPROCS=1")

    # output resolution
    if not flags["r"] and not options["resolution_to_import"]:
        grass.fatal(
            "Use native resolution with the -r flag or specify "
            "'resolution_to_import'.",
        )

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command("g.region", region=f"{orig_region}@{old_mapset}")

    # import DOP tile with original resolution
    grass.message(
        _(f"Started DOP import for key: {tile_key} and URL: {tile_url}"),
    )

    # import DOPs from WMS
    import_dop_from_wms(
        f"{tile_key}@{old_mapset}",
        raster_name,
        {"cir": tile_url, "rgb": tile_url},
        resolution_to_import,
        {"cir": layer_name_cir, "rgb": layer_name_rgb},
        rm_group,
        rm_rast,
        flags["r"],
    )

    # adjust resolution if required
    if resolution_to_import:
        grass.run_command("g.region", res=resolution_to_import, flags="a")
        for band in [1, 2, 3, 4]:
            raster_name_band = f"{raster_name}.{band}"
            grass.run_command(
                "g.rename",
                raster=f"{raster_name_band},{raster_name_band}_tmp1",
            )
            adjust_raster_resolution(
                f"{raster_name_band}_tmp1",
                raster_name_band,
                resolution_to_import,
                type="CELL",
            )
            rm_rast.append(f"{raster_name_band}_tmp1")

    rm_group.append(raster_name)
    grass.message(_(f"Finishing raster import for {raster_name}..."))

    # rescale imported DOPs
    new_rm_rast = rescale_to_1_255("TH", raster_name, extension="num")
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
