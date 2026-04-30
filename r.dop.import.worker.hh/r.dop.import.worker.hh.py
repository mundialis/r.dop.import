#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.hh
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Leon Louwarts
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              in Hamburg
# SPDX-FileCopyrightText: (c) 2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
#############################################################################

# %Module
# % description: Downloads and imports single Digital Orthophotos (DOPs) in Hamburg
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
# % key: tile_urls
# % required: yes
# % label: Input DOP tile URLs
# % description: Comma-separated list of URLs of DOP tiles to import; multiple tiles will be merged before import
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
import sys

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
    from r_dop_import_lib import rescale_to_1_255, import_and_reproject
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

rm_rast = []
rm_group = []

gisdbase = None
TMP_LOC = None
TMP_GISRC = None
original_nprocs = None


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
    """Reset nprocs"""
    if original_nprocs:
        grass.run_command("g.gisenv", set=f"NPROCS={original_nprocs}")


def main():
    """Main function of r.dop.import.worker.hh"""
    global gisdbase, TMP_LOC, TMP_GISRC
    # parser options
    tile_key = options["tile_key"]
    tile_urls = options["tile_urls"].split(",")
    raster_name = options["raster_name"]
    resolution_to_import = None
    if options["resolution_to_import"]:
        resolution_to_import = float(options["resolution_to_import"])
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]
    download_dir = options["download_dir"]
    keep_data = flags["k"]

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

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command("g.region", region=f"{orig_region}@{old_mapset}")
    aoi_map = f"{options['aoi']}@{old_mapset}" if options["aoi"] else None

    # import DOP tile with original resolution
    grass.message(
        _(
            f"Started DOP import for key: {tile_key} with {len(tile_urls)} "
            "URL(s).",
        ),
    )

    imported_rasters = []

    for i, url in enumerate(tile_urls):
        part_name = f"{raster_name}_part{i}"
        gisdbase, TMP_LOC, TMP_GISRC = import_and_reproject(
            url,
            part_name,
            resolution_to_import,
            "HH",
            aoi_map,
            download_dir,
            epsg=25832,
            keep_data=keep_data,
            retries=5,
        )

        # check if Band 1 exists (if no overlap: grass warning)
        test_raster = f"{part_name}.1"
        if grass.find_file(test_raster, element="cell")["name"]:
            imported_rasters.append(part_name)
            grass.message(
                _(
                    f"Successfully imported part {i + 1}/{len(tile_urls)}.",
                ),
            )
        else:
            grass.warning(
                _(
                    f"{url} {i + 1} could not be imported (no overlap).",
                ),
            )

    # merge in case of multiple DOPs
    if len(imported_rasters) > 1:
        grass.message(
            _(
                f"Merging {len(imported_rasters)} of {len(tile_urls)}"
                "raster parts for tile {tile_key}...",
            ),
        )
        for band in range(1, 5):
            band_inputs = [f"{part}.{band}" for part in imported_rasters]
            band_output = f"{raster_name}.{band}"
            grass.run_command(
                "r.patch",
                input=",".join(band_inputs),
                output=band_output,
            )
        # cleanup of single parts
        for part in imported_rasters:
            rm_group.extend(f"{part}.{band}" for band in range(1, 5))
    elif len(imported_rasters) == 1:
        for band in range(1, 5):
            band_input = f"{imported_rasters[0]}.{band}"
            band_output = f"{raster_name}.{band}"
            grass.run_command(
                "g.rename",
                raster=f"{band_input},{band_output}",
            )
    else:
        grass.fatal(
            _(
                f"Tile {tile_key}: None of the {len(tile_urls)} DOP(s)"
                "could be imported successfully.",
            ),
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
    new_rm_rast = rescale_to_1_255("HH", raster_name)
    rm_rast.extend(new_rm_rast)

    # switch back to original location
    switch_back_original_location(gisrc)
    grass.utils.try_remove(newgisrc)
    grass.message(
        _(f"DOP import for key: {tile_key} done!"),
    )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
