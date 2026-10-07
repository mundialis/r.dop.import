#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import.worker.st
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Leon Louwarts
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              in Sachsen-Anhalt
# SPDX-FileCopyrightText: (c) 2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
#############################################################################

# %Module
# % description: Downloads and imports single Digital Orthophotos (DOPs) in Sachsen-Anhalt
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
# % key: item_id
# % required: yes
# % description: Numeric tile ID for the ST mapdownloader
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
# % key: k
# % label: Keep downloaded data in the download directory
# %end

# %rules
# % requires_all: -k,download_dir
# %end


import atexit
import sys
import tempfile
import shutil
import pathlib
import contextlib

import grass.script as grass
from grass.pygrass.utils import get_lib_path

from grass_gis_helpers.cleanup import general_cleanup, cleaning_tmp_location
from grass_gis_helpers.general import test_memory
from grass_gis_helpers.location import switch_back_original_location
from grass_gis_helpers.mapset import switch_to_new_mapset

# import module library
path = get_lib_path(modname="r.dop.import")
if path is None:
    grass.fatal("Unable to find the dop library directory.")
sys.path.append(path)
try:
    from r_dop_import_lib import (
        rescale_to_1_255,
        import_and_reproject,
        download_dop_st,
    )
except Exception as imp_err:
    grass.fatal(f"r.dop.import library could not be imported: {imp_err}")

rm_rast = []
rm_group = []
rm_files = []

gisdbase = None
TMP_LOC = None
TMP_GISRC = None
# pylint: disable=C0103
original_nprocs = None
tmp_download_dir = None
keep_data = False


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
    # Downloaded/extracted ST files: only remove own files, never the
    # whole download_dir, since it's shared with parallel workers
    if not keep_data:
        for f in rm_files:
            with contextlib.suppress(FileNotFoundError):
                pathlib.Path(f).unlink()
            # try:
            #     pathlib.Path(f).unlink()
            # except FileNotFoundError:
            #     pass
    # Remove auto-created temp download dir unless user asked to keep it
    if tmp_download_dir and not keep_data:
        shutil.rmtree(tmp_download_dir, ignore_errors=True)
    # Reset nprocs
    if original_nprocs:
        grass.run_command("g.gisenv", set=f"NPROCS={original_nprocs}")
    else:
        grass.run_command("g.gisenv", unset="NPROCS")


def main():
    """Main function of r.dop.import.worker.st"""
    # pylint: disable=C0301
    global gisdbase, TMP_LOC, TMP_GISRC, original_nprocs, tmp_download_dir, keep_data

    # parser options
    tile_key = options["tile_key"]
    item_id = options["item_id"]
    raster_name = options["raster_name"]
    resolution_to_import = None
    if options["resolution_to_import"]:
        resolution_to_import = float(options["resolution_to_import"])
    orig_region = options["orig_region"]
    new_mapset = options["new_mapset"]
    download_dir = options["download_dir"]
    keep_data = flags["k"]

    # A real download dir is required to store/extract the zip, even without
    # -k; create a temp on if the userr didn't give one
    if not download_dir:
        tmp_download_dir = tempfile.mkdtemp(prefix="rdop_import_st_")
        download_dir = tmp_download_dir

    # set nprocs to 1, write original value in variable
    gisenv = grass.gisenv()
    if "NPROCS" in gisenv:
        original_nprocs = int(gisenv["NPROCS"])
    grass.run_command("g.gisenv", set="NPROCS=1")

    # set memory to input if possible
    options["memory"] = test_memory(options["memory"])

    # switch to new mapset for parallel processing
    gisrc, newgisrc, old_mapset = switch_to_new_mapset(new_mapset)

    # set region
    grass.run_command("g.region", region=f"{orig_region}@{old_mapset}")
    aoi_map = f"{options['aoi']}@{old_mapset}" if options["aoi"] else None

    # import DOP tile with original resolution
    grass.message(
        _(f"Started DOP import for key: {tile_key}, item_id: {item_id}"),
    )

    # download and extract the tile locally (mandatory, not just for -k,
    # since the mapdownloader link needs a session and can't be read
    # directly by r.import/GDAL)
    local_tif, downloaded_files = download_dop_st(item_id, download_dir)
    rm_files.extend(downloaded_files)

    # import and reproject DOP tiles based on tileindex
    gisdbase, TMP_LOC, TMP_GISRC = import_and_reproject(
        local_tif,
        raster_name,
        resolution_to_import,
        "ST",
        aoi_map,
        download_dir,
        epsg=25833,
        keep_data=keep_data,
    )

    test_raster = f"{raster_name}.1"
    if not grass.find_file(test_raster, element="cell")["name"]:
        grass.warning(
            _(f"{item_id} could not be imported (no overlap)."),
        )
        switch_back_original_location(gisrc)
        grass.utils.try_remove(newgisrc)
        return

    sys.stderr.write(f"METADATA_DOP_URL:{item_id}\n")

    rm_group.append(raster_name)
    grass.message(_(f"Finishing raster import for {raster_name}..."))

    # rescale imported DOPs
    new_rm_rast = rescale_to_1_255("ST", raster_name)
    rm_rast.extend(new_rm_rast)

    # switch back to original location
    switch_back_original_location(gisrc)
    grass.utils.try_remove(newgisrc)
    grass.message(
        _(f"DOP import for key: {tile_key}, item_id: {item_id} done!"),
    )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
