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
import pathlib
import contextlib
import os
import zipfile
import requests

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

EPSG = 25833


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
    # Reset nprocs
    if original_nprocs:
        grass.run_command("g.gisenv", set=f"NPROCS={original_nprocs}")
    else:
        grass.run_command("g.gisenv", unset="NPROCS")


def download_dop_st(item_id, download_dir):
    """Download and extract a single ST DOP tile via the two-step
    prepare/download mechanism of the LVermGeo mapdownloader.

    Args:
        item_id (str): Numeric tile ID from the ST tindex
        download_dir (str): Local directory to download/extract into

    Returns:
        tuple: Path to the extracted .tif, list of all created file paths for
               later cleanup
    """
    # Server needs a browser-like User-Agent; default/non-browser UA strings
    # will be blocked with HTTP 503
    pathlib.Path(download_dir).mkdir(exist_ok=True, parents=True)
    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 (X11; Linux x86_64; rv:155.0) "
                "Gecko/20100101 Firefox/155.0"
            ),
        },
    )
    base_url = "https://www.lvermgeo.sachsen-anhalt.de/"
    session.get(base_url)

    prepare_url = (
        f"{base_url}de/mod/4,1962,501/ajax/1/prepare/"
        f"?items={item_id}&format=zip"
    )
    resp = session.get(
        prepare_url,
        headers={"X-Requested-With": "XMLHttpRequest"},
    )
    resp.raise_for_status()
    download_url = resp.text.strip()
    if not download_url.startswith("http"):
        grass.fatal(
            _(
                f"Unexpected prepare response for item {item_id}: "
                f"{download_url}",
            ),
        )

    dl_resp = session.get(download_url)
    dl_resp.raise_for_status()

    zip_path = os.path.join(download_dir, f"dop20_st_{item_id}.zip")
    pathlib.Path(zip_path).write_bytes(dl_resp.content)

    created_files = [zip_path]
    with zipfile.ZipFile(zip_path) as zf:
        tif_names = [n for n in zf.namelist() if n.lower().endswith(".tif")]
        if not tif_names:
            grass.fatal(
                _(f"No .tif found in ZIP for item {item_id}"),
            )
        zf.extractall(download_dir)
        created_files.extend(
            os.path.join(download_dir, n) for n in zf.namelist()
        )

    return os.path.join(download_dir, tif_names[0]), created_files


def main():
    """Main function of r.dop.import.worker.st"""
    global gisdbase, TMP_LOC, TMP_GISRC, original_nprocs, keep_data

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
        EPSG,
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
