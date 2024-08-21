#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r_dop_import_lib
# AUTHOR(S):   Johannes Halbauer
#
# PURPOSE:     Library for r.dop.import
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

import os
from time import sleep
import grass.script as grass

from grass_gis_helpers.general import set_nprocs
from grass_gis_helpers.location import (
    get_current_location,
    create_tmp_location,
)
from grass_gis_helpers.open_geodata_germany.download_data import (
    download_data_using_threadpool,
    extract_compressed_files,
)
from grass_gis_helpers.raster import adjust_raster_resolution, rename_raster

OPEN_DATA_AVAILABILITY = {
    "NO_OPEN_DATA": ["BW", "BY", "HB", "HH", "MV", "RP", "SL", "ST", "SH"],
    "NOT_YET_SUPPORTED": ["NI"],
    "SUPPORTED": ["BE", "BB", "NW", "SN", "TH", "HE"],
}

RETRIES = 30
WAITING_TIME = 10


def setup_parallel_processing(nprocs):
    """Get possible number of workers and modify environment variables
    Args:
        nprocs (int): Number of workers to use
    Returns:
        nprocs (int): Possible number of workers to use
    """
    nprocs = set_nprocs(nprocs)
    # set some common environmental variables, like:
    os.environ.update(
        {
            "GRASS_COMPRESSOR": "LZ4",
            "GRASS_MESSAGE_FORMAT": "plain",
        },
    )
    return nprocs


def rescale_to_1_256(prefix, raster_name, extension="num"):
    """Rescale raster from 0 to 255 to 1 to 256
    Args:
        prefix (str): Name of federal state
        raster_name (str): Name of raster prefix
    """
    rm_rast = []
    if extension == "num":
        band_dict = {
            "red": 1,
            "green": 2,
            "blue": 3,
            "nir": 4,
        }
    for name, num in band_dict.items():
        grass.run_command("g.region", raster=f"{raster_name}.{num}")
        rastername = f"{prefix}_{raster_name}_{name}"
        grass.run_command(
            "r.mapcalc",
            expression=f"{rastername} = {raster_name}.{num} + 1",
            quiet=True,
            region="intersect",
        )
        rm_rast.append(f"{raster_name}.{num}")

    return rm_rast


def create_grid_and_tiles_list(
    ns_res,
    ew_res,
    tile_size,
    grid,
    rm_vectors,
    aoi,
    id,
    fs,
):
    """Check if aoi is smaller than grid tile size and create grid if not.
    Also create a list containing tiles which overlap with the aoi.
    Args:
        ns_res (float): Vertical resolution
        ew_res (float): Horizontal resolution
        tile_size (int): Size of grid tiles to create
        grid (str): Name of grid to create
        rm_vectors (list): List of vectors to remove in cleanup
        aoi (str): Name of aoi
        id (str): id used in GRASS session
        fs (str): Abbreviation of federal state

    Returns:
        rm_vectors (list): Extended list of vectors to remove in cleanup
        nummber_tiles (str): Number of tiles overlapping with aoi
        tiles_list (list): List of tile names overlapping with aoi
    """
    # check if aoi is smaller than tile size
    if ns_res <= float(tile_size) and ew_res <= float(tile_size):
        grass.run_command("v.in.region", output=grid, quiet=True)
        rm_vectors.append(grid)
        grass.run_command(
            "v.db.addtable",
            map=grid,
            columns="cat int",
            quiet=True,
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
    grid_name = f"tmp_grid_area_{id}"

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
            "v.db.select",
            map=grid_name,
            columns="cat",
            flags="c",
            quiet=True,
        ).keys(),
    )
    number_tiles = len(tiles_num_list)

    grass.message(_(f"Number of tiles: {number_tiles}"))
    tiles_list = []
    for tile in tiles_num_list:
        tile_area = f"{fs}_DOP_{tile}"
        grass.run_command(
            "v.extract",
            input=grid_name,
            where=f"cat == {tile}",
            output=tile_area,
            quiet=True,
        )
        tiles_list.append(tile_area)
        rm_vectors.append(tile_area)

    return rm_vectors, number_tiles, tiles_list


def import_dop_from_wms(
    tile_key,
    rastername,
    tile_url,
    resolution_to_import,
    layer_list,
    cir_band,
    wms_layer,
    rm_group,
    rm_rast,
    native_res,
):
    """Import DOPs from WMS

    Args:
        tile_key (str): Key of current tile
        rastername (str): Name of resulting raster
        tile_url (str): WMS URL to get DOPs from
        resolution_to_import (float): Resolution to resample imported raster to
        layer_list (list): List of WMS layers to use
        cir_band (str): Name of CIR band layer in WMS
        wms_layer (str): wms_layer to import from WMS
        rm_group (list): List of elements to remove in cleanup
        rm_rast (list): List of raster maps to remove in cleanup
        native_res (bool): Keep native DOP resolution

    Returns:
        rm_group (list): Extended list of elements to remove in cleanup
        rm_rast (list): Extended list of raster maps to remove in cleanup
    """
    # set region and create variable names
    grass.run_command("g.region", vector=tile_key)
    tile_key = tile_key.split("@")[0]
    for name in layer_list:
        if name == cir_band:
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
                    layer=f"{wms_layer}{name}",
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
                if count > (RETRIES / 2):
                    grass.fatal(f"Download of {tile_url} not working.")
                sleep(10)

        # change band name to band number
        band_nums = {
            "red": 1,
            "green": 2,
            "blue": 3,
        }
        for band in bands:
            oband = "4" if name == cir_band else band_nums[band]

            # create old and new name for adjusting resolution/renaming
            old_name = f"{out_tmp}.{band}"
            new_name = f"{rastername}.{oband}"

            # adjust resolution/rename raster maps
            if not native_res:
                adjust_raster_resolution(
                    old_name,
                    new_name,
                    resolution_to_import,
                )
            else:
                rename_raster(old_name, new_name)
            rm_rast.append(old_name)

    # drop rest of CIR DOP
    grass.run_command(
        "g.remove",
        type="raster",
        pattern=f"{cir_band}_*_tmp*",
        flags="f",
    )


def keep_data_nw(url, download_dir):
    """Download and keep DOPs for NW from url using threadpool

    Args:
        url (string): Url to download data from
        download_dir (str): Path to directory to download data to

    Returns:
        (str): Path to download data
    """
    url = url.replace("/vsicurl/", "")
    basename = os.path.basename(url)
    download_data_using_threadpool([url], download_dir, None)

    return os.path.join(download_dir, basename)


def keep_data_bb_be(url, download_dir):
    """Download and keep DOPs for BB/BE from url using threadpool

    Args:
        url (string): Url to download data from
        download_dir (str): Path to directory to download data to

    Returns:
        (str): Path to download data
    """
    url = url.replace("/vsizip/vsicurl/", "")
    basename = os.path.basename(url)
    url = url.replace(basename, "")[:-1]
    download_data_using_threadpool([url], download_dir, None)
    extract_compressed_files([basename.replace(".tif", ".zip")], download_dir)
    return os.path.join(download_dir, basename)


def keep_data_sn(url, download_dir):
    """Download and keep DOPs for BB/BE from url using threadpool

    Args:
        url (string): Url to download data from
        download_dir (str): Path to directory to download data to

    Returns:
        (str): Path to download data
    """
    basename = os.path.basename(url)
    print(basename)
    url = os.path.dirname(url.replace("/vsizip/vsicurl/", ""))
    download_data_using_threadpool([url], download_dir, None)
    extract_compressed_files([os.path.basename(url)], download_dir)

    return os.path.join(download_dir, basename)


def import_and_reproject(
    url,
    raster_name,
    resolution_to_import,
    fs,
    aoi_map=None,
    download_dir=None,
    epsg=None,
    keep_data=None,
):
    """Import DOPs and reproject them if needed.

    Args:
        url (str): The URL of the DOP to import
        raster_name (str): The prefix name for the output rasters
        resolution_to_import (float): Resolution for region and raster import
        fs (str): Abbreviation of the current federal state
        aoi_map (str): Name of AOI vector map
        download_dir (str): Path to local directory to downlaod DOPs to
        epsg (int): EPSG code which has to be set if the reproduction should be
                    done manually and not by r.import
        keep_data (bool): Download raster data to local directory and keep it

    Returns:
        gisdbase (str): Path to GISDBASE
        tmp_loc (str): Name of temporary location
        tmp_gisrc (str): Path to GISRC file
    """
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

    # create temporary location with EPSG:25832
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

    # define import parameters
    # set memory manually to 1000
    # Process stuck, when memory is too large (100000)
    # GDAL_CACHEMAX will be interpreted as MB only
    kwargs = {
        "input": url,
        "output": raster_name,
        "memory": 1000,
        "quiet": True,
        "extent": "region",
    }
    if fs == "BB_BE":
        kwargs["flags"] = "o"

    # download and keep data to download dir if -k flag is set
    # and change input parameter in kwargs to local path
    if keep_data:
        # call download and keep data function for respective federal state
        function_name = f"keep_data_{fs.lower()}"
        if function_name in globals():
            func = globals()[function_name]
            kwargs["input"] = func(url, download_dir)
        else:
            grass.fatal(
                _(
                    f"Function to download and keep data for {fs} ",
                    "not found in lib.",
                ),
            )

    # import data
    import_sucess = False
    tries = 0
    while not import_sucess:
        tries += 1
        if tries > RETRIES:
            grass.fatal(
                _(
                    f"Importing {kwargs['input']} failed after {RETRIES} "
                    "retries.",
                ),
            )
        try:
            grass.run_command("r.import", **kwargs)
            import_sucess = True
        except Exception:
            sleep(WAITING_TIME)
    if not aoi_map:
        grass.run_command("g.region", raster=f"{raster_name}.1")

    # reproject data
    if resolution_to_import:
        res = resolution_to_import
    else:
        res = float(
            grass.parse_command("r.info", flags="g", map=f"{raster_name}.1")[
                "nsres"
            ],
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

    # return temp location parameters to remove it in cleanup
    return gisdbase, tmp_loc, tmp_gisrc
