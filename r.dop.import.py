#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Anika Weinmann, Julia Haas
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              (currently only for NRW)
# COPYRIGHT:   (C) 2022-2024 by mundialis GmbH & Co. KG and the GRASS
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
#############################################################################

# %Module
# % description: Downloads and imports Digital Othophotos (DOPs) (currently only for NRW)
# % keyword: imagery
# % keyword: download
# % keyword: DOP
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % required: yes
# %end

# %option G_OPT_V_INPUT
# % key: aoi_map
# % required: no
# % description: Vector map to restrict DOPs import to
# %end

# %option
# % key: filepath
# % required: no
# % description: Text file containing federal state to load DOPs for
# %end

# %option
# % key: federal_state
# % multiple: yes
# % required: no
# % description: Federal state to load DOPs for (no alternative to aoi_map; parameter is required/used for getting download-URL only)
# % options: Brandenburg,Berlin,Baden-Württemberg,Bayern,Bremen,Hessen,Hamburg,Mecklenburg-Vorpommern,Niedersachsen,Nordrhein-Westfalen,Rheinland-Pfalz,Schleswig-Holstein,Saarland,Sachsen,Sachsen-Anhalt,Thüringen
# %end

# %option G_OPT_M_DIR
# % key: local_data_dir
# % required: no
# % description: Directory with raster map of DOPs to import (e.g. VRT)
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

# %flag
# % key: r
# % description: use native DOP resolution
# %end

# %rules
# % required: federal_state, filepath
# % excludes: filepath, federal_state
# %end

import atexit
import glob
import os
import wget
import multiprocessing as mp
import grass.script as grass
from grass.pygrass.modules import Module, ParallelModuleQueue
import sys
from grass_gis_helpers.cleanup import general_cleanup
from grass_gis_helpers.general import communicate_grass_command


sys.path.insert(
    1,
    os.path.join(
        os.path.dirname(sys.path[0]),
        "etc",
        "r.dop.import",
    ),
)
from download_urls import URLS
from download_urls import WMS_HE
from download_urls import WMS_TH
from federal_state_info import FS_ABBREVIATION

tmp_dir = None
resolution_to_import = None
rm_vectors = []
rm_rasters = []
rm_groups = []
orig_region = None
mapset_names = []


def cleanup():
    rm_dirs = []
    if tmp_dir:
        rm_dirs.append(tmp_dir)
    general_cleanup(
        rm_vectors=rm_vectors,
        rm_rasters=rm_rasters,
        rm_groups_wo_rasters=rm_groups,
        orig_region=orig_region,
        rm_dirs=rm_dirs,
    )
    # remove mapsets
    for rm_mapset in mapset_names:
        gisenv = grass.gisenv()
        mapset_path = os.path.join(
            gisenv["GISDBASE"], gisenv["LOCATION_NAME"], rm_mapset
        )
        grass.try_rmdir(mapset_path)


def setup_parallel_processing(nprocs):
    if nprocs == -2:
        nprocs = mp.cpu_count() - 1 if mp.cpu_count() > 1 else 1
    else:
        # Test nprocs settings
        nprocs_real = mp.cpu_count()
        if nprocs > nprocs_real:
            grass.warning(
                "Using %d parallel processes but only %d CPUs available."
                % (nprocs, nprocs_real)
            )

    # set some common environmental variables, like:
    os.environ.update(
        dict(
            GRASS_COMPRESSOR="LZ4",
            GRASS_MESSAGE_FORMAT="plain",
        )
    )
    return nprocs


def reset_region(region):
    """Function to set the region to the given region
    Args:
        region (str): the name of the saved region which should be set and
                      deleted
    """
    nulldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nulldev}
    if region:
        if grass.find_file(name=region, element="windows")["file"]:
            grass.run_command("g.region", region=region)
            grass.run_command("g.remove", type="region", name=region, **kwargs)


def create_grid(tile_size, grid_prefix, area):
    """Create a grid for parallelization
    Args:
        tile_size (float): the size for the tiles in map units
        grid_prefix (str): the prefix name for the output grid
        area (str): the name of area for which to create the grid tiles
    Return:
        tiles_list (list): list with the names of the created vector map tiles
    """
    if area is None or area == "":
        area = f"tmp_aoi_{grass.tempname(8)}"
        rm_vectors.append(area)
        grass.run_command("v.in.region", output=area)
    # set region to area
    region = grass.parse_command("g.region", flags="ug", vector=area)
    dist_ns = abs(float(region["n"]) - float(region["s"]))
    dist_ew = abs(float(region["w"]) - float(region["e"]))

    grass.message(_("Creating tiles..."))
    grid = f"tmp_grid_{grass.tempname(8)}"
    # check if region is smaller than tile size
    if dist_ns <= float(tile_size) and dist_ew <= float(tile_size):
        grass.run_command(
            "g.region", vector=area, res=resolution_to_import, flags="a"
        )
        grass.run_command("v.in.region", output=grid, quiet=True)
        grass.run_command(
            "v.db.addtable", map=grid, columns="cat int", quiet=True
        )
    else:
        # set region
        orig_region = f"grid_region_{grass.tempname(8)}"
        grass.run_command("g.region", save=orig_region, quiet=True)
        grass.run_command("g.region", vector=area, quiet=True)
        grass.run_command("g.region", res=tile_size, flags="a", quiet=True)

        # create grid
        grass.run_command(
            "v.mkgrid", map=grid, box=f"{tile_size},{tile_size}", quiet=True
        )
        # reset region
        reset_region(orig_region)
    grid_name = f"tmp_grid_area_{grass.tempname(8)}"
    grass.run_command(
        "v.select",
        ainput=grid,
        binput=area,
        output=grid_name,
        operator="overlap",
        quiet=True,
    )
    if grass.find_file(name=grid_name, element="vector")["file"] == "":
        grass.fatal(
            _(
                f"The set region is not overlapping with {area}. "
                f"Please define another region."
            )
        )

    # create list of tiles
    tiles_num_list = list(
        grass.parse_command(
            "v.db.select", map=grid_name, columns="cat", flags="c", quiet=True
        ).keys()
    )

    number_tiles = len(tiles_num_list)
    grass.message(_(f"Number of tiles is: {number_tiles}"))
    tiles_list = []
    for tile in tiles_num_list:
        tile_area = f"{grid_prefix}_{tile}"
        grass.run_command(
            "v.extract",
            input=grid_name,
            where=f"cat == {tile}",
            output=tile_area,
            quiet=True,
        )
        tiles_list.append(tile_area)
        rm_vectors.append(tile_area)

    # cleanup
    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmv in [grid, grid_name]:
        if grass.find_file(name=rmv, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmv, **kwargs)

    return tiles_list


def get_tindex(tileindex):
    """Download and import tindex
    Args:
        tileindex ... URL to tile index
    Returns:
        vm_import ... Name of the tile index vector map
    """
    # download tindex
    zipname = os.path.basename(tileindex)
    tmp_dir = grass.tempdir()
    download_path = os.path.join(tmp_dir, zipname)
    wget.download(tileindex, download_path, bar=None)

    # unzip tindex
    unzipped_name = zipname.replace(".gz", "")
    unzipped_path = os.path.join(tmp_dir, unzipped_name)
    os.system(f"gunzip {download_path}")

    # import vector map containing URL for each tile
    vm_import = f"vm_import_{grass.tempname(8)}"
    grass.run_command(
        "v.import",
        input=unzipped_path,
        output=vm_import,
        extent="region",
        overwrite=True,
        quiet=True,
    )
    rm_vectors.append(vm_import)
    return vm_import


def download_and_clip_tindex(federal_state, aoi_map=None):
    """Download and clip tindex
    Args:
        aoi_map (str): name of AOI vector map
    Returns:
        (list): list with urls of tiles
    """
    tileindex = URLS[federal_state]
    if tileindex is None:
        grass.warning(_(f"{federal_state} is not yet implemented."))
        return []
    else:
        # if aoi_map given: set region to aoi_map extent
        if aoi_map:
            grass.run_command(
                "g.region",
                vector=aoi_map,
                res=resolution_to_import,
                flags="a",
                quiet=True,
            )

        # download and unzip tile index
        vm_import = get_tindex(tileindex)
        if not grass.find_file(name=vm_import, element="vector")["file"]:
            grass.fatal(
                _(
                    "No tile found in region. Please check if region or aoi "
                    "overlap with federal state"
                )
            )
        if not aoi_map:
            aoi_map = f"aoi_map_{grass.tempname(8)}"
            rm_vectors.append(aoi_map)
            grass.run_command("v.in.region", output=aoi_map)
        # check which tiles are needed for selected AOI
        vm_clip = f"vm_clip_{grass.tempname(8)}"
        grass.run_command(
            "v.clip",
            input=vm_import,
            clip=aoi_map,
            output=vm_clip,
            flags="d",
            overwrite=True,
            quiet=True,
        )
        rm_vectors.append(vm_clip)
        # else:
        #     # if no aoi given, use complete (current set) region:
        #     vm_clip = vm_import

        # import tiles and rename them according to their band
        # and write them in a list
        return grass.vector_db_select(vm_clip, columns="location")[
            "values"
        ].items()


def get_tiles(federal_state, aoi_map=None):
    """Get or create tileindex for federal state
    Args:
        federal_state (str): A string with a federal state
        aoi_map (str): Name of AOI vector map
    Returns:
        tileindex: None if no tileindex exists and one was created
        tiles_list (list): list of tiles (names of the created vector tiles or
                            URLs to download DOPs)
    """
    if federal_state in URLS:
        grass.message(f"Processing {federal_state}...")
        if federal_state in ["Hessen", "Thüringen"]:
            # create grid for wms import
            if federal_state == "Hessen":
                tiles_list = create_grid(1000, "HE_DOP", aoi_map)
            elif federal_state == "Thüringen":
                tiles_list = create_grid(1000, "TH_DOP", aoi_map)
            # no tileindex
            tileindex = None

        else:
            tileindex = URLS[federal_state]
            if tileindex is None:
                grass.warning(_(f"{federal_state} is not yet implemented."))
                tiles_list = []
            else:
                grass.message(_("Import tindex ..."))
                tiles_list = download_and_clip_tindex(federal_state, aoi_map)
    else:
        if options["filepath"]:
            grass.fatal(
                _(
                    "Non valid name of federal state,"
                    " in 'filepath'-option given"
                )
            )
        elif options["federal_state"]:
            grass.fatal(
                _(
                    "Non valid name of federal state,"
                    " in 'federal_states'-option given"
                )
            )
    return tileindex, tiles_list


def adjust_raster_resolution(raster_name, output, res):
    """Resample or inpolate raster to given resolution"""
    res_rast = float(
        grass.parse_command("r.info", map=raster_name, flags="g")["nsres"]
    )
    if res_rast > res:
        grass.run_command(
            "r.resamp.interp",
            input=raster_name,
            output=output,
            overwrite=True,
            quiet=True,
        )
    elif res_rast < res:
        grass.run_command(
            "r.resamp.stats",
            input=raster_name,
            output=output,
            method="median",
            quiet=True,
            overwrite=True,
        )
    else:
        rename_raster(raster_name, output)


def adjust_resolution_for_rasters(raster_base_name):
    """Resample or inpolate raster for all bands"""
    res = resolution_to_import
    # set region to imported raster
    grass.run_command("g.region", raster=f"{raster_base_name}.1", quiet=True)
    grass.run_command("g.region", res=res, flags="a", quiet=True)
    for band in range(1, 5):
        adjust_raster_resolution(
            f"{raster_base_name}.{band}", f"{raster_base_name}.{band}", res
        )


def create_vrt(b_list, out):
    # copy raster maps to current mapset
    for rast in b_list:
        if "@" in rast:
            rast_wo_mapsetname = rast.split("@")[0]
            grass.run_command(
                "g.copy",
                raster=f"{rast},{rast_wo_mapsetname}",
            )
    b_list = [val.split("@")[0] for val in b_list]
    # buildvrt if required + renaming to output name
    if len(b_list) > 1:
        grass.run_command("g.region", raster=b_list)
        grass.run_command(
            "r.buildvrt", input=b_list, output=out, quiet=True, overwrite=True
        )
    else:
        grass.run_command(
            "g.rename", raster=f"{b_list[0]},{out}", quiet=True, overwrite=True
        )


def download_and_import_dops(
    aoi,
    fs,
    federal_state,
    resolution_to_import,
    rm_rasters,
    orig_region,
    mapset_names,
    get_tiles,
    all_raster,
    nprocs,
):
    tileindex, tiles_list = get_tiles(federal_state, aoi)
    number_tiles = len(tiles_list)
    # set number of parallel processes to number of tiles
    # if multiple federal states given,
    # they are not calculated in parallel so far
    if number_tiles < nprocs:
        nprocs = number_tiles
    queue = ParallelModuleQueue(nprocs=nprocs)
    try:
        grass.message(_(f"Importing {len(tiles_list)} DOPs in parallel..."))
        for tile_el in tiles_list:
            if tileindex:
                key = tile_el[0]
                new_mapset = f"tmp_mapset_rdop_import_tile_{key}_{os.getpid()}"
                mapset_names.append(new_mapset)
                b_name = os.path.basename(tile_el[1][0])
                raster_name = (
                    f"{b_name.split('.')[0].replace('-', '_')}"
                    f"_{os.getpid()}"
                )
                for key_rast in all_raster:
                    all_raster[key_rast].append(
                        f"{fs}_{raster_name}_{key_rast}@{new_mapset}"
                    )
                param = {
                    "flags": "t",
                    "tile_key": key,
                    "tile_url": tile_el[1][0],
                    "federal_state": fs,
                    "raster_name": raster_name,
                    "orig_region": orig_region,
                    "memory": 1000,
                    "new_mapset": new_mapset,
                }
            else:
                key = tile_el
                new_mapset = f"tmp_mapset_rdop_import_tile_{key}_{os.getpid()}"
                mapset_names.append(new_mapset)
                raster_name = tile_el
                for key_rast in all_raster:
                    all_raster[key_rast].append(
                        f"{fs}_{raster_name}_{key_rast}@{new_mapset}"
                    )
                if fs == "HE":
                    wms = WMS_HE
                elif fs == "TH":
                    wms = WMS_TH
                param = {
                    "flags": "",
                    "tile_key": key,
                    "tile_url": wms,
                    "federal_state": fs,
                    "raster_name": raster_name,
                    "orig_region": orig_region,
                    "memory": 1000,
                    "new_mapset": new_mapset,
                }
            grass.message(f"raster_name: {raster_name}")
            if aoi:
                param["aoi_map"] = aoi

            if flags["r"]:
                param["flags"] += "r"
            else:
                param["resolution_to_import"] = resolution_to_import
            # add downloaded raster bands to rm_rast
            rm_red = f"{raster_name}_red"
            rm_green = f"{raster_name}_green"
            rm_blue = f"{raster_name}_blue"
            rm_rasters.append(rm_red)
            rm_rasters.append(rm_green)
            rm_rasters.append(rm_blue)
            # grass.run_command(
            r_dop_import_worker = Module(
                "r.dop.import.worker",
                **param,
                run_=False,
            )
            # catch all GRASS outputs to stdout and stderr
            r_dop_import_worker.stdout_ = grass.PIPE
            r_dop_import_worker.stderr_ = grass.PIPE
            queue.put(r_dop_import_worker)
        queue.wait()
    except Exception:
        for proc_num in range(queue.get_num_run_procs()):
            proc = queue.get(proc_num)
            if proc.returncode != 0:
                # save all stderr to a variable and pass it to a GRASS
                # exception
                errmsg = proc.outputs["stderr"].value.strip()
                grass.fatal(
                    _(f"\nERROR by processing <{proc.get_bash()}>: {errmsg}")
                )


def rename_raster(band_name_old, band_name_new):
    """Rename raster map"""
    grass.run_command(
        "g.rename",
        raster=f"{band_name_old},{band_name_new}",
        quiet=True,
        overwrite=True,
    )


def import_local_data(aoi, local_data_dir, fs, all_raster):
    """Import of data from local file path
    Args:
        aoi (str): name of vector map defining AOI
        local_data_dir (str): path to local data
        fs (str): federal state abbreviation
        all_raster (dict): dictionary with bands and output raster list
    Returns:
        imported_local_data (bool): True if local data imported, otherwise False
    """
    grass.message(_("Importing local data..."))
    imported_local_data = False
    # get files (VRT if available otherwise tif)
    dop_files = glob.glob(
        os.path.join(local_data_dir, fs, "**", "*.vrt"),
        recursive=True,
    )
    if not dop_files:
        dop_files = glob.glob(
            os.path.join(local_data_dir, fs, "**", "*.tif"),
            recursive=True,
        )

    # import data for AOI
    # TODO parallize local data import for multi tifs
    for i, dop_file in enumerate(dop_files):
        cur_reg = grass.region()
        ns_res = cur_reg["nsres"]
        ew_res = cur_reg["ewres"]
        if aoi and aoi != "":
            grass.run_command(
                "g.region",
                vector=aoi,
                nsres=ns_res,
                ewres=ew_res,
                flags="a",
                quiet=True,
            )
        name = f"dop_{i}"
        kwargs = {
            "input": dop_file,
            "output": name,
            "extent": "region",
            "quiet": True,
            "overwrite": True,
        }
        # resolution settings: -r native resolution; otherwise from region
        if not flags["r"]:
            kwargs["resolution"] = "region"
        r_import = communicate_grass_command("r.import", **kwargs)
        err_m1 = "Input raster does not overlap current computational region."
        err_m2 = "already exists and will be overwritten"
        if err_m1 in r_import[1].decode():
            continue
        elif err_m2 in r_import[1].decode():
            pass
        elif r_import[1].decode() != "":
            grass.fatal(_(r_import[1].decode()))

        band_dict = {1: "red", 2: "green", 3: "blue", 4: "nir"}
        for band_num, band in band_dict.items():
            band_name_old = f"{name}.{band_num}"
            band_name_new = f"{name}.{band}"
            rm_rasters.append(band_name_old)
            # check resolution and resample if needed
            if not flags["r"]:
                adjust_raster_resolution(band_name_old, band_name_new, ns_res)
            else:
                rename_raster(band_name_old, band_name_new)
            all_raster[band].append(band_name_new)

    # check if all bands have at least one input raster
    if all([len(val) > 0 for band, val in all_raster.items()]):
        imported_local_data = True

    if not imported_local_data and fs in ["BW"]:
        grass.fatal(_("Local data does not overlap with AOI."))
    elif not imported_local_data:
        grass.message(
            _(
                "Local data does not overlap with AOI. Data will be downloaded"
                " from Open Data portal."
            )
        )
    return imported_local_data


def main():
    global tmp_dir, orig_region, rm_vectors, rm_rasters, rm_groups
    global resolution_to_import, mapset_names

    # a vector map, consisting of the DOP-tiles,
    # while each DOP-tile contains its corresponding
    # download link.
    # By overlaying the aoi_map with this vector map
    # the download links of the DOPs within the area
    # of interest are selected.
    # However, when the aoi_map lies exactly at the border
    # of such DOP tile, a problem due to floating point
    # precision limit can occur. This results in a larger
    # overlay-area of aoi and vector map, than it is
    # actually the case.
    # Two approaches for solution:
    # 1. modify the gpkg-file: delete the decimals
    #    (this avoids the floating proint precision problem)
    # 2. do not use the gpkg-file; instead generate the
    #    DOP download links via the UTM-coordinates

    # parser options
    aoi = options["aoi_map"]
    local_data_dir = options["local_data_dir"]
    nprocs = int(options["nprocs"])
    nprocs = setup_parallel_processing(nprocs)

    # check if required addons installed:
    addon = "r.dop.import.worker"
    if not grass.find_program(addon, "--help"):
        msg = (
            f"The '{addon}' module was not found, install  it first:\n"
            f"g.extension {addon}"
        )
        grass.fatal(_(msg))

    # create list for each raster band for building entire raster
    # of all given federal states
    all_raster = {
        "red": [],
        "green": [],
        "blue": [],
        "nir": [],
    }

    # check if overwrite is not activated
    if os.getenv("GRASS_OVERWRITE", "0") != "1":
        # check if output raster maps already exists
        for band in all_raster:
            out = f"{options['output']}_{band}"
            if grass.find_file(name=out, element="raster")["file"]:
                grass.fatal(
                    _(
                        f"output map <{out}> exists. To overwrite, use the "
                        "--overwrite flag"
                    )
                )

    # read federal state(s) from input options
    if options["filepath"]:
        with open(f'{options["filepath"]}') as f:
            federal_states = f.read().strip().split(",")
    else:
        federal_states = options["federal_state"].split(",")

    # get list of local input folders for federal states
    local_fs_list = []
    if local_data_dir and local_data_dir != "":
        local_fs_list = os.listdir(local_data_dir)

    # Berlin and Brandenburg are using the same wms of BB
    if "Berlin" in federal_states and "Brandenburg" in federal_states:
        BE_idx = federal_states.index("Berlin")
        del federal_states[BE_idx]

    # save current region for setting back later in cleanup
    orig_region = f"orig_region_{grass.tempname(8)}"
    grass.run_command("g.region", save=orig_region, overwrite=True, quiet=True)
    reg = grass.region()
    if flags["r"]:
        # TODO: event auf None
        resolution_to_import = 0.2
    else:
        if reg["nsres"] == reg["ewres"]:
            resolution_to_import = float(reg["nsres"])
        else:
            grass.fatal("N/S resolution is not the same as E/W resolution!")

    # loop through federal states and get respective DOPs
    for federal_state in federal_states:
        if federal_state not in FS_ABBREVIATION:
            grass.fatal(_(f"Non valid name of federal state: {federal_state}"))
        fs = FS_ABBREVIATION[federal_state]

        # check if local data for federal state given
        imported_local_data = False
        if fs in local_fs_list:
            # TODO import_local_data with -r flag
            imported_local_data = import_local_data(
                aoi, local_data_dir, fs, all_raster
            )
        elif fs in ["BW"]:
            grass.fatal(
                _(f"No local data for {fs} available. Is the path correct?")
            )

        # import data when local import was not used
        if not imported_local_data:
            # check if federal state is supported
            if fs in ["BE", "BB", "HE", "NW", "SN", "TH"]:
                # import data
                download_and_import_dops(
                    aoi,
                    fs,
                    federal_state,
                    resolution_to_import,
                    rm_rasters,
                    orig_region,
                    mapset_names,
                    get_tiles,
                    all_raster,
                    nprocs,
                )
            else:
                grass.warning(_(f"Support for {fs} is not yet implemented."))

    raster_out = []
    for band, b_list in all_raster.items():
        out = f"{options['output']}_{band}"
        create_vrt(b_list, out)
        raster_out.append(out)

    grass.message(_(f"Generated following raster maps: {raster_out}"))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
