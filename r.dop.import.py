#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer, Lina Krisztian, Anika Weinmann
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
#              (currently only for NRW)
# COPYRIGHT:   (C) 2022-2023 by mundialis GmbH & Co. KG and the GRASS
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
# % description: Federal state to load DOPs for
# % options: Brandenburg,Berlin,Baden-Württemberg,Bayern,Bremen,Hessen,Hamburg,Mecklenburg-Vorpommern,Niedersachsen,Nordrhein-Westfalen,Rheinland-Pfalz,Schleswig-Holstein,Saarland,Sachsen,Sachsen-Anhalt
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
import wget
import os
import shutil
import grass.script as grass
import sys

sys.path.insert(
    1,
    os.path.join(
        os.path.dirname(sys.path[0]),
        "etc",
        "r.dop.import",
    ),
)
from download_urls import URLS
from federal_states import FS

tmp_dir = None
resolution_to_import = None
rm_vec = []
rm_rast = []
rm_group = []
temp_region = None
TMPLOC = None
SRCGISRC = None
TGTGISRC = None
GISDBASE = None


def cleanup():
    if tmp_dir:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmvec in rm_vec:
        if grass.find_file(name=rmvec, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmvec, **kwargs)
    for rmrast in rm_rast:
        if grass.find_file(name=rmrast, element="raster")["file"]:
            grass.run_command("g.remove", type="raster", name=rmrast, **kwargs)
    for rmgroup in rm_group:
        if grass.find_file(name=rmgroup, element="group")["file"]:
            grass.run_command("g.remove", type="group", name=rmgroup, **kwargs)
    if temp_region:
        # set region back and delete saved region:
        grass.run_command("g.region", region=temp_region)
        grass.run_command(
            "g.remove", type="region", name=temp_region, flags="f", quiet=True
        )
    # remove temp location
    if TMPLOC:
        grass.try_rmdir(os.path.join(GISDBASE, TMPLOC))
    if SRCGISRC:
        grass.try_remove(SRCGISRC)


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


def get_actual_location():
    global TGTGISRC, GISDBASE
    # get actual location, mapset, ...
    grassenv = grass.gisenv()
    tgtloc = grassenv["LOCATION_NAME"]
    tgtmapset = grassenv["MAPSET"]
    GISDBASE = grassenv["GISDBASE"]
    TGTGISRC = os.environ["GISRC"]
    return tgtloc, tgtmapset


def createTMPlocation(epsg=4326):
    global TMPLOC, SRCGISRC
    SRCGISRC = grass.tempfile()
    TMPLOC = f"temp_import_location_{grass.tempname(8)}"
    f = open(SRCGISRC, "w")
    f.write("MAPSET: PERMANENT\n")
    f.write("GISDBASE: %s\n" % GISDBASE)
    f.write("LOCATION_NAME: %s\n" % TMPLOC)
    f.write("GUI: text\n")
    f.close()

    proj_test = grass.parse_command("g.proj", flags="g")
    if "epsg" in proj_test:
        epsg_arg = {"epsg": epsg}
    else:
        epsg_arg = {"srid": "EPSG:{}".format(epsg)}
    # create temp location from input without import
    grass.verbose(_("Creating temporary location with EPSG:%d...") % epsg)
    grass.run_command(
        "g.proj", flags="c", location=TMPLOC, quiet=True, **epsg_arg
    )

    # switch to temp location
    os.environ["GISRC"] = str(SRCGISRC)
    proj = grass.parse_command("g.proj", flags="g")
    if "epsg" in proj:
        new_epsg = proj["epsg"]
    else:
        new_epsg = proj["srid"].split("EPSG:")[1]
    if new_epsg != str(epsg):
        grass.fatal(_("Creation of temporary location failed!"))


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
        rm_vec.append(area)
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
        rm_vec.append(tile_area)

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
    rm_vec.append(vm_import)
    return vm_import


def import_and_reproject(url, raster_name, aoi_map=None, epsg=None):
    """Import DOPs and reprojects them if needed. This is needed for the DOPs
    of Brandenburg and Berlin because GDAL (at least smaller 3.6.3) does not
    support the coordinate reference system in the data.

    Args:
        url (str): The URL of the DOP to import
        raster_name (str): The prefix name for the output rasters
        aoi_map (str): Name of AOI vector map
        epsg (int): EPSG code which has to be set if the reproduction should be
                    done manually and not by r.import
    """
    import_flags = ""
    if epsg is not None:
        if not aoi_map:
            aoi_map = f"region_aoi_{grass.tempname(8)}"
            grass.run_command("v.in.region", output=aoi_map, quiet=True)
        # get actual location, mapset, ...
        tgtloc, tgtmapset = get_actual_location()
        # create temporary location with epsg:4326
        createTMPlocation(epsg)
        import_flags = "o"
        # reproject aoi
        if aoi_map:
            grass.run_command(
                "v.proj",
                location=tgtloc,
                mapset=tgtmapset,
                input=aoi_map,
                output=aoi_map,
                quiet=True,
            )
    if aoi_map:
        grass.run_command(
            "g.region", vector=aoi_map, res=resolution_to_import, flags="a"
        )

    # import data
    kwargs = {
        "input": url,
        "output": raster_name,
        "memory": 1000,
        "quiet": True,
        "flags": import_flags,
        "extent": "region",
    }
    grass.run_command("r.import", **kwargs)

    # reproject data
    res = float(
        grass.parse_command("r.info", flags="g", map=f"{raster_name}.1")[
            "nsres"
        ]
    )
    if epsg is not None:
        # switch location
        os.environ["GISRC"] = str(TGTGISRC)
        if aoi_map:
            grass.run_command("g.region", vector=aoi_map, res=res, flags="a")
        else:
            grass.run_command("g.region", res=res, flags="a")
        for i in range(1, 5):
            name = f"{raster_name}.{i}"
            grass.run_command(
                "r.proj",
                location=TMPLOC,
                mapset="PERMANENT",
                input=name,
                output=name,
                resolution=res,
                flags="n",
                quiet=True,
                memory=1000,
            )


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
        if aoi_map:
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
            rm_vec.append(vm_clip)
        else:
            # if no aoi given, use complete (current set) region:
            vm_clip = vm_import

        # import tiles and rename them according to their band
        # and write them in a list
        return grass.vector_db_select(vm_clip, columns="location")[
            "values"
        ].items()


def get_tiles(federal_state, aoi_map=None):
    """Get or create tileindex for federal state
    Args:
        federal_state (str): A string with a federal state
        aoi_mao (str): Name of AOI vector map
    Returns:
        tileindex: None if no tileindex exists and one was created
        tiles_list (list): list of tiles (names of the created vector tiles or
                            URLs to download DOPs)
    """
    if federal_state in URLS:
        grass.message(f"Processing {federal_state}...")
        if federal_state == "Hessen":
            # create grid for wms import
            tiles_list = create_grid(1000, "HE_DOP", aoi_map)
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


def adjust_resolution(raster_name):
    """Resample or inpolate raster"""
    res = resolution_to_import
    # set region to imported tile
    grass.run_command("g.region", raster=f"{raster_name}.1", quiet=True)
    grass.run_command("g.region", res=res, flags="a", quiet=True)

    res_rast = float(
        grass.parse_command("r.info", map=f"{raster_name}.1", flags="g")[
            "nsres"
        ]
    )
    res_region = float(grass.region()["nsres"])
    if res_rast > res_region:
        for band in range(1, 5):
            grass.run_command(
                "r.resamp.interp",
                input=f"{raster_name}.{band}",
                output=f"{raster_name}.{band}",
                overwrite=True,
            )
    elif res_rast < res_region:
        for band in range(1, 5):
            grass.run_command(
                "r.resamp.stats",
                input=f"{raster_name}.{band}",
                output=f"{raster_name}.{band}",
                overwrite=True,
            )


def rescale_to_1_256(prefix, raster_name, all_raster, extension="num"):
    """Rescale raster from 0 to 255 to 1 to 256
    Args:
        raster_name (str): Name of raster prefix
        all_raster (dict): Dict with all rasters which should to be merged
    """
    if extension == "num":
        band_dict = {"red": 1, "green": 2, "blue": 3, "nir": 4}
    else:
        band_dict = {
            "red": "red",
            "green": "green",
            "blue": "blue",
            "nir": "nir",
        }
    for name, num in band_dict.items():
        rastername = f"{prefix}_{raster_name}_{name}"
        grass.run_command(
            "r.mapcalc",
            expression=f"{rastername} = {raster_name}.{num} + 1",
            quiet=True,
            region="intersect",
        )
        rm_rast.append(f"{raster_name}.{num}")
        all_raster[name].append(rastername)


def import_dop_from_he_wms(tile):
    """Import DOP from HE WMS"""
    from download_urls import WMS_HE

    grass.run_command("g.region", vector=tile)
    for name in ["cir", "rgb"]:
        if name == "cir":
            out_tmp = f"{name}_{tile}_tmp"
            bands = ["red"]
        else:
            out_tmp = f"{tile}_tmp"
            bands = ["red", "green", "blue"]
        rm_group.append(out_tmp)
        for band in ["red", "green", "blue"]:
            rm_rast.append(f"{out_tmp}.{band}")
        grass.run_command(
            "r.in.wms",
            url=WMS_HE,
            layer=f"he_dop20_{name}",
            output=out_tmp,
            format="tiff",
            flags="b",
            overwrite=True,
        )
        for band in bands:
            if name == "cir":
                oband = "nir"
            else:
                oband = band
            if not flags["r"]:
                grass.run_command(
                    "r.resample",
                    input=f"{out_tmp}.{band}",
                    output=f"{tile}.{oband}",
                    overwrite=True,
                    quiet=True,
                )
            else:
                grass.run_command(
                    "g.rename",
                    raster=f"{out_tmp}.{band},{tile}.{oband}",
                    overwrite=True,
                    quiet=True,
                )
    # drop rest of CIR DOP
    grass.run_command(
        "g.remove", type="raster", pattern="cir_*_tmp", flags="f"
    )


def create_vrt(b_list, out):
    if len(b_list) > 1:
        grass.run_command("g.region", raster=b_list)
        grass.run_command("r.buildvrt", input=b_list, output=out, quiet=True)
    else:
        grass.run_command("g.rename", raster=f"{b_list[0]},{out}", quiet=True)


def main():
    global tmp_dir, temp_region, rm_vec, rm_rast, rm_group, resolution_to_import

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
    aoi_map = options["aoi_map"]

    # read federal state(s) from input options
    if options["filepath"]:
        with open(f'{options["filepath"]}') as f:
            federal_states = f.read().strip().split(",")
    else:
        federal_states = options["federal_state"].split(",")

    # Berlin and Brandenburg are using the same wms of BB
    if "Berlin" in federal_states and "Brandenburg" in federal_states:
        BE_idx = federal_states.index("Berlin")
        del federal_states[BE_idx]

    # save current region for setting back later in cleanup
    temp_region = f"temp_region_{grass.tempname(8)}"
    grass.run_command("g.region", save=temp_region, overwrite=True, quiet=True)
    reg = grass.region()
    if flags["r"]:
        resolution_to_import = 0.2
    else:
        if reg["nsres"] == reg["ewres"]:
            resolution_to_import = float(reg["nsres"])
        else:
            grass.fatal("N/S resolution is not the same as E/W resolution!")

    # create list for each raster band for building entire raster
    # of all given federal states
    all_raster = {
        "red": [],
        "green": [],
        "blue": [],
        "nir": [],
    }

    # loop through federal states and get respective DOPs
    for federal_state in federal_states:
        tileindex, tiles_list = get_tiles(federal_state, aoi_map)

        if tileindex:
            num_imp = len(tiles_list)
            num = 0
            for item in tiles_list:
                URL_tile = item[1]
                num += 1
                b_name = os.path.basename(URL_tile[0])
                raster_name = (
                    f"{b_name.split('.')[0].replace('-', '_')}"
                    f"_{os.getpid()}"
                )
                grass.message(
                    f"Started raster import {num}/{num_imp} from URL: {b_name}"
                )
                # Process stuck, when memory is too large (100000)
                grass.run_command("g.region", region=temp_region)
                if federal_state in ["Brandenburg", "Berlin"]:
                    import_and_reproject(
                        URL_tile[0], raster_name, aoi_map, epsg=25833
                    )
                else:
                    import_and_reproject(URL_tile[0], raster_name, aoi_map)

                # adjust resolution if required
                if not flags["r"]:
                    adjust_resolution(raster_name)
                rm_group.append(raster_name)

                grass.message(f"Finishing raster import for {raster_name} ...")
                rescale_to_1_256(FS[federal_state], raster_name, all_raster)

        else:
            # no difference between with or without -r flag (HE WMS is res=0.2)
            for tile in tiles_list:
                import_dop_from_he_wms(tile)
                grass.message(f"Finishing raster import for {tile} ...")
                rescale_to_1_256(
                    FS[federal_state], tile, all_raster, extension="str"
                )

    raster_out = []
    for band, b_list in all_raster.items():
        out = f"{options['output']}_{band}"
        create_vrt(b_list, out)
        raster_out.append(out)

    grass.message(_(f"Generated following raster maps: {raster_out}"))

    if tileindex or aoi_map:
        grass.run_command("g.region", region=temp_region)


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
