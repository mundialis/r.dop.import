#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer & Lina Krisztian
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area (currently only for NRW)
# COPYRIGHT:   (C) 2022 by mundialis GmbH & Co. KG and the GRASS Development Team
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
rm_vec = []
rm_rast = []
rm_group = []
temp_region = None
temp_region_2 = None
pid = os.getpid()
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
    if temp_region_2:
        grass.run_command(
            "g.remove",
            type="region",
            name=temp_region_2,
            flags="f",
            quiet=True,
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
        grid_prefix (list): list with the names of the created vector map tiles
        number_tiles (int): Number of created tiles
    """
    # set region to area
    # TODO check if area is set everytime
    grass.run_command("g.region", vector=area, res=0.2, flags="a")
    # check if region is smaller than tile size
    region = grass.region()
    dist_ns = abs(region["n"] - region["s"])
    dist_ew = abs(region["w"] - region["e"])

    grass.message(_("Creating tiles..."))
    grid = f"tmp_grid_{os.getpid()}"
    if dist_ns <= float(tile_size) and dist_ew <= float(tile_size):
        grass.run_command("v.in.region", output=grid, quiet=True)
        grass.run_command(
            "v.db.addtable", map=grid, columns="cat int", quiet=True
        )
    else:
        # set region
        orig_region = f"grid_region_{os.getpid()}"
        grass.run_command("g.region", save=orig_region, quiet=True)
        grass.run_command("g.region", res=tile_size, flags="a", quiet=True)

        # create grid
        grass.run_command(
            "v.mkgrid", map=grid, box=f"{tile_size},{tile_size}", quiet=True
        )
        # reset region
        reset_region(orig_region)
    grid_name = f"tmp_grid_area_{os.getpid()}"
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

    # cleanup
    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmv in [grid, grid_name]:
        if grass.find_file(name=rmv, element="vector")["file"]:
            grass.run_command("g.remove", type="vector", name=rmv, **kwargs)

    return tiles_list  # , number_tiles


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
    vm_import = f"vm_import_{pid}"
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


def import_and_reproject(url, raster_name, epsg=None):
    """Import DOPs and reprojects them if needed. This is needed for the DOPs
    of Brandenburg because GDAL (at least smaller 3.6.3) does not support the
    coordinate reference system in the data.

    Args:
        url (str): The URL of the DOP to import
        raster_name (str): The prefix name for the output rasters
        epsg (int): EPSG code which has to be set if the reproduction should be
                    done manually and not by r.import
    """
    import_flags = ""
    if epsg is not None:
        # get actual location, mapset, ...
        tgtloc, tgtmapset = get_actual_location()
        # create temporary location with epsg:4326
        createTMPlocation(epsg)
        import_flags = "o"
    # import data
    grass.run_command(
        "r.import",
        input=url,
        output=raster_name,
        memory=1000,
        quiet=True,
        flags=import_flags,
    )
    if epsg is not None:
        # switch location
        os.environ["GISRC"] = str(TGTGISRC)
        # reproject data
        for i in range(1, 5):
            name = f"{raster_name}.{i}"
            grass.run_command(
                "r.proj",
                location=TMPLOC,
                mapset="PERMANENT",
                input=name,
                output=name,
                quiet=True,
            )


def main():
    global tmp_dir, temp_region, temp_region_2, rm_vec, rm_rast, rm_group, pid
    # parser options
    aoi_map = options["aoi_map"]
    # Note:
    # the gpkg-file below (URL_ortho_all) contains
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

    # vector maps
    vm_clip = f"vm_clip_{pid}"

    # read federal state(s) from input options
    if options["filepath"]:
        with open(f'{options["filepath"]}') as f:
            federal_states = f.read()
    else:
        federal_states = options["federal_state"]

    # save current region for setting back later in cleanup
    temp_region = f"temp_region_{pid}"
    grass.run_command("g.region", save=temp_region, overwrite=True, quiet=True)

    # create list for each raster band for building entire raster
    # of all given federal states
    all_raster_red = []
    all_raster_green = []
    all_raster_blue = []
    all_raster_NIR = []

    # preventing error caused by new line at the end of federal_states
    if federal_states.endswith("\n"):
        federal_states = federal_states[:-1]

    # loop through federal states and get respective DOPs
    for federal_state in federal_states.split(","):
        if federal_state in URLS:
            grass.message(f"Processing {federal_state}...")
            # create lists for all raster of one band
            raster_red = []
            raster_green = []
            raster_blue = []
            raster_NIR = []

            if federal_state == "Hessen":
                # create grid for wms import
                tiles_list = create_grid(1000, "HE_DOP", options["aoi_map"])
                # no tileindex
                tileindex = None

            else:
                tileindex = URLS[federal_state]
                if tileindex == None:
                    grass.warning(
                        _(f"{federal_state} is not yet implemented.")
                    )

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

        if tileindex:
            # if aoi_map given: set region to aoi_map extens
            if aoi_map:
                # import tileindex only for aoi
                # for this set region to aoi first
                grass.run_command(
                    "g.region",
                    vector=aoi_map,
                    res=0.2,
                    flags="a",
                    quiet=True,
                )

            # download and unzip tile index
            vm_import = get_tindex(tileindex)

            if aoi_map:
                # check which tiles are needed for selected AOI
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
            tiles = grass.vector_db_select(vm_clip, columns="location")[
                "values"
            ].items()

            # safe current region
            if not flags["r"]:
                temp_region_2 = f"temp_region_2_{pid}"
                grass.run_command(
                    "g.region", save=temp_region_2, overwrite=True, quiet=True
                )

                reg = grass.region()
                if reg["nsres"] == reg["ewres"]:
                    res = float(reg["nsres"])
                else:
                    grass.fatal(
                        "N/S resolution is not the same as E/W resolution!"
                    )

            num_imp = len(tiles)
            num = 0
            for key, URL_tile in tiles:
                num += 1
                raster_name = (
                    f"{os.path.basename(URL_tile[0]).split('.')[0].replace('-', '_')}"
                    f"_pid{pid}"
                )
                # import DOP tile with original resolution
                grass.message(
                    f"Started raster import {num}/{num_imp} from URL: {URL_tile[0]}"
                )
                # set memory manually to 1000
                # Process stuck, when memory is too large (100000)
                if federal_state == "Brandenburg":
                    import_and_reproject(URL_tile[0], raster_name, epsg=25833)
                else:
                    import_and_reproject(URL_tile[0], raster_name)

                # adjust resolution if required
                if not flags["r"]:
                    # set region to imported tile
                    grass.run_command(
                        "g.region",
                        raster=f"{raster_name}.1",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.region", res=res, flags="a", quiet=True
                    )

                    res_rast = float(
                        grass.parse_command(
                            "r.info", map=f"{raster_name}.1", flags="g"
                        )["nsres"]
                    )
                    res_region = float(
                        grass.parse_command("g.region", flags="g")["nsres"]
                    )
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
                # when importing DOPs a 'group' element
                # named 'raster_name' is created which should be removed
                # (in most cases already deleted in temp-location of
                # r.import; only for EPSG:25832 r.import reduces
                # to r.in.gdal (since DOPs given in EPGS:25832)
                # and the 'group' element remains in current mapset)
                rm_group.append(raster_name)

                grass.message("Finished raster import")

                grass.run_command("g.region", region=temp_region_2)

                rastername_red = f"{raster_name}_red"
                grass.run_command(
                    "g.rename",
                    raster=f"{raster_name}.1,{rastername_red}",
                    quiet=True,
                )
                raster_red.append(rastername_red)

                rastername_green = f"{raster_name}_green"
                grass.run_command(
                    "g.rename",
                    raster=f"{raster_name}.2,{rastername_green}",
                    quiet=True,
                )
                raster_green.append(rastername_green)

                rastername_blue = f"{raster_name}_blue"
                grass.run_command(
                    "g.rename",
                    raster=f"{raster_name}.3,{rastername_blue}",
                    quiet=True,
                )
                raster_blue.append(rastername_blue)

                rastername_NIR = f"{raster_name}_NIR"
                grass.run_command(
                    "g.rename",
                    raster=f"{raster_name}.4,{rastername_NIR}",
                    quiet=True,
                )
                raster_NIR.append(rastername_NIR)

            grass.run_command("g.region", region=temp_region_2)

            rm_rast += raster_red + raster_green + raster_blue + raster_NIR

        else:
            temp_region_2 = f"temp_region_HE_{pid}"
            # save current region first
            grass.run_command("g.region", save=temp_region_2)

            # no difference between with or without -r flag (HE WMS is res=0.2)
            for tile in tiles_list:
                # set region to respective HE tile
                grass.run_command("g.region", vector=tile)
                # get rgb and cir DOPs through HE WMS
                wms_HE = "https://www.gds-srv.hessen.de/cgi-bin/lika-services/ogc-free-images.ows?language=ger&"
                # get CIR DOP
                grass.run_command(
                    "r.in.wms",
                    url=wms_HE,
                    layer="he_dop20_cir",
                    output=f"cir_{tile}",
                    format="tiff",
                    flags="b",
                    overwrite=True,
                )
                # rename NIR band (called red, but its actually NIR in this case)
                grass.run_command(
                    "g.rename",
                    raster=f"cir_{tile}.red,{tile}.nir",
                    overwrite=True,
                    quiet=True,
                )
                # drop rest of CIR DOP
                grass.run_command(
                    "g.remove", type="raster", pattern="cir_*", flags="f"
                )
                # get RGB DOP
                grass.run_command(
                    "r.in.wms",
                    url=wms_HE,
                    layer="he_dop20_rgb",
                    output=tile,
                    format="tiff",
                    flags="b",
                    overwrite=True,
                )

                rm_group.append(f"cir_{tile}")
                rm_group.append(tile)

                grass.message("Finished raster import")

                rastername_red = f"{tile}_red"
                grass.run_command(
                    "g.rename",
                    raster=f"{tile}.red,{rastername_red}",
                    quiet=True,
                )
                raster_red.append(rastername_red)

                rastername_green = f"{tile}_green"
                grass.run_command(
                    "g.rename",
                    raster=f"{tile}.green,{rastername_green}",
                    quiet=True,
                )
                raster_green.append(rastername_green)

                rastername_blue = f"{tile}_blue"
                grass.run_command(
                    "g.rename",
                    raster=f"{tile}.blue,{rastername_blue}",
                    quiet=True,
                )
                raster_blue.append(rastername_blue)

                rastername_NIR = f"{tile}_NIR"
                grass.run_command(
                    "g.rename",
                    raster=f"{tile}.nir,{rastername_NIR}",
                    quiet=True,
                )
                raster_NIR.append(rastername_NIR)

            grass.run_command("g.region", region=temp_region_2)

            rm_rast += raster_red + raster_green + raster_blue + raster_NIR

        if federal_state in URLS:
            # if multiple tiles, build one raster out of all tiles for each band
            output_R = f'{options["output"]}_{FS[federal_state]}_R'
            rm_rast.append(output_R)
            output_G = f'{options["output"]}_{FS[federal_state]}_G'
            rm_rast.append(output_G)
            output_B = f'{options["output"]}_{FS[federal_state]}_B'
            rm_rast.append(output_B)
            output_NIR = f'{options["output"]}_{FS[federal_state]}_NIR'
            rm_rast.append(output_NIR)
            if tileindex:
                if len(tiles) > 1:
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_red,
                        output=output_R,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_green,
                        output=output_G,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_blue,
                        output=output_B,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_NIR,
                        output=output_NIR,
                        quiet=True,
                    )
                else:
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_red[0]},{output_R}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_green[0]},{output_G}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_blue[0]},{output_B}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_NIR[0]},{output_NIR}",
                        quiet=True,
                    )
            elif federal_state == "Hessen":
                if len(tiles_list) > 1:
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_red,
                        output=output_R,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_green,
                        output=output_G,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_blue,
                        output=output_B,
                        quiet=True,
                    )
                    grass.run_command(
                        "r.buildvrt",
                        input=raster_NIR,
                        output=output_NIR,
                        quiet=True,
                    )
                else:
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_red[0]},{output_R}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_green[0]},{output_G}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_blue[0]},{output_B}",
                        quiet=True,
                    )
                    grass.run_command(
                        "g.rename",
                        raster=f"{raster_NIR[0]},{output_NIR}",
                        quiet=True,
                    )

            # build raster for all given federal state(s) if there are >1 given
            if "," in federal_states:
                all_raster_red.append(output_R)
                all_raster_green.append(output_G)
                all_raster_blue.append(output_B)
                all_raster_NIR.append(output_NIR)

    if "," in federal_states:
        output_R = f'{options["output"]}_R'
        output_G = f'{options["output"]}_G'
        output_B = f'{options["output"]}_B'
        output_NIR = f'{options["output"]}_NIR'

        grass.run_command(
            "r.buildvrt", input=all_raster_red, output=output_R, quiet=True
        )
        grass.run_command(
            "r.buildvrt", input=all_raster_green, output=output_G, quiet=True
        )
        grass.run_command(
            "r.buildvrt", input=all_raster_blue, output=output_B, quiet=True
        )
        grass.run_command(
            "r.buildvrt", input=all_raster_NIR, output=output_NIR, quiet=True
        )

    output_Rp1 = f'{options["output"]}_red'
    output_Gp1 = f'{options["output"]}_green'
    output_Bp1 = f'{options["output"]}_blue'
    output_NIRp1 = f'{options["output"]}_nir'
    grass.run_command(
        "r.mapcalc",
        expression=f"{output_Rp1} = {output_R} + 1",
        quiet=True,
        region="intersect",
    )
    grass.run_command(
        "r.mapcalc",
        expression=f"{output_Gp1} = {output_G} + 1",
        quiet=True,
        region="intersect",
    )
    grass.run_command(
        "r.mapcalc",
        expression=f"{output_Bp1} = {output_B} + 1",
        quiet=True,
        region="intersect",
    )
    grass.run_command(
        "r.mapcalc",
        expression=f"{output_NIRp1} = {output_NIR} + 1",
        quiet=True,
        region="intersect",
    )
    grass.message(
        _(
            "Generated following raster maps:"
            f"{output_Rp1}, {output_Gp1}, "
            f"{output_Bp1}, {output_NIRp1}"
        )
    )

    if tileindex or aoi_map:
        grass.run_command("g.region", region=temp_region)


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
