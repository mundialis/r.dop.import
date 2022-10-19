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
sys.path.insert(1, os.path.join(os.path.dirname(sys.path[0]),
                                'etc',
                                'r.dop.import'))
from download_urls import URLS

tmp_dir = None
rm_vec = []
rm_rast = []
rm_group = []
temp_region = None


def cleanup():
    if tmp_dir:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

    nuldev = open(os.devnull, "w")
    kwargs = {"flags": "f", "quiet": True, "stderr": nuldev}
    for rmvec in rm_vec:
        if grass.find_file(name=rmvec, element="vector")["file"]:
            grass.run_command("g.remove",
                              type="vector",
                              name=rmvec,
                              **kwargs)
    for rmrast in rm_rast:
        if grass.find_file(name=rmrast, element="raster")["file"]:
            grass.run_command("g.remove",
                              type="raster",
                              name=rmrast,
                              **kwargs)
    for rmgroup in rm_group:
        if grass.find_file(name=rmgroup, element="group")["file"]:
            grass.run_command("g.remove",
                              type="group",
                              name=rmgroup,
                              **kwargs)

    if temp_region:
        # set region back and delete saved region:
        grass.run_command("g.region", region=temp_region)
        grass.run_command("g.remove",
                          type='region',
                          name=temp_region,
                          flags='f',
                          quiet=True)


def main():
    pid = os.getpid()
    global tmp_dir, temp_region, rm_vec, rm_rast, rm_group
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
    URL_ortho_all = ("https://github.com/mundialis/openNRW/raw/master/dop/"
                     "openNRW_DOP10_tileindex.gpkg.gz")
    zipname = os.path.basename(URL_ortho_all)
    unzipped_name = zipname.replace(".gz", "")
    tmp_dir = grass.tempdir()
    download_path = os.path.join(tmp_dir, zipname)
    unzipped_path = os.path.join(tmp_dir, unzipped_name)
    # vector maps
    vm_import = f'vm_import_{pid}'
    vm_clip = f'vm_clip_{pid}'

    # read bundesland file
    if options['filepath']:
        with open(f'{options["filepath"]}') as f:
            federal_states = f.read()
    else:
        federal_states = options['federal_state']

    # get URL for corresponding federal state(s)
    URL = None
    # fs = None
    for federal_state in federal_states.split(','):
        if federal_state in URLS:
            if federal_state == 'Nordrhein-Westfalen':
                URL = URLS[federal_state]
                # fs = federal_state
            else:
                grass.warning(_(f"{federal_state} is not yet implemented."))
        else:
            if options['filepath']:
                grass.fatal(_("Non valid name of federal state,"
                              " in 'filepath'-option given"))
            elif options['federal_state']:
                grass.fatal(_("Non valid name of federal state,"
                              " in 'federal_states'-option given"))
    if not URL:
        grass.fatal(_("Given federal state(s) are not yet implemented."))
    # so far, just NRW implemented
    # in case single federal state given, and not NRW:
    #   skips following part
    #   + grass.message: see above
    # in case multiple federal states given, and at least one of them is NRW:
    #   import data only for NRW area

    if URL:
        # if aoi_map given: set region to aoi_map extens
        if aoi_map:
            # import tileindex only for aoi
            # for this first set region to aoi
            # region
            temp_region = f'temp_region_{pid}'
            # save current region for setting back later in cleanup
            grass.run_command("g.region", save=temp_region, quiet=True)
            # set region to aoi_map
            grass.run_command("g.region",
                              vector=aoi_map, flags='a', quiet=True)

        # download and unzip tile index
        wget.download(URL_ortho_all, download_path, bar=None)
        os.system(f"gunzip {download_path}")
        # import vector map containing URL for each tile
        grass.run_command("v.import",
                          input=unzipped_path,
                          output=vm_import,
                          extent='region',
                          quiet=True)
        rm_vec.append(vm_import)

        if aoi_map:
            # check which tiles are needed for selected AOI
            grass.run_command("v.clip",
                              input=vm_import,
                              clip=aoi_map,
                              output=vm_clip,
                              flags='d',
                              quiet=True)
            rm_vec.append(vm_clip)
        else:
            # if no aoi given, use complete (current set) region:
            vm_clip = vm_import

        # create lists for all raster of one band
        raster_red = []
        raster_green = []
        raster_blue = []
        raster_NIR = []
        # import tiles and rename them according to their band
        # and write them in a list
        tiles = grass.vector_db_select(vm_clip,
                                       columns="location"
                                       )["values"].items()
        for key, URL_tile in tiles:
            raster_name = (f"{os.path.basename(URL_tile[0]).split('.')[0]}"
                           f"_pid{pid}")
            # import DOP tile with original resolution
            kwargs = {"extent": "region"}
            print(f"Started raster import for key: {key} and URL: {URL_tile}")
            # set memory manually to 1000
            # Process stuck, when memory is too large (100000)
            grass.run_command("r.import",
                              input=URL_tile[0],
                              output=raster_name,
                              memory=1000,
                              quiet=True,
                              **kwargs)
            # adjust resolution if required
            if not flags["r"]:
                res_rast = float(grass.parse_command("r.info",
                                                     map=f"{raster_name}.1",
                                                     flags='g')['nsres'])
                res_region = float(grass.parse_command("g.region",
                                                       flags='g')['nsres'])
                if res_rast > res_region:
                    for band in range(1, 5):
                        grass.run_command("r.resamp.interp",
                                          input=f"{raster_name}.{band}",
                                          output=f"{raster_name}.{band}",
                                          overwrite=True)
                elif res_rast < res_region:
                    for band in range(1, 5):
                        grass.run_command("r.resamp.stats",
                                          input=f"{raster_name}.{band}",
                                          output=f"{raster_name}.{band}",
                                          overwrite=True)
            # when importing DOPs a 'group' element
            # named 'raster_name' is created which should be removed
            # (in most cases already deleted in temp-location of
            # r.import; only for EPSG:25832 r.import reduces
            # to r.in.gdal (since DOPs given in EPGS:25832)
            # and the 'group' element remains in current mapset)
            rm_group.append(raster_name)

            print("Finished raster import")

            rastername_red = f"{raster_name}_red"
            grass.run_command("g.rename",
                              raster=f"{raster_name}.1,{rastername_red}",
                              quiet=True)
            raster_red.append(rastername_red)

            rastername_green = f"{raster_name}_green"
            grass.run_command("g.rename",
                              raster=f"{raster_name}.2,{rastername_green}",
                              quiet=True)
            raster_green.append(rastername_green)

            rastername_blue = f"{raster_name}_blue"
            grass.run_command("g.rename",
                              raster=f"{raster_name}.3,{rastername_blue}",
                              quiet=True)
            raster_blue.append(rastername_blue)

            rastername_NIR = f"{raster_name}_NIR"
            grass.run_command("g.rename",
                              raster=f"{raster_name}.4,{rastername_NIR}",
                              quiet=True)
            raster_NIR.append(rastername_NIR)

        rm_rast += raster_red+raster_green+raster_blue+raster_NIR

        # if multiple tiles, build one raster out of all tiles for each band
        output_R = f'{options["output"]}_R'
        rm_rast.append(output_R)
        output_G = f'{options["output"]}_G'
        rm_rast.append(output_G)
        output_B = f'{options["output"]}_B'
        rm_rast.append(output_B)
        output_NIR = f'{options["output"]}_NIR'
        rm_rast.append(output_NIR)
        if len(tiles) > 1:
            grass.run_command(
                "r.buildvrt", input=raster_red, output=output_R, quiet=True
            )
            grass.run_command(
                "r.buildvrt", input=raster_green, output=output_G, quiet=True
            )
            grass.run_command(
                "r.buildvrt", input=raster_blue, output=output_B, quiet=True
            )
            grass.run_command(
                "r.buildvrt", input=raster_NIR, output=output_NIR, quiet=True
            )
        else:
            grass.run_command("g.rename",
                              raster=f"{raster_red[0]},{output_R}",
                              quiet=True)
            grass.run_command("g.rename",
                              raster=f"{raster_green[0]},{output_G}",
                              quiet=True)
            grass.run_command("g.rename",
                              raster=f"{raster_blue[0]},{output_B}",
                              quiet=True)
            grass.run_command("g.rename",
                              raster=f"{raster_NIR[0]},{output_NIR}",
                              quiet=True)
        # value range from [0 255] to [1 256]:
        output_Rp1 = f'{options["output"]}_red'
        output_Gp1 = f'{options["output"]}_green'
        output_Bp1 = f'{options["output"]}_blue'
        output_NIRp1 = f'{options["output"]}_nir'
        grass.run_command(
            "r.mapcalc",
            expression=f'{output_Rp1} = {output_R} + 1',
            quiet=True,
            region='intersect'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{output_Gp1} = {output_G} + 1',
            quiet=True,
            region='intersect'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{output_Bp1} = {output_B} + 1',
            quiet=True,
            region='intersect'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{output_NIRp1} = {output_NIR} + 1',
            quiet=True,
            region='intersect'
        )
        grass.message(_("Generated following raster maps:"
                        f"{output_Rp1}, {output_Gp1}, "
                        f"{output_Bp1}, {output_NIRp1}"))


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
