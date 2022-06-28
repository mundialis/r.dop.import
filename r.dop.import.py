#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import
# AUTHOR(S):   Johannes Halbauer
#
# PURPOSE:     Downloads Digital Orthophotos (DOPs) within a specified area
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
# % description: Downloads Digital Othophotos (DOPs) within a specified area.
# % keyword: imagery
# % keyword: download
# % keyword: DOP
# %end

# %option G_OPT_V_INPUT
# % key: input
# % required: yes
# % description: Existing vector map in GRASS as input
# %end

# %option
# % key: filepath
# % type: string
# % required: yes
# % description: Path to .txt file containing the federal state(s) related to the area of interest
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % required: yes
# % description: Name of resulting raster(s)
# %end


import atexit
import wget
import gzip
import os
import requests
from io import BytesIO
import shutil
import grass.script as grass

options, flags = grass.parser()


AOI = options["input"]

tmp_dir = grass.tempdir()

# variables for tile index download
URL = "https://github.com/mundialis/openNRW/raw/master/dop/openNRW_DOP10_tileindex.gpkg.gz"
zipname = os.path.basename(URL)
unzipped_name = zipname.replace(".gz", "")
download_path = os.path.join(tmp_dir, zipname)
unzipped_path = os.path.join(tmp_dir, unzipped_name)


def cleanup():
    if download_path:
        if os.path.exists(download_path):
            os.remove(download_path)

    if unzipped_path:
        if os.path.exists(unzipped_path):
            os.remove(unzipped_path)

    if tmp_dir:
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

    if grass.find_file(name="necessary_tiles", element="vector")["file"]:
        grass.run_command("g.remove", type="vector", name="necessary_tiles", flags="f")

    if grass.find_file(name="tile_index", element="vector")["file"]:
        grass.run_command("g.remove", type="vector", name="tile_index", flags="f")

    if grass.find_file(name=f'{options["output"]}_R', element="raster")["file"]:
        grass.run_command(
            "g.remove", type="raster", name=f'{options["output"]}_R', flags="f"
        )

    if grass.find_file(name=f'{options["output"]}_G', element="raster")["file"]:
        grass.run_command(
            "g.remove", type="raster", name=f'{options["output"]}_G', flags="f"
        )

    if grass.find_file(name=f'{options["output"]}_B', element="raster")["file"]:
        grass.run_command(
            "g.remove", type="raster", name=f'{options["output"]}_B', flags="f"
        )

    if grass.find_file(name=f'{options["output"]}_NIR', element="raster")["file"]:
        grass.run_command(
            "g.remove", type="raster", name=f'{options["output"]}_NIR', flags="f"
        )


def main():

    # read bundesland file
    with open(f'{options["filepath"]}') as f:
        bundesland = f.read()

    # check if AOI is located in more than one federal state
    if "," in bundesland:
        grass.message(_("Area of interest is located in more than one federal state."))
        if "Nordrhein-Westfalen" in bundesland:
            grass.message(
                _(
                    "Only the digital orthophotos of the part within North Rhine-Westphalia will be imported."
                )
            )

    if "Nordrhein-Westfalen" in bundesland:

        # download and unzip tile index
        global URL

        tiles_zipped = wget.download(URL, download_path)
        os.system("gunzip " + download_path)

        # import tileindex
        grass.run_command("v.import", input=unzipped_path, output="tile_index")

        # check which tiles are needed for selected AOI
        grass.run_command(
            "v.overlay",
            ainput=AOI,
            binput="tile_index",
            operator="and",
            output="necessary_tiles",
        )

        # create lists for all raster of one band
        raster_red = []
        raster_green = []
        raster_blue = []
        raster_NIR = []

        # import tiles and rename them according to their band and write them in a list
        for key, URL in grass.vector_db_select(
            ("necessary_tiles"), columns="b_location"
        )["values"].items():

            raster_name = os.path.basename(URL[0]).split(".")[0]

            grass.run_command("r.import", input=(URL[0]), output=raster_name)

            rastername_red = f"{raster_name}_red"
            grass.run_command("g.rename", raster=f"{raster_name}.1,{rastername_red}")
            raster_red.append(rastername_red)

            rastername_green = f"{raster_name}_green"
            grass.run_command("g.rename", raster=f"{raster_name}.2,{rastername_green}")
            raster_green.append(rastername_green)

            rastername_blue = f"{raster_name}_blue"
            grass.run_command("g.rename", raster=f"{raster_name}.3,{rastername_blue}")
            raster_blue.append(rastername_blue)

            rastername_NIR = f"{raster_name}_NIR"
            grass.run_command("g.rename", raster=f"{raster_name}.4,{rastername_NIR}")
            raster_NIR.append(rastername_NIR)

        # build one raster out of all tiles for each band
        grass.run_command(
            "r.buildvrt", input=raster_red, output=f'{options["output"]}_R'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{options["output"]}_R_1_256 = {options["output"]}_R + 1',
        )

        grass.run_command(
            "r.buildvrt", input=raster_green, output=f'{options["output"]}_G'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{options["output"]}_G_1_256 = {options["output"]}_G + 1',
        )

        grass.run_command(
            "r.buildvrt", input=raster_blue, output=f'{options["output"]}_B'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{options["output"]}_B_1_256 = {options["output"]}_B + 1',
        )

        grass.run_command(
            "r.buildvrt", input=raster_NIR, output=f'{options["output"]}_NIR'
        )
        grass.run_command(
            "r.mapcalc",
            expression=f'{options["output"]}_NIR_1_256 = {options["output"]}_NIR + 1',
        )

        grass.message(_("Done."))
    else:
        grass.fatal(
            _(
                "Sorry, area of interest is not located in North Rhine-Westphalia. Digital orthophotos can not be downloaded for other federal states yet."
            )
        )


if __name__ == "__main__":
    options, flags = grass.parser()
    atexit.register(cleanup)
    main()
