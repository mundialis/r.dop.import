#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r_dop_import_lib
# AUTHOR(S):   Johannes Halbauer
#
# PURPOSE:     Library for r.dem.import
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

import grass.script as grass

from grass_gis_helpers.general import set_nprocs

OPEN_DATA_AVAILABILITY = {
    "NO_OPEN_DATA": ["BW", "BY", "HB", "HH", "MV", "RP", "SL", "ST", "SH"],
    "NOT_YET_SUPPORTED": ["TH", "HE", "NI"],
    "SUPPORTED": ["BE", "BB", "NW", "SN"],
}


def setup_parallel_processing(nprocs):
    nprocs = set_nprocs(nprocs)
    # set some common environmental variables, like:
    os.environ.update(
        dict(
            GRASS_COMPRESSOR="LZ4",
            GRASS_MESSAGE_FORMAT="plain",
        )
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
