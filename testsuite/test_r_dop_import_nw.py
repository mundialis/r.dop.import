#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test Nordrhein-Westfalen
# AUTHOR(S):   Lina Krisztian
#
# PURPOSE:     Tests r.dop.import Nordrhein-Westfalen
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

import os
from grass.gunittest.main import test
from grass.gunittest.gmodules import SimpleModule
import grass.script as grass

from r_dop_import_test_base import RDopImportTestBase


class TestRDopImportNW(RDopImportTestBase):
    """Test class for r.dop.import for NW"""

    fs = "NW"
    ref_res = 0.1
    aoi_cells = 424526

    def test_default_settings(self):
        """
        Tests module with default settings:
        importing data for current set region
        with resolution of current set region
        """
        self.default_settings_test()

    def test_extent_aoi_map(self):
        """
        Tests importing data only for AOI given by aoi_map
        """
        self.extent_aoi_map_test()

    def test_dop_resolution(self):
        """
        Tests importing data with original resolution of DOPs
        """
        self.dop_resolution_test()

    def test_single_tile(self):
        """
        Tests importing data, containing only single DOP tile
        """
        cells = 560700
        print("\nTest import of singel tile (NRW) ...")
        aoi_map_single_tile = "aoi_map_single_tile"
        self.rm_vec.append(aoi_map_single_tile)
        grass.run_command(
            "v.import",
            input=os.path.join("data", "beuel_single_tile.gpkg"),
            output=aoi_map_single_tile,
        )
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state="Nordrhein-Westfalen",
            aoi=aoi_map_single_tile,
        )
        self.assertModule(r_check, "Importing single DOP tile fails.")
        # import should have rows=700 & cols=801 ==> cells=560700
        cells_aoi = grass.parse_command(
            "r.info",
            map=self.test_output_all[0],
            flags="g",
        )["cells"]
        self.assertTrue(
            int(cells_aoi) == cells,
            "Data loaded for larger region than AOI (single DOP)",
        )
        print("Test import of singel tile (NRW) successfully finished.\n")


if __name__ == "__main__":
    test()
