#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test Brandenburg and Sachsen
# AUTHOR(S):   Anika Weinmann

# PURPOSE:     Tests r.dop.import Brandenburg and Sachsen
# COPYRIGHT:   (C) 2024 by mundialis GmbH & Co. KG and the GRASS Development
#              Team
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

from grass.gunittest.main import test
import grass.script as grass
from r_dop_import_test_base import RDopImportTestBase


class TestRDopImportBBSN(RDopImportTestBase):
    fs = "BB,SN"
    ref_res = 0.2
    aoi_cells = 312018

    def test_default_settings(self):
        """
        Tests module with default settings:
        importing data for current set region
        with resolution of current set region
        """
        grass.run_command(
            "g.region",
            vector=self.aoi_map,
            res=1,
            flags="a",
        )
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


if __name__ == "__main__":
    test()
