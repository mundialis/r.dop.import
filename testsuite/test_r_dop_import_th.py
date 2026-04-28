#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test Thüringen
# AUTHOR(S):   Lina Krisztian, Anika Weinmann
# PURPOSE:     Tests r.dop.import Thüringen
# SPDX-FileCopyrightText: (c) 2023 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
#############################################################################

from grass.gunittest.main import test
from r_dop_import_test_base import RDopImportTestBase


class TestRDopImportSN(RDopImportTestBase):
    """Test class for r.dop.import for TH"""

    fs = "TH"
    ref_res = 0.2
    aoi_cells = 389351

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


if __name__ == "__main__":
    test()
