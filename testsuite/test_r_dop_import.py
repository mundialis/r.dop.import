#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test
# AUTHOR(S):   Lina Krisztian

# PURPOSE:      Tests r.dop.import
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

import os
import sys

from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.gunittest.gmodules import SimpleModule
import grass.script as grass


class TestRDopImport(TestCase):
    pid = os.getpid()
    rm_vec = []
    region_b_file = os.path.join('data', 'beuel_two_tiles_1.gpkg')
    region_b = f'region_b_{pid}'
    aoi_map_file = os.path.join('data', 'beuel_two_tiles_2.gpkg')
    aoi_map = f'aoi_map_{pid}'
    rm_vec.append(aoi_map)
    rm_vec.append(region_b)
    aoi_map_single_tile_file = os.path.join('data', 'beuel_single_tile.gpkg')
    aoi_map_single_tile = f'aoi_map_single_tile_{pid}'
    rm_vec.append(aoi_map_single_tile)
    # name of output, which r.dop.import creates
    test_output = f'output_{pid}'
    test_output_all = [f"{out}_{band}"
                       for band, out
                       in zip(['red', 'green', 'blue', 'nir'],
                              4*[test_output])]

    @classmethod
    def setUpClass(self):
        """Ensures expected computational region and generated data"""
        # perform test in location with meter
        # needed e.g. for res-check within tests
        proj_unit = grass.parse_command("g.proj",
                                        flags='g')['unit']
        if proj_unit != 'meter':
            sys.exit('WARNING: Tests were skipped. Please perform the test '
                     'in a location with units of meter.')
        proj_epsg = grass.parse_command("g.proj",
                                        flags='g')['srid']
        if proj_epsg != 'EPSG:32632' and proj_epsg != 'EPSG:25832':
            print("WARNING: tests might fail, due to rounding problems "
                  "when importing/reprojecting test data in certain locations"
                  "(e.g. import 701 instead of 700 rows)."
                  "Tests should run successfully (at least) in location: "
                  "EPSG:32632 and EPSG:25832.")
        # set region
        self.runModule("v.import",
                       input=self.region_b_file,
                       output=self.region_b)
        self.runModule("g.region",
                       vector=self.region_b,
                       res=1.,
                       flags='a')
        # import aoi_map for testing
        self.runModule("v.import",
                       input=self.aoi_map_file,
                       output=self.aoi_map)
        # import aoi_map containing single tile DOP
        self.runModule("v.import",
                       input=self.aoi_map_single_tile_file,
                       output=self.aoi_map_single_tile)

    @classmethod
    def tearDownClass(self):
        """Remove the temporary region and generated data"""
        for vec in self.rm_vec:
            self.runModule("g.remove",
                           type="vector",
                           name=vec,
                           flags='f')

    def tearDown(self):
        """Remove the outputs created
        This is executed after each test run.
        """
        self.runModule("g.remove",
                       type="raster",
                       pattern=f"{self.test_output}*",
                       flags="f")

    def test_default_settings(self):
        """
        Tests module with default settings:
        importing data for current set region
        with resolution of current set region
        """
        r_check = SimpleModule("r.dop.import",
                               output=self.test_output,
                               federal_state='Nordrhein-Westfalen')
        self.assertModule(r_check,
                          "Importing data for current set region fails.")
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el,
                                    msg=f"Creation of {el} failed.")
            # test if range of values are correct
            raster_range = grass.parse_command("r.info",
                                               map=el,
                                               flags='r')
            self.assertTrue(((float(raster_range['min']) >= 1)
                            and (float(raster_range['max']) <= 256)),
                            f"Incorrect value range of raster map {el}."
                            "Should be [1 256]")
        # test resolution
        out_res_ns = round(float(
                           grass.parse_command("r.info",
                                               map=self.test_output_all[0],
                                               flags='g')['nsres']))
        out_res_ew = round(float(
                           grass.parse_command("r.info",
                                               map=self.test_output_all[0],
                                               flags='g')['ewres']))
        self.assertTrue(((out_res_ns == 1.) and (out_res_ew == 1.)),
                        "Resolution of DOPs is not as expected."
                        f"Expected res of 1, but got nsres of {out_res_ns}"
                        f"and ewres of {out_res_ew}")

    def test_extent_aoi_map(self):
        """
        Tests import dating only for AOI given by aoi_map
        """
        r_check = SimpleModule("r.dop.import",
                               output=self.test_output,
                               federal_state='Nordrhein-Westfalen',
                               aoi_map=self.aoi_map)
        self.assertModule(r_check,
                          "Importind data for AOI fails.")
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el,
                                    msg=f"Creation of {el} failed.")
        # check if just aoi data loaded ==>
        # should have 424526 cells
        cells_aoi = grass.parse_command("r.info",
                                        map=self.test_output_all[0],
                                        flags='g')['cells']
        self.assertTrue(int(cells_aoi) == 424526,
                        "Data loaded for larger region than AOI")

    def test_single_tile(self):
        """
        Tests importing data, containing only single DOP tile
        """
        r_check = SimpleModule("r.dop.import",
                               output=self.test_output,
                               federal_state='Nordrhein-Westfalen',
                               aoi_map=self.aoi_map_single_tile)
        self.assertModule(r_check,
                          "Importing single DOP tile fails.")
        # import should have rows=700 & cols=801 ==> cells=560700
        cells_aoi = grass.parse_command("r.info",
                                        map=self.test_output_all[0],
                                        flags='g')['cells']
        self.assertTrue(int(cells_aoi) == 560700,
                        "Data loaded for larger region than AOI (single DOP)")

    def test_dop_resolution(self):
        """
        Tests importing data with original resolution of DOPs
        """
        r_check = SimpleModule("r.dop.import",
                               output=self.test_output,
                               aoi_map=self.aoi_map_single_tile,
                               federal_state='Nordrhein-Westfalen',
                               flags='r')
        self.assertModule(r_check,
                          "Importing data with original DOP resolution fails."
                          "Note: could also be a fail of aoi_map import "
                          "(test_extent_aoi_map) or single DOP tile import "
                          "(test_single_tile). Check them first.")
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el,
                                    msg=f"Creation of {el} failed.")
        out_res_ns = round(float(
                           grass.parse_command("r.info",
                                               map=self.test_output_all[0],
                                               flags='g')['nsres']), 2)
        out_res_ew = round(float(
                           grass.parse_command("r.info",
                                               map=self.test_output_all[0],
                                               flags='g')['ewres']), 2)
        self.assertTrue(((out_res_ns == 0.1) and (out_res_ew == 0.1)),
                        "Resolution of DOPs is not as expected."
                        f"Expected res of 0.1, but got nsres of {out_res_ns}"
                        f" and ewres of {out_res_ew}")


if __name__ == "__main__":
    test()
