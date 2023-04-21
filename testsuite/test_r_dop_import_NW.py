#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test Nordrhein-Westfalen
# AUTHOR(S):   Lina Krisztian

# PURPOSE:     Tests r.dop.import Nordrhein-Westfalen
# COPYRIGHT:   (C) 2022 by mundialis GmbH & Co. KG and the GRASS Development
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

import os
import sys

from grass.gunittest.case import TestCase
from grass.gunittest.main import test
from grass.gunittest.gmodules import SimpleModule
import grass.script as grass


class TestRDopImport(TestCase):
    pid = os.getpid()
    rm_vec = []
    region_b_file = os.path.join("data", "beuel_two_tiles_1.gpkg")
    region_b = f"region_b_{pid}"
    aoi_map_file = os.path.join("data", "beuel_two_tiles_2.gpkg")
    aoi_map = f"aoi_map_{pid}"
    rm_vec.append(aoi_map)
    rm_vec.append(region_b)
    aoi_map_single_tile_file = os.path.join("data", "beuel_single_tile.gpkg")
    aoi_map_single_tile = f"aoi_map_single_tile_{pid}"
    rm_vec.append(aoi_map_single_tile)
    # name of output, which r.dop.import creates
    test_output = f"output_{pid}"
    test_output_all = [
        f"{out}_{band}"
        for band, out in zip(
            ["red", "green", "blue", "nir"], 4 * [test_output]
        )
    ]
    GISDBASE = None
    TGTGISRC = None
    TMPLOC = None
    SRCGISRC = None

    @classmethod
    def delete_tmp_location(cls):
        # remove temp location
        if cls.TMPLOC:
            grass.try_rmdir(os.path.join(cls.GISDBASE, cls.TMPLOC))
        if cls.SRCGISRC:
            grass.try_remove(cls.SRCGISRC)

    @classmethod
    def get_actual_location(cls):
        # get actual location, mapset, ...
        grassenv = grass.gisenv()
        tgtloc = grassenv["LOCATION_NAME"]
        tgtmapset = grassenv["MAPSET"]
        cls.GISDBASE = grassenv["GISDBASE"]
        cls.TGTGISRC = os.environ["GISRC"]
        return tgtloc, tgtmapset

    @classmethod
    def createTMPlocation(cls, epsg=4326):
        """Switch location"""
        SRCGISRC = grass.tempfile()
        cls.TMPLOC = f"temp_import_location_{grass.tempname(8)}"
        f = open(SRCGISRC, "w")
        f.write("MAPSET: PERMANENT\n")
        f.write("GISDBASE: %s\n" % cls.GISDBASE)
        f.write("LOCATION_NAME: %s\n" % cls.TMPLOC)
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
            "g.proj", flags="c", location=cls.TMPLOC, quiet=True, **epsg_arg
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

    @classmethod
    def setUpClass(cls):
        """Ensures expected computational region and generated data"""
        # switch location
        cls.get_actual_location()
        cls.createTMPlocation(25832)

        # perform test in location with meter
        # needed e.g. for res-check within tests
        proj_unit = grass.parse_command("g.proj", flags="g")["unit"]
        if proj_unit != "meter":
            sys.exit(
                "WARNING: Tests were skipped. Please perform the test "
                "in a location with units of meter."
            )
        proj_epsg = grass.parse_command("g.proj", flags="g")["srid"]
        if proj_epsg != "EPSG:32632" and proj_epsg != "EPSG:25832":
            print(
                "WARNING: tests might fail, due to rounding problems "
                "when importing/reprojecting test data in certain locations"
                "(e.g. import 701 instead of 700 rows)."
                "Tests should run successfully (at least) in location: "
                "EPSG:32632 and EPSG:25832."
            )
        # set region
        cls.runModule(
            "v.import",
            input=cls.region_b_file,
            output=cls.region_b,
        )
        cls.runModule(
            "g.region",
            vector=cls.region_b,
            res=1.0,
            flags="a",
        )
        # import aoi_map for testing
        cls.runModule(
            "v.import",
            input=cls.aoi_map_file,
            output=cls.aoi_map,
        )
        # import aoi_map containing single tile DOP
        cls.runModule(
            "v.import",
            input=cls.aoi_map_single_tile_file,
            output=cls.aoi_map_single_tile,
        )

    @classmethod
    def tearDownClass(cls):
        """Remove the temporary region and generated data"""
        for vec in cls.rm_vec:
            cls.runModule(
                "g.remove",
                type="vector",
                name=vec,
                flags="f",
            )
        # switch location
        if cls.TGTGISRC:
            os.environ["GISRC"] = str(cls.TGTGISRC)

    def tearDown(self):
        """Remove the outputs created
        This is executed after each test run.
        """
        self.runModule(
            "g.remove",
            type="raster",
            pattern=f"{self.test_output}*",
            flags="f",
        )

    def test_default_settings(self):
        """
        Tests module with default settings:
        importing data for current set region
        with resolution of current set region
        """
        print("\nTest default settings ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state="Nordrhein-Westfalen",
        )
        self.assertModule(
            r_check, "Importing data for current set region fails."
        )
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el, msg=f"Creation of {el} failed.")
            # test if range of values are correct
            raster_range = grass.parse_command("r.info", map=el, flags="r")
            self.assertTrue(
                (
                    (float(raster_range["min"]) >= 1)
                    and (float(raster_range["max"]) <= 256)
                ),
                f"Incorrect value range of raster map {el}."
                "Should be [1 256]",
            )
        # test resolution
        out_res_ns = round(
            float(
                grass.parse_command(
                    "r.info", map=self.test_output_all[0], flags="g"
                )["nsres"]
            )
        )
        out_res_ew = round(
            float(
                grass.parse_command(
                    "r.info", map=self.test_output_all[0], flags="g"
                )["ewres"]
            )
        )
        self.assertTrue(
            ((out_res_ns == 1.0) and (out_res_ew == 1.0)),
            "Resolution of DOPs is not as expected."
            f"Expected res of 1, but got nsres of {out_res_ns}"
            f"and ewres of {out_res_ew}",
        )
        print("Test default settings successfully finished.\n")

    def test_extent_aoi_map(self):
        """
        Tests importing data only for AOI given by aoi_map
        """
        print("\nTest AOI (NRW) ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state="Nordrhein-Westfalen",
            aoi_map=self.aoi_map,
        )
        self.assertModule(r_check, "Importind data for AOI fails.")
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el, msg=f"Creation of {el} failed.")
        # check if just aoi data loaded ==>
        # should have 424526 cells
        cells_aoi = grass.parse_command(
            "r.info", map=self.test_output_all[0], flags="g"
        )["cells"]
        self.assertTrue(
            int(cells_aoi) == 424526, "Data loaded for larger region than AOI"
        )
        print("Test AOI (NRW) successfully finished.\n")

    def test_single_tile(self):
        """
        Tests importing data, containing only single DOP tile
        """
        print("\nTest import of singel tile (NRW) ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state="Nordrhein-Westfalen",
            aoi_map=self.aoi_map_single_tile,
        )
        self.assertModule(r_check, "Importing single DOP tile fails.")
        # import should have rows=700 & cols=801 ==> cells=560700
        cells_aoi = grass.parse_command(
            "r.info", map=self.test_output_all[0], flags="g"
        )["cells"]
        self.assertTrue(
            int(cells_aoi) == 560700,
            "Data loaded for larger region than AOI (single DOP)",
        )
        print("Test import of singel tile (NRW) successfully finished.\n")

    def test_dop_resolution(self):
        """
        Tests importing data with original resolution of DOPs
        """
        print("\nTest resolution (NRW) ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            aoi_map=self.aoi_map_single_tile,
            federal_state="Nordrhein-Westfalen",
            flags="r",
        )
        self.assertModule(
            r_check,
            "Importing data with original DOP resolution fails."
            "Note: could also be a fail of aoi_map import "
            "(test_extent_aoi_map) or single DOP tile import "
            "(test_single_tile). Check them first.",
        )
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el, msg=f"Creation of {el} failed.")
        out_res_ns = round(
            float(
                grass.parse_command(
                    "r.info", map=self.test_output_all[0], flags="g"
                )["nsres"]
            ),
            2,
        )
        out_res_ew = round(
            float(
                grass.parse_command(
                    "r.info", map=self.test_output_all[0], flags="g"
                )["ewres"]
            ),
            2,
        )
        # 10 cm auflösung und AOI
        self.assertTrue(
            ((out_res_ns == 0.1) and (out_res_ew == 0.1)),
            "Resolution of DOPs is not as expected."
            f"Expected res of 0.1, but got nsres of {out_res_ns}"
            f" and ewres of {out_res_ew}",
        )
        print("Test resolution (NRW) successfully finished.\n")


if __name__ == "__main__":
    test()
