#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r.dop.import test base
# AUTHOR(S):   Anika Weinmann
#
# PURPOSE:     Test base for r.dop.import
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

import os
import sys

from grass.gunittest.case import TestCase
from grass.gunittest.gmodules import SimpleModule
import grass.script as grass

from grass_gis_helpers.cleanup import cleaning_tmp_location
from grass_gis_helpers.location import (
    create_tmp_location,
    get_current_location,
)


class RDopImportTestBase(TestCase):
    """Base test class for r.dop.import"""

    fs = ""
    ref_res = None
    aoi_cells = None
    pid = os.getpid()
    aoi_map = f"aoi_map_{pid}"
    orig_region = f"orig_region_{pid}"

    # name of output, which r.dop.import creates
    test_output = f"output_{pid}"
    test_output_all = [
        f"{out}_{band}"
        for band, out in zip(
            ["red", "green", "blue", "nir"],
            4 * [test_output],
        )
    ]
    rm_vec = []
    rm_vec.append(aoi_map)

    ORIG_GISRC = None
    TMP_LOC = None
    GISDBASE = None
    TMP_GISRC = None

    @classmethod
    # pylint: disable=invalid-name
    def setUpClass(cls):
        """Ensures expected computational region and generated data"""
        if cls.fs != "":
            # switch to location with EPSG code 25832
            loc, mapset, cls.GISDBASE, cls.ORIG_GISRC = get_current_location()
            if cls.TMP_LOC is None:
                cls.TMP_LOC, cls.TMP_GISRC = create_tmp_location(epsg=25832)

        # perform test in location with meter
        # needed e.g. for res-check within tests
        proj_unit = grass.parse_command("g.proj", flags="g")["unit"]
        if proj_unit != "meter":
            sys.exit(
                "WARNING: Tests were skipped. Please perform the test "
                "in a location with units of meter.",
            )
        proj_epsg = grass.parse_command("g.proj", flags="g")["srid"]
        if proj_epsg not in ("EPSG:32632", "EPSG:25832"):
            print(
                "WARNING: tests might fail, due to rounding problems "
                "when importing/reprojecting test data in certain locations"
                "(e.g. import 701 instead of 700 rows)."
                "Tests should run successfully (at least) in location: "
                "EPSG:32632 and EPSG:25832.",
            )
        # import aoi_map for testing
        fs = cls.fs.replace(",", "_")
        cls.runModule(
            "v.import",
            input=os.path.join("data", f"test_aoi_{fs}.gpkg"),
            output=cls.aoi_map,
        )
        # set region
        grass.run_command("g.region", save=cls.orig_region)
        grass.run_command("g.region", vector=cls.aoi_map, res=1, flags="a")
        grass.run_command("g.region", n="n+200", s="n-100", w="e-100")

    @classmethod
    # pylint: disable=invalid-name
    def tearDownClass(cls):
        """Remove the temporary region and generated data"""
        for vec in cls.rm_vec:
            cls.runModule(
                "g.remove",
                type="vector",
                name=vec,
                flags="f",
            )
        if cls.fs != "":
            grass.run_command("g.region", region=cls.orig_region)
            grass.run_command(
                "g.remove",
                type="region",
                name=cls.orig_region,
                flags="f",
            )
            # switch location and remove temp location
            cleaning_tmp_location(
                cls.ORIG_GISRC,
                cls.TMP_LOC,
                cls.GISDBASE,
                cls.TMP_GISRC,
            )

    # pylint: disable=invalid-name
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

    def default_settings_test(self):
        """
        Tests module with default settings:
        importing data for current set region
        with resolution of current set region
        """
        print("\nTest default settings ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state=self.fs.split(","),
        )
        self.assertModule(
            r_check,
            "Importing data for current set region fails.",
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
                    "r.info",
                    map=self.test_output_all[0],
                    flags="g",
                )["nsres"],
            ),
        )
        out_res_ew = round(
            float(
                grass.parse_command(
                    "r.info",
                    map=self.test_output_all[0],
                    flags="g",
                )["ewres"],
            ),
        )
        self.assertTrue(
            ((out_res_ns == 1.0) and (out_res_ew == 1.0)),
            "Resolution of DOPs is not as expected."
            f"Expected res of 1, but got nsres of {out_res_ns} "
            f"and ewres of {out_res_ew}",
        )
        print("Test default settings successfully finished.\n")

    def extent_aoi_map_test(self):
        """
        Tests importing data only for AOI given by aoi_map
        """
        print(f"\nTest AOI ({self.fs}) ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            federal_state=self.fs.split(","),
            aoi=self.aoi_map,
        )
        self.assertModule(r_check, "Importing data for AOI fails.")
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el, msg=f"Creation of {el} failed.")
        # check if just aoi data loaded ==>
        # should have aoi_cells cells
        cells_aoi = grass.parse_command(
            "r.info",
            map=self.test_output_all[0],
            flags="g",
        )["cells"]
        self.assertTrue(
            int(cells_aoi) == self.aoi_cells,
            "Data loaded for larger region than AOI",
        )
        print(f"Test AOI ({self.fs}) successfully finished.\n")

    def dop_resolution_test(self):
        """
        Tests importing data with original resolution of DOPs
        """
        print(f"\nTest resolution ({self.fs}) ...")
        r_check = SimpleModule(
            "r.dop.import",
            output=self.test_output,
            aoi=self.aoi_map,
            federal_state=self.fs.split(","),
            flags="r",
        )
        self.assertModule(
            r_check,
            "Importing data with original DOP resolution fails."
            "Note: could also be a fail of aoi_map import "
            "(test_extent_aoi_map). Check this test first.",
        )
        for el in self.test_output_all:
            # test if all four raster maps are created
            self.assertRasterExists(name=el, msg=f"Creation of {el} failed.")
        out_res_ns = round(
            float(
                grass.parse_command(
                    "r.info",
                    map=self.test_output_all[0],
                    flags="g",
                )["nsres"],
            ),
            2,
        )
        out_res_ew = round(
            float(
                grass.parse_command(
                    "r.info",
                    map=self.test_output_all[0],
                    flags="g",
                )["ewres"],
            ),
            2,
        )
        # ref_res Auflösung und AOI
        self.assertTrue(
            ((out_res_ns == self.ref_res) and (out_res_ew == self.ref_res)),
            "Resolution of DOPs is not as expected."
            f"Expected res of {self.ref_res}, but got nsres of {out_res_ns}"
            f" and ewres of {out_res_ew}",
        )
        print(f"Test resolution ({self.fs}) successfully finished.\n")
