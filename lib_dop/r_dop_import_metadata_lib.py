#!/usr/bin/env python3
#
############################################################################
#
# MODULE:      r_dop_import_lib
# AUTHOR(S):   Leon Louwarts
# PURPOSE:     Library for r.dop.import
# SPDX-FileCopyrightText: (c) 2026 by mundialis GmbH & Co. KG and the
#                             GRASS Development Team
# SPDX-License-Identifier: GPL-3.0-or-later.
#
############################################################################

import os
import pathlib

import grass.script as grass


def get_download_urls_and_names(
    metadata_tmpfile,
    keep_data,
    download_dir,
    out_fs,
    fs,
):
    """Read download URLs from tempfile with fallback to download dir and
    raster count.

    Args:
        metadata_tmpfile (str): Path to tempfile written by state-specific
            addon
        keep_data (bool): True if downloaded data is kept in download_dir
        download_dir (str): Path to download directory
        out_fs (str): Prefix of output raster names in current mapset
        fs (str): Federal state abbreviation

    Returns:
        tuple[list, list]: (dop_urls, dop_names)
            dop_urls: list of download URLs (may be empty)
            dop_names: list of filenames or tile count label (may be empty)

    """

    dop_urls = []
    dop_names = []

    # Primary source: tempfile written by the state-specific addon, containing
    # one or more comma-separated URLs per tile (one tile per line)
    if metadata_tmpfile and pathlib.Path(metadata_tmpfile).exists():
        try:
            tile_url_groups = []
            with pathlib.Path(metadata_tmpfile).open(
                "r",
                encoding="utf-8",
            ) as f:
                for line in f:
                    url = line.strip()
                    if url:
                        tile_url_groups.append(url.split(","))

            dop_urls = [group[0] for group in tile_url_groups]
            grass.debug(
                f"Loaded {len(tile_url_groups)} tile groups "
                f"({sum(len(g) for g in tile_url_groups)} "
                "total URLs) from tempfile",
            )
            pathlib.Path(metadata_tmpfile).unlink()
        except Exception as e:
            grass.warning(
                f"Could not read tempfile {metadata_tmpfile} : {e}",
            )
            dop_urls = []

    # Fallback 1: tempfile was missing/empty -> scan the download directory
    # for locally kept files (only possible if -k flag was set)
    if not dop_urls and (
        keep_data and download_dir and pathlib.Path(download_dir).exists()
    ):
        for _root, _dirs, files in os.walk(download_dir):
            dop_names.extend(
                file
                for file in files
                if file.lower().endswith(
                    (".tif", ".tiff", ".jp2", ".jpeg"),
                )
            )

    # Fallback 2: neither URLs nor local files found -> count the imported
    # raster bands already present in the mapset to at least report a number
    # of tiles (used e.g. for WMS-based states where no URLs/files exist)
    if not dop_urls and not dop_names:
        all_rasters = []

        for mapset, rasters in grass.list_grouped(
            "raster",
        ).items():
            matching = [r for r in rasters if r.startswith(out_fs)]
            all_rasters.extend(matching)
            if matching:
                grass.debug(
                    f"Found {len(matching)} rasters in mapset {mapset}",
                )

        # Each DOP tile consists of 4 bands (red, green, blue, nir)
        num_dop_tiles = len(all_rasters) // 4 if len(all_rasters) >= 4 else 0
        grass.debug(
            f"Total rasters: {len(all_rasters)}, DOP tiles: {num_dop_tiles}",
        )

        if num_dop_tiles > 0:
            if fs in ["BW", "BY", "HE", "TH"]:
                dop_names = [
                    (
                        f"{num_dop_tiles} DOP-Kachel(n) "
                        "(via WMS heruntergeladen)"
                    ),
                ]
            else:
                dop_names = [f"{num_dop_tiles} DOP-Kachel(n)"]

    return dop_urls, dop_names
