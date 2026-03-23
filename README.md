# r.dop.import - Toolset for the import of digital orthophotos (DOPs)

It includes import addons for the open geodata digital orthophotos for Germany,
e.g. for the digital orthophotos (DOPs) of Brandenburg.

The r.dop.import toolset consists of the following modules:

- r.dop.import: downloads digital orthophotos (DOPs)
  for specified federal state and area of interest,
  and creates a single file of all downloaded digital orthophotos (DOPs).
- r.dop.import.*: prepares the parallel download of digital orthophotos (DOPs)
  for the respective federal state in the selected AOI
- r.dop.import.worker.\*: downloads a single digital orthophoto (DOP)
  passed by r.dop.import.\*

## Addon coverage for federal states

| Federal state | DOP Import Addon | Tile-index | Data Download | Resolution | Data Source |
| - | - | - | - | - | - |
| BB | &#9745; | &#9745; | data download as .zip and .tif | 20cm | |
| BE | &#9745; | &#9745; | data download as .zip and .tif | 20cm | |
| BW | | | | | [Open GeoData Portal](https://opengeodata.lgl-bw.de/#/) |
| BY | | | | | [Open GeoData](https://geodaten.bayern.de/opengeodata/) |
| HB | | | | | [Geoportal](https://geoportal.bremen.de/geoportal/#) |
| HE | (&#9745;) to be fixed | | | | |
| HH | | | | | [Geoportal](https://geoportal-hamburg.de/) |
| MV | | | | | [Downloadportal](https://laiv.geodaten-mv.de/afgvk/) |
| NI | &#9745; | &#9745; | data download as .tif | 20cm | |
| NW | &#9745; | &#9745; | data download as .jp2 | 10cm | |
| RP | &#9745; | &#9745; (no automatic updating) | | 20cm | |
| SH | | | | | [Downloadportal](https://geodaten.schleswig-holstein.de/gaialight-sh/_apps/dladownload/) |
| SL | | | | | [Geoportal](https://geoportal.saarland.de/) |
| SN | &#9745; | &#9745; | data download as .zip and .tif plus .csv | 20cm | |
| ST | | | | | [Geodatenportal](https://www.lvermgeo.sachsen-anhalt.de/de/gdp-open-data.html) |
| TH | (&#9745;) to be fixed | | | | |
