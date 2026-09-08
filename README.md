<!-- markdownlint-disable MD041 -->

[![image-alt](https://github.com/OSGeo/grass/raw/main/man/grass_logo.png)](https://grass.osgeo.org/grass-stable/manuals/index.html)

______________________________________________________________________

## NAME

***r.dop.import*** - Toolset for the import of digital orthophotos (DOPs). It includes import addons for the open geodata digital orthophotos for Germany and downloads digital orthophotos for a specified area of interest (AOI) and multiple federal states of Germany.

## KEYWORDS

[raster](https://grass.osgeo.org/grass-stable/manuals/keywords.html#raster), [import](https://grass.osgeo.org/grass-stable/manuals/keywords.html#import), [elevation](https://grass.osgeo.org/grass-stable/manuals/keywords.html#elevation)

## DESCRIPTION

### Modules in this toolset

- [r.dop.import](r.dop.import/README.md): downloads digital orthophotos (DOPs) for specified federal state and area of interest, and creates a single file of all downloaded digital orthophotos (DOPs).
- r.dop.import.{fs}: prepares the parallel download of digital orthophotos (DOPs) for the respective federal state in the selected AOI
- r.dop.import.worker.{fs}: downloads a single digital orthophoto (DOP) passed by r.dop.import.{fs}

### Overview of the available digital orthophotos

| Federal state | fs | DOP Import Addon | Tile-index | Data Download | Resolution | Data Source |
| - | - | - | - | - | - | - |
| Baden-Württemberg | BW | &#9745; | (via WMS) | | 20cm | [Open GeoData Portal](https://opengeodata.lgl-bw.de/#/) |
| Bayern | BY | &#9745; | (via WMS) | | 20cm | [Open GeoData](https://geodaten.bayern.de/opengeodata/) |
| Berlin | BE | &#9745; | &#9745; | data download as .zip and .tif | 20cm | |
| Brandenburg | BB | &#9745; | &#9745; | data download as .zip and .tif | 20cm | |
| Bremen | HB | &#9745; | (via WMS) | | 10cm | [Geoportal](https://geoportal.bremen.de/geoportal/#) |
| Hamburg | HH | &#9745; | &#9745; | data download as .zip and .tif | 20cm | [Geoportal](https://geoportal-hamburg.de/) |
| Hessen | HE | &#9745; | (via WMS) | | 20cm | |
| Mecklenburg-Vorpommern | MV | | | | | [Downloadportal](https://laiv.geodaten-mv.de/afgvk/) |
| Niedersachsen | NI | &#9745; | &#9745; | data download as .tif | 20cm | |
| Nordrhein-Westfalen | NW | &#9745; | &#9745; | data download as .jp2 | 10cm | |
| Rheinland-Pfalz | RP | &#9745; | &#9745; (no automatic updating) | | 20cm | |
| Saarland | SL | | | | | [Geoportal](https://geoportal.saarland.de/) |
| Sachsen | SN | &#9745; | &#9745; | data download as .zip and .tif plus .csv | 20cm | |
| Sachsen-Anhalt | ST | | | | | [Geodatenportal](https://www.lvermgeo.sachsen-anhalt.de/de/gdp-open-data.html) |
| Schleswig-Holstein | SH | &#9745; | &#9745; | data download as .tif | 20cm | [Downloadportal](https://geodaten.schleswig-holstein.de/gaialight-sh/_apps/dladownload/) |
| Thüringen | TH | &#9745; | (via WMS) | | 20cm | |

## REQUIREMENTS

[grass-gis-helpers\>=4.0.0](https://pypi.org/project/grass-gis-helpers/)

## SEE ALSO

*[r.dem.import](https://github.com/mundialis/r.dem.import) for import of digital elevation models*

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)  
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)  
Julia Haas, [mundialis GmbH & Co. KG](https://www.mundialis.de/)  
Lina Krisztian, [mundialis GmbH & Co. KG](https://www.mundialis.de/)  
Victoria-Leandra Brunn, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
