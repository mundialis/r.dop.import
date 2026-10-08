<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.sl* downloads and imports [digital orthophotos (DOP)](https://geoportal.saarland.de/app-article/geobasisdatenuebersicht/) for Saarland (SL) and area of interest using the respective [WMS](https://geoportal.saarland.de/freewms/dop2025?).
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland - Namensnennung - Version 2.0,
url: [https://www.govdata.de/dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0),
source: (c) GeoBasis DE/LVGL-SL (2026) ([Landesamt für Vermessung, Geoinformation und Landentwicklung Saarland](https://geoportal.saarland.de/app-article/geobasisdatenuebersicht/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.sl aoi=aoi_SL output=dop_SL -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
