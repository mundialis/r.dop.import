<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.st* downloads and imports [digital orthophotos (DOP)](https://www.lvermgeo.sachsen-anhalt.de/de/gdp-open-data.html) for Sachsen-Anhalt (ST) and area of interest.
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland - Namensnennung - Version 2.0,
url: [https://www.govdata.de/dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0),
source: (c) GeoBasis-DE/LVermGeo ST
([LVermGeo ST](https://www.lvermgeo.sachsen-anhalt.de/de/gdp-open-data.html))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.st aoi=aoi_ST output=dop_ST -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
