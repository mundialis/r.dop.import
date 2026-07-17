<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.bb.be* downloads and imports [digital orthophotos (DOP)](https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/) for Brandenburg and area of interest.
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland Namensnennung 2.0,
url: [https://www.govdata.de/dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0),
source: (c) LGB
([LGB](https://data.geobasis-bb.de/geobasis/daten/dop/rgbi_tif/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.bb aoi=aoi_BB output=dop_BB -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
