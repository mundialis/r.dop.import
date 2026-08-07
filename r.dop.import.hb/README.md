<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.hb* downloads and imports [digital orthophotos (DOP)](https://www.geo.bremen.de/produkte/luftbildprodukte-11703) for Bremen and Bremerhaven (HB) and area of interest using the respective [WMS](https://geodienste.bremen.de/wms_dop_lb).
The data can be used when referencing the source:
id: CC-BY 4.0,
name: Creative Commons Namensnennung 4.0 International,
url: [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/),
source: (c) Landesamt GeoInformation Bremen
([Geo Bremen](https://www.geo.bremen.de/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.hb aoi=aoi_HB output=dop_HB -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
