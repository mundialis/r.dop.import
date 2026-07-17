<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.he* downloads and imports [digital orthophotos (DOP)](https://hvbg.hessen.de/landesvermessung/geotopographie/luftbilder/digitale-orthophotos-true-orthophoto) for Hessen (HE) and area of interest using the respective [WMS](https://www.geoportal.hessen.de/mapbender/php/mod_showMetadata.php?resource=layer&layout=tabs&redirectToMetadataUrl=1&id=52119).
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland - Zero - Version 2.0,
url: [https://www.govdata.de/dl-de/zero-2-0](https://www.govdata.de/dl-de/zero-2-0),
source: (c) Hessische Verwaltung für Bodenmanagement und Geoinformation
([HVBG](https://hvbg.hessen.de/landesvermessung/geotopographie/luftbilder/digitale-orthophotos-true-orthophoto))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.he aoi=aoi_HE output=dop_HE -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
