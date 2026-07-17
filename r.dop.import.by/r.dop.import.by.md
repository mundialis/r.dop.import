<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.by* downloads and imports [digital orthophotos (DOP)](https://geodaten.bayern.de/opengeodata/) for Bayern (BY) and area of interest using the respective [WMS](https://geodaten.bayern.de/opengeodata/OpenDataDetail.html?pn=dop20rgb&active=SERVICE).
The data can be used when referencing the source:
id: CC BY 4.0,
name: Namensnennung 4.0 International,
url: [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/),
source: (c) Bayerische Vermessungsverwaltung ([Bayerische
Vermessungsverwaltung](https://geodaten.bayern.de/opengeodata/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.by aoi=aoi_BY output=dop_BY -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
