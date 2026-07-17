<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.th* downloads and imports [digital orthophotos (DOP)](https://geoportal.thueringen.de/gdi-th/download-offene-geodaten/download-luftbilder-und-orthophotos) for Thüringen (TH) and area of interest using the respective
[WMS](https://geoportal.thueringen.de/gdi-th/download-offene-geodaten/darstellungs-und-downloaddienste).
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland Namensnennung 2.0,
url: [https://www.govdata.de/dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0),
source: (c) GDI-Th, Freistaat Thueringen
([GDI-Th](https://geoportal.thueringen.de/gdi-th/download-offene-geodaten/download-luftbilder-und-orthophotos))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.th aoi=aoi_TH output=dop_TH -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
