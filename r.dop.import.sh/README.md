<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.sh* downloads and imports [digital orthophotos (DOP)](https://geodaten.schleswig-holstein.de/gaialight-sh/_apps/dladownload/dl-dop20.html) for Schleswig-Holstein (SH) and area of interest.
The data can be used when referencing the source:
id: CC BY 4.0,
name: Namensnennung 4.0 International,
url: [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/),
source: (c) Landesamt für Vermessung und Geoinformation Schleswig-Holstein
([LVermGeo SH](https://www.schleswig-holstein.de/DE/landesregierung/ministerien-behoerden/LVERMGEOSH/Service/serviceGeobasisdaten/geodatenService_Geobasisdaten_DOP_digital))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.sh aoi=aoi_SH output=dop_SH -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
