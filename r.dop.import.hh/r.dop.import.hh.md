## DESCRIPTION

*r.dop.import.hh* downloads and imports [digital orthophotos
(DOP)](https://suche.transparenz.hamburg.de/dataset/luftbilder-hamburg-dop-zeitreihe-unbelaubt3)
for Hamburg (HH) and area of interest.
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland Namensnennung 2.0,
url: [https://www.govdata.de/dl-de/by-2-0](https://www.govdata.de/dl-de/by-2-0),
source: (c) LGV-HH
([LGV-HH](https://suche.transparenz.hamburg.de/dataset/luftbilder-hamburg-dop-zeitreihe-unbelaubt3))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.hh aoi=aoi_HH output=dop_HH -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
