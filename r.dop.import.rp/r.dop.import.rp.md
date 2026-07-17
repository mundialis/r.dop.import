## DESCRIPTION

*r.dop.import.rp* downloads and imports [digital orthophotos
(DOP)](https://geobasis-rlp.de/data/dop20rgbi/current) for
Rhineland-Palatine (RP) and area of interest.
The data can be used when referencing the source:
id: dl-by-de/2.0,
name: Datenlizenz Deutschland - Zero - Version 2.0,
url: [https://www.govdata.de/dl-de/zero-2-0](https://www.govdata.de/dl-de/zero-2-0),
source: (c) Landesamt für Vermessung und Geobasisinformation
Rheinland-Pfalz
([LVermGeo](https://geobasis-rlp.de/data//dop20rgbi/current/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.rp aoi=aoi output=dop_RP -r
```

## AUTHORS

Victoria-Leandra Brunn, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)
