## DESCRIPTION

*r.dop.import.sn* downloads and imports [digital orthophotos
(DOP)](https://www.geodaten.sachsen.de/batch-download-4719.html) for
Sachsen (SN) and area of interest.  
The data can be used when referencing the source:  
id: dl-by-de/2.0,  
name: Datenlizenz Deutschland - Zero - Version 2.0,  
url: https://www.govdata.de/dl-de/zero-2-0,  
source: (c) Landesamt für Geobasisinformation Sachsen
([GeoSN](https://www.geodaten.sachsen.de/batch-download-4719.html))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.nw aoi=aoi_SN output=dop_SN -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)  
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
