## DESCRIPTION

*r.dop.import.bw* downloads and imports [digital orthophotos
(DOP)](https://www.lgl-bw.de/Produkte/Open-Data/) for Baden-Württemberg
(BW) and area of interest using the respective
[WMS](https://metadaten.geoportal-bw.de/geonetwork/srv/ger/catalog.search#/metadata/63985d3a-2be1-25bd-4cf8-ed84a1c7f7b1)
and
[WMS](https://metadaten.geoportal-bw.de/geonetwork/srv/ger/catalog.search#/metadata/4927a198-e9ce-da9a-a5ef-c11f21ff173f).  
The data can be used when referencing the source:  
id: dl-by-de/2.0,  
name: Datenlizenz Deutschland Namensnennung 2.0,  
url: https://www.govdata.de/dl-de/by-2-0,  
source: (c) LGL-BW
([LGL-BW](https://www.lgl-bw.de/Produkte//Luftbildprodukte/DOP20/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.bw aoi=aoi_BW output=dop_BW -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)  
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
