## DESCRIPTION

*r.dop.import.nw* downloads and imports [digital orthophotos
(DOP)](https://www.opengeodata.nrw.de/produkte/geobasis/lusat/akt/dop/dop_jp2_f10/)
for Nordrhein-Westfalen (NW) and area of interest.  
The data can be used when referencing the source:  
id: dl-by-de/2.0,  
name: Datenlizenz Deutschland - Zero - Version 2.0,  
url: https://www.govdata.de/dl-de/zero-2-0,  
source: (c) Landesbetrieb Information und Technik Nordrhein-Westfalen
([IT.NRW](https://www.opengeodata.nrw.de/produkte/geobasis/lusat/akt/dop/dop_jp2_f10/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.nw aoi=aoi_NW output=dop_NW -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)  
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
