## DESCRIPTION

*r.dop.import.ni* downloads and imports [digital orthophotos
(DOP)](https://ni-lgln-opengeodata.hub.arcgis.com/apps/lgln-opengeodata::digitales-orthophoto-dop/about)
for Niedersachsen (NI) and area of interest.  
The data can be used when referencing the source:  
id: CC BY 4.0,  
name: Namensnennung 4.0 International,  
url: https://creativecommons.org/licenses/by/4.0/,  
source: (c) LGL-Niedersachsen
([LGLN](https://ni-lgln-opengeodata.hub.arcgis.com/apps/lgln-opengeodata::digitales-orthophoto-dop/about))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.ni aoi=aoi_NI output=dop_NI -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)  
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
