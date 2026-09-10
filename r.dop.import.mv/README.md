<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import.mv* downloads and imports [digital orthophotos (DOP)](https://laiv.geodaten-mv.de/afgvk/Luftbilder/Beschreibung?produkt=DOP) for Mecklenburg-Vorpommern (MV) and area of interest.
The data can be used when referencing the source:
id: CC BY 4.0,
name: Namensnennung 4.0 International,
url: [https://creativecommons.org/licenses/by/4.0/](https://creativecommons.org/licenses/by/4.0/),
source: (c) GeoBasis-DE/M-V
([LAIV-MV](https://www.laiv-mv.de/Geoinformation/))

## EXAMPLES

### Import DOPs

Import DOPs with native resolution:

```sh
r.dop.import.mv aoi=aoi_MV output=dop_MV -r
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Leon Louwarts, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
