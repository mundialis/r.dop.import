<!-- markdownlint-disable MD041 -->
## DESCRIPTION

*r.dop.import* downloads digital orthophotos (DOPs) for specified federal state and area of interest, stores it in a local directory and creates a single file of DOPs in GRASS. Implemented federal state options are:

- Baden-Württemberg (BW)
- Bayern (BY)
- Berlin (BE)
- Brandenburg (BB)
- Hessen (HE)
- Niedersachsen (NI)
- Nordrhein-Westfalen (NW)
- Rheinland-Pfalz (RP)
- Sachsen (SN)
- Thüringen (TH)

## EXAMPLE

```sh
r.dop.import fs=Nordrhein-Westfalen aoi=Polygon_BonnBeuel output=NRW_DOP_output
```

## AUTHORS

Johannes Halbauer, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)
Anika Weinmann, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Lina Krisztian, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Julia Haas, [mundialis GmbH & Co. KG](https://www.mundialis.de/)
Victoria-Leandra Brunn, [mundialis GmbH & Co.
KG](https://www.mundialis.de/)
