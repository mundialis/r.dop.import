# r.dop.import - Toolset for the import of digital orthophotos (DOPs)

It includes import addons for the open geodata digital orthophotos for Germany,
e.g. for the digital orthophotos (DOPs) of Brandenburg.

The r.dop.import toolset consists of the following modules:

-  r.dop.import: downloads digital orthophotos (DOPs) 
  for specified federal state and area of interest,
  and creates a single file of all downloaded digital orthophotos (DOPs).
-  r.dop.import.*: prepares the parallel download of digital orthophotos (DOPs) 
  for the respective federal state in the selected AOI
-  r.dop.import.worker.\*: downloads a single digital orthophoto (DOP) 
  passed by r.dop.import.\*
