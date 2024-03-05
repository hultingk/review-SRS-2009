# review-SRS-2009-Hulting

Contact: Katherine Hulting, <hultingk@msu.edu>

**File name: founder_plant_2009.csv**\
Each row contains data for one structure on one plant (3 structures/plant). 

| Variable          | Description                                                                                                                                                                             |
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| block   | Site name (8, 10, 52, 57, 53N, 53S, 54N, or 54S) |
| patch         | Patch Identifier (A = Center, B = Connected, C = Wing or Rectangle, D = Wing or Rectangle, E = Wing or Rectangle (labeled clock-wise from patch B)    |
| corner              | Corner of patch where plant is located. With back to patch A, I = top left, II = top right, III = bottom left, IV = bottom right     |
| distance              | Categorical distance of plant from edge. a = closest to edge, d = furthest from edge   |
| species              | Genus of focal plant species  |
| stalk         | stalk # of observation (1-3)            |
| plant\_ID        | Unique ID of plant, species.block.patch.corner.distance   |
| ptype    | Patch type (Center, Connected, Rectangle, Winged)      |
| live  | 0 = dead, 1 = alive, NA = not recorded |
| reproductive    | 0 = not reproductive, 1 = reproductive, NA = not recorded or dead plant  |
| no.viable_seeds | Number of fruits that contained a developed seed, NA = non-reproductive plant or not recorded |
| no.nonviable_seeds | Number of fruits that did not contain a developed seed, NA = non-reproductive plant or not recorded |
| notes | Notes on seed data, NA = no notes |
| dispersed_percent | Percent of structure that dispersed prior to collection, NA = 0% dispersal |
| no_basal_fl | Number of flowering stalks originating from the plant base, NA = non-reproductive plant or not recorded |
| no_axil_fl | Number of axilary flowering stalks, NA = non-reproductive plant or species with no axilary flowering stalks (Aristida, Anthaenantia, Sorghastrum) |
| predisp_seedpred | Percent of non viable seeds with signs of pre-dispersal seed predation, NA = non-reproductive plant or not recorded |
| flags | Data is flagged for issues. 0 = no issues, 1 = delete pollination or pre-dispersal seed predation data for this plant because it never flowered, 2 = pollination data not recorded, 3 = pollination and seed predation data not recorded, 4 = pollination data recorded but not being used, 5 = number of flowers not recorded, 6 = missing all data |
| Flag.Notes | Description of flag, NA = no flag or no further notes |
| dist_num | Distance of plant from edge in meters |
