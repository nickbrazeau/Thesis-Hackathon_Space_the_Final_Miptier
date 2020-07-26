# Spatial patterns of Genetic Relatedness in the Democratic Republic of the Congo Driven by Human Movement

The structure of this research compendium loosely follows the framework of an "R-project" described by the [ROpenSci Team](https://github.com/ropensci/rrrpkg) for reproducibility.

An overview of the compendium is below: 
```
├── R                                             # "Helper" R functions  
│   ├── [a-zA-Z0-9_.-].R
├── README.md
├── Space_the_Final_Miptier.Rproj
├── analyses
│   ├── centrality_topology_analysees
│   ├── cluster_inbreeding
│   ├── highly_related_pairs
├── data
│   ├── derived_data
│   ├── distance_data
│   ├── map_bases
│   └── raw_data
├── import_wrangle
│   ├── 00-import_DHS.R
│   ├── 00-import_maps.R
│   ├── 00-scrape_drc_cities.R
│   ├── 00-smpl_mtdt.R
│   ├── 01-MAP_rasters.R
│   ├── 01-import_worldpop.R
│   ├── 01-nightlights_import.R
│   ├── 02-wrangle_urbanicity.R
│   ├── 03-Genetic_Data_Tidy.R
│   ├── 04-import_GreaterCircle.R
│   ├── 05-import_water_distances_v2.R
│   ├── 06-import_road_distances.R
│   ├── 07-import_flight_distances.R
│   └── 08-come_together.R
├── ind_reports
│   ├── 01-Sample_Clsts_Descriptive.Rmd
│   ├── 02-Genetic_Summaries.Rmd
│   ├── 03-IBD_S_distance_decay.Rmd
│   ├── 04-cluster_inbreeding.Rmd
│   ├── 05-network_analyses.Rmd
│   ├── 06-IBD_meioticsibs.Rmd
│   ├── 07-meiotic_permutations.Rmd
│   ├── 999-Prov_IBD.Rmd
│   ├── 999-SNP_autocorr.Rmd
│   └── 999-within_household.Rmd
├── one_ring_to_rull_them_all.Rmd
├── rdhs.json                                         # User-specific API Key to DHS (local not remote)
└── results
    ├── clust_inbd_results                            # Cluster Inbreeding Gradient Descent Results
    ├── figures                                       # Figures for Manuscript (and extras)
    ├── final_clstb_maps
    └── min_cost_inbreedingresult
```
