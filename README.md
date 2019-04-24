# Space_the_Final_Miptier

The repo and files at writing are a series of analysis scripts and some functions
(`R/`) that could be brought back into `MIPanalyzer`. Below is a description of 
the files and who wrote them so people know who to ask questions to etc if stuff 
breaks etc. 

---

## Analysis Scripts

1.  01-coi_adding_fitting.R - Loads COI data and fits a negative binomial to the spread of COI - OJ
2.  01-epi_covariate_import.R - Downloads DHS datasets for creating covariates - Nick
3.  01-genetic_import.R - Read in filtered MIP datasets for tidying and saving - Nick
4.  01-hotosm_CHiRPS_Covariate_import.R - Waterways, OSM and precipitation data collation - Nick
5.  01-spatial_data_import.R - Cleaning DRC province metas in dhs data - Nick
6.  02-data_wrangling_admins.R - Calculating admin weighted measures of DHS covariates for DRC - Nick          
7.  02-data_wrangling_cluster.R - Calculating cluster weighted measures of DHS covariates for DRC - Nick     
8.  03-bring_it_together.R - Script to call functions that calculate genetic distances and spatial distances and group them 
into a long format - OJ/Nick
9.  03-bring_sample_meta.R - Calculate a few possibly useful covariates at the cluster and admin level - OJ
10. 04-XY_plots.R - Script to call functions to plot genetic distances vs spatial distances at different spatial bins and groups 
and also some preliminary analysis looking at how prevalence impacts the decay in IBD relatedness against space - OJ

## R Files

00-base_functions.R - Functions used in the data wrangling file above - Nick
00-basic_map_functions.R - Similarly for getting values from maps - Nick
00-gendist_liftover_functions.R - Creates long formats for genetic distance. Somewhat deprecated by the functions in munging file - Nick  
00-distance_matrix_munging.R - Generic function for converting distance matrix into long format that calculates summary statistics for the distnaces and groups them at the selected level. Also functions that help automate calculating genetic distances - OJ
01-plotting.R- An overly fiddly plotting function that plots scatter plots with exponential smooths at the selected level of grouping - OJ

---

There is still lots to do but with what is here, the next steps could be to bring in the functions in 00-distance_matric_munging.R and 01-plotting.R and adpat them to be brought into `MIPanalyzer`. The fork of `MIPanalyzer` at `OJWatson/MIPanalyzer` also contains one extra distnace function and the normalisation for `get_genomic_distance`. Not sure how good we have been about calling `dplyr::` etc across all the functions. 

After that things to do could be:


1. Bring in the other genetic distnace written by Amy
1. Test this further by treating MOI as an unknown and jumble them to see if the distance measure changes at all
1. Do this jumble using simulated data
1. Test Bob's idea of creating a DRC genoytpe using the average PLAF for DRC, which we could use to test for highly related cutoffs 
1. Sensitivity analysis for whether clusters on admin borders shift the patterns observed if we change their admin allocation
1. Look at other spatial distance other than gc dist, such as the Enquette distances Keith was working on
1. Look at the other genetic distance patterns in terms of covariates and linear modelling (only briefly looked at the IBD)
1. Need to look at the inflection points of each measure and carry out power analysis to see how spatially spaced samples can be to detect the change is genetic distances
