# remote to local
#rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/carbayes_within_prov_models/ /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/carbayes_within_prov_models/
#rsync -avr nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/carbayes_between_prov_models/ /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/carbayes_between_prov_models/
rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/meiotic_null_dist/_rslurm_meiotic_sib_permutations /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/meiotic_null_dist/_rslurm_meiotic_sib_permutations
rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/distance_likelihoods/_rslurm_distance_likelihoods_recursion/ /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/distance_likelihoods/_rslurm_distance_likelihoods_recursion/

# local to remote
rsync -av /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/data/distance_data
rsync -av /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/data/derived_data
rsync -av /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/map_bases/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/data/map_bases
