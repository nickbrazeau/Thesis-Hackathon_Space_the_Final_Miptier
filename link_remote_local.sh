# remote to local
rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/carbayes_sp_dics/_rslurm_CARleroux_DICs_slope/ /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/carbayes_sp_dics/CARleroux_DICs_slope/
rsync -av nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/results/carbayes_sp_dics/_rslurm_CARleroux_DICs_withinIBD/ /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/results/carbayes_sp_dics/CARleroux_DICs_withinIBD/

# local to remote
rsync -av /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/distance_data/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/data/distance_data
rsync -av /Users/nickbrazeau/Documents/GitHub/Space_the_Final_Miptier/data/derived_data/ nfb@longleaf.unc.edu:/proj/ideel/meshnick/users/NickB/Projects/Space_the_Final_Miptier/data/derived_data
