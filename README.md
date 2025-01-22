This is the supplementary data for the Sparse Bayesian Matrix Clustering article. Here is a brief description of the files.


**Domain_Shift_Simuilation.R**

The R file used to generate the simulations seen in in section 4.3 of the article. The results are stored to an .rds file in the local directory.

**MFM_covariance_edited.R**

The file containing the function for Sparse Bayesian Matrix Clustering. The function used is named MFM_MxN_equal_cov_sparse, with hyperparameters discussed here https://arxiv.org/abs/2010.08495.

**U_covariance_update.R**

The file containing the U_update function, used for the covariance update in SBMC algorithm, with details on hyperparameters discussed here: https://www.sciencedirect.com/science/article/abs/pii/S0047259X22000744

**Dauls_Method.R**

The file containing the posterior inference function used on the output of the MFM_MxN_equal_cov_sparse function cluster assignments. 

**austin_energy_array_cleaned.rda**

RDA file contianing the raw data used for energy clustering, without log(1+x) transform. All months are included in this set. 

**austin_house_index_cleaned.rda**

The RDA file distinguishing unique houses matching the index of the austin_energy_array_cleaned.rda file. 
