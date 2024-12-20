# CLL-association-study

## The followig files and scripts have been used to generate the figures related to association analyises in the paper _"A functional prognostic model predicts progression free survival for ibrutinib + venetoclax treated relapsed CLL"_ and are shared to be available for use on other cohorts  

**Associations_DSSs** - R script and input files to perform Wilcoxon Rank Sum tests on drug sensitivities related to binarized mutations and files needed to run the analysis. 

**Associations_pFlow** - R script and input files to perform Wilcoxon Rank Sum tests on phospho flow data related to binarized mutations and files needed to run the analysis.

**DSSs_standardization** - R script for standardizing DSS scores and file containing all combined raw DSS scores for the cohort. 

**Format_clinical_data** - R scripot for formatting clinical data files and all raw files used to make the Mut df 

**Random_forest_complete_data**- R script for performing random forest on single drugs and drug combinations with file used for the analysis 

**Random_forest_test_vs_oob_errors** R scripot to study overfitting in the random forest model using the Out of Bag approach versus a 10-fold testing approach. 

**SHAP_analysis** - Python script used to perform the SHAP analysis and the random forest files (generated in Random_forest_complete_data) to do the analysis

**Permutation_tests_rf** R script to perform permutation tests for the random forest analysis with required files

**pflow_standardization** R script for standardizing and doing noise corrections on phospho flow data with a file containing the combined MFI data for all patients in the cohort. 


