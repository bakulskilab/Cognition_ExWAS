# Exposome-wide association study of cognition among older adults in the National Health and Nutrition Examination Survey

This repository contains R code files to run an exposome-wide association study (ExWAS) of cognition using public National Health and Nutrition Examination Survey (NHANES) data.

Processed NHANES datasets used in this analysis are available at: https://www.kaggle.com/datasets/nguyenvy/nhanes-19882018

 - Files needed: w - nhanes_1988_2018.RData, dictionary_nhanes.csv, chemicals_clean.csv, comments_clean.csv, weights_clean.csv

Original NHANES datasets are available at: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx

## Citation

Middleton LYM, Walker E, Cockell S, et al. Exposome-wide association study of cognition among older adults in the National Health and Nutrition Examination Survey. Exposome. 2025;5(1):osaf002. [doi:10.1093/exposome/osaf002](https://doi.org/10.1093/exposome/osaf002)

## File descriptions

1_Build_and_clean_dataset - data cleaning, participant inclusion and exclusion

2_Multiple_imputation_cognition_ExWAS - perform multiple imputation by chained equations

3_Cognition_ExWAS_descriptive_tables - create descriptive tables

4_Cognition_ExWAS_models - run generalized linear models

5_Cognition_ExWAS_volcano_plots - create volcano plots of results

6_Chemical_exposure_figures - create descriptive figures

Cog_ExWAS_model_functions - functions necessary to run models in code file 4
