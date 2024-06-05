# Environmental chemical-wide association study of cognition among older adults in the National Health and Nutrition Examination Survey

This repository contains R code files to run an environment-wide association study (ExWAS) of cognition using public National Health and Nutrition Examination Survey (NHANES) data.

Processed NHANES datasets used in this analysis are available at: https://www.kaggle.com/datasets/nguyenvy/nhanes-19882018

Original NHANES datasets are available at: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx

## Citation

## Abstract

## File descriptions

1_Build_and_clean_dataset.Rmd - data cleaning, participant inclusion and exclusion

2_Multiple_imputation_cognition_ExWAS.Rmd - perform multiple imputation by chained equations

3_Cognition_ExWAS_descriptive_tables_and_figures.Rmd - create descriptive tables

4_Chemical_supplemental_figures.Rmd - create descriptive figures

5_Cognition_ExWAS_models.Rmd - run generalized linear models

6_Cognition_ExWAS_volcano_plots.Rmd - create volcano plots of results

Cog_ExWAS_model_functions.R - functions necessary to run models in code 5
