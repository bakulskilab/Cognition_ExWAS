

############################## MAIN ANALYSIS ############################## 
run_exwas_model <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(include == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants
  svydata_subset <- subset(svydata,
                           include == 1)
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(DSST ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIAGENDR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}





############################## NORMAL EGFR SUBSET ############################## 
run_exwas_model_normal_egfr <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(egfr_grp == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants
  svydata_subset <- subset(svydata,
                           egfr_grp == 1)
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(DSST ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIAGENDR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}





############################## ADD LIFESTYLE COVARIATES ############################## 
run_exwas_model_lifestyle <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(include == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants
  svydata_subset <- subset(svydata,
                           include == 1)
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(DSST ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIAGENDR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle +
                                             BMXWAIST +
                                             alc_month_freq))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle +
                                               BMXWAIST +
                                               alc_month_freq))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               BMXWAIST +
                                               alc_month_freq))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle +
                                               BMXWAIST +
                                               alc_month_freq))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               BMXWAIST +
                                               alc_month_freq))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle +
                                               BMXWAIST +
                                               alc_month_freq))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               BMXWAIST +
                                               alc_month_freq))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}





############################## POISSON MODEL ############################## 
run_exwas_model_poisson <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(include == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants
  svydata_subset <- subset(svydata,
                           include == 1)
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIAGENDR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle,
                                           family=quasipoisson()))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(exp_est = exp(estimate),
           ci_lower_exp = exp(estimate - (1.96*std.error)),
           ci_upper_exp = exp(estimate + (1.96*std.error)),
           ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}





############################## POISSON MODEL - NORMAL EGFR ############################## 
run_exwas_model_poisson_egfr <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(egfr_grp == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants
  svydata_subset <- subset(svydata,
                           egfr_grp == 1)
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIAGENDR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle,
                                           family=quasipoisson()))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle,
                                             family=quasipoisson()))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(mild_cog_impair ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIAGENDR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp,
                                             family=quasipoisson()))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(exp_est = exp(estimate),
           ci_lower_exp = exp(estimate - (1.96*std.error)),
           ci_upper_exp = exp(estimate + (1.96*std.error)),
           ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}





############################## STRATIFY BY SEX ############################## 
run_exwas_model_males <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(include == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants and males only
  svydata_subset <- subset(svydata,
                           include == 1 & RIAGENDR == 'Male')
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(DSST ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}




run_exwas_model_females <- function(chem_subset, weights, ewas_data)
{
  # pull chemical variable name
  chem_name <- unique(chem_subset$chem_name)
  wt_name <- paste('WT_', chem_name, sep = '')
  
  print(paste('Chemical variable:', chem_name))
  print(paste('Weight variable:', wt_name))
  
  # pull weight variable from weights dataset
  wt_var <- weights %>%
    select(SEQN, all_of(wt_name))
  
  # give the weight var a generic name
  names(wt_var) <- c('SEQN', 'wt_unadj')
  
  # merge chem and weight
  chem_and_wt <- left_join(chem_subset,
                           wt_var,
                           join_by(SEQN)) %>%
    select(-chem_name)
  
  # merge chem_and_wt with ewas_data
  ewas_with_chem <- left_join(ewas_data,
                              chem_and_wt,
                              by = c('SEQN', 'SDDSRVYR'))
  
  # pull the number of cycles with the chem measure for included participants
  incl_only <- ewas_with_chem %>%
    filter(include == 1) %>%
    drop_na(chem_measure)
  
  list_cycles <- unique(incl_only$SDDSRVYR)
  num_cycles <- length(list_cycles)
  
  print(paste('Number of cycles:', num_cycles))
  
  # calculate adjusted weights - because no longer using cycle 2, all are adjusted the same way
  ewas_with_chem <- ewas_with_chem %>%
    mutate(wt_adj = (1/num_cycles)*wt_unadj)
  
  # extract each imputation as a separate dataset and create an imputationList object
  imp1 <- ewas_with_chem %>%
    filter(.imp == 1) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp2 <- ewas_with_chem %>%
    filter(.imp == 2) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp3 <- ewas_with_chem %>%
    filter(.imp == 3) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp4 <- ewas_with_chem %>%
    filter(.imp == 4) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp5 <- ewas_with_chem %>%
    filter(.imp == 5) %>%
    drop_na(wt_adj) %>%
    filter(wt_adj != 0)
  
  imp_list <- imputationList(list(imp1, imp2, imp3, imp4, imp5))
  
  # set up the survey design using the imputationList
  options(survey.lonely.psu = 'adjust')
  
  svydata <- svydesign(ids = ~SDMVPSU, 
                       strata = ~SDMVSTRA,
                       weights = ~wt_adj,
                       nest = TRUE,
                       data = imp_list)
  
  # subset to included participants and females only
  svydata_subset <- subset(svydata,
                           include == 1 & RIAGENDR == 'Female')
  
  # run the linear regression model for each chemical exposure
  # if/else statements control inclusion of cotinine, creatinine, and cycle as covariates:
  # creatinine only included for urinary measurements
  # cotinine excluded for smoking related chems
  # cycle included if num_cycles >1
  
  if(str_detect(chem_name, '^LB')) { # blood measures
    
    print(paste(chem_name, '= blood measure'))
    if(chem_name == 'LBXCOT') { # blood cotinine
      
      print('Blood cotinine model')
      model <- with(svydata_subset, svyglm(DSST ~
                                             scale(log2(chem_measure)) +
                                             RIDAGEYR +
                                             RIDRETH1 +
                                             hs_educ +
                                             smk_status +
                                             total_seafood_grp +
                                             cycle))
      
    } else { # not blood cotinine
      
      if(num_cycles > 1) {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Blood measure, not cotinine: LBXCOT kept in covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) + # use imputed version as covariate
                                               smk_status +
                                               total_seafood_grp))
        
      }
      
    } 
    
  } else { # urine measures
    
    print(paste(chem_name, '= urine measure'))
    if(chem_name %in% c('URXNAL', 'URXCOTT', 'URXHCTT', 'URXSCN')) { # smoking vars
      
      if(num_cycles > 1) {
        
        print('Smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
    } else { # non smoking vars
      
      if(num_cycles > 1) {
        
        print('Not a smoking chem')
        print('Measured in >1 cycles, cycle added to covariates')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp +
                                               cycle))
        
      } else {
        
        print('Not a smoking chem')
        model <- with(svydata_subset, svyglm(DSST ~
                                               scale(log2(chem_measure)) +
                                               RIDAGEYR +
                                               RIDRETH1 +
                                               hs_educ +
                                               scale(log2(LBXCOT_imp)) +
                                               smk_status +
                                               scale(log2(URXUCR)) +
                                               total_seafood_grp))
        
      }
      
      
    }
    
  }
  
  # pool the results and calculate CI
  model_pool <- summary(pool(model,
                             dfcom = model[[1]][['degf.resid']])) %>%
    mutate(ci_lower = (estimate - (1.96*std.error)),
           ci_upper = (estimate + (1.96*std.error)))
  
  model_pool <- data.frame(model_pool)
  
  return(model_pool)
  
}




