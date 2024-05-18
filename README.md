# Firearms and Mental Health Providers

Analysis of the relationship between deaths from firearms, rates of gun ownership, and the ratio of mental health providers to the local population.  This repository contains R scripts for the statistical analyses and the generation of pertinent tables and figures.  Code is in pre-publication finalization status as of 5/17/2024

## Locations

.../R - Scripts

.../data - raw and cleaned data

.../results - tables, charts, figures, imputed data, etc

## Code

Script 1 - import and clean data
- Bins states into quartiles by gun ownership rate
- Generate contents of eTable 1
- Saves data as df.Rds

Script 2 - describe data
- Examine missingness
- Generate Table 1 and eTable 2
- Compare censored vs. uncensored subsets

Script 3 - imputation and modeling
- Pearson's correlations
- Figures 1 and 2
- Perform crude (unadjusted) and intermediate (limited adjustments) models
- Impute values for missing segregation data
- Apply a mixed effects model to the imputed data sets, and pool the coefficients into a final, complete-adjusted model.
- Generate Table 2 (concise and long form versions)

Script 4 - sensitivity analyses
- Compare random forest to primary mean matching, for imputation of segregation missing values
- Generate, visualize, and apply multiple sensitivity models for missing firearms numbers in censored counties

## Tips

To quickly view the model results:
```
load("results/imputation.RData")
summary(mc4)
summary(mi.continuous)
summary(mi.qrt)
summary(pool_pmm) # full model, vs continuous values for firearm fatalities
summary(pool_pmm_qrt) # vs quartiled bins of firearm fatalities
```

To generate sensitivity model coefficients only, run the following parts of script 4
- Part 1
- Part 5, lines to create `sensitivity.list`
- Part 6 (helper functions)
- Part 7 (batch run of models)
