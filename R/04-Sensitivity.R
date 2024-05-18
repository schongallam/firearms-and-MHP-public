# Firearms and mental health providers
# part 4: Supplemental imputation
# 2024-03-07
# Kelly Corbett and Malcolm Schongalla

# main project folder:                     .../firearms-and-MHP
# R scripts folder (assumed working dir):  .../firearms-and-MHP/R
# data folder:                             .../firearms-and-MHP/data

# Script index:
# 01 - (import data) load data from external sources
# 02 - (describe data) inspect, compare, create Table 1 & eTable 3
# 03 - (impute & model) correlations, linear models, and [e]Tables 2
# 04 - (this script) - compares imputation techniques and sensitivity models

# Script 4 sequence:
# 1. Libraries & Load Data
# 2. Random Forest imputation
# 3. Compare pmm and rf imputed values
# 4. Compare pmm and rf models
# 5. Sensitivity Data Set Visualization
# 6. Functions For Processing Sensitivity Analyses
# 7. Sensitivity Modeling
# 8. Generate Forest Plots

########################### ###
# 1. Libraries & Load Data ####
########################### ###

# Pick one or the other:
set.seed(123) # for reproducibility in testing
# set.seed(Sys.time()) # for randomness

library(dplyr)
library(lmerTest) # Satterthwaite methods in lmer (not in lme4)
library(mice)
library(howManyImputations)
library(broom.mixed)
library(ggpubr) # ggscatter. depends ggplot2
library(forestplot)
requireNamespace("tibble") # rownames_to_column

load("results/imputation.RData")
checkfor <- c("df", "pred", "meth_pmm", "imp1_pmm", "fit_pmm", "fit_pmm_qrt",
              "pool_pmm", "pool_pmm_qrt", "model.mi.covariates",
              "formulaStr.mi.continuous", "formulaStr.mi.qrt",
              "model.full.covariates", "model.full.continuous",
              "model.full.qrt")
stopifnot(all(lapply(checkfor, exists) == TRUE))


############################## ###
# 2. Random Forest imputation ####
############################## ###

# TL;DR-
# - Random Forest is slower than PMM.
# - It imputes slightly higher values than PMM.
# - It is not likely to make any significant difference due to the low estimate
#   associated with segregation.
# - There is no reason to prefer this over PMM in any of the main models.

meth_rf <- meth_pmm
meth_rf["segregation"] <- "rf"

imp.ct <- 75
imp1_rf <- mice(df,
                m=imp.ct,
                maxit = 5,
                predictorMatrix = pred,
                method = meth_rf,
                print = TRUE )
stripplot(imp1_rf, segregation ~ .imp, pch=20, cex=2)
#pdf(file="results/mice_stripplot_rf.pdf", width=9, height=9, family="Helvetica")
#  stripplot(imp1_rf, segregation ~ .imp, pch=20, cex=2)
#dev.off()

plot(imp1_rf)
#pdf(file="results/mice_conv_rf.pdf", width=9, height=9, family="Helvetica")
#  plot(imp1_rf,
#       main=sprintf("mice convergence check (rf)\nSegregation Index(0-100), iterations=%d", imp.ct),
#       layout=c(1,2))
#dev.off()

densityplot(imp1_rf, ~segregation | .imp)
#pdf(file="results/mice_dp_rf.pdf", width=9, height=9, family="Helvetica")
#  densityplot(imp1_rf, ~segregation | .imp,
#            main="mice densityplot of imputations (rf)")
#dev.off()

fit_rf <- with(imp1_rf, lmer(as.formula(model.full.continuous)))
how_many_imputations(fit_rf) # 85 with m=20; 62 with m=50; 53 with m=75

# The following code is optional depending on the original imp.ct and
# result of how_many_imputations
imp2_rf <- mice(df,
                m=25, # 20+25 exceeds 39 but ensures total exceeds re-check
                maxit = 3,
                predictorMatrix = pred,
                method = meth_rf,
                print = TRUE )
imp2_rf <- ibind(imp1_rf, imp2_rf)
fit2_rf <- with(imp2_rf, lmer(as.formula(model.full.continuous)))
how_many_imputations(fit2_rf) # #35
#/optional


####################################### ###
# 3. Compare pmm and rf imputed values ####
####################################### ###

suppressMessages(comp_pmm <- complete(imp1_pmm, action="broad")) # ~1,200 lines
segr_pmm <- select(comp_pmm, matches("segregation"))
suppressMessages(comp_rf <- complete(imp1_rf, action="broad"))
segr_rf <- select(comp_rf, matches("segregation"))

segr_mean_pmm <- colMeans(segr_pmm, na.rm=TRUE)
segr_mean_rf <- colMeans(segr_rf, na.rm=TRUE)
t.test(segr_mean_pmm, segr_mean_rf, alternative="two.sided")
# Interpretation: rf imputes slightly higher segregation values than pmm.
# Probably not impactful.

############################### ###
# 4. Compare pmm and rf models ####
############################### ###

pool_rf <- pool(fit_rf)
# pool_rf <- pool(fit2_rf) # use this if you ran the optional extra imputations
summary(pool_pmm)
summary(pool_rf)
# Interpretation comparing pmm to rf: the segregation estimate in both cases
# is low (0.009, 0.010 respectively) and p value > 0.3, so it makes no major
# direct difference.  The impact on gunOwnership and MHProv estimates is
# also miniscule compared to the order of the overall gunOwnership estimate.


######################################## ###
# 5. Sensitivity Data Set Visualization ####
######################################## ###

# Descriptions from script 3
# lo:  fill censored with '1'
# mid: ... with '5'
# hi:  ... with '9'
# proportional: ... with 1:9 proportional to min:max gunOnwership
# invProportional: inverse formula of proportional
# invMHP: ... with 1:9 inversely proportional to min:max MHProv

# Scatter plot can show outliers
ggscatter.sensitivity <- function(X, censored.only=FALSE, ystr = NULL) {
  df %>%
    filter(censored == TRUE | censored == censored.only) %>%
    select(all_of(X), gunOwnership) %>%
    ggscatter(
      x = "gunOwnership", 
      y = X,
      add = "reg.line",
      xlab = "Proportion of Population Identified as Firearm Owners", 
      ylab = sprintf("Firearm Fatality %s\nfor %s", ystr, X),
      title = sprintf("Scatterplot and Pearson's correlation of Firearm Ownership\nand Firearm Fatality Sensitivity Vector; n = %d", nrow(.)),
      caption = "\nNote: Each observation represents a county’s firearm fatality rate, grouped in columns based on the state’s proportion of firearm ownership.\n The state firearm ownership rate was applied equally to each county within the state in the absence of county-specific data"
    ) + stat_cor(method="pearson", p.accuracy = 0.001)
}

# Plots two related scatterplots one after another.
# Plotting them consecutively like this helps with contextual inspection
sidebyside.ggscatter <- function(X, Y) {
  ratestr = "Rate (per 100,000)"
  countstr = "count"
  print(sprintf("X = %s, Y = %s", X,Y))
  print(ggscatter.sensitivity(X, TRUE, countstr))
  print(ggscatter.sensitivity(Y, FALSE, ratestr))
}

# bar plot showing general raw distribution
barplot.sensitivity.raw <- function(X, censored.only=FALSE) {
  df %>%
    filter(censored == TRUE | censored == censored.only) %>%
    select(all_of(X)) %>%
    unlist %>%
    as.vector(mode="integer") %>%
    round() %>%
    factor(levels = 1:9) %>%
    table() %>%
    barplot(main = "Rounded imputations of censored firearms data",
            sub = X,
            xlab = "Rounded value imputed",
            ylab = "Count",
            ylim = range(pretty(c(0, max(.))))
    )
}

# Barplot of the rates of deaths by county instead of the raw numbers
barplot.sensitivity <- function(X, censored.only=FALSE) {
  df %>%
    filter(censored == TRUE | censored == censored.only) %>%
    select(all_of(X)) %>%
    unlist %>%
    as.vector(mode="integer") %>%
    round() %>%
    table() %>%
    barplot(main = "Distribution of firearms deaths per 100,000\nusing sensitivity vector",
            sub = X,
            xlab = "Rate of firearms deaths",
            ylab = "Number of counties",
            ylim = range(pretty(c(0, max(.))))
    )
}

# Similarly to above, aids contextual inspection
sidebyside.bar <- function(X) {
  barplot.sensitivity.raw(X, TRUE)
  barplot.sensitivity(X, FALSE)
}

# This list also used in parts 6 & 7
sensitivity.list <- list('firearms.lo',
                         'firearms.mid',
                         'firearms.hi',
                         'firearms.pop',
                         'firearms.proportional',
                         'firearms.invProportional',
                         'firearms.invMHP',
                         'firearms.MHP')
sensitivity.list.raw <- list('firearms.numerator.lo',
                             'firearms.numerator.mid',
                             'firearms.numerator.hi',
                             'firearms.numerator.pop',
                             'firearms.numerator.proportional',
                             'firearms.numerator.invProportional',
                             'firearms.numerator.invMHP',
                             'firearms.numerator.MHP')

# bar graphs of sensitivity data set distributions
lapply(sensitivity.list.raw, sidebyside.bar)

# scatter plots of sensitivity data set distributions
mapply(sidebyside.ggscatter, sensitivity.list.raw, sensitivity.list)


################################################### ###
# 6. Functions For Processing Sensitivity Analyses ####
################################################### ###

# The intermediate models and the full models are processed differently,
# because an intermediate model is a simple lmer and the full model is a
# product of multiple imputation and needs pooling of estimates.

# lmer was run on the unimputed data set with the full model formula, which
# produces a "complete case" analysis.  TL;DR- inspection of these estimates
# across various sensitivity analyses reveals the dramatic effect of ignoring
# the missingness of 'segregation'!

# Better function documentation is provided because it may be a little less
# obvious what is going on here, and some functions may be re-used later

make.formulas <- function(dependent.var.list,
                          primary.covar,
                          other.covars.str) {
  #' Generate formulas for sensitivity analyses
  #' 
  #' Makes a list of formulas for use with lmer in this script, for batching sensitivity models
  #' 
  #' @param dependent.var.list List of dependent variables of interest
  #' @param primary.covar The covariate of primary interest
  #' @param other.covars.str A string of the remaining covariates in the right-hand side of the desired formula
  #' @usage make.formulas(dependent.var.list, primary.covar, other.covars.str)
  #' @return List of dataframes: first column is the dependent.var.list. Second column is the formatted formula
  #' @note Caution, name conflicts with mice::make.formulas
  data.frame(dep.var = unlist(dependent.var.list)) %>%
    mutate(formulaStr = sprintf("%s ~ %s + %s", dep.var, primary.covar, other.covars.str))
}

get_coefs <- function(model.list) {
  #' Coefficients of a list of lmer objects
  #' 
  #' Takes a list of lmer objects, and returns a list of tables that contain coefficients from those models
  #' 
  #' @param model.list List of lmer objects
  #' @usage get_coefs(model.list)
  #' @return List of dataframes: first column contains dependent variables.  The remaining columns are from stats::coef
  lapply(model.list, \(x) x %>%
           summary %>%
           coef %>%
           as.data.frame %>%
           tibble::rownames_to_column("covar"))
}

get_CIs <- function(model.list) {
  #' Confidence intervals of a list of lmer objects
  #' 
  #' Takes a list of lmer objects, and returns a list of tables that contain 95% confidence intervals.
  #' 
  #' @param model.list List of lmer objects
  #' @usage get_CIs(model.list)
  #' @return List of dataframes: first column is the dependent variables.  The remaining columns are from stats::confint
  lapply(model.list, \(x) x %>%
           confint %>%
           as.data.frame %>%
           tibble::rownames_to_column("covar"))
}

merge_coef_CI_tables <- function(coefs, CIs) {
  #' Merge lists of coefficient tables and CI tables
  #' 
  #' Vectorized left_join on the products of get_coefs and get_CIs
  #' 
  #' @param coefs from get_coefs
  #' @usage CIs from get_CIs
  #' @return List of dataframes: left join of coefs and CIs tables
  mapply(left_join, coefs, CIs,
         MoreArgs = list(by=join_by(covar)),
         SIMPLIFY=FALSE)
}

get_model_results <- function(model_list) {
  #' Model parameters
  #' 
  #' Gets coefficitents and confidence intervals from a list of lmer objects, and returns a list of corresponding tables.
  #'
  #' @param model_list A list of lmer objects
  #' @usage get_model_results(model_list)
  #' @return list of dataframes. First column of each is the dependent variable. Remaining columns are derived from coef and confint.
  merge_coef_CI_tables(get_coefs(model_list),
                       get_CIs(model_list))
}

# main function for getting coefficients, p, CI from a 
batch_models <- function(dependent.var.list,
                         primary.covar,
                         other.covars.str,
                         data) {
  #' Batching sensitivity analyses
  #' 
  #' Runs lmer on given data using formulae generated from the first three arguments.
  #' 
  #' @param dependent.var.list List of dependent variables of interest
  #' @param primary.covar The covariate of primary interest
  #' @param other.covars.str A string of the remaining covariates in the right-hand side of the desired formula
  #' @param data a dataframe to run lmer on
  #' @usage batch_models(dependent.var.list, primary.covar, other.covars.str, data)
  #' @return See get_model_results
  formulas <- make.formulas(dependent.var.list, primary.covar, other.covars.str)
  model_list <- lapply(formulas$formulaStr, \(x) lmer(as.formula(x), data=data))
  names(model_list) <- formulas$dep.var
  get_model_results(model_list = model_list)
}


########################## ###
# 7. Sensitivity Modeling ####
########################## ###

# Intermediate model-- Most useful not in and of itself, but for cross-checking
# the full models, which should be somewhere in the same ballpark.

cont_tables.im <- batch_models(sensitivity.list,
                               "gunOwnership",
                               model.mi.covariates,
                               df)
qrt_tables.im <- batch_models(sensitivity.list,
                              "qrt",
                              model.mi.covariates,
                              df)

# Same method can be applied to the full model on an unimputed (complete case)
# analysis -- see comment at top of this section
cont_tables.full.no_imputation <- batch_models(sensitivity.list,
                                               "gunOwnership",
                                               model.full.covariates,
                                               df)
qrt_tables.full.no_imputation <- batch_models(sensitivity.list,
                                              "qrt",
                                              model.full.covariates,
                                              df)

# individual models done below.  The regular full model is re-integrated
# here and re-done, for ease of comparing estimates to it
sensitivity.list <- append(sensitivity.list, 'firearms', after=0)

formulas.full.cont <- make.formulas(sensitivity.list,
                                    "gunOwnership",
                                    model.full.covariates)
formulas.full.qrt <- make.formulas(sensitivity.list,
                                   "qrt",
                                   model.full.covariates)

# fit both models to the imputation (gunOwnership, qrt)
fit_pmm_list.cont <- lapply(formulas.full.cont$formulaStr, \(x)
                            with(imp1_pmm, lmer(x)))
names(fit_pmm_list.cont) <- formulas.full.cont$dep.var

fit_pmm_list.qrt <- lapply(formulas.full.qrt$formulaStr, \(x)
                           with(imp1_pmm, lmer(x)))
names(fit_pmm_list.qrt) <- formulas.full.cont$dep.var

# pool coefficients of both models
pool_pmm_list.cont <- lapply(fit_pmm_list.cont, pool)
pool_pmm_list.qrt <- lapply(fit_pmm_list.qrt, pool)

# Extract summary data with p and 95% CI values
cont_tables.full <- lapply(pool_pmm_list.cont, summary, conf.int=TRUE)
qrt_tables.full <- lapply(pool_pmm_list.qrt, summary, conf.int=TRUE)

# Output for manual construction of eTable 4
options(width=120)

sink(file="results/sensitivity_models_continuous.txt")
cont_tables.full
sink()

sink(file="results/sensitivity_models_qrt.txt")
qrt_tables.full
sink()

options(width=89)


########################### ###
# 8. Generate Forest Plots ####
########################### ###

# gather and plot gunOwnership estimates
gunOwnership_estimates <- data.frame(
  model_name = names(cont_tables.full),
  mean = sapply(cont_tables.full, \(x) x[x$term=="gunOwnership","estimate"]),
  lower = sapply(cont_tables.full, \(x) x[x$term=="gunOwnership","2.5 %"]),
  upper = sapply(cont_tables.full, \(x) x[x$term=="gunOwnership","97.5 %"]),
  p.value = sapply(cont_tables.full, \(x) x[x$term=="gunOwnership","p.value"])
)

# pretty p-values for the table
gunOwnership_estimates <- gunOwnership_estimates %>%
  mutate(p.replaced = ifelse(p.value >= 0.01,
                             as.character(round(p.value, 2)),
                             ifelse(p.value >= 0.001,
                                    "< 0.01",
                                    "< 0.001")))

# forest plot of gunOwnership estimate across all models
forestplot(x = gunOwnership_estimates,
           labeltext = c(model_name, p.replaced),
           xlog = TRUE,
           boxsize = 0.2,
           vertices = TRUE,
           zero=gunOwnership_estimates[1,2],
           clip = c(10,60),
           xlab = "Estimate",
           xticks = log(c(10, 20, 25, 30, 40, 60)),
           title = "Forest Plot of gunOwnership Estimates\nFull Model vs All Sensitivity Models") |>
  fp_add_lines(h_2 = gpar(lwd = 1), 
               h_3 = gpar(lwd = 1), 
               h_11 = gpar(lwd = 1)) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(model_name = c("Model"),
                p.replaced = c("p value")
                ) |> 
  fp_set_zebra_style("#EFEFEF")

# Gather and plot quartile estimates.
# (It's more complicated than gunOwnership because we are collecting multiple
# values from each data frame in the list of data frames, and we need to track
# which group they came from)
get_quartile_table <- function(data, qrt_name, group_name) {
  data.frame(
    model_name = names(data),
    mean       = sapply(data, \(x) x[x$term==qrt_name,"estimate"]),
    lower      = sapply(data, \(x) x[x$term==qrt_name,"2.5 %"]),
    upper      = sapply(data, \(x) x[x$term==qrt_name,"97.5 %"]),
    p.value    = sapply(data, \(x) x[x$term==qrt_name,"p.value"]),
    qrt_group  = group_name
  )
}

qrtList = c("qrtSecond", "qrtThird", "qrtHighest")
groupList = c("2nd Quartile", "3rd Quartile", "Highest Quartile")

qrt_estimates <- mapply(get_quartile_table,
                        qrtList, groupList,
                        MoreArgs = list(data=qrt_tables.full),
                        SIMPLIFY = FALSE) %>%
  bind_rows %>%
  mutate(p.replaced = ifelse(p.value >= 0.01,
                             as.character(round(p.value, 2)),
                             ifelse(p.value >= 0.001,
                                    "< 0.01",
                                    "< 0.001")),
         .after = p.value)
rownames(qrt_estimates) <- NULL

qrt_est_title <- "Forest Plot of Quartile Estimates\nFull Model vs All Sensitivity Models"
qrtVertLines <- fpShapesGp(grid = list(gpar(lty = 2, col = "blue"),
                                       gpar(lty = 2, col = "green"),
                                       gpar(lty = 2, col = "darkred")))
qrt_estimates |>
  group_by(qrt_group) |>
  forestplot(labeltext = model_name,
             xlab = "Estimate",
             boxsize = 0.2,
             vertices = TRUE,
             xticks = seq(-2, 20, 2),
             title = qrt_est_title,
             grid = c(qrt_estimates[1,2],
                      qrt_estimates[10,2],
                      qrt_estimates[19,2]),
             shapes_gp = qrtVertLines
             ) |>
  fp_set_style(box = c("blue", "green", "darkred"),
               line = c("blue", "green", "darkred")) |>
  fp_add_lines(h_2 = gpar(lwd = 1), 
               h_3 = gpar(lwd = 1), 
               h_11 = gpar(lwd = 1)) |>
  fp_add_header(model_name = c("Model")) |> 
  fp_set_zebra_style("#EFEFEF")

# gather MHProv values from both the gunOwnership and quartiles results
MHProv_estimates.cont <- data.frame(
  model_name = names(cont_tables.full),
  mean       = sapply(cont_tables.full, \(x) x[x$term=="MHProv","estimate"]),
  lower      = sapply(cont_tables.full, \(x) x[x$term=="MHProv","2.5 %"]),
  upper      = sapply(cont_tables.full, \(x) x[x$term=="MHProv","97.5 %"]),
  p.value    = sapply(cont_tables.full, \(x) x[x$term=="MHProv","p.value"]),
  MHP.group  = "Covariate with gunOwnership (Continuous)"
)

MHProv_estimates.qrt <- data.frame(
  model_name = names(qrt_tables.full),
  mean       = sapply(qrt_tables.full, \(x) x[x$term=="MHProv","estimate"]),
  lower      = sapply(qrt_tables.full, \(x) x[x$term=="MHProv","2.5 %"]),
  upper      = sapply(qrt_tables.full, \(x) x[x$term=="MHProv","97.5 %"]),
  p.value    = sapply(qrt_tables.full, \(x) x[x$term=="MHProv","p.value"]),
  MHP.group  = "Covariate with quartiles"
)

MHProv_estimates <- bind_rows(MHProv_estimates.cont,
                              MHProv_estimates.qrt) %>%
  mutate(p.replaced = ifelse(p.value >= 0.01,
                             as.character(round(p.value, 2)),
                             ifelse(p.value >= 0.001,
                                    "< 0.01",
                                    "< 0.001")),
         .after = p.value)
rownames(MHProv_estimates) <- NULL

MHP_est_title <- "MHProv Beta in Continuous & Quartiled Models\nFull Model vs All Sensitivity Models"
MHProv_estimates |>
  group_by(MHP.group) |>
  forestplot(labeltext = model_name,
             xlab = "Estimate",
             boxsize = 0.2,
             vertices = TRUE,
             title = MHP_est_title,
             grid = c(qrt_estimates[1,2],
                      qrt_estimates[10,2],
                      qrt_estimates[19,2]),
             shapes_gp = qrtVertLines) |>
  fp_set_style(box = c("blue", "green"),
               line = c("blue", "green")) |>
  fp_add_lines(h_2 = gpar(lwd = 1), 
               h_3 = gpar(lwd = 1), 
               h_11 = gpar(lwd = 1)) |>
  fp_add_header(model_name = c("Model")) |> 
  fp_set_zebra_style("#EFEFEF")

# The End #