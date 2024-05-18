# Firearms and mental health providers
# part 3: imputation and modeling
# 2024-05-17
# Kelly Corbett and Malcolm Schongalla

# main project folder:                     .../firearms-and-MHP
# R scripts folder (assumed working dir):  .../firearms-and-MHP/R
# data folder:                             .../firearms-and-MHP/data

# Script index:
# 01 - (import data) load data from external sources
# 02 - (describe data) inspect, compare, create Table 1 & eTable 3
# 03 - (this script) correlations, linear models, and [e]Tables 2
# 04 - (sensitivities) - compares imputation techniques and sensitivity models

# Script 3 sequence:
# 1. Libraries & load data
# 2. Pearson's correlation and Figures 1, 2
# 3. Crude & Intermediate linear modeling (Tables 2a, 2b)
# 4. Generate sensitivity vectors
# 5. Initialize imputation and build parameters
# 6. Run imputations and check the products
# 7. Define and run full models
# 8. Complete Tables 2c & 2d and merge

########################### ###
# 1. Libraries & Load Data ####
########################### ###

# Pick one:
set.seed(123)
# set.seed(Sys.time())

library(dplyr)
library(gt)
library(ggpubr) # ggscatter. depends ggplot2
library(broom.mixed) # tidy methods not auto-loaded by all other parts
library(gtsummary) # tbl_summary helper fxns
library(lmerTest) # depends lme4
library(mice)
library(howManyImputations)

requireNamespace("performance") # icc
requireNamespace("car") # vif

df <- readRDS(file="data/initial_import.Rds")

############################################# ###
# 2. Pearson's correlation and Figures 1, 2 #####
############################################# ###

# Pearson's correlations: gunOwnership ~ firearms, and MHProv ~ firearms
cor.test(df$firearms, y=df$MHProv)
cor.test(df$firearms, y=df$gunOwnership)

# Figure 1
cor1 <- df %>% 
  select(firearms, MHProv)

Figure1 <- ggscatter(
  cor1,
  size = 0.5,
  x = "MHProv", 
  y = "firearms",
  add = "reg.line",  #add linear regression line
  xlab = "Ratio of MH Providers in the population per 100,000 pop.", 
  ylab = "Firearm Fatality Rate (per 100,000)",
  title = "Figure 1: Scatterplot and Pearson's correlation of Mental Health\nProvider Density and Firearm Fatality Rate",
  caption = "Note: Each observation represents a county's firearm fatality rate."
) + stat_cor(method="pearson", p.accuracy = 0.001, label.x.npc="center")
Figure1

pdf(file="results/Fig1 Scatterplot.pdf",
    width=7,
    height=7,
    family="Helvetica")
  print(Figure1)
dev.off()

Figure1_notitle <- ggscatter(
  cor1,
  size = 0.5,
  x = "MHProv", 
  y = "firearms",
  add = "reg.line",  #add linear regression line
  xlab = "Ratio of MH Providers in the population per 100,000 pop.", 
  ylab = "Firearm Fatality Rate (per 100,000)",
#  title = "Figure 1: Scatterplot and Pearson's correlation of Mental Health\nProvider Density and Firearm Fatality Rate",
#  caption = "Note: Each observation represents a county's firearm fatality rate."
) + stat_cor(method="pearson", p.accuracy = 0.001, label.x.npc="center")
#Figure1_notitle

pdf(file="results/Fig1 Scatterplot(no title).pdf",
    width=7,
    height=7,
    family="Helvetica")
  print(Figure1_notitle)
dev.off()

# Figure 2
cor2 <- df %>% 
  select(firearms, gunOwnership)

Figure2 <- ggscatter(
  cor2, 
  size = 0.5,
  x = "gunOwnership", 
  y = "firearms",
  add = "reg.line",  #add linear regression line
  xlab = "Proportion of Population Identied as Firearm Owners", 
  ylab = "Firearm Fatality Rate (per 100,000)",
  title = "Figure 2: Scatterplot and Pearson's correlation of Firearm Ownership\nand Firearm Fatality Rate",
  caption = "Note: Each observation represents a county's firearm fatality rate, grouped in columns based on the state's\nproportion of firearm ownership. The state firearm ownership rate was applied equally to each county\nwithin the state in the absence of county-specific data."
) + stat_cor(method="pearson", p.accuracy = 0.001)
Figure2

pdf(file="results/Fig2 Scatterplot.pdf",
    width=7,
    height=7,
    family="Helvetica")
  print(Figure2)
dev.off()

Figure2_notitle <- ggscatter(
  cor2, 
  size = 0.5,
  x = "gunOwnership", 
  y = "firearms",
  add = "reg.line",  #add linear regression line
  xlab = "Proportion of Population Identied as Firearm Owners", 
  ylab = "Firearm Fatality Rate (per 100,000)",
#  title = "Figure 2: Scatterplot and Pearson's correlation of Firearm Ownership\nand Firearm Fatality Rate",
#  caption = "Note: Each observation represents a county's firearm fatality rate, grouped in columns based on the state's\nproportion of firearm ownership. The state firearm ownership rate was applied equally to each county\nwithin the state in the absence of county-specific data."
) + stat_cor(method="pearson", p.accuracy = 0.001)
#Figure2_notitle

pdf(file="results/Fig2 Scatterplot(no title).pdf",
    width=7,
    height=7,
    family="Helvetica")
  print(Figure2_notitle)
dev.off()

# Note: If PNG output is desired, 600x600 is passable, but font size should be
# adjusted in the ggscatter function


########################################################### ###
# 3. Crude & Intermediate linear modeling (Tables 2a, 2b) #####
########################################################### ###

# Nomenclature:
#  'mcX' - crude model, version X
#  'mi' - intermediate model
#  'full' or 'pmm' - full model after multiple imputation by pmm method

df$qrt <- factor(df$qrt, labels=c('Lowest','Second','Third','Highest'))

# Crude analysis: continuous exposure for guns, without mixed effect 
mc1_gunOwnership <- lm(firearms ~ gunOwnership,
                       data=df) 
summary(mc1_gunOwnership)
confint(mc1_gunOwnership)

# Same, but for mental health providers
mc1_MHProv <- lm(firearms ~ MHProv,
                 data=df) 
summary(mc1_MHProv)
confint(mc1_MHProv)

# Crude analysis: continuous exposure for guns, with mixed effect of state
mc2 <- lmer(firearms ~ gunOwnership + (1|state_fips),
            data=df)
summary(mc2)

# Crude analysis: for Table 2, by quartile with mixed effects
mc3 <- lmer(firearms ~ qrt + (1|state_fips),
            data=df)
summary(mc3)

df$injury.denominator.scaled <- df$injury.denominator / 100000

# Final crude model, correcting for county size and mixed effect of state:
mc4 <- lmer(firearms ~ qrt + injury.denominator.scaled + (1|state_fips),
            data=df)
summary(mc4)

table2_mc4.footnote <- "Proportion of state population that owns firearms, categorized by quartiles."
table2_mc4 <- tbl_regression(mc4,
                             exponentiate = FALSE,
                             label = list(qrt ~"Firearm ownership quartiles by state"),
                             include = c(qrt)) %>%
  modify_caption("Table 2a. Unadjusted Model<br>(Corrected for population only)") %>%
  modify_table_styling(columns = label,
                       rows = label == "Firearm ownership quartiles by state",
                       footnote = table2_mc4.footnote)
table2_mc4

# Intermediate model
# It may seem like extra work to separate out the model formulas here, but it
# adds value: it is re-used in script 4, standardizing it yields consistency.

model.mi.covariates <- list("MHProv",
                            "injury.denominator.scaled",
                            "(1 | state_fips)") |>
  paste(collapse=" + ")

formulaStr.mi.continuous <- paste("firearms ~ gunOwnership",
                                  model.mi.covariates,
                                  sep=" + ")
formulaStr.mi.qrt        <- paste("firearms ~ qrt",
                                  model.mi.covariates,
                                  sep=" + ")

# "firearms ~ gunOwnership + MHProv + injury.denominator.scaled + (1 | state_fips)"
mi.continuous <- lmer(as.formula(formulaStr.mi.continuous),
                      data=df)

# "firearms ~ qrt + MHProv + injury.denominator.scaled + (1 | state_fips)"
mi.qrt <- lmer(as.formula(formulaStr.mi.qrt),
               data=df)

summary(mi.continuous)
summary(mi.qrt)

mi <- mi.qrt
table2_mi.main  <- "Table 2b. Intermediate Model with Mental Health Providers"
table2_mi.sub   <- "Corrected for population"
table2_mi.title <- paste(table2_mi.main, table2_mi.sub, sep="<br>")
table2_mi.footnote.qrt <- "Proportion of state population that owns firearms, categorized by quartiles."
table2_mi.footnote.MHProv <- "Ratio of primary care physicians per 1000 people."
table2_mi <-
  tbl_regression(mi,
                 exponentiate = FALSE,
                 label = list(qrt ~ "Firearm ownership quartiles by state",
                              MHProv ~ "Mental health providers ratio"),
                 include = c(qrt, MHProv)) %>%
#                              injury.denominator.scaled ~ "County Population")) %>%
  modify_caption(table2_mi.title) %>%
  modify_table_styling(columns = label,
                       rows = label == "Firearm ownership quartiles by state", # use 'lable ==' here to keep the footnote off of "Lowest...Highest" rows
                       footnote = table2_mi.footnote.qrt) %>%
  modify_table_styling(columns = label,
                       rows = variable == "MHProv",
                       footnote = table2_mi.footnote.MHProv)
  
table2_mi # to be stacked later

performance::icc(mi.continuous)
performance::icc(mi.qrt)

car::vif(mi.continuous)
car::vif(mi.qrt) # values around 1, reassures against problematic collinearity

################################## ###
# 4. Generate sensitivity vectors ####
################################## ###

# Add sensitivity analysis columns now, so that they are present in the final
# mids product.  Allows for running the sensitivity analysis on the full model.

# 3/8/2024 note: After inspecting the various sensitivity sets, it is apparent
# that at least one county (FIPS 48301) becomes a major outlier.  This includes
# at least one county with a low population count that is censored-- i.e., at
# least 1 firearm death. This one small county having a minimum 1 death adds a
# firearm death rate to the dataset that is far greater than (~ twice) the next
# highest county.  It raises the question - should these small counties have the
# same impact that a big county has?  Statistically, no.

# Overiew of each sensitivity set:
# - pop (population-proportional: replace NAs with a range of 1:9 depending on
#   injury.denominator)
# - low (replace NAs with 1)
# - mid (replace NAs with 5) -- limited utility
# - hi (replace NAs with 9)
# - proportional (replace NAs with a range of 1:9 depending on gunOwnership)
# - invProportional (inverse of proportional)
# - invMHP (replace NAs with a range of 1:9 depending inversely on MHProv)
#   -- i.e. assumes MHProv is inversely correlated with firearms deaths

# Building the 'population-proportional' sensitivity model: Assumes the largest
# censored county would have 9 firearms deaths, which is reasonable by
# inspection of the data.

popMax <- max(df[is.na(df$firearms.numerator),]$injury.denominator) # 243457
popMin <- min(df[is.na(df$firearms.numerator),]$injury.denominator) # 749

df$firearms.numerator.pop <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  (8*(df$injury.denominator-popMin)/(popMax-popMin))+1)
df$firearms.pop <- (df$firearms.numerator.pop / df$injury.denominator)*100000

df$firearms.numerator.lo <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  1)
df$firearms.lo <- (df$firearms.numerator.lo / df$injury.denominator)*100000

df$firearms.numerator.mid <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  5)
df$firearms.mid <- (df$firearms.numerator.mid / df$injury.denominator)*100000

df$firearms.numerator.hi <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  9)
df$firearms.hi <- (df$firearms.numerator.hi / df$injury.denominator)*100000

# Build 'proportional' set
gOmax <- max(df$gunOwnership)
gOmin <- min(df$gunOwnership)

df$firearms.numerator.proportional <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  (8*(df$gunOwnership-gOmin)/(gOmax-gOmin))+1)
df$firearms.proportional <-
  (df$firearms.numerator.proportional / df$injury.denominator)*100000

# Build 'invProportional' set
df$firearms.numerator.invProportional <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  10-((8*(df$gunOwnership-gOmin)/(gOmax-gOmin))+1))
df$firearms.invProportional <-
  (df$firearms.numerator.invProportional / df$injury.denominator)*100000

# Build 'MHP' and 'invMHP' sets
MHPmax <- max(df[is.na(df$firearms),]$MHProv)
MHPmin <- min(df[is.na(df$firearms),]$MHProv)

df$firearms.numerator.invMHP <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  10-(8*(df$MHProv-MHPmin)/(MHPmax-MHPmin)+1))
df$firearms.invMHP <- (df$firearms.numerator.invMHP / df$injury.denominator)*100000

df$firearms.numerator.MHP <- ifelse(
  !is.na(df$firearms.numerator),
  df$firearms.numerator,
  8*(df$MHProv-MHPmin)/(MHPmax-MHPmin)+1)
df$firearms.MHP <- (df$firearms.numerator.MHP / df$injury.denominator)*100000


################################################ ###
# 5. Initialize imputation and build parameters ####
################################################ ###

# Due to high multicollinearity, imputing NAs in all columns is not feasible.
# Using car::vif and stats:alias (early versions of this work, not shown here)
# on samples of linear models shows many instances of severe correlation &
# linear dependency.

# Based on missing_table (script 2), the following have low rates of
# missingness.  We consider these rates to be negligible and allowable.  No
# imputation done for these categories.
#   rural (0.2%)
#   urban (0.2%) (removed from final analysis as redundant)
#   poverty (0.2%)
#   income (<0.1%)
#   PCP (<0.1%)
#   uninsured (<0.1%)
#   unemployed (<0.1%)

# The following categories have high missingness rates and are not imputed,
# either being outcome measures, or, are ignored as covariates in the subsequent
# modeling:
#   homicide (57.8%)
#   firearms.numerator and firearms(26.0%) (outcome measure)
#   suicide (22.6%)

# Only one category remains, and it supports our modeling needs:
#   segregation (33.7% missing)

# The subsets (uncensored vs. censored) are statistically different
# across many dimensions.  Therefore, imputation of missings in each subset is
# best undertaken by a nonparametric method.

# Imputation steps:
#   Factorize appropriate variables
#   Scale down large order variables as needed (for floating point limitations)
#   For complete data set, for uncensored subset, and for censored subset (each
#   separately):
#     -initialize and customize imputation matrices
#     -impute with method(s) -> pmm and rf are used (see script 4 for rf)
#     -check for convergence and gross mis-imputations
#   compare the imputed values between the complete dataset and the subsets

df$poverty <- factor(df$poverty, labels=c("No","Yes"))

# rank should not be used due to computational cost.  qrt can be used instead,
# in keeping with crude and intermediate models. Uncomment to include.

#df$rank <- as.ordered(df2$rank)

# Scaling down the large-value columns avoids floating-point precision issues,
# (not necessarily a problem for random forest, but it allows for flexibility
# in statistical design)

#df$injury.denominator <- df$injury.denominator/100000 # already done
df$income.scaled <- df$income/10000

# initialize mice, check logged events
imp0 <- mice(df, maxit=0)
imp0$loggedEvents
# note collinear firearms.denominator, plus constant covariates
# this means injury.denominator also needs to be excluded

# pmm was selected for imputation.  See comparison vs rf in script 4.
meth <- imp0$method
meth[names(meth)] <- ""

meth_pmm <- meth
meth_pmm["segregation"] <- "pmm"

# Set up the prediction matrix.  quickpred sets up the matrix, with explicit
# exclusions, then we just set everything not the "segregation" row to zero.
non_predictors <- c("state_fips", # constant
                    "state_abbrev", # constant
                    "county_name", # constant
                    "county_fips", # constant
                    "firearms.numerator", # primary outcome measure
                    "firearms.numerator.pop", # sensitivity analysis
                    "firearms.numerator.lo", # sensitivity analysis
                    "firearms.numerator.hi", # sensitivity analysis
                    "firearms.numerator.mid", # sensitivity analysis
                    "firearms.numerator.invProportional", # sensitivity analysis
                    "firearms.numerator.proportional", # sensitivity analysis
                    "firearms.numerator.invMHP", # sensitivity analysis
                    "firearms.numerator.MHP", # sensitivity analysis
                    "firearms.denominator", # collinear
                    "injury.denominator", # collinear
                    "injury.denominator.scaled", #collinear
                    "income", # in favor of income.scaled
                    "MHProv", # bias reinforcement
                    "censored", # not an observed value
                    "firearms", # primary outcome statistic
                    "firearms.pop", # sensitivity analysis
                    "firearms.lo", # primary outcome statistic
                    "firearms.hi", # primary outcome statistic
                    "firearms.mid", # primary outcome statistic
                    "firearms.invProportional", # sensitivity analysis
                    "firearms.proportional", # sensitivity analysis
                    "firearms.invMHP", # sensitivity analysis
                    "firearms.MHP", # sensitivity analysis
                    "rank", # collinearity risk
                    "gunOwnership", # collinearity risk
                    "qrt") # not an observed value

pred <- quickpred(df, exclude=non_predictors)

#zero out predictors for all rows except segregation
pred[rownames(pred) != "segregation",] <- 0

# remaining predictors for segregation:
# income.scaled, PCP, homicide, suicide, uninsured, unemployed, grad, black, white,
# rural, totalGuns, poverty
#
# Empirically, selectively excluding (homicide | suicide | grad) values impacts
# the final NA count after imputation.  Specifically, removing homicide and
# suicide from prediction matrix has the biggest benefit in improving the
# completion rate in imputation.  Comparing the six possible combinations of
# removing homicide, suicide, and grad: 1059 NA became ->
#-> 529 NA, with '0' for homicide
#-> 1010 NA, with '0' for grad
#-> 1009 NA, with '0' for suicide
#-> 524 NA, with '0' for homicide & grad
#-> 1009 NA, with '0' for suicide & grad
#-> 16 NA, with '0' for homicide & suicide
#-> 16 NA, with '0' for grad, homicide, and suicide
pred["segregation",c("homicide","suicide")] <- 0


############################################ ###
# 6. Run imputations and check the products ####
############################################ ###

# Run pmm imputation
imp.ct <- 33
imp1_pmm <- mice(df,
                 m=imp.ct,
                 maxit = 5,
                 predictorMatrix = pred, 
                 method = meth_pmm,
                 print = TRUE)

imp1_pmm.stripplot <-stripplot(imp1_pmm, segregation ~ .imp, pch=20, cex=2)
imp1_pmm.stripplot

pdf(file="results/mice_stripplot_pmm.pdf",
    width=9,
    height=9,
    family="Helvetica")
  print(imp1_pmm.stripplot)
dev.off()

# check for convergence
conv.plot.title <- "mice convergence check (pmm)"
conv.plot.sub <- sprintf("Segregation Index(0-100), iterations=%d", imp.ct)
imp1_pmm.convCheck <- plot(imp1_pmm,
     main = conv.plot.title,
     sub = conv.plot.sub,
     layout=c(1,2))
imp1_pmm.convCheck

pdf(file="results/mice_conv_pmm.pdf", width=9, height=9, family="Helvetica")
  print(imp1_pmm.convCheck)
dev.off()

imp1_pmm.densityPlot <- densityplot(
  imp1_pmm, ~segregation | .imp,
  main="mice densityplot of imputations (pmm)")
imp1_pmm.densityPlot

pdf(file="results/mice_dp_pmm.pdf", width=9, height=9, family="Helvetica")
  print(imp1_pmm.densityPlot)
dev.off()


################################ ###
# 7. Define and run full models ####
################################ ###

# Build formula from a defined list in the interest of consistency and avoiding
# typo-driven bugs.  This adds value because the formula was built progressively
# during study development, and the model is run several times on different data
# sets, across at least two scripts.

model.full.covariates <- list("MHProv",
                              "rural",
                              "PCP",
                              "unemployed",
                              "uninsured",
                              "grad",
                              "segregation",
                              "poverty",
                              "injury.denominator.scaled",
                              "(1 | state_fips)") |>
  paste(collapse=" + ")

model.full.continuous <- paste("firearms ~ gunOwnership",
                               model.full.covariates,
                               sep=" + ")
model.full.qrt        <- paste("firearms ~ qrt",
                               model.full.covariates,
                               sep=" + ")

fit_pmm <- with(imp1_pmm, lmer(as.formula(model.full.continuous)))

how_many_imputations(fit_pmm) #24 if run on imp(m=33), but #111 if m=5
# Note, cycling imputation iterations and how_many_imputations will not
# asymptote toward zero.

fit_pmm_qrt <- with(imp1_pmm, lmer(as.formula(model.full.qrt)))
how_many_imputations(fit_pmm_qrt) # less than m

pool_pmm <- pool(fit_pmm)
summary(pool_pmm)

pool_pmm_qrt <- pool(fit_pmm_qrt)
summary(pool_pmm_qrt)

# Check icc and vif, which we did on the intermediate model.
# It's not reasonable to do this on each and every imputed dataset, so just
# extract and model the first imputed dataset and check on that.
df.imputed <- complete(imp1_pmm)
mf.1.cont <- lmer(as.formula(model.full.continuous),
                  data=df.imputed)
mf.1.qrt <- lmer(as.formula(model.full.qrt),
                 data=df.imputed)
performance::icc(mf.1.cont)
performance::icc(mf.1.qrt)
car::vif(mf.1.cont)
car::vif(mf.1.qrt)
# ok

# Creates a large file (~75mb).  Generate for script 4 "as needed."
save(df, pred, meth_pmm, imp1_pmm, fit_pmm, fit_pmm_qrt, pool_pmm, pool_pmm_qrt,
     model.mi.covariates, formulaStr.mi.continuous, formulaStr.mi.qrt,
     model.full.covariates, model.full.continuous, model.full.qrt,
     file="results/imputation.RData")


###################################### ###
# 8. Complete Table 2c & 2d and merge ####
###################################### ###

# Table 2c lists all covariates. Table 2d is an alternate version that replaces
# multiple rows with an explanatory footnote, for brevity.  The final Table 2
# incorportates 2a, 2b, and 2d.  Table 2c is preserved for eTable 2 (comprised
# of tables 2a, 2b, 2c)

tab2c_include <- c('qrt', # gets footnote
                   'MHProv', # gets footnote
                   'rural',
                   'PCP', # gets footnote
                   'unemployed',
                   'uninsured',
                   'grad',
                   'segregation', # gets footnote
                   'poverty',
                   'injury.denominator.scaled') # gets footnote
tab2c_rownames <- c("Firearm ownership quartiles by state",
                    "Mental health providers ratio",
                    "Proportion of population living in rural areas (%)",
                    "Primary care physicians ratio",
                    "Unemployed (%)",
                    "Without health insurance (%)",
                    "High school graduation or equivalent (%)",
                    "Residential segregation index",
                    "Poverty",
                    "County Population")

tab2c_labels <- mapply(function(x, y) as.formula(sprintf("%s ~ \'%s\'",x, y)),
                        tab2c_include, tab2c_rownames, USE.NAMES=FALSE)

table2_mf <- tbl_regression(fit_pmm_qrt,
                            label = tab2c_labels[1:9],
                            include = tab2c_include[1:9]) %>%
  modify_caption("Table 2c. Full Model<br>Corrected for Population") %>%
  modify_table_styling(columns = label,
                       rows = label == "Firearm ownership quartiles by state",
                       footnote = "Proportion of state population that owns firearms, categorized by quartiles.") %>%
  modify_table_styling(columns = label,
                       rows = variable == "MHProv",
                       footnote = "Ratio of primary care physicians per 1000 people.") %>%
  modify_table_styling(columns = label,
                       rows = variable == "PCP",
                       footnote = "Ratio of primary care physicians per 1000 people.") %>%
  modify_table_styling(columns = label,
                       rows = variable == "segregation",
                       footnote = "Index of dissimilarity, where 0 (complete integration) to 100 (complete segregation)") %>%
  modify_table_styling(columns = label,
                       rows = label == "Poverty",
                       footnote = "USDA Poverty Area Measures metric 'HiPov1519'")
table2_mf


#Concise version of table2_mf (footnote instead of showing all estimates)
table2_mf.concise <- tbl_regression(fit_pmm_qrt,
                                    label = tab2c_labels[1:2],
                                    include = tab2c_include[1:2]) %>%
  modify_caption("Table 2d. Full Model (Concise)") %>%
  
  modify_table_styling(columns = label,
                       rows = label == "Firearm ownership quartiles by state",
                       footnote = "Proportion of state population that owns firearms, categorized by quartiles.") %>%
  modify_table_styling(columns = label,
                       rows = variable == "MHProv",
                       footnote = "Ratio of primary care physicians per 1000 people.")
# Summary footnote: corrections for rurality, PCPs per 1,000 people,
# unemployment, (lack of) health insurance, high school education, residential
# segregation, poverty (USDA Poverty Area Measures metric ‘HiPov1519’)
table2_mf.concise


# table2 is tables a,b,c merged -> "eTable 2" in manuscript
table2 <- tbl_merge(list(table2_mc4, table2_mi, table2_mf),
                    tab_spanner = c("Unadjusted Model",
                                    "Intermediate Model",
                                    "Full Model")) %>%
  modify_caption("eTable 2: Mixed effect linear regression model showing a significant association between the Proportion of Firearm Ownership by State and Firearm Fatality Rate (per 100,000)")

# table2.concise is tables a,b,d -> "Table 2" in manuscript
table2.concise <- tbl_merge(list(table2_mc4, table2_mi, table2_mf.concise),
                            tab_spanner = c("Unadjusted Model",
                                            "Intermediate Model",
                                            "Full Model")) %>%
  modify_caption("Table 2: Mixed effect linear regression model showing a significant association between the Proportion of Firearm Ownership by State and Firearm Fatality Rate (per 100,000)")

table2.footnote1 <- "All models corrected for county population."
table2_gt <- as_gt(table2) %>%
  tab_footnote(table2.footnote1)
table2_gt

table2.footnote2 <- "Full model corrected for rurality, ratio of primary care physicians per 1,000 people, unemployment %, without health insurance %, high school graduation %, Residential segregation index, Poverty indicator (USDA PAM HiPov1519)"
table2_gt.concise <- as_gt(table2.concise) %>%
  tab_footnote(paste(table2.footnote1, table2.footnote2))
table2_gt.concise

gt::gtsave(table2_gt, filename='eTable2.rtf', path='results/')
gt::gtsave(table2_gt, filename='eTable2.pdf', path='results/')
gt::gtsave(table2_gt.concise, filename='Table2.rtf', path='results/')
gt::gtsave(table2_gt.concise, filename='Table2.pdf', path='results/')

# The End #