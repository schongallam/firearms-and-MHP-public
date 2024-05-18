# Firearms and mental health providers
# part 2: preliminary inspection, comparison
# 2024-05-17
# Kelly Corbett and Malcolm Schongalla

# main project folder:                     .../firearms-and-MHP
# R scripts folder (assumed working dir):  .../firearms-and-MHP/R
# data folder:                             .../firearms-and-MHP/data

# Script index:
# 01 - (import data) load data from external sources
# 02 - (this script) inspect, compare, create Table 1 & eTable 3
# 03 - (impute & model) correlations, linear models, and [e]Tables 2
# 04 - (sensitivities) - compares imputation techniques and sensitivity models

# Script 02 sequence:
# 1. PRG Seed, Libraries & Load Data
# 2. Examine missingness
# 3. Generate and save Table 1
# 4. Generate and save eTable 3 - censored vs uncensored


##################################### ###
# 1. PRG Seed, Libraries & Load Data ####
##################################### ###

# Pick one or the other:
set.seed(123) # for reproducibility in testing
# set.seed(Sys.time()) # for randomness

library(finalfit)
library(gtsummary)
library(dplyr)

# bridge from prior-generated data from part 1
df <- readRDS(file="data/initial_import.Rds")

######################### ###
# 2. Examine missingness ####
######################### ###

sapply(df, function(x) sum(is.na(x)))
# Observe these missingness counts.  Tend to be large, > 700 NAs, or small, < 8.

missing_ratio <- sapply(df, function(x) sum(is.na(x)))/nrow(df)
missing_table <- sort(missing_ratio[missing_ratio > 0], decreasing = TRUE)
# Notable missing results:
#  firearms (primary outcome) is 26% missing;
#  segregation index -> 33%.  This is a candidate for imputation.
#  homicide ~> 58%.  Homicides are outside of scope, and missing too many for
#   imputation.
#  suicide ~> 23%.  Could be imputed, but also outside of modeling scope.
#
# The "small missing" categories could be addressed with either imputation or 
# discarding.  We favor discarding in these cases, because there is unlikely to
# be a large bias introduction from discarding this relatively small number.
# On the other hand, imputation risks introducing noise and other biases, and
# adds complexity.  See script 4 for relevant sensitivity analyses.

############################# ###
# 3. Generate & save Table 1 ####
############################# ###

# order of elements in tab1_include and tab1_rownames correspond to each other.
tab1_include <- c('firearms',
                  'gunOwnership',
                  'poverty',
                  'homicide',
                  'suicide',
                  'unemployed',
                  'grad',
                  'segregation',
                  'income',
                  'rural',
                  'uninsured',
                  'MHProv',
                  'PCP',
                  'white',
                  'black',
                  'hisp')

tab1_rownames <- c("Firearm fatality rate (per 100,000)",
                   "Proportion of firearm ownership",
                   "Counties with > 20% poverty",
                   "Homicide rate (per 100,000)",
                   "Suicide rate (per 100,000)",
                   "Unemployed (%)",
                   "High school graduation or equivalent (%)",
                   "Residential segregation index",
                   "Median household income ($)",
                   "Proportion of population living in rural areas (%)",
                   "Without health insurance (%)",
                   "Mental health providers ratio",
                   "Primary care physicians ratio",
                   "Non-Hispanic White",
                   "Non-Hispanic Black",
                   "Hispanic")

tab1_vars <- c("group1","group2","group3","group4","group5")
tab1_var_labels <- c("Community Safety",
                     "Health Factors",
                     "Socioeconomic Factors",
                     "Racial & Ethnic Distribution (%)",
                     "Firearms Per State (2016)")

tab1_roworder <- c('group1',
                   'firearms',
                   'homicide',
                   'suicide',
                   'group2',
                   'MHProv',
                   'PCP',
                   'uninsured',
                   'group3',
                   'unemployed',
                   'grad',
                   'income',
                   'poverty',
                   'segregation',
                   'rural',
                   'group4',
                   'white',
                   'black',
                   'hisp',
                   'group5',
                   'gunOwnership')

#generate formulas for tbl_summary to use to label each row
tab1_labels <- mapply(function(x, y) as.formula(sprintf("%s ~ \'%s\'",x, y)),
                      tab1_include, tab1_rownames, USE.NAMES=FALSE)

table1 <- df %>%
  tbl_summary(include = all_of(tab1_include),
              label = tab1_labels,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c(
                "- {N_nonmiss}",
                "- {mean} ({sd})",
                "- {median} ({p25}, {p75})",
                "- {min}, {max}"),
              missing = "no") %>%
  modify_caption("**Table 1: Baseline Characteristics of US Counties & States**") %>%
  modify_header(label="**Baseline Characteristics per County in 2023**",
                stat_0="**Summary:<br>N={N}**") %>%
  modify_table_body(~ .x %>%
                      add_row(variable=tab1_vars,
                              var_label=tab1_var_labels,
                              row_type="label",
                              label=tab1_var_labels) %>%
                      arrange(factor(variable, levels=tab1_roworder))) %>%
  modify_footnote(update=label~'Source: Country Health Rankings &amp; Roadmaps^14^') %>%
  modify_table_styling(columns = label,
                       rows = variable=='segregation',
                       footnote = "Index of 0 (complete integration) to 100 (complete segregation)") %>%
  modify_table_styling(columns = label,
                     rows = variable %in% c('MHProv','PCP'),
                     footnote = "Providers per 1000 residents") %>%
  modify_table_styling(columns = label,
                       rows = variable=='group5',
                       footnote = 'Source: “State-level estimates of household firearm ownership”^16^') %>%
  modify_column_indent(columns = label,
                       rows = !(variable %in% tab1_vars))

# The next line fails to bold the rows if part of the previous pipe
table1 <- modify_table_styling(table1, columns = label,
                               rows = .data$variable %in% tab1_vars,
                               text_format = "bold")

table1_gt <- as_gt(table1)
table1_gt

gt::gtsave(table1_gt, filename='table1.rtf', path='results/')
gt::gtsave(table1_gt, filename='table1.pdf', path='results/')


######################################################### ###
# 4. Generate and save eTable 3 - censored vs uncensored ####
######################################################### ###

cont_candidates <- c('injury.denominator',
                     'income',
                     'segregation',
                     'MHProv',
                     'PCP',
                     'uninsured',
                     'unemployed',
                     'grad',
                     'hisp',
                     'black',
                     'white',
                     'rural')

chisq_candidates <- c('poverty_c')

eTab3_include <- c(cont_candidates, chisq_candidates)

eTab3_rownames <- c("County population",
                   "Income ($)",
                   "Residential segregation index",
                   "Mental health providers ratio",
                   "Primary care physicians ratio",
                   "Uninsured rate (%)",
                   "Unemployed rate (%)",
                   "High school grad. or equiv. (%)",
                   "Hispanic (%)",
                   "Non-Hispanic Black (%)",
                   "Non-Hispanic White (%)",
                   "Rurality (%)",
                   "Poverty index")

# Generate the associations between row variables and row names (labels)
eTab3_labels <- mapply(function(x, y) as.formula(sprintf("%s ~ \'%s\'",x, y)),
                      eTab3_include, eTab3_rownames, USE.NAMES=FALSE)

# creating these columns makes the tbl_summary product look nicer
df <- df %>%
  mutate(subset=ifelse(censored==TRUE,'Censored','Uncensored')) %>%
  mutate(poverty_c=ifelse(poverty==1,'>20%','<20%'))

# Helps with column & row order in tbl_summary
df$subset <- factor(df$subset, levels=c('Uncensored','Censored'), ordered=TRUE)
df$poverty_c <- factor(df$poverty_c, levels=c(">20%","<20%"), ordered=TRUE)

eTab3_vars <- c('group1', 'group2', 'group3', 'group4')
eTab3_var_labels <- c('General',
                     'Health Factors',
                     'Socioeconomic Factors',
                     'Racial & Ethnic Distribution (%)')
eTab3_roworder <- c('group1', # 'General'
                   'injury.denominator',
                   'group2', # 'Health Factors'
                   'MHProv',
                   'PCP',
                   'uninsured',
                   'group3', # 'Socioeconomic Factors'
                   'unemployed',
                   'grad',
                   'income',
                   'poverty_c', # note different from table 1
                   'segregation',
                   'rural',
                   'group4', # 'Racial & Ethnic Distribution (%)'
                   'white',
                   'black',
                   'hisp')

eTable3 <- df %>%
  tbl_summary(
    include = all_of(eTab3_include),
    label = eTab3_labels,
    #statistic = list(
    #  all_continuous() ~ "{median} ({p25}, {p75})",
    #  all_categorical() ~ "{n} / {N} ({p}%)"),
    type = all_continuous() ~ "continuous2",
    statistic = all_continuous() ~ c(
      "- {N_nonmiss}",
      "- {mean} ({sd})",
      "- {median} ({p25}, {p75})",
      "- {min}, {max}"),
    by=subset,
    digits=list(injury.denominator ~ c(0,2),
                income ~ c(0,2)),
    missing='no') %>%
  add_p(test = list(all_continuous() ~ "t.test"),
        pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  modify_caption("**eTable 3: Comparing Uncensored vs. Censored Subsets of Firearms Deaths**") %>%
  modify_footnote(update=c('stat_1','stat_2')~'Median (IQR); n / N (%); "Censored" are counties with 1 to 9 firearms deaths.  "Uncensored" are counties with 0 or >9 firearm deaths.') %>%
  modify_table_styling(columns = label,
                       rows = variable=='segregation',
                       footnote = "Index of 0 (complete integration) to 100 (complete segregation)") %>%
  modify_table_styling(columns = label,
                       rows = variable %in% c('MHProv','PCP'),
                       footnote = "Providers per 1000 residents") %>%
  modify_table_body(~ .x %>% # thanks to Guilherme, https://stackoverflow.com/questions/65665465/grouping-rows-in-gtsummary
                      add_row(variable=eTab3_vars,
                              var_label=eTab3_var_labels,
                              row_type='label',
                              label=eTab3_var_labels) %>%
                      arrange(factor(variable, levels=eTab3_roworder))) %>%
  modify_column_indent(columns = label,
                       rows = !(variable %in% eTab3_vars)) %>%
  modify_column_indent(columns = label,
                       rows = (row_type == "level"),
                       double_indent=TRUE)

# The next line fails to bold the rows if part of the pipe
eTable3 <- modify_table_styling(eTable3, columns = label,
                               rows = .data$variable %in% eTab3_vars,
                               text_format = "bold")

eTable3_gt <- as_gt(eTable3)

eTable3_gt

gt::gtsave(eTable3_gt, filename='eTable3.rtf', path='results/')
gt::gtsave(eTable3_gt, filename='eTable3.pdf', path='results/')

# The End #