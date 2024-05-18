# Firearms and Mental Health Providers
# Part 1, data extraction and joining
# County Health Rankings + RAND gun data + usda PAM
# 2024-03-07
# Kelly Corbett and Malcolm Schongalla

# main project folder:                     .../firearms-and-MHP
# R scripts folder (assumed working dir):  .../firearms-and-MHP/R
# data folder:                             .../firearms-and-MHP/data

# Script index:
# 01 - (this script) load data from external sources
# 02 - (describe data) inspect, compare, create Table 1 & eTable 2
# 03 - (impute & model) correlations, linear models, and Table 2
# 04 - (sensitivities) - compares imputation techniques and sensitivity models

# Script 01 sequence:
# 1. Libraries
# 2. Load RAND data on gun ownership by state
# 3. Load USDA poverty area measures
# 4. Load CDC WONDER firearm death numbers by county
# 5. Load County Health Rankings and Roadmaps (CHR&R) dataset
# 6. Apply CDC WONDER data on 0-death counties
# 7. Join, finalize columns, save

############### ###
# 1. Libraries ####
############### ###

library(dplyr)
library(stringr)

############################################## ###
# 2. Load RAND data on gun ownership by state ####
############################################## ###

# RAND gun ownership data by state. Missing Washington DC.
guns0 <- read.csv("data/data.csv") %>%
  select (fips, gunOwnership, totalGuns, rank) %>%
  rename(state_fips = fips) %>%
  mutate(state_fips = str_pad(as.character(state_fips),
                              width=2,
                              side="left",
                              pad="0"))
# guns0 ready to left join

##################################### ###
# 3. Load USDA poverty area measures ####
##################################### ###

# the USDA dataset contains many different measures of poverty.
# HiPov1519 is the most current poverty measure and the most sensitive for
# poverty.
pam0 <- read.csv("data/PovertyAreaMeasures2022.csv") %>%
  select(fips, HiPov1519_c) %>%
  distinct(fips, .keep_all=TRUE) %>% # exclude duplicates
  rename(poverty = HiPov1519_c, county_fips = fips) %>%
  mutate(county_fips=str_pad(as.character(county_fips),
                             width=5,
                             side="left",
                             pad="0"))

table(pam0$poverty) # 7 were NA but coded as -1.  Re-code to NA
pam0$poverty <- replace(pam0$poverty, pam0$poverty<0, NA)
summary(pam0$poverty) # now 7 NA's

### pam0 ready to join

##################################################### ###
# 4. Load CDC WONDER firearm death numbers by county ####
##################################################### ###

# The next file contains zero-death counties that are excluded from CHR&R data.
# Counties with 1-9 deaths are still suppressed.
# Note: Five counties are "Missing" all data for the row.  These can be
# distinguished because "Deaths" and "Population" both are left as NA,
# whereas the suppressed rows are NA "Deaths" but have a value for "Population"
county_deaths <- read.csv("data/CDC WONDER Underlying Cause of Death, 1999-2020.txt",
                          sep="\t",
                          header=TRUE,
                          na.strings=c("Suppressed","Missing"),
                          col.names=c("Notes",
                                      "County",
                                      "County_fips",
                                      "Deaths",
                                      "Population",
                                      "Crude.Rate"),                          
                          colClasses=c("Notes"      ="NULL",
                                       "County"     ="NULL",
                                       "County_fips"="character",
                                       "Deaths"     ="numeric",
                                       "Population"  ="numeric",
                                       "Crude.Rate" ="NULL")) %>%
  filter(nchar(County_fips)==5) %>%
  filter(!is.na(Population))
# county_deaths ready to join

############################################################## ###
# 5. Load County Health Rankings and Roadmaps (CHR&R) dataset ####
############################################################## ###

# County Health Rankings & Roadmaps, injury data by county including firearms
# deaths.  This data set masks both "zero-death" counties, as well as counties
# with suppressed death statistics in the range of 1-9 (for privacy).
# CDC WONDER does actually provide the zero-death counties and we bring
# that in further below.

# set up vectors for reading efficiently
county_cols <- c("State.FIPS.Code", # our columns of interest.  Order counts!
                 "State.Abbreviation", 
                 "Name", 
                 "X5.digit.FIPS.Code",
                 # the following need to be as.numeric-ized:
                 "Firearm.Fatalities.numerator",
                 "Firearm.Fatalities.denominator",       
                 "Injury.Deaths.denominator",
                 "Median.Household.Income.raw.value",
                 "Residential.Segregation...Black.White.raw.value",
                 "Mental.Health.Providers.numerator",
                 "Mental.Health.Providers.denominator",
                 "Primary.Care.Physicians.numerator",
                 "Primary.Care.Physicians.denominator",
                 "Homicides.numerator",
                 "Homicides.denominator",
                 "Suicides.numerator",
                 "Suicides.denominator",
                 "Uninsured.numerator",
                 "Uninsured.denominator",
                 "Unemployment.numerator",
                 "Unemployment.denominator",
                 "High.School.Completion.numerator",
                 "High.School.Completion.denominator",
                 "X..Hispanic.numerator",
                 "X..Hispanic.denominator",
                 "X..Non.Hispanic.Black.numerator",
                 "X..Non.Hispanic.Black.denominator",
                 "X..Non.Hispanic.White.numerator",
                 "X..Non.Hispanic.White.denominator",
                 "X..Rural.numerator",
                 "X..Rural.denominator")
county_cols_to_numeric <- county_cols[5:length(county_cols)]
county_cols_new_names <- c("state_fips", # new names for first 6 columns needed
                           "state_abbrev", # order counts!
                           "county_name",
                           "county_fips",
                           "firearms.numerator",
                           "firearms.denominator",
                           "injury.denominator",
                           "income",
                           "segregation")
county_cols_old_names <- county_cols[1:9]
names(county_cols_old_names) <- county_cols_new_names

# read CHR&R data, stripping header, stripping state-only FIPS, select columns
# of interest, rename, and create dependent columns.
county0 <- read.csv("data/analytic_data2023_0.csv", header=TRUE)[3:3195,] %>%
  filter(County.FIPS.Code != "000") %>%
  select(all_of(county_cols)) %>%
  mutate(across(all_of(county_cols_to_numeric), as.numeric)) %>%
  rename(all_of(county_cols_old_names)) %>%
  mutate(MHProv = 1000*Mental.Health.Providers.numerator/Mental.Health.Providers.denominator,
         PCP = 1000*Primary.Care.Physicians.numerator/Primary.Care.Physicians.denominator,
         homicide = 100000*Homicides.numerator/Homicides.denominator,
         suicide = 100000*Suicides.numerator/Suicides.denominator,
         uninsured = 100*Uninsured.numerator/Uninsured.denominator,
         unemployed = 100*Unemployment.numerator/Unemployment.denominator,
         grad = 100*High.School.Completion.numerator / High.School.Completion.denominator,
         hisp = 100*X..Hispanic.numerator/X..Hispanic.denominator,
         black = 100*X..Non.Hispanic.Black.numerator/X..Non.Hispanic.Black.denominator,
         white = 100*X..Non.Hispanic.White.numerator/X..Non.Hispanic.White.denominator,
         rural = 100*X..Rural.numerator/X..Rural.denominator,
         urban = ifelse(rural <51, 1,  #where 1 = urban
                        ifelse(rural >=51, 2, NA)), #where 2 = rural
         .keep="unused") %>%
  mutate(censored = ifelse(is.na(firearms.numerator), TRUE, FALSE))

# Note 1: column 'firearms' (ratio) is created not here but in a later script
#         instead, after we perform some corrections (note 2, and part 6).
# Note 2: firearm fatality denominator is provided, but it is not as complete
#         as the injury death denominator.  After we show the injury death
#         denominator is at least as complete, and otherwise equal to, the
#         firearm denominator, we use the injury death denominator instead.
#         It may be preferable to just discard these rows however, at the cost
#         of higher missingness but the benefit of less imputation.
# Now, Confirm the following:
# 1) If there is an injury.denominator value, there is also a
#    firearms.denominator value
# 2) wherever there are both values, they are the same
#    If these are both true, we can use injury.denominator in lieu of
#    firearms.denominator going forward, because injury.denominator is more
#    complete.

county0 <- county0 %>%
  mutate(
    missing_check = ifelse(is.numeric(firearms.denominator),
                           ifelse(is.na(injury.denominator),
                                  TRUE,
                                  FALSE),
                           FALSE ),
    denoms_same = ifelse(!is.na(firearms.denominator) & !is.na(injury.denominator),
                         firearms.denominator == injury.denominator,
                         TRUE) )

sum(county0$missing_check) # expected == 0. i.e., if there is a firearms denominator, there is also an injury denominator
sum(county0$denoms_same==FALSE) # expected == 0. i.e., if there is a number for both denominators, they are the same

county0 <- county0 %>% select(-missing_check, -denoms_same)

sum(is.na(county0$injury.denominator)) # still missing 105 denominators.
# More on that, later in part 6

summary(county0$homicide)
summary(county0$suicide)
summary(county0$uninsured)
summary(county0$unemployed)
summary(county0$grad)
summary(county0$hisp)
summary(county0$black)
summary(county0$white)
summary(county0$rural)
table(county0$urban)
table(county0$censored)

# county0 preparation for joining is complete.

############################################### ###
# 6. Apply CDC WONDER data on 0-death counties ####
############################################### ###

# Lets see if we can complete missings with the CDC WONDER data on firearm
# deaths.  County Health Measures suppresses deaths < 10 at the county level,
# but CDC WONDER will provide counties with 0 deaths.  This can be re-incorporated
# without compromising privacy, IF the numbers match up where available:
# - The 0's in the CDC WONDER set can be put directly into the df.
# - The Population denominator in CDC WONDER can be used to help complete the
# missing injury.denominator values

tryCatch({
          df <- left_join(county0, county_deaths, by=("county_fips" = "County_fips"))
          print("left_join successful")
         },
         error=function(e) { message("left_join error"); print(e)},
         warning=function(e) { message("left_join warning"); print(e)} )
# FAILS because county_deaths is missing some County_fips values that exist in
# county0.  Let's explore what the differences are:

county0 %>%
  select(county_fips, firearms.numerator, injury.denominator) %>%
  filter(!(county_fips %in% county_deaths$County_fips))
county_deaths %>% filter(!(County_fips %in% county0$county_fips))
# notice that:
# -County_fips '02270' has same numerator and denominator as county_fips '02158'
# -County_fips '46113' has same numerator and denominator as county_fips '46102'
# According to
# https://www.cdc.gov/nchs/data/data_acces_files/County-Geography.pdf
# the codes are equivalent but County Health Measures appears to be using the
# values preferred since 2014.  We can manually change county_deaths to allow
# for left_join to work

county_deaths[county_deaths$County_fips=="02270","County_fips"] <- "02158"
county_deaths[county_deaths$County_fips=="46113","County_fips"] <- "46102"

tryCatch({
          df <- left_join(county0, county_deaths, by=c("county_fips" = "County_fips"))
          print("left_join successful")
         },
         error=function(e) { message("left_join error"); print(e)},
         warning=function(e) { message("left_join warning"); print(e)} ) # works now

# Make sure the firearm fatalities are the same (ignoring 0s and NAs),
# then incorporate 0s
check.deaths <- df %>%
  select(county_fips, firearms.numerator, Deaths) %>%
  filter(!is.na(firearms.numerator))
identical(check.deaths$firearms.numerator, check.deaths$Deaths) #TRUE

sum(is.na(df$firearms.numerator)) #871, ~28%
df <- df %>%
  mutate(firearms.numerator = ifelse(!is.na(firearms.numerator),
                                     firearms.numerator, Deaths)) %>%
  select(-Deaths) # doing it this way preserves column order
sum(is.na(df$firearms.numerator)) #818, ~26%

# Now, we can determine if something is censored or uncensored:
df <- df %>% mutate(censored = ifelse(is.na(firearms.numerator), TRUE, FALSE))

# Now make sure the population sizes are the same (ignoring NAs)
check.population <- df %>%
  select(county_fips, injury.denominator, Population) %>%
  filter(!is.na(injury.denominator))

identical(check.population$injury.denominator, check.population$Population) #FALSE
which(check.population$injury.denominator != check.population$Population) #88 in check.population is row 93 in df
check.population[which(check.population$injury.denominator != check.population$Population),] #hmm...
# The County Health Measures denominator is 46345.  The Population size is 36999
# According to the CDC WONDER data a/o 10/2/2023.
# CHM denominator is "denominator is the aggregate annual population over the
# 5-year period."
# Examination of the other CHM dataset columns does not present an obvious
# explanation.
# 10/3/2023 - KC emailed CHR&R maintainers to query this situation
# 10/4/2023 - KC received explanation that the CHR&R value was derived from
# examination of this county_fips code before and after it was divided into
# two other counties and then abolished. This county appears to be an outlier
# situation, and we have opted to use the CHR&R value for consistency.
# "CDC County Geography Changes: 1990-present" demonstrates this county is an
# isolated incident and this discrepancy is does not present again with occult
# discrepancies between masked CHR&R 'injury.denominator' values and CDC WONDER
# 'Population' values.

# prioritize the CHR&R data over CDC WONDER for the single instance of conflict.

df <- df %>%
  mutate(df, injury.denominator = ifelse(!is.na(injury.denominator),
                                         injury.denominator,
                                         Population)) %>%
  select(-Population)

sum(is.na(df$injury.denominator)) #0

################################## ###
# 7. Join, finalize columns, save ####
################################## ###

df$firearms <- (df$firearms.numerator / df$injury.denominator)*100000
summary(df$firearms)

# remaining joins: guns0, pam0, remove fips for Washington DC
# (because Washington DC is missing the primary exposure)
df <- df |>
  # join guns0, pam0, and remove District of Columbia (missing primary exposure)
  left_join(guns0, by=("state_fips" = "state_fips")) |>
  left_join(pam0, by=("county_fips" = "county_fips")) |>
  subset(state_fips !="11")

# Missingness check of final analysis df
summary(df)

# Create a quartile column for gun ownership rate.  Rank comes from RAND.
# This will be used in the modeling portion.
df$qrt <- ifelse(df$rank %in% 1:12, "4",
                  ifelse(df$rank %in% 13:25, "3",
                         ifelse(df$rank %in% 26:37, "2",
                                ifelse(df$rank %in% 38:50,"1", NA))))
df$qrt <- as.factor(df$qrt)
bins <- df %>%
  select(state_abbrev, qrt) %>%
  distinct(state_abbrev, .keep_all = TRUE) %>%
  arrange(qrt)

# eTable 1
bins[1:13,]  # states in the lowest (1st) quartile of gun ownership
bins[14:25,] # second quartile
bins[26:38,] # third quartile
bins[39:50,] # states in the highest (4th) quartile of gun ownership

# Saving as RDS file preserves all class information:
saveRDS(df, file="data/initial_import.Rds")

# The End part 1 #