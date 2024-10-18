#' This is primarily extraction of code from the Admiral vignette for creating
#' ADSL: https://pharmaverse.github.io/admiral/articles/adsl.html

library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)

# We can read XPT files
dm_xpt <- haven::read_xpt("./dm.xpt")

# Or the new Dataset JSON format! 
library(datasetjson)
dm_json <- read_dataset_json("./dm.json")

# Use some test data for the demo
data("dm")
data("ds")
data("ex")
data("ae")
data("lb")

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)

adsl <- dm %>%
  select(-DOMAIN) %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM)


# impute start and end time of exposure to first and last respectively,
# do not impute date
ex_ext <- ex %>%
  # Derive EXSTDTM and EXSTTMF
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  # Derive EXENDTM and EXENTMF
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

# Derive TRTSDTM and TRTSTMF
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Derive TRTEDTM and TRTETMF
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )

# convert datetime variables to dates
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

# And now calculate durations
adsl <- adsl %>%
  derive_var_trtdurd()

# More date variables from DS
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )

# Disposition status
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )

# Reason for discontinuation
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      DSDECOD %notin% c("SCREEN FAILURE", "COMPLETED", NA)
  ) %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD %in% "OTHER"
  )

# Cause of Death
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS)
  )


# Grouping Variables
format_agegr1 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 64) ~ "18-64",
    var_input > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}

format_region1 <- function(var_input) {
  case_when(
    var_input %in% c("CAN", "USA") ~ "North America",
    !is.na(var_input) ~ "Rest of the World",
    TRUE ~ "Missing"
  )
}

adsl <- adsl %>%
  mutate(
    AGEGR1 = format_agegr1(AGE),
    REGION1 = format_region1(COUNTRY)
  )


# Population flags
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  )


# Get a template program to start from:
admiral::use_ad_template("ADSL")




# Load some pre-saved metadata, based on xportr data
var_spec <- readRDS("var_spec.Rds")
dataset_spec <- readRDS("dataset_spec.Rds")

# Write data using xportr
library(xportr)
xportr_options(
  xportr.df_domain_name = "Dataset",
  xportr.df_label = "Description",
  xportr.variable_name = "variable",
  xportr.label = "Label",
  xportr.type_name = "type",
  xportr.format_name = "Format",
  xportr.length = "Length",
  xportr.order_name = "Order"
)

adsl %>%
  xportr_type(var_spec, "ADSL", "message") %>%
  xportr_length(var_spec, "ADSL", "message") %>%
  xportr_label(var_spec, "ADSL", "message") %>%
  xportr_order(var_spec, "ADSL", "message") %>%
  xportr_format(var_spec, "ADSL") %>%
  xportr_df_label(dataset_spec, "ADSL") %>%
  xportr_write("adsl.xpt")

# Or export using datasetjson
adsl_items <- readRDS("./adsl_items.Rds") # Load prepared metadata
ds_json <- dataset_json(adsl, "IG.ADSL", "ADSL", "Subject-Level Analysis Dataset", adsl_items)
write_dataset_json(ds_json, file="adsl.json")
