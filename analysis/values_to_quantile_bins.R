# This function takes a list of variants annotated with pathogenicity
# scores as input and converts those annotations to "bins", where each
# bin is based on the observed distribution of the scores. As a
# result, the values in each column become integers, but the number of
# variants remains unchanged.

library(dplyr)
library(magrittr)

x <- vroom::vroom(snakemake@input[["In"]])

if (snakemake@params[["variants_input"]] == TRUE) {
  x <- select(x, -locus, -alleles) %>%
    mutate(variant_count = 1, singleton_count = ifelse(AC == 1, 1, 0))
}

if (!is.null(snakemake@params[["scores_to_include"]])) {
  scores <- c(grep("^sai_", colnames(x), value = TRUE), grep("^ssm_", colnames(x), value = TRUE))
  scores_to_remove <- setdiff(scores, snakemake@params[["scores_to_include"]])
  x <- select(x, -all_of(scores_to_remove))
}

########################################
# TODO: this is just for testing!
is.na(x$sai_xgb) %>% summary()
x <- na.omit(x)
is.na(x$sai_xgb) %>% summary()
########################################

x %>%
  mutate_at(vars(starts_with("sai_"), starts_with("ssm_")), ~ {
    if (any(is.na(.))) {
      stop("NAs found in input scores")
    } else {
      quantiles <- quantile(.,
        probs = snakemake@params[["quantiles"]],
      )
      if (length(snakemake@params[["quantiles"]]) != 16) stop("Wrong quantiles")
      case_when(
        (. >= quantiles[16]) ~ 16,
        (. >= quantiles[15]) ~ 15,
        (. >= quantiles[14]) ~ 14,
        (. >= quantiles[13]) ~ 13,
        (. >= quantiles[12]) ~ 12,
        (. >= quantiles[11]) ~ 11,
        (. >= quantiles[10]) ~ 10,
        (. >= quantiles[9]) ~ 9,
        (. >= quantiles[8]) ~ 8,
        (. >= quantiles[7]) ~ 7,
        (. >= quantiles[6]) ~ 6,
        (. >= quantiles[5]) ~ 5,
        (. >= quantiles[4]) ~ 4,
        (. >= quantiles[3]) ~ 3,
        (. >= quantiles[2]) ~ 2,
        (. >= quantiles[1]) ~ 1,
        TRUE ~ NA_real_
      )
    }
  }) -> x

x %>%
  # All scores should end with "_score", with the only exception being CADD
  mutate_at(vars(ends_with("_score"), starts_with("CADD_")), ~ {
    if (any(is.na(.))) {
      stop("NAs found in quantile bins")
    } else {
      .
    }
  }) -> x

write.table(x, file = snakemake@output[["Out"]], quote = FALSE, row.names = FALSE, sep = "\t")
