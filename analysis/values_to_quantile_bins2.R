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

create_top_half_mask <- function(len, max_groups = NULL) {
  mask <- c()
  iterations <- log2(len)
  for (i in 0:iterations) {
    mask <- c(mask, rep(i + 1, len / (2^i) - floor(len / (2^(i + 1)))))
  }
  if (!is.null(max_groups)) mask[mask > max_groups] <- max_groups
  return(mask)
}

########################################
# TODO: this is just for testing!
is.na(x$sai_xgb) %>% summary()
is.na(x$ssm_anr) %>% summary()
x <- na.omit(x)
is.na(x$sai_xgb) %>% summary()
is.na(x$ssm_anr) %>% summary()
########################################

x %>%
  mutate_at(vars(starts_with("sai_"), starts_with("ssm_")), ~ {
    if (any(is.na(.))) {
      stop("NAs found in input scores")
    } else {
      # Record the original ordering
      indices <- order(.)
      # Re-order the mask based on the original ordering
      create_top_half_mask(length(indices), snakemake@params[["max_group"]])[order(indices)]
    }
  }) -> x

write.table(x, file = snakemake@output[["Out"]], quote = FALSE, row.names = FALSE, sep = "\t")
