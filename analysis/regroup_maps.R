#' Recalculates MAPS for a different set of variables
library(magrittr)
library(dplyr)
regroup_maps <- function(maps_table, maps_vars) {
  maps_table <- maps_table %>%
    group_by_at(vars(maps_vars)) %>%
    summarize(
      singleton_count = sum(singleton_count),
      expected_singletons = sum(expected_singletons),
      variant_count = sum(variant_count),
      proportion_singletons = singleton_count / variant_count,
      maps = (singleton_count - expected_singletons) / variant_count,
      maps_sem = sqrt(
        proportion_singletons *
          (1 - proportion_singletons) / variant_count
      ),
      maps_confint_upper = maps + 1.96 * maps_sem,
      maps_confint_lower = maps - 1.96 * maps_sem
    ) %>%
    ungroup()
  return(maps_table)
}
