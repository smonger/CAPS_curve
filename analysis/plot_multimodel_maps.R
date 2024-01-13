library(dplyr)
library(ggplot2)

# TODO: there are bugs here!

# TODO: make a parameter?
source("regroup_maps.R")

df <- data.frame()

for (i in 1:length(snakemake@input[["maps_dataframes"]])) {
  x <- read.table(snakemake@input[["maps_dataframes"]][i], header = TRUE)
  # Add a column with the name of the varible from which the "Top
  # N%" groups were calculated
  x <- mutate(x, model = rep(colnames(x)[1], dim(x)[1]))
  colnames(x)[1] <- "group"
  df <- rbind(df, x)
}

if (!is.null(snakemake@params[["groups"]])) df <- filter(df, group %in% snakemake@params[["groups"]])

# Rename models
if (!is.null(snakemake@params[["model_labels"]]) | !is.null(snakemake@params[["new_model_labels"]])) {
  if (!is.null(snakemake@params[["model_labels"]]) & !is.null(snakemake@params[["new_model_labels"]])) {
    df[["variable"]] <- factor(df[["variable"]],
      levels = snakemake@params[["model_labels"]],
      labels = snakemake@params[["new_model_labels"]]
    )
  } else {
    stop("Either original or new labels were not provided")
  }
}

maps_table <- regroup_maps(df, c("group", "variable"))
pdf(snakemake@output[["plot"]])
ggplot(maps_table) +
  aes(
    # TODO: unstable?
    x = reorder(group, maps),
    y = maps,
    color = variable
  ) +
  {
    if (snakemake@params[["separated"]]) {
      geom_errorbar(
        aes(ymin = maps_confint_lower, ymax = maps_confint_upper),
        width = 0,
        position = position_dodge(width = 0.8)
      )
    } else {
      geom_errorbar(
        aes(ymin = maps_confint_lower, ymax = maps_confint_upper),
        width = 0
      )
    }
  } +
  {
    if (snakemake@params[["separated"]]) geom_point(size = 4, position = position_dodge(width = 0.8)) else geom_point(size = 4)
  } +
  # ylim(
  #   ifelse(is.null(snakemake@params[["ylim_min"]]),
  #     min(maps_table$maps_confint_lower),
  #     snakemake@params[["ylim_min"]]
  #   ),
  #   ifelse(is.null(snakemake@params[["ylim_max"]]),
  #     max(maps_table$maps_confint_upper),
  #     snakemake@params[["ylim_max"]]
  #   )
  # ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      vjust = snakemake@params[["xlab_vjust"]],
      hjust = snakemake@params[["xlab_hjust"]],
      angle = ifelse(is.null(snakemake@params[["xlab_angle"]]), 0, snakemake@params[["xlab_angle"]])
    ),
    text = element_text(size = 30),
    plot.margin = margin(0, 5.5, 0, 5.5)
  ) +
  scale_color_manual(snakemake@params[["legend_title"]],
    values = c("red", "green", "blue", "orange", "magenta")
  ) +
  ylab(snakemake@params[["ylab"]]) +
  xlab(snakemake@params[["xlab"]]) +
  ggtitle(snakemake@params[["title"]]) +
  {
    if (!is.null(snakemake@params[["missense_level"]])) geom_hline(aes(yintercept = snakemake@params[["missense_level"]]), linetype = "dashed")
  } +
  {
    if (!is.null(snakemake@params[["splicing_level"]])) geom_hline(aes(yintercept = snakemake@params[["splicing_level"]]), linetype = "dashed")
  } +
    #TODO: the code below does nothing
  {
    if (!is.null(snakemake@params[["nonsense_level"]])) geom_hline(aes(yintercept = snakemake@params[["nonsense_level"]]), linetype = "dashed")
  }
dev.off()
