# strategy
# each run for a given K is represented by a q_tbl
# a q_tbl has columns name, cluster, prob, and population
# multiple runs are then assembled into a q_multi_tbl
# which has metadata *k* and *run* and *q_tbl* column which is a list
# of q_tbl's
# We can then create a function that plots an q_tbl
# we then need a strategy to combine multiple q_tbl for a given k
# we finally want to align multiple k

# we need a function to read Q files
# import_q <- function(files, indiv_meta)
# where files is a either a vector of names
# and data is a grouped tibble with *name* as the sample name column
# first create a list of q matrices



# modified from https://luisdva.github.io/rstats/model-cluster-plots/
# plotting

library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)

k2plot <-
  ggplot(q_tbl, aes(x = factor(sampleID), y = prob, fill = factor(popGroup))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(loc), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 1)) +
  theme(
    panel.spacing.x = unit(0.1, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_gdocs(guide = FALSE)
k2plot



#### for dapc
q_tbl <- object$posterior %>%
  as_tibble() %>% rename_with(~paste0("Q",.x)) %>%
  # add the pops data for plotting
  mutate(name = rownames(object$posterior),
         group = object$grp) %>%
  tidyr::pivot_longer(cols = starts_with("Q"), names_to = "q", values_to = "prob")

ggplot(q_tbl, aes(x = name, y = prob, fill = q)) +
  geom_col(color = "gray", size = 0.1)+
  facet_grid(~group, switch = "x", scales = "free", space = "free")+
  theme_minimal() + labs(x = "", y = "") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expand_scale(add = 0.7)) +
  theme(
    panel.spacing.x = unit(0.01, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) + guides(fill = "none")
