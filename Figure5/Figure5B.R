suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidygraph)
  library(ggraph)
  library(scales)
})

## -------------------------------
## 0) Input check
## -------------------------------
required_cols <- c("ref", "sample", "variable", "value")
stopifnot(exists("obj"),
          all(required_cols %in% colnames(obj)))

edge_min <- 0.01   # min edge weight
set.seed(123)

## -------------------------------
## 1) Pre-processing
## -------------------------------
obj <- obj %>%
  mutate(group2 = paste0(sample, "~", as.character(variable)))

## -------------------------------
## 2) Aggregate neighbor signal per reference
## -------------------------------
list_dat <- vector("list", length = length(unique(obj$ref)))
names(list_dat) <- sort(unique(obj$ref))

for (i in names(list_dat)) {
  mean_net_score <- obj %>%
    filter(ref == i, variable != i) %>%
    group_by(group2) %>%
    summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(sample  = str_split(group2, "~", simplify = TRUE)[, 1],
           cluster = str_split(group2, "~", simplify = TRUE)[, 2])
  
  dat_net_order <- mean_net_score %>%
    group_by(cluster) %>%
    summarize(value = median(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(ref = i)
  
  list_dat[[i]] <- dat_net_order
}

## -------------------------------
## 3) Combine edges and filter
## -------------------------------
aa <- bind_rows(list_dat)
aa_filt <- aa %>% filter(value > edge_min)

graph_edges <- tibble(
  from   = aa_filt$ref,
  to     = aa_filt$cluster,
  weight = aa_filt$value
)

## -------------------------------
## 4) Build network
## -------------------------------
graph <- tbl_graph(edges = graph_edges, directed = FALSE)

graph <- graph %>%
  activate(nodes) %>%
  mutate(deg = centrality_degree())

## Colors based on node names
node_names <- graph %>% activate(nodes) %>% pull(name)
pal <- setNames(hue_pal()(length(node_names)), node_names)

## -------------------------------
## 5) Plot
## -------------------------------
p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "grey70") +
  geom_node_point(aes(size = deg, color = name), show.legend = FALSE) +
  geom_node_text(aes(label = name),
                 vjust = 0.5, color = "black") +
  scale_size(range = c(2, 18)) +
  scale_edge_width(limits = c(edge_min, max(graph_edges$weight, na.rm = TRUE)),
                   range = c(0.2, 2)) +
  scale_color_manual(values = pal) +
  theme_void()

print(p)
