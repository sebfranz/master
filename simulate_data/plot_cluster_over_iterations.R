library(ggplot2)
library(ggalluvial)
library(reshape2)

d <- cell_cluster_history
rownames(d) <- paste0("cell", 1:nrow(d))
d <- melt(d)
colnames(d) <- c("cell", "iteration", "cluster")
d['cluster'] <- as.factor(d[, 'cluster'])

# Plotting it
# Slow. But keeps track of individual cells
# ggplot(d, aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
#   scale_fill_brewer(type = "qual", palette = "Set2") +
#   geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
#   geom_stratum() +
#   theme(legend.position = "bottom") +
#   ggtitle("Treatment across observation period")

# Doesn't keep track of individual cells
ggplot(d, aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow() +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Treatment across observation period")
