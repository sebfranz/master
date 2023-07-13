library(ggplot2)
library(ggalluvial)
library(reshape2)

plot_cluster_history <- function(cell_cluster_history){

  d <- cell_cluster_history
  d <- d[ , colSums(is.na(d))==0]
  d <- melt(d, id.vars="Cell ID")
  colnames(d) <- c("cell", "iteration", "cluster")
  d['cluster'] <- as.factor(d[, 'cluster'])

  # Plotting it
  # Slow. But keeps track of individual cells
  # ggplot(d, aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
  #   scale_fill_brewer(type = "qual", palette = "Set2") +
  #   geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgray") +
  #   geom_stratum() +
  #   theme(legend.position = "bottom") +
  #   ggtitle("Cluster allocation for each iteration")

  # Doesn't keep track of individual cells
  ggplot(d, aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow() +
    geom_stratum() +
    ylab("Cells") +
    xlab("Iteration") +
    labs(fill="Cluster") +
    theme(legend.position = "bottom") +
    ggtitle("Log of cluster allocation")

}
