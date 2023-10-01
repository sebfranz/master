# library(ggplot2)
# library(ggalluvial)
# library(reshape2)

if (!require(aricode)) install.packages('aricode')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(dplyr)) install.packages('dplyr')
if (!require(ggalluvial)) install.packages('ggalluvial')
if (!require(reshape2)) install.packages('reshape2')
if (!require(ggsankey)){
  if (!require(remotes)) install.packages('remotes')
  library(remotes)
  remotes::install_github("davidsjoberg/ggsankey")
}
library(ggsankey)
library(ggalluvial)
library(ggplot2)
library(reshape2)
library(dplyr)
library(aricode)  # To calculate rand index

plot_cluster_history <- function(cell_cluster_history){

  d <- cell_cluster_history
  d <- d[ , colSums(is.na(d))==0]  # Remove NA data

  new_colnames <- colnames(d)

  rand_ind <- vector(length = (ncol(cell_cluster_history)-1))
  for(i in 2:ncol(cell_cluster_history)){
    rand_ind[i-1] <- round(RI(cell_cluster_history[,2],
                   cell_cluster_history[,i]), 2)
    new_colnames[i] <- paste0(new_colnames[i], "\nRI:", rand_ind[i-1])
  }



  colnames(d) <- new_colnames
  #
  # # This is some toy data to play around with the sankey plot
  # a <- matrix(c(1,10,11,4,2,10,3,8,2,10,11,12), nrow=3, ncol=4, byrow = T)
  # df <- data.frame(a)
  # colnames(df) <- c("a","b","c","d")
  # DF <- df
  # for(icol in 2:(ncol(df)-1)){
  #   for(irow in 1:nrow(df)){
  #     DF[irow, icol] <- paste0(df[irow,icol], "f", df[irow,icol-1], "t", df[irow,icol+1])
  #   }
  # }
  #
  # df <- DF
  # df <- make_long(.df=df, colnames(df))
  # ggplot(df, aes(x = x,
  #                next_x = next_x,
  #                node = node,
  #                next_node = next_node,
  #                fill = factor(node),
  #                label = node)) +
  #   geom_sankey(flow.alpha = 0.5, node.color = 1) +
  #   geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  #   scale_fill_viridis_d() +
  #   theme_sankey(base_size = 16) +
  #   guides(fill = guide_legend(title = "Cluster")) +
  #   labs(x = "") +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   ggtitle("Sankey diagram tracking cluster allocation in each iteration")
  #
  # ggsave('plot.png', dpi = 300, height = 12, width = 20, unit = 'in')

  # d <- df
  # d['Cell ID'] <- 1:3
  d <- melt(d, id.vars="Cell ID")
  colnames(d) <- c("cell", "iteration", "cluster")
  d['cluster'] <- as.factor(d[, 'cluster'])


  # Plotting it
  # Slow. But keeps track of individual cells
  p <- ggplot(d, aes(x = iteration,
                stratum = cluster,
                alluvium = cell,
                fill = cluster,
                label = cluster)) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_flow(stat = "alluvium", lode.guidance = "frontback") +
    geom_stratum(alpha=0.5) +
    geom_text(stat = "stratum", size = 3) +
    theme(legend.position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("Alluvial diagram tracking cluster allocation in each iteration")

  # ggsave('cell_cluster_history.png', plot=p, dpi = 300, height = 6, width = 12, unit = 'in')
  # Doesn't keep track of individual cells
  # ggplot(d, aes(x = iteration, stratum = cluster, alluvium = cell, fill = cluster, label = cluster)) +
  #   scale_fill_brewer(type = "qual", palette = "Set2") +
  #   geom_flow() +
  #   geom_stratum() +
  #   ylab("Cells") +
  #   xlab("Iteration") +
  #   labs(fill="Cluster") +
  #   theme(legend.position = "bottom") +
  #   ggtitle(paste0("Log of cluster allocation\nRand index of true vs final: ",round(rand_ind,2)))




}
