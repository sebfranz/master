library(scatterplot3d)
library(fields)  # For colors

plot_lm_planes <- function(dat, current_cell_cluster_allocation, ind_reggenes, ind_targetgenes, n_cell_clusters, title){

  library(fields)  # For colors
  par(mfrow=c(1,2))
  colorSet <- two.colors(2)
  for(cell_cluster in 1:n_cell_clusters){
    cluster_rows <- which(current_cell_cluster_allocation == cell_cluster)
    colors <- colorSet[true_cell_clust_sample[cluster_rows]]
    colors = alpha(colors, 0.09)
    x <- dat[cluster_rows, ind_reggenes[1]]
    y <- dat[cluster_rows, ind_reggenes[2]]
    p <- plot(x, y, col=colors, pch=18, xlab='Regulator gene 1', ylab='Regulator gene 2')
    legend(x=min(x),
           col= colorSet,
           bg="white", lty=c(1,1), lwd=2, yjust=0,
           legend = c("True cluster 1", "True cluster 2"), cex = 1.1)
    plot(dat[cluster_rows, ind_targetgenes[1]], dat[cluster_rows, ind_targetgenes[2]], col=colors, pch=18, xlab='Target gene 1', ylab='Target gene 2')
    mtext(paste("Cell cluster", cell_cluster), side = 3, line = - 2, outer = TRUE)
  }


}

