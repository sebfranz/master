# Plot histograms of r2
r2plot <- function(iteration, r2, prev_cell_clust, n_cell_clusters = 2){
  if(T){
    par(mfrow=c(length(unique(prev_cell_clust)),n_cell_clusters))
    for(i_cells_from_cell_cluster in 1:length(unique(prev_cell_clust))){
      for(i_fits_into_cell_cluster in 1:length(unique(prev_cell_clust))){
        print(paste(i_cells_from_cell_cluster, i_fits_into_cell_cluster))
        ind_for_cell_cluster = which(rep(1:n_cell_clusters, n_target_gene_clusters)==i_fits_into_cell_cluster)
        hist(r2[ prev_cell_clust==i_cells_from_cell_cluster, ind_for_cell_cluster], breaks=100, main=paste("Cells from cell cluster", i_cells_from_cell_cluster, "\nfits into cell cluster", i_fits_into_cell_cluster, "with r2:"))
      }
    }
    mtext(paste("Iteration", iteration), side = 3, line = -1, outer = TRUE)
  }
}
