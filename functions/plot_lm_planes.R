library(scatterplot3d)
library(fields)  # For colors

plot_lm_planes <- function(dat, current_cell_cluster_allocation, ind_reggenes, ind_targetgenes, n_cell_clusters, title){
  data_for_plot <- data.frame(matrix(0, nrow = length(current_cell_cluster_allocation)*length(ind_reggenes), ncol = 4))
  colnames(data_for_plot) <- c('x','y','z','color')
  color <- 0
  max_row <- 0
  previous_max_row <- max_row +1
  for(cell_cluster in 1:n_cell_clusters){
    for(current_target_gene in ind_targetgenes){
      color <- color + 1

      cluster_rows <- which(current_cell_cluster_allocation == cell_cluster)
      previous_max_row <- max_row + 1
      max_row <- max_row + length(cluster_rows)
      ind <- previous_max_row:max_row
      # print(paste(previous_max_row, max_row, color))
      data_for_plot[ind,'x'] <- dat[cluster_rows, ind_reggenes[1]]  # /max(dat[cluster_rows, ind_reggenes[1]])
      data_for_plot[ind,'y'] <- dat[cluster_rows, ind_reggenes[2]]  # /max(dat[cluster_rows, ind_reggenes[2]])
      data_for_plot[ind,'z'] <- dat[cluster_rows, current_target_gene]  # /max(dat[cluster_rows, current_target_gene])
      data_for_plot[ind,'color'] <- color
    }
  }

  # print(data.frame(table(data_for_plot[,'color'])))
  n_unique_colors <- length(unique(data_for_plot[,'color']))
  colorSet <- tim.colors(n_unique_colors)
  colors <- colorSet[data_for_plot[,'color']]
  colors <- adjustcolor(colors, alpha.f = 0.22)
  p <- scatterplot3d(x=data_for_plot[,'x'],
                     y=data_for_plot[,'y'],
                     z=data_for_plot[,'z'],
                     color=colors,
                     pch = 19,
                     xlab = "Regulator gene 1",
                     ylab = "Regulator gene 2",
                     zlab = "Target gene",
                     main = title)


  color <- 0
  for(cell_cluster in 1:n_cell_clusters){
    for(current_target_gene in 1:length(ind_targetgenes)){
      color <- color + 1
      p$plane3d(models[[cell_cluster]]$coefficients[,current_target_gene], col = colorSet[color])
    }
  }
  legend(p$xyz.convert(-0.5, -0.5, -0.5),
         col= c(colorSet[1],colorSet[2], colorSet[3], colorSet[4]),
         bg="white", lty=c(1,1), lwd=2, yjust=0,
         legend = c("C1T1", "C1T2", "C2T1", "C2T2"), cex = 1.1)




}

