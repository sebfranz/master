if (!require(testthat)) install.packages('testthat')
library(testthat)
context("Testing data generation")
filename <- "./simulate_data/functions/generate_dummy_data_for_scregclust.R"
if(!file.exists(filename)){
  filename <- "generate_dummy_data_for_scregclust.R"
}
source(filename);

test_that("Check if list is returned", {

  res <- generate_dummy_data_for_scregclust()

  # Test that the result is numeric
  expect_that( is.list(res), equals(TRUE) )

  # Test that the result is the correct value
  expect_that( length(res), equals(5) )
})


# Generate many test inputs ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


n_target_genes <- c(2,10)
n_regulator_genes <- c(2, 50)
n_cells <- c(100,1000)
n_target_gene_clusters <- c(2, 3)
regulator_mean <- seq(0.1,0.5,0.1)


test_data <- expand.grid(n_target_genes, n_regulator_genes, n_cells, n_target_gene_clusters, regulator_mean)
colnames(test_data) <- c("n_target_genes", "n_regulator_genes", "n_cells", "n_target_gene_clusters", "regulator_mean")

n_target_genes <- as.vector(test_data["n_target_genes"][[1]])
n_regulator_genes <- as.vector(test_data["n_regulator_genes"][[1]])
n_cells <- as.vector(test_data["n_cells"][[1]])
n_target_gene_clusters <- as.vector(test_data["n_target_gene_clusters"][[1]])
regulator_mean <- as.vector(test_data["regulator_mean"][[1]])

for(i in 1:length(regulator_mean)){
  coefficient_mean <- seq(1, 100, length.out= n_target_gene_clusters[i])
  description = paste("\n n_target_genes", n_target_genes[i], "\n",
                      "n_regulator_genes", n_regulator_genes[i], "\n",
                      "n_cells", n_cells[i], "\n",
                      "n_target_gene_clusters", n_target_gene_clusters[i], "\n",
                      "regulator_mean", regulator_mean[i], "\n",
                      "coefficient_mean", paste(coefficient_mean, collapse = " "))

  test_that(desc = description, {



    res <- generate_dummy_data_for_scregclust(n_target_genes = n_target_genes[i],              #number of target genes
                                              n_regulator_genes = n_regulator_genes[i],               #number of regulator genes
                                              n_cells  = n_cells[i],             #number of cells
                                              n_target_gene_clusters  = n_target_gene_clusters[i],               #Number of target gene clusters
                                              regulator_mean   = regulator_mean[i], #mean expression of regulator genes
                                              coefficient_mean = coefficient_mean #mean coefficients in true  model, length n_target_gene_clusters
    )

    # Test that the result is numeric
    expect_that( is.list(res), equals(TRUE) )

    # Test that the result is the correct value
    expect_that( length(res), equals(5) )

  })

}

