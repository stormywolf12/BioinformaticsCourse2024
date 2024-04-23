install.packages("scLink")
library(scLink)

install.packages('Seurat')
library(Seurat)

expression_matrix <- ReadMtx(
  mtx = "wt_matrix.mtx", features = "wt_features.tsv",
  cells = "wt_barcodes.tsv"
)
#seurat_object <- CreateSeuratObject(counts = expression_matrix)

geneexpressionmatrix <- sclink_norm(expression_matrix, filter.genes = TRUE)
genecoexpressionnetworks <- sclink_net(geneexpressionmatrix, ncores = 1)
correlationmatrix <- sclink_cor(geneexpressionmatrix, ncores = 1)

thresholded_cor <- correlationmatrix
thresholded_cor[abs(correlationmatrix) < 0.7] <- 0
library(igraph)
network <- graph_from_adjacency_matrix(as.matrix(thresholded_cor), weighted = TRUE, mode = "undirected")
plot(network, layout = layout_with_fr)



