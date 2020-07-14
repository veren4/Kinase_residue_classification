
###########################################################
#                                                         #
#  Clustering.3: Trying out different clustering methods  #
#                                                         #
###########################################################

library(tidyverse)
library(pheatmap)
library(viridis)

save_pheatmap_png = function(x, filepath, width = 600, height = 1000, res = 750){
  png(filepath, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

annotation = read_rds("result_2019_08_16_clustering_row_annotation")
seven_aa_groups_separated_values = read_rds("result_2019_08_19_seven_aa_groups_separated_values")
seven_aa_groups_separated_values_scored = read_rds("result_2019_08_19_seven_aa_groups_separated_values_scored")

###
scored_dataset = seven_aa_groups_separated_values_scored
###

# https://www.r-bloggers.com/k-means-clustering-in-r/

set.seed(42)

my_kmeans_clustering = kmeans(x = seven_aa_groups_separated_values_scored,
                              centers = 4,
                              nstart = 20)

my_kmeans_clustering


# https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/

library(factoextra)



fviz_cluster()












classification_annotation_colors = list(
  app_class = c(key = "#009900",   # green
                potency = "#0066ff",  # blue
                scaffold = "#8c8c8c",    # grey
                selectivity = "#ff6600"),  # orange
  Steffi_manual_class =  c(key = "#009900",   # green
                           potency = "#0066ff",  # blue
                           scaffold = "#8c8c8c",    # grey
                           selectivity = "#ff6600")  # orange
)


my_heatmap = pheatmap(mat = scored_dataset,
                      cluster_rows = T,
                      cluster_cols = F,
                      clustering_distance_rows = "euclidean",    # distance measure
                      clustering_method = "ward.D2",             # agglomeration method of the hierarchical clustering (as in hclust)
                      
                      annotation_row = annotation,
                      # cutree_rows = 4,
                      annotation_colors = classification_annotation_colors,
                      cellwidth = 10,
                      treeheight_row = 100,
                      show_rownames = F,
                      
                      color = viridis(200),
                      angle_col = 45,
                      fontsize = 7,
                      
                      kmeans_k = 4
)

save_pheatmap_png(x = my_heatmap,
                  filepath = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/plots/M.png",
                  width = 700,
                  height = 1200,
                  res = 150)