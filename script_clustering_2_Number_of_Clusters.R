
###########################################################
#                                                         #
#  Clustering.2: Determining the best number of Clusters  #
#                Better PCA plotting                      #
#                                                         #
###########################################################

# https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/

library(tidyverse)
library(factoextra)

# ## input upon original creation of this file:
# annotation = read_rds("result_2019_08_16_clustering_row_annotation")
# # seven_aa_groups = read_rds("result_2019_08_16_clustering_seven_aa_groups")
# # seven_aa_groups_scored = read_rds("result_2019_08_16_clustering_seven_aa_groups_scored")  # matrix
# 
# seven_aa_groups_separated_values = read_rds("result_2019_08_19_seven_aa_groups_separated_values")
# seven_aa_groups_separated_values_scored = read_rds("result_2019_08_19_seven_aa_groups_separated_values_scored")
# 
# scored_dataset = seven_aa_groups_separated_values_scored

# input definition on the 2.9.19:
annotation = row_annotation

scored_dataset = scored

# average silhouette width
fviz_nbclust(x = scored_dataset,
             FUNcluster = hcut,
             method = "silhouette",
             k.max = 40) +
  labs(title = "Silhouette method, k.max = 40") +
  xlab("\nNumber of clusters k") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_text(size = rel(1.3)))

# gap statistics                                          # Ueber Nacht laufen lassen!!!!!!!!!
set.seed(123)
fviz_nbclust(x = scored_dataset,
             FUNcluster = hcut,
             method = "gap_stat",
             nboot = 500,   # number of bootstrap samples
             k.max = 30     # maximum number of clusters to consider
             ) +
  labs(title = "Gap statistics, nboot = 500, k.max = 30")

# total within sum of square
fviz_nbclust(x = scored_dataset,
             FUNcluster = hcut,
             method = "wss",
             k.max = 40) +
  labs(title = "Ellbow method, k.max = 40") +
  xlab("\nNumber of clusters k") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        title = element_text(size = rel(1.3)))



# another try with a different package: NbClust -------------------------------

library(NbClust)

my_NbClust = NbClust(data = scored_dataset,
                     distance = "euclidean",
                     min.nc = 1,
                     max.nc = 30,
                     method = "ward.D2",
                     index = "gap"  # alternativen: "gap", "silhouette
)  

NbClust_plot_input = as.data.frame(my_NbClust$All.CriticalValues) %>% 
  dplyr::rename("critical_values" = "my_NbClust$All.CriticalValues") %>% 
  tibble::rownames_to_column(var = "k")
NbClust_plot_input$k = as.numeric(NbClust_plot_input$k)

NbClust_plot = ggplot(data = NbClust_plot_input, mapping = aes(x = k, y = critical_values)) +
  geom_line(stat = "identity", color = "#3385ff", size = 1.1) +
  labs(title = "NbClust method combination, k.max = 30") +
  xlab("\nNumber of clusters k") +
  ylab("Critical values") +
  theme(axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        title = element_text(size = rel(1.3)),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

NbClust_plot


# PCA copied from script_clustering.R -----------------------------------------

my_pca = prcomp(scored_dataset,
                scale. = F,
                center = T,
                retx = T)

# str(my_pca)
# plot(my_pca)
# plot(my_pca$x[,1], my_pca$x[,2])

# pca_labels = as.matrix(df_Steffi_classes)
# my_pca$x = cbind(my_pca$x, pca_labels)

# pca plot according to http://huboqiang.cn/2016/03/03/RscatterPlotPCA

library(grid)
library(gridExtra)

my_pca_df = as.data.frame(my_pca$x)
my_pca_df$labellls = rownames(my_pca$x)

temp = select(annotation, app_class) %>% 
  rownames_to_column()

my_pca_df = inner_join(my_pca_df, temp, by = c("labellls" = "rowname"))



percentage <- round(my_pca$sdev / sum(my_pca$sdev) * 100, 2) #standard deviation
percentage <- paste( colnames(my_pca_df), "(", paste( as.character(percentage), "%", " )", sep = "") )

p = ggplot(data = my_pca_df, mapping = aes(x = PC1, y = PC2, color = app_class)) +
  geom_point() +
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")) +
  scale_color_manual(values = c("#009900", "#0066ff", "#8c8c8c", "#ff6600")) +
  xlab(percentage[1]) +
  ylab(percentage[2])

p


# plot features which contribute to classification
my_pca_df_r = as.data.frame(my_pca$rotation)
my_pca_df_r$feature = row.names(my_pca_df_r)

contributing_variables = ggplot(data = my_pca_df_r, mapping = aes(x=PC1,y=PC2,label=feature,color=feature )) +
  geom_point() +
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")) +
  geom_text(size = 3)

contributing_variables









