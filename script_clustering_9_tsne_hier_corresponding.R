
###
### Do the clusters in tsne correspond to HC branches?
###

library(tidyverse)


# get the t-SNE groups --------------------------------------------------------

tsne_input = read_rds("result_2019_10_07_tsne_input")
my_tsne = read_rds("result_2019_10_07_my_tsne")   # the order of rows remains the same

### t-SNE plot
rownames = read_rds("result_2019_10_11_rownames")
tsne_ggplot_input = data.frame(x = my_tsne$Y[,1], y = my_tsne$Y[,2], col = rownames$`automatic classification`, sample = rownames$rowname)

tsne_plot = ggplot(data = tsne_ggplot_input, mapping = aes(x, y, colour = col)) +
  geom_point() +
  scale_color_manual(values = c("key"="#009900", # green
                                "potency"="#0066ff", # blue
                                "scaffold"="#8c8c8c", # grey
                                "selectivity"="#ff6600")) +
  
  theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),
        axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"), plot.margin=unit(c(1,1,1,1),"line"),
        legend.key = element_blank(),
        panel.grid = element_line(color = "black")) +
  scale_y_continuous(breaks = seq(-30, 30, by = 10),
                     minor_breaks = seq(-30, 30, by = 5)) +
  scale_x_continuous(breaks = seq(-40, 40, by = 10),
                     minor_breaks = seq(-40, 40, by = 5))

tsne_plot
###


tsne_coordinates = my_tsne$Y %>% 
  as.data.frame() %>% 
  dplyr::rename(x = "V1", y = "V2")

tsne_cluster_identifier = tsne_coordinates %>% 
  dplyr::mutate(cluster_1 = if_else(condition = (x > -14 & x < -5 &  y > 10 & y < 20),
                                    true = 1, false = 0),
                cluster_2 = if_else(condition = (x > -30 & x < -15 & y > 7 & y < 16),
                                    true = 1, false = 0),
                cluster_3 = if_else(condition = (x > -37 & x < -29 & y > -5 & y < 0),
                                    true = 1, false = 0),
                cluster_4 = if_else(condition = (x > -27 & x < -15 & y > -24 & y < -8),
                                    true = 1, false = 0)
                )

rownames(tsne_cluster_identifier) = rownames(tsne_input)

tsne_cluster_1 = tsne_cluster_identifier %>% rownames_to_column() %>% dplyr::filter(cluster_1 == 1) %>% dplyr::pull(rowname)
tsne_cluster_2 = tsne_cluster_identifier %>% rownames_to_column() %>% dplyr::filter(cluster_2 == 1) %>% dplyr::pull(rowname)
tsne_cluster_3 = tsne_cluster_identifier %>% rownames_to_column() %>% dplyr::filter(cluster_3 == 1) %>% dplyr::pull(rowname)
tsne_cluster_4 = tsne_cluster_identifier %>% rownames_to_column() %>% dplyr::filter(cluster_4 == 1) %>% dplyr::pull(rowname)

rm(my_tsne, tsne_cluster_identifier, tsne_coordinates, tsne_input, rownames, tsne_plot)


# get the hierarchical clustering groups --------------------------------------

scored = read_rds("result_2019_10_07_scored")

my_hclust = hclust(dist(scored,
                              method = "euclidean"),
                         method = "ward.D2")
hclust_groups  = stats::cutree(tree = my_hclust, k = 4)

hclust_clusters = scored %>%
  as.data.frame()

hclust_clusters$cluster = hclust_groups

hclust_clusters = hclust_clusters %>% 
  dplyr::select(cluster) %>% 
  rownames_to_column()
                                                                                        # Which cluster is which:
hc_cluster_1 = hclust_clusters %>% dplyr::filter(cluster == 1)# %>% dplyr::pull(rowname) # 1647
hc_cluster_2 = hclust_clusters %>% dplyr::filter(cluster == 2)# %>% dplyr::pull(rowname) # 1674
hc_cluster_3 = hclust_clusters %>% dplyr::filter(cluster == 3)# %>% dplyr::pull(rowname) # 791            D
hc_cluster_4 = hclust_clusters %>% dplyr::filter(cluster == 4)# %>% dplyr::pull(rowname) # 256            C

# From C and D, separate blue and grey:

annotation = read_rds("result_2019_10_07_row_annotation") %>% 
  dplyr::select(`automatic classification`) %>% 
  rownames_to_column()

hc_cluster3_potency = hc_cluster_3 %>% 
  left_join(annotation, by = c("rowname" = "rowname")) %>% 
  dplyr::filter(`automatic classification` == "potency") %>% 
  pull(rowname)

hc_cluster3_scaffold = hc_cluster_3 %>% 
  left_join(annotation, by = c("rowname" = "rowname")) %>% 
  dplyr::filter(`automatic classification` == "scaffold") %>% 
  pull(rowname)

hc_cluster4_potency = hc_cluster_4 %>% 
  left_join(annotation, by = c("rowname" = "rowname")) %>% 
  dplyr::filter(`automatic classification` == "potency") %>% 
  pull(rowname)

hc_cluster4_scaffold = hc_cluster_4 %>% 
  left_join(annotation, by = c("rowname" = "rowname")) %>% 
  dplyr::filter(`automatic classification` == "scaffold") %>% 
  pull(rowname)

rm(annotation, hc_cluster_1, hc_cluster_2, hc_cluster_3, hc_cluster_4,
   hclust_groups, scored, my_hclust)



# die beiden Clusterings vergleichen ------------------------------------------

# hc 3 p <-> tsne 4B
length(hc_cluster3_potency) # 414
length(which(hc_cluster3_potency %in% tsne_cluster_4)) # 262
length(which(tsne_cluster_4 %in% hc_cluster3_potency)) # 262

# hc 3 p <-> tsne 3A
tsne_cluster_3 %in% hc_cluster3_potency # F

# hc 4 p <-> tsne 4B
length(hc_cluster4_potency) # 177
length(which(hc_cluster4_potency %in% tsne_cluster_4)) # 0
length(which(tsne_cluster_4 %in% hc_cluster4_potency)) # 0

# hc 3 sc <-> tsne 2D
tsne_cluster_2 %in% hc_cluster3_scaffold # T


# hc 4 sc <-> tsne 1
length(hc_cluster4_scaffold) # 68
length(which(tsne_cluster_1 %in% hc_cluster4_scaffold)) # 0
length(which(hc_cluster4_scaffold %in% tsne_cluster_1)) # 0

# hc 4 sc <-> tsne 2D
tsne_cluster_2 %in% hc_cluster4_scaffold # F

rm(hc_cluster3_potency, hc_cluster3_scaffold, hc_cluster4_potency, hc_cluster4_scaffold,
   tsne_cluster_1, tsne_cluster_2, tsne_cluster_3, tsne_cluster_4)


# Advanced plot ---------------------------------------------------------------

advanced_tsne_plot_input = dplyr::inner_join(tsne_ggplot_input, hclust_clusters, by = c("sample" = "rowname")) %>% 
  dplyr::rename("hc_cluster" = "cluster")

tsne_plot = ggplot(data = advanced_tsne_plot_input, mapping = aes(x, y, colour = col)) +
  geom_point() +
  scale_color_manual(values = c("key"="#009900", # green
                                "potency"="#0066ff", # blue
                                "scaffold"="#8c8c8c", # grey
                                "selectivity"="#ff6600")) +
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
        legend.key = element_blank(), legend.title = element_blank()) +
  
  geom_point(data = subset(advanced_tsne_plot_input, hc_cluster == 4),
             colour = "black", size = 2.1, shape = 21)
  

tsne_plot




