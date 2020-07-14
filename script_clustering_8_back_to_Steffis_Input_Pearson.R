
####################################################################################
#                                                                                  #
#  Clustering.8: Back to the roots                                                 #
#                Weiterentwicklung von Clustering.4:                               #
#                Input genau wie bei Steffi                                        #
#                Pearson Distanz statt Euklidisch (Winkel statt absolute Distanz)  #
#                                                                                  #
####################################################################################

require(tidyverse)
require(pheatmap)
require(dendextend)
library(viridis)
library(grid)
library(ggbiplot)

# scaling function
scale_this_column_range_01 = function(col){
  my_range = range(col)
  col = col - my_range[1]
  col = col / (my_range[2] - my_range[1])
  return(col)
}

# setwd("C:/Users/vburger/code/partitioned_residue_classification")
# base = read_csv("Z:/users_files/Verena Burger/11_PCA_and_Clustering/Kinase_Residue_Classification_2019-09-29.csv")

# use Classification with Mathias' Alignment:
base = read_csv("Z:/users_files/Verena Burger/14_Mathias_Alignment_as_App_MSA/Kinase_Residue_Classification_2019-10-04.csv")


Steffi_supplementary_tables = read_rds("result_2019_08_27_Seffi_supplementary_tables")


clustering_input = dplyr::select(base,
                                 gene_name,
                                 residue_kinase,
                                 kinome_wide_conservation,
                                 bb,
                                 sc, 
                                 bbsc_total_observations,
                                 targetable_inert,
                                 target_wide_conservation,
                                 median_pkDapp_M,
                                 weighted_median_pkDapp_M_different_aa,
                                 functional_class,
                                 kinase_position,
                                 conservation_overrepresentation_factor,
                                 increased_affinity) %>%
  dplyr::filter(gene_name %in% c("EPHA2", "ABL1", "MELK")) %>%
  dplyr::mutate(seen_sc = if_else(condition = bbsc_total_observations != 0,
                                  true = if_else(condition = sc >= 1,
                                                 true = 1,
                                                 false = 0),
                                  false = 0)) %>% 
  dplyr::mutate(aa_group = if_else(condition = targetable_inert == "aliphatic",
                                   true = "aliphatic",
                                   false = "other")
  ) %>% 
  dplyr::select(-targetable_inert)


# add unique rownames for clustering ------------------------------------------

clustering_input$unique_case_name = NA
for(i in 1:nrow(clustering_input)){
  clustering_input[i, "unique_case_name"] = paste0(i, "_", clustering_input[i, "gene_name"])
}
rm(i)

clustering_input = column_to_rownames(clustering_input, var = "unique_case_name")


# annotation ------------------------------------------------------------------

### from my app

functional_class_myApp = dplyr::select(clustering_input,
                                       gene_name,
                                       functional_class,
                                       conservation_overrepresentation_factor,
                                       increased_affinity,
                                       bbsc_total_observations,
                                       sc,
                                       kinome_wide_conservation,
                                       aa_group
) %>% 
  
  # shiny app classification
  dplyr::rename("automatic_classification" = "functional_class",
  ) %>%
  
  # overrepresentation in ts?
  mutate(overrep_in_ts = if_else(condition = conservation_overrepresentation_factor > 1.5,
                                 true = "yes", false = "no")) %>% 
  select(-conservation_overrepresentation_factor) %>% 
  
  # conveys increased affinity?
  mutate(conv_increased_affinity = if_else(condition = increased_affinity >= 0.5,
                                           true = "yes", false = "no")) %>% 
  select(-increased_affinity) %>% 
  
  # backbone exposed or chemically inert?
  mutate(backbone_exposed = if_else(condition = bbsc_total_observations == 0,
                                    true = "no information",
                                    false = if_else(condition = sc >= 1,
                                                    true = "no backbone",
                                                    false = "backbone exposed"))) %>% 
  select(-bbsc_total_observations, -sc) %>% 
  mutate(exposed_bb_OR_chem_inert = if_else(condition = (backbone_exposed == "backbone exposed" | aa_group == "aliphatic"),
                                            true = "yes", false = "no")) %>% 
  select(-backbone_exposed, -aa_group) %>% 
  
  # kw conservation level >= 50%?
  mutate(kw_cons_over_50 = if_else(condition = kinome_wide_conservation >= 50,
                                   true = "yes", false = "no")) %>% 
  select(-kinome_wide_conservation) %>% 
  
  rownames_to_column()

functional_class_myApp$automatic_classification = factor(functional_class_myApp$automatic_classification)

### from Steffi's manual classification

functional_class_Steffi_manual = dplyr::select(clustering_input, gene_name, kinase_position) %>% 
  tibble::rownames_to_column() %>% 
  dplyr::inner_join(y = Steffi_supplementary_tables, by = c("gene_name" = "gene_name",
                                                            "kinase_position" = "kinase_position")) %>% 
  dplyr::select(rowname, Steffi_manual_class) %>% 
  dplyr::rename("manual_classification" = "Steffi_manual_class")
functional_class_Steffi_manual$manual_classification = factor(functional_class_Steffi_manual$manual_classification)


row_annotation = full_join(functional_class_myApp, functional_class_Steffi_manual,
                           by = c("rowname" = "rowname")) %>% 
  column_to_rownames(var = "rowname") %>% 
  select(manual_classification, automatic_classification, kw_cons_over_50, exposed_bb_OR_chem_inert, conv_increased_affinity, overrep_in_ts, gene_name) %>% 
  dplyr::rename("manual classification" = "manual_classification",
                "automatic classification" = "automatic_classification",
                "6 kw cons over 50" = "kw_cons_over_50",
                "5 exposed bb OR aliphatic" = "exposed_bb_OR_chem_inert",
                "3 conv increased affinity" = "conv_increased_affinity",
                "3 overrep in ts" = "overrep_in_ts",
                "gene name" = "gene_name")

### annotation colors

classification_annotation_colors = list(
  `automatic classification` = c(key = "#009900",   # green
                                 potency = "#0066ff",  # blue
                                 scaffold = "#8c8c8c",    # grey
                                 selectivity = "#ff6600"),  # orange
  `manual classification` =  c(key = "#009900",   # green
                               potency = "#0066ff",  # blue
                               scaffold = "#8c8c8c",    # grey
                               selectivity = "#ff6600"),  # orange
  `gene name` = c(ABL1 = "#FFD700", # gelb
                  EPHA2 = "#FF1493",  # pink
                  MELK = "#663399"),   # purple
  `3 overrep in ts` = c(yes = "green", no = "red"),
  `3 conv increased affinity` = c(yes = "green", no = "red"),
  `5 exposed bb OR aliphatic` = c(yes = "green", no = "red"),
  `6 kw cons over 50` = c(yes = "green", no = "red")
)


# scaling ---------------------------------------------------------------------

scored = dplyr::select(clustering_input,
                       -gene_name, -residue_kinase,
                       -bbsc_total_observations, -bb, -sc,
                       -functional_class, -kinase_position,
                       -conservation_overrepresentation_factor, -increased_affinity) %>% 
  dplyr::rename("affinity same aa" = "median_pkDapp_M",
                "affinity diff aa" = "weighted_median_pkDapp_M_different_aa",
                "exposed backbone" = "seen_sc",
                "kinome wide conservation" = "kinome_wide_conservation",
                "target wide conservation" = "target_wide_conservation") %>% 
  rownames_to_column() %>% 
  dplyr::mutate(`aa group` = if_else(condition = aa_group == "aliphatic",
                                     true = 0.5,
                                     false = 0)) %>% 
  column_to_rownames() %>% 
  dplyr::select(-aa_group)

scored$`kinome wide conservation` = scale_this_column_range_01(scored$`kinome wide conservation`)
scored$`target wide conservation` = scale_this_column_range_01(scored$`target wide conservation`)

scored$`affinity same aa` = scale_this_column_range_01(scored$`affinity same aa`)
scored$`affinity diff aa` = scale_this_column_range_01(scored$`affinity diff aa`)

scored$`exposed backbone` = scored$`exposed backbone`/2

scored = as.matrix(scored)

# plot ------------------------------------------------------------------------

my_heatmap = pheatmap(mat = scored,
                      cluster_rows = T,
                      cluster_cols = F,
                      clustering_distance_rows = "euclidean",    # distance measure: "correlation"/ "euclidean"
                      clustering_method = "ward.D2",             # agglomeration method of the hierarchical clustering (as in hclust)
                      
                      annotation_row = row_annotation,
                      cutree_rows = 7,
                      annotation_colors = classification_annotation_colors,
                      cellwidth = 10,
                      treeheight_row = 100,
                      show_rownames = F,
                      
                      color = viridis(200),
                      fontsize = 7
)



# PCA -------------------------------------------------------------------------

my_pca = prcomp(scored,
                scale. = F,
                center = T,
                retx = T)

group_colors = row_annotation$`automatic classification`


a_gg_plot = ggbiplot(my_pca,
                     ellipse = F,
                     # ellipse.prob = 0.95,
                     choices = c(1,2),
                     groups = group_colors,
                     var.axes = T
) + scale_color_manual(name = "groups", values = c("key"="#009900", # green
                                                   "potency"="#0066ff", # blue
                                                   "scaffold"="#8c8c8c", # grey
                                                   "selectivity"="#ff6600")) + # orange
  theme_light() +
  theme(legend.title = element_blank())
# theme(legend.title = "Automatic Classification")
# theme(legend.title = element_text())
# labs(groups = "Automatic Classification")
# legend.text = element_text()

a_gg_plot

# save_png = function(x, filepath, width = 500, height = 480, res = 750){
#   png(filepath, width = width, height = height, res = res)
#   grid::grid.newpage()
#   grid::grid.draw(x$gtable)
#   dev.off()
# }
# save_png(x = a_gg_plot,
#          filepath = "Z:/users_files/Verena Burger/1_presentations/Abschluss_Praesentation/A.png",
#          width = 1000,
#          height = 100,
#          res = 150)
# dev.off()

# umap ------------------------------------------------------------------------
library(umap)

set.seed(123)
my_umap = umap(d = scored)


source("lib_umap_plot.R")

plot_umap(x = my_umap,
          labels = row_annotation$`automatic classification`,
          # main = "",
          colors = c("#009900", "#0066ff", "#8c8c8c", "#ff6600")
)


# t-SNE -----------------------------------------------------------------------
library(Rtsne)

# scored = read_rds("result_2019_10_10_scored")

tsne_input = scored %>% 
  as.data.frame() %>% 
  rownames_to_column() %>%
  distinct(`kinome wide conservation`,
           `target wide conservation`,
           `affinity same aa`,
           `affinity diff aa`,
           `exposed backbone`,
           `aa group`,
           .keep_all = T) %>% 
  column_to_rownames() %>% 
  as.matrix() %>% 
  normalize_input()

rownames = row_annotation %>% 
  select(`automatic classification`) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% rownames(tsne_input))



set.seed(90)
my_tsne = Rtsne(X = tsne_input,
                check_duplicates = F,
                pca_center = F,
                pca_scale = F,
                normalize = F,
                eta = 100,   # default: 200
                perplexity = 80) # default: 30  # MAX: 99

# Maria adapted:
# perplexity 2, 5, 30, 50, 100
# eta (learning) 10, 20, 50, 100, 1000

tsne_ggplot_input = data.frame(x = my_tsne$Y[,1], y = my_tsne$Y[,2], col = rownames$`automatic classification`, sample = rownames$rowname)



tsne_plot = ggplot(data = tsne_ggplot_input, mapping = aes(x, y, colour = col)) +
  geom_point() +
  scale_color_manual(values = c("key"="#009900", # green
                                "potency"="#0066ff", # blue
                                "scaffold"="#8c8c8c", # grey
                                "selectivity"="#ff6600")) +
  
  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
        legend.key = element_blank())

tsne_plot

library(plotly)
ggplotly(tsne_plot)



# Overlap manual <-> automatic classification ---------------------------------

manual = dplyr::pull(row_annotation, `manual classification`)
manual = manual[!is.na(manual)]
manual = as.data.frame(table(manual))
manual$classification = "manual classification"
manual = dplyr::rename(manual, functional_class = manual)


automatic = dplyr::pull(row_annotation, `automatic classification`)
automatic = as.data.frame(table(automatic))
automatic$classification = "automatic classification"
automatic = dplyr::rename(automatic, functional_class = automatic)


df_comparison = dplyr::bind_rows(manual, automatic, .id = NULL)
df_comparison$Freq = as.numeric(df_comparison$Freq)
df_comparison$classification = base::factor(df_comparison$classification, levels = c("manual classification", "automatic classification"))

sum_manual = dplyr::filter(df_comparison, classification == "manual classification")
sum_manual = sum(sum_manual$Freq)

sum_automatic = dplyr::filter(df_comparison, classification == "automatic classification")
sum_automatic = sum(sum_automatic$Freq)

df_comparison$percentage = 0

for(i in 1:4){
  df_comparison[i, "percentage"] = (df_comparison[i, "Freq"]/sum_manual)*100
}

for(i in 5:8){
  df_comparison[i, "percentage"] = (df_comparison[i, "Freq"]/sum_automatic)*100
}

df_comparison$percentage = round(df_comparison$percentage, digits = 0)

df_comparison$functional_class = factor(df_comparison$functional_class, levels = c("key", "scaffold", "potency", "selectivity"))




classification_comparison_plot = ggplot(data = df_comparison, mapping = aes(x = functional_class, y = Freq, fill = functional_class, alpha = classification)) +
  geom_col(position = "dodge2") +  # geom_col uses stat = "identity" per default
  # facet_grid(rows = vars(classification)) +
  theme_light() +
  scale_fill_manual(values = c("key"="#009900", # green
                               "potency"="#0066ff", # blue
                               "scaffold"="#8c8c8c", # grey
                               "selectivity"="#ff6600"),
                    guide = "none") +
  xlab(element_blank()) +
  ylab("Frequency") +
  geom_text(mapping = aes(label = paste0(percentage, "%")),
            color = "white",
            vjust = 1.8,
            position = position_dodge2(width = 0.9)) +
  scale_alpha_manual(values = c("manual classification" = 0.7,
                                "automatic classification" = 1),
                     guide = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(color = "black", size = rel(1.3)),
        axis.text.y = element_text(color = "black", size = rel(1.3)),
        axis.line = element_line(colour = "black"), axis.ticks = element_line(colour = "black"))
# geom_text(mapping = aes(label = label),
#           color = "white",
#           # vjust = "top",
#           vjust = 3.8,
#           position = position_dodge2(width = 0.9))



classification_comparison_plot

# position = "dodge


### OR: alluvial/ Sanky plot:
library(ggforce)

my_Sankey_data = dplyr::select(row_annotation, `manual classification`, `automatic classification`) %>% 
  as.data.frame()

levels(my_Sankey_data$`manual classification`) = c(levels(my_Sankey_data$`manual classification`), "NA")

my_Sankey_data[is.na(my_Sankey_data$`manual classification`), "manual classification"] = "NA"

my_Sankey_data = my_Sankey_data %>% 
  dplyr::group_by(`manual classification`, `automatic classification`) %>% 
  dplyr::count(name = "value") %>% 
  as.data.frame()



my_Sankey_data = ggforce::gather_set_data(my_Sankey_data, 1:2)

my_Sankey_data$x = factor(my_Sankey_data$x, levels = c("manual classification", "automatic classification"))


TEXT_size = 8

# als svg speichern, damit ich es ggf. noch nachtraeglich aendern kann
# grDevices::pdf(file = "Z:/users_files/Verena Burger/2_thesis/figures/Sankey_classification_comparison.pdf",
#                height = 3.9,
#                width = 3)

ggplot(data = my_Sankey_data, mapping = aes(x, id = id, split = y, value = value)) +
  
  geom_parallel_sets(aes(fill = `manual classification`), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white') +
  
  scale_fill_manual(values = c("key"="#009900", # green
                               "potency"="#0066ff", # blue
                               "scaffold"="#707070", # grey
                               "selectivity"="#ff6600",
                               "NA" = "#D3D3D3"),
                    guide = "none") +
  
  theme_classic(base_size = TEXT_size, base_family = "") +
  theme(
    # plot.title = element_text(size = rel(3)),
    # plot.subtitle = element_text(size = rel(1.5)),
    axis.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(color = "black", size = rel(2)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    # axis.text.x=element_blank(),
    # axis.text.y=element_blank(),
    # axis.ticks.length = unit(1, "cm"),
    # legend.text = element_text(size = 100),
    # legend.key.size = unit(3, "line"),
  )

# dev.off()






