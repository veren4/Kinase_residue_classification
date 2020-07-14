
###########################################################
#                                                         #
#  Clustering.1: Creating the clustering Input            #
#                Hierarchical Clustering + Plotting       #
#                                                         #
###########################################################


# setwd("C:/Users/vburger/code/partitioned_residue_classification")
require(tidyverse)
require(pheatmap)
require(dendextend)
library(viridis)
library(grid)
library(ggbiplot)
library(umap)

######### functions #################
scale_this_column_range_01 = function(col){
  
  my_range = range(col)  # lowest and highest value of the column
  
  col = col - my_range[1]    # so dass der niedrigste Wert bei 0 liegt
  
  col = col / (my_range[2] - my_range[1])  # jeden Wert durch die originale "Spannweite" teilen
  # man bekommt quasi einen "Prozentwert", wo der Wert liegt. => das muss
  # also immer zwischen 0 und 1 liegen.
  return(col)
  
}

save_pheatmap_png = function(x, filepath, width = 600, height = 1000, res = 750){
  png(filepath, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#####################################

base = read_csv("Z:/users_files/Verena Burger/11_PCA_and_Clustering/classified_Kinobead_residues.csv")

set.seed(42)

base_subset = select(base,
                     gene_name,
                     residue_kinase,
                     kinome_wide_conservation,
                     bb,
                     sc, 
                     bbsc_total_observations,
                     targetable_inert,
                     conservation_overrepresentation_factor, 
                     increased_affinity,
                     functional_class,
                     kinase_position) %>%
  filter(gene_name %in% c("EPHA2", "ABL1", "MELK")) %>%
  sample_n(size = 300)


# reference data from Steffi --------------------------------------------------
Steffi_supplementary_tables = tribble(~gene_name, ~kinase_position, ~Steffi_manual_class,
                                      "ABL1",     248,              "scaffold",
                                      "ABL1",     249,              "scaffold",
                                      "ABL1",     250,              "scaffold",
                                      "ABL1",     251,              "scaffold",
                                      "ABL1",     252,              "selectivity",
                                      "ABL1",     253,              "selectivity",
                                      "ABL1",     256,              "scaffold",
                                      "ABL1",     269,              "scaffold",
                                      "ABL1",     270,	             "scaffold",
                                      "ABL1",     271,	             "potency",
                                      "ABL1",     282,              "selectivity",
                                      "ABL1",     285,              "scaffold",
                                      "ABL1",     286,              "potency",
                                      "ABL1",     289,              "scaffold",
                                      "ABL1",     290,	             "key",
                                      "ABL1",     293,              "scaffold",
                                      "ABL1",     298,              "scaffold",
                                      "ABL1",     299,	             "scaffold",
                                      "ABL1",     313,	             "key",
                                      "ABL1",     315,	             "key",
                                      "ABL1",     316,	             "scaffold",
                                      "ABL1",     317,	             "selectivity",
                                      "ABL1",     318,              "key",
                                      "ABL1",     319,	             "scaffold",
                                      "ABL1",     320,              "key",
                                      "ABL1",     321,	             "scaffold",
                                      "ABL1",     322,              "selectivity",
                                      "ABL1",     325,              "selectivity",
                                      "ABL1",     329,	             "selectivity",
                                      "ABL1",     354,	             "scaffold",
                                      "ABL1",     359,              "key",
                                      "ABL1",     360,	             "scaffold",
                                      "ABL1",     361,	             "potency",
                                      "ABL1",     362,	             "potency",
                                      "ABL1",     367,              "scaffold",
                                      "ABL1",     368,              "potency",
                                      "ABL1",     370,	             "scaffold",
                                      "ABL1",     379,	             "scaffold",
                                      "ABL1",     380,	             "selectivity",                          # Steffi: reverse selectivity
                                      "ABL1",     381,	             "potency",
                                      "ABL1",     382,              "potency",
                                      "ABL1",     383,	             "scaffold",
                                      
                                      "EPHA2",    617,	             "selectivity",
                                      "EPHA2",    619,	             "scaffold",
                                      "EPHA2",    620,	"scaffold",
                                      "EPHA2",    621,	"scaffold",
                                      "EPHA2",    622,	"scaffold",
                                      "EPHA2",    623,	"selectivity",
                                      "EPHA2",    627,	"scaffold",
                                      "EPHA2",    644,	"scaffold",
                                      "EPHA2",    646,	"potency",
                                      "EPHA2",    663,	"potency",
                                      "EPHA2",    666,	"scaffold",
                                      "EPHA2",    667,	"key",
                                      "EPHA2",    670,	"selectivity",
                                      "EPHA2",    675,	"scaffold",
                                      "EPHA2",    676,	"scaffold",
                                      "EPHA2",    690,	"key",
                                      "EPHA2",    692,	"key",
                                      "EPHA2",    693,	"scaffold",
                                      "EPHA2",    694,	"selectivity",
                                      "EPHA2",    695,	"key",
                                      "EPHA2",    696,	"scaffold",
                                      "EPHA2",    697,	"key",
                                      "EPHA2",    698,	"scaffold",
                                      "EPHA2",    699,	"selectivity",   # Steffi: reverse selectivity
                                      "EPHA2",    702,	"selectivity",
                                      "EPHA2",    706,	"selectivity",
                                      "EPHA2",    730,	"scaffold",
                                      "EPHA2",    735,	"selectivity",
                                      "EPHA2",    737,	"potency",
                                      "EPHA2",    743,	"scaffold",
                                      "EPHA2",    744,	"potency",
                                      "EPHA2",    746,	"scaffold",
                                      "EPHA2",    755,	"scaffold",
                                      "EPHA2",    756,	"selectivity",
                                      "EPHA2",    757,	"potency",
                                      "EPHA2",    758,	"potency",
                                      
                                      "MELK",     15,	"selectivity",
                                      "MELK",     17, "scaffold",
                                      "MELK",     18,	"scaffold",
                                      "MELK",     19,	"scaffold",
                                      "MELK",     20,	"scaffold",
                                      "MELK",     22,	"scaffold",
                                      "MELK",     23,	"scaffold",
                                      "MELK",     25,	"scaffold",
                                      "MELK",     38,	"scaffold",
                                      "MELK",     40,	"potency",
                                      "MELK",     57,	"potency",
                                      "MELK",     61,	"scaffold",
                                      "MELK",     70,	"selectivity",
                                      "MELK",     86,	"selectivity",  # Steffi: reverse selectivity
                                      "MELK",     87,	"scaffold",
                                      "MELK",     88,	"selectivity",
                                      "MELK",     89,	"selectivity",
                                      "MELK",     90,	"scaffold",
                                      "MELK",     92,	"scaffold",
                                      "MELK",     93,	"selectivity",
                                      "MELK",     136,	"scaffold",
                                      "MELK",     137,	"potency",
                                      "MELK",     139,	"scaffold",
                                      "MELK",     149,	"selectivity",  # Steffi: reverse selectivity
                                      "MELK",     150,	"potency",
                                      "MELK",     151,	"potency"
)


#############################################################################################################################

# Clustering with 7 aa groups (with the delta + factor not splitted)

set.seed(42)

base_subset = select(base,
                     gene_name,
                     residue_kinase,
                     kinome_wide_conservation,
                     bb,
                     sc, 
                     bbsc_total_observations,
                     targetable_inert,
                     conservation_overrepresentation_factor, 
                     increased_affinity,
                     functional_class,
                     kinase_position) %>%
  filter(gene_name %in% c("EPHA2", "ABL1", "MELK")) %>%
  sample_n(size = 300)



seven_aa_groups = base_subset %>% 
  mutate(seen_sc = if_else(condition = bbsc_total_observations != 0,
                           true = if_else(condition = sc >= 1,
                                          true = 1,
                                          false = 0),
                           false = 0)) %>% 
  mutate(aa_group = case_when(residue_kinase %in% c("A", "G", "I", "L", "P", "V", "M") ~ "aliphatic",
                              residue_kinase %in% c("F", "W", "Y")                     ~ "aromatic",
                              residue_kinase %in% c("D", "E")                          ~ "acidic",
                              residue_kinase %in% c("R", "H", "K")                     ~ "basic",
                              residue_kinase %in% c("S", "T")                          ~ "hydroxylic",
                              residue_kinase == "C"                                    ~ "cysteine",
                              residue_kinase %in% c("N", "Q")                          ~ "amidic"))

# add dummy encoding
aa_groups_factor = factor(seven_aa_groups$aa_group)
aa_groups_dummy_coding = model.matrix(~ -1 + aa_groups_factor) %>% 
  as.data.frame()
aa_groups_dummy_coding = aa_groups_dummy_coding/3

seven_aa_groups = bind_cols(seven_aa_groups, aa_groups_dummy_coding) %>% 
  select(-aa_group)

# add unique rownames for clustering
seven_aa_groups$unique_case_name = NA
for(i in 1:nrow(seven_aa_groups)){
  seven_aa_groups[i, "unique_case_name"] = paste0(i, "_", seven_aa_groups[i, "gene_name"])
}
rm(i)

seven_aa_groups = column_to_rownames(seven_aa_groups, var = "unique_case_name") %>% 
  dplyr::rename("acidic" = "aa_groups_factoracidic",
         "aliphatic" = "aa_groups_factoraliphatic",
         "amidic" = "aa_groups_factoramidic",
         "aromatic" = "aa_groups_factoraromatic",
         "basic" = "aa_groups_factorbasic",
         "cysteine" = "aa_groups_factorcysteine",
         "hydroxylic" = "aa_groups_factorhydroxylic")


# annotation ------------------------------------------------------------------

functional_class_myApp = select(seven_aa_groups, 
                                gene_name,
                                functional_class,
                                conservation_overrepresentation_factor,
                                increased_affinity,
                                bbsc_total_observations,
                                sc,
                                kinome_wide_conservation,
                                targetable_inert
                                ) %>% 
  
  # shiny app classification
  dplyr::rename("app_class" = "functional_class") %>%
  
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
  mutate(exposed_bb_OR_chem_inert = if_else(condition = (backbone_exposed == "backbone exposed" | targetable_inert == "aliphatic"),
                                            true = "yes", false = "no")) %>% 
  select(-backbone_exposed, -targetable_inert) %>% 
  
  # kw conservation level >= 50%?
  mutate(kw_cons_over_50 = if_else(condition = kinome_wide_conservation >= 50,
                                   true = "yes", false = "no")) %>% 
  select(-kinome_wide_conservation) %>% 
  
  rownames_to_column()

functional_class_myApp$app_class = factor(functional_class_myApp$app_class)


# Steffi manual classification
functional_class_Steffi_manual = select(seven_aa_groups, gene_name, kinase_position) %>% 
  rownames_to_column() %>% 
  inner_join(y = Steffi_supplementary_tables, by = c("gene_name" = "gene_name",
                                                     "kinase_position" = "kinase_position")) %>% 
  select(rowname, Steffi_manual_class)
functional_class_Steffi_manual$Steffi_manual_class = factor(functional_class_Steffi_manual$Steffi_manual_class)


row_annotation = full_join(functional_class_myApp, functional_class_Steffi_manual,
                           by = c("rowname" = "rowname")) %>% 
  column_to_rownames(var = "rowname") %>% 
  select(Steffi_manual_class, app_class, kw_cons_over_50, exposed_bb_OR_chem_inert, conv_increased_affinity, overrep_in_ts, gene_name)

classification_annotation_colors = list(
  app_class = c(key = "#009900",   # green
                potency = "#0066ff",  # blue
                scaffold = "#8c8c8c",    # grey
                selectivity = "#ff6600"),  # orange
  Steffi_manual_class =  c(key = "#009900",   # green
                           potency = "#0066ff",  # blue
                           scaffold = "#8c8c8c",    # grey
                           selectivity = "#ff6600"),  # orange
  gene_name = c(ABL1 = "#FF4500",
                EPHA2 = "#FFA500",
                MELK = "#6B8E23"),
  overrep_in_ts = c(yes = "green", no = "red"),
  conv_increased_affinity = c(yes = "green", no = "red"),
  exposed_bb_OR_chem_inert = c(yes = "green", no = "red"),
  kw_cons_over_50 = c(yes = "green", no = "red"),
  delta_factor_clusters = c("cluster_1" = "#FF0000FF",
                            "cluster_2" = "#FFFF00FF",
                            "cluster_3" = "#00FF00FF",
                            "cluster_4" = "#00FFFFFF",
                            "cluster_5" = "#0000FFFF",
                            "cluster_6" = "#FF00FFFF")
)


# scaling ---------------------------------------------------------------------

scored = select(seven_aa_groups, -bbsc_total_observations, -bb, -sc, -targetable_inert, -gene_name, -residue_kinase, -functional_class, -kinase_position)

scored$kinome_wide_conservation = scale_this_column_range_01(scored$kinome_wide_conservation)
scored$conservation_overrepresentation_factor = scale_this_column_range_01(scored$conservation_overrepresentation_factor)
scored$increased_affinity = scale_this_column_range_01(scored$increased_affinity)

scored$seen_sc = scored$seen_sc/3

scored = as.matrix(scored)

# plot ------------------------------------------------------------------------

my_heatmap = pheatmap(mat = scored,
                      cluster_rows = T,
                      cluster_cols = F,
                      clustering_distance_rows = "euclidean",    # distance measure
                      clustering_method = "ward.D2",             # agglomeration method of the hierarchical clustering (as in hclust)
                      
                      annotation_row = row_annotation_with_clusters,
                      cutree_rows = 6,
                      annotation_colors = classification_annotation_colors,
                      cellwidth = 10,
                      treeheight_row = 100,
                      show_rownames = F,
                      
                      color = viridis(200),
                      angle_col = 45,
                      fontsize = 7
)


# save_pheatmap_png(x = my_heatmap,
#                   filepath = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/plots/K.png",
#                   width = 600,
#                   height = 1200,
#                   res = 150)


# save the clusters for annotation in later plot ------------------------------

my_hclust_annotation = hclust(dist(scored), method = "ward.D2")

# as.dendrogram(my_hclust_annotation) %>% 
#   plot(horiz = T)

my_delta_cluster_annotation = cutree(tree = as.dendrogram(my_hclust_annotation), k = 6)




########### PCA ########################

my_pca = prcomp(scored,
                scale. = F,
                center = T,
                retx = T)

# str(my_pca)
# plot(my_pca)

# pca_labels = as.matrix(df_Steffi_classes)
# my_pca$x = cbind(my_pca$x, pca_labels)
# 
# str(my_pca)

group_colors = row_annotation$app_class

a_gg_plot = ggbiplot(my_pca,
                     ellipse = T,
                     choices = c(1,2),
                     groups = group_colors,
                     var.axes = T
)

a_gg_plot + scale_color_manual(name = "groups", values = c("key"="green", "potency"="blue", "scaffold"="grey", "selectivity"="orange"))
a_gg_plot + scale_color_manual(name = "groups", values = c("key"="green", "potency"="blue", "scaffold"="grey", "selectivity"="orange"))
a_gg_plot + scale_color_discrete(values = c("green", "blue", "grey", "orange"))
# nrow(my_pca$x)


############### UMAP ###########################################################################

# umap.defaults

my_umap = umap(d = scored)

head(my_umap)

source("lib_umap_plot.R")

plot_umap(x = my_umap,
          labels = row_annotation$app_class,
          main = "7 aa groups, dummy encoding, z-normalized umap",
          colors = c("#009900", "#0066ff", "#8c8c8c", "#ff6600")
)

# c(key = "#009900",   # green
#   potency = "#0066ff",  # blue
#   scaffold = "#8c8c8c",    # grey
#   selectivity = "#ff6600")  # orange

################ t-SNE ###################################################
library(Rtsne)

nrow(scored)

tsne_input = scored %>% 
  as.data.frame() %>% 
  distinct() %>% 
  as.matrix() %>% 
  normalize_input()


set.seed(90)
my_tsne = Rtsne(X = tsne_input,
                check_duplicates = F,
                pca_center = F,
                pca_scale = F,
                normalize = F,
                eta = 100,   # default: 200
                perplexity = 80) # default: 30  # MAX: 99

# Maria adapted:
# perpelxity 2, 5, 30, 50, 100
# eta (learning) 10, 20, 50, 100, 1000

plot(my_tsne$Y, col = row_annotation$app_class, asp = 1)




############ two values instead of delta/ factor ##########################

# Clustering with 7 aa groups and with two values (instead of delta/ factor)

set.seed(42)

seven_aa_groups_separated_values = select(base,
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
                                          kinase_position) %>%
  filter(gene_name %in% c("EPHA2", "ABL1", "MELK")) %>%
  sample_n(size = 300) %>% 
  mutate(seen_sc = if_else(condition = bbsc_total_observations != 0,
                           true = if_else(condition = sc >= 1,
                                          true = 1,
                                          false = 0),
                           false = 0)) %>% 
  mutate(aa_group = case_when(residue_kinase %in% c("A", "G", "I", "L", "P", "V", "M") ~ "aliphatic",
                              residue_kinase %in% c("F", "W", "Y")                     ~ "aromatic",
                              residue_kinase %in% c("D", "E")                          ~ "acidic",
                              residue_kinase %in% c("R", "H", "K")                     ~ "basic",
                              residue_kinase %in% c("S", "T")                          ~ "hydroxylic",
                              residue_kinase == "C"                                    ~ "cysteine",
                              residue_kinase %in% c("N", "Q")                          ~ "amidic"))

# add dummy encoding
aa_groups_factor = factor(seven_aa_groups_separated_values$aa_group)
aa_groups_dummy_coding = model.matrix(~ -1 + aa_groups_factor) %>% 
  as.data.frame()
aa_groups_dummy_coding = aa_groups_dummy_coding/3

seven_aa_groups_separated_values = bind_cols(seven_aa_groups_separated_values, aa_groups_dummy_coding) %>% 
  select(-aa_group)

# add unique rownames for clustering
seven_aa_groups_separated_values$unique_case_name = NA
for(i in 1:nrow(seven_aa_groups_separated_values)){
  seven_aa_groups_separated_values[i, "unique_case_name"] = paste0(i, "_", seven_aa_groups_separated_values[i, "gene_name"])
}
rm(i)

seven_aa_groups_separated_values = column_to_rownames(seven_aa_groups_separated_values, var = "unique_case_name") %>% 
  dplyr::rename("acidic" = "aa_groups_factoracidic",
                "aliphatic" = "aa_groups_factoraliphatic",
                "amidic" = "aa_groups_factoramidic",
                "aromatic" = "aa_groups_factoraromatic",
                "basic" = "aa_groups_factorbasic",
                "cysteine" = "aa_groups_factorcysteine",
                "hydroxylic" = "aa_groups_factorhydroxylic")


# scaling ---------------------------------------------------------------------

scored = select(seven_aa_groups_separated_values,
                -gene_name, -residue_kinase,
                -bbsc_total_observations, -bb, -sc,
                -targetable_inert, -functional_class, -kinase_position)

scored$kinome_wide_conservation = scale_this_column_range_01(scored$kinome_wide_conservation)
scored$target_wide_conservation = scale_this_column_range_01(scored$target_wide_conservation)

scored$median_pkDapp_M = scale_this_column_range_01(scored$median_pkDapp_M)
scored$weighted_median_pkDapp_M_different_aa = scale_this_column_range_01(scored$weighted_median_pkDapp_M_different_aa)

scored$seen_sc = scored$seen_sc/3

scored = as.matrix(scored)

# plot ------------------------------------------------------------------------

row_annotation = read_rds("result_2019_08_16_clustering_row_annotation")

my_delta_cluster_annotation = as.data.frame(my_delta_cluster_annotation)
my_delta_cluster_annotation$helper = "cluster_"
my_delta_cluster_annotation = unite(my_delta_cluster_annotation, helper, my_delta_cluster_annotation, sep = "", col = "delta_factor_clusters")

row_annotation = rownames_to_column(row_annotation)
my_delta_cluster_annotation = rownames_to_column(my_delta_cluster_annotation)
row_annotation_with_clusters = dplyr::inner_join(x = row_annotation, y = my_delta_cluster_annotation,
                                                 by = c("rowname" = "rowname")) %>% 
  column_to_rownames()

my_heatmap = pheatmap(mat = scored,
                      cluster_rows = T,
                      cluster_cols = F,
                      clustering_distance_rows = "euclidean",    # distance measure
                      clustering_method = "ward.D2",             # agglomeration method of the hierarchical clustering (as in hclust)
                      
                      annotation_row = row_annotation_with_clusters,
                      cutree_rows = 6,
                      annotation_colors = classification_annotation_colors,
                      cellwidth = 10,
                      treeheight_row = 100,
                      show_rownames = F,
                      
                      color = viridis(200),
                      angle_col = 45,
                      fontsize = 7
)

save_pheatmap_png(x = my_heatmap,
                  filepath = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/plots/M.png",
                  width = 700,
                  height = 1200,
                  res = 150)


