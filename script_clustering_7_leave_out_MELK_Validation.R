
##########################################################################
#                                                                        #
#  Clustering.7: Clustering with 2 Kinases, then add MELK as Validation  #
#                                                                        #
##########################################################################

# Dieses Skript hangelt sich entlang an: http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/#supplementary-individuals

library(tidyverse)


# Prapare dataset -------------------------------------------------------------

base = read_csv("Z:/users_files/Verena Burger/11_PCA_and_Clustering/Kinase_Residue_Classification_2019-09-08.csv")

df_classification = select(base,
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
  filter(gene_name %in% c("JAK1", "MELK", "EGFR", "MAPK1", "EPHA2")) %>%
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
aa_groups_factor = factor(df_classification$aa_group)
aa_groups_dummy_coding = model.matrix(~ -1 + aa_groups_factor) %>% 
  as.data.frame()

aa_groups_dummy_coding = aa_groups_dummy_coding/2

df_classification = bind_cols(df_classification, aa_groups_dummy_coding) %>% 
  select(-aa_group)

rm(aa_groups_factor, aa_groups_dummy_coding)

# add unique rownames for clustering
df_classification$unique_case_name = NA
for(i in 1:nrow(df_classification)){
  df_classification[i, "unique_case_name"] = paste0(i, "_", df_classification[i, "gene_name"])
}
rm(i)

df_classification = df_classification %>% 
  
  column_to_rownames(var = "unique_case_name") %>% 
  
  dplyr::rename("acidic" = "aa_groups_factoracidic",
                "aliphatic" = "aa_groups_factoraliphatic",
                "amidic" = "aa_groups_factoramidic",
                "aromatic" = "aa_groups_factoraromatic",
                "basic" = "aa_groups_factorbasic",
                "cysteine" = "aa_groups_factorcysteine",
                "hydroxylic" = "aa_groups_factorhydroxylic")

# Annotation ------------------------------------------------------------------

row_annotation = df_classification %>% 
  dplyr::select(functional_class)

row_annotation$functional_class = factor(row_annotation$functional_class)


# Finish the data preparation -------------------------------------------------

df_classification = df_classification %>% 
  
  dplyr::select(-gene_name, -residue_kinase,
                -bbsc_total_observations, -bb, -sc,
                -targetable_inert, -functional_class, -kinase_position,
                -conservation_overrepresentation_factor, -increased_affinity) %>% 
  
  dplyr::rename("affinity_same_aa" = "median_pkDapp_M",
                "affinity_diff_aa" = "weighted_median_pkDapp_M_different_aa",
                "exposed_backbone" = "seen_sc")


# Scaling ---------------------------------------------------------------------

scale_this_column_range_01 = function(col){
  my_range = range(col)  # lowest and highest value of the column
  col = col - my_range[1]    # so dass der niedrigste Wert bei 0 liegt
  col = col / (my_range[2] - my_range[1])  # jeden Wert durch die originale "Spannweite" teilen
  # man bekommt quasi einen "Prozentwert", wo der Wert liegt. => das muss
  # also immer zwischen 0 und 1 liegen.
  return(col)
}

scored = df_classification

scored$kinome_wide_conservation = scale_this_column_range_01(scored$kinome_wide_conservation)
scored$target_wide_conservation = scale_this_column_range_01(scored$target_wide_conservation)

scored$affinity_same_aa = scale_this_column_range_01(scored$affinity_same_aa)
scored$affinity_diff_aa = scale_this_column_range_01(scored$affinity_diff_aa)

scored$exposed_backbone = scored$exposed_backbone/2

scored = as.matrix(scored)


# Splitting up into active and supplementary ----------------------------------

# "active individuals": JAK1, MELK, EGFR, MAPK1
train = str_detect(string = rownames(scored),
                   pattern = "EPHA2",              
                   negate = T)                            

# supplementary individuals: EPHA2
validate = str_detect(string = rownames(scored),
                      pattern = "EPHA2",
                      negate = F)


# perform the PCA -------------------------------------------------------------
my_pca_object = prcomp(scored[train,],
                       scale. = F,
                       center = T,
                       retx = T)

library(factoextra)

# Visualize eigenvalues (scree plot)
fviz_eig(my_pca_object)

# variables
fviz_pca_var(my_pca_object,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


# predict supplementary individuals: EPHA2 ------------------------------------

suppl_predicted_coordinates = stats::predict(my_pca_object, newdata = scored[validate,])


# plot it all -----------------------------------------------------------------
# http://huboqiang.cn/2016/03/03/RscatterPlotPCA

# library(grid)
# library(gridExtra)

# add variables for grouping

pca_train = my_pca_object$x[, c("PC1", "PC2")] %>% 
  as.data.frame() %>% 
  add_column(origin = "train") %>% 
  rownames_to_column()


pca_validate = suppl_predicted_coordinates[, c("PC1", "PC2")] %>% 
  as.data.frame() %>% 
  add_column(origin = "validate") %>% 
  rownames_to_column()

row_annotation = rownames_to_column(row_annotation)

pca_complete = bind_rows(pca_train, pca_validate) %>% 
  left_join(y = row_annotation, by = c("rowname" = "rowname")) %>% 
  column_to_rownames() %>% 
  unite(col = color_separator, origin, functional_class, sep = "_", remove = F)


explained_variance = round(my_pca_object$sdev / sum(my_pca_object$sdev) * 100, 1)
explained_variance = paste0(colnames(as.data.frame(my_pca_object$x)),
                           " (",
                           paste0(as.character(explained_variance),
                                 "% explained variance)"
                                 )
                           )



knorke_plot = ggplot(data = pca_complete, mapping = aes(x = PC1, y = PC2)) +
  
  geom_point(mapping = aes(alpha = origin,
                           color = color_separator,
                           size = origin,
                           shape = origin)) +
  
  scale_alpha_manual(values = c("train" = 0.15, "validate" = 0.7)) +
  
  scale_color_manual(values = c("train_key" = "#00e600",
                                "train_potency" = "#80b3ff",
                                "train_scaffold" = "#cccccc",
                                "train_selectivity" = "#ffa366",
                                "validate_key"="#009900", # green
                                "validate_potency"="#0066ff", # blue
                                "validate_scaffold"="#8c8c8c", # grey
                                "validate_selectivity"="#ff6600", # orange
                                "key"="#009900", # green
                                "potency"="#0066ff", # blue
                                "scaffold"="#8c8c8c", # grey
                                "selectivity"="#ff6600" # orange
                                )) +  
  
  scale_size_manual(values = c("train" = 7, "validate" = 2)) +

  theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"),
        legend.key = element_blank()) +
  
  xlab(explained_variance[1]) +  ylab(explained_variance[2]) +
  
  stat_ellipse(mapping = aes(color = functional_class))

knorke_plot










  
















