
#############################################
#                                           #
#  Clustering.5: Punkte in PCA verschieben  #
#                                           #
#############################################

library(tidyverse)
library(ggbiplot)
scored = read_rds("result_2019_09_02_scored")
annotation = read_rds("result_2019_09_02_row_annotation")


my_pca = prcomp(scored,
                scale. = F,
                center = T,
                retx = T)

group_colors = annotation$app_class

my_pca_plot = ggbiplot(my_pca,
                     ellipse = T,
                     choices = c(1,2),
                     groups = group_colors,
                     var.axes = T
) + scale_color_manual(name = "groups", values = c("key"="#009900", # green
                                                   "potency"="#0066ff", # blue
                                                   "scaffold"="#8c8c8c", # grey
                                                   "selectivity"="#ff6600", # orange
                                                   "focus" = "red")) + 
  theme(legend.title = element_blank())

my_pca_plot



# pick a point and change its functional class to "red"
# and subsequently expand the annotation colors to color the functional class "red" red.

# which(rownames(annotation) == "28_ABL1")   # 28
levels(annotation[, "app_class"]) = c(levels(annotation[, "app_class"]), "focus")
annotation[28, "app_class"] = "focus"


# scored[28, "aliphatic"] = 0.5
# scored[28, "aromatic"] = 0.5
scored[28, "target_wide_conservation"] = 0.2

View(scored[28,])


