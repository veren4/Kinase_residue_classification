
##########
#  MELK  #
##########

library(tidyverse)

supplementary_table = readxl::read_xlsx(path = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/Steffi_Classifications/MELK/aan4368_Tables1-11/aan4368_Table_S11.xlsx",
                                        sheet = "Residue classification")

# drugs = colnames(supplementary_table)[10:14]

Steffis_class_MELK = select(supplementary_table, Residue, Position, Annotation) %>% 
  unite(Residue, Position, col = "residue", sep = "") %>% 
  
  group_by(Annotation) %>%     # workaround for spread
  mutate(grouped_id = row_number()) %>% 
  
  spread(key = Annotation,
         value = residue,
         drop = F) %>% 
  select(-grouped_id)
  

my_app_classification = read_csv(file = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/Steffi_Classifications/MELK/my_app_classification_filteredfor_MELK.csv")


source("lib_validation_formattable_color_Steffis_classification.R")

colored_table_MELK = color_Steffis_classification_according_to_Veris(df_Steffis_classification = Steffis_class_MELK,
                                                                     filepath_download_from_App_filtered = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/Steffi_Classifications/MELK/my_app_classification_filteredfor_MELK.csv")
colored_table_MELK
