library(tidyverse)
source("lib_validation_formattable_color_Steffis_classification.R")




#####################

classification_steffi_EPHA2_Dasatinib = data.frame(
  
  key = c("M667", "I690", "T692",
          "M695", "N697",
          NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  scaffold = c("I619", "G620", "A621", "G622",
               "V627", "A644", "I666", "I675",
               "I676", "E693", "E696", "G698",
               "L730", "R743", "L746", "V755"),
  
  potency = c("K646", "E663", "H737",
              "N744", "D757", "F758",
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  
  selectivity = c("K617", "E623", "F670",
                  "Y694", "A699", "K702",
                  "E706", "Y735", "S756",
                  NA, NA, NA, NA, NA, NA, NA)
  
)


###############################################################################

# Filters" EPHA2, Dasatinib!
colored_table_EPHA2_Dasatinib = color_Steffis_classification_according_to_Veris(df_Steffis_classification = classification_steffi_EPHA2_Dasatinib,
                                                                                filepath_download_from_App_filtered = "Z:/users_files/Verena Burger/11_PCA_and_Clustering/Steffi_Classifications/EPHA2/my_app_classification_EPHA2_Dasatinib_umgestelltes_class.csv")
colored_table


