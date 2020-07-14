
######################################
#                                    #
#  Generating the Kinobead Call Set  #
#                                    #
######################################

# pKD input matrix/ data frame
require(tidyverse)
require(readxl)
paper_path = "Z:/users_files/Verena Burger/4_datasets/current/Kinobeads/The_target_landscape_of_clinical_kinase_drugs/aan4368_Tables1-11/aan4368_Table_S2.xlsx"

paper = read_excel(path = paper_path,
                   sheet = "Kinobeads") %>% 
  select("Drug", "Gene Name", "Apparent Kd") %>%          # long format, Mathias had it in wide format
  drop_na("Apparent Kd")

paper = rename(paper, "apparent_Kd_nM" = "Apparent Kd",
               "compound" = "Drug",
               "gene_name" = "Gene Name")

# write_csv(paper, "result_2019_07_31_affinity_upload_KB.csv")
write_csv(paper, "C:/Users/vburger/code/partitioned_residue_classification/Residue_classification_Shiny_App/www/AFFINITY_UPLOAD.csv")

