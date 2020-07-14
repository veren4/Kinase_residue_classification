
#################
#  Validierung  #
#################

require(tidyverse)
require(formattable)

App_data_Dasatinib_EPHA2 = read_csv("C:/Users/vburger/Desktop/Kinase_Residue_Classification_2019-07-31.csv")

my_download = read_csv("C:/Users/vburger/Desktop/Kinase_Residue_Classification_2019-07-31.csv") %>% 
  select(residue_kinase, kinase_position, functional_class) %>% 
  unite(residue_kinase, kinase_position,
        col = "residue",
        sep = "",
        remove = T) %>% 
  
  group_by(functional_class) %>%    # workaround for spread
  mutate(grouped_id = row_number()) %>%   # also workaround
  
  spread(key = functional_class,
         value = residue,
         drop = F) %>% 
  select(-grouped_id)


Steffi_EPHA2_key = c("M667", "I690", "T692", "M695", "N697")
Steffi_EPHA2_scaffold = c("I619", "G620", "A621", "G622", "V627", "A644", "I666", "I675", "I676", "E693", "E696", "G698", "L730", "R743", "L746", "V755")
Steffi_EPHA2_potency = c("K646", "E663", "H737", "N744", "D757", "F758")
Steffi_EPHA2_selectivity = c("K617", "E623", "F670", "Y694", "A699", "K702", "E706", "Y735", "S756")

all_of_Steffis_EPHA2_cases = c(Steffi_EPHA2_key, Steffi_EPHA2_scaffold, Steffi_EPHA2_potency, Steffi_EPHA2_selectivity)

my_additional_key = drop_na(my_download, key) %>% pull(key)
my_additional_key = my_additional_key[!(my_additional_key %in% Steffi_EPHA2_key)]
my_additional_key = my_additional_key[!(my_additional_key %in% all_of_Steffis_EPHA2_cases)]

my_additional_scaffold = drop_na(my_download, scaffold) %>% pull(scaffold)
my_additional_scaffold = my_additional_scaffold[!(my_additional_scaffold %in% Steffi_EPHA2_scaffold)]
my_additional_scaffold = my_additional_scaffold[!(my_additional_scaffold %in% all_of_Steffis_EPHA2_cases)]

my_additional_potency = drop_na(my_download, potency) %>% pull(potency)
my_additional_potency = my_additional_potency[!(my_additional_potency %in% Steffi_EPHA2_potency)]
my_additional_potency = my_additional_potency[!(my_additional_potency %in% all_of_Steffis_EPHA2_cases)]

my_additional_selectivity = drop_na(my_download, selectivity) %>% pull(selectivity)
my_additional_selectivity = my_additional_selectivity[!(my_additional_selectivity %in% Steffi_EPHA2_selectivity)]
my_additional_selectivity = my_additional_selectivity[!(my_additional_selectivity %in% all_of_Steffis_EPHA2_cases)]








