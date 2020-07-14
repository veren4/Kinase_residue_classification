
####################################################
#                                                  #
#  Clustering.6: comparing the clustering outcome  #
#                                                  #
####################################################


library(tidyverse)
# annotation = read_rds("result_2019_09_02_row_annotation")
annotation = row_annotation
annotation = dplyr::select(annotation, `manual classification`, `automatic classification`)

# colnames(my_table) = c("manual", "automatic")
my_table = annotation %>% 
  tidyr::drop_na() %>% 
  dplyr::rename(manual = `manual classification`, automatic = `automatic classification`) %>%
  tidyr::gather(manual, automatic, key = "classification_method", value = "functional_class") %>% 
  dplyr::group_by(classification_method, functional_class) %>% 
  dplyr::count(name = "total_observations") %>% 
  dplyr::ungroup() %>% 
  
  # group_sums = my_table %>%
  #   group_by(classification_method) %>%
  #   summarise(sumsi = sum(total_observations)) # 229
  
  mutate(percentage = round((total_observations/3317)*100, digits = 0))


my_table$classification_method = factor(my_table$classification_method)
my_table$functional_class = factor(my_table$functional_class, levels = c("key", "scaffold", "potency", "selectivity"))

facet_labels = c("automatic" = "automatic classification", "manual" = "manual classification")

cont_plot = ggplot(data = my_table, mapping = aes(x = functional_class, y = total_observations, fill = functional_class)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(classification_method),
             labeller = labeller(classification_method = facet_labels)) +
  theme_light() +      # oder: theme_minimal()
  scale_fill_manual(values = c("key"="#009900", # green
                               "potency"="#0066ff", # blue
                               "scaffold"="#8c8c8c", # grey
                               "selectivity"="#ff6600")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  xlab(element_blank()) +
  ylab("Frequency") + 
  labs(title = "Classification comparison R/S") +
  geom_text(mapping = aes(label = paste0(percentage, "%")),
            position = position_dodge(width = 0.9),
            vjust = 1.5, size = 3.5, color = "white")

cont_plot




# write.csv(my_table, file = "todays_table")
# write_rds(my_table, "Z:/users_files/Verena Burger/2_thesis/figures/classification_comparison_data_Robject_2")



