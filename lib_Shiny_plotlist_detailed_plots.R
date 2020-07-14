require(tidyverse)
require(ggrepel)


make_plotlist_detailed_plots = function(information_about_all_kinases,
                                        current_filtered_information,
                                        kinome_wide_conservation_helper_df,
                                        target_space_conservation_helper_df){
  
  
  all_kinases_information = information_about_all_kinases

  my_subset = current_filtered_information

  kinome_wide_conservation = kinome_wide_conservation_helper_df

  target_space_conservation = target_space_conservation_helper_df
  

  aas <- c("G", "A", "V", "L", "I", "M", "F", "W", "P", "S", "T", "C", "Y", "N", "Q", "D", "E", "K", "R", "H", "-")
  
  plotlist_detailed = list()
  
  
  withProgress(message = "Plotting...",
               expr = {
                 
                 for(ff in 1:nrow(my_subset)){
                   incProgress()
                   
                   # make shared y label for this row -------------------------
                   
                   current_row_functional_class = as.character(my_subset[ff, "functional_class"])
                  
                   class_color = case_when(current_row_functional_class == "key"         ~ "#009900", # green
                                           current_row_functional_class == "potency"     ~ "#0066ff", # blue
                                           current_row_functional_class == "scaffold"    ~ "#8c8c8c", # grey
                                           current_row_functional_class == "selectivity" ~ "#ff6600") # orange
                   
                   y_label_plot = ggplot() + 
                     annotate("text", x = 4, y = 25, size = 6, label = paste(current_row_functional_class, "residue"),
                              colour = class_color, angle = 90) + 
                     theme_void() +
                     theme(panel.grid.major=element_blank(),
                           panel.grid.minor=element_blank())
                   
                   
                   # make df for conservation plot ---------------------------------------------
                   
                   kw = dplyr::filter(kinome_wide_conservation, alignment_position == my_subset[ff, "alignment_position"])
                   
                   df_conservation_plot_helper = data.frame(aa = aas, stringsAsFactors = F) %>% 
                     left_join(kw, by = c("aa" = "residue_kinase"))
                   
                   current_row_drug = as.character(my_subset[ff, "compound"])
                   
                   row_kinase_alignment_position = as.numeric(my_subset[ff, "alignment_position"])
                   
                   current_row_kinase = as.character(my_subset[ff, "gene_name"])
                   
                   
                   if(round(sum(df_conservation_plot_helper$kinome_wide_conservation, na.rm = T), digits = 3) == 100){   ### ist die kw conservation summe 100?
                   
                     if(sum(is.na(df_conservation_plot_helper$kinome_wide_conservation)) > 0){ ### falls ja: gibt es NAs?
                       
                       df_conservation_plot_helper[which(is.na(df_conservation_plot_helper$kinome_wide_conservation)), "kinome_wide_conservation"] = 0  ### falls ja, auffuellen
                     
                     } # ggf. NAs in kw cons. aufgefuellt
                     
                     # ggf. alignment positon auffuellen
                     df_conservation_plot_helper[which(is.na(df_conservation_plot_helper$alignment_position)), "alignment_position"] = my_subset[ff, "alignment_position"]
                       
                     ts = dplyr::filter(target_space_conservation,
                                        compound == current_row_drug & alignment_position == row_kinase_alignment_position)
                     
                     df_conservation_plot_helper = left_join(df_conservation_plot_helper, ts, by = c('alignment_position' = 'alignment_position',
                                                                                                     'aa' = 'residue_kinase'))
                     
                     if(round(sum(df_conservation_plot_helper$target_wide_conservation, na.rm = T), digits = 3) == 100){  ### ist die ts cons. summe 100?
                     
                       if(sum(is.na(df_conservation_plot_helper$target_wide_conservation)) > 0){  #### is iwo in der ts conservation ein NA?
                     
                         df_conservation_plot_helper[which(is.na(df_conservation_plot_helper$target_wide_conservation)), "target_wide_conservation"] = 0 ### ggf. NAs auffuellen
                           
                       }
                       
                       # hier normal weiter
                       
                       my_aa = which(df_conservation_plot_helper$aa == as.character(my_subset[ff, "residue_kinase"]))
                       
                       df_conservation_plot_helper$color = "a"
                       df_conservation_plot_helper[my_aa, "color"] = "b"
                       
                       df_conservation_plot_helper[my_aa, "my_aa_kw_label"] = paste0(" ", round(df_conservation_plot_helper[my_aa, "kinome_wide_conservation"]), "%")
                       df_conservation_plot_helper[my_aa, "my_aa_ts_label"] = paste0(" ", round(df_conservation_plot_helper[my_aa, "target_wide_conservation"]), "%")
                       
                       df_conservation_plot_helper$aa = factor(df_conservation_plot_helper$aa, levels = aas)
                       
                       # make the conservation plot ------------------------------------------------
                       
                       conservation_plot = ggplot(data = df_conservation_plot_helper) +                # data = conservation_per_position[[wo]]
                         geom_col(mapping = aes(x = aa, y = target_wide_conservation, fill = color)) +
                         geom_col(mapping = aes(x = aa, y = -kinome_wide_conservation, fill = color)) +
                         xlab("aa or gap") +
                         ylab("conservation level [%]\n\n        kinome-wide                   target space-wide") +
                         ggtitle(label = paste("\nObserved aas in kinase position", as.numeric(my_subset[ff, "kinase_position"])),
                                 subtitle = paste("Coloring based on amino acid in kinase", current_row_kinase)) +
                         coord_cartesian(ylim = c(-100,100)) +
                         scale_y_continuous(labels = function(x)abs(x)) +
                         geom_hline(yintercept = 0) +
                         scale_fill_manual(values = c("a" = "grey",
                                                      "b" = "royalblue"),
                                           guide = F) +
                         geom_text(mapping = aes(x = aa, y = target_wide_conservation, label = my_aa_ts_label),
                                   vjust = -0.5) +
                         geom_text(mapping = aes(x = aa, y = -kinome_wide_conservation, label = my_aa_kw_label),
                                   vjust = 1.4)
                       
                       
                     }else{ # wenn die ts summe nicht 100 ist:
                       
                       my_text = paste("\n   The conservation plot cannot be calculated\n",
                                       "       for this position because\n",
                                       "       the sum of the target space-wide conservations of\nthe aas in this position is not 100 %.")
                       
                       conservation_plot = ggplot() + 
                         annotate("text", x = 4, y = 25, size = 5, label = my_text) + 
                         theme_void() +
                         theme(panel.grid.major=element_blank(),
                               panel.grid.minor=element_blank())
                       
                       }
                       
                       
                     }else{ # wenn die kw summe nicht 100 ist:
                       
                       my_text = paste("\n   The conservation plot cannot be calculated\n",
                                       "       for this position because\n",
                                       "       the sum of the kinome-wide conservations of\nthe aas in this position is not 100 %.")
                       
                       conservation_plot = ggplot() + 
                         annotate("text", x = 4, y = 25, size = 5, label = my_text) + 
                         theme_void() +
                         theme(panel.grid.major=element_blank(),
                               panel.grid.minor=element_blank())
                       
                     }
                       
                
                   # make affinity plot helper -------------------------------------------------
                   
                   current_row_residue = as.character(my_subset[ff, "residue_kinase"])
                   
                   df_affinity_plot_helper = select(all_kinases_information, compound, gene_name, pKDapp_M, residue_kinase,
                                                    alignment_position, kinase_position) %>% 
                     distinct() %>% 
                     dplyr::filter(compound == current_row_drug & alignment_position == row_kinase_alignment_position) %>% 
                     # von der Kinase der row die kinase position, von allen anderen kinasen
                     # die alignment_position
                     
                     dplyr::mutate(color = case_when(
                       
                       gene_name == current_row_kinase       ~ "c",  # the current kinase
                       
                       residue_kinase == current_row_residue ~ "b",  # kinases with the same aa
                       
                       residue_kinase != current_row_residue ~ "a"   # kinases with a different aa
                       
                     )) %>% 
                     
                     mutate(repel = if_else(condition = gene_name == current_row_kinase,  # labels for ggrepel
                                            true = current_row_kinase,
                                            false = ""
                     )
                     )
                   
                   # to have all 21 cases of aa:
                   df_affinity_plot_helper$residue_kinase = factor(df_affinity_plot_helper$residue_kinase, levels = aas)
                   
                   
                   # make the affinity plot ----------------------------------------------------

                   affinity_plot = ggplot(data = df_affinity_plot_helper, mapping = aes(x = residue_kinase, y = pKDapp_M)) +
                     scale_x_discrete(drop = F) +
                     geom_dotplot(mapping = aes(fill = color, width = 1),
                                  binaxis = "y",
                                  stackdir = "center",
                                  binwidth = 0.01,
                                  stackratio = 1,
                                  position = position_jitter(width = 0.03, height = 0.04, seed = 12),
                                  dotsize = 10, # previously 5.5
                                  color = NA) +
                     xlab("aa or gap") +
                     ylab(bquote(pK[D]^app~.(current_row_drug)~"[M]")) +    # can't add a \n for distance to plot
                     geom_hline(yintercept = 6,
                                linetype = "dashed") +
                     scale_fill_manual(values = c("a" = "#999999",
                                                  "b" = "royalblue",
                                                  "c" = "#CC0000"),
                                       guide = F) +
                     geom_text_repel(mapping = aes(x = residue_kinase, y = pKDapp_M, label = repel),
                                     box.padding = 0.2,
                                     segment.alpha = 0,
                                     color = "#CC0000") +
                     ggtitle(label = paste("\nAffinity of", current_row_drug, "to diff. kinases, separated by aa in position", my_subset[ff, "kinase_position"]),
                             subtitle = paste("Coloring based on amino acid in kinase", current_row_kinase)) +
                     stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.55) +
                     coord_cartesian(ylim = c(4, 10))
                   
                   
                   # add the 3 plots to a plotlist ----------------------------------------------
                   
                   plotlist_detailed[[length(plotlist_detailed)+1]] = y_label_plot
                   
                   plotlist_detailed[[length(plotlist_detailed)+1]] = conservation_plot
                   
                   plotlist_detailed[[length(plotlist_detailed)+1]] = affinity_plot
                   
                   
                   df_conservation_plot_helper = NULL
                   df_affinity_plot_helper = NULL
                   
                   
                 } # for(ff in 1:nrow(my_subset))
                 
               } # withProgress expr
    
  ) # withProress
  
  
  return(plotlist_detailed)

}




