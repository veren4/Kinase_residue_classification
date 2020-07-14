require(formattable)





color_Steffis_classification_according_to_Veris = function(df_Steffis_classification,
                                                           filepath_download_from_App_filtered){
  
  
  app_download = read_csv(filepath_download_from_App_filtered) %>% 
    select(residue_kinase, kinase_position, functional_class) %>% 
    unite(residue_kinase, kinase_position,
          col = "residue",
          sep = "",
          remove = T)
    
  my_classification_key = filter(app_download, functional_class == "key") %>% 
    pull(residue)
  
  my_classification_scaffold = filter(app_download, functional_class == "scaffold") %>% 
    pull(residue)
  
  my_classification_potency = filter(app_download, functional_class == "potency") %>% 
    pull(residue)
  
  my_classification_selectivity = filter(app_download, functional_class == "selectivity") %>% 
    pull(residue)
  
  all_my_cases = c(my_classification_key, my_classification_scaffold, my_classification_selectivity, my_classification_potency)
  
  color_picker_na = function(z){
    
    if(is.na(z)){return("white")}
    
  }
  
  
  color_picker_key = function(z){
    
    if(is.na(z)){return("white")}                                # NA
    else if(z %in% my_classification_key){return("#CCFF99")}     # identically classified
    else if(!(z %in% all_my_cases)){return("#FFFF99")}              # missing
    else{return("#FF9999")}  # differently classified
    
  }
  
  color_picker_scaffold = function(z){
    
    if(is.na(z)){return("white")}                                # NA
    else if(z %in% my_classification_scaffold){return("#CCFF99")}     # identically classified
    else if(!(z %in% all_my_cases)){return("#FFFF99")}              # missing
    else{return("#FF9999")}  # differently classified
    
  }
  
  color_picker_potency = function(z){
    
    if(is.na(z)){return("white")}                                # NA
    else if(z %in% my_classification_potency){return("#CCFF99")}     # identically classified
    else if(!(z %in% all_my_cases)){return("#FFFF99")}              # missing
    else{return("#FF9999")}  # differently classified
    
  }
  
  color_picker_selectivity = function(z){
    
    if(is.na(z)){return("white")}                                # NA
    else if(z %in% my_classification_selectivity){return("#CCFF99")}     # identically classified
    else if(!(z %in% all_my_cases)){return("#FFFF99")}              # missing
    else{return("#FF9999")}  # differently classified
    
  }
  
  result = formattable(x = df_Steffis_classification,
              list(
                
                key = formatter("span",
                                style = x ~ style(display = "block",
                                                  "border-radius" = "4px",
                                                  "padding-right" = "4px",
                                                  color = sapply(x, color_picker_na),
                                                  "background-color" = sapply(x, color_picker_key))
                ),
                
                scaffold = formatter("span",
                                     style = x ~ style(display = "block",
                                                       "border-radius" = "4px",
                                                       "padding-right" = "4px",
                                                       color = sapply(x, color_picker_na),
                                                       "background-color" = sapply(x, color_picker_scaffold))
                ),
                
                potency = formatter("span",
                                    style = x ~ style(display = "block",
                                                      "border-radius" = "4px",
                                                      "padding-right" = "4px",
                                                      color = sapply(x, color_picker_na),
                                                      "background-color" = sapply(x, color_picker_potency))
                ),
                
                selectivity = formatter("span",
                                        style = x ~ style(display = "block",
                                                          "border-radius" = "4px",
                                                          "padding-right" = "4px",
                                                          color = sapply(x, color_picker_na),
                                                          "background-color" = sapply(x, color_picker_selectivity))
                )
                
                
                
              ) # list
  ) # formattable
  
  return(result)
  
}