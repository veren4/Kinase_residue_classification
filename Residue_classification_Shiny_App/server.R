
############
#          #
#  server  #
#          #
############

server <- function(input, output, session){
  
  session$userData$all_kinases_information_interacting <- data.frame(gene_name = character(), compound= character(), kinase_position= integer(), functional_class = character(), kinase_family = character())
  
  output$detailed_plots = renderPlot(
    placeholder_detailed_plots,
    bg = "transparent",
    width = 1200,
    height = 300
  )
  
  output$heatmap_plot = renderPlot(
    placeholder_overview_plot,
    bg = "transparent",
    width = 1200,
    height = 300
  )
  
  r_upload <- reactive({
    
    if(is.null(input$affinity_upload)){
      
      return(data.frame(gene_name = character(), compound = character(), kinase_position = integer(), functional_class = character(), kinase_family = character()) )
      
    }else{
      
      # Check the input -------------------------------------------------------
      
      # Is it a csv file?
      if(file_ext(input$affinity_upload$datapath) != "csv"){
        
        showModal(modalDialog(title = "Input error",
                              "The selected file doesn't have a .csv file extension and is therefore not considered a csv file.\n
                              Please select a suitable input file. For input requirements, see Help panel.",
                              easyClose = T,
                              fade = T,
                              footer = NULL
                              
        ))
        return(data.frame(gene_name = character(), compound= character(), kinase_position= integer()) )
      }
      
      # are there column headers, and are they correct?
      if(readLines(input$affinity_upload$datapath, n = 1) != "compound,gene_name,apparent_Kd_nM"){
        
        showModal(modalDialog(title = "Input error",
                              "The uploaded file appears not to have the required headers.\n
                              These are: \"compound\", \"gene_name\" and \"apparent_Kd_nM\" (see Help panel).",
                              easyClose = T,
                              fade = T,
                              footer = NULL
                              
        ))
        return(data.frame(gene_name = character(), compound= character(), kinase_position= integer()) )
      }
      
      # Assignment
      df_affinity = read_csv(file = input$affinity_upload$datapath) %>% as_tibble()
      # df_affinity = read_csv("C:/Users/vburger/code/partitioned_residue_classification/Residue_classification_Shiny_App/www/AFFINITY_UPLOAD.csv") %>% as_tibble()
      
      # check: column number + names + data types
      if(ncol(df_affinity) != 3){
        
        showModal(modalDialog(title = "Input error",
                              "The uploaded file appears not to have 3 columns.\n
                              Please select a suitable input file. For input requirements, see Help panel.",
                              easyClose = T,
                              fade = T,
                              footer = NULL
                              
        ))
        return(data.frame(gene_name = character(), compound= character(), kinase_position= integer()) )
      }
      
      if(class(df_affinity$compound) != "character" | class(df_affinity$gene_name) != "character" | class(df_affinity$apparent_Kd_nM) != "numeric"){
        
        showModal(modalDialog(title = "Input error",
                              "The columns in the uploaded files have unexpected data types.\n
                              For input requirements, see Help panel",
                              easyClose = T,
                              fade = T,
                              footer = NULL
                              
        ))
        return(data.frame(gene_name = character(), compound= character(), kinase_position= integer()) )
      }
      
      
      # Calculate residue information -----------------------------------------
      
      df_affinity$pKDapp_M = -log10(df_affinity$apparent_Kd_nM*1E-9)
      # start = Sys.time()
      # all_kinases_information_list <- calculate_information_for_all_kinases(multiple_sequence_alignment = MSA,
      #                                                                       generic_structural_information = df_structure_information,
      #                                                                       targetable_inert_helper_df = df_targetable_inert_helper,
      #                                                                       uploaded_affinities = df_affinity)
      # end = Sys.time()
      # cat("Calculating the kinase information for all kinases took ", end-start, "\n")
      all_kinases_information_list = read_rds("Z:/users_files/Verena Burger/1_presentations/Abschluss_Praesentation/all_kinases_information_list")
      
      session$userData$all_kinases_information_unfiltered <- all_kinases_information_list$kinase_information_matrix
      session$userData$all_kinases_information_interacting <- all_kinases_information_list$kinase_information_matrix %>% 
        dplyr::filter(interacting_residue == 1)
      session$userData$kinome_wide_conservation <- all_kinases_information_list$kw_conservation_helper
      session$userData$target_space_conservation <- all_kinases_information_list$ts_conservation_helper
      
      updateCheckboxGroupInput(session = session,
                               inputId = "display_columns_selection",
                               selected = c("gene_name", "compound", "pKDapp_M", "kinase_position", "residue_kinase", "functional_class"))
      
      enable(id = "display_columns_selection")
      enable(id = "download_table")
      enable(id = "overview_plots")
      
      # Reset the plot spaces in case a new dataset gets uploaded:
      output$detailed_plots = renderPlot(
        placeholder_detailed_plots,
        bg = "transparent",
        width = 1200,
        height = 300
      )
      output$heatmap_plot = renderPlot(
        placeholder_overview_plot,
        bg = "transparent",
        width = 1200,
        height = 300
      )
      
      return(session$userData$all_kinases_information_interacting)
      
      } # else
    
  })# r_upload
  
  

  r_selection_kinase <- reactive({
    if (is.null(input$kinase_selection)) {
      return(!logical(length = nrow(session$userData$all_kinases_information_interacting)))
    } else{
      return(session$userData$all_kinases_information_interacting$gene_name %in% input$kinase_selection)
    }
  })
  
  
  r_selection_compound <- reactive({
    if (is.null(input$compound_selection)) {
      return(!logical(length = nrow(session$userData$all_kinases_information_interacting)))
    } else{
      return(session$userData$all_kinases_information_interacting$compound %in% input$compound_selection)
    }
  })
  
  r_selection_position <- reactive({
    if (is.null(input$position_selection)) {
      return(!logical(length = nrow(session$userData$all_kinases_information_interacting)))
    } else{
      return(session$userData$all_kinases_information_interacting$kinase_position %in% input$position_selection)
    }
  })
  
  
  
  
  
  
  ##### update the selection options ##########################################
  
  output$kinase_selection = renderUI({ 

    if(nrow(r_upload()) > 0){
      
      choices_list = r_upload() %>%               # grouping of kinases according to kinase family
        filter(r_selection_position() & r_selection_compound()) %>%
        select(gene_name, kinase_family) %>%
        distinct()
    
      choices_list = base::split(x = choices_list$gene_name, f = choices_list$kinase_family)
      choices_list = map(.x = choices_list, .f = base::sort)
      choices_list = map(.x = choices_list, .f = as.list)
      
    }else{
      
      choices_list = NULL
      
    }
    
    selectizeInput(inputId = "kinase_selection",
                   label = "Filter for kinases",
                   choices = choices_list,
                   selected = input$kinase_selection,
                   multiple = TRUE,
                   options = list(placeholder = "all", allowEmptyOption = T, selectOnTab = F)
    )
    
    })



  output$position_selection = renderUI({
    
    if(nrow(r_upload()) > 0){
      
      if((length(input$kinase_selection) == 1) & (length(input$compound_selection))){
        
        choices_list = r_upload() %>%               # grouping of positions according to classification
          filter(r_selection_kinase() & r_selection_compound()) %>%
          select(kinase_position, functional_class) %>%
          distinct()
        
        choices_list = base::split(x = choices_list$kinase_position, f = choices_list$functional_class)
        choices_list = map(.x = choices_list, .f = base::sort)
        choices_list = map(.x = choices_list, .f = as.list)
        
      }else{
        
        choices_list = r_upload() %>%
          filter(r_selection_kinase() & r_selection_compound()) %>% 
          pull(kinase_position)%>%
          unique() %>%
          sort()
        
      } # inner ifelse
      
    }else{
      
      choices_list = NULL
      
    } # outer ifelse
  
    selectizeInput(inputId = "position_selection",
                   label = "Filter for positions",
                   choices = choices_list,
                   selected = input$position_selection,
                   multiple = TRUE,
                   options = list(placeholder = "all", allowEmptyOption = T, selectOnTab = F)
    )

    })


  output$compound_selection = renderUI({
    
    selectizeInput(inputId = "compound_selection",
                   label = "Filter for compounds",
                   choices = r_upload() %>%
                     filter(r_selection_kinase() & r_selection_position()) %>% 
                     pull(compound) %>%
                     unique() %>%
                     str_sort(),
                   selected = input$compound_selection,
                   multiple = TRUE,
                   options = list(placeholder = "all", allowEmptyOption = T, selectOnTab = F)
    )

    })
  
  
  
  r_filtered_table = reactive({  ##### definition of crucial reactive value: filtered table #####
      
      if(!is.null(input$affinity_upload)){
        
        my_data = r_upload()

        
        if(!is.null(input$kinase_selection)){
          
          if(input$kinase_selection[1] != ""){
            
            my_data = filter(my_data, gene_name %in% input$kinase_selection)     #  an empty selection is "", not NULL
            
          }
        } 
        
        if(!is.null(input$position_selection)){
          
          if(input$position_selection[1] != ""){
            
            my_data = filter(my_data, kinase_position %in% input$position_selection)
            
          }
        }
        
        if(!is.null(input$compound_selection)){
          
          if(input$compound_selection[1] != ""){
            
            my_data = filter(my_data, compound %in% input$compound_selection)
            
          }
        }
        
        return(my_data) 
        
      }
    
  })
  

  
  
  observe({   ############################################ update overview table #######################################################################
    
    overview_table = NULL

      if(is.null(r_filtered_table())){
        
        overview_table = session$userData$all_kinases_information_interacting
        
        if(!is.null(input$display_columns_selection)){
          
          overview_table = select(overview_table, input$display_columns_selection)
          
        }
        
      }else{  # filterted table != NULL
        
        if(is.null(input$display_columns_selection)){
          
          overview_table = r_filtered_table()
          
        }else{
          
          overview_table = select(r_filtered_table(), input$display_columns_selection)
          
        }
        
      }
    
    if(!is.null(overview_table)){
      
      if("pKDapp_M" %in% colnames(overview_table))                {overview_table$pKDapp_M                 = round(overview_table$pKDapp_M, digits = 2)}
      
      if("apparent_Kd_nM" %in% colnames(overview_table))          {overview_table$apparent_Kd_nM           = round(overview_table$apparent_Kd_nM, digits = 2)}
      
      if("kinome_wide_conservation" %in% colnames(overview_table)){overview_table$kinome_wide_conservation = round(overview_table$kinome_wide_conservation, digits = 2)}
      
      if("target_wide_conservation" %in% colnames(overview_table)){overview_table$target_wide_conservation = round(overview_table$target_wide_conservation, digits = 2)}
      
      if("conservation_overrepresentation_factor" %in% colnames(overview_table)){overview_table$conservation_overrepresentation_factor = round(overview_table$conservation_overrepresentation_factor, digits = 2)}
      
      if("median_pkDapp_M" %in% colnames(overview_table))         {overview_table$median_pkDapp_M          = round(overview_table$median_pkDapp_M, digits = 2)}
      
      if("weighted_median_pkDapp_M_different_aa" %in% colnames(overview_table)){overview_table$weighted_median_pkDapp_M_different_aa = round(overview_table$weighted_median_pkDapp_M_different_aa, digits = 2)}
      
      if("increased_affinity" %in% colnames(overview_table))      {overview_table$increased_affinity       = round(overview_table$increased_affinity, digits = 2)}
      
      filtered_descriptions = descriptions[which(names(descriptions) %in% colnames(overview_table))]
      
      col_desc = sprintf("var tips = [%s], header = table.columns().header(); for (var i = 0; i < tips.length; i++) { $(header[i]).attr('title', tips[i]); } ",
                         paste(paste0("'", filtered_descriptions, "'"), collapse = " , "))
      
      renaming_key = c(kinase_family = "Kinase family",
                       gene_name = "Kinase",
                       compound = "Compound",
                       pKDapp_M = "Apparent pKD [M]",
                       apparent_Kd_nM = "Apparent KD [nM]",
                       kinase_position = "Position",
                       residue_kinase = "Amino acid",
                       alignment_position = "Alignment position",
                       targetable_inert = "Chemical aa property",
                       functional_class = "Functional class",
                       kinome_wide_conservation = "Kinome-wide conservation [%]",
                       target_wide_conservation = "Target-wide conservation [%]",
                       conservation_overrepresentation_factor = "Conservation overrepresentation",
                       bb = "Backbone interaction observations",
                       sc = "Sidechain interaction observations",
                       bbsc_total_observations = "Total bb/sc observations",
                       median_pkDapp_M = "Median pKDapp same aa [M]" ,
                       weighted_median_pkDapp_M_different_aa = "Median pKDapp diff. aa [M]",
                       increased_affinity = "Delta affinity same/ diff. aa [M]"
      )
      
      names(overview_table) = renaming_key[colnames(overview_table)]
      
      if("Kinase" %in% colnames(overview_table)){overview_table$Kinase = factor(overview_table$Kinase)} # damit kein range slider
      if("Compound" %in% colnames(overview_table)){overview_table$Compound = factor(overview_table$Compound)}
      if("Position" %in% colnames(overview_table)){overview_table$Position = factor(overview_table$Position)}
      if("Amino acid" %in% colnames(overview_table)){overview_table$`Amino acid` = factor(overview_table$`Amino acid`)}
      if("Alignment position" %in% colnames(overview_table)){overview_table$`Alignment position` = factor(overview_table$`Alignment position`)}
      if("Functional class" %in% colnames(overview_table)){overview_table$`Functional class` = factor(overview_table$`Functional class`)}
      if("Chemical aa property" %in% colnames(overview_table)){overview_table$`Chemical aa property` = factor(overview_table$`Chemical aa property`)}
      
      dt_overview_table = datatable(data = overview_table,
                                    rownames = F,
                                    # style = 'bootstrap',
                                    filter = list(position = 'top', clear = TRUE, plain = FALSE),
                                    options = list(pageLength = 25,
                                                   searching = T,
                                                   columnDefs = list(list(className = 'dt-center', targets = "_all"))
                                    ),
                                    callback = JS(col_desc)
      )
      
      output$overview_table = DT::renderDataTable(expr = dt_overview_table)
      
    }
    
  })
  

  observeEvent(eventExpr = r_filtered_table(),  ############### warning message when too many rows for plotting detailed plots #######################
               handlerExpr = {
                
                 if(nrow(r_filtered_table()) <= 30){
                   
                   enable(id = "plot_detailed_plots")
                   
                   hide(id = "too_many_rows_to_plot")
                   
                 }else{
                   
                   disable(id = "plot_detailed_plots")
                   
                   output$too_many_rows_to_plot = renderText(paste("The current selection is too extensive to plot. Please apply
                                                                   more filtering criteria. Max. number of rows: 30. Number of selected rows:",
                                                                   nrow(r_filtered_table())))
                   
                   show(id = "too_many_rows_to_plot")
                   
                   
                   
                 } # else
               } # handlerexpr
               ) # observeEvent
  

  observeEvent(eventExpr = input$plot_detailed_plots,  ################ detailed plots ####################################
               handlerExpr = {
                 
                 plotlist_detailed_plots = make_plotlist_detailed_plots(information_about_all_kinases = session$userData$all_kinases_information_unfiltered,
                                                                        current_filtered_information = r_filtered_table(),
                                                                        kinome_wide_conservation_helper_df = session$userData$kinome_wide_conservation,
                                                                        target_space_conservation_helper_df = session$userData$target_space_conservation)
                
                 dynamic_height = (length(plotlist_detailed_plots)/3)*375
                 
                 output$detailed_plots = renderPlot(

                   # expr = do.call("grid.arrange", c(plotlist_detailed_plots, ncol = 3)),
                   expr = plot_grid(plotlist = plotlist_detailed_plots, ncol = 3,
                                    rel_widths = c(1, 5, 5)),
                   width = 1200,
                   height = dynamic_height,
                   
                   bg = "transparent"

                 ) # render plot
                 
               }
    
  ) # making plots

  
  output$example_input_download = downloadHandler(
    
    filename = "example_affinity_upload.csv",
    
    content = function(file){file.copy('www/AFFINITY_UPLOAD.csv', file)}
    
  )
  
  output$example_input_download_front = downloadHandler(
    
    filename = "example_affinity_upload.csv",
    
    content = function(file){file.copy('www/AFFINITY_UPLOAD.csv', file)}
    
  )

  
  output$download_table = downloadHandler(
    
    filename = function(){paste("Kinase_Residue_Classification_", Sys.Date(), ".csv", sep = "")},
    
    content = function(file){
      
      write.csv(r_filtered_table(), file, row.names = F)
      
    }
    
  )
  
  output$kinase_list_download = downloadHandler(
    
    filename = "Residue_Classification_Kinases.csv",
    
    content = function(file){file.copy('www/Residue_Classification_Kinases.csv', file)}
    
  )
  
  output$structure_list_download = downloadHandler(
    
    filename = "Residue_Classification_PDB_structure_IDs.csv",
    
    content = function(file){file.copy('www/Residue_Classification_PDB_structure_IDs.csv', file)}
    
  )
  
  output$alignment_download = downloadHandler(
    
    filename = "Residue_Classification_Alignment.clustal",
    
    content = function(file){file.copy('www/Clustal_Omega_Alignment_from_Website.clustal', file)}
    
  )
  
  ############################# warning when too many rows selected for direct detailed plots ###############
  
  observeEvent(eventExpr = input$overview_table_rows_selected,
               handlerExpr = {
                 
                 if(length(input$overview_table_rows_selected) <= 30){
                   
                   enable(id = "make_selected_rows_detailed_plots")
                   hide(id = "too_many_rows_selected")
                   
                 }else{
                   
                   disable(id = "make_selected_rows_detailed_plots")
                   
                   output$too_many_rows_selected = renderText(paste("The current selection is too extensive to plot. 
                                                                    For generating detailed plots, you can select a maximum of 
                                                                    30 rows. Number of selected rows:", length(input$overview_table_rows_selected)))
                   
                   show(id = "too_many_rows_selected")
                   
                 } # else
                 
               })
  
  
  ####################### Making detailed plots directly from selected rows ###############################
  
  observeEvent(eventExpr = input$make_selected_rows_detailed_plots,
               handlerExpr = {
                 
                 plotlist_detailed_plots = make_plotlist_detailed_plots(information_about_all_kinases = session$userData$all_kinases_information_interacting,
                                                                        current_filtered_information = r_filtered_table()[sort(input$overview_table_rows_selected),],
                                                                        kinome_wide_conservation_helper_df = session$userData$kinome_wide_conservation,
                                                                        target_space_conservation_helper_df = session$userData$target_space_conservation)

                 dynamic_height = (length(plotlist_detailed_plots)/3)*375

                 output$detailed_plots = renderPlot(
                   
                   expr = plot_grid(plotlist = plotlist_detailed_plots, ncol = 3,
                                    rel_widths = c(1, 5, 5)),
                   width = 1200,
                   height = dynamic_height,
                   
                   bg = "transparent"
                 )
                 
               }
               ) # observeEvent

  
  
  observeEvent(eventExpr = input$overview_plots,     ########### Making the heatmap residue overview plot #################################################
               handlerExpr = {
                 
                 withProgress(message = "Plotting...",
                              expr = {
                                
                                overview_plot_helper_1 = r_filtered_table() %>%
                                  select(gene_name, kinase_position, residue_kinase, functional_class, alignment_position) %>% 
                                  group_by(gene_name, alignment_position, functional_class) %>% 
                                  count() %>% 
                                  spread(key = functional_class,
                                         value = n)
                                
                                if(!("key" %in% colnames(overview_plot_helper_1))){overview_plot_helper_1$key = NA}
                                if(!("selectivity" %in% colnames(overview_plot_helper_1))){overview_plot_helper_1$selectivity = NA}
                                if(!("potency" %in% colnames(overview_plot_helper_1))){overview_plot_helper_1$potency = NA}
                                if(!("scaffold" %in% colnames(overview_plot_helper_1))){overview_plot_helper_1$scaffold = NA}
                                
                                overview_plot_helper_1 = mutate(overview_plot_helper_1, coloring_class = case_when(
                                    !is.na(key) ~ "key",
                                    !is.na(selectivity) ~ "selectivity",
                                    !is.na(potency) ~ "potency",
                                    !is.na(scaffold) ~ "scaffold"
                                  )
                                  )
                                
                                incProgress()
                                
                                overview_plot_helper_2 = select(r_filtered_table(), gene_name, kinase_position, residue_kinase, alignment_position) %>%
                                  distinct()
                                
                                incProgress()
                                
                                overview_plot_input = inner_join(overview_plot_helper_1, overview_plot_helper_2,
                                                                 by = c("gene_name" = "gene_name", "alignment_position" = "alignment_position")) %>% 
                                  unite(col = "display_residue", residue_kinase, kinase_position, sep = "", remove = F)
                                
                                overview_plot_input$alignment_position = as.factor(overview_plot_input$alignment_position)
                                overview_plot_input$gene_name          = as.factor(overview_plot_input$gene_name)
                                
                                incProgress()
                                
                                residue_heatmap = ggplot(data = overview_plot_input) +
                                  geom_raster(mapping = aes(x = alignment_position, y = gene_name, fill = coloring_class)) +
                                  scale_fill_manual(values = c("key" = "#009900",   # green
                                                               "potency" = "#0066ff",  # blue
                                                               "scaffold" = "#8c8c8c",    # grey
                                                               "selectivity" = "#ff6600")  # orange
                                  ) +
                                  scale_y_discrete(name = "Kinases\n") + 
                                  scale_x_discrete(name = "Position in the multiple sequence alignment of the kinase domains of all human kinases (not consecutive)\n",
                                                   position = "top",
                                                   labels = NULL) +
                                  theme_minimal() +
                                  theme(legend.position = "top",
                                        legend.title = element_blank(),
                                        axis.text.x = element_text(angle = 90))
                           
                                incProgress()
                                
                                my_width = length(unique(overview_plot_input$alignment_position))*11.4
                                my_height = nrow(overview_plot_input)*0.7
                                
                                session$userData$heatmap_input <- overview_plot_input
                                
                                incProgress()
                                
                                output$heatmap_plot = renderPlot(residue_heatmap,
                                                                 width = 80 + my_width,    # 3000
                                                                 height = 95 + my_height,   # 2000
                                                                 bg = "transparent"
                                                                 )
                                
                              }) # withProgress
                 
               }) # observeEvent overview Plots
  
  output$hover_info <- renderUI({       ############################### heatmap hover information #####################################################

    if(!is.null(input$heatmap_hover)){
      
      hover <- input$heatmap_hover
      point <- nearPoints(df = session$userData$heatmap_input, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
      if (nrow(point) == 0) return(NULL)
      
      # calculate point position INSIDE the image as percent of total dimensions
      # from left (horizontal) and from top (vertical)
      left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
      top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
      
      # calculate distance from left and bottom side of the picture in pixels
      left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
      top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
      
      # create style property fot tooltip
      # background color is set so tooltip is a bit transparent
      # z-index is set so we are sure are tooltip will be on top
      style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.9); ",
                      "left:", left_px + 2, "px; top:", top_px + 2, "px;")
      
      # actual tooltip created as wellPanel
      wellPanel(
        style = style,
        p(HTML(paste0("<b> Kinase: </b>", point$gene_name, "<br/>",
                      "<b> Residue: </b>", point$display_residue, "<br/>",
                      "<b> Classifications with different drugs: <br/>",
                      "<b> key: </b>", point$key, "<br/>",
                      "<b> potency: </b>", point$potency, "<br/>",
                      "<b> scaffold: </b>", point$scaffold, "<br/>",
                      "<b> selectivity: </b>", point$selectivity)
        ) #paste0
        ) # p()
      ) # wellPanel
      
    } # if(!is.null hover)
    
  })
  
  
} # server

