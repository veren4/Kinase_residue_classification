
####################
#                  #
#  user interface  #
#                  #
####################

ui <- fluidPage(

  useShinyjs(),
  
  theme = shinytheme("simplex"),
  
  titlePanel(
    
    windowTitle = "Kinase Residue Classification",
    
    div("Kinase Residue Classification",
         actionLink(inputId = "gitlab",
                    label = tags$img(src = "logo_GitLab.svg", width = "40px", height = "40px"),
                    style = "background-size:cover; background-position:center; position:absolute; right:42px;",
                    onclick = "window.open('https://gitlab.lrz.de/proteomics/students/verenaburger/partitioned_residue_classification')")
         
    )
    
    
  ), # titlePanel
  
  sidebarLayout(
    
    sidebarPanel(
      
      tags$head(
        
        tags$style(HTML(".optgroup[data-group=\"key\"] .optgroup-header { color: #009900 !important; }
                         .optgroup[data-group=\"potency\"] .optgroup-header { color: #0066ff !important; }
                         .optgroup[data-group=\"scaffold\"] .optgroup-header { color: #8c8c8c !important; }
                         .optgroup[data-group=\"selectivity\"] .optgroup-header { color: #ff6600 !important; }
                         .optgroup .optgroup-header { color: #909090  !important; }
                        ") # HTML         # unterste Zeile (group header) rot: #D80000
                   ) # tags$style
      ), # tags$head
      
      width = 2,
      
      span(downloadLink(outputId = "example_input_download_front", label = "For an input example, click here!")),
      
      tipify(
        
        fileInput(inputId = "affinity_upload",                              # affinity data upload
                  label = "Upload drug profiles",
                  placeholder = "upload here",
                  multiple = F,
                  accept = ".csv"),
        title = "For input requirements, see Help panel",
        placement = "right"
      ),
      

        uiOutput(outputId = "kinase_selection"),

        uiOutput(outputId = "position_selection"),

        uiOutput(outputId = "compound_selection"),

      
      conditionalPanel(condition = 'input.tabs_panel === "Graphic overview"',
                       
                       # Suppression of the hover error message to the user:
                       tags$style(type="text/css",
                                  ".shiny-output-error { visibility: hidden; }",
                                  ".shiny-output-error:before { visibility: hidden; }"
                       ),
                       
                       disabled(
                         tipify(actionButton(inputId = "overview_plots", label = "Plot selection overview!"),
                                title = "For this plot, it is advisable to apply few or no filters.",
                                placement = "right"
                         )
                         
                       ) # disabled
      ),
      
      conditionalPanel(condition = 'input.tabs_panel === "Detailed plots"',
                       
                       disabled(
                       actionButton(inputId = "plot_detailed_plots",
                                    label = "Plot selection details!")
                       ), br(), br(),
                       
                       hidden(
                         
                         #br(),
                         
                         # fluidRow(column(width = 12,
                                         
                                         # tags$style(".fa-exclamation {color: red;}"),   #  display: inline;
                                         tags$head(tags$style("#too_many_rows_to_plot {color: red;}")),
                                         # icon("exclamation"),
                                         textOutput(outputId = "too_many_rows_to_plot")       
                                         
                         # ) # width
                         # ) # column
                       ) # hidden
                       
      ),
      
      conditionalPanel(condition = 'input.tabs_panel === "Tabular overview"',
                       
                       disabled(
                         
                         checkboxGroupInput(inputId = "display_columns_selection",
                                            label = "Choose columns to display",
                                            choices = list("Kinase family"    = "kinase_family",
                                                           "Kinase"           = "gene_name",
                                                           "Compound"         = "compound",
                                                           "pKDapp [M]"       = "pKDapp_M",
                                                           "Apparent KD [nM]" = "apparent_Kd_nM",
                                                           "Position"         = "kinase_position",
                                                           "Amino acid"       = "residue_kinase",
                                                           "Alignment position" = "alignment_position",
                                                           "Chemical aa property" = "targetable_inert",
                                                           "Functional class" = "functional_class",
                                                           "Kinome-wide conservation" = "kinome_wide_conservation",
                                                           "Target-wide conservation" = "target_wide_conservation",
                                                           "Conservation overrepresentation" = "conservation_overrepresentation_factor",
                                                           "Backbone interaction observations" = "bb",
                                                           "Sidechain interaction observations" = "sc",
                                                           "Total bb/sc observations" = "bbsc_total_observations",
                                                           "Median pKDapp same aa" = "median_pkDapp_M",
                                                           "Median pKDapp diff. aa" = "weighted_median_pkDapp_M_different_aa",
                                                           "Difference in affinity" = "increased_affinity"
                                                           )
                                            
                         ), # checkboxgroupinput
                         
                         br(),
                         
                         downloadButton(outputId = "download_table",
                                      label = "Download the table"
                         ),
                         
                         helpText("The download file contains all available columns (it is independent from the column selection), and the selected filters 
                                  are applied.")
                         
                       ) # disabled
      ) # conditionalPanel
      
      
    ),  # sidePanel
    
    
    mainPanel(
      width=10,
      
      tabsetPanel(
        
        id = 'tabs_panel',
        
        tabPanel(title = "Tabular overview",
                 br(),
                 
                 fluidRow(column(width = 12, align = "left",
                                 "The listed kinase residues are positions which engage in interactions with a ligand in at least one crystal of this kinase in PDB.", br(),
                                 "Therefore, they don't necessarily interact with the respective drug. For further information, see Help panel."
                                 )),
                 br(),
                 
                 tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.selected {background-color: #999999 !important; color: white !important;}')), # I didn't yet manage to change the font color to white
                 # #ff8080
                 # #ff3333
                 DT::dataTableOutput("overview_table"),
                 br(),
                 fluidRow(
                   
                   column(width = 6, align = "center", offset = 3,
                          disabled(
                            actionButton(inputId = "make_selected_rows_detailed_plots",
                                         label = "Make detailed plots of the selected rows in the \"Detailed plots\" panel!")
                                   ) #disabled
                          ) # column
                 ), # fluidrow
                   
                 br(),
                 
                 fluidRow(column(width = 6, align = "center", offset = 3, 
                                 
                                 hidden(
                                   tags$head(tags$style("#too_many_rows_selected {color: red;}")),
                                   textOutput(outputId = "too_many_rows_selected")       
                                 ) # hidden
                                 
                 )), # fluidrow
                 br(),
                 bsCollapse(multiple = T,
                            
                            bsCollapsePanel(title = "How are the residues classified?",
                                            "The classification system used in this application is based on",
                                            a(href = "https://pubs.acs.org/doi/10.1021/acschembio.6b00709", "Heinzlmeir et al. 2016", target = "_blank"),
                                            ":", br(), br(),
                                            tags$img(src = 'my_decision_tree.svg', width = "900", height = "900")),
                            
                            bsCollapsePanel(title = "Which amino acids are in the groups of chemical amino acid properties?",
                                            "Amino acid groups:", br(),
                                            "aliphatic: Alanine (A), Glycine (G), Isoleucine (I), Leucine (L), Proline (P), Valine (V), Methionine (M)", br(),
                                            "aromatic: Phenylalanine (F), Tryptophan (W), Tyrosine (Y)", br(),
                                            "acid base: Aspartic acid (D), Glutamic acid (E), Arginine (R), Histidine (H), Lysine (K)", br(),
                                            "polar and else: Serine (S), Threonine (T), Cysteine (C), Asparagine (N), Glutamine (Q)"),
                            
                            bsCollapsePanel(title = "How is it determined which kinase positions are listed?",
                                            "The underlying dataset contains crystal structures of the kinases. With the online tool",
                                            a(href = "https://plip.biotec.tu-dresden.de/plip-web/plip/index", "PLIP", target = "_blank"),
                                            ", it is determined which amino acid positions of a given kinase interact with a ligand at least once in the 
                                            crystal structures of this kinase. These positions are considered as \"interacting\" from then on. Therefore,
                                            it may occur that a residue is classified as e.g. \"key\" with a compound which doesn't interact with this 
                                            particular kinase residue. However, from the perspective of medicinal chemistry and drug discovery, it is still 
                                            beneficial to know about the potential role of this kinase position (E.g. in the exemple above with a key residue 
                                            classification, it can be deduced that if that position is to be addressed by an inhibitor, it would take the 
                                            role of a key residue.).")
                 ) # bsCollapse
                 
        ),
        
  
        tabPanel(title = "Graphic overview",
                 br(),
                 
                 div(
                   style = "position:relative",
                   plotOutput("heatmap_plot", 
                              hover = hoverOpts("heatmap_hover", delay = 100, delayType = "debounce"),
                              inline = T),
                   uiOutput("hover_info")
                 ),
                 br(),
                 bsCollapse(multiple = T,
                            
                            bsCollapsePanel(title = "If the residue of a kinase is classified differently with different compounds, how is the corresponding tile colored in the overview plot?",
                                            "The tiles are colored according to a hierarchy: If there is at least one key classification for this kinase residue
                                             with this inhibitor, the color is green (key). This continues analogously in the order key > selectivity > potency > scaffold. This 
                                            hierarchy favours classifications interesting for drug discovery."),
                            
                            bsCollapsePanel(title = "What does the x-axis represent?",
                                            "The x-axis shows positions of the alignment of kinase domain protein sequences. Therefore, an x-axis value represents a 
                                            functional position in the kinase domain. The absolute protein sequence positions of the depcted residues can be seen by 
                                            hovering over the tiles.")
                 ) # bsCollapse

        ),
        
        tabPanel(title = "Detailed plots",
                 br(),
                 
                 "The kinase domain sequences were aligned in a Clustal Omega Multiple Sequence Alignment. In the plots, the selected position is the absolute position of the selected kinase (blue),", br(),
                 "while the positions of the compared kinases are their respective alignment positions. This ensures the functional comparability of the kinase positions.", br(),
                 
                 br(),
                 
                 plotOutput("detailed_plots",
                            inline = T),
                 br(),
                 bsCollapse(multiple = T,
                            
                            bsCollapsePanel(title = "What does blue coloring signify?",
                                            "Blue coloring signifies the amino acid in the kinase at this position."),
                            
                            bsCollapsePanel(title = "What is represented by a point in the affinity plot (right)?",
                                            "Each point represents a kinase."),
                            
                            bsCollapsePanel(title = "What does the dashed linie at 6 molar in the detailed affinity plots mean?",
                                            "This line represents an arbitrarily chosen threshold for which interaction can be considered an
                                            actual binding event and which one is random. Values with an affinity lower than 6 M can be
                                            considered random interactions. This threshold has no effect on the assignment of the functional
                                            class and is only there for visual guidance."),
                            
                            bsCollapsePanel(title = "What do the horizontal bars in the affinity plot (right) signify?",
                                            "The short horizontal bars indicate the median of the pkDapp values shown for this amino acid.")
                            ) # bsCollapse
        ),
        
        tabPanel(title = "Help",
                 br(),
                 
                 "The underlying dataset comprises information from", strong("3201 X-ray crystal structures"), "of kinases bound by ligands from PDB.", HTML('&emsp;'), downloadLink(outputId = "structure_list_download", label = list(icon("download"), "Click here to download a list of the PDB IDs of the structures!")), br(), br(),
                 
                 "These structures represent", strong("238 kinases"), "from nine kinase families and", strong("2337 ligands"), ".", HTML('&emsp;'), downloadLink(outputId = "kinase_list_download", label = list(icon("download"), "Click here to download a list of the gene names of the kinases!")), br(), br(),
                 
                 "Furthermore, there are the", strong("protein sequences of 496 kinase domains"), "from the human reference proteome from UniProt.", HTML('&emsp;'), downloadLink(outputId = "alignment_download", label = list(icon("download"), "Click here to download the alignment in clustal format!")), br(), br(),
                 
                 br(),
                
                 bsCollapse(multiple = T,  # FAQ panel
                            
                            bsCollapsePanel(title = "What format does my input dataset have to be in?",
                                            "The affinity data (/ \"drug profiles\") has to be in csv format and have 3 columns: \"compound\", \"gene_name\" and \"apparent_Kd_nM\". 
                                            The information in the \"gene_name\" column has to be the canonical gene names of the kinases, as used e.g. in",
                                            a(href = "https://www.uniprot.org/", "UniProt", target="_blank"), ". The delimiter (digit group separators) between 
                                            the number and the decimal points has to be a dot (.).",
                                            br(), br(),
                                            "Here is an ", downloadLink(outputId = "example_input_download", label = "input example"), "based on",
                                            a(href = "https://science.sciencemag.org/content/358/6367/eaan4368", "Klaeger et al. 2017", target="_blank"), "."),
                            
                            bsCollapsePanel(title = "How are the residues classified?",
                                            "The classification system used in this application is based on",
                                            a(href = "https://pubs.acs.org/doi/10.1021/acschembio.6b00709", "Heinzlmeir et al. 2016", target = "_blank"),
                                            ":", br(), br(),
                                            tags$img(src = 'my_decision_tree.svg', width = "900", height = "900")
                                            ),
                          
                            bsCollapsePanel(title = "Why is a kinase from my uploaded affinity dataset missing in the gene name selection options?",
                                            "If a gene name which is included in your uploaded affinity dataset is missing from the kinase selection options, 
                                            that can be because there is no structure for this kinase in the underlying PDB dataset (more likely) or it can be 
                                            because this gene isn't included in the current reference kinome dataset (pkinfam.txt from Uniprot) (unlikely). If the 
                                            reason is the first one, to solve the issue, there would have to be added an X-Ray crystal structure of this 
                                            kinase bound by an inhibitor to PDB, and the underlying dataset of this app would have to be updated."),
                            
                            bsCollapsePanel(title = "What is an apparent Kd?",
                                            "An apparent Kd is determined to take into account the depletion of protein in pulldown- or titration assays. 
                                            To calculate it, a correction factor (cf) needs to be determined by performing a second pulldown with the 
                                            same lysate and fresh resin. The amount of protein captured in the two pulldowns is then compared/ quantified. 
                                            By multiplying the EC50 value with the correction factor, the EC50 value can be converted into an apparent 
                                            binding constant Kdapp.", br(), br(),
                                            "Kdapp = EC50 * cf", br(), br(),
                                            "Adapted from", a(href = "https://science.sciencemag.org/content/358/6367/eaan4368", "Klaeger et al. 2017", target="_blank"), "."
                                            )
                          
                            ) # bsCollapse
            
        )
        
      )
      
    )
  ),
  
  hr(),
  div(
    width = 10,
    a(href = "https://proteomics.wzw.tum.de/index.php?id=2", "TU Munich Chair of Proteomics and Bioanalytics", target="_blank"),
        
        actionLink(inputId = "WZW",
                   label = tags$img(src = "logo_WZW.svg", width = "30px", height = "30px"),
                   style = "background-size:cover; background-position:center; position:absolute; right:47px;",
                   onclick = "window.open('https://www.wzw.tum.de/index.php?id=2&L=1')")
        
  ),
  
  br()
  
)

