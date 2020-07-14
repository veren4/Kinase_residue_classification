
############
#          #
#  global  #
#          #
############


require(tidyverse)
require(shinythemes)
require(shinyjs)
require(shinyBS)
require(DT)
require(cowplot)
require(tools)         # for checking the file extension of the input file

source("C:/Users/vburger/code/partitioned_residue_classification/lib_Shiny_calculate_information_for_all_kinases.R")
source("C:/Users/vburger/code/partitioned_residue_classification/lib_Shiny_plotlist_detailed_plots.R")



# Load the MSA file -----------------------------------------------------------

# MSA = read_rds("C:/Users/vburger/code/partitioned_residue_classification/result_2019_09_27_df_alignment_melted")

# use Mathias' Alignment as underlying MSA instead of line 24:
MSA = read_rds("Z:/users_files/Verena Burger/14_Mathias_Alignment_as_App_MSA/result_2019_10_04_Mathias_Alignment_finished_adaptation")


# Load the generic structural information -------------------------------------

df_structure_information = read_rds("C:/Users/vburger/code/partitioned_residue_classification/result_2019_08_09_mkmd")
# explanation for column "MSA_aa_cases":
# 1 = in this alignment position is just 1 case: 1 aa case (100% conservation)
# 2 = in this alignment position are >= 2 cases: >= 2 aa cases and maybe gaps


# Load the targetable/ inert helper df ----------------------------------------

df_targetable_inert_helper = read_rds("C:/Users/vburger/code/partitioned_residue_classification/result_2019_09_05_targetable_inert")


# column descriptions for overview table --------------------------------------

descriptions = c(
  
  kinase_family = "The kinase_family as defined by UniProt (pkinfam.txt)",
  
  gene_name = "The gene name, as defined in UniProt",
  
  compound = "Compound name as defined in upload file",
  
  pKDapp_M = "Apparent pKD in Molar",
  
  apparent_Kd_nM = "Apparent KD in nanoMolar",
  
  kinase_position = "Position in the kinase protein sequence",
  
  residue_kinase = "Amino acid",
  
  alignment_position = "Position in the multiple sequence alignment of the kinase domains of the complete human kinome",
  
  targetable_inert = "See Help panel for associated amino acids",
  
  functional_class = "Kinase residue classification based on Heinzlmeir et al. 2016",
  
  kinome_wide_conservation = "Conservation of this amino acid amongst all kinases",
  
  target_wide_conservation = "Conservation of this amino acid amongst the kinases in the target space of the drug",
  
  bb = "Number of observations in crystal structures, in which this residue is enaged in an interaction with the compound with its backbone; depends exclusively on the kinase position",
  
  sc = "Number of observations in crystal structures, in which this residue is enaged in an interaction with the compound with its sidechain; depends exclusively on the kinase position",
  
  bbsc_total_observations = "Sum of all backbone- and sidechain interaction observations; depends exclusively on the kinase position",
  
  median_pkDapp_M = "Median of the affinities of the drug to kinases with the same amino acid in this position in Molar",
  
  weighted_median_pkDapp_M_different_aa = "Weighted median of the affinities of the drug to kinases with a different (or no) amino acid in this position in Molar",
  
  conservation_overrepresentation_factor = "The factor by which this amino acid is overrepresented in the target space of the compound compared to its kinome-wide conservation",
  
  increased_affinity = "The difference between the median apparent pKD of this compound to kinases with the same amino acid in this position and to kinases with a different (or no) amino acid in this position in Molar"
  
  )

# plot placeholders -----------------------------------------------------------

placeholder_detailed_plots = ggplot() +
  annotate("text", x = 4, y = 25, size = 5, label = "Click \"Plot selection details!\" in the left panel to create this plot\nor select rows in the tabular overview and click on the button below the table.", colour = "#909090") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


placeholder_overview_plot = ggplot() +
  annotate("text", x = 4, y = 25, size = 5, label = "Click \"Plot selection overview!\" in the left panel to create this plot.", colour = "#909090") +
  theme_void() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


