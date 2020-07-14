
## 
## In this script, the input files for the Shiny app aka 'the underlying dataset" are created.
##
## ATTENTION: Don't just run the script as a whole, instead, go line by line and read ALL the 
## comments. There are occasions where you have to run scripts outside of RStudio to
## continue (e.g. running a pyhton script).
##


# PDB_annotated_kinases_ligandreport_csv_filepath:
# PDB > Search > Browse by annotation > Enzyme Classification > Transferases > Transferring phosphorus-containing groups > filter:
# in tab "sequence entities: Organism: Homo sapiens; Experimental Method: X-ray
# (PDB Search Parameter is now: Ligands associated with structures from: 'EnzymeClassificationTree Search for 2.7: Transferring phosphorus-containing groups and Taxonomy is just Homo sapiens (human) and Experimental Method is X-RAY ')
# in Dropdown "Reports:" select "Ligands" > Download csv
PDB_annotated_kinases_ligandreport_csv_filepath = "Z:/users_files/Verena Burger/4_datasets/current/PDB/PDB_annotated_kinases_Ligandreport.csv"

altered_pkinfam_filepath = "Z:/users_files/Verena Burger/4_datasets/current/UniProt/pkinfam_altered.txt"

# The output from parseFlatFile.py. For more details, search this script for the variable name => see below
# kinase_domain_sequences_fasta = "Z:/users_files/Verena Burger/4_datasets/current/UniProt/domains.fasta"






require(tidyverse)
require(msa)
require(seqinr)
require(reshape2)
require(bio3d)

source("lib_PDB_ID_retrieval.R")
source("lib_parser_Uniprot_pkinfam.R")
source("lib_generate_mkmd.R")
source("lib_PLIP_XML_parser.R")


# filter for Ligands of a molecular weight between 300 and 600 Dalton. ------

var <- c("Ligand MW")
cond <- c(300, 600)
df_PDB_structures = readr::read_csv(PDB_annotated_kinases_ligandreport_csv_filepath) %>% 
  dplyr::select("PDB ID", "Ligand ID", "Ligand MW", "Ligand Name") %>% 
  dplyr::distinct() %>% 
  dplyr::filter(.data[[var[[1]]]] > cond[[1]],
                .data[[var[[1]]]] < cond[[2]])
rm(var, cond)

# now, df_PDB_structures is a dataframe which contains the PDB structure- and Ligand IDs from PDB structures
# in which there is a protein which PDB classifies as "Transferase which transfers phopsh.-containing groups"
# from Homo sapiens, generated with X-ray.


# PLIP Download -------------------------------------------------------------

# PLIP_input = df_PDB_structures %>%
#   pull(`PDB ID`) %>%
#   unique()

# write_lines(x = PLIP_input, path = "Z:/users_files/Verena Burger/4_datasets/current/PLIP/PLIP_input_PDB_IDs.txt")
# this is the input for the PLIP download

# Now happens the PLIP Download with the 2 respective scripts for that written by Tobi:
# script_PLIP_download.R
# script_PLIP_download_curl.R

# [Tobi downloads the files for me. (this takes several hours)]


# retrieve gene names ---------------------------------------------------------

pkinfam = data.table::fread(altered_pkinfam_filepath,
                stringsAsFactors = F,
                fill = T,
                header = F)
pkinfam_gene_names = unique(pkinfam$V1)       # all kinase gene names from the pkinfam file

# Now: searching PDB with the pkinfam kinase names as input to get the structure IDs
# This section is taken from pdb2kb > PDB_structure_retieval.R

PDB_IDs_gene_names = data.frame()
k_ = length(pkinfam_gene_names)

start_time = Sys.time()
for (i in 1:k_) {            # takes ca. 6.5 minutes

  cat("PDB ID query is at " , base::round(i / k_, 4) * 100, " %\n")

  df_result = get_PDB_IDs_for_genename(gene_name = pkinfam_gene_names[i])

  if (dim(df_result)[1] > 0) {
    df_result$gene_name = pkinfam_gene_names[i]
    PDB_IDs_gene_names = dplyr::bind_rows(PDB_IDs_gene_names, df_result)

  }else{
    message("Nothing found for ", pkinfam_gene_names[i])
  }
}
end_time = Sys.time()
end_time - start_time
rm(end_time, start_time, k_, i)
# In the above loop, I don"t find structures for a lot of pkinfam gene names, including
# many CDKs (i.e. CDK3) and many MAPKs (i.e. MAP4K5).
# If nothing is found in the loop, it means that there is no PDB entry associated with
# the pkinfam gene name as an entity.


# add the retrieved gene_names to df_PDB_structures ---------------------------

PDB_IDs_gene_names = PDB_IDs_gene_names %>% 
  dplyr::select(PDB_ID, gene_name) %>% 
  dplyr::distinct()

df_PDB_structures = dplyr::inner_join(df_PDB_structures, PDB_IDs_gene_names, by = c("PDB ID" = "PDB_ID"))
# vorher habe ich 4468 rows = 4468 structures, in denen PDB sagt, dass eine Kinase drin ist.
# nach dem inner_join habe ich noch 3462 rows, ich verliere also 1006 structures.
# Theoretisch sollte ich keine rows verlieren. Dass das trotzdem passiert, kann an der unperfekten PDB query liegen.
# -> Das heisst, wahrscheinlich finde ich ueber den gene name nicht immer alle Strukturen


# exclude structures which are associated with 2 kinase names -----------------

two_kinases_checker = df_PDB_structures %>% 
  dplyr::select(`PDB ID`, gene_name) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(`PDB ID`) %>% 
  dplyr::summarise(counts = dplyr::n()) %>% 
  dplyr::filter(counts == 1) %>% 
  dplyr::select(`PDB ID`)

df_PDB_structures = dplyr::inner_join(df_PDB_structures, two_kinases_checker, by = c("PDB ID" = "PDB ID"))
rm(two_kinases_checker)

# add kinase families from pkinfam --------------------------------------------

Uniprot_pkinfam = parse_pkinfam(pkinfam_altered_filepath = altered_pkinfam_filepath)   ## I could combine this with the use of pkinfam weiter oben!!!!!!
Uniprot_pkinfam = Uniprot_pkinfam[, c(1,3)]

df_PDB_structures = dplyr::left_join(df_PDB_structures, Uniprot_pkinfam, by = c("gene_name" = "kinase name"))


# Exclude kinases with less than 4 structures ---------------------------------

# minimal number of structures of a kinase: 4
# (Wie ist das mit mehreren Liganden in einer Struktur? => Die Anzahl der PDB IDs zaehlt, weil ich nicht weiss, ob die Liganden wirklich an unterschiedlichen
# ATP-Bindestellen binden.)
# backup = df_PDB_structures

structure_number_checker = df_PDB_structures %>% 
  dplyr::select(`PDB ID`, gene_name) %>% 
  dplyr::group_by(gene_name) %>% 
  dplyr::count() %>% 
  dplyr::filter(n >= 4) %>% 
  dplyr::select(gene_name)

df_PDB_structures = dplyr::inner_join(df_PDB_structures, structure_number_checker, by = c("gene_name" = "gene_name"))



# Now, df_PDB_structures is the generic PDB data, enriched for gene_name and kinase_family.
# Following, we will add the interaction information from PDB.
# Therefore, the PDB crystal structures are analyzed with the PLIP tool.

# Check for empty PLIP outputs: go to the file folder containing the XML files.
# sort the XML files by size. If there are files of size 0, check if they are 
# really empty and then delete them.

df_PDB_structures = df_PDB_structures %>% 
  dplyr::rename(PDB_ID = `PDB ID`, ligand_ID = `Ligand ID`, ligand_molecular_weight = `Ligand MW`, ligand_name = `Ligand Name`) %>% 
  dplyr::filter(PDB_ID != "5EBZ")   # empty XML file (see OneNote 10.6.)


# parse the PLIP XML output files ---------------------------------------------

# write_rds(df_PDB_structures, "result_2019_09_09_df_PDB_structures_after_row_170")
# df_PDB_structures = read_rds("result_2019_07_17_df_PDB_structures_row_146")

mkmd = generate_mkmd(generic_data = df_PDB_structures,   # takes a long time to run
                     XML_directory = "Z:/users_files/Verena Burger/4_datasets/current/PLIP/Tobi_Download/xml")

mkmd$resnr = as.numeric(mkmd$resnr)

# write_rds(mkmd, "result_2019_08_09_mkmd")


# Make the Multiple Sequence Alignment ----------------------------------------

# Now, I take the human reference proteome download file from Uniprot -> contains
# the sequences for all human proteins and the information, at what positions the
# kinase domain begins and ends.
# I use Mathias' script_parseFlatFile.py to filter the reference proteome for
# the kinases and then extract just the kinase domain sequences.
# The output file is stored at kinase_domain_sequences_fasta.

# Before: I did the alignment in this Skript. But it turns out, that the
# function msaClustalOmega is not deterministic. Therefore, I conduct the
# MSA in the Clustal Omega Website (EMBL EBI), save it as a clustal file and
# parse it here (see below).

# On the website, I set the output format to clustal, perform the alignment, go
# to "Alignment Summary" (or so) and download the clustal file.
# the downloaded file is parsed here:
clustal_omega_alignment = seqinr::read.alignment(file = "Z:/users_files/Verena Burger/4_datasets/current/UniProt/Clustal_Omega_Alignment_from_Website.clustal",
                                                 format = "clustal", forceToLower = F)

matricized_msa = seqinr::as.matrix.alignment(clustal_omega_alignment) %>% 
  as.data.frame() %>% 
  rownames_to_column()

matricized_msa = mutate(matricized_msa, gene_name = gsub(pattern = "^([a-zA-Z0-9-]+?)_[a-zA-Z0-9-]+?_[a-zA-Z0-9-]+?_\\d+?-\\d+?_.+",
                                                         replacement = "\\1",
                                                         x = matricized_msa$rowname,
                                                         perl = T)) %>% 
  mutate(start     = gsub(pattern = "^[a-zA-Z0-9-]+?_[a-zA-Z0-9-]+?_[a-zA-Z0-9-]+?_(\\d+?)-\\d+?_.+",       
                          replacement = "\\1",
                          x = matricized_msa$rowname,
                          perl = T)) %>% 
  mutate(end       = gsub(pattern = "^[a-zA-Z0-9-]+?_[a-zA-Z0-9-]+?_[a-zA-Z0-9-]+?_\\d+?-(\\d+?)_.+",
                          replacement = "\\1",
                          x = matricized_msa$rowname,
                          perl = T)) %>% 
  select(rowname, gene_name, start, end, everything()) %>%   # reorder the columns
  dplyr::rename(sequence_name = rowname)


# Prepare Shiny input files, including matricized_msa and mkmd ----------------

# Be aware: Some kinases in the reference proteome (and therefore in the MSA) have
# more than 1 kinase domain (so far: max. 2). In this case, if you look just at
# the gene_names (as opposed to at the sequence_names), there will be duplicates.

df_alignment = matricized_msa  # the MSA with Sequence- and gene names
# Factors to character
df_alignment[] <- lapply(df_alignment, as.character)
df_alignment$start = as.numeric(df_alignment$start)
df_alignment$end = as.numeric(df_alignment$end)
df_alignment$end = as.numeric(df_alignment$end)

# transform the alignment from a wide to a long format for more efficient computing further down the code
df_alignment_melted = reshape2::melt(df_alignment, id = c("gene_name", "start", "end", "sequence_name")) %>% 
  dplyr::rename(alignment_position = "variable") %>%
  dplyr::rename(residue_kinase = "value")

df_alignment_melted$alignment_position = as.numeric(df_alignment_melted$alignment_position)

# df_kinome_wide_conservation can just be used for calculation in
# inbetween start and end. Otherwise multiple domains can be considered

# add a column with the information, if there is >= 1 kinase with this aa in this alignment
# position in the alignment:

counter = df_alignment_melted %>% 
  select(alignment_position, residue_kinase) %>% 
  distinct()

# 1 = in this alignment position is just 1 case: 1 aa case (100% conservation)
# 2 = in this alignment position are >= 2 cases: >= 2 aa cases and maybe gaps

counter$MSA_aa_cases = 2

for(i in 1:length(unique(counter$alignment_position))){
  
  idx = which(counter$alignment_position == i)

  
  if(length(idx) == 1){
    
    counter[idx, "MSA_aa_cases"] = 1
    
  }
  
}

counter = select(counter, -residue_kinase) %>% 
  distinct()

df_alignment_melted = left_join(df_alignment_melted, counter, by = c("alignment_position" = "alignment_position"))

write_rds(df_alignment_melted, "result_2019_09_27_df_alignment_melted")

# In the Shiny app, mkmd is loaded and named "df_structure_information.




# for preparation of the Kinobead affinity upload file, see separate R script:
# script_prepare_Kinobead_call_set.R

