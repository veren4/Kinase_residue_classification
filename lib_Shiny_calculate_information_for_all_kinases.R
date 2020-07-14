

calculate_information_for_all_kinases = function(multiple_sequence_alignment,
                                                 generic_structural_information,
                                                 targetable_inert_helper_df,
                                                 uploaded_affinities){
  
  
  withProgress(message = "Processing affinity input...",
               detail = "Depending on the extent of the input, this may take several minutes.",
               min = 0,
               max = 2,
               expr = {
                 
                 df_alignment_melted        = multiple_sequence_alignment
                 
                 df_structure_information   = generic_structural_information
                 
                 df_targetable_inert_helper = targetable_inert_helper_df
                 
                 df_affinity                = uploaded_affinities
                 
                 
                 # Create the dataframe we'll keep adding information to: df_final -------------
                 
                 df_final = df_affinity %>%                                 # wir gehen aus vom hochgeladenen user-affinity-file
                   dplyr::inner_join(df_alignment_melted, by = "gene_name") %>%    # das MSA dranjoinen.
                   # erledigt: checken, ob, wenn 1 kinase 2 kindomaenen hat, diese Zeile dann 2x drin ist mit jeweils einer
                   # der Alignment-Zeilen
                   dplyr::group_by(sequence_name) %>%
                   dplyr::mutate(kinase_position = start) %>%                      # diese + folgende Zeilen: add absolute kinase position
                   dplyr::ungroup() %>%
                   dplyr::mutate(help = ifelse(residue_kinase == "-", 0, 1)) %>%
                   dplyr::group_by(sequence_name, compound) %>% 
                   dplyr::mutate(help2 = cumsum(help)) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(kinase_position = kinase_position + help2 - 1) %>%
                   dplyr::select(-help2, -help) %>%
                   dplyr::mutate(kinase_position = ifelse(residue_kinase == "-", -1, kinase_position))
                 
                 
                 # help df for kw conservation -------------------------------------------------
                 incProgress()
                 df_kinome_wide_conservation = df_alignment_melted %>%
                   dplyr::select(alignment_position, residue_kinase) %>% 
                   dplyr::group_by(alignment_position, residue_kinase) %>% 
                   dplyr::count() %>% 
                   dplyr::rename("kinome_wide_conservation" = "n") %>%
                   dplyr::ungroup() %>% 
                   dplyr::group_by(alignment_position) %>%
                   dplyr::mutate(kinome_wide_conservation = kinome_wide_conservation/ sum(kinome_wide_conservation) * 100)
                 
                 
                 # join the kw_conservation help df to df_final --------------------------------
                 incProgress()
                 df_final = df_final %>% dplyr::inner_join(df_kinome_wide_conservation, by = c("alignment_position", "residue_kinase"))
                 
                 
                 # help df for ts conservation -------------------------------------------------
                 incProgress()
                 df_target_space_conservation = df_final %>%
                   dplyr::select(compound, gene_name, alignment_position, residue_kinase) %>%
                   dplyr::distinct() %>%
                   dplyr::group_by(compound, alignment_position, residue_kinase) %>%
                   dplyr::count() %>% 
                   dplyr::rename("target_wide_conservation" = "n") %>%
                   dplyr::ungroup() %>%
                   dplyr::group_by(compound, alignment_position) %>%
                   dplyr::mutate(target_wide_conservation = target_wide_conservation/ sum(target_wide_conservation) * 100)
                 
                 
                 # join the help df for ts conservation to df_final ----------------------------
                 incProgress()
                 df_final = df_final %>% dplyr::left_join(df_target_space_conservation, by =c("compound", "alignment_position", "residue_kinase"))
                 
                 
                 # Calculation of the 2 affinity columns ----------------------------------
                 
                 # the affinity columns are TARGET-SPACE-specific! (aka drug specific)
                 incProgress()
                 df_affinity_median = df_final %>%                                                 # I will use this later as the source for affinity same aa, and I
                   dplyr::select(compound, alignment_position, pKDapp_M, residue_kinase) %>%       # also use it in the subsequent command as a helper
                   dplyr::group_by(compound, alignment_position, residue_kinase)  %>%
                   dplyr::summarise(median_pkDapp_M = median(pKDapp_M), number_of_points = dplyr::n()) 
                 incProgress()
                 df_aff_diff_aa = df_affinity_median %>%                                # calculating the affinity of the drug to kinases which have a diff. aa in this place
                   dplyr::inner_join(df_affinity_median, by = c("compound" = "compound",  "alignment_position" = "alignment_position")) %>%
                   dplyr::filter(residue_kinase.x != residue_kinase.y) %>%
                   dplyr::mutate(weighted_median_pkDapp_M.y_different_aa = median_pkDapp_M.y * number_of_points.y) %>%
                   dplyr::group_by(compound, alignment_position, residue_kinase.x) %>%
                   dplyr::summarise(weighted_median_pkDapp_M_different_aa = sum(weighted_median_pkDapp_M.y_different_aa) / sum(number_of_points.y)) %>%
                   dplyr::rename("residue_kinase" = "residue_kinase.x")
                 
                 
                 # Calculation of the backbone/ sidechain information --------------------------
                 
                 ###### BE AWARE - INNER JOIN WILL LOOSE US ANYTHING WHICH IS NOT INTERACTING   (which is what we want :) )
                 
                 ###
                 ### Folgenden Abschnitt ggf. anpassen.
                 ### 
                 incProgress()
                 df_bbsc_helper1 = df_final %>%
                   dplyr::select(gene_name, kinase_position, residue_kinase, alignment_position) %>% # Kinase_position muss mit
                   # drinbleiben, weil ich es fuer das joinen mit df_structure_info brauche
                   dplyr::distinct() %>%
                   dplyr::inner_join(df_structure_information, by = c("gene_name" = "gene_name", "kinase_position" = "resnr", "residue_kinase" = "restype")) %>%
                   # I have to join on kinase_position, as I don't have the 
                   # alignment_position in df_structure_info.
                   # There is a column called alignment_position in df_structure_inf., but
                   # that is old and probably wrong. -> Don't use it!
                   
                   dplyr::group_by(gene_name, alignment_position, residue_kinase, bb_or_sc) %>%  # I changed this to alignment_position here
                   dplyr::summarise(n = n())
                 incProgress()
                 df_bbsc_helper2 =  df_final %>%
                   dplyr::filter(kinase_position != -1) %>%     # helper2 habe ich, weil ich diese bb/sc nur wissen muss, wenn die
                   # kinase Position auch besetzt ist, wenn sie also nicht -1 ist
                   dplyr::select(gene_name, kinase_position, residue_kinase, alignment_position) %>%
                   dplyr::distinct()
                 # Ich denke, dass ich die kinase_position hier eigentlich nicht mitnehmen muss. Ich lasse sie jetzt aber mal drin,
                 # weil es so funktioniert und korrekt zaehlt.
                 incProgress()
                 df_bb_sc = df_bbsc_helper1 %>%
                   dplyr::inner_join(df_bbsc_helper2, by = c("gene_name" = "gene_name", "residue_kinase" = "residue_kinase", "alignment_position" = "alignment_position")) %>%
                   # Jetzt habe ich nur noch besetzte kinasepositionen drin.
                   
                   dplyr::group_by(residue_kinase, bb_or_sc, alignment_position) %>%
                   dplyr::summarise(n = sum(n)) %>% 
                   tidyr::spread(key = "bb_or_sc", value = "n") %>%   # hier kommt noch eine NA Spalte dazu, da NA von spread neben bb und
                   # sc als 3. Fall gesehen wird. Diese Zahlen wollen wir ignorieren, deshalb select.
                   dplyr::select(residue_kinase, alignment_position, bb, sc)                                            
                 
                 df_bb_sc$bb[is.na(df_bb_sc$bb)] = 0
                 df_bb_sc$sc[is.na(df_bb_sc$sc)] = 0
                 incProgress()
                 df_bb_sc = dplyr::mutate(df_bb_sc, bbsc_total_observations = bb + sc)
                 
                 
                 # join bbsc to df_final -------------------------------------------------------
                 incProgress()
                 df_final = df_final %>% dplyr::left_join(df_bb_sc, by = c("residue_kinase", "alignment_position"))
                 
                 
                 # join same_aa_affinity to df_final ---------------------------------------
                 incProgress()
                 df_final = df_final %>% dplyr::inner_join(df_affinity_median, by = c("compound", "alignment_position", "residue_kinase")) %>%
                   dplyr::select(-number_of_points)
                 
                 
                 # join diff_aa_affinity to df_final -------------------------------------------
                 incProgress()
                 df_final = df_final %>% dplyr::left_join(df_aff_diff_aa, by = c("compound", "alignment_position", "residue_kinase"))
                 

                 # mark interacting residues + kinase_family -----------------------------------
                 incProgress()
                 df_interacting_residues = df_structure_information %>% 
                   dplyr::select(kinase_family, gene_name, restype, resnr) %>% 
                   dplyr::distinct() %>% 
                   tibble::add_column(interacting_residue = 1)
                 incProgress()
                 df_final = dplyr::left_join(df_final, df_interacting_residues, by = c("gene_name" = "gene_name", "residue_kinase" = "restype", "kinase_position" = "resnr"))
                 
                 
                 # add targetable/ inert -------------------------------------------------------
                 incProgress()
                 df_final = dplyr::left_join(df_final, df_targetable_inert_helper, by = c("residue_kinase" = "aa")) # Performance-Optimierungsmoeglichkeit: 
                 # Evtl. kann ich das schon am Anfang joinen, so dass es nicht in der
                 # App gemacht werden muss.

                 
                 # set the rows where I have no affinity information for the drug with kinases
                 # with a different aa in this position to 0:
                 df_final[which(is.na(df_final$weighted_median_pkDapp_M_different_aa)), "weighted_median_pkDapp_M_different_aa"] = 0
                
                
                 # add overrepresentation factor + incr. affinity delta ------------------------
                 incProgress()
                 df_final = dplyr::mutate(df_final,
                                          conservation_overrepresentation_factor = if_else(condition = interacting_residue == 1,   # calculation necessary only for the interacting residues
                                                                                             true = target_wide_conservation / kinome_wide_conservation,
                                                                                             false = 0),
                                          increased_affinity = if_else(condition = interacting_residue == 1,
                                                                         true = median_pkDapp_M - weighted_median_pkDapp_M_different_aa,
                                                                         false = 0))
                 
                 
                 # add classification ----------------------------------------------------------
                 # add filter: classify only the interacting residues
                 df_final = dplyr::mutate(df_final,
                                          functional_class = if_else(condition = interacting_residue == 1,
                                                                     true = case_when(
                                                                                      MSA_aa_cases!=1 &  conservation_overrepresentation_factor>1.5 & increased_affinity>=0.5 ~ "key",
                                                                                      sc==0 | targetable_inert=="aliphatic"     ~ "scaffold",
                                                                                      is.na(sc) | is.na(targetable_inert)       ~ "-",
                                                                                      kinome_wide_conservation>=50              ~ "potency",
                                                                                      kinome_wide_conservation<50               ~ "selectivity"),
                                                                     false = "-")
                 )
                 
                 
                 result = list(kinase_information_matrix = df_final,
                               kw_conservation_helper = df_kinome_wide_conservation,
                               ts_conservation_helper = df_target_space_conservation)
                 
               } # with Progress expr
               
  ) # with progress
  
  
  
  return(result)
  
}
