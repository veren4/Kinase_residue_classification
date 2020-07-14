require(XML)
# require(dplyr)
# require(tibble)
# require(bio3d)
# require(tidyr)
# require(data.table)


parse_PLIP_XML = function(ligandId, structureId, xml_directory){

  XML_filepath = base::paste(xml_directory, "/", structureId, ".xml", sep = "")
  
  XML_file = XML::xmlParse(XML_filepath)
  u = c("//water_bridge[@id]", "//hydrogen_bond[@id]", "//hydrophobic_interaction[@id]", "//salt_bridge[@id]", "//halogen_bond[@id]", "//pi_stack[@id]", "//pi_cation_interaction[@id]")
  w = c("water_bridge", "hydrogen_bond", "hydrophobic_interaction", "salt_bridge", "halogen_bond", "pi_stack", "pi_cation_interaction")
  
  df = data.frame()
  
  for(ia in 1:length(u)){
    nodes = XML::getNodeSet(XML_file, u[ia])
    xml_content = XML::xmlToDataFrame(doc = XML_file, stringsAsFactors = F, nodes = nodes)
    
    if(all(dim(xml_content) == 0)){
      }else{
      cont = dplyr::select(xml_content, resnr, restype, restype_lig, protcoo, ligcoo)
      cont = tibble::add_column(cont, interaction_type = w[ia])
      df = base::rbind(df, cont)
      }
  }
    
  df = dplyr::filter(df, restype_lig == ligandId)
  
  
  if(nrow(df) == 0){
    return(NULL)
  }else{
    
    df$restype = bio3d::aa321(df$restype)
    
    
    ###########################################
    ##  add sc/ bb information with bio3d    ##
    ###########################################
    
    pdb = suppressWarnings(read.pdb(structureId))
    bb = c("N", "CA", "C", "O")     # Soll ich bb lieber ausserhalb der Funktion definieren, und es darunter rm(bb)?
    
    df2 = pdb$atom
    
    df2$x = base::formatC(df2$x, digits = 3, format = "f")         # or: sprintf("%.3f", b)
    df2$y = base::formatC(df2$y, digits = 3, format = "f")
    df2$z = base::formatC(df2$z, digits = 3, format = "f")
    df2 = dplyr::mutate(df2, COO = paste(x,y,z, sep = ""))

    df2 = df2 %>% 
      dplyr::mutate(bb_or_sc = if_else(elety %in% bb, "bb", "sc"))
    
    df3 = dplyr::select(df2, COO, resid, elety, bb_or_sc) %>% 
      dplyr::rename(kinase_resid = resid, kinase_elety = elety)

    df = dplyr::left_join(df, df3, by = c("protcoo" = "COO"))
    
    df4 = dplyr::select(df2, COO, resid, elety) %>% 
      dplyr::rename(drug_resid = resid, drug_elety = elety)
    
    df = dplyr::left_join(df, df4, by = c("ligcoo" = "COO"))
    
    # the columns protcoo, ligcoo, drug_resid can be removed (and added again for debugging)
    df = dplyr::select(df, -protcoo, -ligcoo, -drug_resid, -kinase_resid) %>% 
      dplyr::distinct() %>% 
      dplyr::rename(drug_PDBid = restype_lig)
    
    # reorder the columns
    df = df[c("drug_PDBid",  "restype", "resnr", "interaction_type", "kinase_elety", "bb_or_sc", "drug_elety")]
    
    return(df)
    
   }
  
  }
