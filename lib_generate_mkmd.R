require(progress)

generate_mkmd = function(generic_data, XML_directory){
  
  xml_files_directory = XML_directory
  
  mkmd = data.frame()    # mkmd = many kinases many drugs
  # k_ = nrow(generic_data)
  # start_time = Sys.time()
  
  
  pb <- progress_bar$new(format = "  parsing XMLs [:bar] :percent remaining time: :eta", total = nrow(generic_data), clear = F, width= 80)
  
  # add new line \n
  
  
  for(i in 1:nrow(generic_data)){
    
    pb$tick()
    
    PLIP_output = parse_PLIP_XML(
      ligandId = generic_data$ligand_ID[i],
      structureId = generic_data$PDB_ID[i],
      xml_directory = xml_files_directory
    )
    
    PLIP_output$gene_name = generic_data$gene_name[i]
    PLIP_output$kinase_family = generic_data$kinase_family[i]
    
    mkmd = dplyr::bind_rows(mkmd, PLIP_output)  
    
    # end_time = Sys.time()
    # cat("time since start:", end_time - start_time, "\n")
    
  }
  
  return(mkmd)
  
}

