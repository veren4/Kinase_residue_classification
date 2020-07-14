
# require(data.table)
require(stringr)
require(RCurl)

get_PDB_IDs_for_genename = function(gene_name) {
  url1 <- 'http://www.rcsb.org/pdb/rest/search'
  
  xml_base <- '<orgPdbQuery>
  <queryType>org.pdb.query.simple.UniprotGeneNameQuery</queryType>
  <description>UniProt Gene Name:  HBB1</description>
  <query>${gene_name}</query>
  </orgPdbQuery>'
  xml_text = stringr::str_interp(xml_base)           # xml-Text zusammensetzen (interpolieren) (= Genename einsetzen)
  
  h = basicTextGatherer()
  httpheader = c(Accept = "*/*",                         # general expression for POSTING data with R
                 "Content-Type" = "application/x-www-form-urlencoded")
  
  result <- RCurl::curlPerform(
    url = url1,
    httpheader = httpheader,
    postfields = xml_text,
    writefunction = h$update,
    verbose = F
  )
  
  if (h$value() == "") {
    return(data.frame())
  } else {
    df = data.table::fread(h$value(), header = F, sep = ":") %>% 
      dplyr::rename(PDB_ID = V1, entity_number = V2)
    
    # base::colnames(df) = c("PDB_ID", "entity_number")
    return(df)
  }
  
}