require(data.table)

parse_pkinfam = function(pkinfam_altered_filepath){
  pkinfam = fread(pkinfam_altered_filepath,
                  stringsAsFactors = F,
                  fill = T,
                  header = F)
  
  pkinfam[, "acc human" := paste0(V3,V4)]
  pkinfam[, V3 := NULL]
  pkinfam[, V4 := NULL]
  pkinfam[, "group" := ""]
  pkinfam[, "group name" := ""]
  setcolorder(pkinfam, c("group", "group name", "V1", "V2", "acc human", "V5", "V6"))
  colnames(pkinfam) = c("group", "group name", "kinase name", "entry human","acc human","entry mouse","acc mouse")
  
  
  pkinfam$group[1:58] = "AGC"       # AKT1 - STK38L
  pkinfam$group[59:142] = "CAMK"    # "" - Tssk5
  pkinfam$group[143:154] = "CK1"    # CSNK1A1 - VRK3
  pkinfam$group[155:216] = "CMGC"   # CDK1 - SRPK3
  pkinfam$group[217:227] = "NEK"    # NEK1 - NEK9
  pkinfam$group[228:234] = "RGC"    # GUCY2C - NPR2
  pkinfam$group[235:290] = "STE"    # MAP2K1 - TNIK
  pkinfam$group[291:324] = "TKL"     # ACVR1 - TNNI3K
  pkinfam$group[325:416] = "Tyr"     # AATK - ZAP70
  pkinfam$group[417:496] = "Other"     # AAK1 - WNK4
  
  pkinfam$`group name`[which(pkinfam$group == "AGC")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "CAMK")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "CK1")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "CMGC")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "NEK")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "RGC")] = "kinase: adenylyl cyclase class-4/guanylyl cyclase"
  pkinfam$`group name`[which(pkinfam$group == "STE")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "TKL")] = "Ser/Thr protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "Tyr")] = "protein kinase family"
  pkinfam$`group name`[which(pkinfam$group == "Other")] = ""
  
  
  for(i in 1:nrow(pkinfam)){
    pkinfam$kinase_family[i] = paste(pkinfam$group[i], pkinfam$`group name`[i])
  }

  pkinfam = select(pkinfam, c("kinase name", "entry human", "kinase_family"))
  
  pkinfam = pkinfam[!grepl("MOUSE", pkinfam$`kinase name`),]
  pkinfam = pkinfam[!grepl("MOUSE", pkinfam$`entry human`),]
  
  return(pkinfam)
}
