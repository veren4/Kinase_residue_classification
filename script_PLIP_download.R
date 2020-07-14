
require(data.table)
PATH_CSV = "/media/kusterlab/users_files/Verena Burger/4_datasets/current/PLIP/PLIP_input_PDB_IDs.txt"    # 11.07.19: I updated this path

dt_verena = fread(PATH_CSV, data.table = T, header = F)

dt = fread("/home/tschmidt/Sync/2_WORK/verena/ids_2.txt", data.table = T, header = F)
dt_verena[,url:=dt$V1]

download.file(url, destfile, method, quiet = FALSE, mode = "w",
              cacheOK = TRUE,
              extra = getOption("download.file.extra"),
              headers = NULL, ...)

download.file(dt[1], "/tmp/test.txt")

f = function(x){
 return(unlist(strsplit(x, "/"))[7]) 
}

dt_verena$id = unlist(lapply(dt_verena$url, f))

dt_verena[url == ""]
dt_verena[,dl_link := paste0("https://projects.biotec.tu-dresden.de/plip-web/plip/download/", dt_verena$id, "?filePath=outputs%2Freport.xml")]
require(progress)
POSITION = "/home/tschmidt/Sync/2_WORK/verena/dls"
total = nrow(dt_verena)
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = total, clear = FALSE, width= 60)
for (i in 1:total) {
  name = paste0(dt_verena[i]$V1,".xml")
  
  pb$tick()
  if (dt_verena[i]$url != ""){
  download.file(dt_verena[i]$dl_link, file.path(POSITION, name), quiet = T)
  }
}
PATH_CSV_OUT = "/media/kusterlab/users_files/Verena Burger/4_datasets/current/PLIP/PDB_IDs_dl.txt"
fwrite(dt_verena, file = PATH_CSV_OUT)

#https://projects.biotec.tu-dresden.de/plip-web/plip/download/ad7fcd33-8ec0-452b-aac0-cc6f358db47c?filePath=outputs%2Freport.xml