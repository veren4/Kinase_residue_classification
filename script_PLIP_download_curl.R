
require(data.table)
require(pbapply)
PATH_CSV = "/media/kusterlab/users_files/Verena Burger/4_datasets/current/PLIP/PDB_IDs_Kinases_for_PLIP.txt"

dt = fread(PATH_CSV, data.table =T, header = F)
id = "4USE"
dl = function(id){
#  id = "5I9Y"
curl =paste0("curl -v \'https://projects.biotec.tu-dresden.de/plip-web/plip/submitCalculationJob\'  -H \'Origin: https://projects.biotec.tu-dresden.de\' -H \'Upgrade-Insecure-Requests: 1\' -H \'DNT: 1\' -H \'Content-Type: multipart/form-data; boundary=----WebKitFormBoundaryaM1s1iyyEBAh9kg4\' -H \'Referer: https://projects.biotec.tu-dresden.de/plip-web/plip/index\' --data-binary $\'------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"select-pdb\"\r\n\r\nby-id\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"pdbId\"\r\n\r\n", 
             id,
"\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"enableReplaceObsolete\"\r\n\r\n1\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"showAdvancedOptions\"\r\n\r\n\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"name\"\r\n\r\n\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"userEmail\"\r\n\r\n\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4\r\nContent-Disposition: form-data; name=\"submit\"\r\n\r\nRun analysis\r\n------WebKitFormBoundaryaM1s1iyyEBAh9kg4--\r\n\' 2>&1 | grep Location")


x = system(curl, intern = TRUE)
Sys.sleep(1)
return(x)
#xx = strsplit(x, "/")
#x = xx[[1]][length(xx[[1]])]
# x = gsub('.{1}$', '', x)

return(x)
}

y = pblapply(dt$V1, FUN = dl)
z =unlist(lapply(y, function(c) length(c) == 0)) 
y[z] =list(c("not_found"))
fwrite(data.table(plip_id = unlist(y)), "~/Sync/2_WORK/verena/ids.txt")
