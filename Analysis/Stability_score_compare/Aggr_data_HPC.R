## Aggregate data in HPC parallel results

setwd("~/ConsensusClustering/Data")

k = 6
Kopt = c()
for (sd_max in 1:10)
  for (rep in 0:9){
    File = paste0("temp/g_hc", k, sd_max, rep, "_n600.rds")
    if (file.exists(File)){
      itt = readRDS(File)
      itt["sd_max"] = sd_max
      itt["rep"] = rep
      Kopt = rbind(Kopt, itt)
    }else
      print(File)
  }
Kopt = data.frame(Kopt)
dim(Kopt)

x = Kopt
sort(apply(x, 2, function(x){sum(x==k)}), decreasing = TRUE)


saveRDS(Kopt, file = paste0("results/run3_g_hc", k , "_500.rds"))
