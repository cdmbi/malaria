
setwd("~/Documents/malaria")
library(rcdk)
library(data.table)


target_data <- fread("ChEMBL_21_MWTlt900_standardized.csv")
query_data <- read.csv("GAMO_PFdata_200115.csv", stringsAsFactors = FALSE)
query_smiles <- query_data$smiles[1:13403]
target_smiles <- target_data$smiles[110001:120000]

sample <- fread("results_10000.csv")

##calculating all 
options(java.parameters = "-Xmx31000m")
target.mols <- parse.smiles(target_smiles)
#library(parallel)
#cl <- makeCluster(24)
#clusterExport(cl = cl, varlist = "target.mols")
#target.fps <- parLapply(cl = cl,
#                        target.mols,
#                        get.fingerprint, type = "circular")
target.fps <- lapply(target.mols, get.fingerprint, type = "circular")

saveRDS(target.fps, "target_fps_120000.Rds")



query_mols <- parse.smiles(query_smiles)[[1]]
target.mols <- parse.smiles(target_smiles)

query.fp <- get.fingerprint(query_mols, type = "circular")
#target.fps <- lapply(target.mols, get.fingerprint, type = "circular")
#tanimoto_similarity <- unlist(lapply(target.fps, distance, fp2 = query.fp,
#                                     method = 'tanimoto'))
#results <- as.data.frame(tanimoto_similarity)


#my_results <- data.frame()
setwd("~/Documents/malaria")
target_data <- readRDS("target_fps_30000.Rds")
query_data <- readRDS("query_fp_GAMPO.Rds")

library(parallel)
library(doSNOW)
cl <- makeCluster(24)
#clusterExport(cl, "target.fps")
registerDoSNOW(cl)
my_results <- list(1:13403)
my_results <- foreach(i = 1:13403, .packages = 'rcdk') %dopar% {
  #query.mol <- parse.smiles(query_smiles)[[i]]
  #target.mols <- parse.smiles(target_smiles)
  #query.fp <- get.fingerprint(query.mol, type = "circular")
  query.fp <- query_data[[i]]
  #target.fps <- lapply(target.mols, get.fingerprint, type = "circular")
  #target.fps <- readRDS("target_fps_10000.Rds")
  target.fps <- target_data
  tanimoto_similarity <- unlist(lapply(target.fps, distance,
                                       fp2 = query.fp,
                                       method = "tanimoto"))
  rm(query.mol)
  rm(target.mols)
  rm(query.fp)
  rm(target.fps)
  my_results[[i]] <- tanimoto_similarity
                      }

my_results_df <- as.data.frame(do.call("rbind", my_results))

write.csv(my_results_df, file = "results_300000.csv")
