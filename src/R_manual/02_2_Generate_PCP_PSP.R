library(seqinr)
library(parallel)
library(parallelly)
# library(foreach)
source("/home/yhorokh/Snakemake/SPI-snake/src/functions.R")

#################################################################
# (1) read input file and extract info
#################################################################

# Nmers
Nmers = as.numeric(8:15)

# max intervening sequence length
MiSl = as.numeric(25)

# CPUs
# Ncpu = availableCores(27)
Ncpu = availableCores()
cl <- parallel::makeForkCluster()gc()
# cl <- parallelly::makeClusterPSOCK(Ncpu)
cl <- parallelly::autoStopCluster(cl)

# database
dat = read.fasta(file="/home/yhorokh/Snakemake/SPI-snake/results/DB_exhaustive/inputSequence_100422_100551.fasta", 
                 seqtype="AA", as.string = TRUE)
# Output dir
directory = "/home/yhorokh/Snakemake/SPI-snake/results/exhaustive_DB"
suppressWarnings(
  dir.create(directory)
)

#################################################################
# (2) Generate cleavage and spliced peptides
#################################################################

# 
# Save_prot_chunk_PSP <- function(prot, Nmers=Nmers, nF=nF, MiSl=MiSl, directory=directory){
#   for(i in 1:length(Nmers)){
#     for(j in 1:nF){
#       start=Pi[j]+1
#       end=Pi[j+1]
#       for(i in start:end) {
#         pept=i
#         fileNum=j
#         nmer=Nmers[i]
#         
#         index = which(attributes(proteome)$names==orderedProteomeEntries[pept])
#         inputSequence = proteome[[index]]
#         SP1s = CutAndPaste(inputSequence,nmer=nmer,MiSl=MiSl)
#         PCP = list()
#         PSP = list()
#         PCP[[1]] = SP1s[[3]]
#         PCP[[2]] = SP1s[[4]]
#         PSP[[1]] = SP1s[[1]]
#         PSP[[2]] = SP1s[[2]]
#         save(PSP,file=paste(directory, "/splicedPeptides_", nmer, "_", nF,'.RData',sep=''))
#         save(PCP,file=paste(directory, "/nonsplicedPeptides_", nmer,  "_", nF,'.RData',sep=''))
#       }
#     }
#   }
# }

inputSequence = dat[[130]]
nmer=8
MiSl=25


a <- CutAndPaste_seq(dat[[130]],nmer=8,MiSl=25)
a <- dat[20000:20050]
a <- mclapply(a, CutAndPaste_seq, nmer=8,MiSl=25, mc.cores = Ncpu)
