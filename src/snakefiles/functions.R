### ---------------------------- Proteome pre-processing ----------------------------

splitWithOverlap <- function(seq, max_length, overlap_length) {
  start = seq(1, length(seq), by=max_length-overlap_length)
  end   = start + max_length - 1
  end[end > length(seq)] = length(seq)
  
  out <- lapply(1:length(start), function(i){
    paste(seq[start[i]:end[i]], sep="", collapse="")
  })
  names(out) <- paste(start, end, sep="-")
  return(out)
}

Split_list_max_length <- function(String_list, max_length=2000, overlap_length=MiSl*2){
  out <- lapply(seq_along(String_list), function(i){
    
    # Extract sequence attributes
    protein_name <- attr(String_list[[i]], "name")
    protein_seq <- getSequence(String_list[[i]])
    
    # Split input into overlaping chunks
    protein_chunks <- splitWithOverlap(seq = protein_seq, 
                                       max_length = max_length, 
                                       overlap_length = overlap_length)
    
    # Save chunk properties in sequence names
    names(protein_chunks) <- paste0(protein_name, "|chunk:", names(protein_chunks))
    return(protein_chunks)
  })
  # Convert to incoming list format
  unl <- unlist(out)
  out <- as.list(unl)
  names(out) <- names(unl)
  return(out)
}

Split_list_max_length_parallel <- function(String_list, max_length=2000, overlap_length=MiSl*2){
  out <- mclapply(seq_along(String_list), mc.cores = Ncpu, function(i){
    
    # Extract sequence attributes
    protein_name <- attr(String_list[[i]], "name")
    protein_seq <- getSequence(String_list[[i]])
    
    # Split input into overlaping chunks
    protein_chunks <- splitWithOverlap(seq = protein_seq, 
                                       max_length = max_length, 
                                       overlap_length = overlap_length)
    
    # Save chunk properties in sequence names
    names(protein_chunks) <- paste0(protein_name, "|chunk:", names(protein_chunks))
    return(protein_chunks)
  })
  # Convert to incoming list format
  unl <- unlist(out)
  out <- as.list(unl)
  names(out) <- names(unl)
  return(out)
}

Save_prot_chunk <- function(dat, 
                            nF=nF, 
                            orderedProteomeEntries=orderedProteomeEntries, 
                            directory=directory){
  for(j in 1:nF){
    start=Pi[j]+1
    end=Pi[j+1]
    
    inputSequence = dat[orderedProteomeEntries[start:end]]
    seqinr::write.fasta(sequences = inputSequence, names = names(inputSequence), 
                        file.out = paste0(directory, "/", proteome, "_", start, "_", end, ".fasta"))
  }
}

### ---------------------------- PCP/PSP generation ----------------------------

CutAndPaste_seq <- function(inputSequence,nmer,MiSl){
  
  results <- list()
  peptide = strsplit(inputSequence,"")[[1]]
  L = length(peptide)
  
  if(L>nmer){
    
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L,nmer)
    # print("compute all PCP with length")
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1)==nmer)
    cpNmer = cp[index,]
    # print("got all Nmer PCP")
    
    CPseq = translateCP(cpNmer,peptide)
    # print("CP translated")
    
    if(length(cp[,1])>1){
      
      sp = computeSPcomplete(cp,maxL=nmer,minL=nmer,MiSl=MiSl)
      # print(paste("all SP:",length(sp[,1])))
      
      SPseq = translateSP(sp,peptide)
      # print("All SPs translated")
      
      CPseqClean = unique(CPseq)
      # print(paste("CP without doubles:",length(CPseqClean)))
      
      SPseqClean = unique(SPseq)
      # print(paste("SP without doubles:",length(SPseqClean)))
      
      x = removeCPfromSP_seq(SPseqClean,CPseqClean)
      SPseqClean = x[[1]]
      # print(paste("SP without CP:",length(SPseqClean)))
      
      prot_stats <- data.frame(protein = attr(inputSequence, "name"),
                               all_PCP = length(CPseq),
                               all_PSP = length(SPseq),
                               unique_PCP = length(CPseqClean),
                               unique_PSP = length(SPseqClean),
                               unique_PSP_noPCP = length(SPseqClean)
                               )
    }
    else{
      SPseqClean = as.character("")
      CPseqClean = as.character("")
    }
    results[[1]] = SPseqClean
    results[[2]] = CPseqClean
    results[[3]] = prot_stats
    
  }
  else{
    results[[1]] = as.character("")
    results[[2]] = as.character(inputSequence)
    results[[3]] = data.frame()
    
    if(L==nmer){
      results[[1]] = as.character("")
      results[[2]] = as.character(inputSequence)
      results[[3]] = data.frame()
    }
  }
  return(results)
  rm(results, peptide, L, cp, index, cpNmer, CPseq, sp, SPseq, CPseqClean, SPseqClean, x, prot_stats)
}

# compute all PCP with length <=nmer
computeCPomplete <- function(L,nmer){
  
  maxL = nmer+1
  CP <- numeric()
  
  for(i in 1:L){
    CP = rbind(CP, cbind(rep(i,length(c(i:min(L,(i+maxL-2))))),
                         c(i:min(L,(i+maxL-2)))))
  }
  return(CP)
  rm(maxL, CP)
}

# translate PCP
translateCP <- function(CP,peptide){
  
  CPseq = rep(NA,dim(CP)[1])
  for(i in 1:dim(CP)[1]){
    CPseq[i] = paste(peptide[CP[i,1]:CP[i,2]],sep="",collapse="")
  }
  return(CPseq)
  rm(CPseq)
}

# compute all PSP with length == nmer
computeSPcomplete <- function(cp,maxL,minL,MiSl){
  
  SP <- numeric()
  N = dim(cp)[1]
  NN = 5 * (10**6)
  
  SP = matrix(NA,NN,4)
  
  a = 1
  for(i in 1:N){
    
    temp = cbind(rep(cp[i,1],N),
                 rep(cp[i,2],N),
                 cp[,1],
                 cp[,2])
    
    L = temp[,4]-temp[,3]+temp[,2]-temp[,1]+2
    
    ind = which(((temp[,3]-temp[,2])==1)|(L>maxL)|(L<minL)|((temp[,3]-temp[,2])>(MiSl+1))|((temp[,1]-temp[,4])>(MiSl+1))|((temp[,3]<=temp[,2])&(temp[,4]>=temp[,1])))
    
    if(length(ind)>0){
      temp = temp[-ind,]
    }
    
    if(!is.null(dim(temp)[1])){
      if((a+dim(temp)[1]-1)>NN){
        SP = rbind(SP,matrix(NA,NN,4))
      }
      
      if(dim(temp)[1]>0){
        SP[c(a:(a+dim(temp)[1]-1)),] = temp
      }
      a = a+dim(temp)[1]
    }
  }
  # remove all empty lines from SP
  if(a<NN){
    SP = SP[-c(a:NN),]
  }
  return(SP)
  rm(temp, N, NN, L, ind, SP, a)
}

# translate PSP
translateSP <- function(SP,peptide){
  
  SPseq = rep(NA,dim(SP)[1])
  for(i in 1:dim(SP)[1]){
    SPseq[i] = paste(peptide[c(SP[i,1]:SP[i,2],SP[i,3]:SP[i,4])],sep="",collapse="")
  }
  return(SPseq)
  rm(SPseq)
}

removeCPfromSP_seq <- function(x,z){
  
  result <- list()
  result[[1]] = x[-which(x%in%z)]
  
  return(result)
  rm(result)
}

### ---------------------------- Legacy ----------------------------

#################################################################
# create bash scripts
#################################################################

createBash <- function(nmer,cpuName,cpuNum,start,end,directory){

    # run files

    FILE = paste(directory,"/scripts/",nmer,"mers/bashFiles/runDatabaseComputation",cpuNum,".sh",sep="")

    x1 = paste("cd ",directory,"/scripts/",nmer,"mers/RFiles/",sep="")
    x2 = paste("/cluster/soft/linux64/bin/R CMD BATCH runDatabaseComputation",cpuNum,".r",sep="")
    x3 = paste("/cluster/soft/linux64/bin/R CMD BATCH runMWcomputation",cpuNum,".r",sep="")

    write(x1,file=FILE,ncolumns=1)
    write(x2,file=FILE,ncolumns=1,append=TRUE)
    write(x3,file=FILE,ncolumns=1,append=TRUE)
    
    
    
    # qsub file
    
    FILE = paste(directory,"/scripts/queue/QSUB_CPU.sh",sep="")
    
    x1 = paste("qsub -keo -q long -lnodes=",cpuName," ",directory,"/scripts/",nmer,"mers/bashFiles/runDatabaseComputation",cpuNum,".sh",sep="")
    write(x1,file=FILE,ncolumns=1,append=TRUE)


}


#################################################################
# create R scripts database coputation
#################################################################

createRscriptsDB <- function(nmer,cpuName,cpuNum,start,end,directory,MiSl){

    x = rep(NA,21)
    
    x[1] = "library(seqinr)"
    x[2] = paste("source('",directory,"/TDF/code/loadFunctions.r')",sep="")
    x[3] = paste("load('",directory,"/scripts/others/orderedProteomeEntries.RData')",sep="")
    x[4] = "acc = orderedProteomeEntries"
    x[5] = paste("proteins = read.fasta(file='",directory,"/",DB,"', seqtype='AA',as.string=TRUE)",sep="")
    
    x[6] = paste("start =",start)
    x[7] = paste("end =",end)
    
    x[8] = paste("for(pept in ",start,":",end,"){",sep="")
    x[9] = "index = which(attributes(proteins)$names==acc[pept])"
    x[10] = "inputSequence = proteins[[index]]"
    x[11] = paste("SP1s = CutAndPaste(inputSequence,nmer=",nmer,",MiSl=",MiSl,")",sep="")
    x[12] = "PCP = list()"
    x[13] = "PSP = list()"
    x[14] = "PCP[[1]] = SP1s[[3]]"
    x[15] = "PCP[[2]] = SP1s[[4]]"
    x[16] = "PSP[[1]] = SP1s[[1]]"
    x[17] = "PSP[[2]] = SP1s[[2]]"
    
    f1 = paste(directory,"/databases/",nmer,"mers/splicedPeptides/",sep="")
    x[18] = paste("save(PSP,file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")
    
    f1 = paste(directory,"/databases/",nmer,"mers/nonsplicedPeptides/",sep="")
    x[19] = paste("save(PCP,file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")

    # write entry into status
    FILE = paste(directory,"/status/Nmers/",nmer,"mers/runStatus",cpuNum,".txt",sep="")
    
    x[20] =  paste("write(paste('Database computation ',pept-start+1,'; protein ',pept,' ; ',100*(pept-start+1)/(end-start+1),'%',sep=''),file='",FILE,"',append=TRUE,ncolumns=1)",sep="")

    x[21] = "}"


    FILE = paste(directory,"/scripts/",nmer,"mers/RFiles/runDatabaseComputation",cpuNum,".r",sep="")

    write(x[1],file=FILE,ncolumns=1)
    for(i in 2:21){
        write(x[i],file=FILE,ncolumns=1,append=TRUE)
    }
    


}



#################################################################
# create R scripts MW computation
#################################################################

createRscriptsMW <- function(nmer,cpuName,cpuNum,start,end,directory,MiSl){
    
    x = rep(NA,21)
    

    x[1] = paste("source('",directory,"/TDF/code/loadFunctions.r')",sep="")
    x[2] = paste("load('",directory,"/scripts/others/orderedProteomeEntries.RData')",sep="")
    x[3] = "acc = orderedProteomeEntries"

    x[4] = paste("start =",start)
    x[5] = paste("end =",end)
    x[6] = paste("for(pept in ",start,":",end,"){",sep="")
    
    
    f1 = paste(directory,"/databases/",nmer,"mers/splicedPeptides/",sep="")
    x[7] = paste("load(file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")
    x[8] = "MW = computeMZ(PSP[[2]])"
    f1 = paste(directory,"/databases/",nmer,"mers/MWsplicedPeptides/",sep="")
    x[9] = paste("save(MW,file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")
    
    
    
    f1 = paste(directory,"/databases/",nmer,"mers/nonsplicedPeptides/",sep="")
    x[10] = paste("load(file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")
    x[11] = "MW = computeMZ(PCP[[2]])"
    f1 = paste(directory,"/databases/",nmer,"mers/MWnonsplicedPeptides/",sep="")
    x[12] = paste("save(MW,file=paste('",f1,"',acc[pept],'.RData',sep=''))",sep="")
    
    
    
    # write entry into status
    FILE = paste(directory,"/status/Nmers/",nmer,"mers/runStatus",cpuNum,".txt",sep="")
    
    x[13] =  paste("write(paste('Molecular weight computation ',pept-start+1,'; protein ',pept,' ; ',100*(pept-start+1)/(end-start+1),'%',sep=''),file='",FILE,"',append=TRUE,ncolumns=1)",sep="")
    
    
    x[14] = "}"
    
    
    
    
    FILE = paste(directory,"/scripts/",nmer,"mers/RFiles/runMWComputation",cpuNum,".r",sep="")
    
    write(x[1],file=FILE,ncolumns=1)
    for(i in 2:14){
        write(x[i],file=FILE,ncolumns=1,append=TRUE)
    }
    
    
    
}






#################################################################
# database functions
#################################################################




CutAndPaste <- function(inputSequence,nmer,MiSl){
    
    results <- list()
    peptide = strsplit(inputSequence,"")[[1]]
    L = length(peptide)
    
    if(L>nmer){
        
        # compute all PCP with length <=nmer
        cp = computeCPomplete(L,nmer)
        print("compute all PCP with length")
       
    # get all PCP with length == nmer
        index = which((cp[,2]-cp[,1]+1)==nmer)
        cpNmer = cp[index,]
        print("got all Nmer PCP")

        CPseq = translateCP(cpNmer,peptide)
        print("CP translated")
        
        if(length(cp[,1])>1){
            
            sp = computeSPcomplete(cp,maxL=nmer,minL=nmer,MiSl=MiSl)
            print(paste("all SP:",length(sp[,1])))
            
            SPseq = translateSP(sp,peptide)
            print("All SPs translated")

           # x = removeDoubles(CPseq,cpNmer)
            CPseqClean = CPseq
            cpClean = cpNmer
            print(paste("CP without doubles:",length(CPseqClean)))
            
           # x = removeDoubles(SPseq,sp)
            SPseqClean = unique(SPseq)
            spClean = sp
            print(paste("SP without doubles:",length(SPseqClean)))

            x = removeCPfromSP(SPseqClean,spClean,CPseqClean)
            SPseqClean = x[[1]]
            spClean = x[[2]]
            print(paste("SP without CP:",length(SPseqClean)))
        }
        else{
            spClean = numeric()
            SPseqClean = numeric()
            cpClean = numeric()
            CPseqClean = numeric()
        }
        results[[1]] = spClean
        results[[2]] = SPseqClean
        results[[3]] = cpClean
        results[[4]] = CPseqClean
        
    }
    else{
        results[[1]] = numeric()
        results[[2]] = numeric()
        results[[3]] = numeric()
        results[[4]] = numeric()
        
        if(L==nmer){
            results[[3]] = matrix(c(1,nmer),1,2)
            results[[4]] = inputSequence
        }
    }
    return(results)
}









# remove double sequences
removeDoubles <- function(x,y){
    
    result <- list()
    
    temp = unique(x)
    x2 = rep(NA,length(temp))
    y2 = matrix(NA,length(temp),dim(y)[2])
    for(i in 1:length(temp)){
        ind = which(x==temp[i])
        x2[i] = x[ind[1]]
        y2[i,] = y[ind[1],]
    }
    
    result[[1]] = x2
    result[[2]] = y2
    return(result)
    
}

# remove PCP from PSP list
removeCPfromSP <- function(x,y,z){
  
  result <- list()
  index = which(x%in%z)
  
  result[[1]] = x[-index]
  result[[2]] = y[-index,]
  
  return(result)
}





#################################################################
# MZ computation
#################################################################


computeMZ <- function(seq){
    
    
    if(length(seq)<1){
            mz = numeric()
    }
    if(length(seq)>0){

        mz = rep(NA,length(seq))
        
        for(i in 1:length(seq)){
            
            aa = matrix(NA,2,20)
            aa1 = c("A",    "R",    "N",    "D",    "C",    "E",    "Q",    "G",    "H",    "I",    "L",    "K",    "M",    "F",    "P",    "S",    "T",    "W",    "Y",    "V")
            aa2 = c(71.037114,  156.101111, 114.042927, 115.026943, 103.009185, 129.042593, 128.058578, 57.021464,  137.058912, 113.084064, 113.084064, 128.094963, 131.040485, 147.068414, 97.052764,  87.032028,  101.047679,     186.079313, 163.06332,  99.068414)
            
            t = strsplit(seq[i],"")[[1]]
            MW = 18.01528
            
            for(k in 1:length(t)){
                ind = which(aa1==t[k])
                if(length(ind)>0){
                    MW = MW + aa2[ind]
                }
                else{
                    MW = NA
                }
            }
            
            mz[i] = round(MW, digits = 5)
            
        }
    }
    
    return(mz)
}




