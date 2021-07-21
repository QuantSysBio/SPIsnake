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
