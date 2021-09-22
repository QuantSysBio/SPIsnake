### ---------------------------- Utils ---------------------------------------------
# https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

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
                            directory=directory,
                            maxE=maxE){
  for(j in 1:nF){
    start=Pi[j]+1
    end=Pi[j+1]
    
    inputSequence = dat[orderedProteomeEntries[start:end]]
    seqinr::write.fasta(sequences = inputSequence, names = names(inputSequence), 
                        file.out = paste0(directory, "/", maxE, "_", proteome, "_", start, "_", end, ".fasta"))
  }
}


### ---------------------------- PCP/PSP generation ----------------------------

CutAndPaste_seq_return_sp <- function(inputSequence,nmer,MiSl){
  # Makes PCP, cis-PSP, revcis-PCP from a given input
  names(inputSequence)="test"
  
  results <- list()
  peptide = strsplit(inputSequence,"")[[1]]
  #print(peptide)
  L = length(peptide)
  
  if(L>nmer){
    #print("L>nmer")
    
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L,nmer)
    # print("compute all PCP with length")
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1)==nmer)
    cpNmer = cp[index,]
    # print("got all Nmer PCP")
    print("cp, cpNmer")
    #print(cp)
    #print(cpNmer)
    CPseq = translateCP(cpNmer,peptide)
    # print("CP translated")
    
    if(length(cp[,1])>1){
      #print("length(cp[,1])>1")
      
      ## these are the indicese
      sp = computeSPcomplete(cp,maxL=nmer,minL=nmer,MiSl=MiSl)
      return(sp)
    }
  }
}


translateSP_fast <- function(SP,peptide){
  peptide <- paste0(peptide,collapse="")
  peptide <- AAString(peptide)
  ## apply SP (coordinate) to peptie using biostrings
  left <- IRanges(start = SP[,1],end = SP[,2])
  right <- IRanges(start = SP[,3],end = SP[,4])
  rez <- stringi::stri_join(extractAt(peptide,left),extractAt(peptide,right))
  return(rez)
  rm(rez)
}

## this is the fast implementation of Juliane's function
CutAndPaste_seq_from_big_sp_fast <- function(inputSequence,big_sp_input,nmer,MiSl){
  # we take the precalculated sp object for ALL protein lengths, select the subset we need
  # the subset we need is selected via length(peptide) as index
  # big_sp_input then obviously must be ordered from 1:N
  names(inputSequence)="test"
  
  results <- list()
  peptide = strsplit(inputSequence,"")[[1]]
  #print(peptide)
  L = length(peptide)
  
  sp_input <- big_sp_input[[L]]
  rm(big_sp_input)
  
  if(L>nmer){
    #print("L>nmer")
    
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L,nmer)
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1)==nmer)
    cpNmer = cp[index,]
    
    CPseq = translateCP(cpNmer,peptide)
    # print("CP translated")
    
    if(length(cp[,1])>1){
      #print("length(cp[,1])>1")
      
      sp = sp_input
      
      SPseq = translateSP_fast2(sp,peptide)
      #print("All SPs translated")
      
      CPseqClean = unique(CPseq)
      #print(paste("CP without doubles:",length(CPseqClean)))
      
      SPseqClean = unique(SPseq)
      #print(paste("SP without doubles:",length(SPseqClean)))
      
      #print("SPseqClean")
      #print(SPseqClean)
      x = removeCPfromSP_seq(SPseqClean,CPseqClean)
      SPseqClean = x[[1]]
      #print(paste("SP without CP:",length(SPseqClean)))
      
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
    ## NB:this returns PCP and PSP sequences that are unique on the protein level
    
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


CutAndPaste_seq <- function(inputSequence,nmer,MiSl){
  # Makes PCP, cis-PSP, revcis-PCP from a given input
  
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

CutAndPaste_seq_PCP <- function(inputSequence,nmer){
  # Make PCP sequences from a given input
  
  results <- list()
  peptide = strsplit(inputSequence,"")[[1]]
  L = length(peptide)
  
  if(L>nmer){
    
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L,nmer)
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1)==nmer)
    cpNmer = cp[index,]
    
    CPseq = translateCP(cpNmer,peptide)
    
    if(length(cp[,1])>1){
      
      CPseqClean = unique(CPseq)
      
      prot_stats <- data.frame(protein = attr(inputSequence, "name"),
                               all_PCP = length(CPseq),
                               all_PSP = 0,
                               unique_PCP = length(CPseqClean),
                               unique_PSP = 0,
                               unique_PSP_noPCP = 0
      )
    }
    else{
      SPseqClean = as.character("")
      CPseqClean = as.character("")
    }
    results[[1]] = ""
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
  rm(results, peptide, L, cp, index, cpNmer, CPseq, CPseqClean, prot_stats)
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
  ## aded L!=maxL
  SP <- numeric()
  N = dim(cp)[1]
  NN = 5 * (10**6)
  
  SP = matrix(NA,NN,4)
  #print(SP)
  a = 1
  # repeat as many times as you have cp
  for(i in 1:N){
    temp1 <- rep(cp[i,1],N)
    temp2 <- rep(cp[i,2],N)
    temp3 <- cp[,1]
    temp4 <- cp[,2]
    
    L = temp4-temp3+temp2-temp1+2
    
    ind = which(((temp3-temp2)==1)|L!=maxL|((temp3-temp2)>(MiSl+1))|((temp1-temp4)>(MiSl+1))|((temp3<=temp2)&(temp4>=temp1)))
    
    if(length(ind)>0){
      temp1 = temp1[-ind]
      temp2 = temp2[-ind]
      temp3 = temp3[-ind]
      temp4 = temp4[-ind]
    }
    
    #
    if(!is.null(length(temp1))){
      if((a+length(temp1)-1)>NN){
        ## this looks slow but is rarely called so it is fine
        SP = rbind(SP,matrix(NA,NN,4))
      }
      
      if(length(temp1)>0){
        SP[c(a:(a+length(temp1)-1)),1] = temp1
        SP[c(a:(a+length(temp1)-1)),2] = temp2
        SP[c(a:(a+length(temp1)-1)),3] = temp3
        SP[c(a:(a+length(temp1)-1)),4] = temp4
      }
      a = a+length(temp1)
      
      #return("1")
    }
    #
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

### ---------------------------- MW and PTMs ----------------------------

mixANDmatch3 <- function(mzMin, mzMax, MW0){
  ## object recycling is intentional
  X = data.table(a=MW0, b=mzMin, c=mzMax)
  index = which(data.table::inrange(X$a, X$b, X$c, incbounds=TRUE))
  return(index)
}

computeMZ_biostrings <- function(seq){
  seq <- AAStringSet(seq)
 
  if(length(seq)<1){
    mz = numeric()
  }
  if(length(seq)>0){
   
    mz = rep(NA,length(seq))
   
    aa = matrix(NA,2,20)
      aa1 = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    aa2 = c(71.037114, 156.101111, 114.042927, 115.026943, 103.009185, 129.042593, 128.058578, 57.021464, 137.058912, 113.084064, 113.084064, 128.094963, 131.040485, 147.068414, 97.052764, 87.032028, 101.047679, 186.079313, 163.06332, 99.068414)
   
   
    aa3 <-  matrix(NA, nrow=20, ncol=length(seq))
    for(n in 1:length(aa1)){
      count <- Biostrings::vcountPattern(aa1[n],seq)
      aa3[n,] <- count*aa2[n]
    }
    MW=colSums(aa3)+18.01528
    ## add NA for zero values
    MW=round(MW, digits = 5)
    MW[MW==18.01528] <- NA
  }
  return(MW)
}
