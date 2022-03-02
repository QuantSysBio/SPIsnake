#################################################################
# database functions
#################################################################

### ---------------------------- Utils ---------------------------------------------
# https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stop_quietly <- function() {
  opt = options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

seq_list_to_dt <- function(seq_list){
  rbindlist(lapply(seq_list, as.data.table), idcol = "id")
}

splitWithOverlap <- function(seq, max_length, overlap_length) {
  start = seq(1, length(seq), by=max_length-overlap_length)
  end   = start + max_length - 1
  end[end > length(seq)] = length(seq)
  
  out <- lapply(1:length(start), function(i){
    stri_c(seq[start[i]:end[i]], sep="", collapse="")
  })
  names(out) = stri_c(start, end, sep="-")
  return(out)
}

Split_list_max_length <- function(String_list, max_length=2000, overlap_length=MiSl*2){
  out = lapply(seq_along(String_list), function(i){
    
    # Extract sequence attributes
    protein_name = attr(String_list[[i]], "name")
    protein_seq = getSequence(String_list[[i]])
    
    # Split input into overlaping chunks
    protein_chunks = splitWithOverlap(seq = protein_seq, 
                                      max_length = max_length, 
                                      overlap_length = overlap_length)
    
    # Save chunk properties in sequence names
    names(protein_chunks) = paste0(protein_name, "|chunk:", names(protein_chunks))
    return(protein_chunks)
  })
  # Convert to incoming list format
  unl = unlist(out)
  out = as.list(unl)
  names(out) = names(unl)
  return(out)
}

Split_list_max_length_parallel <- function(String_list, max_length=2000, overlap_length=MiSl*2){
  out = bettermc::mclapply(seq_along(String_list), mc.cores = Ncpu, function(i){
    
    # Extract sequence attributes
    protein_name = attr(String_list[[i]], "name")
    protein_seq = getSequence(String_list[[i]])
    
    # Split input into overlaping chunks
    protein_chunks = splitWithOverlap(seq = protein_seq, 
                                      max_length = max_length, 
                                      overlap_length = overlap_length)
    
    # Save chunk properties in sequence names
    names(protein_chunks) = paste0(protein_name, "|chunk:", names(protein_chunks))
    return(protein_chunks)
  })
  # Convert to incoming list format
  unl = unlist(out)
  out = as.list(unl)
  names(out) = names(unl)
  return(out)
}

Split_max_length2 <- function(XStringSet, max_length=500, overlap_length=MiSl*2){
  out <- as.list(seq_along(XStringSet))
  for (i in seq_along(XStringSet)) {
    # Extract coordinates
    at <- breakInChunks(totalsize = width(XStringSet[i]), chunksize = max_length)
    at <- successiveIRanges(width = width(at), gapwidth = -overlap_length)
    
    # Extract sequence
    protein_chunks <- unlist(extractAt(XStringSet[i], at))
    
    # Update names
    names(protein_chunks) <- paste0(names(XStringSet[i]), "|chunk:", at@start, "-", at@start + at@width - 1)
    out[[i]] <- protein_chunks
  }
  return(unlist(AAStringSetList(out)))
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
            SPseqClean = SPseq
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


# compute all PCP with length <=nmer
computeCPomplete <- function(L,nmer){
    
    maxL = nmer+1
    CP <- numeric()
    
    for(i in 1:L){
        temp2 = c(i:min(L,(i+maxL-2)))
        temp1 = rep(i,length(temp2))
        temp3 = cbind(temp1,temp2)
        CP = rbind(CP,temp3)
    }
    
    return(CP)
}


# translate PCP
translateCP <- function(CP,peptide){
    
    CPseq = rep(NA,dim(CP)[1])
    for(i in 1:dim(CP)[1]){
        CPseq[i] = paste(peptide[CP[i,1]:CP[i,2]],sep="",collapse="")
    }
    return(CPseq)
    
}


# compute all PSP with length == nmer
computeSPcomplete <- function(cp,maxL,minL,MiSl){
    
    SP <- numeric()
    N = dim(cp)[1]
    NN = 10**8
    
    SP = matrix(NA,NN,4)
    
    a = 1
    for(i in 1:N){
        
    print(i)
        temp1 = rep(cp[i,1],N)
        temp2 = rep(cp[i,2],N)
        temp3 = cp[,1]
        temp4 = cp[,2]
        
        temp = cbind(temp1,temp2,temp3,temp4)
        L = temp[,4]-temp[,3]+temp[,2]-temp[,1]+2

        ind = which(((temp[,3]-temp[,2])==1)|(L>maxL)|(L<minL)|((temp[,3]-temp[,2])>(MiSl+1))|((temp[,1]-temp[,4])>(MiSl+1))|((temp[,3]<=temp[,2])&(temp[,4]>=temp[,1])))

        if(length(ind)>0){
            temp = temp[-ind,]
        }

        if((a+dim(temp)[1]-1)>NN){
            SP = rbind(SP,matrix(NA,NN,4))
        }

        if(dim(temp)[1]>0){
            SP[c(a:(a+dim(temp)[1]-1)),] = temp
        }

        a = a+dim(temp)[1]
        
        
    }
    # remove all empty lines from SP
    if(a<NN){
        SP = SP[-c(a:NN),]
    }
    

    
    return(SP)
}


# translate PSP
translateSP <- function(SP,peptide){
    
    SPseq = rep(NA,dim(SP)[1])
    for(i in 1:dim(SP)[1]){
        SPseq[i] = paste(peptide[c(SP[i,1]:SP[i,2],SP[i,3]:SP[i,4])],sep="",collapse="")
    }
    
    return(SPseq)
}

removeCPfromSP_seq <- function(x,z){
  result = list()
  result[[1]] = x[-which(x%in%z)]
  
  return(result)
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



translateSP_fast <- function(SP,peptide){
  peptide = stri_c(peptide,collapse="")
  peptide = AAString(peptide)
  ## apply SP (coordinate) to peptie using biostrings
  left = IRanges(start = SP[,1],end = SP[,2])
  right = IRanges(start = SP[,3],end = SP[,4])
  rez = stri_join(extractAt(peptide,left), extractAt(peptide,right))
  return(rez)
}


Generate_PSP <- function(protein_inputs, nmer, MiSl){
  inputSequence <- protein_inputs[[1]]
  peptide = strsplit(inputSequence,"")[[1]]
  L = length(peptide)
  
  results = vector(mode = "list", length = 3)
  
  if(L>nmer){
    # compute all PCP with length <=nmer
    cp = computeCPomplete(L,nmer)
    
    # get all PCP with length == nmer
    index = which((cp[,2]-cp[,1]+1)==nmer)
    cpNmer = cp[index,]
    
    CPseq = translateCP(cpNmer,peptide)
    # print("CP translated")
    
    if(length(cp[,1])>1){
      #print("length(cp[,1])>1")
      sp = protein_inputs[[2]]
      
      SPseq = translateSP_fast(sp,peptide)
      #print("All SPs translated")
      
      CPseqClean = unique(CPseq)
      #print(paste("CP without doubles:",length(CPseqClean)))
      
      SPseqClean = unique(SPseq)
      l_SPseqClean = length(SPseqClean)
      #print(paste("SP without doubles:",length(SPseqClean)))
      #print("SPseqClean")
      #print(SPseqClean)
      #x = removeCPfromSP_seq(SPseqClean,CPseqClean)
      #SPseqClean = x[[1]]
      #print(paste("SP without CP:",length(SPseqClean)))
      
      prot_stats = data.frame(protein = attr(inputSequence, "name"),
                              all_PCP = length(CPseq),
                              all_PSP = length(SPseq),
                              unique_PCP = length(CPseqClean),
                              unique_PSP = l_SPseqClean,
                              unique_PSP_noPCP = length(SPseqClean)
      )
    } else {
      SPseqClean = as.character("")
      CPseqClean = as.character("")
    }
    results[[1]] = SPseqClean
    results[[2]] = CPseqClean
    results[[3]] = prot_stats
  } else {
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
}

CutAndPaste_seq <- function(inputSequence,nmer,MiSl){
  # Makes PCP, cis-PSP, revcis-PCP from a given input
  
  results = list()
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
      
      prot_stats = data.frame(protein = attr(inputSequence, "name"),
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

CutAndPaste_seq_PCP <- function(inputSequence, nmer){
  # Make PCP sequences from a given input
  results = list()
  L = nchar(inputSequence)
  
  if(L>nmer){
    # compute all PCP with length <=nmer
    cp <- computeCPomplete(L, nmer)
    
    # # get all PCP with length == nmer
    # index = which((cp[,2]-cp[,1]+1)==nmer)
    # cpNmer = cp[index,]
    CPseq <- str_sub(inputSequence, cp[,1], cp[,2])
    
    if(nrow(cp) > 1){
      
      CPseqClean = unique(CPseq)
      prot_stats = data.frame(protein = attr(inputSequence, "name"),
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


getPTMcombinations_fast <- function(input,NmaxMod,mods_input=mods){
  s = as.vector(input[1])
  m = as.numeric(as.vector(input[2]))
  
  aa = strsplit(s,split="")[[1]]
  
  kn = which(mods_input$Site=="N-term" & mods_input$Position=="Any N-term")
  kc = which(mods_input$Site=="C-term" & mods_input$Position=="Any C-term")
  
  ka = list()
  for(i in 1:length(aa)){
    ka[[i]] = which(mods_input$Site%in%aa[i] & mods_input$Position=="Anywhere")
  }
  ka = unlist(ka)
  
  modIndex = c(kn,kc,ka)
  modId = mods_input$Id[modIndex]
  modSite = mods_input$Site[modIndex]
  modPos = mods_input$Position[modIndex]
  modDelta = as.numeric(as.vector(mods_input$MonoMass[modIndex]))
  
  # get all combinations of NmaxMod modifications
  combi = list()
  combi[[1]] = c(1:length(modIndex))
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod,length(modIndex))){
      combi[[i]] <- arrangements::combinations(c(1:length(modIndex)),i)
    }
  }
  
  deltaMass = list()
  IDs = list()
  deltaMass[[1]] = modDelta[combi[[1]]]
  IDs[[1]] = modId[combi[[1]]]
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod,length(modIndex))){
      
      matx <- matrix(modDelta[combi[[i]]],dim(combi[[i]])[1],i)
      deltaMass[[i]] <- rowSums(matx)
      
      IDs[[i]] = matrix(modId[combi[[i]]],dim(combi[[i]])[1],i)
      IDs[[i]] <- do.call(stringi::stri_join, c(input=as.data.frame(IDs[[i]], stringsAsFactors = F), sep=";"))
    }
  }
  
  # generate final mod sequences with delta Masses
  peptide = rep(s,sum(unlist(lapply(deltaMass,length)))) ## length of the peptide and
  ids = unlist(IDs)
  deltaMW = unlist(deltaMass)
  finalMW = m+deltaMW
  
  res = data.table(peptide,ids,MW=finalMW)
  return(res)
}

## this is the fast implementation of Juliane's function
CutAndPaste_seq_from_big_sp_fast <- function(inputSequence,big_sp_input,nmer,MiSl){
  ## we take the precalculated sp object for ALL protein lengths, select the subset we need
  ## the subset we need is selected via length(peptide) as index
  ## big_sp_input then obviously must be ordered from 1:N
  # names(inputSequence)="test"
  
  results = list()
  peptide = strsplit(inputSequence,"")[[1]]
  #print(peptide)
  L = length(peptide)
  
  sp_input = big_sp_input[[L]]
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
      
      SPseq = translateSP_fast(sp,peptide)
      #print("All SPs translated")
      
      CPseqClean = unique(CPseq)
      #print(paste("CP without doubles:",length(CPseqClean)))
      
      SPseqClean = unique(SPseq)
      l_SPseqClean = length(SPseqClean)
      #print(paste("SP without doubles:",length(SPseqClean)))
      
      #print("SPseqClean")
      #print(SPseqClean)
      x = removeCPfromSP_seq(SPseqClean,CPseqClean)
      SPseqClean = x[[1]]
      #print(paste("SP without CP:",length(SPseqClean)))
      
      prot_stats = data.frame(protein = attr(inputSequence, "name"),
                              all_PCP = length(CPseq),
                              all_PSP = length(SPseq),
                              unique_PCP = length(CPseqClean),
                              unique_PSP = l_SPseqClean,
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

# compute all PCP with length <=nmer
computeCPomplete <- function(L,nmer){
  
  maxL = nmer+1
  CP = numeric()
  
  for(i in 1:L){
    CP = rbind(CP, cbind(rep(i,length(c(i:min(L,(i+maxL-2))))),
                         c(i:min(L,(i+maxL-2)))))
  }
  return(CP)
  rm(maxL, CP)
}

CutAndPaste_seq_PCP <- function(inputSequence,nmer){
  # Make PCP sequences from a given input
  
  results = list()
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
      
      prot_stats = data.frame(protein = attr(inputSequence, "name"),
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
}



getPTMcombinations_fast <- function(s = peptide, m = MW, NmaxMod=max_variable_PTM, mods_input=mods){
  aa = strsplit(s,split="")[[1]]
  
  kn = which(mods_input$Site=="N-term" & mods_input$Position=="Any N-term")
  kc = which(mods_input$Site=="C-term" & mods_input$Position=="Any C-term")
  
  ka = list()
  for(i in 1:length(aa)){
    ka[[i]] = which(mods_input$Site%in%aa[i] & mods_input$Position=="Anywhere")
  }
  ka = unlist(ka)
  
  modIndex = c(kn,kc,ka)
  modId = mods_input$Id[modIndex]
  modSite = mods_input$Site[modIndex]
  modPos = mods_input$Position[modIndex]
  modDelta = as.numeric(as.vector(mods_input$MonoMass[modIndex]))
  
  # get all combinations of NmaxMod modifications
  combi = list()
  combi[[1]] = c(1:length(modIndex))
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod,length(modIndex))){
      combi[[i]] <- arrangements::combinations(c(1:length(modIndex)),i)
    }
  }
  
  deltaMass = list()
  IDs = list()
  deltaMass[[1]] = modDelta[combi[[1]]]
  IDs[[1]] = modId[combi[[1]]]
  
  if(NmaxMod>1 & length(modIndex)>1){
    for(i in 2:min(NmaxMod,length(modIndex))){
      deltaMass[[i]] <- rowSums(matrix(modDelta[combi[[i]]],dim(combi[[i]])[1],i))
      
      IDs[[i]] = matrix(modId[combi[[i]]],dim(combi[[i]])[1],i)
      IDs[[i]] <- do.call(stringi::stri_join, c(input=as.data.frame(IDs[[i]], stringsAsFactors = F), sep=";"))
    }
  }
  
  # generate final mod sequences with delta Masses
  PTMcombinations = data.table(peptide = rep(s,sum(unlist(lapply(deltaMass,length)))),
                               ids = unlist(IDs),
                               MW = m+unlist(deltaMass))
  return(PTMcombinations)
}

getPTMcombinations_fast_vec <- Vectorize(getPTMcombinations_fast, vectorize.args = c("s", "m"), SIMPLIFY = F)
