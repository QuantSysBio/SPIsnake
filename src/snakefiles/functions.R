### ---------------------------- Proteome pre-processing ----------------------------
Split_max_length <- function(x, max_length=500, overlap_length=MiSl*2){
  lapply(x, function(x){
    
    if (nchar(x) <= max_length) {
      protein_chunks <- x
    } else {
      # Extract coordinates
      totalsize = 200*nchar(x) 
      
      at <- breakInChunks(totalsize = totalsize, chunksize = max_length)
      at <- successiveIRanges(width = width(at), gapwidth = -overlap_length)
      
      # Subset last relevant IRanges
      keep <- which((at@start + at@width) <= nchar(x))
      keep <- c(keep, max(keep) + 1)
      at <- at[keep]
      
      ### Correct the last coordinate
      at@width[length(at@width)] <- as.integer(nchar(x) - at@start[length(at@start)] + 1)
      
      # Extract sequence
      protein_chunks <- extractAt(x, at)
      
      # Update names
      names(protein_chunks) <- paste0("|chunk:", at@start, "-", at@start + at@width - 1)
    }
    return(protein_chunks)
  })
}

Save_prot_chunk_biostrings <- function(dat, 
                                       nF=nF, 
                                       Pi=Pi,
                                       orderedProteomeEntries=orderedProteomeEntries, 
                                       directory=directory,
                                       proteome_name=proteome_name,
                                       protein_counter_start=protein_counter_start,
                                       protein_counter_end=protein_counter_end,
                                       maxE=maxE){
  for(j in 1:nF){
    start=Pi[j]+1
    end=Pi[j+1]
    
    writeXStringSet(x = dat[orderedProteomeEntries[start:end]], 
                    filepath = paste0(directory, "/", maxE, "_", proteome_name, "_", protein_counter_start, "_", protein_counter_end, "_", j,".fasta"))
  }
}

### ---------------------------- PCP/PSP generation ----------------------------
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
  SP = numeric()
  N = dim(cp)[1]
  NN = 5 * (10**6)
  
  SP = matrix(NA,NN,4)
  #print(SP)
  a = 1
  # repeat as many times as you have cp
  for(i in 1:N){
    temp1 = rep(cp[i,1],N)
    temp2 = rep(cp[i,2],N)
    temp3 = cp[,1]
    temp4 = cp[,2]
    
    L = temp4-temp3+temp2-temp1+2
    
    # L!=maxL
    ind = which(((temp3-temp2)==1)|(L>maxL)|(L<minL)|((temp3-temp2)>(MiSl+1))|((temp1-temp4)>(MiSl+1))|((temp3<=temp2)&(temp4>=temp1)))
    
    #print("ind")    
    #print(ind)
    
    if(length(ind)>0){
      temp1 = temp1[-ind]
      temp2 = temp2[-ind]
      temp3 = temp3[-ind]
      temp4 = temp4[-ind]
    }
    
    #
    if(length(temp1)!=1){
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

CutAndPaste_seq_return_sp <- function(inputSequence,nmer,MiSl){
  ### Makes PCP, cis-PSP, revcis-PCP from a given input
  # names(inputSequence)="test"
  
  results = list()
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

Generate_PSP_2 <- function(protein_inputs){
  stri_join(extractAt(protein_inputs[[1]], IRanges(start = protein_inputs[[2]][,1], 
                                                   end =   protein_inputs[[2]][,2])), 
            extractAt(protein_inputs[[1]], IRanges(start = protein_inputs[[2]][,3], 
                                                   end =   protein_inputs[[2]][,4])))
}

### ---------------------------- MW and PTMs ----------------------------

read_MW_file <- function(file, num_threads){
  vroom(file = file, 
        delim = " ", num_threads = num_threads, 
        col_names = c("Precursor_mass","RT"),
        show_col_types = FALSE) %>%
    lazy_dt() %>%
    mutate(MW_Min = Precursor_mass - Precursor_mass * tolerance * 10 ** (-6)) %>%
    mutate(MW_Max = Precursor_mass + Precursor_mass * tolerance * 10 ** (-6)) %>%
    # Min/sec for RT
    # mutate(RT = RT / 60) %>%
    mutate(RT_Min = RT - RT_tolerance) %>%
    mutate(RT_Max = RT + RT_tolerance) %>%
    select(-Precursor_mass)  %>%
    unique() %>%
    as.data.table()
}

computeMZ_biostrings <- function(seq){
  seq = AAStringSet(seq)
  
  if(length(seq) < 1){
    MW = NA
  }
  if(length(seq) > 0){
    aa1 = c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    aa2 = c(71.037114, 156.101111, 114.042927, 115.026943, 103.009185, 129.042593, 128.058578, 57.021464, 137.058912, 113.084064, 113.084064, 128.094963, 131.040485, 147.068414, 97.052764, 87.032028, 101.047679, 186.079313, 163.06332, 99.068414)
    aa3 = letterFrequency(seq, letters = aa1) %*% diag(aa2)
    MW=rowSums(aa3)+18.01528
    ## add NA for zero values
    # MW=round(MW, digits = 5)
    MW[MW==18.01528] = NA
  }
  return(MW)
}

### ---------------------------- RT filtering ----------------------------
regression_stats <- function(obs, pred){
  # lm
  data = data.frame(pred = pred, 
                    obs = obs)
  pred.lm = lm(pred ~ obs, data = data)
  
  # metrics
  # correlation coefficients
  pc = cor(obs, pred, method = "pearson")
  sm = cor(obs, pred, method = "spearman")
  
  # mean squared error
  mse =  round(mean((obs - pred)^2), 4)
  # root mean squared error
  rmse = round(sqrt(mse), 4)
  # mean absolute deviation
  mae =  round(mean(abs((obs - pred))), 4)
  
  # sumarize
  all.metrics = c(summary(pred.lm)$r.squared, pc, mse, rmse, mae)
  names(all.metrics) = c("Rsquared", "PCC", "MSE", "RMSE", "MAE")
  
  return(all.metrics)
}


### ---------------------------- PTM generation ----------------------------

getPTMcombinations_fast <- function(s, m, NmaxMod,mods_input=mods){
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
