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




