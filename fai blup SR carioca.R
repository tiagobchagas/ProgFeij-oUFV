#FAI BLUP 

getwd()
setwd("/Users/tiagobchagas/Desktop/Qualificação/seleção recorrente selecionar")
#file_paths=list.files(pattern = "*.csv", full.names = T)


library(dplyr)
library(fs)
library(tidyverse)

#rm(target)
#ff
#getwd()
file_paths = fs::dir_ls("/Users/tiagobchagas/Desktop/Qualificação/seleção recorrente selecionar")

list_of_file <- file_paths %>%
  map(function (path) {
    read_csv(path)
  })

names(list_of_file) = c("Prod2019", "Ag2020", "Prod2020", "Ag2021", "Arq2021", "Prod2021", "Ag2022", "Arq2022", "Prod2022") 

#names(list_of_file) = "/Users/tiagobchagas/Desktop/Qualificação/seleção recorrente selecionar/Effect_Mix_lme4-"="")

# convert each list into a tibble
list_of_file%>% map(as_tibble) %>% 
  
# set each list item to an environment variable
  list2env(envir = .GlobalEnv)

# merge two data frames by ID
#total <- merge(data frameA,data frameB,by="ID")


#rename cols

colnames(Prod2019)[2] = "prod2019"
colnames(Prod2020)[2] = "prod2020"
colnames(Ag2020)[2] = "ag2020"
colnames(Prod2021)[2] = "prod2021"
colnames(Ag2021)[2] = "ag2021"
colnames(Arq2021)[2] = "arq2021"
colnames(Prod2022)[2] = "prod2022"
colnames(Ag2022)[2] = "ag2022"
colnames(Arq2022)[2] = "arq2022"


#Fatores
Prod2022 <- Prod2022 %>% 
  mutate_at(vars(Line), as.factor)

str(Prod2022)
#filtrar as famílias anteriormente selecionadas nas safras anteriores
target <- c(Prod2022[,1])
#target <- paste(target, collapse = ", ")
#target == Prod2022["Line"]
target=unlist(target)
#rm(MeanProd2019)
MeanProd2019= select(filter(Prod2019, Line %in% c(target)), c(Line,prod2019))
MeanProd2020= select(filter(Prod2020, Line %in% c(target)), c(Line,prod2020))
MeanAg2020= select(filter(Ag2020, Line %in% c(target)), c(Line,ag2020))
MeanProd2021= select(filter(Prod2021, Line %in% c(target)), c(Line,prod2021))
MeanAg2021= select(filter(Ag2021, Line %in% c(target)), c(Line,ag2021))
MeanArq2021= select(filter(Arq2021, Line %in% c(target)), c(Line,arq2021))
MeanProd2022= select(filter(Prod2022, Line %in% c(target)), c(Line,prod2022))
MeanAg2022= select(filter(Ag2022, Line %in% c(target)), c(Line,ag2022))
MeanArq2022= select(filter(Arq2022, Line %in% c(target)), c(Line,arq2022))



cariocatraits=cbind(MeanProd2019,MeanProd2020[,2],MeanProd2021[,2],MeanProd2022[,2],MeanAg2020[,2],MeanAg2021[,2],MeanAg2022[,2], MeanArq2021[,2], MeanArq2022[,2])
cariocatraits

####################################################################################################################
################################################## FAI-BLUP index ##################################################
####################################################################################################################
##### Arguments 
## Data: Data file lies: genotypes in the first column and BLUP means (genetic values) in other columns.
## Show: Shows some components of factors analysis.
## D.ideotype: Desirable ideotype (max, min, mean or numeric for each trait).
## U.ideotype: Undesirable ideotype (max, min, mean or numeric for each trait).
## SN: Number of genotypes to be selected.
## eigen.value.min: Criterion for choosing the number of factors 
####################################################################################################################
## Created: 06-October-2016.
## By: Tiago de Souza MarÁal e Jo„o Romero do Amaral Santos de Carvalho Rocha
####################################################################################################################
require(MASS) ## The MASS package is required to execute the routine in R program. 
####################################################################################################################
fai.blup <- function(data, show = TRUE, ideotype.D, ideotype.U, SN = NULL, eigen.value.min = 1){
  means <- data[,2:ncol(data)]
  rownames(means) <- data[,1]
  normalize.means <- scale(means, center = FALSE, scale = apply(means, 2, sd))
  cor.means <- cor(normalize.means)
  eigen.decomposition <- eigen(cor.means)
  eigen.values <- eigen.decomposition$values 
  eigen.vectors <- eigen.decomposition$vectors
  colnames(eigen.vectors) <- paste("PC",1:ncol(cor.means),sep="")
  rownames(eigen.vectors) <- colnames(means)
  if(length(eigen.values[eigen.values >= eigen.value.min]) == 1){
    eigen.values.factors <- as.vector(c(as.matrix(sqrt(eigen.values[eigen.values >= eigen.value.min]))))
    initial.loadings <- cbind(eigen.vectors[, eigen.values >= eigen.value.min]*eigen.values.factors)
    finish.loadings <- initial.loadings 
  }
  if(length(eigen.values[eigen.values >= eigen.value.min]) > 1){
    eigen.values.factors <- t(replicate(ncol(cor.means), c(as.matrix(sqrt(eigen.values[eigen.values >= eigen.value.min])))))
    initial.loadings <- eigen.vectors[, eigen.values >= eigen.value.min]*eigen.values.factors
    finish.loadings <- varimax(initial.loadings)[[1]][]
  }  
  colnames(finish.loadings) <- paste("FA",1:ncol(initial.loadings),sep="")
  rownames(finish.loadings) <- colnames(means)
  comunalits <- rowSums(finish.loadings^2)
  cumulative.var <- cumsum(eigen.values/sum(eigen.values))*100
  pca <- cbind(eigen.values,cumulative.var)
  rownames(pca) <- paste("PC",1:ncol(means),sep="") 
  fa <- cbind(finish.loadings,comunalits)
  canonical.loadings <- ginv(finish.loadings%*%t(finish.loadings))%*%finish.loadings
  rownames(canonical.loadings) <- colnames(means)    
  scores <- t(t(canonical.loadings)%*%t(normalize.means))
  colnames(scores) <- paste("SC",1:ncol(scores),sep="")
  rownames(scores) <- data[,1]
  IN <- 2^ncol(finish.loadings)
  pos.var.factor <- which(abs(finish.loadings) == apply(abs(finish.loadings),1,max) , arr.ind = T)
  var.factor <- lapply(1:ncol(finish.loadings),function(i){rownames(pos.var.factor)[pos.var.factor[,2] == i]})
  names(var.factor) <- paste("FA",1:ncol(finish.loadings),sep="")
  names.pos.var.factor <- rownames(pos.var.factor)
  names(ideotype.D) <- colnames(means)
  names(ideotype.U) <- colnames(means)
  ideotype.D.test <- as.numeric(gsub("[^0-9]","",x = ideotype.D))
  ideotype.U.test <- as.numeric(gsub("[^0-9]","",x = ideotype.U))
  names(ideotype.D.test) <- colnames(means)
  names(ideotype.U.test) <- colnames(means)
  ideotype.D.test <- ideotype.D.test[names.pos.var.factor]  
  ideotype.U.test <- ideotype.U.test[names.pos.var.factor]
  canonical.loadings.factor <- canonical.loadings[names.pos.var.factor,] 
  ideotype.factor.D <- ideotype.D[names.pos.var.factor]
  ideotype.factor.U <- ideotype.U[names.pos.var.factor]             
  id.D <- rev(paste("D",1:ncol(finish.loadings),sep=""))
  id.U <- rev(paste("U",1:ncol(finish.loadings),sep=""))
  D.U <- rbind(id.D,id.U)
  groups.factor <- lapply(1:ncol(finish.loadings),function(i){D.U[,i]})
  construction.ideotypes <- as.matrix(rev(expand.grid(groups.factor)))
  colnames(construction.ideotypes) <- paste("Factor",1:ncol(construction.ideotypes),sep="")
  D <- numeric(0)
  U <- numeric(0)
  normalize.means.factor <- normalize.means[,names.pos.var.factor]    
  for(i in 1:ncol(normalize.means)){
    if(is.na(ideotype.D.test[i])){
      if(ideotype.factor.D[i] == "max"){
        D <- c(D, max(normalize.means.factor[,i]))
      }
      if(ideotype.factor.D[i] == "min"){
        D <- c(D, min(normalize.means.factor[,i]))
      }
      if(ideotype.factor.D[i] == "mean"){
        D <- c(D, mean(normalize.means.factor[,i]))
      }                                     
    }
    if(!is.na(ideotype.D.test[i])){
      D <- c(D, as.numeric(ideotype.factor.D[i]))                  
    }
    if(is.na(ideotype.U.test[i])){
      if(ideotype.factor.U[i] == "max"){
        U <- c(U, max(normalize.means.factor[,i]))
      }
      if(ideotype.factor.U[i] == "min"){
        U <- c(U, min(normalize.means.factor[,i]))
      }
      if(ideotype.factor.U[i] == "mean"){
        U <- c(U, mean(normalize.means.factor[,i]))
      }                  
    }
    if(!is.na(ideotype.U.test[i])){
      U <- c(U, as.numeric(ideotype.factor.U[i]))                  
    }               
  }
  names(D) <- names(ideotype.factor.D)
  names(U) <- names(ideotype.factor.U)
  Di <- lapply(1:ncol(finish.loadings),function(i){D[pos.var.factor[,2] == i]})
  Ui <- lapply(1:ncol(finish.loadings),function(i){U[pos.var.factor[,2] == i]})
  names(Di) <- paste("D",1:ncol(finish.loadings),sep="")
  names(Ui) <- paste("U",1:ncol(finish.loadings),sep="")
  comb.U.D <- c(Di,Ui)
  ideotypes.matrix <- matrix(0,IN,ncol(means))
  for(i in 1:IN){
    ideotypes.matrix[i,] <- unlist(comb.U.D[construction.ideotypes[i,]])
  }
  rownames(ideotypes.matrix) <- paste("ID",1:IN, sep = "")
  colnames(ideotypes.matrix) <- colnames(normalize.means.factor)
  ideotypes.scores <- ideotypes.matrix%*%canonical.loadings.factor
  sd.scores <- scale(rbind(scores,ideotypes.scores), center = FALSE, scale = apply(rbind(scores,ideotypes.scores), 2, sd))
  DE <- dist(sd.scores)
  DEM <- as.matrix(sqrt( (1/ncol(scores))*((DE)^2) ))
  GID <- DEM[1:nrow(scores), (nrow(scores) + 1):nrow(sd.scores)]
  spatial.prob <- (1/GID)/(replicate(IN, c(as.numeric(apply((1/GID),1,sum)))))
  ideotype.rank <- lapply(1:IN, function(i){sort(spatial.prob[,i],decreasing = TRUE)})
  names(ideotype.rank) <- paste("ID",1:IN,sep="")
  means.factor <- means[,names.pos.var.factor] 
  if(!is.null(SN)){
    genetic.gain <- lapply(1:IN, function(i){cbind(pos.var.factor[,2], ((colMeans(means.factor[names(ideotype.rank[[i]])[1:SN],]) - colMeans(means.factor))/colMeans(means.factor))*100)})
    for(i in 1:IN){
      colnames(genetic.gain[[i]]) <- c("Factor","Genetic Gain (%)")
    }
    names(genetic.gain) <- paste("ID",1:IN,sep="")
  }
  if(is.null(SN)){
    genetic.gain <- NULL
  }    
  if(show){
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nPrincipal Component Analysis\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(pca)              
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nFactor Analysis\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(fa)
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nComunalit Mean:",mean(comunalits),"\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nIdeotype Numbers:",IN,"\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nIdeotype Matrix\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(data.frame(construction.ideotypes))
    cat("\n\n")
    print(var.factor)
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nFAI-BLUP Index\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(list(ID1 = ideotype.rank$ID1))
    cat("\n-----------------------------------------------------------------------------------\n")
    if(!is.null(SN)){
      cat("\n Genetic Gain\n")
      cat("\n-----------------------------------------------------------------------------------\n")
      print(list(ID1 = genetic.gain$ID1))
      cat("\n\n")
      print(list("Selected Genotypes" = names(ideotype.rank[[1]])[1:SN]))
      cat("\n-----------------------------------------------------------------------------------\n")
    }              
  }
  if(!show){
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nIdeotype Matrix\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(data.frame(construction.ideotypes))
    cat("\n\n")
    print(var.factor)                           
    cat("\n-----------------------------------------------------------------------------------\n")
    cat("\nFAI-BLUP Index\n")
    cat("\n-----------------------------------------------------------------------------------\n")
    print(list(ID1 = ideotype.rank$ID1))
    cat("\n-----------------------------------------------------------------------------------\n")
    if(!is.null(SN)){
      cat("\nGenetic Gain\n")
      cat("\n-----------------------------------------------------------------------------------\n")
      print(list(ID1 = genetic.gain$ID1))
      cat("\n\n")
      print(list("Selected Genotypes" = names(ideotype.rank[[1]])[1:SN]))
      cat("\n-----------------------------------------------------------------------------------\n")
    }
  }        
  output <- list(IN = IN, comunalits = comunalits, finish.loadings = finish.loadings, canonical.loadings = canonical.loadings, scores = scores,construction.ideotypes = data.frame(construction.ideotypes), final.fai.blup.rank = ideotype.rank, genetic.gain = genetic.gain) 
}

####################################################################################################################

#setwd("C:\\Users\\t_ago\\Desktop\\DS UFV\\an·lise-conjunta\\")

#data=read.table("C:\\Users\\t_ago\\Desktop\\DS UFV\\an·lise-conjunta\\data.txt",h=T) # Read data

data=cariocatraits
D.ideotype <- c("max","max","max","max","min","min","min","min","min") ## Derirable ideotype
U.ideotype <- c("min","min","min","min","max","max","max","max","max") ## Underirable ideotype

names(data)
(data[,-1])

ID=fai.blup(data = data, show = TRUE, ideotype.D = D.ideotype, ideotype.U = U.ideotype, SN = 121, eigen.value.min = 0.7) # FAI-BLUP index output

ID$final.fai.blup.rank$ID1
ID$final.fai.blup.rank$ID2
ID$final.fai.blup.rank$ID3
ID$scores

genotipos = data.frame(ID$final.fai.blup.rank$ID1)
genotipos <- data.frame(names = row.names(genotipos), genotipos)
#rm(genotipos)
str(ID)
str(genotipos)
genotipos$RANK<-rank(genotipos$ID.final.fai.blup.rank.ID1)

colnames(genotipos)[1] = "Line"
dados = merge(data,genotipos, by= "Line" )

#rm(dados)
#rm(dados)
#write.table(genotipos, "/Users/tiagobchagas/Desktop/Qualificação/seleção recorrente selecionar"      ,sep=";", col.names=TRUE, row.names=TRUE, quote=FALSE, na="NA")
library(openxlsx)
write.xlsx(dados, "CariocaGenotipos.xlsx", sheetName = "MEANS", 
           colNames = TRUE, rowNnames = TRUE, append = FALSE)
getwd()
#setwd("/Users/tiagobchagas/Desktop/HTP CREA 2022/")
write.xlsx(cariocameans, "CariocaMeans.xlsx", sheetName = "MEANS", 
           colNames = TRUE, rowNnames = F, append = FALSE)