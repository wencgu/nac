#type in the paths of your choice
setwd("C:/Users/S174729/Google Drive/Data/LCMS/NAC")
OutputFile <- "C:/Users/S174729/Google Drive/Data/LCMS/NAC/Aslan 2019-06-18 20190618-AT_20190430_Tracing Results.xlsx"
RawData<-read.csv("C:/Users/S174729/Google Drive/Data/LCMS/NAC/Aslan 2019-06-18 20190618-AT_20190430_Tracing.csv")
#SampleInfo<-read.csv("C:/Users/S174729/Google Drive/Data/LCMS/NAC/Aslan 2019-05-08 sam.csv")
MetaboliteList <- read.csv("C:/Users/S174729/Google Drive/Data/LCMS/NAC/nam.csv",header=T,check.names=F)


require(gsubfn)
require(nnls)
require(tidyverse)
require(gdata)
require(xlsx) #might need to manually install x64 java

#convert table data to basic data 
head(RawData)
colnames(RawData)[1]<-"Compound.Name"
RawData<-
RawData %>% 
  gather(Sample.Raw.File.Name,Peak.Area,-Compound.Name)

#adjust problematic metabolites
#RawData$Compound.Name<-as.character(RawData$Compound.Name)
#RawData$Compound.Name[RawData$Compound.Name=="2-/3-Phosophoglycerate"]<-"2-/3-Phosophoglycerate M+0"

#check duplication
#get rid of the duplicates of glucose data
head(RawData)
unique(RawData$Compound.Name)
dim(RawData)
RawData<-RawData %>% distinct(Compound.Name,Sample.Raw.File.Name,.keep_all = T)
dim(RawData)


#check
unique(RawData$Sample.Raw.File.Name)
head(RawData)
class(RawData$Sample.Raw.File.Name)

#necessary
RawData$Sample.Raw.File.Name<-as.character(RawData$Sample.Raw.File.Name)

#necessary
#SampleInfo$Original_Name<-as.character(SampleInfo$Original_Name)

RawData<-RawData %>% dplyr::select(Compound.Name,Sample.Raw.File.Name,Peak.Area)%>%
    dplyr::mutate(Metabolite=gsub(" M.*$", "",RawData$Compound.Name)) %>% 
    dplyr::mutate(Isotopomer=gsub("^.* M.","",RawData$Compound.Name))

#check to see whether there is negative or NA values in peak area
summary(RawData$Peak.Area)
mean(!is.na(RawData$Peak.Area))
RawData<- RawData %>%
  select(Sample.Raw.File.Name,Peak.Area,Isotopomer,Metabolite)

RawData$Peak.Area<-as.character(RawData$Peak.Area)
RawData$Peak.Area[RawData$Peak.Area=="N/F"]<-"0"
RawData$Peak.Area<-as.numeric(RawData$Peak.Area)
RawData$Peak.Area[RawData$Peak.Area<0]<-0


CarbonNaturalAbundance <- c(0.9893, 0.0107)
HydrogenNaturalAbundance <- c(0.999885, 0.000115)
NitrogenNaturalAbundance <- c(0.99636, 0.00364)
OxygenNaturalAbundance <- c(0.99757, 0.00038, 0.00205)
SulfurNaturalAbundance <- c(0.9493, 0.00762, 0.0429)

C13Purity <- 0.99
ReportPoolSize <- TRUE

Resolution <- 60000 #60,000--TOM
ResDefAt <- 200

#check to see whether all metabolites from raw data are included in the metabolites list, vise versa
#if not, start over to change to proper names, or add metabolites to the Metabolist
mean(MetaboliteList$Compound %in% RawData$Metabolite)
mean(RawData$Metabolite %in% MetaboliteList$Compound )
unique(RawData$Metabolite)
unique(RawData$Metabolite)%in% MetaboliteList$Compound
count(RawData,Metabolite)
length(unique(RawData$Sample.Raw.File.Name))

#check
head(RawData)
unique(RawData$Isotopomer)

#further adjustment
#RawData$Metabolite[RawData$Metabolite=="6-Phosphogluconolactone neg"]<-"6-Phosphogluconolactone"

OutputUncorrectedMatrix<-           matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)
OutputUncorrectedPercentageMatrix<- matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)

OutputMatrix <-           matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)
OutputPercentageMatrix <- matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)
OutputPoolBefore <-       matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)
OutputPoolAfter <-        matrix(0, ncol=length(unique(RawData$Sample.Raw.File.Name)), nrow=0)

OutputCompound <- NULL
OutputLabel <- NULL
OutputPoolCompound <- NULL

if.not.null <- function(x) if(!is.null(x)) x else 0 #customized function for later extracting the numbers of atoms

IsotopeCorrection <- function(formula, datamatrix, label) {
  
  AtomNumber <- rep(0,6)
  names(AtomNumber) <- c("C","H","N","O","P","S")
  MassDifference <- abs(c((2.0141-1.00783),(15.00011-14.00307),(16.99913-15.99491),
                          (17.99916-15.99491),(32.97146-31.97207),(33.96787-31.97207))-
                          ((13.00335-12)*c(1,1,1,2,1,2)))
  names(MassDifference) <- c("H2","N15","O17","O18","S33","S34")
  CorrectionLimit <- rep(0,6)
  names(CorrectionLimit) <- c("H2","N15","O17","O18","S33","S34")
  
  
  AtomNumber["C"] <- if.not.null(unlist(strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["H"] <- if.not.null(unlist(strapply(formula, "(H)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["N"] <- if.not.null(unlist(strapply(formula, "(N)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["O"] <- if.not.null(unlist(strapply(formula, "(O)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["P"] <- if.not.null(unlist(strapply(formula, "(P)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  AtomNumber["S"] <- if.not.null(unlist(strapply(formula, "(S)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
  
  MolecularWeight <- sum(AtomNumber*c(12,1,14,16,31,32))
  
  CorrectionLimit <- floor(MolecularWeight^(3/2)*1.66/(Resolution*sqrt(ResDefAt))/MassDifference)
  
  ExpMatrix <-        matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["C"]+1)
  CorrectedMatrix <-  matrix(0, ncol=ncol(datamatrix), nrow=AtomNumber["C"]+1)
  
 for (i in 1:length(label)) {
    ExpMatrix[label[i]+1,] <- datamatrix[i,] 
  }
  
  
  PurityMatrix <- diag(AtomNumber["C"]+1)
  CarbonMatrix <- diag(AtomNumber["C"]+1)
  NitrogenMatrix <- diag(AtomNumber["C"]+1)
  HydrogenMatrix <- diag(AtomNumber["C"]+1)
  OxygenMatrix <- matrix(0,ncol=(AtomNumber["C"]+1),nrow=(AtomNumber["C"]+1))
  SulfurMatrix <- matrix(0,ncol=(AtomNumber["C"]+1),nrow=(AtomNumber["C"]+1))
  
  
  for(i in 1:(AtomNumber["C"]+1)){
    PurityMatrix[i,] <- sapply(0:(AtomNumber["C"]), function(x) dbinom(x-i+1, x , (1-C13Purity)))
  }
  
  for(i in 1:(AtomNumber["C"]+1)){
    CarbonMatrix[,i] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, AtomNumber["C"]-i+1 , CarbonNaturalAbundance[2]))
  }
  
  for(j in 0:min(AtomNumber["N"], CorrectionLimit["N15"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      NitrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["N"], NitrogenNaturalAbundance[2])
    }
  
  for(j in 0:min(AtomNumber["H"], CorrectionLimit["H2"], AtomNumber["C"]))
    for(i in 1:(AtomNumber["C"]-j+1)){
      HydrogenMatrix[(i+j),i] <- dbinom(j, AtomNumber["H"], HydrogenNaturalAbundance[2])
    }
  
  for(i in 0:min(AtomNumber["O"],CorrectionLimit["O17"])) {
    for(j in 0:min(AtomNumber["O"],CorrectionLimit["O18"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["O"]|k>AtomNumber["C"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["C"]-k+1)) {
          OxygenMatrix[(m+k),m] <- OxygenMatrix[(m+k),m] + dmultinom(c((AtomNumber["O"]-i-j),i,j), AtomNumber["O"], OxygenNaturalAbundance)      
        }
      }
    }
  }
  
  for(i in 0:min(AtomNumber["S"],CorrectionLimit["S33"])) {
    for(j in 0:min(AtomNumber["S"],CorrectionLimit["S34"])){
      k<-(i+j*2)
      if ((i+j)>AtomNumber["S"]|k>AtomNumber["C"]) {
        break
      }
      else {
        for (m in 1:(AtomNumber["C"]-k+1)) {
          SulfurMatrix[(m+k),m] <- SulfurMatrix[(m+k),m] + dmultinom(c((AtomNumber["S"]-i-j),i,j), AtomNumber["S"], SulfurNaturalAbundance)      
        }
      }
    }
  }
  
  for(i in 1:ncol(datamatrix)) {
    CorrectedMatrix[,i] <- coef(nnls(SulfurMatrix %*% OxygenMatrix %*% NitrogenMatrix %*% 
                                       HydrogenMatrix %*% CarbonMatrix %*% PurityMatrix, ExpMatrix[,i]))
  }
  
  return(CorrectedMatrix)
  
}


FillUncorrectedMatrix <- function(formula, datamatrix, label) {
  
  CNumber <- rep(0,1)
  
  names(CNumber) <- c("C")
  
  CNumber["C"] <- if.not.null(unlist(strapply(formula, "(C)(\\d*)", ~ as.numeric(if (nchar(..2)) ..2 else 1))))
   
  UncorrectedMatrix <-        matrix(0, ncol=ncol(datamatrix), nrow=CNumber["C"]+1)
   
  for (i in 1:length(label)) {
    UncorrectedMatrix[label[i]+1,] <- datamatrix[i,] 
  }
  
  return(UncorrectedMatrix)
  
}


for (i in unique(RawData$Metabolite)) {
  Formula=as.character(MetaboliteList$Formula[MetaboliteList$Compound==i])
   if(length(Formula)==0) {
    print(paste("The formula of",i,"is unknown",sep=" "))
    break
  }
  CurrentMetabolite <- RawData %>% 
    filter(Metabolite==i) %>%
    spread(Sample.Raw.File.Name,Peak.Area)
  DataMatrix <- data.matrix(CurrentMetabolite[,3:ncol(CurrentMetabolite)])
      
      
  DataMatrix[is.na(DataMatrix)] <- 0
  Corrected <- IsotopeCorrection(Formula, DataMatrix, as.numeric(CurrentMetabolite$Isotopomer))
  CorrectedPercentage <- scale(Corrected,scale=colSums(Corrected),center=FALSE)
  
  Uncorrected<-FillUncorrectedMatrix(Formula,DataMatrix, as.numeric(CurrentMetabolite$Isotopomer))
  OutputUncorrectedPercentage<- scale(Uncorrected,scale=colSums(Uncorrected),center=FALSE)
  
  OutputMatrix <- rbind(OutputMatrix, Corrected)
  OutputPercentageMatrix <- rbind(OutputPercentageMatrix, CorrectedPercentage)
  OutputPoolBefore <- rbind(OutputPoolBefore, colSums(DataMatrix))
  OutputPoolAfter <- rbind(OutputPoolAfter, colSums(Corrected))
  OutputCompound <- append(OutputCompound, rep(i,nrow(Corrected)))
  OutputLabel <- append(OutputLabel, c(0:(nrow(Corrected)-1)))
  OutputPoolCompound <- append(OutputPoolCompound, i)
  
  OutputUncorrectedMatrix <- rbind(OutputUncorrectedMatrix,Uncorrected)
  OutputUncorrectedPercentageMatrix <- rbind(OutputUncorrectedPercentageMatrix, OutputUncorrectedPercentage)
} 

OutputDF <- data.frame(OutputCompound, OutputLabel, OutputMatrix)
OutputPercentageDF <- data.frame(OutputCompound, OutputLabel, OutputPercentageMatrix)
OutputPoolBeforeDF <- data.frame(OutputPoolCompound, OutputPoolBefore)
OutputPoolAfterDF <- data.frame(OutputPoolCompound, OutputPoolAfter)

OutputUncorrectedDF<- data.frame(OutputCompound, OutputLabel, OutputUncorrectedMatrix)
OutputUncorrectedPercentageDF <- data.frame(OutputCompound, OutputLabel, OutputUncorrectedPercentageMatrix)

#adjust "names" sources for output purpose
names(OutputDF) <- c("Compound", "Label", unique(RawData$Sample.Raw.File.Name))
names(OutputPercentageDF) <- names(OutputDF)
names(OutputPoolBeforeDF) <- c("Compound", unique(RawData$Sample.Raw.File.Name))
names(OutputPoolAfterDF) <- c("Compound", unique(RawData$Sample.Raw.File.Name))

names(OutputUncorrectedDF) <- c("Compound", "Label", unique(RawData$Sample.Raw.File.Name))
names(OutputUncorrectedPercentageDF) <- names(OutputUncorrectedDF)


write.xlsx2(OutputDF, file=OutputFile, sheetName = paste("Corrected",sep="_"), row.names=FALSE, append=TRUE)
write.xlsx2(OutputPercentageDF, file=OutputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPoolAfterDF, file=OutputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)

write.xlsx2(OutputUncorrectedDF, file=OutputFile, sheetName = paste("Uncorrected",sep="_"), row.names=FALSE, append=TRUE)
write.xlsx2(OutputUncorrectedPercentageDF, file=OutputFile, sheetName = "Uncorrected Normalized", row.names=FALSE, append=TRUE)

#STOP
#poolsize normalize to tic area
SampleInfo
class(SampleInfo$TIC)
mean(SampleInfo$TIC)
SampleInfo<-SampleInfo %>% mutate(TIC.Relative=TIC/mean(SampleInfo$TIC))

OutputPoolAfterDF
#df.normalize to tic area
df.ntt<-OutputPoolAfterDF %>% gather(Number,Peak.Area,-Compound)
class(SampleInfo$Number)
SampleInfo$Number<-as.character(SampleInfo$Number)
class(df.ntt$Number)
df.ntt<-df.ntt %>% left_join(SampleInfo,by="Number") %>% mutate(Peak.Area.NTT=Peak.Area/TIC.Relative) %>%
  mutate(Peak.Area.NTT.SM=Peak.Area/TIC) %>%
  mutate(Peak.Area.NTT.LOG=log(Peak.Area.NTT+1))

write.csv(df.ntt,"temp.csv",row.names = F)


df.ntt.log<-df.ntt %>% select(Compound,Number,Peak.Area.NTT.LOG) %>% spread(Number,Peak.Area.NTT.LOG)
df.ntt.original<-df.ntt %>% select(Compound,Number,Peak.Area.NTT) %>% spread(Number,Peak.Area.NTT)
#df.ntt.sm<-df.ntt %>% select(Compound,Number,Peak.Area.NTT.SM) %>% spread(Number,Peak.Area.NTT.SM)

write.xlsx2(df.ntt.original, file=OutputFile, sheetName = paste("Pool2TIC",sep="_"), row.names=FALSE, append=TRUE)
write.xlsx2(df.ntt.log, file=OutputFile, sheetName = paste("Pool2TIC.LOG",sep="_"), row.names=FALSE, append=TRUE)
#write.xlsx2(df.ntt.sm, file=OutputFile, sheetName = paste("PoolNormalizedSmall",sep="_"), row.names=FALSE, append=TRUE)
df.ntt

#Corrected Normalized to TIC area
OutputUncorrectedDF
df.cnt<-
  OutputUncorrectedDF %>% 
    gather(Number,Peak.Area,-Compound,-Label)%>% 
    left_join(SampleInfo,by="Number") %>% 
    mutate(Peak.Area.CNT=Peak.Area/TIC.Relative) %>%
    mutate(Peak.Area.CNT.LOG=log(Peak.Area.CNT+1)) 

df.cnt.original<-
  df.cnt%>%
  select(Compound,Label,Number,Peak.Area.CNT) %>% 
  spread(Number,Peak.Area.CNT)

df.cnt.log<-
    df.cnt%>%
    select(Compound,Label,Number,Peak.Area.CNT.LOG) %>% 
    spread(Number,Peak.Area.CNT.LOG)

write.xlsx2(df.cnt.original, file=OutputFile, sheetName = paste("Corrected2TIC",sep="_"), row.names=FALSE, append=TRUE)
write.xlsx2(df.cnt.log, file=OutputFile, sheetName = paste("Corrected2TIC.LOG",sep="_"), row.names=FALSE, append=TRUE)
  




