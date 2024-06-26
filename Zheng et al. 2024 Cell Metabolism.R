
# library(devtools)
# devtools::install_github("XiaoyangSu/AccuCor")

library(accucor)

require(gsubfn)
require(nnls)
require(tidyverse)
require(gdata)
require(openxlsx) 

setwd("G:/My Drive/Data/LCMS/NAC/")

#type in the name of the data file WITHOUT the .xlsx extension
inputFileName<-"Zheng 2023-09-21 2023-09-20 15N Gln GAR re-run"

#code below should read your data file and 'metabolist' file;

#also generated is the name of the results file which is your original file name plus 'Results' in the end

rawData<-read.xlsx(paste0(getwd(),"/",inputFileName,".xlsx"), sheet = 1,
                   skipEmptyRows = T, skipEmptyCols = T, 
                   check.names = F, sep.names = " ") # ! sep.names = " "

metaboliteList <- read.csv(paste0(getwd(),"/nam.csv"),header=T,check.names=F)

outputFile <- paste0(getwd(),"/",inputFileName," Results.xlsx")

#take a look at the data

head(rawData)
dim(rawData)
glimpse(rawData)


# keep columns that are not all NAs
dim(rawData)
rawData<-
  rawData[,!apply(is.na(rawData),2,all)] 
dim(rawData)

# # clean data UNIQUE TO THIS SET

dim(rawData)
rawData<-
  rawData %>%
  # slice_head(n = 5) %>% ## ONLY FIRST FIVE ROWS IN THIS CASE!!
  filter(X1!="Negative TIC") %>%
  filter(X1!="-TIC") %>% 
  filter(X1!="+TIC") %>% 
  filter(X1!="TIC Area") %>% 
  filter(!grepl('ISTD',X1))  %>% 
  filter(!grepl('Citrate \\[M\\+0\\]\\/Isocitrate \\[M\\+0\\]',X1))  %>%
  filter(!grepl('Citrate 13C\\[M',X1))



dim(rawData)
## clean data END

# create a vector that contains original column names for later SELECT use
vecOriginalColumns<-colnames(rawData)[-1] 


# convert to such data format: rows are observations(i.e. metabolites or isotopomers), 
# columns are variables (such as compoundName, fileName, peakArea)

rawData<-
  rawData%>%
  rename(compoundName=1) %>%
  gather(fileName,peakArea,-compoundName)

head(rawData)

# check if there are any duplicate observations
# if the dim(rawData) returns the same numbers before and after the distinct() line, then there was no duplication
# if the first returned number is smaller after distinct(), 
# then there were duplications and only the first observation was kept, all others were deleted

unique(rawData$compoundName)
dim(rawData)
rawData<-
  rawData %>% 
  distinct(compoundName,fileName,.keep_all = T)
dim(rawData)


# show all the file names, please note that no two file names should be the same
# convert the fileName column into character if it is not a character vector already
unique(rawData$fileName)
class(rawData$fileName)
rawData$fileName<-as.character(rawData$fileName)

# show all the compound names, please note that if it is tracing data, naming format should be:
# Metabolite M+N
# an example: Lactate M+3
# convert the compoundName column into character if it is not a character vector already
unique(rawData$compoundName)
class(rawData$compoundName)
rawData$compoundName<-as.character(rawData$compoundName)




## interface with accucor basic data format
head(rawData)

rawData$Compound<-trimws(trimws(rawData$compoundName,"right","\\S"),"right","\\s")

rawData<-
rawData %>% 
  mutate(
    Isotope = case_when(
      grepl("\\[M\\+0\\]", compoundName) ~ "C12 PARENT",
      grepl("13C\\[M\\+", compoundName) ~ "C13-label-",
      grepl("15N\\[M\\+", compoundName) ~ "N15-label-",
      grepl("2H\\[M\\+", compoundName) ~ "D-label-",
            TRUE ~ "Other"
    )
  )

rawData$Label<-
substring(rawData$compoundName,regexpr("\\[M\\+",rawData$compoundName) + 3, 
          regexpr("\\]",rawData$compoundName) - 1)
rawData<-
rawData %>% 
  mutate(
    LabelTemp = case_when(
      Label == "0" ~ "",
      TRUE ~ Label
    )
  ) %>% 
  unite("IsotopeLabel",c(Isotope,LabelTemp), sep = "", remove = F) %>% 
  select(-LabelTemp)

# check to see whether there are negative or NA values in peak area
summary(rawData$peakArea)
mean(!is.na(rawData$peakArea))

# convert any "N/F" in the peak area into 0
rawData$peakArea<-as.character(rawData$peakArea)
rawData$peakArea[rawData$peakArea=="N/F"]<-"0"
rawData$peakArea[rawData$peakArea=="NaN"]<-"0"
rawData$peakArea<-as.numeric(rawData$peakArea)
rawData$peakArea[rawData$peakArea<0]<-0


# check to see whether every compound is documented in the metaboliteList file (named nam.csv)
# fraction of the compounds that are documented in the metaboliteList file
mean(rawData$Compound %in% metaboliteList$Compound )

# if not, below lines will tell you which metabolite(s) is not on file, add them to file 
# all the compounds in your data:
(metabolitesInData<-unique(rawData$Compound))
# true or false they are on the list:
(logicMetabolitesInDataDocumented<-unique(rawData$Compound)%in% metaboliteList$Compound)

# the one(s) that is not on the list:
metabolitesInData[!logicMetabolitesInDataDocumented]




# if necessary, make adjustment to nam.csv, start the script over 





# add the Formula column to the rawData dataframe

rawData<-
rawData %>% 
  left_join(metaboliteList, by = c("Compound"="Compound")) 

dfAllSampleInput<-
rawData %>% 
  select(Compound,Formula,IsotopeLabel,fileName,peakArea) %>% 
  spread(key = fileName, value = peakArea) %>% 
  select(Compound, Formula, IsotopeLabel, all_of(vecOriginalColumns)) 

dfCarbonSampleInput<-
rawData %>% 
  filter(Isotope %in% c("C12 PARENT","C13-label-")) %>% 
  select(Compound,Formula,IsotopeLabel,fileName,peakArea) %>% 
  spread(key = fileName, value = peakArea) %>% 
  select(Compound, Formula, IsotopeLabel, all_of(vecOriginalColumns)) 

dfNitrogenSampleInput<-
  rawData %>% 
  filter(Isotope %in% c("C12 PARENT","N15-label-")) %>% 
  select(Compound,Formula,IsotopeLabel,fileName,peakArea) %>% 
  spread(key = fileName, value = peakArea) %>% 
  select(Compound, Formula, IsotopeLabel, all_of(vecOriginalColumns))

dfDeuteriumSampleInput<-
  rawData %>% 
  filter(Isotope %in% c("C12 PARENT","D-label-")) %>% 
  select(Compound,Formula,IsotopeLabel,fileName,peakArea) %>% 
  spread(key = fileName, value = peakArea) %>% 
  select(Compound, Formula, IsotopeLabel, all_of(vecOriginalColumns))

## unique to this dataset!
# dfCarbonSampleInput<-
#   dfCarbonSampleInput %>% 
#   select(!matches("N_"))
# 
# dfNitrogenSampleInput<-
#   dfNitrogenSampleInput %>% 
#   select(!matches("C_")) 


## unique to this dataset end!


write.csv(x = dfCarbonSampleInput, file = "CarbonSampleInput.csv", row.names = F)
write.csv(x = dfNitrogenSampleInput, file = "NitrogenSampleInput.csv", row.names = F)
write.csv(x = dfDeuteriumSampleInput, file = "DeuteriumSampleInput.csv", row.names = F)

# Input file (example file included)
# Or use your own: carbon_input_file <- "/path/to/my/datafile.csv"

# carbon_input_file <- system.file("extdata", "C_Sample_Input_Simple.csv", package = "accucor")


# Output is written to [input_file]_corrected.xlsx by default
# Be sure to specify the appropriate resolution.
# For Exactive, the resolution is 100000, defined at 200 Mw

carbon_corrected <- natural_abundance_correction(
  path = paste0(getwd(),"/","CarbonSampleInput.csv"),
  resolution = 240000,
  purity = 0.99)


# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "PoolBeforeDF", "PoolAfterDF"

carbon_corrected

# get not corrected (before correction) fractional enrichment 
dfNotCorrected<-
carbon_corrected$Original

dfNotCorrected <-   
  dfNotCorrected %>% 
  select(-Formula)  

dfNotCorrected$IsotopeLabel<-gsub("C12 PARENT","0",dfNotCorrected$IsotopeLabel) 
dfNotCorrected$IsotopeLabel<-gsub("C13-label-","",dfNotCorrected$IsotopeLabel)
dfNotCorrected$IsotopeLabel<-as.numeric(dfNotCorrected$IsotopeLabel)
colnames(dfNotCorrected)[2]<-"Label"

listNotCorrected <- split(dfNotCorrected,dfNotCorrected$Compound) 

scale_in_dfNotCorrected <- 
  function(x){ cbind(x[,1:2],scale(x[,-(1:2)], scale=colSums(x[,-(1:2)]),center=FALSE))}
listNotCorrectedPercentage<-lapply(listNotCorrected, scale_in_dfNotCorrected)

dfNotCorrectedPercentage<-bind_rows(listNotCorrectedPercentage)

listSheets <- list("correctedPeakArea" = carbon_corrected$Corrected,
                   "correctedFractionalEnrichment" = carbon_corrected$Normalized,
                   "correctedPoolSize" = carbon_corrected$PoolAfterDF,
                   "rawPeakArea" = carbon_corrected$Original,
                   "rawFractionalEnrichment" = dfNotCorrectedPercentage)
write.xlsx(listSheets, file = outputFile)






## normalize to tic,

## unique to this dataset ##
ticDF<-
   read.xlsx(paste0(getwd(),"/",inputFileName,".xlsx"), sheet = 2, skipEmptyRows = T,
             skipEmptyCols = T, check.names = F)   %>% ## CHANGE TAB ACCORDINGLY
     filter(X1 %in% c("Negative TIC","TIC Area")) %>%
     gather(fileName,TIC, -X1) %>%
     select(-X1)

colnames(ticDF)<-c("fileName", "TIC")


## STOP HERE FOR CARBON CORRECTION


# Output is written to [input_file]_corrected.xlsx by default
# Be sure to specify the appropriate resolution.
# For Exactive, the resolution is 100000, defined at 200 Mw

nitrogen_corrected <- natural_abundance_correction(
  path = paste0(getwd(),"/","NitrogenSampleInput.csv"),
  resolution = 400000,
  purity = 0.98)

nitrogen_corrected

# get not corrected (before correction) fractional enrichment 
dfNotCorrected<-
  nitrogen_corrected$Original

dfNotCorrected <-   
  dfNotCorrected %>% 
  select(-Formula)  

dfNotCorrected$IsotopeLabel<-gsub("C12 PARENT","0",dfNotCorrected$IsotopeLabel) 
dfNotCorrected$IsotopeLabel<-gsub("N15-label-","",dfNotCorrected$IsotopeLabel)
dfNotCorrected$IsotopeLabel<-as.numeric(dfNotCorrected$IsotopeLabel)
colnames(dfNotCorrected)[2]<-"Label"

listNotCorrected <- split(dfNotCorrected,dfNotCorrected$Compound) 

scale_in_dfNotCorrected <- 
  function(x){ cbind(x[,1:2],scale(x[,-(1:2)], scale=colSums(x[,-(1:2)]),center=FALSE))}
listNotCorrectedPercentage<-lapply(listNotCorrected, scale_in_dfNotCorrected)

dfNotCorrectedPercentage<-bind_rows(listNotCorrectedPercentage)


##

outputFile <- paste0(getwd(),"/",inputFileName," Results 15N Res400k.xlsx")

##

listSheets <- list("correctedPeakArea" = nitrogen_corrected$Corrected,
                   "correctedFractionalEnrichment" = nitrogen_corrected$Normalized,
                   "correctedPoolSize" = nitrogen_corrected$PoolAfterDF,
                   "rawPeakArea" = nitrogen_corrected$Original,
                   "rawFractionalEnrichment" = dfNotCorrectedPercentage)
write.xlsx(listSheets, file = outputFile)




# Deuterium below:


deuterium_corrected <- natural_abundance_correction(
  path = paste0(getwd(),"/","DeuteriumSampleInput.csv"),
  resolution = 60000, # !!!
  purity = 0.99)


# The results are also returned as a named list of dataframes for further processing in R
# "Original", "Corrected", "Normalized", "PoolBeforeDF", "PoolAfterDF"

deuterium_corrected

# get not corrected (before correction) fractional enrichment 
dfNotCorrected<-
  deuterium_corrected$Original

dfNotCorrected <-   
  dfNotCorrected %>% 
  select(-Formula)  

dfNotCorrected$IsotopeLabel<-gsub("C12 PARENT","0",dfNotCorrected$IsotopeLabel) 
dfNotCorrected$IsotopeLabel<-gsub("D-label-","",dfNotCorrected$IsotopeLabel)
dfNotCorrected$IsotopeLabel<-as.numeric(dfNotCorrected$IsotopeLabel)
colnames(dfNotCorrected)[2]<-"Label"

listNotCorrected <- split(dfNotCorrected,dfNotCorrected$Compound) 

scale_in_dfNotCorrected <- 
  function(x){ cbind(x[,1:2],scale(x[,-(1:2)], scale=colSums(x[,-(1:2)]),center=FALSE))}
listNotCorrectedPercentage<-lapply(listNotCorrected, scale_in_dfNotCorrected)

dfNotCorrectedPercentage<-bind_rows(listNotCorrectedPercentage)

listSheets <- list("correctedPeakArea" = deuterium_corrected$Corrected,
                   "correctedFractionalEnrichment" = deuterium_corrected$Normalized,
                   "correctedPoolSize" = deuterium_corrected$PoolAfterDF,
                   "rawPeakArea" = deuterium_corrected$Original,
                   "rawFractionalEnrichment" = dfNotCorrectedPercentage)

outputFile <- paste0(getwd(),"/",inputFileName," Results 2H 60k.xlsx")

write.xlsx(listSheets, file = outputFile)
