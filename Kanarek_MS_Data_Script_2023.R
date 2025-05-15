# # R-Script to combine datasets into one Excel (extracting peak areas), normalize and filter metabos by Rsq and CV; updated Oct 18, 2022

options(scipen = 100) #displays decimals in numeric, NOT sig-figs
remove(list = ls()) #clean-up environment first


library(dplyr) #load dependencies
library(purrr)
#library(rmarkdown)

##### DESCRIPTION:
# This script functions on R versions 4.1.2 and 4.2.0, and on both MacOS (Catalina) and Windows 10 operating systems
#
# How to use:
# (In top toolbar:) Code > Run Region > Run All
# Keyboard shortcut on Macs: (Cmd + Option + R)
#
# Description:
# This script takes input of all the raw CSVs from TraceFinder, and outputs the following files (in order of normalization steps):
#
# CSV outputs:
# "MSdataframe.csv" --> raw data of combined CSVs from TraceFinder
# "standards.csv" --> CSV of just the isoAAs, QRESS and aminopterin standards
# "normalization_vector_fromstandards.csv" --> normalization factor from standards
# "sampledata_normalizedtostandards.csv" --> raw data normalized to standards
# "sampledata_normalized_withRsq_and_CV.csv" --> standard-normalized raw data with Rsq and CV values
# "RsqCVfiltered_normalized_sampledata.csv" --> filtered metabos by Rsq and CV values
# "normalizationvector_fromtotalpolars.csv" --> normalization factor from metabolites
# "finaldata.csv" --> final data (normalized to polars) with Rsq and CV values
#
# Sample names have to match across all CSVs (so you may have to run the first section of code for your C18 and HILIC runs separately, then combine manually)
#
#
#####

print("Where can I find your Tracefinder .csv's?  (e.g.: /Documents/MS Data/202201208 (can use quotation marks)")
wd <- readline(prompt="Working Directory/Folder:")
wd <- gsub('[\"]', '', wd) #remove any slashes in the wrong direction
wd <- gsub('[\\"]', '/', wd) #replace doubled slashes (from Windows to R input) to single slashes
wd <- gsub('~', '', wd) #remove any tildas at the beginning
if (substr(wd,nchar(wd), nchar(wd)) != "/") # look if there is a backslash at the end of the working directory input, if not then add one
{
  wd = paste(wd,"/",sep="")
}
setwd(wd)
outputpath <- paste(wd, "RscriptOutput/", sep = "") #makes output folder = Rscriptoutput
if (file.exists(outputpath) == FALSE) #If the output folder doesn't exist, make one
{
  dir.create(outputpath)
}


# # Combine raw metabolite CSVs from Tracefinder ---------------------------------------------
i <- 1
filenames <- list.files(pattern = "*.csv");
sample <- read.csv(filenames[i], stringsAsFactor=FALSE)
df <- matrix(nrow = nrow(sample), ncol = length(filenames))
row.names(df) <- sample$Filename
columnnames <- matrix(nrow=length(filenames),1)

for (i in 1:length(filenames))
{
  data <- read.csv(filenames[i],stringsAsFactor=FALSE)
  df[,i] <- data$Area
  columnnames[i,] <- data[1,1]
}

colnames(df) <- columnnames[,1]

View(df)

print("combined data is now displayed.")
input <- readline(prompt="should I remove any rows from subsequent analysis? (y/n):")
if (input == "y")
{
  i <- 0
  while(i == 0)
  {
    rowstoremove <- unlist(strsplit(readline(prompt="which rows should be removed?  e.g. abc def (space separated list of unique identifiers:"), " "))
    for (i in 1:length(rowstoremove))
    {
      temp <- grepl(rowstoremove[i],rownames(df))
      if (i == 1)
      {
        rowlogical <- temp
      } else {
        rowlogical <- rowlogical + temp
      }
    }
    rows <- df[as.logical(rowlogical)]

    print(rownames(df)[as.logical(rowlogical)])
    yn <- readline(prompt="are these all the rows to be removed?  (y/n):")
    if (yn == "y")
    {
      i <- 1
    } else {
      i <- 0
    }
  }

  df <- df[!rowlogical,]

} else {

}

input <- readline(prompt="Combined data matrix is now displayed.  continue? (y/n):")
if (input == "y")
{

} else {
  stop()
}

print("Saving data as MSdataframe.csv...")
write.csv(df,file=paste(outputpath,"MSdataframe.csv", sep=""))


# Extract and normalize raw data and pools/dilutions to standards ----------------------------------------------------------------

print("reading MSdataframe.csv...")
rawdata <- read.csv(paste(outputpath,"MSdataframe.csv", sep=""), row.names = 1)
rawdata[rawdata == "N/F" ] <- NA
for (i in 1:length(rawdata)) #loop that converts to numeric
{
  suppressWarnings(class(rawdata[,i]) <- "numeric") #this step can generate errors when N/F is coerced to NA
}

# Find 1x pool(s)
i <- 0
while(i == 0)
{
  print("Input text to find names of 1x pool injections for CV calculation.  Space-separated list (e.g.: pool1xa pool1xb pool1xc), at least 3 names")
  pool1x <- unlist(strsplit(readline(prompt="pool1x names:"), " "))
  rm(rowlogical)
  for (i in 1:length(pool1x))
  {
    temp <- grepl(pool1x[i],rownames(rawdata))
    if (i == 1)
    {
      rowlogical <- temp
    } else {
      rowlogical <- rowlogical + temp
    }
  }
  pool1xs <- rownames(rawdata)[as.logical(rowlogical)]

  print(pool1xs)
  yn <- readline(prompt="are these all the pools?  (y/n):")
  if (yn == "y")
  {
    i <- 1
  } else {
    i <- 0
  }
}

# Tell me what text can identify the pool dilutions you want to use for Rsq calculations?
# Whatever text you pick must select out ONLY 1 of each dilution (e.g. all "1xa" injections).

curve <- as.numeric(unlist(strsplit(readline(prompt="What pool dilutions were included in the experiment?  (e.g. 1 0.33 0.1) as a space-separated list:"), " ")))

i <- 0
# invalid curve:
if (length(curve) <3|length(curve)>4)
{
  print("invalid standard curve")
  stop()
}

# 3 point curve:
while (i == 0)
{
  if (length(curve) == 3)
  {
    print("What text should I use to find dilution injections?  (e.g. pool1x pool033x pool01x) space-separated list")
    dilutions <- unlist(strsplit(readline(prompt="dilution names:"), " "))
    print("What text should I use to find injections for Rsq dilutions?  (e.g. pool1xa pool033xa pool01xa) space-separated list")
    Rsqdilutions <- unlist(strsplit(readline(prompt="dilution names:"), " "))
    poolmatrix <- matrix(0,nrow = 3, ncol = 3)
    poolmatrix[,1] <- curve
    poolmatrix[,2] <- dilutions
    poolmatrix[,3] <- Rsqdilutions
    colnames(poolmatrix) <- c("curve", "names","RsqDilutionNames")
    print("Pool curve dilutions:")
    print(poolmatrix)

    #pull the names of the Rsq dilutions
    rm(rowlogical)
    for (i in 1:length(poolmatrix[,3]))
    {
      temp <- grepl(poolmatrix[i,3],rownames(rawdata))
      if (i == 1)
      {
        rowlogical <- temp
      } else {
        rowlogical <- rowlogical + temp
      }
    }
    Rsqpools <- rownames(rawdata)[as.logical(rowlogical)]

    print("Rows to use for Rsq:")
    print(Rsqpools)

    input <- readline(prompt="are curve names and values correct? (y/n):")
    if (input == "y")
    {
      i <- 1
    }
    else
    {
      i <- 0
    }
  }
}

# 4 point curve:
while (i == 0)
{
  if (length(curve) == 4)
  {
    print("What text should I use to find dilution injections?  (e.g. pool1x pool033x pool02x pool01x) space-separated list")
    dilutions <- unlist(strsplit(readline(prompt="dilution names:"), " "))
    print("What text should I use to find injections for Rsq dilutions?  (e.g. pool1xa pool033xa pool01xa) space-separated list")
    Rsqdilutions <- unlist(strsplit(readline(prompt="dilution names:"), " "))
    poolmatrix <- matrix(0,nrow = 4, ncol = 3)
    poolmatrix[,1] <- curve
    poolmatrix[,2] <- dilutions
    poolmatrix[,3] <- Rsqdilutions
    colnames(poolmatrix) <- c("curve", "names","RsqDilutionNames")
    print("Pool curve dilutions:")
    print(poolmatrix)

    input <- readline(prompt="are curve names and values correct? (y/n):")
    if (input == "y")
    {
      i <- 1
    }
  }
}

write.csv(poolmatrix, file=paste(outputpath,"poolmatrix.csv", sep=""))

standardcols <- grepl("13C",colnames(rawdata))|grepl("15N",colnames(rawdata))|grepl("D3",colnames(rawdata))|grepl("D4",colnames(rawdata))|grepl("aminopterin",colnames(rawdata))
standards <- rawdata[,standardcols] #pull columns matching conditions above
rownames(standards) <- rownames(rawdata)
print("Identified standards metabolites in samples:")

blankpoolrows <- grepl("blank",rownames(standards))|grepl("uAA",rownames(standards))|grepl("QRESS",rownames(standards)) #identify blanks
standards <- standards[!blankpoolrows,] #remove blanks

print("removing any standards metabolites with >10% NA...")
filteredstandards <- standards[, which(colMeans(!is.na(standards)) > 0.90)] #remove any standards with >10% NA

print("writing standards into standards.csv...")
write.csv(filteredstandards,file=paste(outputpath,"filteredstandards.csv", sep=""))
View(filteredstandards)
input <- readline(prompt="Filtered Standards are now displayed (prior to dilution normalization.  continue? (y/n):")
if (input == "y")
{

} else {
  stop()
}

filteredstandards <- read.csv(paste(outputpath,"filteredstandards.csv", sep=""), row.names = 1)
#multiply standards by dilution factor in pool dilution rows
for (i in 2:length(poolmatrix[,2])) #start in 2nd column
{
  templogical <- grepl(poolmatrix[i,2], rownames(filteredstandards))
  filteredstandards[templogical, ] <- filteredstandards[templogical, ]*(1/as.numeric(poolmatrix[i,1]))
}

View(filteredstandards)
input <- readline(prompt="Filtered Standards are now displayed (POST-dilution normalization.  continue? (y/n):")
if (input == "y")
{
} else {
  stop()
}


colaverages <- colSums(filteredstandards, na.rm = TRUE)/nrow(filteredstandards) #get column averages
norm_standards <- sweep(filteredstandards, 2, colaverages, FUN = '/') #divide standard matrix by average
norm_vect <- rowSums(norm_standards, na.rm = TRUE)/ncol(filteredstandards) #calculate normalization factor into norm_vect

write.csv(norm_vect,file=paste(outputpath,"normalization_vector_fromstandards.csv", sep=""))
View(norm_vect)


if (sum(norm_vect>1.25|norm_vect<0.75) > 0)
{
  print("********")
  print("********")
  print("WARNING: the following rows have abnormal norm factors (<0.75 or >1.25):")
  print(rownames(filteredstandards)[norm_vect>1.2|norm_vect<0.8])
}
input <- readline(prompt="Normalization vector based on internal standards is now displayed.  continue? (y/n):")
if (input == "y")
{

} else {
  stop()
}

sampledata <- rawdata #pull specified rows for analysis
sampledata_normtostandards <- sampledata

####

standard_norm_index <- match(names(norm_vect),rownames(sampledata_normtostandards)) #matches norm_vect names to the sampledata names
for (j in 1:length(standard_norm_index))
{
    sampledata_normtostandards[standard_norm_index[j],] <- sampledata[standard_norm_index[j],]/norm_vect[j] #normalizes each row in sample data to the norm_vect valuem, based on which row is indicated in standard_norm_index
}

print("writing normalized data to sampledata_normalizedtostandards.csv")

write.csv(sampledata_normtostandards,file=paste(outputpath,"sampledata_normalizedtostandards.csv", sep=""))

View(sampledata_normtostandards)
input <- readline(prompt="Normalized data is now displayed.  continue? (y/n):")
if (input == "y")
{

} else {
  stop()
}

# CV and RSq calculation --------------------------------------------------

sampledata_normtostandards <- read.csv(paste(outputpath,"sampledata_normalizedtostandards.csv", sep=""), row.names = 1)
sampledata_normtostandards[sampledata_normtostandards == "N/F" ] <- NA
sampledata_normtostandards[sampledata_normtostandards == "#VALUE!" ] <- NA

#CV - to calculate, find stdev of columns/average of columns
cv <- numeric(length(sampledata_normtostandards))

#find the 1x pool rows based on user inputs
rm(rowlogical)
for (i in 1:length(pool1x))
{
  temp <- grepl(pool1x[i],rownames(sampledata_normtostandards))
  if (i == 1)
  {
    rowlogical <- temp
  } else {
    rowlogical <- rowlogical + temp
  }
}
pool1xdata <- sampledata_normtostandards[as.logical(rowlogical),]


#calculate cvs for each column
for (i in 1:length(pool1xdata))
{
  if (sum(!is.na(pool1xdata[,i])) > 0)
  {
    cv[i] <- sd(pool1xdata[,i],na.rm = TRUE)/mean(pool1xdata[,i], na.rm = TRUE) #calculate cvs for each column by pools
  }
}
cv[is.na(cv)] <- 0 #replace values with only 1 pool value

#Rsq - pull user specified columns

Rsq <- numeric(length(sampledata_normtostandards))

rm(rowlogical)
for (i in 1:length(poolmatrix[,3]))
{
  temp <- grepl(poolmatrix[i,3],rownames(sampledata_normtostandards))
  if (i == 1)
  {
    rowlogical <- temp
  } else {
    rowlogical <- rowlogical + temp
  }
}
#poolcurvedata <- sampledata_normtostandards[as.logical(rowlogical),] #this old code is problematic, it takes the rows without sorting them
#new code to take the pool data by order:
temp <- sampledata_normtostandards[as.logical(rowlogical),]
poolcurvedata <- temp

for (i in 1:length(poolmatrix[,3])) #this for loop pulls the pool curve injections for Rsq calculation
{
  poolcurvedata[i,] <- temp[grepl(poolmatrix[i,3],rownames(temp)),]
  rownames(poolcurvedata)[i] <- poolmatrix[i,3]
}

for (i in 1:length(poolcurvedata)) #calculate the Rsq values
{
  rsq <- function(x, y) summary(lm(y~x))$r.squared
  if (sum(is.na(poolcurvedata[,i])) == 0)
  {
    Rsq[i] <- rsq(curve, poolcurvedata[,i])
  }
}

# save calculations
temp <- sampledata_normtostandards
temp[nrow(temp)+1,] <- t(Rsq)
temp[nrow(temp)+1,] <- t(cv)
names <- rownames(temp)
names[length(names)-1] <- "Rsq"
names[length(names)] <- "CV"
rownames(temp) <- names
write.csv(temp,file=paste(outputpath,"sampledata_normalized_withRsq_and_CV.csv", sep=""))

View(temp)

input <- readline(prompt="Data is now displayed with Rsq and CV.  continue? (y/n):")
if (input == "y")
{

} else {
  stop()
}

# Filter out metabolites by Rsq and CV --------------------------------------
#keep if Rsq >0.95; CV <0.3

Rsq <- readline(prompt="What cut-off should be used for Rsq? (e.g. 0.95):")
CV <- readline(prompt="What cut-off should be used for CV? (e.g. 0.30):")

write.csv(c(Rsq, CV),file=paste(outputpath,"Rsq_CV.csv", sep=""))

filtered_sampledata <- read.csv(paste(outputpath,"sampledata_normalized_withRsq_and_CV.csv", sep=""), row.names = 1)
for (i in 1:ncol(filtered_sampledata))
{
  if ((filtered_sampledata[(nrow(filtered_sampledata)-1),i]<Rsq)|(filtered_sampledata[(nrow(filtered_sampledata)-0),i]>CV)) #filter out if Rsq <0.949 or CV > 0.301
  {
    filtered_sampledata[,i] <- NA #assign NA to all values if didn't pass Rsq and CV filtering
  }
}

#subset those that were marked OK by last step (take rows that have at least some good values)
filtered_sampledata <- filtered_sampledata[,colSums(is.na(filtered_sampledata)) != nrow(filtered_sampledata)] 

#remove standards
standardcols <- grepl("13C",colnames(filtered_sampledata))|grepl("15N",colnames(filtered_sampledata))|grepl("D3",colnames(filtered_sampledata))|grepl("D4",colnames(filtered_sampledata))|grepl("aminopterin",colnames(filtered_sampledata))
filtered_sampledata <- filtered_sampledata[,!standardcols]

write.csv(filtered_sampledata,file=paste(outputpath,"RsqCVfiltered_normalized_sampledata.csv", sep=""))

View(filtered_sampledata)
input <- readline(prompt="Data filtered by Rsq and CV cut-offs is now displayed.  continue? (y/n):")
if (input == "y")
{
  
} else {
  stop()
}

# Normalization to total metabolites -------------------------------------------

finaldata <- filtered_sampledata[1:(nrow(filtered_sampledata)-2),] #remove Rsq and CV rows
finaldata <- finaldata[!grepl("pool",rownames(finaldata)),] #pull only rows that are sample or pool
averages <- colSums(finaldata, na.rm = TRUE)/nrow(finaldata) #get column averages
norm_data <- sweep(finaldata, 2, averages, FUN = '/') #divide sample matrix by averages
data_norm_vect <- rowSums(norm_data, na.rm= TRUE)/ncol(norm_data) #calculate normalization factor into data_norm_vect

write.csv(data_norm_vect,file=paste(outputpath,"normalizationvector_fromtotalpolars.csv", sep=""))

View(data_norm_vect)
if (sum(data_norm_vect>1.25|data_norm_vect<0.75) > 0)
{
  print("********")
  print("********")
  print("WARNING: the following rows have abnormal norm factors (<0.75 or >1.25):")
  print(rownames(finaldata)[data_norm_vect>1.25|data_norm_vect<0.75])
}
input <- readline(prompt="Norm vector to total polars now displayed.  continue? (y/n):")
if (input == "y")
{
  
} else {
  stop()
}

finalnormdata <- finaldata
polars_norm_index <- match(names(data_norm_vect),rownames(finaldata)) #matches norm_vect names to the sampledata names
for (j in 1:length(polars_norm_index))
{
  finalnormdata[polars_norm_index[j],] <- finaldata[polars_norm_index[j],]/data_norm_vect[j] #normalizes each row in sample data to the data_norm_vect valuem, based on which row is indicated in standard_norm_index
}

#for (i in 1:nrow(finaldata)) #this for loop normalizes to standards in norm_vector
#{  
#  finalnormdata[i,] <- finaldata[i,]/data_norm_vect[i]
#}

#add back Rsq and CV data
finalnormdata[nrow(finalnormdata)+1,] <- filtered_sampledata[nrow(filtered_sampledata)-1,]
finalnormdata[nrow(finalnormdata)+1,] <- filtered_sampledata[nrow(filtered_sampledata),]
write.csv(finalnormdata,file=paste(outputpath,"finaldata.csv", sep=""))

View(finalnormdata)
print("Final data is displayed")
