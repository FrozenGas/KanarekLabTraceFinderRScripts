# KanarekLab TraceFinder RScript

**Introduction**

The R script included here streamlines the normalization and analysis process for data outputted by Tracefinder.  It takes the output .csv's from Tracefinder and performs: data aggregation, normalization to internal standards, and normalization to total signal, with user input.

**System Requirements**

Software:
1. R and RStudio: https://posit.co/download/rstudio-desktop/
  I recommend using R version 4.2.0
  The R script works well in our lab on both Windows 10 and Mac OS Catalina

**Installation and Usage Guide**
1. Download the R Script to your computer
2. Open the RScript in RStudio
3. To start the script, click: "Code > Run Region > Run All"
   - This is critical to ensure the user inputs are processed correctly
4. Once the script begins executing, it will prompt user input in the Console.  See demo folder for an example.

**Intput**
1. Tracefinder ".csv" files containing peak areas.  Each ".csv" represents a different metabolite, with individual rows within each file for each sample.
   
**Outputs**

.csv outputs (will be in a folder called "RScriptOutput" in the working directory):

- "MSdataframe.csv" --> raw data of combined CSVs from TraceFinder
- "standards.csv" --> CSV of just the isoAAs, QRESS and aminopterin standards
- "normalization_vector_fromstandards.csv" --> mean-centred normalization factor from standards
- "sampledata_normalizedtostandards.csv" --> raw data normalized to normalization factor from standards
- "sampledata_normalized_withRsq_and_CV.csv" --> standard-normalized raw data with Rsq and CV values
- "RsqCVfiltered_normalized_sampledata.csv" --> filtered metabolites based on user-inputted thresholds for Rsq and CV values
- "normalizationvector_fromtotalpolars.csv" --> mean-centred normalization factor from total signal of all metabolites
- "finaldata.csv" --> final normalized data with Rsq and CV values

Expected run-time on a computer: <5 minutes

**Demo**

The included demo.zip file includes demo ".csv" files from Tracefinder, sample output, as well as the example R console text from a successful run of the code.  

Written by Alan Wong
Contact: alan.wong@childrens.harvard.edu
