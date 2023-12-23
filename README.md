# TCellModel
Reconstruction and analysis metabolic network mode for T cell subtypes

# Expression data preprocessing

## Affymetrix Data ##
Normalization and probe calling of affymetrix expression data, using Affy package and MAS5.0
Probes with _p_-value<0.04, 0.04< _p_-value<0.06, and _p_-value>0.06 are called _Present_, _Marginal_, and _Absent_, respectively

Here we considered Marginal and Present calls as expressed (1) and Absent calls as unexpressed (0)
```R
# loading required Affy library
library(Affy)

# use dataset's raw reads file
untar("GSE71566_RAW.tar")
eset = ReadAffy()
pData(eset)

# MAS5.0 normalization and visualization
eset.mas5 = mas5(eset,normalize = TRUE)
expset = exprs(eset.mas5)
expset.log2 = log(expset, 2)
boxplot(expset.log2, las=2)

# MAS5.0 Calling
calls.mas5 = mas5calls(eset)
calls = data.frame(exprs(calls.mas5))
pvals = data.frame(calls.mas5@assayData$se.exprs)
calls[calls=="P" | calls=="M"] = 1
calls[calls=="A"] = 0
head(calls)

# Writing csv output
write.csv(calls, file = "Calls_GSE71566.csv")
```

## Agilent Data ##
Agilent data mainly consist of images from the microarray chip, containing the information on th brightness of each probe for each sample. Data is analyzed using limma library.

Here we consider probes, present in at least half of the samples, and at least 40% brighter than negative control dark corners, as expressed (1) and the rest as unexpressed (0)
```R
# loading requried libraries
library(limma)
library(biomaRt)

# creating the mapping from agilent probe ID to gene ID, using BiomaRt library
mart = useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'agilent_sureprint_g3_ge_8x60k',
    'wikigene_description',
    'ensembl_gene_id',
    'gene_biotype',
    'external_gene_name'))

write.table(
  annotLookup,
  paste0('Human_agilent_sureprint_g3_ge_8x60k', gsub("-", "_", as.character(Sys.Date())), '.tsv'),
  sep = '\t',
  row.names = FALSE,
  quote = FALSE)

# Targets.txt file is a tab-delimited text file containing FileName, and if required, other columns for other sample properties
targets = readTargets('Targets.txt', sep = '\t')
project = read.maimages(targets, source = 'agilent', green.only = T, other.columns = 'gIsWellAboveBG')
project.bgcorrect <- limma::backgroundCorrect(project, method = 'normexp')
project.bgcorrect.norm <- normalizeBetweenArrays(project.bgcorrect, method = 'quantile')
Control <- project.bgcorrect.norm$genes$ControlType==1L

# using dark corners of the chip as the negative control
NegControl <- project.bgcorrect.norm$genes$ProbeName=='DarkCorner'
project.bgcorrect.norm.ctrl <- project.bgcorrect.norm[NegControl,]
NegValue <- colMeans(project.bgcorrect.norm.ctrl$E)
NoSymbol <- is.na(project.bgcorrect.norm$genes$external_gene_name)

# selecting probes, expressed in at least half of the samples (here this dataset had 8 samples)
IsExpr <- rowSums(project.bgcorrect.norm$other$gIsWellAboveBG > 0) >= 4
project.bgcorrect.norm.ctrl <- project.bgcorrect.norm[Control,]
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & IsExpr, ]

# probes with intensity at least 40% higher than the negative control's intensity
Expr40 <- project.bgcorrect.norm.filt$E > 1.4*NegValue
project.bgcorrect.norm.expr40 <- project.bgcorrect.norm[!Control & Expr40, ]

# generating and saving output calls
Expr40 = data.frame(Expr40)
Expr40$Gene = project.bgcorrect.norm.filt$genes$GeneName
Expr40$Probe_ID = project.bgcorrect.norm.filt$genes$ProbeName
head(Expr40)
Expr40[Expr40==T] = 1
Expr40[Expr40==F] = 0
head(Expr40)
write.csv(Expr40,file = "Calls_GSE60678.csv")
```

## Illumina Beadchip Data ##
This section uses limma library to read Illumina Beadchip data and conduct detection process to call _Present/Absent_ calls

## Clariom Data ##
This section uses oligo library to read Calriom chip and conduct base calling based on _paCalls_ function

```R
library(oligo)

celfiles <- list.files("CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
calls <- oligo::paCalls(rawData)
normData <- rma(rawData)
boxplot(normData)
expset = exprs(normData)
expset.log2 = log(expset, 2)
boxplot(expset.log2, las=2)
calls.mas5 = mas5calls(rawData)
calls = data.frame(exprs(calls.mas5))
pvals = data.frame(calls.mas5@assayData$se.exprs)
calls[calls=="P" | calls=="M"] = 1
calls[calls=="A"] = 0
head(calls)
write.csv(calls, file = "Calls_GSE71566.csv")
CALLS = paCalls(object=rawData)
rawData@assayData
```
