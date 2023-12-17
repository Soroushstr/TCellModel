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
