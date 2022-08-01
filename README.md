# readcell
The optional parameters are single sample or multiple copies. In a single sample, file_site corresponds to the file name, and multiple copies correspond to the folder. Sample =T indicates that the file is a single sample.
F is multi-sample, data_type parameter is read data type, you can set 'TXT' 'RDS' 'H5' '10X' 'CSV' 'TSV' 'loom' 'H5ad', note: 'loom' 'H5ad' does not have multiple books option
```{r}
suppressMessages(suppressWarnings({
    library(tidyverse)
    library(dplyr)
    library(stringr)
    library(Seurat)
    library(SummarizedExperiment)
    library(loomR)
    library(devtools)
}))  ###loading package
install_github("yang1997434/readcell")
library(readcell)
file_site='/opt/test/yangpeng/test'
sample=F
data_type='txt' ###set parameter
data<-load_data(file_site,sample=sample,data_type=data_type)
```
