# readcell
The file can be read using the load_data function with three parameters
```{r}
data<-load_data(file_site,sample=sample,data_type=data_type)
```
File_site: indicates the file name of a single sample file and the folder name of multiple copies
Sample =sample: sample=T indicates that the file is a single sample and F is a multi-sample file
The data_type parameter is the read data type. You can set 'TXT' 'RDS' 'H5' '10X' 'CSV' 'TSV' 'Loom' 'H5AD'
Note: 'loom' 'H5ad' does not have multiple copy options
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
It is used to draw scenes of Featureplot and Vlnplot of Marker gene, simplify code blocks, and automatically plan the size of the picture with the number of marker gene
```{r}
Jd(c(“CD79A”,”PECAM”))
```
![caption]("file:///C:/Users/DELL/Downloads/plot.pdf")
