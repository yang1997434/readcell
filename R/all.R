#' read file
#'
#' @param file_site File name or folder
#' @param sample T or F
#' @param data_type single cell data type,such as:rds h5 loom 10X
#'
#' @return seurat object
#' @export
#'
#' @examples read_file(file_site,sample=T,data_type='rds')
read_file<-  function(file_site,
                      sample=T,
                      data_type='rds'){
  if (sample==T){
    if (data_type=='txt'){
      cat("====read txt(txt.gz) file start==== \n")
      seaurat <- read.table(gzfile(file_site),h = T,sep='\t',row.names = 1,check.names = F)
      seaurat <- CreateSeuratObject(seaurat,min.cell=3,min.features=200)
      cat("====seaurta file manufacture finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='rds'){
      cat("====read rds file start==== \n")
      seaurat <- readRDS(file_site)
      cat("====seaurta file loading finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='h5'){
      cat("====read 10X h5 file start====\n")
      seaurat <-Read10X_h5(file_site, use.names = TRUE, unique.features = TRUE)
      cat("====read 10X h5 file loading finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='10X'){
      cat("====read 10X file start==== \n")
      seaurat <- Read10X(data.dir = paste0(file_site,'/'))
      cat("====read 10X file loading finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='csv'){
      cat("====read csv file start==== \n")
      seaurat <- read.table(file_site,sep=',',header = T,row.names = 1)
      seaurat <- CreateSeuratObject(seaurat,min.cell=3,min.features=200)
      cat("====seaurta file manufacture finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='tsv'){
      cat("====read tsv file start==== \n")
      seaurat <- read.table(file_site,sep='\t',header = T,row.names = 1)
      seaurat <- CreateSeuratObject(seaurat,min.cell=3,min.features=200)
      cat("====seaurta file manufacture finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='h5ad'){
      cat("====read h5ad file start==== \n")
      seaurat <- readH5AD(file_site,reader = "R")
      cat("====loading h5ad file finsish==== \n")
      return (seaurat)
    }
    else if (data_type=='loom'){
      cat("====read loom file start==== \n")
      seaurat <- connect(filename = file_site, mode = "r+")
      cat("====loading loom file finsish==== \n")
      return (seaurat)
    }
    else cat("====data type erro==== \n")
  }
  else if (sample==F){
    seaurat<- list()
    file_name=as.list(list.files(path=file_site))
    for (i in 1:length(file_name)){
      if (data_type=='txt'){
        cat("====read txt(txt.gz) file start==== \n")
        seaurat[[i]] <- read.table(gzfile(paste0(file_site,'/',file_name[[i]])),h = T,sep='\t',row.names = 1,check.names = F)
        seaurat[[i]] <- CreateSeuratObject(seaurat[[i]],min.cell=3,min.features=200)
        cat("====seaurta file manufacture finsish==== \n")
      }
      else if (data_type=='rds'){
        cat("====read rds file start==== \n")
        seaurat[[i]] <- readRDS(paste0(file_site,'/',file_name[[i]]))
        cat("====seaurta file loading finsish==== \n")
      }
      else if (data_type=='h5'){
        cat("====read 10X h5 file start==== \n")
        seaurat[[i]] <-Read10X_h5(paste0(file_site,'/',file_name[[i]]), use.names = TRUE, unique.features = TRUE)
        cat("====read 10X h5 file loading finsish==== \n")
      }
      else if (data_type=='10X'){
        cat("====read 10X file start==== \n")
        seaurat[[i]] <- Read10X(data.dir = paste0(file_site,'/',file_name[[i]],'/'))
        cat("====read 10X file loading finsish==== \n")
      }
      else if (data_type=='csv'){
        cat("====read csv file start==== \n")
        seaurat[[i]] <- read.table(paste0(file_site,'/',file_name[[i]]),sep=',',header = T,row.names = 1)
        seaurat[[i]] <- CreateSeuratObject(seaurat[[i]],min.cell=3,min.features=200)
        cat("====seaurta file manufacture finsish==== \n")
        return (seaurat)
      }
      else if (data_type=='tsv'){
        cat("====read tsv file start==== \n")
        seaurat[[i]] <- read.table(paste0(file_site,'/',file_name[[i]]),sep='\t',header = T,row.names = 1)
        seaurat[[i]] <- CreateSeuratObject(seaurat[[i]],min.cell=3,min.features=200)
        cat("====seaurta file manufacture finsish==== \n")
        return (seaurat)
      }
      else cat("====data type erro==== \n")
    }
    return (seaurat)
  }
  else cat("====sample type is undefined==== \n")
}
#' plolt
#'
#' @param genes marker gene
#'
#' @return vlnplot Featureplot
#' @export
#'
#' @examples jd(c("CD47","PECAM"))
jd<-function(genes){
  res = length(genes)
  if (res < 5) {
    if (res<=3){
      options(repr.plot.width=14, repr.plot.height=5*1, repr.plot.res = 200)
      print(FeaturePlot(dataobj,features = genes,pt.size = pt.size,order = T,ncol=3))
      options(repr.plot.width=14, repr.plot.height=5*1, repr.plot.res = 200)
      print(VlnPlot(dataobj,features = genes,pt.size = 0,ncol=1))
    }
    else {
      options(repr.plot.width=14, repr.plot.height=5*1, repr.plot.res = 200)
      print(FeaturePlot(dataobj,features = genes,pt.size = pt.size,order = T,ncol=3))
      options(repr.plot.width=14, repr.plot.height=5*2, repr.plot.res = 200)
      print(VlnPlot(dataobj,features = genes,pt.size = 0,ncol=1))
    }

  }
  else {
    options(repr.plot.width=14, repr.plot.height=5*1, repr.plot.res = 200)
    print(FeaturePlot(dataobj,features = genes,pt.size = pt.size,order = T,ncol=3))
    options(repr.plot.width=14, repr.plot.height=5*(res-3), repr.plot.res = 200)
    print(VlnPlot(dataobj,features = genes,pt.size = 0,ncol=1))
  }
}

#' file merge
#'
#' @param seaurat object
#' @param file_site  folder
#'
#' @return seurat object
#' @export
#'
#' @examples file_merge(data,file_site)
file_merge<-function(seaurat,file_site){
  file_name=as.list(list.files(path=file_site))
  for (i in 1:length(file_name)){
    seaurat[[i]]@meta.data$orig.ident <- file_name[[i]]
    names(file_name) = seq(from=1, to=length(file_name), by=1)
    seaurat[[i]]@meta.data$sample_lable <- names(file_name)[i]
  }
  j=1
  while (j <= (length(seaurat)-1)){
    if (j==1){
      seaurat_merge<-merge( seaurat[[j]],seaurat[[j+1]])
      j=j+1
    }
    else {
      seaurat_merge<-merge(seaurat_merge,seaurat[[j+1]])
      j=j+1
    }
  }
  return (seaurat_merge)
}
#' load data
#'
#' @param file_site File name or folder
#' @param sample T or F
#' @param data_type single cell data type,such as:rds h5 loom 10X
#'
#' @return seurat object
#' @export
#'
#' @examples load_data(file_site,sample=T,data_type='rds')
load_data<-function(file_site,
                    sample=T,
                    data_type='rds'){
  if (sample==T){
    data<-read_file(file_site,sample=T,data_type=data_type)
    return (data)
  }
  if (sample==F){
    data<-read_file(file_site,sample=F,data_type=data_type)
    data<-file_merge(data,file_site=file_site)
    return (data)
  }
}
