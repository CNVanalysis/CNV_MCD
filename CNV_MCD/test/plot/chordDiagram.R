rm(list = ls())
library(circlize)
"plotchord" <- function()
{
    file = '/.. /chord_data_file'
    data = read.table(file,header=TRUE)
    data = as.matrix(data)
    rownames(data) = c("CNV_MCD","FREEC","CNVnator","CNV_IFTV")
    colnames(data) = c("chr1",2:22)
    grid.col = c(CNV_MCD = "red", dpCNV="yellow",FREEC="Green",CNVnator="orange",CNV_IFTV="purple",chr1="grey","2"="grey","3"="grey","4"="grey","5"="grey","6"="grey","7"="grey","8"="grey","9"="grey","10"="grey","11"="grey","12"="grey","13"="grey","14"="grey","15"="grey","16"="grey","17"="grey","18"="grey","19"="grey","20"="grey","21"="grey","22"="grey")
    circos.par(start.degree = 180,clock.wise = TRUE)
    chordDiagram(data,grid.col = grid.col)
    circos.clear()

}


