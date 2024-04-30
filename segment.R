library(cghFLasso)

"segment" <- function(data)
{
   seg<-cghFLasso(data)

   segdata<-seg$Esti.CopyN

   out.file="seg.txt"

   write.table(segdata, file=out.file, sep="\t")

}