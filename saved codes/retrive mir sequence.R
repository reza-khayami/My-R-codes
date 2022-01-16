setwd("E:/university files/MSc/Articles/Dr Kerachian/RS96/results/RRI/mirDB")
library(data.table)
files <- list.files(pattern = ".csv")
cnDT <- lapply(files, fread)
names(cnDT) <- files
setdiff(cnDT$C3total.csv$`miRNA Name`,cnDT$C2total.csv$`miRNA Name`)

miRNANames = cnDT$C1total.csv$`miRNA Name`
version=checkMiRNAVersion(miRNANames, verbose = TRUE)

names <- lapply(cnDT, function(x){x$`miRNA Name`})
result1 <-  lapply(names, function(x){miRNA_NameToAccession(x,version = "v22")})

NAs <- lapply(result1, function(x){which(is.na(x$Accession))})

result2 <-  lapply(result1, function(x){getMiRNASequence(x$Accession,targetVersion = "v22")})
NAs2 <- lapply(result2, function(x){which(is.na(x$miRNASequence_v22))})
NAs2

result3 <- result2

for (i in 1:6){
  result3[[i]]$miRNAName_v22 <- result1[[i]]$miRNAName_v22
}



result4 <- lapply(result3, function(x){x$fasta <- paste0(">", x$Accession,":", x$miRNAName_v22,  "\n", x$miRNASequence_v22)})

A <- function(x){write.table(x,paste0("A", names(x), ".txt"), quote = F, row.names = F, col.names = F)}
for (i in 1:6){
  A(result4[i])
}
