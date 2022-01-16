################ All models##################
mySNP <- data.frame(id = c(1:623),
                    casco = c(rep(1, 448), rep(0, 175)),
                    SNP1 = factor(c(rep(1, 95), rep(2, 165), rep(3, 188),
                                    rep(1, 62), rep(2, 52), rep(3, 61)),
                                  levels = c(1, 2, 3), labels =
                                    c("CC", "CT", "TT")))
library("SNPassoc")
SNPanlys <- setupSNP(mySNP, 3, sep = "")
aaa <- association(casco ~ SNP1, SNPanlys, model = c("all"))
aaa
write.csv(aaa,file = 'teh8.csv', row.names = TRUE, col.names = TRUE, sep = ',')



###########multiplicative model##############

Data <- data.frame(SNP = "SNP",
                   C_EX_Case = (95*2)+165,
                   C_EX_Cont = (62*2+52),
                   T_NEN_Case = 188*2 + 165,
                   T_NEN_Cont = 61*2 + 52)
#X = table or list
OddsRatio <- function(Data, X, r){
  #Data tabulation
  a <- Data[, 2]
  b <- Data[, 3]
  c <- Data[, 4]
  d <- Data[, 5]
  AF <- data.frame(c(a/(a+c), (b/(b+d)), (c/(a+c)), (d/(b+d))))
  #calculating OR
  exposed <- (a) / (b)
  notexposed <- (c) / (d)
  OR1 <- round(exposed / notexposed, r)
  ORUCI1 <-  round(exp(log(OR1) + 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)), r)
  ORLCI1 <-  round(exp(log(OR1) - 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)), r)
  OR1status <- ifelse(ORUCI1 > 1 & ORLCI1 > 1 | ORUCI1 < 1 & ORLCI1 < 1,
                      paste("Significant"), paste("Not Significant"))
  
  OR2 <- round(notexposed / exposed, r)
  ORUCI2 <-  round(exp(log(OR2) + 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)), r)
  ORLCI2 <-  round(exp(log(OR2) - 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)), r)
  OR2status <- ifelse(ORUCI2 > 1 & ORLCI2 > 1 | ORUCI2 < 1 & ORLCI2 < 1,
                      paste("Significant"), paste("Not Significant"))
  
  #Calculating RR
  total1 <- ((a) + (b))
  total2 <- ((c) + (d))
  A1 <- ((b) / (a))
  A2 <-  ((d) / (c))
  RR1 <- round((a / total1) / (c / total2), r)
  D <- (A1 / total1) + (A2 / total2)
  RRUCI1 <- round(exp(log(RR1) + 1.96 * sqrt(D)), r)
  RRLCI1 <- round(exp(log(RR1) - 1.96 * sqrt(D)), r)
  RR1status <- ifelse(RRUCI1 > 1 & RRLCI1 > 1 | RRUCI1 < 1 & RRLCI1 < 1,
                      paste("Significant"), paste("Not Significant"))
  
  RR2 <- round((c / total2) / (a / total1), r)
  RRUCI2 <- round(exp(log(RR2) + 1.96 * sqrt(D)), r)
  RRLCI2 <- round(exp(log(RR2) - 1.96 * sqrt(D)), r)
  RR2status <- ifelse(RRUCI2 > 1 & RRLCI2 > 1 | RRUCI2 < 1 & RRLCI2 < 1,
                      paste("Significant"), paste("Not Significant"))
  
  #chi and fisher
  matrix1 <- matrix(c(a, c, b, d), nrow=2)
  matrix2 <- matrix(c(b, d, a, c), nrow=2)
  fisher1 <- fisher.test(matrix1)
  fisher2 <- fisher.test(matrix2)
  chiX <- chisq.test(matrix1)
  
  #Returning the Values
  if(X == "list"){
    results_as_list <- list(SNP = "SNP", OR_Disease = c(`Odds Ratio D` = OR1,
                                                        CI_95_ORD = paste0(ORLCI1, "-", ORUCI1),
                                                        Status_ORD = OR1status),
                            
                            OR_Normal = c(`Odds Ratio N` = OR2,
                                          CI_95_ORN = paste0(ORLCI2, "-", ORUCI2),
                                          Status_ORN = OR2status),
                            chiX2 = chiX,
                            fisher1 = fisher1,
                            fisher2 = fisher2,
                            
                            RR_Disease = c(`Relative Risk D` = RR1,
                                           CI_95_RRD = paste0(ORLCI1, "-", ORUCI1),
                                           Status_RRD = RR1status),
                            
                            RR_Normal = c(`Relative Risk N` = RR2,
                                          CI_95_ORN = paste0(ORLCI2, "-", ORUCI2),
                                          Status_RRN = RR2status))
    return(results_as_list)
  }
  if(X == "table"){
    results_as_table <- data.frame(Entry = c("C/EX", "T/NEX"),
                                   Normal = c(b, d),
                                   N_Percent = c(AF[2, ], AF[4, ]),
                                   Cancer = c(a, c),
                                   C_Percent = c(AF[1, ], AF[3, ]),
                                   OR = c(round(fisher1$estimate, r),
                                          round(fisher2$estimate, r)),
                                   lower = c(round(fisher1$conf.int[1], r),
                                             round(fisher2$conf.int[1], r)),
                                   upper = c(round(fisher1$conf.int[2], r),
                                             round(fisher2$conf.int[2], r)),
                                   p_value = c(fisher1$p.value,
                                               fisher2$p.value))
    
    return(results_as_table)
  }
  
  
}
SNP <- OddsRatio(Data,"table", 2)
write.csv(SNP,file = 'teh12.csv', row.names = TRUE, col.names = TRUE, sep = ',')


########## Hardy-Weinberg######
library("HardyWeinberg")
GenotypeCounts <- data.frame(AA = c(95, 62), AB = c(165,52), BB = c(188,61))
rownames(GenotypeCounts) <- c('cancer', 'normal')
colnames(GenotypeCounts) <- c('AA','AB','BB')
s <- apply(GenotypeCounts, 1,  HWChisq)
if(s$cancer$pval < 0.05) {print('case notHW')} else {print('HW')}
if(s$normal$pval < 0.05) {print('control notHW')} else {print('HW')}

