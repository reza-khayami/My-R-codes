OddsRatio <- function(a,b,c,d, r){
  # 
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

  return(data.frame(Type = c("Odds Ratio 1","Odds Ratio 2",
                             "Relative Risk 1", "Relative Risk 2" ),
                    Value = c(OR1, OR2, RR1, RR2),
                    CI_95 = c(paste0(ORLCI1, "-", ORUCI1),
                            paste0(ORLCI2, "-", ORUCI2),
                           paste0(RRLCI1, "-", RRUCI1),
                           paste0(RRLCI2, "-", RRUCI2)),
                    Status = c(OR1status, OR2status,
                               RR1status, RR2status)))
}


