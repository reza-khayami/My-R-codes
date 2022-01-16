mySNP <- function(table, AA,BB,AB, SNPimp , matrixn, iteration, seed, f1 ){
  
  #table : table containing SNPs
  #AA,BB,AB: Genotypes (ref first)
  #SNPimp: set TRUE for imputation
  #matrixn: number of datasets in mice
  #iteratuin: number of iterations for imputation
  #seed: set seed
  #f1: formula for regression, should have genotype
  #example: "Outcome~Genotype + ... "
  #Two columns named Genotype and Status are needed
  
  AB.BB <- paste0(AB,"-",BB)
  AA.AB <- paste0(AA,"-",AB)
  AA.BB <- paste0(AA,"-",BB)
  
  
  #prepare tables------------------------
  #Dominant
  table$Dominant <- ifelse(table$Genotype == AA, AA, AB.BB)
  table$Dominant <- factor(table$Dominant)
  table$Dominant <- relevel(table$Dominant, ref = AA)
  #Recessive
  table$Recessive <- ifelse(table$Genotype == BB, BB, AA.AB)
  table$Recessive <- factor(table$Recessive)
  table$Recessive <- relevel(table$Recessive, ref = AA.AB)
  
  #Overdominant
  table$Overdominant <- ifelse(table$Genotype == AB, AB, AA.BB)
  table$Overdominant <- factor(table$Overdominant)
  table$Overdominant <- relevel(table$Overdominant, ref = AA.BB)
  
  #log additive
  table$Log <- ifelse(table$Genotype == AA, 0, 
                      ifelse(table$Genotype == AB, 1,2))
  
  #Build models--------------------
  #univariate
  uniglm <- list(
    glmadd = glm(Status~Genotype, family = binomial(link = "logit"), data = table),
    glmdom = glm(Status~Dominant, family = binomial(link = "logit"), data = table),
    glmrec = glm(Status~Recessive, family = binomial(link = "logit"), data = table),
    glmovr = glm(Status~Overdominant, family = binomial(link = "logit"), data = table),
    glmlog = glm(Status~Log, family = binomial(link = "logit"), data = table),
    nmod = glm(Status~1, family = 'binomial', data = table) ##"null" mod
  )
  
  if(SNPimp != TRUE){impglm <- list(
    Additive = glm(as.formula(f1),family = binomial(link = "logit"), data = table),
    Dominant = glm(as.formula(gsub("Genotype","Dominant",f1)),family = binomial(link = "logit"), data = table),
    Recessive = glm(as.formula(gsub("Genotype","Recessive",f1)),family = binomial(link = "logit"), data = table),
    Overdominant = glm(as.formula(gsub("Genotype","Overdominant",f1)),family = binomial(link = "logit"), data = table),
    Log = glm(as.formula(gsub("Genotype","Log",f1)),family = binomial(link = "logit"), data = table)
  )} else  if(SNPimp == TRUE){
    #imputation-----------------------
    imp <- mice(table,
                m = matrixn,
                maxit = iteration,
                seed=seed)
    
    impglm <- list(
      Additive = with(data=imp,glm(as.formula(f1),family = binomial)),
      Dominant = with(data=imp,glm(as.formula(gsub("Genotype","Dominant",f1)),family = binomial)),
      Recessive = with(data=imp,glm(as.formula(gsub("Genotype","Recessive",f1)),family = binomial)),
      Overdominant = with(data=imp,glm(as.formula(gsub("Genotype","Overdominant",f1)),family = binomial)),
      Log = with(data=imp,glm(as.formula(gsub("Genotype","Log",f1)),family = binomial))
    )
  
    #summary of multi
    sumimp <-  lapply(impglm, function(x){summary(pool(x),conf.int = TRUE,conf.level = 0.95)})
    
    for (i in 1:length(sumimp)){
      rownames(sumimp[[i]]) <- sumimp[[i]]$term
    }
  }
    
    
  #get colsums for percentage in table-------------------------
  Sums <- colSums(table(table$Genotype,table$Status))
  tables <- list(
    Additive = table(table$Genotype,table$Status),
    Dominant = table(table$Dominant,table$Status),
    Recessive = table(table$Recessive,table$Status),
    Overdominant = table(table$Overdominant,table$Status),
    Log = data.frame(Control = table(table$Status)[1], Case = table(table$Status)[2]))
  
    models <- list(uniglm,impglm)
    names(models) <- c("uniglm","impglm")
    
    
    data <- list()
    for ( i in 1:5){
      data[[i]] <-  data.frame(Model = c(names(tables)[i],rep("",nrow(tables[[i]])-1)),
                               Genotype = if(i==5){"0,1,2"}else{rownames(tables[[i]])},
                               Cases = c(
                                 paste0(tables[[i]][,"Case"],
                                        " ",
                                        "(", round(tables[[i]][,"Case"]/Sums["Case"],2)*100,
                                        "%)")),
                               
                               Controls = c(
                                 paste0(tables[[i]][,"Control"],
                                        " ",
                                        "(", round(tables[[i]][,"Control"]/Sums["Control"],2)*100,
                                        "%)")),
                               
                               `OR (95%CI)` = if(i ==5){
                                 c(paste0(round(exp(coef(models[["uniglm"]][[i]]))[-1],2),
                                          " (",
                                          round(exp(confint(models[["uniglm"]][[i]])[-1,"2.5 %"]),2),
                                          "-",
                                          round(exp(confint(models[["uniglm"]][[i]])[-1,"97.5 %"]),2),
                                          ")"))}
                               else{
                                 c(1
                                   ,
                                   paste0(round(exp(coef(models[["uniglm"]][[i]]))[-1],2),
                                          " (",
                                          round(exp(confint(models[["uniglm"]][[i]])[-1,"2.5 %"]),2),
                                          "-",
                                          round(exp(confint(models[["uniglm"]][[i]])[-1,"97.5 %"]),2),
                                          ")"))},
                               P = if(nrow(tables[[i]])<2){round(summary(models[["uniglm"]][[i]])$coefficients[-1,"Pr(>|z|)"],5)
                               }else {c("",
                                        round(summary(models[["uniglm"]][[i]])$coefficients[-1,"Pr(>|z|)"],5))} ,
                               
                               `OR adj (95%CI)` = 
                                 if(nrow(tables[[i]])>2){
                                   c(1,paste0(round(exp(sumimp[[i]][2:3,"estimate"]),2),
                                              " (",
                                              round(exp(sumimp[[i]][2:3, "2.5 %"]),2),
                                              "-",
                                              round(exp(sumimp[[i]][2:3, "97.5 %"]),2),
                                              ")"))
                                 }else if(nrow(tables[[i]])==2){
                                   c(1,
                                     paste0(round(exp(sumimp[[i]][2,"estimate"]),2),
                                            " (",
                                            round(exp(sumimp[[i]][2, "2.5 %"]),2),
                                            "-",
                                            round(exp(sumimp[[i]][2, "97.5 %"]),2),")"))
                                 }else if(nrow(tables[[i]])<2){
                                   paste0(round(exp(sumimp[[i]][2,"estimate"]),2),
                                          " (",
                                          round(exp(sumimp[[i]][2, "2.5 %"]),2),
                                          "-",
                                          round(exp(sumimp[[i]][2, "97.5 %"]),2),")")
                                 },
                               
                               
                               Padj = if(nrow(tables[[i]])>2){
                                 c("",
                                   paste0(round(sumimp[[i]][2:3,"p.value"],5)))
                               }else if(nrow(tables[[i]])==2){
                                 c("",
                                   paste0(round(sumimp[[i]][2,"p.value"],5)))
                               }else if(nrow(tables[[i]]<2)){
                                 paste0(round(sumimp[[i]][2,"p.value"],5))
                               }
                               
                               
                               
                               
                               
                               
                               
      )
      
      
      
      
    }
    Results <- list(
      univariate= uniglm,
      multiple = impglm,
      OR.table = do.call(rbind.data.frame, data),
      imputation_table = imp)
    return(Results)
      
  }
  
  
  
  
  
  
  
  save(mySNP, file = "codes/mySNPfunction.rdata")
  
  
 
  
  