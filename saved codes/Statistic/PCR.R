install.packages('pcr')

library(pcr)
file <- "D:/university files/MSc/Thesis/SLC30A10-3/results/expression.xlsx"

ct1 <- readxl::read_excel(file,
                          sheet = "SLC10")
ct1 <- ct1[,-3]
group_var <- rep(c('OVR', 'control'), each = 2)
res <- pcr_analyze(ct1,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'control')


ct2 <- readxl::read_excel(file,
                          sheet = "SLC3")
ct2 <- ct2[,-3]
res2 <- pcr_analyze(ct2,
                   group_var = group_var,
                   reference_gene = 'GAPDH',
                   reference_group = 'control')

ct3 <- readxl::read_excel(file,
                          sheet = "103")
ct3 <- ct3[,-4]
res3 <- pcr_analyze(ct3[,c(1,3)],
                    group_var = group_var,
                    reference_gene = 'GAPDH',
                    reference_group = 'control')

