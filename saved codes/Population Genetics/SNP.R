
install.packages("GWASExactHW")
library("GWASExactHW")
#GWASExactHW
GenotypeCounts <- data.frame(nAA = c(112, 189), nAa = c(195,157), naa = c(223,184))
#hardy wienberg
GenotypeCounts <- as.vector(c(112, 195, 223))
myTest <- HWExact(GenotypeCounts) #dataframe
names(myTest) <- rownames(GenotypeCounts)
myTest
plot(-log10(myTest), type="b", ylab="-log10(p-value)", main="HWE
p-values for SNPS")
abline(h=-log10(0.05), col="red")
sum(myTest<0.05)
names(myTest)[which(myTest<0.05)]
library("HardyWeinberg")
a <- HWChisq(GenotypeCounts)
rownames(GenotypeCounts) <- c('cancer', 'normal')
s <- apply(GenotypeCounts, 1,  HWChisq)
if(s$cancer$pval < 0.05) {print('notHW')} else {print('HW')}
if(s$normal$pval < 0.05) {print('notHW')} else {print('HW')}


GenotypeCounts[1]
chisq.test(GenotypeCounts)
?cite
citation("HardyWeinberg")
