########    1. Load the packages    ########
library(flowCore)
library(flowStats)
library(flowViz)
library(ggcyto)
library(flowAI)


setwd("E:/university files/MSc/Thesis/SLC30A10-3/results/Flowcytometry/FlowCyto/FCS")
fs <- read.flowSet(path = "E:/university files/MSc/Thesis/SLC30A10-3/results/Flowcytometry/FlowCyto/FCS", pattern = ".fcs")

#####    2.Compensation    ####

spillover(fs)
fs_comp <-compensate(fs, spillover(fs)$SPILL)

####    3.Cleaning    ####

fs_clean <- flow_auto_qc(fs,second_fractionFR = 1)

####    4.Transformation    ####

trans.function <- function(x){
  trans <- estimateLogicle(x, colnames(x[,c(3, 5:9)]), m = 4.737)
  fs_clean_trans <- transform(x, trans)
}
fs_clean_trans <- fsApply(fs_clean, trans.function)

####    5.Visualise the results    ####

autoplot(fs_clean[[1]])
autoplot(fs_clean_trans[[1]])
autoplot(fs_clean_trans, x = "FSC-A", y = "SSC-H")
autoplot(fs_clean_trans, x="BL1-A", y="BL2-A", bins = 256)
autoplot(fs_clean_trans, x="Time", y="FSC-A", bins = 128)

#Automatic gating
#create the empty gating set
auto_gs<-GatingSet(HT29)

#cell gate
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC-A","SSC-A")))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC-A", y="SSC-A", "noneDebris_gate", bins=700) + cleanup

#singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC-A', y = 'FSC-H', "singlets", bins = 700) + cleanup

#Quad gate
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
BGquad_gate <- fsApply(fs_data, function(fr) openCyto:::.quadGate.seq(fr, gFunc="mindensity", channels =c('BL1-A', 'BL2-A')))
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("1","2","3","4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs, x = 'BL1-A', y = 'BL2-A', gs_get_pop_paths(auto_gs)[4], bins = 256)

fs_data <- gs_pop_get_data(auto_gs, "defaultEllipsoidGate") #get parent data
BGquad_gate <- fsApply(fs_data, function(fr) openCyto:::.quadGate.seq(fr, gFunc="mindensity",channels =c('BL1-A', 'BL2-A')))
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("1","2","3","4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)

p<-ggcyto(auto_gs,aes(x = 'BL1-A', y = 'BL2-A'), subset="singlets", arrange = FALSE)
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7])
p<- p + theme(strip.text = element_text(size = 7))
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p<- p  + cleanup
p

pdf('ss.pdf')
par(mfrow = c(3,4))
for (i in 1:11){
  plotDens(fs_clean_trans[[i]], channels = c(3,7),density.overlay = c(T,T),main = keyword(fs_clean_trans[[i]])$GUID.original )
}
dev.off()


#combined

combined_gs<-GatingSet(fs_clean_trans)

#cell gate
fs_data<- gs_pop_get_data(combined_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC-A","SSC-A")))
gs_pop_add(combined_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(combined_gs)
autoplot(combined_gs, x="FSC-A", y="SSC-A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(combined_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-H")))
gs_pop_add(combined_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(combined_gs)


qg <- quadGate("BL1-A"= 1.4, "BL2-A"= 1.7)
gs_pop_add(combined_gs, qg, parent = "singlets", name = c("-BL1+BL2", "+BL1+BL2", "+BL1-BL2", "-BL1-BL2"))
gs_get_pop_paths(combined_gs)
recompute(combined_gs)
autoplot(combined_gs, x = 'BL1-A', y = 'BL2-A', gs_get_pop_paths(combined_gs)[4:7], bins = 256)
