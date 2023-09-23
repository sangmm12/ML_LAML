
setwd('~/machine_learning/linkET')

library(linkET)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)

rf_risk <- read.csv('~/diff/ml_laml/rfrisk.csv')
rf_risk <- rf_risk[,-c(ncol(rf_risk))]
rf_group <- rf_risk[which(rf_risk$risk=='High'),]$SampleName



riskscore_geneset <- c("ALDOC","KLF9","PPM1F","CSTB","CYYR1","DDIT4","DGAT2","IFITM3","PMM1","SASH1","SH2D3C","SUSD3","ANXA4","ARAP1","CYC1","DHRS9","ECHDC3","FBXO6","FKBP5","HIP1","HPCAL1","IL6R","LTB4R","PLCB2","RAB3D","RGL2","RPS6KA1","S100A13","TRIB1","ADGRE5","PELO","AVPI1","CYP19A1","HLX","PKN1","PPM1N","SHARPIN","SLA","DGAT1","GAS6","KRT5","AK1","CBLN3","CCND3","ECE1","EHBP1L1","GADD45A","GLTP","GNG11","ICAM3","IGF2R","ISG20","LSP1","MAD2L1BP","MAST3","MYL6","MZT2A","NEDD9","NFKBIL1","OPTN","OSGIN1","PCGF5","PF4","PIM1","PRDX5","RELB","RHOC","RPL3L","SEC14L1","SESN2","SH3BP5","SH3TC1","SIAH2","SLC10A3","SRXN1","TFE3","TMEM63C","TREML2","TUBA4A","TWIST2")


##################expression
CIBER <- read.csv('~/machine_learning/data/tcga/survival_out.csv',header=T,check.names = F)


ssgsva <- load('SSGSVA.Rds')



mantel <- mantel_test(ssgsva, out,
                      spec_'High'ect = list(ssGSVA = 1:1)) %>% 
mutate(rd = cut(r, breaks = c(-1, 0.2, 0.4, 1),
                labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
       pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


pdf(paste('High','.pdf',sep=''),height=10,width=10)
qcorrplot(correlate(out), type = "lower", diag = FALSE) 
  geom_square() 
  geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu")) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "R", order = 3))
dev.off()
