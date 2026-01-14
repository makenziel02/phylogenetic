### raydisc commands for ancestral reconstruction of morphological characters

### load libraries ###
library(ape)
library(evobiR)
library(picante)
library(corHMM)
library(phytools)
library(geiger)
library("phangorn")

### set directory ###
setwd("~/Desktop/reu_ancestral")

### load in table with habitat type data ###
Oph_Morph_Char1 <- read.csv(file="Comparative Morphology - Scoring - Character1.csv")
Oph_Morph_Char2 <- read.csv(file="Comparative Morphology - Scoring - Character2.csv")
Oph_Morph_Char3 <- read.csv(file="Comparative Morphology - Scoring - Character3.csv")
Oph_Morph_Char4 <- read.csv(file="Comparative Morphology - Scoring - Character4.csv")
Oph_Morph_Char5 <- read.csv(file="Comparative Morphology - Scoring - Character5.csv")
Oph_Morph_Char6 <- read.csv(file="Comparative Morphology - Scoring - Character6.csv")
Oph_Morph_Char7 <- read.csv(file="Comparative Morphology - Scoring - Character7.csv")
Oph_Morph_Char8 <- read.csv(file="Comparative Morphology - Scoring - Character8.csv")

### read in tree file in nexus format and then plot ###
OphTree <- read.nexus(file="Ophio_concat_probe_20kb_woAspersus.nexus")
plot(OphTree)

### ladderize the tree; done so outgroups were at the top (my preferred aesthetic)
OphTree <- ladderize(OphTree, right=T)

### make sure tree is ultrametric if not already ###
OphTree <- chronopl(OphTree, lambda = 0.1)

### rayDISC command with equal rates model

recon_ER_Char1 <- rayDISC(OphTree, Oph_Morph_Char1, model="ER", node.states="marginal")
pdf(file="recon_ER_Char1.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char1$states,cex=0.5,pie.cex=1, width=10, height=90)
dev.off()

recon_ER_Char2 <- rayDISC(OphTree, Oph_Morph_Char2, model="ER", node.states="marginal")
pdf(file="recon_ER_Char2.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char2$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ER_Char3 <- rayDISC(OphTree, Oph_Morph_Char3, model="ER", node.states="marginal")
pdf(file="recon_ER_Char3.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char3$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ER_Char4 <- rayDISC(OphTree, Oph_Morph_Char4, model="ER", node.states="marginal")
pdf(file="recon_ER_Char4.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char4$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ER_Char5 <- rayDISC(OphTree, Oph_Morph_Char5, model="ER", node.states="marginal")
pdf(file="recon_ER_Char5.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char5$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ER_Char6 <- rayDISC(OphTree, Oph_Morph_Char6, model="ER", node.states="marginal")
pdf(file="recon_ER_Char6.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char6$states,cex=0.5,pie.cex=1, width=10, height=90)
dev.off()

recon_ER_Char7 <- rayDISC(OphTree, Oph_Morph_Char7, model="ER", node.states="marginal")
pdf(file="recon_ER_Char7.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char7$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ER_Char8 <- rayDISC(OphTree, Oph_Morph_Char8, model="ER", node.states="marginal")
pdf(file="recon_ER_Char8.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ER_Char8$states,cex=0.5,pie.cex=1, width=10, height=90)
dev.off()


### rayDISC command with ARD model

recon_ARD_Char1 <- rayDISC(OphTree, Oph_Morph_Char1, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char1.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char1$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char2 <- rayDISC(OphTree, Oph_Morph_Char2, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char2.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char2$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char3 <- rayDISC(OphTree, Oph_Morph_Char3, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char3.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char3$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char4 <- rayDISC(OphTree, Oph_Morph_Char4, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char4.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char4$states,cex=0.5,pie.cex=1, width=10, height=90)
dev.off()

recon_ARD_Char5 <- rayDISC(OphTree, Oph_Morph_Char5, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char5.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char5$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char6 <- rayDISC(OphTree, Oph_Morph_Char6, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char6.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char6$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char7 <- rayDISC(OphTree, Oph_Morph_Char7, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char7.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char7$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

recon_ARD_Char8 <- rayDISC(OphTree, Oph_Morph_Char8, model="ARD", node.states="marginal")
pdf(file="recon_ARD_Char8.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_ARD_Char8$states,cex=0.5,pie.cex=0.25, width=10, height=90)
dev.off()

### rayDISC command with SYM model

recon_SYM_Char1 <- rayDISC(OphTree, Oph_Morph_Char1, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char1.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char1$states,cex=0.5,pie.cex=1,width=10, height=20)
dev.off()

recon_SYM_Char2 <- rayDISC(OphTree, Oph_Morph_Char2, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char2.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char2$states,cex=0.5,pie.cex=1,width=10, height=20)
dev.off()

recon_SYM_Char3 <- rayDISC(OphTree, Oph_Morph_Char3, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char3.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char3$states,cex=0.5,pie.cex=1,width=10, height=20)
dev.off()

recon_SYM_Char4 <- rayDISC(OphTree, Oph_Morph_Char4, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char4.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char4$states,cex=0.5,pie.cex=0.25,width=10, height=20)
dev.off()

recon_SYM_Char5 <- rayDISC(OphTree, Oph_Morph_Char5, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char5.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char5$states,cex=0.5,pie.cex=0.25,width=10, height=20)
dev.off()

recon_SYM_Char6 <- rayDISC(OphTree, Oph_Morph_Char6, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char6.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char6$states,cex=0.5,pie.cex=0.25,width=10, height=20)
dev.off()

recon_SYM_Char7 <- rayDISC(OphTree, Oph_Morph_Char7, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char7.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char7$states,cex=0.5,pie.cex=0.25,width=10, height=20)
dev.off()

recon_SYM_Char8 <- rayDISC(OphTree, Oph_Morph_Char8, model="SYM", node.states="marginal")
pdf(file="recon_SYM_Char8.pdf", width = 10, height = 20)
plotRECON(OphTree,recon_SYM_Char8$states,cex=0.5,pie.cex=0.25,width=10, height=20)
dev.off()

### get model selection statistics; pick one with lowest AICc as best model
recon_ER_Char8$AIC
recon_ER_Char8$AICc

recon_ARD_Char8$AIC
recon_ARD_Char8$AICc

recon_SYM_Char8$AIC
recon_SYM_Char8$AICc



