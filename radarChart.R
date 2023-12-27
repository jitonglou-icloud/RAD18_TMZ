rm(list = ls())
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

library(dplyr)

outPath <- "./out"
dir.create(outPath, showWarnings = FALSE)

### read data

data.wtd <- readxl::read_xlsx("../../PinAPL/N3_DMSO_PD20_avg_0.01_Sidak_sgRNAList.xlsx")

data.wtt <- readxl::read_xlsx("../../PinAPL/N3_TMZ_PD20_avg_0.01_Sidak_sgRNAList.xlsx")

data.kod <- readxl::read_xlsx("../../PinAPL/S18_DMSO_PD20_avg_0.01_Sidak_sgRNAList.xlsx")

data.kot <- readxl::read_xlsx("../../PinAPL/S18_TMZ_PD20_avg_0.01_Sidak_sgRNAList.xlsx")

pathways <- readxl::read_xlsx("../Short_Pathway_List.xlsx")


SetTest <- function(data.t, data.d, pathways) {
  # .t treatment
  # .d control
  
  # control guides
  data.t.control <- data.t %>% filter(sub("\\_.*", "", gene) == "NonTargeting")
  data.d.control <- data.d %>% filter(sub("\\_.*", "", gene) == "NonTargeting")
  
  data.t$Scaled_Log2_TD_norm <- log2(data.t$`fold change`)
  data.d$Scaled_Log2_TD_norm <- log2(data.d$`fold change`)
  
  ### construct Z-stat for each gene in each group

  data.t.gene <- data.t %>% filter(sub("\\_.*", "", gene) != "NonTargeting")
  data.d.gene <- data.d %>% filter(sub("\\_.*", "", gene) != "NonTargeting")
  
  data.t.gene <- data.t.gene %>% 
    group_by(gene) %>% 
    summarise(muhat = mean(Scaled_Log2_TD_norm),
              S2 = var(Scaled_Log2_TD_norm)) %>% 
    ungroup()
  
  data.d.gene <- data.d.gene %>% 
    group_by(gene) %>% 
    summarise(muhat = mean(Scaled_Log2_TD_norm),
              S2 = var(Scaled_Log2_TD_norm)) %>% 
    ungroup()
  
  ### construct two-sample Z-stat for each gene
  
  data.gene <- inner_join(data.t.gene, data.d.gene, by = "gene", suffix = c(".t", ".d"))
  
  data.gene <- data.gene %>% 
    mutate(muhat.g = muhat.t - muhat.d,
           S2.g = S2.t/10 + S2.d/10)
  
  ### For each pathway
  
  data.stats <- data.frame(pathway = colnames(pathways), muhat.p = NA, S2.p = NA, mu_c = NA)
  
  for (i in 1:ncol(pathways) ) {
    
    pathway <- pathways[, i] %>% unlist()
    pathway <- pathway[!is.na(pathway)]
    
    data.pathway <- data.gene %>% 
      filter(gene %in% pathway)
    
    ### combine Z-stats for the genes in a pathway
    
    muhat.pathway <- mean(data.pathway$muhat.g)
    S2.pathway <- mean(data.pathway$S2.g)/nrow(data.pathway)
    
    data.stats$muhat.p[i] <- muhat.pathway
    data.stats$S2.p[i] <- S2.pathway
   }
  
  return(data.stats)
}

`%notin%` <- Negate(`%in%`)

### two-sided p-value

pval <- function(mu, v2) {
  pnorm(-abs(mu/sqrt(v2)))*2
}

### one-sided p-value (alternative: "less")

pvall <- function(mu, v2) {
  pnorm((mu/sqrt(v2)))
}

########
data.stats.ko <- SetTest(data.t = data.kot, data.d = data.kod, pathways = pathways)
data.stats.wt <- SetTest(data.t = data.wtt, data.d = data.wtd, pathways = pathways)

data.stats.ko <- data.stats.ko %>% mutate(class = "KO PD20\nTMZ vs DMSO")
data.stats.wt <- data.stats.wt %>% mutate(class = "WT PD20\nTMZ vs DMSO")

### self-contained  test

data.selfcontained.ko <- data.stats.ko %>% 
  mutate(p.val = pval(muhat.p, S2.p))
data.selfcontained.wt <- data.stats.wt %>% 
  mutate(p.val = pval(muhat.p, S2.p))

### re-order the pathways
order0 <- c("BER", "TLS", "FA", "HR", "NHEJ", "NER", "Checkpoint Signaling",
            "Mitosis/Spindle Assembly Checkpoint", "Poly ADP Ribose Polymerases",
            "Nucleotide Metabolism","Template Switch"
            )

data.selfcontained.ko <- data.selfcontained.ko %>% arrange(factor(pathway, levels = order0))
data.selfcontained.wt <- data.selfcontained.wt %>% arrange(factor(pathway, levels = order0))



### radar plot

library(fmsb)

# sort the pathways

# radar plot of pvals

pathway.names <- data.selfcontained.ko$pathway

data.pval <- as.data.frame( -log10(matrix( c(data.selfcontained.ko$p.val, data.selfcontained.wt$p.val), nrow = 2, byrow = T )) )
colnames(data.pval) <- pathway.names
rownames(data.pval) <- c("KO PD20\nTMZ vs DMSO", "WT PD20\nTMZ vs DMSO")

np <- ncol(data.pval)
max(-log10(data.selfcontained$p.val))
min(-log10(data.selfcontained$p.val))
data.pval <- rbind(rep(80,np), rep(0,np), data.pval)

rgba <- function(...) {
  rgb(..., maxColorValue = 255)
}

colors_border=c( rgba(191,97,165,255), rgba(61,88,167,255) )
colors_in=c( rgba(191,97,165,100), rgba(61,88,167,100) )

########
setEPS()
postscript(file.path(outPath, "pvalRadar.eps"), width = 14, height = 9)

radarchart( data.pval  , axistype=1 ,
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey50", cglty=1, axislabcol="black", caxislabels=seq(0,80,20), cglwd=0.8,
            #custom labels
            vlcex=2.5 ,
            vlabels = pathway.names, calcex = 2
)

# Add a legend
legend.names <- c("RAD18-/- PD20\nTMZ vs DMSO", "RAD18+/+ PD20\nTMZ vs DMSO")
legend(x=1.1, y=1.1, legend = legend.names, bty = "n", pch=20 , col=colors_border ,
       text.col = "black", cex=2, pt.cex=3, y.intersp=1.8)
title("-log10 p-value", adj = 0.5, cex.main = 3)
########
dev.off()


