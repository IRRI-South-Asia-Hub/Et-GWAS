extract_irisID <- function(trait, infile){
  library(crayon)
  library(ggplot2)
  library(crayon)
  library(tidyr)
  library(ggpubr)
  
  theme_cus <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    strip.background=element_blank(),axis.text.x=element_text(colour="black"),
                    axis.text.y=element_text(colour="black"),
                    axis.ticks=element_line(colour="black"),
                    plot.margin=unit(c(1,1,1,1),"line"),
                    legend.title = element_text(color = "black", size = 10),
                    legend.text = element_text(color = "red"),
                    legend.background = element_rect(fill = "lightgray"),
                    legend.key = element_rect(fill = "white", color = NA))

# Step1: Phenotypic extract -----------------------------------------------

  outlist <- paste0(trait,"_list.txt")
  hist_plot <- paste0(trait,"_histogram.png")
  box_plot <- paste0(trait,"_boxplot.png")
  
  phe <- read.csv(infile, header = T, "\t")
  colnames(phe) <- c("Accessions","trait")
  
  p1 <- ggplot(phe) + aes(y = trait) +
    geom_boxplot(fill = "#0c4c8a") + theme_cus + xlab("") + ylab(trait)
  ggsave(p1, filename = paste0("www/",box_plot),width = 10,height = 6)

  p <- ggplot(data = phe, aes(x = trait)) + geom_histogram()+
     theme_cus + xlab("")
  
  ggsave(p,filename = paste0("www/",hist_plot),width = 10,height = 6)

  out1 <- phe[,c(1,1)]
  write.table(out1, file = outlist, col.names = F, row.names = F,
              quote = F)
  

# Step 2: genotypic data extract------------------------------------------------
  infam <- paste0(trait,"_geno")
  system(paste0("plink --bfile data/pruned_v2.1 --keep ",
                trait,"_list.txt --noweb --make-bed --out ",infam))
  
  a <- read.delim(paste0(infam,".fam"), header = F,sep = "")
  names(a)[1:6] <- c("Fam_ID","Ind_ID","Paternal","Maternal","Sex","Phenotype")
  b <- phe
  names(b) <- c("IRIS.ID","Phenotype")
  comb <- merge(a,b, by.x = "Ind_ID", by.y = "IRIS.ID", all = F)
  comb$Phenotype.x <- comb$Phenotype.y
  out <- comb[!is.na(comb$Fam_ID),]
  out <- out[,c(1:6)]
  out$Paternal <- 0
  out$Maternal <- 0
  out$Sex <- 0
  out$Phenotype.x[is.na(out$Phenotype.x)] <- -9
  
  write.table(out, file = paste0(infam,".fam"), col.names = F, row.names = F,
              quote = F, sep = " ")
  
}
