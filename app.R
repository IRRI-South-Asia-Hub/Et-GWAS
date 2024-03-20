library(shiny)
library(shinyjs)
library(shinythemes)
library(crayon)
library(ggplot2)
library(ggridges)
library(dplyr)
library(ggfortify)
library(ggpubr)
library(extrafont)
library(data.table)
library(CMplot)
library(foreach)
library(doParallel)
library(MASS)
library(zip)
library(tidyr)
library(plyr)
library(tidyverse)
library(haplotypes)
library(agricolae)
library(readxl)
library(RColorBrewer)


#setwd("/srv/shiny-server/pooled-GWAS/")
source("scirpts/theme_Publication.R")
source("scirpts/association_main.R")
#source("scirpts/haplo_pheno.R")

#Function to extract the genomic data
extract_irisID <- function(trait, infile, perc){
  
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
  phenofile <- infile
  # Step1: Phenotypic extract -----------------------------------------------
  
  outlist <- paste0(trait,"_list.txt")
  hist_plot <- paste0(trait,"_histogram.png")
  box_plot <- paste0(trait,"_boxplot.png")
  
  phe <- read.delim(infile, header = T)
  colnames(phe) <- c("Accessions","trait")
  
  box <- ggplot(phe) + aes(y = trait) +
    geom_boxplot(fill = "#0c4c8a") + theme_cus + xlab("") + ylab(trait)
  ggsave(box, filename = paste0("www/",box_plot),width = 10,height = 6)
  
  his <- ggplot(data = phe, aes(x = trait)) + geom_histogram()+
    theme_cus + xlab("")
  
  ggsave(his,filename = paste0("www/",hist_plot),width = 10,height = 6)
  
  out1 <- phe[,c(1,1)]
  write.table(out1, file = outlist, col.names = F, row.names = F,
              quote = F)
  
  
  # Step 2: genotypic data extract------------------------------------------------
  infam <- paste0(trait,"_geno")
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile data/server2 --keep ",
                trait,"_list.txt --noweb --make-bed --out ",infam))
  
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile ",
                infam," --noweb --recode --out ",infam))
  
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
  
  
  # pca_plink ---------------------------------------------------------------
  pheno <- "data/IRIS_pop_all.txt"
  theme_cus <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    strip.background=element_blank(),
                    axis.text.y=element_text(colour="black", face = "bold"),
                    axis.text.x=element_text(colour="black", face = "bold",angle = 90),
                    axis.ticks=element_line(colour="black"),
                    strip.text.x = element_text(colour="black", face = "bold"),
                    title = element_text(color = "black", size = 12, family = "Helvitica",
                                         face = "bold"),
                    plot.margin=unit(c(1,1,1,1),"line"),
                    legend.title = element_text(color = "black", size = 12, family = "Helvitica",
                                                face = "bold"),
                    legend.text = element_text(colour = "black", size = 12,family = "Helvitica"),
                    legend.key = element_rect(fill = "white", color = NA))
  
  
  infam2 <- paste0(trait,"_geno2")
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile data/server2 --keep ",
                trait,"_list.txt --noweb --make-bed --out ",infam2))
  
  system(paste0("softwares_external/plink2 --bfile ",infam2,
                " --pca 5 --nonfounders --out ",infam,"-pc"))
  
  pca_plot <- paste0(trait,"_pca.png")
  pov_plot <- paste0(trait, "_pov.png")
  
  
  #PLINK
  pca_file <- paste0(infam, "-pc.eigenvec")
  val_file <- paste0(infam, "-pc.eigenval")
  m <- read.table(pca_file, header = F)
  cat(blue(paste("Loading",pca_file,"...\n")))
  names(m)[1:5] <- c("FID","IID","PC1","PC2","PC3")
  
  phe <- read.delim(pheno, header = T)
  cat(blue(paste("Loading",pheno,"...\n")))
  names(phe) <- c("Name","IRIS.ID","SUBPOPULATION","COUNTRY")
  
  comb_p <- merge(m, phe, by.x = "FID", by.y = "IRIS.ID", all.x = T)
  
  p1 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC2, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                  "aus","japx","temp",
                                  "trop","subtrop","admix","aro" ),
                       values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352",
                                   "#00AFBB", "#000000", "#E69F00","#F0E442",
                                   "#D55E00", "#CC79A7","#999999", "#52854C"))
  p1 <- p1 + theme_Publication()
  
  ggsave(paste0(trait,"_pca.jpeg"), p1, width = 8, height = 8)
  
  p2 <- ggplot(comb_p) + geom_point(aes(x=PC1,y=PC3, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                  "aus","japx","temp",
                                  "trop","subtrop","admix","aro" ),
                       values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352",
                                   "#00AFBB", "#000000", "#E69F00","#F0E442",
                                   "#D55E00", "#CC79A7","#999999", "#52854C"))
  p3 <- ggplot(comb_p) + geom_point(aes(x=PC2,y=PC3, color=SUBPOPULATION)) +
    theme_cus + theme(legend.position = "none")+
    scale_color_manual(breaks = c("indx","ind2","ind1B","ind1A","ind3",
                                  "aus","japx","temp",
                                  "trop","subtrop","admix","aro" ),
                       values = c( "#56B4E9", "#56B4E9","#0072B2", "#293352","#00AFBB",
                                   "#000000", "#E69F00","#F0E442",
                                   "#D55E00", "#CC79A7","#999999", "#52854C"))
  
  ar <- ggarrange(p1,ggarrange(p2,p3,labels = c("B", "C"),ncol = 2,nrow = 1),
                  ncol = 1,nrow = 2,labels = "A", 
                  common.legend = T, legend = "bottom")
  print(ar)
  cat(green(paste0("Saving PCA plot for first 4 components in ",pca_plot, "...\n")))
  ggsave(pca_plot, ar, width = 8, height = 6)
  
  value <- read.table(val_file, header = F)
  var <- (value/sum(value))*100
  var$V2 <- c(1:nrow(value))
  
  pov <- ggplot(var) + geom_line(aes(x = V2, y = V1), color="red") +
    geom_point(aes(x = V2, y = V1),shape = 8) + 
    theme_cus + xlab("Principal Components") + ylab("PoV")
  
  cat(green(paste0("Saving PCA PoV of first 4 components in pca_pov.pdf", 
                   pov_plot,"...\n")))
  ggsave(pov_plot, pov, width = 8, height = 6)
  
  #convert to ped/map and then to genotype file
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile ",
                infam2 ," --noweb --recode --out ",infam2))
  
  system(paste0("perl softwares_external/PopLDdecay/bin/mis/plink2genotype.pl -inPED ",
                infam2,".ped -inMAP ",infam2,".map -outGenotype ", infam2,".genotype"))
  
  #LD calculation
  system(paste0("softwares_external/PopLDdecay/bin/PopLDdecay -InGenotype ",
                infam2,".genotype -OutStat ", infam2,"_LDdecay"))
  
  # #LD plot
  system(paste0("perl softwares_external/PopLDdecay/bin/Plot_OnePop.pl -inFile ",
                infam2,"_LDdecay.stat.gz -output ",trait,"_ldplot"))
  
  # pheno_dist_xp -----------------------------------------------------------
  cat(green("Pheno_dist\n"))
  
  popfile <- "data/IRIS_pop_all.txt"
  
  Type <- c("indx","ind2","aus","ind1B","ind1A",
            "ind3","japx","temp","trop","subtrop","admix","aro")
  
  # Initial distribution ----------------------------------------------------
  b <- read.delim(phenofile, header = T)
  names(b) <- c("Designation","phenotype")
  pop <- read.delim(popfile, header = T)
  names(pop) <- c("Name","Designation","Subpopulation","COUNTRY")
  
  val <- ceiling((perc*nrow(b))/100)
  
  bb <- merge(b, pop, by = "Designation", all = F)
  b <- bb[!is.na(bb$phenotype),]
  
  b$Type <- b$Subpopulation
  b$Type <- gsub(pattern = "indx",replacement = "ind",b$Type)
  b$Type <- gsub(pattern = "ind2",replacement = "ind",b$Type)
  b$Type <- gsub(pattern = "ind1B",replacement = "ind",b$Type)
  b$Type <- gsub(pattern = "ind1A",replacement = "ind",b$Type)
  b$Type <- gsub(pattern = "ind3",replacement = "ind",b$Type)
  
  low <- max(head(b[order(b$phenotype, decreasing=F), ], val)[,2])
  high <- min(head(b[order(b$phenotype, decreasing=T), ], val)[,2])
  
  b$group <- ifelse(b$phenotype <= low, 1, ifelse(b$phenotype >= high, 3, 2))
  supp.labs <- c("Low bulk", "Mixed","High bulk")
  names(supp.labs) <- c("1", "2", "3")
  
  image_file <- paste0(trait,"_",perc,"_dist1.png")
  #png(image_file, width = 1000, height = 650)
  dist <- ggplot(data=b, aes(phenotype)) + theme_bw() + 
    theme(legend.position="none",
          panel.background = element_blank(),panel.border=element_rect(fill=NA),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          panel.spacing = unit(0.3, "lines"),
          strip.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    geom_bar(data=subset(b,group==1),col="black", fill="red", 
             alpha = 0.5)+
    geom_bar(data=subset(b,group==2),col="black", fill="green", 
             alpha = 0.5) +
    geom_bar(data=subset(b,group==3),col="black", fill="blue", 
             alpha = 0.5, )+
    facet_grid(. ~ group, scales = "free", labeller = labeller(group = supp.labs))
  #dev.off()
  ggsave(image_file, dist, width = 8, height = 6)
  
  image_file <- paste0(trait,"_",perc,"_subdist1.png")
  #png(image_file, width = 1000, height = 650)
  dist2 <- ggplot(b, aes(x = phenotype, y = Subpopulation, fill=Subpopulation)) + 
    geom_density_ridges2(alpha=0.5)+
    xlab("Yield") +
    theme(legend.position="none",
          panel.background = element_blank(),panel.border=element_rect(fill=NA),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          panel.spacing = unit(0.3, "lines"),
          strip.text.x = element_text(size = 12),
          axis.title.y = element_blank(),
          axis.ticks = element_line(colour = "black"),
          axis.text.x=element_text(colour="black", size = 12),
          axis.text.y=element_text(colour="black", size = 12)) +
    geom_vline(xintercept = c(low,high),linetype="dashed", color = "red")
  #dev.off()
  ggsave(image_file, dist2, width = 8, height = 6)
  
  pheno_out <- paste0(trait,"_",perc,"_XPpheno.txt")
  write.table(b, file = pheno_out, sep = "\t", row.names = F, 
              col.names = T, quote = F)
  
  #Extra
  low <- b[b$group==1,c(1,1)]
  write.table(low, file = paste0(trait,"_low",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  high <- b[b$group==3,c(1,1)]
  write.table(high, file = paste0(trait,"_high",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  rand <- b[b$group==2,c(1,1)]
  rand2 <- rand[as.integer(runif(min(nrow(low),nrow(high)),1,nrow(rand))),]
  write.table(rand2, file = paste0(trait,"_rand",perc,".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  cat(green("Completed the pheno analysis\n"))
  
  
  # Pooling genotype --------------------------------------------------------
  infam <- paste0(trait,"_geno")
  hmpfile <- paste0(infam,".hmp.txt")
  
  low_in <- paste0(trait,"_low",perc,".txt")
  rand_in <- paste0(trait,"_rand",perc,".txt")
  high_in <- paste0(trait,"_high",perc,".txt")
  
  low_out <- paste0(trait,"_low",perc)
  rand_out <- paste0(trait,"_rand",perc)
  high_out <- paste0(trait,"_high",perc)
  
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile ",
                infam," --noweb --keep ",low_in," --recode --out ", low_out))
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile ",
                infam," --noweb --keep ",rand_in," --recode --out ",rand_out))
  system(paste0("softwares_external/plink-1.07-x86_64/./plink --bfile ",
                infam," --noweb --keep ",high_in," --recode --out ",high_out))
  
  system(paste0("bash scirpts/fileconversion.sh -h ",hmpfile," -p ",
                infam," -t ped2hap"))
  
  system(paste0("bash scirpts/pooling_snp.sh -h ",infam))
  
  bulks <- c("high","low","rand")
  
  
  start_time <- Sys.time()
  pooling_snp <- function(i){
    
    infile <- paste0(trait,"_",bulks[i],perc,".ped")
    c <- read.delim(infile, header = F, sep = " ",colClasses = c("character"))
    d <- c[,c(7:ncol(c))]
    out <- apply(d, 2, function(x){
      run = rle(sort(x[x!=0]))
      if(length(run$lengths) < 2){
        run$lengths[2] <- 0
        run$values[2] <- 0
      } 
      g <- stack(run)[1]
    })
    
    lstData <- Map(as.data.frame, out)
    dfrData <- rbindlist(lstData)
    
    count_cols <- sort(c(seq(1,nrow(dfrData),4),seq(2,nrow(dfrData),4)))
    all_cols <- sort(c(seq(3,nrow(dfrData),4),seq(4,nrow(dfrData),4)))
    df <- dfrData[count_cols]
    df$all <- dfrData[all_cols]
    
    colhead <- data.frame(rep(seq(1,(nrow(df)/4)), each=4))
    colnames(colhead) <- "num"
    
    comb <- cbind(colhead,df)
    comb$values <- as.numeric(comb$values)
    comb$alle <- paste(comb$num, comb$all, sep = "_")
    comb <- comb[comb$all != 0,]
    out <- comb %>% dplyr::group_by(alle) %>% 
      dplyr::summarize(tot = sum(values), nums = unique(num), .groups = 'drop')
    
    outfile <- paste0("temp_",bulks[i],".csv")
    write.csv(out, outfile, row.names = F)
  }
  cl<- detectCores()
  registerDoParallel(cl)
  result <- foreach (i=1:3) %dopar% {
    pooling_snp(i)
  }
  stopImplicitCluster()
  
  out_low <- read.csv("temp_low.csv")
  out_high <- read.csv("temp_high.csv")
  out_rand <- read.csv("temp_rand.csv")
  
  temp <- merge(out_high,out_low, by = "alle", all = T)
  comb <- merge(temp,out_rand, by = "alle", all = T)
  out <- comb[order(comb$nums.x),c(1,2,4,6,7)]
  colnames(out) <- c("RefAlt", "high" ,"low","random","num")
  end_time <- Sys.time()
  print(paste0("for pooling: ", (end_time - start_time)))
  
  # Ref and Alt alleles From HapMap -----------------------------------------
  hmpfile <- paste0(infam,"_hmp_allele.txt")
  hmp <- read.delim(hmpfile,header = T,sep = "\t")
  hmp$num <- c(1:nrow(hmp))
  hmp$Ref1 <- paste(hmp$num,hmp$Ref,sep = "_")
  hmp$Alt1 <- paste(hmp$num,hmp$Alt,sep = "_")
  ref <- hmp[,c(3,4)]
  alt <- hmp[,c(3,5)]
  
  ref_comb <- merge(ref,out,by.x = "Ref1", by.y = "RefAlt")
  alt_comb <- merge(alt,out,by.x = "Alt1", by.y = "RefAlt")
  names(alt_comb)[3:5] <- c("high_alt","low_alt","random_alt")
  names(ref_comb)[3:5] <- c("high_ref","low_ref","random_ref")
  
  comb <- merge(ref_comb,alt_comb, by = "num.x", all = T)
  final <- comb[,c(2,7,3,8,4,9,5,10)]
  colnames(final) <- c("Ref","Alt","high_ref","high_alt","low_ref",
                       "low_alt","random_ref","random_alt")
  
  final_out <- paste0(trait,"_pgwas_input",perc,".txt")
  write.table(final, file = final_out, sep = "\t",col.names = T, row.names = F,
              quote = F)
  print("the input file has been created")
  
  # Association -------------------------------------------------------------
  mapfile <- paste0(trait,"_geno.map")
  
  outfile <- paste(trait,"_pgwas_out",perc,".txt",sep = "")
  outqtls <- paste(trait,"_qtls",perc,".txt",sep = "")
  outsnps <- paste(trait,"_snp",perc,".txt",sep = "")
  
  input <- final
  input[is.na(input)]=0
  
  cols = c(3:8)    
  input[,cols] = apply(input[,cols], 2, function(x) as.numeric(as.character(x)));
  map <- read.delim(mapfile, header = F)
  input2 <- cbind(map[,c(2,1,4)],input[,c(3:8)])
  input <- input2
  names(input)[1] <- "snpid"
  names(input)[2] <- "chr"
  names(input)[3] <- "pos"
  
  xpgwas_modified <- function(input, filter=50,  plotlambda=TRUE){
    statout <- get_chistat(snps=input, filter=filter)
    
    lam <- estlambda(statout$stat, plot = plotlambda, proportion = 1, 
                     method = "regression",
                     filter = TRUE, main="Before genomic control")
    
    outqval <- xpgwas_qval(stat=statout, lambda=lam[['estimate']])
    return(outqval)  
  }
  
  qval <- xpgwas_modified(input, filter=50,  plotlambda=TRUE)
  
  b <- qval[,c(1:3,5)]
  names(b) <- c("SNP","Chromosome","Position","trait1")
  
  CMplot(b,plot.type="m",threshold=c(0.01,0.05)/nrow(b),
         threshold.col=c('red','orange'),
         multracks=FALSE, chr.den.col=NULL, file.output=FALSE)
  
  write.table(qval, file = outfile, col.names = T, row.names = F, 
              quote = F,sep = "\t")
  
  
  #xpplot_mine
  
  # Findign significant snps ------------------------------------------------
  suggestiveline = -log10(5e-05)
  genomewideline = -log10(5e-08)
  
  qval$log10p <- -log10(qval$pval)
  SNPset <- qval
  qtltable <- SNPset[SNPset$log10p >= suggestiveline,]
  write.table(qtltable$snpid, file = outsnps, col.names = F, row.names = F, 
              quote = F, sep = "\t")
  
  
  qtltable <- SNPset %>% dplyr::mutate(passThresh = log10p >= suggestiveline) %>%
    dplyr::group_by(chr, run = {
      run = rle(passThresh)
      rep(seq_along(run$lengths), run$lengths)
    }) %>% dplyr::filter(passThresh == TRUE) %>% dplyr::ungroup()
  
  if(nrow(qtltable) > 0){
    qtltable <- qtltable %>% dplyr::ungroup() %>% dplyr::group_by(chr) %>%
      dplyr::group_by(chr, qtl = {
        qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
      }) %>% dplyr::select(-run) %>%
      dplyr::summarize(start = min(pos), end = max(pos), length = end - start + 1,
                       nSNPs = length(pos),
                       avgSNPs_Mb = round(length(pos)/(max(pos) - min(pos)) * 1e+06),
                       .groups = 'drop')
  }
  names(qtltable)[1] <- "CHR"
  write.table(qtltable, file = outqtls, col.names = T, row.names = F, 
              quote = F, sep = "\t")
  print("everything is done")
  
  list(plot1 = box, plot2 = his,plot3 = ar, plot4 = pov, plot5 = dist,
       plot6 = dist2)
}

haplo_pheno <- function(can_genes, pheno, path){
  
  filenames <- read.csv(can_genes, header = F)
  colnames(filenames) <- "V1"
  #list.files(path = path, pattern = ".csv")
  
  #filenames <- gsub(".csv","",filenames)
  print(head(filenames))
  
  pheno <- read_xlsx(pheno)
  names(pheno) <- c("Designation","trait")
  print(paste0("size:",nrow(filenames)))
  
  for (n  in c(1:nrow(filenames))) {
    skip_to_next <<- FALSE
    i <- filenames$V1[n]
    print(i)
    a <- read.csv(paste0(path,"/",i,".csv"), header = TRUE, stringsAsFactors = FALSE)
    
    dx <- paste0("Now running: ",i)
    if(ncol(a) > 6) {
      a <- data.frame(a [-c(2:2),-c(2:5)])
      a[is.na(a)] <- "-"
      b <- a %>% filter(across(everything(a), ~ !grepl("/",.)))
      #extract----
      lst <- data.frame(pheno$Designation)
      bb <- merge(b, lst, by.x = "JAPONICA.NIPPONBARE.POSITIONS", by.y = "pheno.Designation")
      b <- bb
      #------------
      c <- gsub(' ', '_', b$JAPONICA.NIPPONBARE.POSITIONS)
      c <- cbind(c,b[,-c(1:1)])
      c <-data.frame(c)
      colnames(c)<-colnames(b)
      d<- unite(c[,-c(1:1)], col = "seq", remove =  TRUE, sep = "", na.rm = TRUE)
      e <- data.frame(cbind(c$JAPONICA.NIPPONBARE.POSITIONS,d$seq))
      e<-e[!(str_count(e$X2) < 1),]
      
      dir.create(paste0(path,"/",i,""))
      sink(file.path((paste0(path,"/",i,"")), "Sequence.fas"))
      for (j in 1:nrow(e)){
        name = paste0(">",e[j,1])
        sequence = paste0(e[j,2])
        cat(name,sep="\n")
        cat(sequence,sep="\n")
      }
      sink()
      x<- read.fas(file.path((paste0(path,"/",i,"")), "Sequence.fas"))
      haplo<- haplotype(x)
      haploy<- data.frame(haplo@sequence)
      col <- dimnames(c [-c(1:1)])
      colnames(haploy)<- col[[2]]
      row<- paste("H",seq(1:nrow(haploy)) , sep = "")
      row.names(haploy)<- row
      haploy$counts <- haplo@freq
      haploy$freq <- haplo@freq/nrow(e)
      popul <- paste("pop", seq(1:haplo@nhap), sep = "")
      par<-grouping(haplo, popul)
      groups <- paste0("H",par[["hapvec"]])
      variety<- data.frame(e [,1:1])
      variety$group <- groups
      var_pheno <- variety
      colnames(var_pheno) <- c("Accession", "group")
      genotype <- data.frame(rev((haploy))[-c(1:2)])
      genotype <- data.frame(rev((genotype)))
      variety<- merge(variety,c, by.x = 'e...1.1.',by.y = 'JAPONICA.NIPPONBARE.POSITIONS')
      write.csv(haploy,file.path((paste0(path,"/",i,"")),file = "Haplotype.csv"), row.names = F)
      write.csv(variety,file.path((paste0(path,"/",i,"")), file = "Variety.csv"), row.names = F)
      write.table(genotype,file.path((paste0(path,"/",i,"")),file = "Flapjack.genotype"),sep="\t", quote = F)
      var <- gsub('_', ' ', var_pheno$Accession)
      var_pheno$names <- var
      
      haplopheno <- merge(var_pheno, pheno, by.x = 'names', by.y = 'Designation')
      haplopheno$group<- as.factor(haplopheno$group)
      haplopheno<-haplopheno[order(haplopheno$group),]
      hp <- count(haplopheno$group)
      hp <- hp[!(hp$freq <2),]
      sd<- merge(hp,haplopheno, by.x = 'x', by.y = 'group')
      sd<- sd  %>% dplyr::select(1, 5)
      
      colu<- colnames(sd)
      colnames(sd)<- c("Haplotype",colu[2])
      colu <- colnames(sd)
      sink(file.path((paste0(path,"/",i,"")), "anova.txt"))
      aooov<- aov(formula((paste0(colnames(sd)[2], "~", colnames(sd)[1]))), data = sd)
      print(summary(aooov))
      sink()
      sink(file.path((paste0(path,"/",i,"")), "duncan.txt"))
      tryCatch(print(duncan<- duncan.test(aooov,trt = "Haplotype",console = TRUE, group = T)), 
               error = function(e) { skip_to_next <<- TRUE})
      if(skip_to_next) { 
        cat(red(paste0("There is error in: ",i,". so skipping to next")))
        sink()
      }
      if(skip_to_next) { 
        next 
      }     
      #print(duncan<- duncan.test(aooov,trt = "Haplotype",console = TRUE, group = T))
      sink()
      sink(file.path((paste0(path,"/",i,"")), "duncan_pval.txt"))
      print(duncan.pval<- duncan.test(aooov,trt = "Haplotype",console = TRUE, group = F))
      sink()
      summary<- data.frame(duncan$means)
      summary<- cbind(hp,summary)
      
      write.csv(summary,file.path((paste0(path,"/",i,"")),file = "summary.stat.csv"))
      labels <- data.frame(duncan$groups)
      labels$Haplotype<- row.names(labels)
      labels<- labels[order(labels$Haplotype),]
      specify<- cbind(hp,labels)
      
      specify2 <- specify
      specify2$final <- paste0(specify2$x,"-", round(specify2[,3], digits=2),"(", specify2$groups,")")
      write.csv(specify2,file.path((paste0(path,"/",i,"")),file =  "group.csv"))
      
      mycolors = c(brewer.pal(name="BuGn", n = 9), brewer.pal(name="OrRd", n = 9),brewer.pal(name = "YlOrBr", n=9),brewer.pal(name = "Pastel1",n=9))
      jpeg(file.path((paste0(path,"/",i,"")), "boxplot.jpg"), height = 1200,width = 1400, res = 300)
      p<-ggplot(sd, aes(x=Haplotype, y=sd[,2])) +labs(y=paste0(colnames(sd)[2]))+ geom_boxplot(aes(fill = factor(Haplotype)), show.legend = F)+scale_color_manual(values = mycolors, aesthetics = "fill")+ theme(axis.text.x = element_text(angle = 90, hjust= 1.0))
      q<-p+geom_text(data = specify, aes(x,Inf, label= paste0("n=",freq,"\n",groups)), vjust = 1,size= 1.5)
      r<- q+theme_bw()+theme(text = element_text(size=10, face="bold"),axis.text.x = element_text(angle = 90, hjust= 1.0))
      s<- r + geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5,color = "red")
      t<- s+ggtitle(label = paste0(i,""))
      print(s)
      dev.off()
      var_300 <- gsub('_', ' ', variety[,1])
      variety_300 <- variety
      variety_300$names <- var_300
      
      variety_300<- merge(variety_300, pheno, by.x = 'names', by.y = 'Designation')
      variety_300<- variety_300[order (variety_300$group),]
      write.csv(variety_300,file.path((paste0(path,"/",i,"")),file = "variety_300.csv"))
      count <- count(haplopheno$group)
      colnames(count)<- c("x", "count_300")
      count$freq_3k <- count$count_300/nrow(e)
      count$freq_300 <- count$count_300/nrow(haplopheno)
      hap<- cbind(row.names(haploy), haploy)
      haplotype_300<- merge(hap, count, by.x = 'row.names(haploy)', by.y = 'x' )
      write.csv(haplotype_300,file.path((paste0(path,"/",i,"")),file = "haplotype_300.csv"))
    }
  }
  
  return(dx)
}

ui = fluidPage(
  tabsetPanel(
    tabPanel("Map", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(selectInput(inputId = "Bulk_size",
                                        label = "Bulk size:",
                                        choices = c("10", "15", "20")),
                            
                            textInput("Trait", "Trait" , "Ex: SPY "),
                            verbatimTextOutput("value"),
                            
                            fileInput("upload", "Phenotype"),
                            
                            actionButton("runButton", "Run"),
                            
                            downloadButton("downloadData", label = "Download MTAs"),
                            conditionalPanel(
                              condition = "input.runButton > 0",
                              p("~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"),
                              p("Association analysis started"),
                              p("~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"),
                              verbatimTextOutput("progress_output")
                            )
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Phenotype", plotOutput("plot1", width = "100%", 
                                                    height = "300px"),
                            plotOutput("plot2", width = "100%", height = "300px")),
                   tabPanel("Population structure", plotOutput("plot3", width = "100%", 
                                                               height = "300px"),
                            plotOutput("plot4", width = "100%", height = "300px")),
                   tabPanel("Phenotypic distribution", plotOutput("plot5", width = "100%", 
                                                                  height = "300px"),
                            plotOutput("plot6", width = "100%", height = "300px")),
                   tabPanel("Association", plotOutput("plot7", width = "100%", 
                                                      height = "300px"))),
                 h6("Sample download", align = "center"),
                 verbatimTextOutput("result_output")
               )
             )
    ),
    tabPanel("Haplo-Pheno Analysis", fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 fileInput("upload2", "Phenotype"),
                 fileInput("upload3", "Candidates"),
                 textInput("input_path", "Enter a Path:", value = ""),
                 actionButton("runHapPheno", "Submit"),
                 
                 downloadButton("downloadHaps", label = "Download Haps"),
               ),
               mainPanel(fluidRow(
                 textOutput("Hap_progress"),
                 textOutput("Hap_final")
               )
               )
             )
    )
  )
)

server <- function(input, output, session) {
  
  plotData <- reactiveVal(NULL)
  processingResult <- reactiveVal(NULL)
  
  observeEvent(input$runButton, {
    
    processingResult(NULL)  
    system.time({
      plotData(extract_irisID(input$Trait,input$upload$datapath, as.numeric(input$Bulk_size)))
    })
    processingResult()
  })
  
  output$resultText <- renderPrint({
    processingResult()
  })
  
  output$plot1 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot1
  })
  output$plot2 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot2
  })
  output$plot3 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot3
  })
  output$plot4 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot4
  })
  output$plot5 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot5
  })
  output$plot6 <- renderPlot({
    if (!is.null(plotData()))
      plotData()$plot6
  })
  output$plot7 <- renderPlot({
    
    qval <- read.delim(paste(input$Trait,"_pgwas_out",as.numeric(input$Bulk_size),".txt",sep = ""))
    b <- qval[,c(1:3,5)]
    names(b) <- c("SNP","Chromosome","Position","trait1")
    
    CMplot(b,plot.type="m",threshold=c(1e-4,1e-6),
           threshold.col=c('red','orange'),
           multracks=FALSE, chr.den.col=NULL, file.output=FALSE)
  })
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("output", "zip", sep=".")
    },
    content <- function(file) {
      dir.create("out")
      system(paste0("cp ",trait,"*.png ",
                    trait,"*.jpeg ",
                    trait,"_pgwas* ",
                    trait,"*.pdf ",
                    "out"))
      files2zip <- dir('out', full.names = TRUE)
      zip(zipfile = 'testZip', files = files2zip)
      file.copy("testZip", file)
    },
    contentType = "application/zip"
  )
  
  #for the second panel
  hapData <- reactiveVal(NULL)
  running <- reactiveVal(FALSE)
  
  observeEvent(input$runHapPheno, {
    running(TRUE)
    hapData(haplo_pheno(input$upload3$datapath,input$upload2$datapath,
                        input$input_path))
  })
  running(FALSE)
  
  observe({
    if (running()) {
      output$Hap_final <- renderText("Haplo Pheno Analysis done. You can stop the App")
    } else {
      output$Hap_final <- renderText("Analysis is still running...")
    }
  })

}

shinyApp(ui, server)
