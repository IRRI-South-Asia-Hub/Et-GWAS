library(ggplot2)
library(ggpubr)
library(forcats)
library(ggthemes)
library(tidyr)
library(crayon)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  cat(green("Give the following files \n"))
  cat(green("1.Trait name (this will be in the plot names) \n"))
  cat(green("2.CV file (K values and corresponding cv values) \n"))
  cat(green("3.Pheno final file (IRIS ID and phenotype) \n"))
  cat(red("*************************************************************\n"))
  cat(blue("Rscript --vanilla <traitname> <cv_file> \n"))
  cat(red("*************************************************************\n"))
  stop("one or many missing files \n", call.=FALSE)
}

trait <- args[1]
cv_file <- args[2]
phe_file <- args[3]

pop_all <- read.delim("/home/niranjani/Documents/data/IRIS_pop_all.txt", header = T)
phe_list <- read.delim(phe_file, header = F, sep = "")
colnames(pop_all) <- c("Name","IRIS.ID","Subpopulation","COUNTRY")
colnames(phe_list) <- c("IRIS.ID","phe")

pop <- merge(phe_list,pop_all, by = "IRIS.ID", all.x = T,sort = F)

cluster_png <- paste0(trait,"_cluster.png")

theme_cus <-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),
                  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                  strip.background=element_blank(),
                  axis.text.x=element_text(colour="black"),
                  axis.text.y=element_text(colour="black"),
                  axis.ticks=element_line(colour="black"),
                  plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none")

# cv_error ----------------------------------------------------------------
cv <- read.delim(cv_file, header = F)
names(cv) <- c("K","CV")
ggplot(cv) + geom_line(aes(x = K, y = CV)) + 
  geom_point(aes(x = K, y = CV),shape = 8, colour = "red", size = 2)+
  scale_x_continuous(breaks = cv$K) + ylab("Cross-Validation error")+
  theme(panel.background = element_rect(colour = "black"),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black", size = 12),
        axis.text.y=element_text(colour="black", size = 12),
        axis.ticks=element_line(colour="black"),
        axis.title.y = element_text(colour="black", size = 12),
        plot.margin=unit(c(1,1,1,1),"line"),legend.position = "none")


# cluster_plot ------------------------------------------------------------
gvec<-vector("list",length=nrow(cv))
gvec<-vector("list",length=3)

for (i in c(1:nrow(cv))) {
  filename=paste0(trait,".",i,".Q")
  tbl=read.table(filename, header = F)

  tbl$sub <- pop$Subpopulation
  tbl <- tbl %>% mutate(id = rownames(tbl))
  df <- gather(tbl, key, value, -id, -sub)
  df <- df[!is.na(df$sub),]
  df <- df[order(df$sub),]
  kplot <-
    ggplot(df, aes(x = id, y = value, fill = factor(key))) +
    geom_col(color = "gray", size = 0.1) +
    theme_minimal() + labs(x ="", y = "Ancestry") +
    facet_grid(~fct_inorder(sub), switch = "x", scales = "free", space = "free") + 
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    theme(
      panel.spacing.x = unit(0.001, "lines"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(colour = "black", size = 12),
      axis.text.x = element_blank(),
      strip.text.x = element_blank()) +
    scale_fill_gdocs(guide = FALSE)
  if(i==nrow(cv)){
    kplot <-
      ggplot(df, aes(x = id, y = value, fill = factor(key))) +
      geom_col(color = "gray", size = 0.1) +
      theme_minimal() + labs(x ="", y = "Ancestry") +
      facet_grid(~fct_inorder(sub), switch = "x", scales = "free", space = "free") + 
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = expand_scale(add = 1)) +
      theme(
        panel.spacing.x = unit(0.001, "lines"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(colour = "black", size = 16, family = "Arial", face = "bold"),
        axis.text.x = element_blank(),
        strip.text.x = element_text(color = "black", size = 16, angle = 90, 
                                    family = "Arial", face = "bold", hjust = 1)) +
      scale_fill_gdocs(guide = FALSE)
    
  }
  kplot
  gvec[[i]]<-kplot
}
out[[3]] <- gvec[[3]]

ggsave("kcluster.jpeg", kplot, width = 25, height = 5, dpi = 300)

png(cluster_png,width = 2000, height = 1314)
ggarrange(plotlist = list(gvec[[1]]:gvec[[8]]), ncol = 1, nrow = 8)
dev.off()


# using igraph ------------------------------------------------------------
library(igraph)

plot(x, axes = FALSE, add = FALSE, xlim = c(-1, 1),
     ylim = c(-1, 1), mark.groups = list(), mark.shape = 1/2,
     mark.col = rainbow(length(mark.groups), alpha = 0.3),
     mark.border = rainbow(length(mark.groups), alpha = 1),
     mark.expand = 15)
