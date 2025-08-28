##############################################################################
# Analysis of SOLA (2009-2014) time series: KEGG and gene tables
##############################################################################

# packages
library(readr)
library(data.table)
library(stringr)
library(lomb)
library(ggplot2)
library(hydroTSM)
library("pgirmess")
library(tidyr)

# 1/ Import, format, and explore table
# -------------------------------------

# KEGG table
keggs <- read_tsv("~/Desktop/PETRIMED/DATA/METAG/KEGG_SOLA_2009_2014/BBMOSOLA-common-dates-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl.hierarchy.txt",quote = "")
keggs <- keggs[,colSums(is.na(keggs))<nrow(keggs)]
keggs<-keggs[,-grep("BL",colnames(keggs))]
keggs<-as.data.frame(keggs[rowSums(keggs[, grep("SO", colnames(keggs))])>0,])
rownames(keggs)<-keggs$annot
keggs.ab<-keggs[, grep("SO", colnames(keggs))]

# Functional Annotation of KEGGs
keggs.funct<-keggs[,grep("SO|..162", colnames(keggs), invert = TRUE)]
functs <- read_tsv("~/Desktop/PETRIMED/DATA/METAG/KEGG_SOLA_2009_2014/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl",quote = "")
functs<-as.data.frame(cbind(functs$annot,str_split_fixed(functs$`Gene name; KEGG annotation`, pattern = "; ", n = 2)))
functs<-as.data.frame(cbind(functs,str_split_fixed(functs$V2, pattern = ", ", n=9)))
row.names(functs)<-functs$V1
functs<-merge(functs, keggs.funct[,-2], by.x = "V1", by.y = "annot")
colnames(functs)<-c("kegg", "Gene names", "kegg annotation", paste0("GN", 1:9), paste0("KH",1:31 ) )
functs<-merge(functs,data.frame(tot.ab = rowSums(keggs.ab)), by.x = "kegg", by.y = "row.names")

# Info Samples
info<-data.frame(ID = colnames(keggs.ab),Date = as.Date(gsub("SO","",colnames(keggs.ab)), format = "%y%m%d") )
info$DATE<-as.POSIXct(info$Date)
info$nDate<-as.numeric(info$Date)
info$Year<-format(info$Date, "%Y")
info$Month<-format(info$Date, "%b")
info$Season<-factor(str_to_title(time2season(info$Date, out.fmt = "seasons")), levels = c("Winter", "Spring", "Summer", "Autumn"), ordered = TRUE)

# Gene table
genes <- fread("~/Desktop/PETRIMED/DATA/METAG/KEGG_SOLA_2009_2014/SOLA_2009_2014_gene_annot_filtered.tbl", fill = 0,)
genes.annot<-as.data.frame(genes[,c(1,2)])
genes<-as.data.frame(genes[,-c(1,2)])
rownames(genes)<-genes.annot$gene
genes.annot <- as.data.frame(separate_rows(genes.annot, annotation, sep = ","))

# info gene info to KEGG table
functs<-merge(functs,as.data.frame(table(genes.annot$annotation)) ,by.x = "kegg", by.y = "Var1", all.x = TRUE);colnames(functs)[ncol(functs)]<-"nb.genes"
functs[is.na(functs$nb.genes), "nb.genes"]<-0

# 2/ Analysis of KEGG rhythmicity
# -------------------------------------

## a. Analysis of rhythmicity
# prepare data
lsp.keggs<-merge(info, t(keggs.ab), by.x = "ID", by.y = "row.names" )

# loop the command over all Keggs
#r=NULL
#for (i in 9:ncol(lsp.keggs)){
#  rl<-summary(randlsp(x = lsp.keggs[, c(4, i)], repeats = 99, type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
#  r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[i], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[13,]) ))
#  #rl<-summary(lsp(x = lsp.keggs[, c(4, i)],type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
#  #r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[i], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[12,]) ))
#  print(paste(i, ncol(lsp.keggs), sep = " / "))
#}
#saveRDS(r, "Desktop/PETRIMED/ANALYSES/RHYTHM/lsp_keggs_r99.RData")
r<-readRDS("Desktop/PETRIMED/ANALYSES/RHYTHM/lsp_keggs_r99.RData")
r<-merge(r,functs, by ="kegg" )

# Filter rhytmic KEGGs
r.keggs<-r[r$PNmax >= 0.1 & r$p.value < 0.01,]
nr.keggs<-r[!r$kegg %in% r.keggs$kegg ,]

## b. Analysis of specific KEGG rhythmicity

# Define periods to highlight
highlight_periods <- data.frame(
  xmin = seq(as.POSIXct("2009/1/1"), as.POSIXct("2015/1/1"), "years"),
  xmax = seq(as.POSIXct("2009/12/31"), as.POSIXct("2015/12/31"), "years"),
  ymin = -Inf,
  ymax = Inf,
  fill = c("A","B","A","B", "A", "B", "A") )

## explore patterns over the time-series: #1 most abundant KEGG
# identify most rhythmic KEGG
r.keggs[r.keggs$PNmax == max(r.keggs$PNmax),]
# test correlation between KEGG and Total gene abundance of the KEGG
gene.tab<-data.frame(gene.tab = colSums(genes[rownames(genes) %in% genes.annot[grep("K06287", genes.annot$annotation) , "gene" ],]))
gene.tab<-merge(gene.tab, lsp.keggs, by.x = "row.names", by.y = "ID")
plot(gene.tab$K06287, gene.tab$gene.tab) # not the same values, order, but overall correlation

# plot
functs[functs$kegg == "K06287",]
ggplot(lsp.keggs, aes(x = DATE, y = K06287))+
  #facet_grid(.~Year, space = "free", scale = "free")+
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y", expand = c(0,0))+
  geom_rect(data = highlight_periods, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_manual(values = c("#ccede2", "#ccd9ed"), guide = "none")+
  geom_line(linewidth = 0.9, col = "#184873")+
  theme_minimal()+
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 12))

# plot genes
m.genes<-reshape2::melt(as.matrix(genes[rownames(genes) %in% genes.annot[grep("K06287", genes.annot$annotation) , "gene" ],]))
m.genes<-merge(info, m.genes, by.x = "ID", by.y = "Var2")
#m.genes<-m.genes[m.genes$value>0,]

ggplot(m.genes, aes(x = DATE, y = value, fill = Var1))+
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y", expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_area(position = "stack", colour = "transparent", show.legend = FALSE)+
  #geom_rect(data = highlight_periods, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), alpha = 0.5, inherit.aes = FALSE) +
  #scale_fill_manual(values = c("#ccede2", "#ccd9ed"), guide = "none")+
  theme_minimal()+
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 12))

## c. Other properties of rhythmicity
# Are abundant KEGGs more rhythmic?
ggplot(r, aes(x =tot.ab, y =PNmax))+
  scale_x_log10()+
  geom_point(col = "#184873")+
  geom_smooth(method = "lm", colour = "coral3", se = FALSE)+
  theme_minimal()+
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.title  = element_text(size = 12))

# Are KEGG 1st level hierarchy categories rhythmically distinct?
kruskal.test(r$PNmax, r$KH1)
kruskalmc(r$PNmax, r$KH1)
ggplot(r, aes(x =KH1, y =PNmax))+
  geom_boxplot(col = "#184873")+
  theme_minimal()+
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title  = element_text(size = 12))

# 3/ Analysis of genes rhythmicity
# -------------------------------------

## a. Analysis of rhythmicity
# prepare data
lsp.genes<-merge(info, t(genes), by.x = "ID", by.y = "row.names" )
colnames(lsp.genes)
# loop the command over all Keggs
r=NULL
for (i in 8:ncol(lsp.genes)){
  #rl<-summary(randlsp(x = lsp.genes[, c(4, i)], repeats = 99, type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
  #r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[i], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[13,]) ))
  rl<-summary(lsp(x = lsp.genes[, c(4, i)],type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
  r<-rbind(r,data.frame(kegg = colnames(lsp.genes)[i], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[12,]) ))
  print(paste(i, ncol(lsp.genes), sep = " / "))
}
saveRDS(r, "Desktop/PETRIMED/ANALYSES/RHYTHM/lsp_genes_base.RData")
r<-readRDS("Desktop/PETRIMED/ANALYSES/RHYTHM/lsp_genes_base.RData")


