##############################################################################
# Analysis of SOLA (2009-2014) time series: KEGG and gene tables
##############################################################################

# packages
library("ggsci")
library(RColorBrewer)
library(readr)
library(data.table)
library(stringr)
library(lomb)
library(ggplot2)
library(hydroTSM)
#library("pgirmess")
library(tidyr)
library(seas)
library(dplyr)
library(lubridate)
library(ggpubr)



# 1/ Import, format, and explore KEGG table
# -----------------------------------------

# KEGG table
keggs <- read_tsv("~/Desktop/PETRIMED/DATA/METAG/KEGG_SOLA_2009_2014/BBMOSOLA-common-dates-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl.hierarchy.txt",quote = "")
#keggs <- read_tsv("BBMOSOLA-common-dates-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl.hierarchy.txt",quote = "")
keggs <- keggs[,colSums(is.na(keggs))<nrow(keggs)]
keggs<-keggs[,-grep("BL",colnames(keggs))]
keggs<-as.data.frame(keggs[rowSums(keggs[, grep("SO", colnames(keggs))])>0,])
rownames(keggs)<-keggs$annot
keggs.ab<-keggs[, grep("SO", colnames(keggs))]

# Functional Annotation of KEGGs
keggs.funct<-keggs[,grep("SO|..162", colnames(keggs), invert = TRUE)]
functs <- read_tsv("~/Desktop/PETRIMED/DATA/METAG/KEGG_SOLA_2009_2014/BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl",quote = "")
#functs <- read_tsv("BBMOSOLA-GC_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.gene_name-kegg_annotation.tbl",quote = "")
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
info$JD<-as.numeric(format(info$Date, "%j"))

# 2/ Create randomization list/matrix
# -----------------------------------------

## a) Create year combinations
# Initialize a list to store data.frames
Years=unique(info$Year)
nYears=length(unique(info$Year))
yl<-year.length(Years)

# Create all unique combinations of years from 2 years to 7
comb_list <- list()
for (k in 2:nYears) {
    tmp <- as.data.frame(t(combn(Years, k)))
    colnames(tmp) <- paste0("Year", 1:k)
    comb_list[[as.character(k)]] <- tmp
}

# Combine all into a single data.frame (with varying number of columns)
all_combs <- bind_rows(comb_list, .id = "size" )[,-1]
all_combs$NB <- apply(all_combs, 1, function(x) max(which(!is.na(x))))

## b) save info for each combinations and create a continuous Julian Days column for each timeseries
# Function to check leap year
is_leap <- function(year) {(year %% 4 == 0 & year %% 100 != 0) | (year %% 400 == 0)}
info <- info[order(info$Date), ]

# loop to save info with new Julian Days column
list_comb_info<-list()
for(j in seq_len(nrow(all_combs))){
  # Initialize cumulative offset
  info.tmp<-info[info$Year %in% all_combs[j,],]
  year_offset <- 0
  info.tmp$JDD <- NA
  prev_year <- NA
  
  # loop Julian days computation
  for (i in seq_len(nrow(info.tmp))) {
    this_year <- info.tmp$Year[i]
  
    # If new year, update offset
    if (!identical(this_year, prev_year)) {
      if (!is.na(prev_year)) {
        year_offset <- year_offset + ifelse(is_leap(as.numeric(prev_year)), 366, 365)}
      prev_year <- this_year
    }
    # Add cumulative offset to JD
    info.tmp$JDD[i] <- info.tmp$JD[i] + year_offset
  }
  # save it for each combination
  list_comb_info[[paste(j)]]<-info.tmp
}

# 3/ Compute lomb-scargle periodogram over all year combinations
# --------------------------------------------------------------

# list_comb_lsp<-list()
# for (cb in seq_len(length(list_comb_info))){
# 
#   # merge info on the time series with KEGGs abundances
#   lsp.keggs<-merge(list_comb_info[[cb]], t(keggs.ab), by.x = "ID", by.y = "row.names" )
# 
#   # loop the command over all Keggs
#   r=NULL
#   for (k in 11:ncol(lsp.keggs)){
#     rl<-summary(randlsp(x = lsp.keggs[, c(9, k)], repeats = 99, type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
#     r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[k], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[13,]), Phase = format(lsp.keggs[which.max(lsp.keggs[,k]),"Date"], "%b") ))
#     print(paste(cb, k, ncol(lsp.keggs), sep = " / "))
#   }
#   list_comb_lsp[[paste(cb)]]<-r
# }
# saveRDS(list_comb_lsp, "list_comb_lsp.RData")

# 4/ Analyze best combinations per KEGG
# --------------------------------------------------------------

# results
list_comb_lsp<-readRDS("~/Desktop/PETRIMED/ANALYSES/ts_rhythm/list_comb_lsp.RData")

## format results
# empty files
pnmax=NULL
period=NULL
pv=NULL
phase=NULL

# extract and merge results per coefficient
for (i in 1:120){
  pnmax<-cbind(pnmax,list_comb_lsp[[i]][,2]);colnames(pnmax)[ncol(pnmax)]<-i
  period<-cbind(period,list_comb_lsp[[i]][,3]);colnames(period)[ncol(period)]<-i
  pv<-cbind(pv,list_comb_lsp[[i]][,4]);colnames(pv)[ncol(pv)]<-i
  phase<-cbind(phase,list_comb_lsp[[i]][,5]);colnames(phase)[ncol(phase)]<-i
}

# format rownmaes of result table
rownames(pnmax)<-list_comb_lsp[[i]][,1]
rownames(period)<-list_comb_lsp[[i]][,1]
rownames(pv)<-list_comb_lsp[[i]][,1]
rownames(phase)<-list_comb_lsp[[i]][,1]

# which combination is the best for each KEGG in terms of PNmax
comb.PN<-data.frame(PNmax =  apply(pnmax,1, which.max))

# which Year's removal generate the highest PNmax?
missing_years<-list()
for (y in 1:120){missing_years[[y]]<-Years[!Years %in% list_comb_info[[y]]$Year]}
missing_years

all.KO=NULL
for (k in rownames(comb.PN) ){all.KO<-c(all.KO,missing_years[[comb.PN[k,]]])}
table(all.KO)

all.years=NULL
for (y in 1:120 ){all.years<-rbind(all.years,unlist(all_combs[y,-ncol(all_combs)],use.names=FALSE))}
table(all.years)

# filter only significant results
list_comb_lsp_filtered<-list()
for (cb in 1:120){list_comb_lsp_filtered[[cb]]<-list_comb_lsp[[cb]][list_comb_lsp[[cb]]$PNmax>0.1 &list_comb_lsp[[cb]]$p.value<=0.01 & list_comb_lsp[[cb]]$Period>100,]}
comb.kr<-cbind(all_combs,data.frame(nb.rhytmic.keggs = unlist(lapply(list_comb_lsp_filtered,nrow))))
comb.kr$string=apply(comb.kr[,1:7], 1, function(x) paste(x[!is.na(x)], collapse = " ") )

ggplot(data = comb.kr, aes(y = string, x = nb.rhytmic.keggs))+
  facet_grid(NB~., space = "free", scale = "free")+
  scale_x_continuous(expand = c(0,0))+
  geom_bar(stat = "identity", fill = "steelblue4" )+
  labs(y = "Years included in the\ncomputation of Rythmicity", x = "# of significantly rhythmic KEGGs\n(over 8374 KEGGs)")+
  theme_minimal()+
  theme(strip.text.y = element_text(angle= 0, hjust = 0))

# reshape to long format
m.comb.kr<-na.omit(reshape2::melt(comb.kr[,-ncol(comb.kr)], id.vars = c("NB", "nb.rhytmic.keggs")))
table(m.comb.kr$value)
ag.comb.kr<-aggregate(nb.rhytmic.keggs ~ value+NB, data = m.comb.kr, FUN = mean)

# NB of year combinations
library(pgirmess)
nb=3
kruskalmc(m.comb.kr[m.comb.kr$NB == nb,"nb.rhytmic.keggs"], m.comb.kr[m.comb.kr$NB == nb,"value"], alpha = 0.05)
my_comparisons <- list( c("2010", "2013"), c("2010", "2014"), c("2010", "2015"))

ggplot(data = m.comb.kr[m.comb.kr$NB == nb,], aes(x =value , y = nb.rhytmic.keggs) )+
  geom_boxplot()+
  stat_compare_means(label.y.npc = "bottom", label.x.npc = "left")+ # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons )+ # Add pairwise comparisons p-value
  labs(y = "# of significantly rhytmic KEGG per combination", x = "Year included", title = paste0("For combinations of ", nb, " years"))+
  theme_minimal()+
  theme(legend.position = "none",panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

ggplot(data = m.comb.kr, aes(x =value , y = nb.rhytmic.keggs) )+
  facet_grid(NB~., scale = "free")+
  geom_boxplot()+
  labs(y = "# of significantly rhytmic KEGG per combination", x = "Year included")+
  theme_get()+
  theme(legend.position = "none",panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),strip.text.y = element_text(angle = 0, hjust = 0))

# Let's explore the patterns of the environment across these Years
env.somlit<-read.csv("~/Desktop/PETRIMED/DATA/ENV/SOLA_hydro_2001-2024.csv")
env<-env.somlit[,grep("q|DN|DC", colnames(env.somlit),invert = TRUE)]
env<-env[, c(match(c("DATE","month", "year"), colnames(env)), 11:23)]
env<-env[env$year %in% Years,]
menv<-reshape2::melt(env, id.vars = c("DATE", "month", "year") )
menv$Date<-as.Date(menv$DATE, format = "%Y-%m-%d")
menv$JD<-as.numeric(format(menv$Date, "%j"))
menv$Year<-factor(as.character(menv$year), levels = rev(Years), ordered = TRUE)

ggplot(data = menv, aes(x = JD, y = value, color = Year, group = Year) )+
  facet_grid(variable~., scale = "free")+
  geom_line(linewidth = 0.6)+
  labs(x = "", y="")+
  scale_x_continuous(expand = c(0,0))+
  scale_color_tron()+
  theme_minimal()+
  theme(strip.text.y = element_text(angle= 0, hjust = 0))





