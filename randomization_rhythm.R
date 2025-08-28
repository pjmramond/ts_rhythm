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
library(seas)

# 1/ Import, format, and explore KEGG table
# -----------------------------------------

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
library(dplyr)
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

list_comb_lsp<-list()
for (cb in seq_len(length(list_comb_info))){
  
  # merge info on the time series with KEGGs abundances
  lsp.keggs<-merge(list_comb_info[[cb]], t(keggs.ab), by.x = "ID", by.y = "row.names" )
  
  # loop the command over all Keggs
  r=NULL
  for (k in 11:ncol(lsp.keggs)){
    rl<-summary(randlsp(x = lsp.keggs[, c(9, k)], repeats = 99, type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
    r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[k], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[13,]), Phase = format(lsp.keggs[which.max(lsp.keggs[,k]),"Date"], "%b") ))
    print(paste(c, k, ncol(lsp.keggs), sep = " / "))
  }
  list_comb_lsp[[paste(cb)]]<-r
}



