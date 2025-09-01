# ts_rhythm

OG files from Lidia were pre-treated in bash before R import:

```
# Just gene-KEGG annotation
awk '{print $1, $2}' kegg_BBMOSOLA-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts.tbl > gene_kegg.tbl

# Just SOLA abundance table
awk '{for(i=87;i<=176;i++) printf "%s%s", $i, (i==176?RS:FS)}' kegg_BBMOSOLA-GC_250bp_gene.lengthNorm.SingleCopyGeneNorm.counts.tbl > SOLA_2009_2014_gene.tbl

# stitch them back together
paste -d ' ' gene_kegg.tbl SOLA_2009_2014_gene.tbl > SOLA_2009_2014_gene_annot.tbl

# filtered out empty rows (absent genes)
awk 'NR==1 {print; next} {sum=0; for(i=3;i<=NF;i++) sum+=$i; if(sum!=0) print}' SOLA_2009_2014_gene_annot.tbl > SOLA_2009_2014_gene_annot_filtered.tbl

# convert delimiters to tab sep
awk -v OFS='\t' '{$1=$1}1' SOLA_2009_2014_gene_annot_filtered.tbl > SOLA_2009_2014_gene.tbl

# remove transition files
rm SOLA_2009_2014_gene_annot.tbl
rm mOLA_2009_2014_gene_annot_filtered.tbl
rm gene_kegg.tbl
```

### 1st analyses with script: sola_2009_2014_kegg.R
### 1) Based on three objects, this script first imports and formats:
- 1 sample info table (with dates)
- 1 KEGG Abundance table
- 1 KEGG Functional Annotation table
- 1 Genes Abundance table
- 1 Genes-KEGG annotation table

### 2) Lomb-Scargle Periodogram computation and analysis of results:
Results table for the 8375 KEGGs looks like this:
```
      kegg    PNmax  Period p.value Phase
1   K00525 0.207330 362.000   0.010   Mar
2   K02014 0.140910 181.000   0.212   Oct
3   K02335 0.186610 362.000   0.061   May
4   K03046 0.284880 362.000   0.000   Apr
5   K03043 0.223030 362.000   0.000   Apr
6   K02337 0.088779  38.394   0.939   May
7   K01915 0.180820 362.000   0.020   Apr
8   K00605 0.251340 362.000   0.000   Apr
9   K01952 0.124520  31.675   0.424   May
10  K01652 0.111540 362.000   0.535   Apr
```
2857 Keggs were found to be rhythmic (PNmax > 0.1 and p.value < 0.01)
5517 Keggs were not significantly rhythmic (over the 7 years time period)

As expected, rhythmicity was well correlated to KEGG's total abundance, with rarer KEGGs having sparser distribution (potentially absent for long time periods):
<img width="444" height="384" alt="Screenshot 2025-08-28 at 17 12 20" src="https://github.com/user-attachments/assets/11f87fab-a307-41f9-934c-836d679f230a" />

We also can explore the differences in rhythmicity across KEGGs functional categories:
<img width="908" height="641" alt="Screenshot 2025-08-28 at 17 13 52" src="https://github.com/user-attachments/assets/763b4c2f-9ad7-4606-be39-2ee9abad08ba" />

We also investigated the most rhythmic KEGG's patterns over the time series:
<img width="621" height="326" alt="Screenshot 2025-08-28 at 17 03 31" src="https://github.com/user-attachments/assets/6ebbcd29-f16c-4502-a4c7-3897da3bca19" />

Its functional annotation is not very revealing:
```
egg Gene names          kegg annotation GN1 GN2 GN3 GN4 GN5 GN6 GN7 GN8 GN9                                    KH1
3292 K06287        maf septum formation protein maf                                 09190 Not Included in Pathway or Brite
                                                      KH2               KH3  KH4  KH5  KH6  KH7  KH8  KH9 KH10 KH11 KH12 KH13 KH14 KH15 KH16
3292 09193 Unclassified: signaling and cellular processes 99978 Cell growth <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA>
     KH17 KH18 KH19 KH20 KH21 KH22 KH23 KH24 KH25 KH26 KH27 KH28 KH29 KH30 KH31   tot.ab nb.genes
3292 <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> <NA> 48.81385     5067
```
Next we check the patterns of the genes annotated to this KEGG (the KEGG's composition).
5067 genes were annotated to this specific KEGG, here are their cumulated abundance over the time-series:
<img width="889" height="414" alt="Screenshot 2025-08-28 at 17 06 34" src="https://github.com/user-attachments/assets/42562468-3b78-4b1a-8167-599254764534" />

We note that the abundances of genes annotated to this KEGG is not equal to the abundance of the KEGG (see here across samples):
<img width="368" height="347" alt="Screenshot 2025-08-28 at 17 02 19" src="https://github.com/user-attachments/assets/98c795ad-6e2d-4c50-9226-e5a9dd8633a4" />

### Randomize rhythmicity with script: randomization_rhythm.R
### 1) Create list of Year combination
The idea of this script is to compute LSP rhythmicity for all KEGGs based on all possible random combinations of Years.
First we create a matrix of years combinations, from 2-Years to 7-Years:
```
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
```
Which gives us the following table with 120 combinations, in which each year is represented equitatively (63 times):
```
      Year1  Year2  Year3  Year4  Year5  Year6  Year7    NB
  1:   2009   2010   <NA>   <NA>   <NA>   <NA>   <NA>     2
  2:   2009   2011   <NA>   <NA>   <NA>   <NA>   <NA>     2
  3:   2009   2012   <NA>   <NA>   <NA>   <NA>   <NA>     2
  4:   2009   2013   <NA>   <NA>   <NA>   <NA>   <NA>     2
  5:   2009   2014   <NA>   <NA>   <NA>   <NA>   <NA>     2
 ---                                                       
116:   2009   2010   2011   2013   2014   2015   <NA>     6
117:   2009   2010   2012   2013   2014   2015   <NA>     6
118:   2009   2011   2012   2013   2014   2015   <NA>     6
119:   2010   2011   2012   2013   2014   2015   <NA>     6
120:   2009   2010   2011   2012   2013   2014   2015     7
```
### 2) Run LSP for all combinations
We then create 120 subdatasets using only each of these Years combinations, we create a fake consecutive "Julian Days" column that run over the year-gaps.
We can then run the LSP algorithm and save the results for all 8374 KEGGs over the 120 Year-combinations:
```
list_comb_lsp<-list()
for (cb in seq_len(length(list_comb_info))){
  
  # merge info on the time series with KEGGs abundances
  lsp.keggs<-merge(list_comb_info[[cb]], t(keggs.ab), by.x = "ID", by.y = "row.names" )
  
  # loop the command over all Keggs
  r=NULL
  for (k in 11:ncol(lsp.keggs)){
    rl<-summary(randlsp(x = lsp.keggs[, c(9, k)], repeats = 99, type = "period", from = 30, to = 1000, trace = FALSE, plot = FALSE ))
    r<-rbind(r,data.frame(kegg = colnames(lsp.keggs)[k], PNmax = as.numeric(rl[9,]), Period = as.numeric(rl[10,]), p.value = as.numeric(rl[13,]), Phase = format(lsp.keggs[which.max(lsp.keggs[,k]),"Date"], "%b") ))
    print(paste(cb, k, ncol(lsp.keggs), sep = " / "))
  }
  list_comb_lsp[[paste(cb)]]<-r
}
saveRDS(list_comb_lsp, "list_comb_lsp.RData")
```
### 3) Analysis of the results:
```
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
```
We studied which combinations yielded the highest PNmax (rhytmicity) for each KEGG.
We can also study which Years' removal generated the highest PNmax over all KEGGs:
```
# which Year's removal generate the highest PNmax?
missing_years<-list()
for (y in 1:120){missing_years[[y]]<-Years[!Years %in% list_comb_info[[y]]$Year]}

all.KO=NULL
for (k in rownames(comb.PN) ){all.KO<-c(all.KO,missing_years[[comb.PN[k,]]])}
table(all.KO)
```
which gives us:
```
table(all.KO)
all.KO
2009 2010 2011 2012 2013 2014 2015 
5599 5030 5615 5911 6691 6364 6494 
```
> It seems that the absence of the years 2013 and 2015 generated the highest PNmax for most KEGGs (6691 and 6494 over 8734).
> Could these years harbor more perturbations than other years causing a lack of rythmicity for most KEGGs? 

Let's check the variability of environmental variables across years (source: SOMLIT SOLA):
<img width="682" height="889" alt="Screenshot 2025-09-01 at 15 20 20" src="https://github.com/user-attachments/assets/57f1beeb-11df-4e79-a454-6aca330061ed" />


