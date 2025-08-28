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


