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

###1st analyses with script: sola_2009_2014_kegg.R
###1) Based on three objects, this script first imports and formats:
- 1 sample info table (with dates)
- 1 KEGG Abundance table
- 1 KEGG Functional Annotation table
- 1 Genes Abundance table
- 1 Genes-KEGG annotation table

###2) Lomb-Scargle Periodogram computation and analysis of results



