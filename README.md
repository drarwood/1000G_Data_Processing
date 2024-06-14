# 1000G Data Processing
Repository containing data and code for processing 1000G data
---
### Allele Frequency Differences
Differences in allele frequencies across the 5 super-ancestries and 26 sub-ancestries in the 1000G Project have been calculated 
for 78,1222,255 autosomal variants in 2,584 indivuals. 
Data is stored in a Zenodo repo here: https://zenodo.org/doi/10.5281/zenodo.11635380. 
The data is split by chromosome with the file extensions `.pop`, `.var`, `.bin`. 
An example Python script and SNP list for lookup have been provided for use. For example, to lookup frequency differences of the subset 
of 12,111 height variants published in Yengo et al on chromosome 1 across ancestral groups in the 1000G project, the following command
can be used:

```
  python GetFreqDiffs.py \
    --freq-diffs 1000g_build38_pop_freq_diffs_1 \
    --vars 12111_height_snps.txt \
    --out 12111_height_snps_1000g_freq_diffs.txt
```

