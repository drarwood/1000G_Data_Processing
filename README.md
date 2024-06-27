# 1000G Data Processing
Repository containing data and code for processing 1000G data
---
### Association vs Allele Frequency differences (build 38)
A script [`Assoc_P_vs_AF_Diffs.R`](assoc_p_vs_af_diffs/Assoc_P_vs_AF_Diffs.R) and dataset containing sparse matrices of allele frequency differences within the
5 broad genetic ancestral groups (AFR, AMR, EAS, EUR, SAS).   

** Ancestry-specific matrices yet to be provided ** 

R libraries required: `data.table`, `r2r`, `ggplot2`, `optparse`.   

Help using this R script call be called using `-h` or `--help`

```
Rscript Assoc_P_vs_AF_Diffs.R -h

Usage: Assoc_P_vs_AF_Diffs.R [options]

Options:
	-s STUDY, --study=STUDY
		study name (no spaces) (required)

	-g GWAS, --gwas=GWAS
		REGENIE results filename (required)

	-a ANCESTRY, --ancestry=ANCESTRY
		Genetic ancestry of study: AFR|AMR|EAS|EUR|SAS (required)

	-r REFDIR, --refdir=REFDIR
		Main directory holding 1000G frequency data (required)

	-o OUTDIR, --outdir=OUTDIR
		Output directory for results. Default: ./

	-h, --help
		Show this help message and exit
```
Example command:
```
Rscript Assoc_P_vs_AF_Diffs.R \
  --study     MY_STUDY \
  --gwas      gwas_results.txt \
  --ancestry  EUR \
  -r          /path/to/1000g_build38_pop_freq_diffs_by_subancestry/ \
  -o          /path/to/output/directory/
```
---


### Allele Frequency Differences (build 38)
Differences in allele frequencies across the 5 super-ancestries and 26 sub-ancestries in the 1000G Project have been calculated 
for 78,122,255 autosomal variants in 2,584 indivuals. <br/><br/>
Data is stored in a Zenodo repo here: https://zenodo.org/doi/10.5281/zenodo.11635380
<br/><br/>
The data is split by chromosome with the file extensions `.pop`, `.var`, `.bin`. 
The file specification can be found [here](https://github.com/drarwood/1000G_Data_Processing/blob/master/1000G_pop_freq_diffs_file_format.pdf).
<br/><br/>
An example Python script and SNP list for lookup have been provided for use. For example, to lookup frequency differences of the subset 
of 12,111 height variants published in Yengo et al on chromosome 1 across ancestral groups in the 1000G project, the following command
can be used:

```
  python GetFreqDiffs.py \
    --freq-diffs 1000g_build38_pop_freq_diffs_1 \
    --vars 12111_height_snps.txt \
    --out 12111_height_snps_1000g_freq_diffs_1.txt
```
---


