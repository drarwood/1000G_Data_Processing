# Allele Frequency Differences vs P-value bin
# Andrew Wood (a.r.wood@exeter.ac.uk)
# The purpose of this script is to obtain allele frequency differences between 
# ancestry groups in 1000G within one of the five super-populations

# Import required libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(r2r))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))

### Functions ###

CheckArgs <- function(o) {
  
  # This function checks arguments input by user
  # Parameters:
  #  o (list of options) : provides access to options from command line
  # Returns:
  #  Not applicable.
  
  # Return command line arguments to user
  print("Study_P_vs_1000G_Freq_Diffs.R")
  print(paste0("-s  --study     ", o$s))
  print(paste0("-g  --gwas      ", o$g))
  print(paste0("-a  --ancestry  ", o$a))
  print(paste0("-r  --refdir    ", o$r))
  print(paste0("-o  --outdir    ", o$o))

  # Quit is any values NA
  if (is.na(opt$s)||is.na(opt$g)||is.na(opt$a)||is.na(opt$r)) {
    stop("Please ensure arguments required are provided. Use -h to see options")
  }
  
  # Check valid ancestry code provided
  if (!(opt$a=="AFR"||opt$a=="AMR"||opt$a=="EAS"||opt$a=="EUR"||opt$a=="SAS")) {
    stop("Please provide a value ancestry code [AFR|AMR|EAS|EUR|SAS]")
  }
  
}

LoadRegenieFile <- function(f) {

  # This function loads in a REGENIE results file
  # Parameters:
  #   f (character): file name of REGENIE results file 
  # Returns:
  #   x (df): data frame based on REGENIE GWAS data
  # Current explicit assumptions:
  #   Results from all autosomal chromosomes present in REGENIE file

  # read in the REGENIE file - assumed columns and order: 
  # CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTR
  x <- fread(f, data.table = FALSE, header = TRUE, sep=" ",
             colClasses=c("character", "character", "character", 
                          "character", "character", "numeric", 
                          "numeric", "numeric", "character",
                          "numeric", "numeric", "numeric", 
                          "numeric", "character"))

  # filter out any results with INFO < 0.3 and MAF < 0.1%
  x <- subset(x, INFO>=0.3 & A1FREQ >= 0.001 & A1FREQ <= 0.999)

  # determine trait raising and lowering alleles
  x$ALLELE_RAISING<-ifelse(x$BETA >= 0, x$ALLELE1, x$ALLELE0)
  x$ALLELE_LOWERING<-ifelse(x$ALLELE1 == x$ALLELE_RAISING, x$ALLELE0, x$ALLELE1)

  x$ID1<-paste0(x$CHROM,":",x$GENPOS,":",x$ALLELE_RAISING,":",x$ALLELE_LOWERING)
  x$ID2<-paste0(x$CHROM,":",x$GENPOS,":",x$ALLELE_LOWERING,":",x$ALLELE_RAISING)
  
  # convert -LOG10(P) to P=value ceiling to nearest 0.025
  x$PBIN <- as.character(ceiling((10^-x$LOG10P) * 40)/40)
    
  # keep only the columns we need
  x <- x[,c("ID1","ID2","ALLELE_RAISING","ALLELE_LOWERING","PBIN")]
  
  # return DF
  return (x)
  
}

Load1000GVarData <- function(d, p, c) {
  
  # This function loads the 1000G frequency differences file
  # Parameters:
  #   d (string) : directory holding variant data
  #   p (string) : broad ancestral group
  #   c (int)    : chromomsome
  # Returns:
  #   x (df)      : df holding all variant data across chromosomes
  
  print(paste0("Loading 1000G variants for chromosome ", as.character(c)))
  x <- fread(paste0(d,"/",p,"/1000g_build38_pop_freq_diffs_",p,"_",c), 
             data.table=FALSE)
  x$id <- paste(x$chr, x$pos, x$a1_freq_allele, x$a2_freq_allele, sep=":")

  # Return DF
  return (x)
  
}

### FUNCTIONS END ###


### MAIN ###

### 1.  Define and obtain arguments from command line
option_list <- list( 
  make_option(c("-s", "--study"), type="character", default=NA, 
              help="study name (no spaces) (required)"),
  make_option(c("-g", "--gwas"), type="character", default=NA,
              help="REGENIE results filename (required)"),
  make_option(c("-a", "--ancestry"), type="character", default=NA,
              help="Genetic ancestry of study: AFR|AMR|EAS|EUR|SAS (required)"),
  make_option(c("-r", "--refdir"), type="character", default=NA, 
              help="Main directory holding 1000G frequency data (required)"),
  make_option(c("-o", "--outdir"), type="character", default="./",
              help="Output directory for results. Default: ./"))

opt <- parse_args(OptionParser(option_list=option_list))


### 2. Check command line arguments
CheckArgs(opt)

 
### 3. Load and extract required GWAS data
print("Loading GWAS results from REGENIE")
gwas <- LoadRegenieFile(opt$g)


### 4. Load population labels using chromosome 1 file only
h <- fread(paste0(opt$refdir, "/", opt$ancestry, 
                  "/1000g_build38_pop_freq_diffs_EUR_21"), 
                  header=T, nrows=0, data.table=FALSE)
pops_1000G = colnames(h)[5:ncol(h)]


### 5. Get allele  frequencies from the 1000G for variants in the 1000G data set

## 5a. Create a hashmap for PBIN -> column indices in arrays and matrix
# 0.025 = 1, 0.05  = 2, ..., 0.995 = 39, 1 = 40
h <- hashmap()
i <- 0
for (bin in seq(0.025,1, 0.025)) {
  i <- i+1
  h[[as.character(bin)]] <- i
}

## 5b. Create matrix to store total of all trait raising allele frequency diffs
# across each of the P-value bins for each ancestry pair
total_diff_matrix  <- matrix(0, nrow=40, ncol=length(pops_1000G))
total_count_vector <- numeric(40)

# 5c. Cycle through the chromosomes of the allele frequency diff. dataset
for (chr in 1:22) {

  # Load variants for which we have allele frequency diffs in 1000G
  vars_1000G <- Load1000GVarData(opt$refdir, opt$ancestry, chr)
  
  # Now merge ID1 and ID2 in gwas to the var_1000G DF
  merge1 <- merge(gwas, vars_1000G, by.x = "ID1", by.y = "id")
  merge2 <- merge(gwas, vars_1000G, by.x = "ID2", by.y = "id")
  x <- rbind(merge1, merge2)
  
  # Determine if frequency differences need to be flipped. Set as multiplier
  x$FLIP_DIFF <- ifelse(x$ALLELE_RAISING == x$a1_freq_allele, 1, -1)

  # cycle through the p-value bin
  for (bin in seq(0.025,1, 0.025)) {
    
    # Subset DF based on current P-value bin
    xs = subset(x, PBIN == as.character(bin))
    
    # Generate freq. diff matrix for all variants on chr
    # Get column 10 until column n-1 - these columns hold the freq diffs in x
    M = data.matrix(xs[10:(ncol(xs)-1)])
    # Flip as required
    M = M*xs$FLIP_DIFF
    
    # Get totals of freq differences for the difference ancestry pairs
    total_freq_diffs = colSums(M)
    
    # Update total_diff_matrix used to capture all chromosomes
    total_diff_matrix[h[[as.character(bin)]], ] = 
      total_diff_matrix[h[[as.character(bin)]], ] + total_freq_diffs

    # Update total_count_matrix used to capture all chromosomes
    total_count_vector[h[[as.character(bin)]]] = nrow(M)
    
  }
  
}

## 5d. Divide columns of total diffs matrix by vector of total counts & store DF 
avg_freq_diffs_by_p_value_bin <- data.frame(seq(0.025,1,0.025),
                                           total_diff_matrix/total_count_vector)

# Update names of df
names(avg_freq_diffs_by_p_value_bin) <- c("pbin", toupper(pops_1000G))


### 6. Output average frequency differences data to file
print("Writing mean frequency differences across P-value bins to file")
write.table(avg_freq_diffs_by_p_value_bin, 
            file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study, 
                        "_avg_freq_diffs_by_p_value_bin.txt"),
            quote = F, sep = "\t", row.names=F)


### 7. Output simple linear regression results of average freq. diff ~ P bin
print("Writing statistics for frequency difference ~ P-value bin to file")

cat(paste("ancestry_pair","effect", "se", "tval", "pval", sep="\t"), "\n", 
    file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study, 
                "_pval_vs_freq_diffs_stats.txt"))

for (i in 2:ncol(avg_freq_diffs_by_p_value_bin)) {
  
  pair   <- names(avg_freq_diffs_by_p_value_bin[i])
  model  <- lm(avg_freq_diffs_by_p_value_bin[,i] ~
               avg_freq_diffs_by_p_value_bin$pbin)
  
  effect <- summary(model)$coefficients[2, 1]
  se     <- summary(model)$coefficients[2, 2]
  t      <- summary(model)$coefficients[2, 3]
  p      <- summary(model)$coefficients[2, 4]
  
  cat(paste(pair,effect, se, t, p, sep="\t"), "\n",
      file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study,
                  "_pval_vs_freq_diffs_stats.txt"), append=T)
}


### 7. Create plots showing correlations
print("Plotting frequency differences against strength of associations")
# Combined:
df <- melt(as.data.table(avg_freq_diffs_by_p_value_bin), id.vars="pbin")
p <- ggplot(df, aes(pbin,value, col=variable)) +
     geom_point() +
     geom_smooth(method='lm', formula= y~x, se = TRUE) +
     ggtitle(paste0("Trait raising allele frequency difference between all 
                    1000G sub-ancesties"))+
     labs(x = "P-value bin", y ="Trait Raising Allele Freq. Difference") +
     theme(plot.title = element_text(hjust = 0.5))
# output plot
suppressMessages(
  ggsave(filename = paste0(opt$outdir, "/", opt$ancestry, "_", opt$study, 
                         "_1000G_ALL_freq_diffs.pdf"), p)
)
       
# Each sub-ancestry pair
for (i in 2:ncol(avg_freq_diffs_by_p_value_bin)) {
  assign("pair", colnames(avg_freq_diffs_by_p_value_bin)[i])
  p <- ggplot(avg_freq_diffs_by_p_value_bin, aes(pbin, .data[[pair]])) +
    geom_point() +
    geom_smooth(method='lm', formula= y~x, se = TRUE) +
    ggtitle(paste0("Trait raising allele frequency difference between ", pair))+
    labs(x = "P-value bin", y ="Trait Raising Allele Freq. Difference") +
    theme(plot.title = element_text(hjust = 0.5))
  # output plot
  suppressMessages(
    ggsave(filename = paste0(opt$outdir, "/", opt$ancestry, "_", opt$study,
                             "_1000G_", pair, 
                             "_trait_raising_allele_freq_diffs.pdf"), p)
  )

}

# Final user feedback
print("Analysis finished")
