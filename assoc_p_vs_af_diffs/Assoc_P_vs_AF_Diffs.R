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
  x <- x[,c("ID1","ID2","CHROM","GENPOS","ALLELE_RAISING","ALLELE_LOWERING",
            "PBIN")]
  
  # return DF
  return (x)
  
}


LoadPopPairLabelsFile <- function(f) {
  
  # This function loads a .pop file
  # Parameters:
  #   f (string) : filename of .pop
  # Returns:
  #   x (df) : df containing data in .pop file
  
  x <- fread(f, data.table=FALSE, 
             colClasses=c("integer", "character"),
             col.names=c("index", "population_pair"))
  
  # Return DF
  return (x)
  
}


Load1000GVarData <- function(d, p, c) {
  
  # This function loads the 1000G variant data associated with .bin file
  # Parameters:
  #   d (string) : directory holding variant data
  #   p (string) : broad ancestral group
  #   c (int)    : chromomsome
  # Returns:
  #   x (df) : df holding all variant data across chromosomes
  
  print(paste0("Loading 1000G variants for chromosome ", as.character(c)))
  x <- fread(paste0(d,"/",p,"/1000g_build38_pop_freq_diffs_",p,"_",c,".var"), 
             col.names=c("chr", "pos", "primary_allele", "other_allele", 
                         "start_byte", "non_zero_diffs"), data.table=FALSE, 
             colClasses=c("character", "numeric", "character","character",
                          "numeric", "numeric"))
  x$id <- paste(x$chr, x$pos, x$primary_allele, x$other_allele, sep=":")
  
  # Return DF
  return (x)
  
}


ExtractAlleleFreqDiffs <- function(d, p, c, g, n, h) {
  
  # This function builds and returns a matrix (40 x N ancestry-pairs)
  # storing the total difference of allele frequencies between the N-ancestry
  # pairs for the 40 P-value bins for a given chromosome.
  # Parameters:
  #  d (string)  : reference data directory sotring .bin file
  #  p (string)  : population [AFR|AMR|EAS|EUR|SAS]
  #  c (int)     : chromosome
  #  g (df)      : DF holding variant byte addresses etc (based on .var file)
  #  n (int)     : number of population-pairs for given ancestry
  #  h (hashmap) : p_bin -> array / matrix column index
  # Returns:
  #  M (matrix)  : total trait raising frequency differences across 40 p-value
  #                bins for each of the ancestry pairs
  
  # Create matrix of zeros where rows represent the p-value bins and columns 
  # represent the population pairs. So columns updated based on p-value bin
  M <- matrix(0, nrow=40, ncol=n)
  
  # Open file
  bin <- file(paste0(d,"/",p,"/1000g_build38_pop_freq_diffs_",p,"_",c,".bin"), 
              "rb")
  
  # check magic number and format. Stop if not as expected
  magic_number <- readBin(bin, "raw", 6)
  format       <- readBin(bin, "int", 1, size=1)
  stopifnot(rawToChar(magic_number) == "ARWOOD" && as.numeric(format) == 1)

  
  # cycle through the rows of the DF storting the start byte & N non-zeros
  for (i in 1:nrow(g)) {
    
    # create vector of zeros of length = n (number of population pairs)
    v <- numeric(n)
    
    # check non-zero diffs exist
    if (g$non_zero_diffs[i] > 0) {
      
      # Go to start byte of this variant in diffs file
      seek(bin,g$start_byte[i])
      
      # Read block of 6 x N non-zero diffs bytes into vector called buffer
      buffer = readBin(bin, "raw", n=6*g$non_zero_diffs[i], size=1, endian="little")
      
      offset = 1
      # Go through the index and diff value pairs
      for (j in 1:g$non_zero_diffs[i]) {
        # Obtain the ancestry-pair index the non-zero diff associates with
        index <- readBin(buffer[offset:(offset+1)], "int", n=1, size=2, endian="little")
        # Obtain the ancesty-pair allele frequency difference
        diff  <- readBin(buffer[(offset+2):(offset+5)], "double", n=1, size=4, endian="little")
        # add 1 as population pairs zero-indexed but arrays start at 1
        v[index+1] <- diff
        offset <- offset+6
      }

      # check whether we need to flip the freqs if trait raising allele not
      # the same as the allele the frequency differences are based on
      if (g$FLIP_DIFF[i]) {
        v <- v*-1
      }
      
      # Update matrix by adding the values of the vector to the row values in
      # the matrix, where the row represents the P-value bin of association.
      pbin_index <- as.numeric(h[[as.character(g$PBIN[i])]])
      M[pbin_index,] <- M[pbin_index,] + v
      
    }
    
  }
  
  # close the file
  close(bin)
  
  # return the matrix
  return(M)
  
}

### FUNCTIONS END ###



### "MAIN" ---

# 1.  Define and obtain arguments from command line
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

# Return command line arguments to user
print("Study_P_vs_1000G_Freq_Diffs.R")
print(paste0("-s  --study     ", opt$s))
print(paste0("-g  --gwas      ", opt$g))
print(paste0("-a  --ancestry  ", opt$a))
print(paste0("-r  --refdir    ", opt$r))
print(paste0("-o  --outdir    ", opt$o))

# Quit is any values NA
if (is.na(opt$s)||is.na(opt$g)||is.na(opt$a)||is.na(opt$r)) {
  stop("Please ensure arguments required are provided. Use -h to see options")
}
# Check valid ancestry code provided
if (!(opt$a=="AFR"||opt$a=="AMR"||opt$a=="EAS"||opt$a=="EUR"||opt$a=="SAS")) {
  stop("Please provide a value ancestry code [AFR|AMR|EAS|EUR|SAS]")
}

# Go to directory for output of data
setwd(opt$o)

### 1. Load and extract required GWAS data
gwas <- LoadRegenieFile(opt$g)

### 2. Load population labels using chromosome 1 file only
pops_1000G <- LoadPopPairLabelsFile(paste0(opt$refdir, "/", opt$ancestry,
                                          "/1000g_build38_pop_freq_diffs_",
                                           opt$ancestry, "_1.pop"))

### 3. Get allele  frequencies from the 1000G for variants in the 1000G data set

# 3a. Create a hashmap for PBIN -> column indices in arrays and matrix
# 0.025 = 1, 0.05  = 2, ..., 0.995 = 39, 1 = 40
h <- hashmap()
i <- 0
for (bin in seq(0.025,1, 0.025)) {
  i <- i+1
  h[[as.character(bin)]] <- i
}

# 3b. Create matrix to store total of all trait raising allele frequency diffs
# across each of the P-value bins for each ancestry pair
total_diff_matrix  <- matrix(0, nrow=40, ncol=nrow(pops_1000G))
total_count_vector <- numeric(40)


# 3c. Cycle through the chromosomes of the allele frequency diff. dataset
for (chr in 21:22) {

  # Load variants for which we have allele frequency diffs in 1000G
  vars_1000G <- Load1000GVarData(opt$refdir, opt$ancestry, chr)

  # Now merge ID1 and ID2 in gwas to the var_1000G DF
  merge1 <- merge(gwas, vars_1000G, by.x = "ID1", by.y = "id")
  merge2 <- merge(gwas, vars_1000G, by.x = "ID2", by.y = "id")
  x <- rbind(merge1, merge2)

  # Determine if frequency differences need to be flipped
  x$FLIP_DIFF <- ifelse(x$ALLELE_RAISING == x$primary_allele, F, T)

  # Now cycle through the variants in X, extracting the allele freqs
  c_diff_matrix <- ExtractAlleleFreqDiffs(opt$refdir, opt$ancestry, chr,
                                         x, nrow(pops_1000G), h)

  # Create a count of the number of variants in each P-value bin for this chr
  chr_count_vector <- numeric(40)
  for (i in seq(0.025,1, 0.025)) {
    chr_count_vector[h[[as.character(i)]]] <-
      nrow(subset(x, PBIN == as.character(i)))
  }

  # Update the total_diff_matrix by adding the chromosome-specific diff matrix
  total_diff_matrix <- total_diff_matrix+c_diff_matrix

  # Update the total_count_vector by adding the chromosome-specific counts of
  # variants in each P-value bin
  total_count_vector <- total_count_vector+chr_count_vector

}

### 4. Divide columns of total diffs matrix by vector of total counts & store DF 
avg_freq_diffs_by_p_value_bin <- data.frame(seq(0.025,1,0.025),
                                           total_diff_matrix/total_count_vector)
# Update names of df
names(avg_freq_diffs_by_p_value_bin) <- c("pbin",
                                          toupper(pops_1000G$population_pair))

### 5. Output average frequency differences data to file
write.table(avg_freq_diffs_by_p_value_bin, 
            file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study, 
                        "_avg_freq_diffs_by_p_value_bin.txt"),
            quote = F, sep = "\t", row.names=F)


### 6. Output simple linear regression results of average freq. diff ~ P bin
cat(paste("ancestry_pair","effect", "se", "tval", "pval", sep="\t"), "\n", 
    file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study, 
                "_pval_vs_freq_diffs_stats.txt"))

for (i in 2:ncol(avg_freq_diffs_by_p_value_bin)) {
  pair  <- names(avg_freq_diffs_by_p_value_bin[i])
  model <- lm(avg_freq_diffs_by_p_value_bin[,i] ~
                avg_freq_diffs_by_p_value_bin$pbin)
  effect <- summary(model)$coefficients[2, 1]
  se <- summary(model)$coefficients[2, 2]
  t <- summary(model)$coefficients[2, 3]
  p <- summary(model)$coefficients[2, 4]
  cat(paste(pair,effect, se, t, p, sep="\t"), "\n",
      file=paste0(opt$outdir,"/", opt$ancestry, "_", opt$study,
                  "_pval_vs_freq_diffs_stats.txt"), 
      append=T)
}


### 7. Create plot showing correlations

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
       
# each sub-ancestry pair
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

### END SCRIPT ###