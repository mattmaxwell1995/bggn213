Class 13: Genomics part I
================

Load MXL SNP analysis.

``` r
snp_data <- read.csv("MXL_snp_analysis.csv")
summary(snp_data)
```

    ##  Sample..Male.Female.Unknown. Genotype..forward.strand.
    ##  NA19648 (F): 1               A|A:22                   
    ##  NA19649 (M): 1               A|G:21                   
    ##  NA19651 (F): 1               G|A:12                   
    ##  NA19652 (M): 1               G|G: 9                   
    ##  NA19654 (F): 1                                        
    ##  NA19655 (M): 1                                        
    ##  (Other)    :58                                        
    ##        Population.s. Father Mother
    ##  ALL, AMR, MXL:64    -:64   -:64  
    ##                                   
    ##                                   
    ##                                   
    ##                                   
    ##                                   
    ## 

Determine percentage of GG variants from snp\_data. Use table to sort into bins and calculate percentage.

``` r
table(snp_data$Genotype..forward.strand.) / nrow(snp_data) * 100
```

    ## 
    ##     A|A     A|G     G|A     G|G 
    ## 34.3750 32.8125 18.7500 14.0625

FASTQ Quality Scores
--------------------

``` r
#install.packages("seqinr")
#install.packages("gtools")

library(seqinr)
library(gtools)
chars <- s2c("DDDDCDEDCDDDDBBDDDCC@")
chars
```

    ##  [1] "D" "D" "D" "D" "C" "D" "E" "D" "C" "D" "D" "D" "D" "B" "B" "D" "D"
    ## [18] "D" "C" "C" "@"

``` r
phred <- asc( chars ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

Load RNA seq files of 'GG' SNP patient into our environment. Let's see if this SNP changes gene expression.

``` r
seq1 <- read.csv("HG00109_1.fastq")
seq2 <- read.csv("HG00109_2.fastq")
#tail(seq1, col=10, rows=7)
summary(seq1)
```

    ##                                                       X.HW.ST546.136.D0HWFACXX.5.1101.10675.54538
    ##  +                                                                          : 3863               
    ##  TTGGAGTCTGCAGAGGAACGGCGTGAGCGAGAACAGCAGGACTTGGAGTTTGCCAAGGAGATGGCAGAAGATGAT:   13               
    ##  ACAAGGACTTGGAGTCTGCAGAGGAACGGCGTGAGCGAGAACAGCAGGACTTGGAGTTTGCCAAGGAGATGGCAG:   11               
    ##  CGAAGATGCAGAGTTCATTGTTGCCAAGGCCATCCGGGATGGTGTCATTGAGGCCAGCATCAACCACGAGAAGGG:   11               
    ##  CGATGACCAACGCCCTTCGCAAGGCCCCTCAGCACACAGCTGTCGGCTTCAAACAGACGGTGCACAAGCTTCTCA:   11               
    ##                                                                             :   10               
    ##  (Other)                                                                    :11679

Genotype Based Expression Levels
--------------------------------

Population scale analysis of SNP for Asthma. Load data from url with read.table(), looks better than read.csv

``` r
expr <- read.table("https://bioboot.github.io/bggn213_W19/class-material/rs8067378_ENSG00000172057.6.txt")
head(expr)
```

    ##    sample geno      exp
    ## 1 HG00367  A/G 28.96038
    ## 2 NA20768  A/G 20.24449
    ## 3 HG00361  A/A 31.32628
    ## 4 HG00135  A/A 34.11169
    ## 5 NA18870  G/G 18.25141
    ## 6 NA11993  A/A 32.89721

``` r
table(expr$geno)
```

    ## 
    ## A/A A/G G/G 
    ## 108 233 121

``` r
inds.gg <- expr$geno == "G/G"
summary( expr[inds.gg,"exp"] )
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
inds.ag <- expr$geno == "A/G"
summary( expr[inds.ag,"exp"] )
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds.aa <- expr$geno == "A/A"
summary( expr[inds.aa,"exp"] )
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
boxplot(exp ~ geno, data = expr, notch = TRUE)
```

![](class13_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# Boxplot with the data shown
#ggplot(expr, aes(geno, exp, fill=geno)) + 
 # geom_boxplot(notch=TRUE, outlier.shape = NA) + 
  #geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.4)
```
