# SomaGene 
  
  some explanation and logo ans stuffs //TODO

  
  
## Installation 
  
### Requirements
The following R packages:
  * [data.table](https://github.com/Rdatatable/data.table)
  * [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
  * [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)
  
  

Clone the package from the github page in the directory DIR
```
git clone https://github.com/bcb-sut/SomaGene.git DIR

```

Then install the SomaGene using R function `install.packages()`

```r
# Install SomaGene from source
install.packages("DIR/SomaGene_src", repos = NULL, type = "source")
```


## Usage

First you should import the SomaGene library

```r
library(SomaGene)
```

Then you should read your input data (for example your gene list)
There's a small test data in github.

```r
table <- fread("DIR/test_data/test_table.tsv")
```

The dataset should contains four columns, ID, chromosome, start, and end.
For example the forrmat of test_table.tsv is like below:

| GeneID  | chr | start | end |
| ------------- | ------------- | ------------- | ------------- |
| ENSG00000223764	  | chr1  | 845735 | 864616 |


**The values of chromosome column should have the format of ['chr' + number of chromosomoe] like chr2 or chrX**




* `annotateBinary(input_id, input_chr, input_start, input_end, annot_id = NULL, annot_chr, annot_start, annot_end)`: This function annotates the input genomic regions with a given "binary" annotation. A "binary" annotation is simply a set of genomic regions without any extra attribute (e.g. the set of enhancers/promoters or the set of eQTLs).

                
  __inputs__:
  * **input_id**: Character vector defining the name of input genomic regions (e.g. gene id)
  * **input_chr**: 	Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
  * **annot_id**: Optional (NULL by default): Character vector defining the name or ID of each annotation entry. If NULL, the row numbers (indexes) are used as IDs.
  * **annot_chr**: Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **annot_start**: Numeric vector specifying the starting position of annotation entries.
  * **annot_end**: Numeric vector specifying the ending position of annotation entries.
  
   __outputs__:
   * **input_ID**: Character vector defining the name of input genomic regions (e.g. gene id)
   * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
   * **input_start**: Numeric vector specifying the starting position of input genomic regions.
   * **input_end**: Numeric vector specifying the ending position of input genomic regions.
   * **overlap_count**: Numeric vector specifying the count of overlapped annotations with that geneID.
   * **overlapping_annot_IDs**: Character vector specifying the IDs of overlapped annotations with that geneID (separated by comma).
   
   
There is a 'Fantom5_enhancer_phase1and2.txt' file for testing.
   
```r
annot <- fread("DIR/test_data/Fantom5_enhancer_phase1and2.txt")

# Then you run the annotateBinary function
result <- annotateBinary(input_id = table$geneID, input_chr =  table$chr, input_start =  table$start, input_end =  table$end, annot_id = rownames(annot), annot_chr = annot$chr, annot_start = annot$start, annot_end = annot$end)

```

The result is like below

| GeneID  | chr | start | end | overlap_count | overlapping_annot_IDs
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| ENSG00000223764	  | chr1  | 845735 | 864616 | 3 | 4,5,6 |




  
* `annotateMultiScore(input_id, input_chr, input_start, input_end, annot_chr, annot_start, annot_end, annot_sub_id, annot_sub_score)`: This function annotates the input genomic regions with a given "multi-score" annotation. A "multi-score" annotation specifies a set of genomic regions and assigns a set of sub-IDs and their corresponding sub-scores to each of them (e.g. the annotation of DNAse hypersensitive clusters by ENCODE).
                
  __inputs__:
  * **input_id**: Character vector defining the name of input genomic regions (e.g. gene id)
  * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
  * **annot_chr**: Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **annot_start**: Numeric vector specifying the starting position of annotation entries.
  * **annot_end**: Numeric vector specifying the ending position of annotation entries.
  * **annot_sub_ids**: Character vector defining sets of sub-IDs for annotation entries (separated by comma).
  * **annot_sub_score**: Character vector defining sets of sub-scores for annotation entries (separated by comma). 
  
   __outputs__:
   * **input_ID**: Character vector specifying the operlapping annotation ids of output (separated by comma).
    * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
   * **overlap_sub_ids**: The sub-IDs of annotation entries that overlap with each input region (separated by comma).
   * **overlap_sub_scores**: CThe average score of annotation for each overlapping sub-ID over each input region (separated by comma).
   * **overlap_sub_percentages**: The percentage of overlap between each input region and its overlapping sub-IDs (separated by comma).
   
   __Note__:
   * For each input genomic region, the average score of annotation for each sub-ID is measured as the weighted average of scores of overlapping annotation entries where the percentages of overlaps of annotation entries with the input region are taken as weights (this is measured for each sub-ID separately). The percentage of overlap between an input genomic region with a specific annotation sub-ID is calculated as the total portion of the input region that is covered by that annotation sub-ID divided by the length of the region.


   
There is a 'dnase.tsv' file for testing.
   
```r
dnase <- fread("DIR/test_data/dnase.tsv")

# Then you run the annotateBinary function
result <- annotateMultiScore(input_id = table$geneID, input_chr = table$chr, input_start = table$start, input_end = table$end, annot_chr = dnase$chrom, annot_start = dnase$chromStart, annot_end = dnase$chromEnd, annot_sub_id = dnase$sourceIds, annot_sub_score = dnase$sourceScores)


```

The result is like below

| GeneID  | chr | start | end | overlapping_sub_ids | overlap_sub_scores | overlap_percentages |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| ENSG00000223764	  | chr1  | 845735 | 864616 | 90,115,89 | 433.94021101993,560.202813599062,315 | 7.34838042729152,7.34838042729152,2.42935906271537 |





* `annotateSingleScore(input_id, input_chr, input_start, input_end, annot_chr, annot_start, annot_end, annot_score)`: This function annotates the input genomic regions with a given "single-score" annotation. A "single-score" annotation specifies a set of genomic regions and assigns a numeric score to each of them (e.g. the annotation of histone modification peaks by ENCODE).


                
  __inputs__:
  * **input_id**: Character vector defining the name of input genomic regions (e.g. gene id)
  * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
  * **annot_chr**: Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **annot_start**: Numeric vector specifying the starting position of annotation entries.
  * **annot_end**: Numeric vector specifying the ending position of annotation entries.
  * **annot_score**: Numeric vector specifying scores for annotation entries (separated by comma). 
  
   __outputs__:
   * **input_ID**: Character vector specifying the operlapping annotation ids of output (separated by comma).
    * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
   * **overlap_score**: The average score of annotation over each input region.
   * **overlap_percentage**: The percentage of overlap between each input region and annotation entries.
  
  __Note__:
  * For each input genomic region, the average score of annotation is measured as the weighted average of scores of overlapping annotation entries where the percentages of overlaps of annotation entries with the input region are taken as weights. The percentage of overlap between an input genomic region with annotation entries is calculated as the total portion of the input region that is covered by the annotation entries divided by the length of the region.
  
  

There is a 'ctcf.tsv' file for testing.
   
```r
ctcf <- fread("DIR/test_data/ctcf.tsv")

# Then you run the annotateBinary function

result <- annotateSingleScore(input_id = table$geneID, input_chr = table$chr, input_start = table$start, input_end = table$end, annot_chr = ctcf$chr, annot_start = ctcf$start, annot_end = ctcf$end, annot_score = ctcf$score)


```

The result is like below

| GeneID  | chr | start | end | overlap_score | overlap_percentages |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| ENSG00000223764	  | chr1  | 845735 | 864616 | 435.9652 | 0.2223493 |


  
* `annotateCategorical(input_id, input_chr, input_start, input_end, annot_chr, annot_start, annot_end, annot_category)`: This function annotates the input genomic regions with a given "categorical" annotation. A "categorical" annotation specifies a set of genomic regions and assigns a category to each of them (e.g. chromHMM annotation which assigns a chromatin state to any segment of the genome).
                
  __inputs__:
  * **input_id**: Character vector defining the name of input genomic regions (e.g. gene id)
  * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
  * **annot_chr**: Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **annot_start**: Numeric vector specifying the starting position of annotation entries.
  * **annot_end**: Numeric vector specifying the ending position of annotation entries.
  * **annot_category**: Character vector specifying categories for annotation entries.
  
   __outputs__:
   * **input_ID**: Character vector specifying the operlapping annotation ids of output (separated by comma).
    * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
   * **overlapping_categories**: The categories that overlap with each input region (separated by comma).
   * **overlapping_percentage**: The percentages of overlaps between each input region and its overlapping categories (separated by comma).
   
   __Note__:
   * The percentage of overlap between an input genomic region with a specific category is calculated as the total portion of the input region that is covered by that category devided by the length of the region.


There is a 'hmm.tsv' file for testing.

```r

hmm <- fread("DIR/test_data/hmm.tsv")

# Then you run the annotateBinary function

result <- annotateCategorical(input_id = table$geneID, input_chr = table$chr, input_start = table$start, input_end = table$end, annot_chr = hmm$chr, annot_start = hmm$start, annot_end = hmm$end, annot_category = hmm$name)

```

The result is like below

| GeneID  | chr | start | end | overlapping_categories | overlap_percentages |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| ENSG00000223764	  | chr1  | 845735 | 864616 | 4_Strong_Enhancer | 100 |



   

* `findHiCinteractions(input_id, input_chr, input_start, input_end, hic_f1_id, hic_f1_chr, hic_f1_start, hic_f1_end, hic_f2_id, hic_f2_chr, hic_f2_start, hic_f2_end, target_id, target_chr, target_start, target_end)`: This function finds interactions between input and target genomic regions based on a given set of Hi-C interactions. It basically finds the subset of given Hi-C interactions in which one fragment overlaps with input genomic regions and the other fragment overlaps with target genomic regions.
 
                
  __inputs__:
  * **input_id**: Character vector defining the name of input genomic regions (e.g. gene id)
  * **input_chr**: Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **input_start**: Numeric vector specifying the starting position of input genomic regions.
  * **input_end**: Numeric vector specifying the ending position of input genomic regions.
  * **hic_f1_id**: Character defining the ID for left fragments of Hi-C interactions.
  * **hic_f1_chr**: Character vector defining the name of the chromosome for left fragments of Hi-C interactions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **hic_f1_start**: Numeric vector specifying the starting position of left fragments of Hi-C interactions.
  * **hic_f1_end**: Numeric vector specifying the ending position of left fragments of Hi-C interactions.
  * **hic_f2_id**: Character defining the ID for left fragments of Hi-C interactions.
  * **hic_f2_chr**: Character vector defining the name of the chromosome for left fragments of Hi-C interactions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **hic_f2_start**: Numeric vector specifying the starting position of left fragments of Hi-C interactions.
  * **hic_f2_end**: Numeric vector specifying the ending position of left fragments of Hi-C interactions.
  * **target_id**: Character vector defining the name of target genomic regions (e.g. gene id)
  * **target_chr**: Character vector defining the name of the chromosome for target genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
  * **target_start**: Numeric vector specifying the starting position of target genomic regions.
  * **target_end**: Numeric vector specifying the ending position of target genomic regions.
  
   __outputs__:
   A list of two data.frames:

  The first data.frame contains the input genomic regions and 3 extra columns which specify:

  1. The number of target genomic regions that interact with each input genomic region through Hi-C interactions.

  2. The IDs of target genomic regions that interact with each input genomic region through Hi-C interactions (separated by comma).

  3. The IDs of Hi-C interactions that connect each input genomic region with target genomic regions (separated by comma).

  The second data.frame provides further information about each Hi-C interaction. It specifies the genomic regions of the left and right Hi-C fragments. Furthermore, the IDs of input genomic regions and target genomic regions that interact through these Hi-C interactions are listed in each row (separated by comma).
  
  __Note__:
  * The ID of a Hi-C interaction is defined by concatenating the left fragment ID and the right fragment ID (separated by dash).
  

* `testRegions(ROIs, mutational_catalogs)`: Given a set of input genomic regions and a list of mutational catalogues in various cancers, a p-value for each genomic region (in each cancer separately) is calculated which can be used as a measure of how significantly each genomic region is mutated in each cancer (see [SomaGene paper]() for more details).


  __inputs__:
  
    * **ROIs**:	A data.frame containing input genomic regions. The required columns are: region_id - Defines the name of the genomic regions (e.g. gene id).
      * chr - The name of the chromosome (one of chr1, chr2, ..., chrX, chrY or chrM).
      * start - The starting position of the region in the chromosome.
      * end - The ending position of the region in the chromosome.

    * **mutational_catalogs**:	A list of data.frames. Each data.frame contains the mutational catalogue in a specific cancer and must have the following columns: 
      * subject_id - Defines the id of the subject for which the mutation is recorded (e.g. ICGC donor-id or sample-id).
      * chr - The name of the chromosome (one of chr1, chr2, ..., chrX, chrY or chrM).
      * pos - The position of the mutation in the chromosome.
      * cancer - The cancer type of the mutational catalogue.
      
  __outputs__:
  * A data.frame containing the input genomic regions, plus the calculated p-values for each genomic region and each cancer.

