#' Apply a statistical test to identify significantly mutated genomic regions in each cancer.
#'
#' Given a set of input genomic regions and a list of mutational catalogues in various cancers,
#' a p-value for each genomic region (in each cancer separately) is calculated which can be used
#' as a measure of how significantly each genomic region is mutated in each cancer
#' (see References for more details).
#'
#' @param ROIs A data.frame containing input genomic regions. The required columns are: \cr
#' \itemize{
#'   \item region_id - Defines the name of the genomic regions (e.g. gene id).
#'   \item chr - The name of the chromosome (one of chr1, chr2, ..., chrX, chrY or chrM).
#'   \item start - The starting position of the region in the chromosome.
#'   \item end - The ending position of the region in the chromosome.
#' }
#' @param mutational_catalogs A list of data.frames. Each data.frame contains the mutational catalogue in a specific cancer
#' and must have the following columns: \cr
#' \itemize{
#'   \item subject_id - Defines the id of the subject for which the mutation is recorded (e.g. ICGC donor-id or sample-id).
#'   \item chr - The name of the chromosome (one of chr1, chr2, ..., chrX, chrY or chrM).
#'   \item pos - The position of the mutation in the chromosome.
#'   \item cancer - The cancer type of the mutational catalogue.
#' }
#'
#' @return A data.frame containing the input genomic regions, plus the calculated p-values for each genomic region and each cancer.
#'
#' @references asdfas
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table
#'
#' @export
testRegions <- function(ROIs, mutational_catalogs){

    if( ! inherits(ROIs,'data.frame')){
        stop('ROIs must be a data.frame')
    }
    if(class(mutational_catalogs) != 'list'){
        stop('mutational_catalogs must be a list of data.frames')
    }

    # ROIs <- '../../data/data_for_pipeline/ROIs_for_pipeline.tsv'

    # check for required columns in ROIs
    if(!all(c('region_id','chr','start','end') %in% colnames(ROIs))){
        stop('ROIs must contain the columns {region_id, chr, start, end}')
    }
    if(length(setdiff( unique(ROIs$chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('the values of "chr" column in ROIs file must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }




    # mutations_data_path_list <- list.files('../../data/mutations_data',full.names = T)
    # for(p in mutations_data_path_list){
    #     s <- strsplit(strsplit(p,'/')[[1]][5],'_')[[1]]
    #     ss <- s[-c(length(s)-1,length(s))]
    #     cancer <- paste(ss,collapse = '-')
    #     cat(cancer,'\n')
    #     a <- readRDS(p)
    #     b <- data.table(subject_id = a$donor_id, chr = a$chr, pos = a$pos, cancer = cancer)
    #     fwrite(b, paste0('../../data/data_for_pipeline/mutations_data_for_pipeline/',cancer,'_mutations_data.tsv'), sep = '\t', col.names = T, row.names = F)
    # }

    # mutations_data_path_list <- list.files('../../data/data_for_pipeline/mutations_data_for_pipeline', full.names = T)

    cancers <- c()

    for(i in 1:length(mutational_catalogs)){

        cat(paste0('processing mutational catalog ',i,' (out of ',length(mutational_catalogs),')\n'))

        mutations_loci_samples <- mutational_catalogs[[i]]

        if( ! inherits(mutations_loci_samples,'data.frame')){
            stop(paste0('mutational_catalogs[[',i,']] must be a data.frame (mutational_catalogs must be a list of data.frames)'))
        }

        if(!all(c('subject_id','chr','pos','cancer') %in% colnames(mutations_loci_samples))){
            stop(paste0('mutational_catalogs[[',i,']] must contain the columns {subject_id, chr, pos, cancer}'))
        }
        if(length(setdiff( unique(mutations_loci_samples$chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
            stop(paste0('The values of "chr" column in mutational_catalogs[[',i,']] must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}'))
        }
        cancer <- unique(mutations_loci_samples$cancer)
        if(length(cancer) != 1){
            stop(paste0('The "cancer" column in mutational_catalogs[[',i,']] must have unique value'))
        }

        cancers <- c(cancers, cancer)

        mutations_loci_samples$pos <- as.numeric(mutations_loci_samples$pos)

        cat('counting mutated samples for each ROI for',cancer,'cancer\n')
        ovlp_mut_samp <-
            findOverlaps(
                GRanges(
                    seqnames=ROIs$chr,
                    ranges=IRanges(ROIs$start,ROIs$end)
                ),
                GRanges(
                    seqnames=mutations_loci_samples$chr,
                    ranges=IRanges(mutations_loci_samples$pos,mutations_loci_samples$pos)
                )
            )
        ovlp_samp_dt <- data.table(roi_indx = queryHits(ovlp_mut_samp),
                                               ovlp_samp = mutations_loci_samples$subject_id[subjectHits(ovlp_mut_samp)])

        aggregated_ovlp_samp_dt <-  ovlp_samp_dt[,length(unique(ovlp_samp)), by=roi_indx]

        in_cancer_mutated <- paste0('in_',cancer,'_mutated')
        ROIs[[in_cancer_mutated]] <- 0
        ROIs[[in_cancer_mutated]][aggregated_ovlp_samp_dt$roi_indx] <- aggregated_ovlp_samp_dt$V1

        in_cancer_not_mutated <- paste0('in_',cancer,'_not_mutated')
        ROIs[[in_cancer_not_mutated]] <- length(unique(mutations_loci_samples$subject_id)) - ROIs[[in_cancer_mutated]]
    }

    cat('creating remaining columns\n')
    for(this_cancer in cancers){
        cols <- paste0('in_',setdiff(cancers,this_cancer),'_mutated')
        this_col <- paste0('not_in_',this_cancer,'_mutated')
        ROIs[[this_col]] <- rowSums(ROIs[,..cols])

        cols <- paste0('in_',setdiff(cancers,this_cancer),'_not_mutated')
        this_col <- paste0('not_in_',this_cancer,'_not_mutated')
        ROIs[[this_col]] <- rowSums(ROIs[,..cols])

        mutated_col <- paste0('in_',this_cancer,'_mutated')
        not_mutated_col <- paste0('in_',this_cancer,'_not_mutated')

        cancer_total_samples_col <- paste0(this_cancer,'_total_samples')
        ROIs[[cancer_total_samples_col]] <- ROIs[[mutated_col]] + ROIs[[not_mutated_col]]
    }

    for(this_cancer in cancers){
        cat('performing fisher exact test for',this_cancer,'cancer\n')
        cols <- c(paste0('in_',this_cancer,'_mutated'),paste0('in_',this_cancer,'_not_mutated'),
                  paste0('not_in_',this_cancer,'_mutated'),paste0('not_in_',this_cancer,'_not_mutated'))
        this_ROIs <- ROIs[,..cols]
        this_col <- paste0('pvalue_in_',this_cancer)
        ROIs[[this_col]] <- 10
        for(i in 1:dim(this_ROIs)[1]){
            ROIs[[this_col]][i] <- fisher.test(matrix(unlist(this_ROIs[i,]),2,2), alternative="greater")$p.value
        }
    }

    ROIs <- ROIs[,c(1:4,grep('pvalue',colnames(ROIs))), with = F]

    return(ROIs)
}


