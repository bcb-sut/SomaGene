#' Annotate the input genomic regions with annotation of type "binary".
#'
#' This function annotates the input genomic regions with a given "binary" annotation. A "binary" annotation is simply
#' a set of genomic regions without any extra attribute (e.g. the set of enhancers/promoters or the set of eQTLs).
#'
#' @usage
#' annotateBinary(input_id, input_chr, input_start, input_end,
#'                 annot_id = NULL, annot_chr, annot_start, annot_end)
#'
#' @param input_id Character vector defining the name of input genomic regions (e.g. gene id)
#' @param input_chr Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param input_start Numeric vector specifying the starting position of input genomic regions.
#' @param input_end Numeric vector specifying the ending position of input genomic regions.
#'
#' @param annot_id Optional (NULL by default): Character vector defining the name or ID of each annotation entry. If NULL, the row numbers (indexes) are used as IDs.
#' @param annot_chr Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param annot_start Numeric vector specifying the starting position of annotation entries.
#' @param annot_end Numeric vector specifying the ending position of annotation entries.
#'
#' @return A data.frame containing the input genomic regions and 2 extra columns which specify:
#' \enumerate{
#'   \item The number of annotation entries that overlap with each input region.
#'   \item The IDs (names) of annotation entries that overlap with each input region (separated by comma).
#' }
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table as.data.table
#'
#' @export
annotateBinary <- function(
    input_id,
    input_chr,
    input_start,
    input_end,

    annot_id = NULL,
    annot_chr,
    annot_start,
    annot_end) {



    if(length(setdiff( unique(input_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for input regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(annot_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for annotation regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }

    temp <- length(unique(length(input_id),length(input_chr),length(input_start),length(input_end)))
    if(temp == 0){ stop('NULL vectors for input regions') }
    if(temp > 1){ stop('the lengths of vectors for input regions are unequal') }

    temp <- length(unique(length(annot_chr),length(annot_start),length(annot_end)))
    if(temp == 0){ stop('NULL vectors for annotation regions') }
    if(temp > 1){ stop('the lengths of vectors for annotation regions are unequal') }


    cat('step 1/2','\n')

    input_ranges <- GRanges(seqnames = input_chr, ranges = IRanges(input_start, input_end))
    annot_ranges <- GRanges(seqnames = annot_chr, ranges = IRanges(annot_start, annot_end))

    ovlp_hits_input_annot <- suppressWarnings(findOverlaps(input_ranges, annot_ranges))
    ovlp_input_annot <- as.data.table(ovlp_hits_input_annot)
    colnames(ovlp_input_annot) <- c('inputHits', 'annotHits')


    if(is.null(annot_id)){
        annot_id <- as.character(1:length(annot_chr))
    }

    ovlp_input_annot$ovlp_ID <- annot_id[ovlp_input_annot$annotHits]
    ovlp_agg <- ovlp_input_annot[,
                                 list(overlap_count = length(annotHits),
                                      overlapping_annot_IDs = paste(ovlp_ID, collapse = ',')),
                                 by = inputHits]

    output <- data.table(input_ID = input_id, input_chr = input_chr, input_start = input_start, input_end = input_end)
    output$overlap_count <- 0
    output$overlap_count[ovlp_agg$inputHits] <- ovlp_agg$overlap_count
    output$overlapping_annot_IDs <- ''
    output$overlapping_annot_IDs[ovlp_agg$inputHits] <- ovlp_agg$overlapping_annot_IDs

    cat('step 2/2','\n')

    return(output)
}



