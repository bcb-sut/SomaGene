#' Annotate the input genomic regions with annotation of type "categorical".
#'
#' This function annotates the input genomic regions with a given "categorical" annotation. A "categorical" annotation specifies
#' a set of genomic regions and assigns a category to each of them (e.g. chromHMM annotation which assigns a chromatin state to
#' any segment of the genome).
#'
#' @param input_id Character vector defining the name of input genomic regions (e.g. gene id)
#' @param input_chr Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param input_start Numeric vector specifying the starting position of input genomic regions.
#' @param input_end Numeric vector specifying the ending position of input genomic regions.
#'
#' @param annot_chr Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param annot_start Numeric vector specifying the starting position of annotation entries.
#' @param annot_end Numeric vector specifying the ending position of annotation entries.
#' @param annot_category Character vector specifying categories for annotation entries.
#'
#' @return A data.frame containing the input genomic regions and 2 extra columns which specify:
#' \enumerate{
#'   \item The categories that overlap with each input region (separated by comma).
#'   \item The percentages of overlaps between each input region and its overlapping categories (separated by comma).
#' }
#'
#' @note The percentage of overlap between an input genomic region with a specific category is calculated as the
#' total portion of the input region that is covered by that category devided by the length of the region.
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table as.data.table
#'
#' @export
annotateCategorical <- function(
    input_id,
    input_chr,
    input_start,
    input_end,

    annot_chr,
    annot_start,
    annot_end,
    annot_category) {



    if(length(setdiff( unique(input_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for input regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(annot_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for annotation regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }

    temp <- length(unique(length(input_id),length(input_chr),length(input_start),length(input_end)))
    if(temp == 0){ stop('NULL vectors for input regions') }
    if(temp > 1){ stop('the lengths of vectors for input regions are unequal') }

    temp <- length(unique(length(annot_category),length(annot_chr),length(annot_start),length(annot_end)))
    if(temp == 0){ stop('NULL vectors for annotation regions') }
    if(temp > 1){ stop('the lengths of vectors for annotation regions are unequal') }


    cat('step 1/3','\n')

    input_ranges <- GRanges(seqnames = input_chr, ranges = IRanges(input_start, input_end))
    annot_ranges <- GRanges(seqnames = annot_chr, ranges = IRanges(annot_start, annot_end))

    ovlp_hits_input_annot <- suppressWarnings(findOverlaps(input_ranges, annot_ranges))
    ovlp_input_annot <- as.data.table(ovlp_hits_input_annot)
    colnames(ovlp_input_annot) <- c('inputHits', 'annotHits')

    if(length(annot_category) != length(annot_chr)){
        stop('err: the lengths of vectors for annotation categories and annotation entries do not match')

    }
    widths_of_ovlps <- width(overlapsRanges(IRanges(input_start, input_end),
                                            IRanges(annot_start, annot_end),
                                            ovlp_hits_input_annot))
    ovlp_input_annot$ovlp_width <- widths_of_ovlps

    ovlp_input_annot$category <- annot_category[ovlp_input_annot$annotHits]
    ovlp_input_annot <- ovlp_input_annot[,
                                         list(ovlp_width = sum(ovlp_width)),
                                         by = list(inputHits, category)]

    ovlp_input_annot$input_width <- width(input_ranges[ovlp_input_annot$inputHits])

    ovlp_input_annot$ovlp_prct <- 100 * (ovlp_input_annot$ovlp_width / ovlp_input_annot$input_width)
    ovlp_input_annot$ovlp_prct[which(ovlp_input_annot$ovlp_prct > 100)] <- 100

    cat('step 2/3','\n')

    ovlp_agg <- ovlp_input_annot[,
                                 list(overlapping_categories = paste0(category, collapse = ','),
                                      ovlp_prcts = paste0(ovlp_prct, collapse = ',')),
                                 by = inputHits]

    output <- data.table(input_ID = input_id, input_chr = input_chr, input_start = input_start, input_end = input_end)
    output$overlapping_categories <- ''
    output$overlapping_categories[ovlp_agg$inputHits] <- ovlp_agg$overlapping_categories
    output$overlap_percentages <- ''
    output$overlap_percentages[ovlp_agg$inputHits] <- ovlp_agg$ovlp_prcts

    cat('step 3/3','\n')

    return(output)
}



