#' Annotate the input genomic regions with annotation of type "single-score".
#'
#' This function annotates the input genomic regions with a given "single-score" annotation. A "single-score" annotation specifies
#' a set of genomic regions and assigns a numeric score to each of them (e.g. the annotation of histone modification peaks by ENCODE).
#'
#' @param input_id Character vector defining the name of input genomic regions (e.g. gene id)
#' @param input_chr Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param input_start Numeric vector specifying the starting position of input genomic regions.
#' @param input_end Numeric vector specifying the ending position of input genomic regions.
#'
#' @param annot_chr Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param annot_start Numeric vector specifying the starting position of annotation entries.
#' @param annot_end Numeric vector specifying the ending position of annotation entries.
#' @param annot_score Numeric vector specifying scores for annotation entries.
#'
#' @return A data.frame containing the input genomic regions and 2 extra columns which specify:
#' \enumerate{
#'   \item The average score of annotation over each input region.
#'   \item The percentage of overlap between each input region and annotation entries.
#' }
#'
#' @note For each input genomic region, the average score of annotation is measured as the weighted average of scores of
#' overlapping annotation entries where the percentages of overlaps of annotation entries with the input region are taken
#' as weights.
#' The percentage of overlap between an input genomic region with annotation entries is calculated as the
#' total portion of the input region that is covered by the annotation entries divided by the length of the region.
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table as.data.table
#'
#' @export
annotateSingleScore <- function(
    input_id,
    input_chr,
    input_start,
    input_end,

    annot_chr,
    annot_start,
    annot_end,
    annot_score) {



    if(length(setdiff( unique(input_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for input regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(annot_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for annotation regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }

    temp <- length(unique(length(input_id),length(input_chr),length(input_start),length(input_end)))
    if(temp == 0){ stop('NULL vectors for input regions') }
    if(temp > 1){ stop('the lengths of vectors for input regions are unequal') }

    temp <- length(unique(length(annot_chr),length(annot_start),length(annot_end),length(annot_score)))
    if(temp == 0){ stop('NULL vectors for annotation regions') }
    if(temp > 1){ stop('the lengths of vectors for annotation regions are unequal') }


    cat('Resolving...')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Resolve the overlapping ranges of the annotation  - - - - - - - -
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    annot_ranges <- GRanges(seqnames = annot_chr, ranges = IRanges(annot_start, annot_end))

    annot_disjoin_ranges <- disjoin(annot_ranges, with.revmap = T)

    resolved_annot <- data.table(chr = as.character(seqnames(annot_disjoin_ranges)),
                                 start = start(annot_disjoin_ranges),
                                 end = end(annot_disjoin_ranges))

    revmap <- mcols(annot_disjoin_ranges)$revmap

    rm(annot_disjoin_ranges)

    if( any(sapply(revmap, function(s) length(s)) > 1) ){

        cat('found overlapping ranges...')

        reps_of_revmap <- sapply(revmap, length)

        unfolded_revmap <- unlist(revmap)

        rm(revmap)

        by_indicator <- rep(x = 1:length(reps_of_revmap), times = reps_of_revmap)

        rm(reps_of_revmap)

        unfolded_annot_metadata <- data.table(signalValue = annot_score[unfolded_revmap])
        unfolded_annot_metadata$by_indicator <- by_indicator

        rm(unfolded_revmap)

        aggregated_annot_metadata <- unfolded_annot_metadata[, list(signalValue = mean(signalValue)),
                                                             by = by_indicator]

        rm(unfolded_annot_metadata)
        rm(by_indicator)

        resolved_annot <- cbind(resolved_annot, aggregated_annot_metadata)

        rm(aggregated_annot_metadata)

        annot_chr <- resolved_annot$chr
        annot_start <- resolved_annot$start
        annot_end <- resolved_annot$end
        annot_score <- resolved_annot$signalValue

        rm(resolved_annot)
    } else {
        rm(revmap)
        rm(resolved_annot)
    }
    cat('Ok\n')


    cat(1,'\n')

    input_ranges <- GRanges(seqnames = input_chr, ranges = IRanges(input_start, input_end))
    annot_ranges <- GRanges(seqnames = annot_chr, ranges = IRanges(annot_start, annot_end))

    ovlp_hits_input_annot <- suppressWarnings(findOverlaps(input_ranges, annot_ranges))
    ovlp_input_annot <- as.data.table(ovlp_hits_input_annot)
    colnames(ovlp_input_annot) <- c('inputHits', 'annotHits')


    widths_of_ovlps <- width(overlapsRanges(IRanges(input_start, input_end),
                                            IRanges(annot_start, annot_end),
                                            ovlp_hits_input_annot))
    ovlp_input_annot$ovlp_width <- widths_of_ovlps

    ovlp_input_annot$score <- annot_score[ovlp_input_annot$annotHits]

    ovlp_input_annot$input_width <- width(input_ranges[ovlp_input_annot$inputHits])

    ovlp_input_annot$ovlp_prct <- 100 * (ovlp_input_annot$ovlp_width / ovlp_input_annot$input_width)
    ovlp_input_annot$ovlp_prct[which(ovlp_input_annot$ovlp_prct > 100)] <- 100

    cat(2,'\n')

    ovlp_agg <- ovlp_input_annot[,
                                 list(score = sum(as.numeric(score)*ovlp_prct) / sum(ovlp_prct),
                                      ovlp_prct = sum(ovlp_prct)),
                                 by = inputHits]
    ovlp_agg$ovlp_prct[which(ovlp_agg$ovlp_prct > 100)] <- 100


    output <- data.table(input_ID = input_id, input_chr = input_chr, input_start = input_start, input_end = input_end)
    output$overlap_score <- ''
    output$overlap_score[ovlp_agg$inputHits] <- ovlp_agg$score
    output$overlap_percentage <- ''
    output$overlap_percentage[ovlp_agg$inputHits] <- ovlp_agg$ovlp_prct

    cat(3,'\n')

    return(output)
}



