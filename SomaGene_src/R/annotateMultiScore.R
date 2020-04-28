#' Annotate the input genomic regions with annotation of type "multi-score".
#'
#' This function annotates the input genomic regions with a given "multi-score" annotation. A "multi-score" annotation specifies
#' a set of genomic regions and assigns a set of sub-IDs and their corresponding sub-scores to each of them
#' (e.g. the annotation of DNAse hypersensitive clusters by ENCODE).
#'
#' @param input_id Character vector defining the name of input genomic regions (e.g. gene id)
#' @param input_chr Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param input_start Numeric vector specifying the starting position of input genomic regions.
#' @param input_end Numeric vector specifying the ending position of input genomic regions.
#'
#' @param annot_chr Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param annot_start Numeric vector specifying the starting position of annotation entries.
#' @param annot_end Numeric vector specifying the ending position of annotation entries.
#' @param annot_sub_id Character vector defining sets of sub-IDs for annotation entries (separated by comma).
#' @param annot_sub_score Character vector defining sets of sub-scores for annotation entries (separated by comma).
#'
#' @return A data.frame containing the input genomic regions and 3 extra columns which specify:
#' \enumerate{
#'   \item The sub-IDs of annotation entries that overlap with each input region (separated by comma).
#'   \item The average score of annotation for each overlapping sub-ID over each input region.
#'   \item The percentage of overlap between each input region and its overlapping sub-IDs.
#' }
#'
#' @note For each input genomic region, the average score of annotation for each sub-ID is measured as the weighted average of scores of
#' overlapping annotation entries where the percentages of overlaps of annotation entries with the input region are taken
#' as weights (this is measured for each sub-ID separately).
#' The percentage of overlap between an input genomic region with a specific annotation sub-ID is calculated as the
#' total portion of the input region that is covered by that annotation sub-ID divided by the length of the region.
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table as.data.table
#'
#' @export
annotateMultiScore <- function(
    input_id,
    input_chr,
    input_start,
    input_end,

    annot_chr,
    annot_start,
    annot_end,
    annot_sub_id,
    annot_sub_score) {



    if(length(setdiff( unique(input_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for input regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(annot_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for annotation regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }

    temp <- length(unique(length(input_id),length(input_chr),length(input_start),length(input_end)))
    if(temp == 0){ stop('NULL vectors for input regions') }
    if(temp > 1){ stop('the lengths of vectors for input regions are unequal') }

    temp <- length(unique(length(annot_chr),length(annot_start),length(annot_end),length(annot_sub_id),length(annot_sub_score)))
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

        unfolded_annot_metadata <- data.table(sub_id = annot_sub_id[unfolded_revmap], sub_score = annot_sub_score[unfolded_revmap])
        unfolded_annot_metadata$by_indicator <- by_indicator

        rm(unfolded_revmap)

        aggregated_annot_metadata <- unfolded_annot_metadata[, list(sub_id = paste0(sub_id, collapse = ','),
                                                                   sub_score = paste0(sub_score, collapse = ',')),
                                                            by = by_indicator]
        rm(unfolded_annot_metadata)
        rm(by_indicator)

        resolved_annot <- cbind(resolved_annot, aggregated_annot_metadata)

        rm(aggregated_annot_metadata)

        annot_chr <- resolved_annot$chr
        annot_start <- resolved_annot$start
        annot_end <- resolved_annot$end
        annot_sub_id <- resolved_annot$sub_id
        annot_sub_score <- resolved_annot$sub_score

        rm(resolved_annot)
    } else {
        rm(revmap)
        rm(resolved_annot)
    }
    cat('Ok\n')


    input_ranges <- GRanges(seqnames = input_chr, ranges = IRanges(input_start, input_end))
    annot_ranges <- GRanges(seqnames = annot_chr, ranges = IRanges(annot_start, annot_end))

    ovlp_hits_input_annot <- suppressWarnings(findOverlaps(input_ranges, annot_ranges))
    ovlp_input_annot <- as.data.table(ovlp_hits_input_annot)
    colnames(ovlp_input_annot) <- c('inputHits', 'annotHits')



    cat(1,'\n')

    widths_of_ovlps <- width(overlapsRanges(IRanges(input_start, input_end),
                                            IRanges(annot_start, annot_end),
                                            ovlp_hits_input_annot))
    ovlp_input_annot$ovlp_width <- widths_of_ovlps

    ovlp_input_annot$sub_id <- annot_sub_id[ovlp_input_annot$annotHits]
    ovlp_input_annot$sub_score <- annot_sub_score[ovlp_input_annot$annotHits]
    ovlp_input_annot$num_sub_ids <- sapply(strsplit(annot_sub_id, ','), length)[ovlp_input_annot$annotHits]

    ovlp_input_annot$input_width <- width(input_ranges[ovlp_input_annot$inputHits])

    ovlp_input_annot$ovlp_prct <- 100 * (ovlp_input_annot$ovlp_width / ovlp_input_annot$input_width)
    ovlp_input_annot$ovlp_prct[which(ovlp_input_annot$ovlp_prct > 100)] <- 100


    cat(2,'\n')
    ovlp_input_annot$ovlp_prct_str <- lapply(c(1:dim(ovlp_input_annot)[1]),
                                             function(i) {
                                                 paste(rep(ovlp_input_annot$ovlp_prct[i], ovlp_input_annot$num_sub_ids[i]),collapse = ',')
                                             })

    cat(3,'\n')
    ovlp_input_annot$annotHits_str <- lapply(c(1:dim(ovlp_input_annot)[1]),
                                             function(i) {
                                                 paste(rep(ovlp_input_annot$annotHits[i], ovlp_input_annot$num_sub_ids[i]),collapse = ',')
                                             })

    cat(4,'\n')
    ovlp_agg <- ovlp_input_annot[,
                                 list(sub_id = paste0(sub_id, collapse = ','),
                                      sub_score = paste0(sub_score, collapse = ','),
                                      ovlp_prct_str = paste0(ovlp_prct_str, collapse = ','),
                                      annotHits_str = paste0(annotHits_str, collapse = ',')),
                                 by = inputHits]

    cat(5,'\n')
    for(row_indx in 1:dim(ovlp_agg)[1]) {
        dt <- data.table(sub_id = unlist(strsplit(ovlp_agg[row_indx]$sub_id,',')),
                         sub_score = as.numeric(unlist(strsplit(ovlp_agg[row_indx]$sub_score,','))),
                         ovlp_prct = as.numeric(unlist(strsplit(ovlp_agg[row_indx]$ovlp_prct_str,','))),
                         annotHits = as.numeric(unlist(strsplit(ovlp_agg[row_indx]$annotHits_str,','))))

        # because in some cases, there are duplicated IDs in a single row of annot data table
        dt <- dt[, {sub_score = mean(sub_score); ovlp_prct = mean(ovlp_prct); list(sub_id=sub_id,
                                                                                   sub_score=sub_score,
                                                                                   ovlp_prct=ovlp_prct,
                                                                                   annotHits=annotHits)}, by = list(sub_id,annotHits)]

        # to resolve multiple scores and overlaps existing for a single ID when merging the overlapping intervals
        dt <- dt[, list(sub_score = sum(sub_score*ovlp_prct) / sum(ovlp_prct), ovlp_prct = sum(ovlp_prct)), by = sub_id]

        ovlp_agg[row_indx]$sub_id <- paste(dt$sub_id,collapse=',')
        ovlp_agg[row_indx]$sub_score <- paste(dt$sub_score,collapse=',')
        ovlp_agg[row_indx]$ovlp_prct_str <- paste(dt$ovlp_prct,collapse=',')
    }

    cat(6,'\n')
    output <- data.table(input_ID = input_id, input_chr = input_chr, input_start = input_start, input_end = input_end)
    output$overlapping_sub_ids <- ''
    output$overlapping_sub_ids[ovlp_agg$inputHits] <- ovlp_agg$sub_id
    output$overlap_sub_scores <- ''
    output$overlap_sub_scores[ovlp_agg$inputHits] <- ovlp_agg$sub_score
    output$overlap_percentages <- ''
    output$overlap_percentages[ovlp_agg$inputHits] <- ovlp_agg$ovlp_prct_str

    return(output)
}



