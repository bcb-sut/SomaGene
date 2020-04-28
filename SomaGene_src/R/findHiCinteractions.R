#' Find interactions between two sets of genomic regions through Hi-C interactions.
#'
#' This function finds interactions between input and target genomic regions based on a given set of Hi-C interactions.
#' It basically finds the subset of given Hi-C interactions in which one fragment overlaps with input genomic regions
#' and the other fragment overlaps with target genomic regions.
#'
#' @param input_id Character vector defining the name of input genomic regions (e.g. gene id)
#' @param input_chr Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param input_start Numeric vector specifying the starting position of input genomic regions.
#' @param input_end Numeric vector specifying the ending position of input genomic regions.
#'
#' @param hic_f1_id Character defining the ID for left fragments of Hi-C interactions.
#' @param hic_f1_chr Character vector defining the name of the chromosome for left fragments of Hi-C interactions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param hic_f1_start Numeric vector specifying the starting position of left fragments of Hi-C interactions.
#' @param hic_f1_end Numeric vector specifying the ending position of left fragments of Hi-C interactions.
#'
#' @param hic_f2_id Character defining the ID for right fragments of Hi-C interactions.
#' @param hic_f2_chr Character vector defining the name of the chromosome for right fragments of Hi-C interactions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param hic_f2_start Numeric vector specifying the starting position of right fragments of Hi-C interactions.
#' @param hic_f2_end Numeric vector specifying the ending position of right fragments of Hi-C interactions.
#'
#' @param target_id Character vector defining the name of target genomic regions (e.g. gene id)
#' @param target_chr Character vector defining the name of the chromosome for target genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).
#' @param target_start Numeric vector specifying the starting position of target genomic regions.
#' @param target_end Numeric vector specifying the ending position of target genomic regions.
#'
#' @return A list of two data.frames:
#' \itemize{
#'   \item {
#'     The first data.frame contains the input genomic regions and 3 extra columns which specify:
#'     \enumerate{
#'       \item The number of target genomic regions that interact with each input genomic region through Hi-C interactions.
#'       \item The IDs of target genomic regions that interact with each input genomic region through Hi-C interactions (separated by comma).
#'       \item The IDs of Hi-C interactions that connect each input genomic region with target genomic regions (separated by comma).
#'     }
#'   }
#'   \item {
#'     The second data.frame provides further information about each Hi-C interaction. It specifies the genomic regions of
#'     the left and right Hi-C fragments. Furthermore, the IDs of input genomic regions and target genomic regions that
#'     interact through these Hi-C interactions are listed in each row (separated by comma).
#'   }
#' }
#'
#' @note The ID of a Hi-C interaction is defined by concatenating the left fragment ID and the right fragment ID (separated by dash).
#'
#' @examples
#'
#' @importFrom GenomicRanges findOverlaps GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom data.table data.table as.data.table
#'
#' @export
findHiCinteractions <- function(
    input_id,
    input_chr,
    input_start,
    input_end,

    hic_f1_id,
    hic_f1_chr,
    hic_f1_start,
    hic_f1_end,

    hic_f2_id,
    hic_f2_chr,
    hic_f2_start,
    hic_f2_end,

    target_id,
    target_chr,
    target_start,
    target_end) {






    if(length(setdiff( unique(input_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for input regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(hic_f1_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for HiC-F1 must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(hic_f2_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for HiC-F2 must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }
    if(length(setdiff( unique(target_chr), paste0('chr',c(1:22,'X','Y','M')))) != 0){
        stop('"chromosome" vector for target regions must be a subset of {chr1,chr2,...,chr22,chrX,chrY,chrM}')
    }

    temp <- length(unique(length(input_id),length(input_chr),length(input_start),length(input_end)))
    if(temp == 0){ stop('NULL vectors for input regions') }
    if(temp > 1){ stop('the lengths of vectors for input regions are unequal') }

    temp <- length(unique(length(hic_f1_id),length(hic_f1_chr),length(hic_f1_start),length(hic_f1_end)))
    if(temp == 0){ stop('NULL vectors for HiC (fragment 1)') }
    if(temp > 1){ stop('the lengths of vectors for HiC (fragment 1) are unequal') }

    temp <- length(unique(length(hic_f2_id),length(hic_f2_chr),length(hic_f2_start),length(hic_f2_end)))
    if(temp == 0){ stop('NULL vectors for HiC (fragment 2)') }
    if(temp > 1){ stop('the lengths of vectors for HiC (fragment 2) are unequal') }

    temp <- length(unique(length(target_id),length(target_chr),length(target_start),length(target_end)))
    if(temp == 0){ stop('NULL vectors for target regions') }
    if(temp > 1){ stop('the lengths of vectors for target regions are unequal') }


    cat(1,'\n')


    # the Hi-C data is expected to have these columns: {F1_ID, F1_chr, F1_start, F1_end, F2_ID, F2_chr, F2_start, F2_end}
    HiC_all_sig_interactions <- data.table(F1_ID = hic_f1_id, F1_chr = hic_f1_chr, F1_start = hic_f1_start, F1_end = hic_f1_end,
                                           F2_ID = hic_f2_id, F2_chr = hic_f2_chr, F2_start = hic_f2_start, F2_end = hic_f2_end)

    HiC_all_sig_interactions$uni_ID <- paste0(HiC_all_sig_interactions$F1_ID, '-', HiC_all_sig_interactions$F2_ID)



    # the targets data is expected to have these columns: {target_ID, chr, start, end}
    targets <- data.table(target_ID = target_id, chr = target_chr, start = target_start, end = target_end)




    input_ranges <- GRanges(seqnames=input_chr, ranges=IRanges(input_start, input_end))

    HiC_F1_ranges <- GRanges(seqnames=HiC_all_sig_interactions$F1_chr,
                             ranges=IRanges(HiC_all_sig_interactions$F1_start, HiC_all_sig_interactions$F1_end))
    HiC_F2_ranges <- GRanges(seqnames=HiC_all_sig_interactions$F2_chr,
                             ranges=IRanges(HiC_all_sig_interactions$F2_start, HiC_all_sig_interactions$F2_end))

    targets_ranges <- GRanges(seqnames=targets$chr,
                              ranges=IRanges(targets$start, targets$end))

    cat(2,'\n')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    input_ovlp_with_HiC <- rbind(
        data.table(as.data.table(suppressWarnings(findOverlaps(input_ranges, HiC_F1_ranges))), fragment = 'L'),
        data.table(as.data.table(suppressWarnings(findOverlaps(input_ranges, HiC_F2_ranges))), fragment = 'R')
    )
    input_ovlp_with_HiC <- input_ovlp_with_HiC[ which ( ! (
        duplicated(input_ovlp_with_HiC, by = c('queryHits', 'subjectHits'), fromLast = F) |
            duplicated(input_ovlp_with_HiC, by = c('queryHits', 'subjectHits'), fromLast = T) ) ) ]

    colnames(input_ovlp_with_HiC) <- c('input_indx', 'HiC_ovlp_with_input', 'fragment_ovlp_with_input')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    cat(3,'\n')

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    targets_ovlp_with_HiC <- rbind(
        data.table(as.data.table(suppressWarnings(findOverlaps(targets_ranges, HiC_F1_ranges))), fragment = 'L'),
        data.table(as.data.table(suppressWarnings(findOverlaps(targets_ranges, HiC_F2_ranges))), fragment = 'R')
    )
    targets_ovlp_with_HiC <- targets_ovlp_with_HiC[ which ( ! (
        duplicated(targets_ovlp_with_HiC, by = c('queryHits', 'subjectHits'), fromLast = F) |
            duplicated(targets_ovlp_with_HiC, by = c('queryHits', 'subjectHits'), fromLast = T) ) ) ]

    targets_ovlp_with_HiC <- targets_ovlp_with_HiC[, list(overlapping_targets = paste0(targets$target_ID[queryHits], collapse = ',')), by = list(subjectHits, fragment)]
    colnames(targets_ovlp_with_HiC)[c(1,2)] <- c('HiC_ovlp_with_target', 'fragment_ovlp_with_target')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    cat(4,'\n')

    target_fragment <- c('L','R')
    names(target_fragment) <- c('R','L')
    targets_ovlp_with_HiC$target_fragment_for_target <- as.character(target_fragment[targets_ovlp_with_HiC$fragment_ovlp_with_target])


    result_ovlps <-
        cbind(input_ovlp_with_HiC,
              targets_ovlp_with_HiC[,c('HiC_ovlp_with_target', 'fragment_ovlp_with_target', 'overlapping_targets')] [ match(
                  paste(input_ovlp_with_HiC$HiC_ovlp_with_input, input_ovlp_with_HiC$fragment_ovlp_with_input),
                  paste(targets_ovlp_with_HiC$HiC_ovlp_with_target, targets_ovlp_with_HiC$target_fragment_for_target) ) ] )
    result_ovlps <- result_ovlps[ ! is.na(result_ovlps$HiC_ovlp_with_target)]

    result_ovlps <- result_ovlps[,-4]
    colnames(result_ovlps)[2] <- 'HiC_indx'

    result_ovlps$input_ID <- input_id[result_ovlps$input_indx]
    result_ovlps$HiC_uni_ID <- HiC_all_sig_interactions$uni_ID[result_ovlps$HiC_indx]


    cat(5,'\n')


    # unfolding result_ovlps based on overlapping_targets
    reps <- sapply(strsplit(result_ovlps$overlapping_targets, ','), length)
    unfolded_indxs <-  unlist(lapply(1:length(reps), function(i) rep(i, reps[i])))
    unfolded_overlapping_targets <- unlist(strsplit(result_ovlps$overlapping_targets, ','))
    result_ovlps <- result_ovlps[unfolded_indxs]
    result_ovlps$overlapping_targets <- unfolded_overlapping_targets
    colnames(result_ovlps)[5] <- 'target_ID'
    # - - - - - - - - - - - - - - - - - - - - - - - - -

    cat(6,'\n')

    summary_ovlps <- result_ovlps[ ,
                                   list(HiC_IDs = paste0(HiC_uni_ID, collapse = ',')),
                                   by = list(input_indx, target_ID)] [ ,
                                                                      list(HiC_interacting_target_IDs = paste0(target_ID, collapse = ';'),
                                                                           HiC_interaction_IDs = paste0(HiC_IDs, collapse = ';')),
                                                                      by = input_indx]




    cat(7,'\n')


    ############################################################################################################
    # create "output"

    output <- data.table(input_ID = input_id, input_chr = input_chr, input_start = input_start, input_end = input_end)

    output[['number_of_target_IDs_interacting_through_HiC']] <- 0
    output[['number_of_target_IDs_interacting_through_HiC']][summary_ovlps$input_indx] <- sapply(strsplit(summary_ovlps$HiC_interacting_target_IDs, ';'),length)

    output[['target_IDs_interacting_through_HiC']] <- ''
    output[['target_IDs_interacting_through_HiC']][summary_ovlps$input_indx] <- summary_ovlps$HiC_interacting_target_IDs

    output[['HiC_interaction_IDs']] <- ''
    output[['HiC_interaction_IDs']][summary_ovlps$input_indx] <- summary_ovlps$HiC_interaction_IDs
    ############################################################################################################


    cat(8,'\n')


    ############################################################################################################
    # create "HiC_interactions_report"

    result_ovlps$F1_input_ID <- result_ovlps$input_ID
    result_ovlps$F1_input_ID[which(result_ovlps$fragment_ovlp_with_input != 'L')] <- ''

    result_ovlps$F2_input_ID <- result_ovlps$input_ID
    result_ovlps$F2_input_ID[which(result_ovlps$fragment_ovlp_with_input != 'R')] <- ''

    result_ovlps$F1_target_ID <- result_ovlps$target_ID
    result_ovlps$F1_target_ID[which(result_ovlps$fragment_ovlp_with_target != 'L')] <- ''

    result_ovlps$F2_target_ID <- result_ovlps$target_ID
    result_ovlps$F2_target_ID[which(result_ovlps$fragment_ovlp_with_target != 'R')] <- ''



    summary_HiC <- result_ovlps[, list(F1_input_IDs = paste0(unique(F1_input_ID[!(F1_input_ID == '')]), collapse = ','),
                                       F2_input_IDs = paste0(unique(F2_input_ID[!(F2_input_ID == '')]), collapse = ','),
                                       F1_target_IDs = paste0(unique(F1_target_ID[!(F1_target_ID == '')]), collapse = ','),
                                       F2_target_IDs = paste0(unique(F2_target_ID[!(F2_target_ID == '')]), collapse = ',')), by = list(HiC_indx, HiC_uni_ID)]


    HiC_indx <- summary_HiC$HiC_indx
    HiC_interactions_report <- cbind(summary_HiC[,'HiC_uni_ID'],
                                HiC_all_sig_interactions[,c('F1_ID','F1_chr','F1_start','F1_end')][HiC_indx],
                                summary_HiC[,c('F1_input_IDs','F1_target_IDs')],
                                HiC_all_sig_interactions[,c('F2_ID','F2_chr','F2_start','F2_end')][HiC_indx],
                                summary_HiC[,c('F2_input_IDs','F2_target_IDs')])


    colnames(HiC_interactions_report) <- c('HiC_interaction_ID',
                                      'F1_ID','F1_chr','F1_start','F1_end',
                                      'input_IDs_overlapping_with_F1','target_IDs_overlapping_with_F1',
                                      'F2_ID','F2_chr','F2_start','F2_end',
                                      'input_IDs_overlapping_with_F2','target_IDs_overlapping_with_F2')

    ############################################################################################################


    cat(9,'\n')


    return(list(output=output, HiC_interactions_report=HiC_interactions_report))
}



