% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateSingleScore.R
\name{annotateSingleScore}
\alias{annotateSingleScore}
\title{Annotate the input genomic regions with annotation of type "single-score".}
\usage{
annotateSingleScore(input_id, input_chr, input_start, input_end, annot_chr,
  annot_start, annot_end, annot_score)
}
\arguments{
\item{input_id}{Character vector defining the name of input genomic regions (e.g. gene id)}

\item{input_chr}{Character vector defining the name of the chromosome for input genomic regions (one of chr1, chr2, ..., chrX, chrY or chrM).}

\item{input_start}{Numeric vector specifying the starting position of input genomic regions.}

\item{input_end}{Numeric vector specifying the ending position of input genomic regions.}

\item{annot_chr}{Character vector defining the name of the chromosome for annotation entries (one of chr1, chr2, ..., chrX, chrY or chrM).}

\item{annot_start}{Numeric vector specifying the starting position of annotation entries.}

\item{annot_end}{Numeric vector specifying the ending position of annotation entries.}

\item{annot_score}{Numeric vector specifying scores for annotation entries.}
}
\value{
A data.frame containing the input genomic regions and 2 extra columns which specify:
\enumerate{
  \item The average score of annotation over each input region.
  \item The percentage of overlap between each input region and annotation entries.
}
}
\description{
This function annotates the input genomic regions with a given "single-score" annotation. A "single-score" annotation specifies
a set of genomic regions and assigns a numeric score to each of them (e.g. the annotation of histone modification peaks by ENCODE).
}
\note{
For each input genomic region, the average score of annotation is measured as the weighted average of scores of
overlapping annotation entries where the percentages of overlaps of annotation entries with the input region are taken
as weights.
The percentage of overlap between an input genomic region with annotation entries is calculated as the
total portion of the input region that is covered by the annotation entries divided by the length of the region.
}
\examples{

}
