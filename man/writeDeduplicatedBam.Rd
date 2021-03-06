% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deduplicate.R
\name{writeDeduplicatedBam}
\alias{writeDeduplicatedBam}
\title{Deduplicate a single bam file}
\usage{
writeDeduplicatedBam(bamfilename, destination.prefix = "../DEDUPLICATED_BAM/",
  chunksize = 5e+06, trim.len = 50, stats.path = NULL, write = TRUE,
  debug = 0, ...)
}
\arguments{
\item{bamfilename}{character giving filename to bam}

\item{destination.prefix}{prefix, relative to the bamfilename path, to write output}

\item{chunksize}{How many records to read in at a time?}

\item{trim.len}{How long to trim reads for the purposes of finding duplicates?}

\item{stats.path}{If not NULL, the path to which  prefix_dedup_stats.rds files will be written.}

\item{write}{Should the deduplicated .bam be written (vs only calculate stats)}

\item{debug}{higher values report more diagnostics.}

\item{...}{arguments passed to getUniqueQname}
}
\value{
a list containing `ndiscard` (number of reads discarded), `nkeep` (number of reads kept), `dest` (the stem of the destination filename).
}
\description{
Parse a file and test for duplicate reads by comparing the rnames (reference name) and sequence identity for equivalency.
At the moment, this involves streaming over the file twice (once to mark duplicate rnaes
}

