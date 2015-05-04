## trimm.R v0.16 15-04-2015
## Script for finding UMI/WellID in paired reads, trimming, and partitioning into new fastq files
## Jingyuan Deng and Andrew McDavid
## Updates: Check putative wellid matches for approximate match (hamming distance of 1) to expected in initial search for bridge
## Trim based on fastq qualities
## Only consider TTTAG to identify bridge (but restrict to first 30 nucleotides)
## Throw away UMI/wellID if qualities are too low and flag that read in the output.
## Place UMI/wellID before spaces in the qname


## maximum index we'll accept for the adaptor
MAX_BRIDGE_INDEX <- 30
## Adaptor sequence
ADAPTER <- 'TTTAG'
## Length of WellIDSeq
N_WELLID <- 7
## Length of UMI
N_UMI <- 5
## How much will we trim (beyond first adaptor nucleotide) when we find adaptor?
TRIM_OFFSET <- nchar(ADAPTER)+N_WELLID+N_UMI+5 # 5 = TGGGG due to template switching behavior from polymerase
## How many records to read from fastq files at once
YIELD_SIZE <- 1.5e6


if(getRversion() >= "2.15.1") globalVariables(c('wellIDecc',
                  'wellIDecc.fq1',
                  'wellIDecc.fq2',
                  'WellIDSeq',
                  'hadSuffix.fq1',
                  'hadSuffix.fq2',
                  'lowqual.fq1',
                  'lowqual.fq2',
                  'nsuffix',
                  'lowqual',
                  'hadBridge',
                  'UMI.fq1',
                  'UMI.fq2',
                  'exact.fq1',
                  'exact.fq2',
                  'reindex',
                  'N',
                  'N.x',
                  'N.y',
                  'fq1',
                  'fq2'))

##' Look for bridge adapter, UMI and wellID + TGGG.
##'
##' Report matching start positions for bridge, and matching UMI and wellID.
##' @param fq fastq
##' @param adapter \code{character} adapter sequence
##' @param adapterTol number of mismatches tolerated
##' @param suffix \code{character} expected terminating sequence
##' @param minQualityTol \code{integer} this number of bases below \code{minAlphaQual} will cause a putative adaptor match to be ignored
##' @param minAlphaQual \code{character} quality threshold as a phred \code{character}
##' @param wellIDTable \code{data.table} in wellIDsequences are taken
##' @return \code{data.table}
##' @import data.table
##' @import ShortRead
##' @import parallel
##' @import stringi
##' @import Biostrings
findBridge <- function(fq, adapter, adapterTol=0, suffix='', minQualityTol, minAlphaQual, wellIDTable){
    ## gives R a sad if we set > 1 when using mclapply
    nthreads <- max(getOption('mc.cores'), 2)
    seqs <- sread(fq)
    nsuffix <- nchar(suffix)
    nadapter <- nchar(adapter)
    ntotal <- nadapter+N_UMI + N_WELLID
    pos <- vmatchPattern(adapter, subseq(seqs, 1, MAX_BRIDGE_INDEX), max.mismatch=adapterTol, fixed=TRUE) ## find the positions of the matches
    ## if we match twice, we take first match.  seems sensible.
    start.pos.na <- sapply(startIndex(pos), function(el) if(is.null(el)) NA else el[1])  ## start positions of first match or NA if no match
    ## found adaptor (not NA) AND position less than MAX_BRIDGE_INDEX
    indices_bridge_filter <- which(start.pos.na < MAX_BRIDGE_INDEX)
    start.pos <- start.pos.na[indices_bridge_filter]
    ## don't attempt to use UMI/wellID from low quality sequence
    bridge_qual <-narrow(quality(fq[indices_bridge_filter]), start=start.pos + nadapter, end=start.pos+ntotal+nsuffix-1)
    ## gives R another sad when we use mclapply and i have _no_ idea why
    lowqual <- width(trimTails(bridge_qual, k=minQualityTol, a=minAlphaQual))- (ntotal+nsuffix  - nadapter)
    ##indices_bridge_filter <- setdiff(indices_bridge_filter, lowqual)
    start.pos <- start.pos.na[indices_bridge_filter]
    TGGG <- DNAStringSet(seqs[indices_bridge_filter], start=start.pos+ntotal, end=start.pos+ntotal+nsuffix-1)
    hadSuffix <- TGGG == suffix
    UMI <- substr(seqs[indices_bridge_filter], start.pos+nadapter, start.pos+nadapter+N_UMI-1)
    wellID <- substr(seqs[indices_bridge_filter], start.pos+nadapter+N_UMI, start.pos+ntotal-1)
    dt <- data.table(indices_bridge_filter, hadSuffix, start.pos, UMI, wellID, lowqual)
    dt[, wellIDecc:=wellIDTable[stringdist::amatch(wellID, wellIDTable, maxDist=1, method='hamming', nthread=nthreads)], keyby=wellID]
    dt[,':='(hadSuffix=hadSuffix & !is.na(wellIDecc),
             exact=vTRUE(wellID == wellIDecc))]
    dt
}

## !is.na && TRUE 
vTRUE <- function(x){
    x[is.na(x)] <- FALSE
    x
}

##' Use a sampleSheet to partition fastq by wellID
##' 
##' Each fastq (fq1 and fq2) will be divided by wellID.
##' Output file names are determined by column sampleID + _1.fastq and _2.fastq.
##' Samples that appear repeatedly in the sample sheet are concatenated.
##' Reading fastq files, and writing output is parallelized. Set \code{options(mc.cores=NCORES)} to utilize.
##' We can't easily parallelize over input fastq if we might be concatenating across input fastq files.
##' @param sampleSheet \code{data.table giving the map between wellID and input fastq file}
##' @param out.dir directory that partitioned fastq will be written.  If NULL, then will use the FASTQ directory from the project.
##' @param minQualityTol \code{integer} number of bases below minQuality that will be tolerated in bridge-umi-wellID-suffix
##' @param minQuality \code{numeric} \code{minQualityTol} or more bases below this numeric quality in the bridge-umi-wellID-suffix will result in the read being discarded
##' @return Returns a \code{data.table} providing summary statistics (how many codes could be determined unambiguously by wellID and input fq1, as well as how many reads couldn't be matched by code or were ambiguous).
##' @export
partitionFastqFiles <- function(sampleSheet, out.dir=NULL, minQualityTol=4, minQuality=15) {
    out.dir <- getConfig()[["subdirs"]][["FASTQ"]]
    uniqFiles <- unique(sampleSheet[,list(fq1, fq2)])
    stats <- vector('list', length=nrow(uniqFiles))
    cur.files <- list.files(out.dir, pattern='*.fastq')
    if(length(cur.files)>0) stop("`out.dir` must be empty")
    for(i in seq_len(nrow(uniqFiles))){
        f1 <- FastqStreamer(uniqFiles[i,fq1], n=YIELD_SIZE)
        f2 <- FastqStreamer(uniqFiles[i,fq2], n=YIELD_SIZE)
        ## a list holds forward and reverse fastq
        fq <- list(fq1=NA, fq2=NA)
        this.stat <- data.table(WellIDSeq=character(), N=integer())
        chunk <- 1-YIELD_SIZE
        ## Portion of sampleSheet corresponding to this input file
        thisSS <- sampleSheet[fq1==uniqFiles[i,fq1]]
        thisSSwid <- unique(thisSS$WellIDSeq)
        minAlphaQual <- NULL
        while(length(fq$fq1 <- yield(f1))){
            if(is.null(minAlphaQual)){
                E <- encoding(quality(fq[[1]]))
                minAlphaQual <- names(E)[which(E==minQuality)]
            }
            setkey(this.stat, WellIDSeq)
            chunk <- chunk+YIELD_SIZE
            message('File ', uniqFiles[i,fq1], ' chunk ', chunk)
            fq$fq2 <- yield(f2)
            message('Finding bridges...', appendLF=FALSE)
            indices <- lapply(fq, findBridge, adapter=ADAPTER, minQualityTol=minQualityTol, minAlphaQual=minAlphaQual, wellIDTable=thisSSwid)
            message('Found...', appendLF=FALSE)
            indices_all <- merge(indices$fq1, indices$fq2, by=c("indices_bridge_filter"), all=TRUE, suffixes=c(".fq1", ".fq2"))
            indices_all[ ,':='(nsuffix=vTRUE(hadSuffix.fq1)+ vTRUE(hadSuffix.fq2),
                               lowqual=pmin(lowqual.fq1*vTRUE(hadSuffix.fq1), lowqual.fq2*vTRUE(hadSuffix.fq2), na.rm=TRUE))]
            indices_good <- indices_all[nsuffix==1 & lowqual>=0, ]
            indices_good[,hadBridge:=ifelse(vTRUE(hadSuffix.fq1), 1, 2)]
            indices_good[,':='(UMI=UMI.fq1,
                               WellIDSeq=wellIDecc.fq1,
                               exact=exact.fq1,
                               lowqual=lowqual.fq1)]
            indices_good[hadBridge==2,':='(UMI=UMI.fq2,
                             WellIDSeq=wellIDecc.fq2,
                             exact=exact.fq2,
                             lowqual=lowqual.fq2)]
            ## NAs are mates without adapter--fill with -TRIM_OFFSET + 1 (so that trim starts at 1, ie, no trimming done)
            indices_good[is.na(indices_good)] <- -TRIM_OFFSET + 1
            indices_good <- merge(indices_good, thisSS, by='WellIDSeq')
            #rebuild trimmed fastqs
            message('Trimming...', appendLF=FALSE)
            mateTrim <- setNames(vector(mode='list', length=length(fq)), names(fq))
            setkey(indices_good, sampleID)
            indices_good[, reindex:=1:nrow(indices_good)]
            mateTrim <- mclapply(c('fq1', 'fq2'), function(j){
                thisfq <- fq[[j]][indices_good$indices_bridge_filter]
                start <- indices_good[[paste0('start.pos.', j)]]+TRIM_OFFSET
                seq_trim <- DNAStringSet(sread(thisfq), start=start)
                qual_trim <- BStringSet(quality(quality(thisfq)), start=start)
                qual_trim <- SFastqQuality(qual_trim) # reapply quality score type
                ## bowtie truncates the qnames when it encounters spaces
                ids_split <- stri_split_fixed(as.character(ShortRead::id(thisfq)), pattern=' ', n=2, simplify=TRUE)
                ids_annot <- paste(paste(ids_split[,1],  indices_good$WellIDSeq, indices_good$UMI, sep='!'), ids_split[,2])
                ids_trim <- BStringSet(ids_annot)
                ShortReadQ(seq_trim, quality=qual_trim, id=ids_trim)  # Rebuild reads object with trimmed sequences and quality scores
                #mateTrim[[j]] <- do.call(trimTails, c(list(object=mateTrim[[j]]), trim.control))
})
            stopifnot(length(mateTrim[[1]])==length(mateTrim[[2]]))
            
            sampleID <- unique(sampleSheet[fq1==uniqFiles[i,fq1]]$sampleID)
            message('Writing to disk...', appendLF=FALSE)
            mclapply(sampleID, function(s){
                for(j in seq_along(mateTrim))
                    writeFastq(mateTrim[[j]][indices_good[s,reindex,nomatch=0]], paste0(out.dir, s, "_", j, ".fastq"), mode='a', compress=FALSE)
            })
            yield.stat <- indices_good[,list(.N),key=list(WellIDSeq)]
            badSuffix <- indices_all[,.N,keyby=nsuffix][nsuffix != 1]
            nlowqual <- nrow(indices_all[nsuffix==1 & lowqual < 0,])
            yield.stat <- rbind(yield.stat,
                               data.table(WellIDSeq=c('UNSUFFIXED', 'AMBIGUOUS', 'LOWQUAL', 'NONE'),
                                          N=c(badSuffix[,N],
                                              nlowqual,
                                              length(fq[[1]])-nrow(indices_all))))
            setkey(yield.stat, WellIDSeq)
            this.stat <- merge(this.stat, yield.stat, all=TRUE, by='WellIDSeq')
            this.stat[is.na(this.stat)] <- 0
            this.stat <- this.stat[,list(WellIDSeq=WellIDSeq,            
                                         N=N.x+N.y)]
         message('Done.')   
        }
        this.stat[,fq1:=uniqFiles[i,'fq1', with=FALSE]]
        stats[[i]] <- this.stat
        close(f1)
        close(f2)
    }
    saveRDS(file.path(getConfig()[["subdirs"]][["STATS"]], 'partition_stats.rds'))
    invisible(rbindlist(stats))
}

## script to test output
testFastq <- function(prefix, expectedwellID){
    ## check that both files contain same fastq id (aside from :1 :2 differences)
    ## check for expectedwellID
    ## check that there's no adapter found
    ## (but occasionally, there's two bridge sequences found in a read)
}

