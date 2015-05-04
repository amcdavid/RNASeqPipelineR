## deduplicate version .3
## Keep highest quality read (based on sum of qscores)
## Counts multiplicities more carefully
## TODO:
## -Analyze/Report sequence but not UMI duplicates
## -abstract data collection
## Andrew McDavid


if(getRversion() >= "2.15.1") globalVariables(c('qname',
                  'rname',
                  'pos',
                  'squal',
                  'umi',
                  'umidup',
                  'isDup',
                  'keep',
                  'Dist',
                  '.'))

stripExtension <- function(x){
    stri_replace_last(x,  '', regex='(\\.transcript)?\\.bam$')
}

##' Deduplicate a folder of .bams by using UMIs
##'
##' Each read is compared to an \emph{equivalence class} of reads that have the same mapping coordinates and UMI.
##' If the sequences match within a certain edit distance (see \code{getUniqueQname}), then all except one read (the highest quality read) are removed.
##' Note that for paired end reads, only the left-most (3 prime) mate will typically be compared in this process
##' @param ncores \code{integer} number of parallel processes to use
##' @param bam.path If provided, \code{character} giving the folder of input .bams.  If missing, then use the \code{RSEM} folder.
##' @param destination.prefix If provided, a \code{character} giving the folder of output .bams.  If missing, then use '../DEDUPLICATED_BAM' (evaluated relative to bam.path)
##' @param ... additional arguments passed to \code{writeDeduplicatedBam}
##' @return data.table with statistics
##' @import ShortRead
##' @import stringi
##' @import data.table
##' @import Rsamtools
##' @export
deduplicateBam <- function(ncores=1, bam.path, destination.prefix='../DEDUPLICATED_BAM/', ...){
    if(missing(bam.path)){
        bam.path <- getConfig()[["subdirs"]][["RSEM"]]
    }
    bamin <- list.files(path=bam.path, pattern='*transcript.bam$', full.names=TRUE)

    isRelative <- stri_detect(destination.prefix, regex='^/')
    dest <- if(isRelative) file.path(dirname(bamin[1]), destination.prefix) else destination.prefix
    if(!dir.exists(dest)) {
        message('Making directory ', dest)
        system(paste0('mkdir -p ', dest))
    }

    bamout <- list.files(dest, pattern='*.(transcript.)?bam$', full.names=TRUE)
    no <-  stripExtension(basename(bamin)) %in% stripExtension(basename(bamout))
    bamyes <- bamin[!no]
    if(length(bamyes)>0){
        message("Deduplicating ", length(bamyes), " files to ", dest)
        out <- parallel::mclapply(bamyes, function(x) writeDeduplicatedBam(x, destination.prefix=dest, ...), mc.cores=ncores)
        out <- lapply(bamyes[1:5], function(x) writeDeduplicatedBam(x, destination.prefix=dest, ...))
        err <- sapply(out, inherits, 'try-error')
        if(any(err)) warning('Problems deduplicating file(s) ', paste(bamyes[err], collapse=','))
        save(out, file=file.path(getConfig()[["subdirs"]][["STATS"]], 'deduplicateStats.RData'))
        invisible(out)
    } else{
        message("Nothing to deduplicate")
    }
}


##bamfilename: character giving filename to bam
##destination.prefix: string (maybe directory name) to prepend to output
## index=FALSE: is the bam filed indexed by its coordinate system?  Leave FALSE.
## chunksize: how many records to read in at once?
## trim.len: how far should the reads be trimmed for the purposes of finding duplicates
## return.stats: should statistics regarding the deduplication be returned
## write: should the deduplicated files be written (or just stats returned).
writeDeduplicatedBam <- function(bamfilename, destination.prefix='../DEDUPLICATED_BAM/', chunksize=5e6, trim.len=50, return.stats=FALSE,write=TRUE,...){
    what <- c("qname", "rname", "pos", "seq", "strand", "qual")
    bf <- open(BamFile(bamfilename, yieldSize=chunksize))
    uniq <- list()
    i <- 1
    repeat{
        bam <- scanBam(bf, param=ScanBamParam(what=what), )[[1]]
        if(length(bam$qname)==0) break
        ## squal is negative sum of qualities (so that sorting in ascending order gives highest quality reads first)
        bamdt <- data.table(qname = bam$qname, rname=bam$rname, pos=bam$pos, idx=seq_along(bam$qname), seq=substr(as.character(bam$seq), 1, trim.len), strand=bam$strand, squal=-alphabetScore(bam$qual))
        rm(bam)
        uniq[[i]] <- getUniqueQname(bamdt, ...)
        i <- i+1
    }
    close(bf)
    ## in case we had to process the file in chunks
    badqname <- unique(do.call(c, lapply(uniq, '[[', 'badqnames')))
    goodqnametable <- unique(rbindlist(lapply(uniq, '[[', 'goodqnames')))
    goodqname <- goodqnametable$qname
    message('In file ', bamfilename, ' ', length(goodqname)+length(badqname), 
' reads total')     
    message('Keeping ', nrow(goodqnametable), ' of which ', 
         nrow(goodqnametable[multiplicity>0]), ' are mapped.')
    FR <- IRanges::FilterRules(function(DF) !(DF$qname %in% badqname))
    destname <- file.path(destination.prefix, paste0(stripExtension(basename(bamfilename)), '.bam'))
    dest <- stats <- NULL
    if(write){
        dest <- filterBam(bamfilename, destination=destname, index=character(0), filter=FR, params=ScanBamParam(what='qname'), indexDestination=FALSE)
    }
    if(return.stats){
        multiplicityUnmapped <- do.call(c, lapply(uniq, '[[', 'multiplicity'))      
        stats <- list(potentialDupED=do.call(c, lapply(uniq, '[[', 'potentialDupED')),
                      multiplicity=multiplicityUnmapped[multiplicityUnmapped>0])
    }
    
    list(ndiscard=length(badqname), nkeep=length(setdiff(goodqname, badqname)), dest=stripExtension(basename(bamfilename)), stats=stats)
}


## bamdt: a data.table with qname, rname, pos and seq
## umipattern: regular expression to extract UMI from qname
## max.edit.dist: distance below which we declare two candidate to be duplicates
## debug: larger values print more diagnostics
## returns a list with elements badqname (a character vector of qnames to discard)
## goodqname (data.table of qnames to keep)
## ed (sample of edit distances of potential duplicates, useful to tune max.edit.dist)
## multiplicity (number of times each duplicate appeared, of interest to estimate fragment dropout rate)
getUniqueQname <- function(bamdt,  umipattern='[ACGT]+$', max.edit.dist=8, debug=1){ 
    setkey(bamdt, qname, rname, pos, squal)
    #list of qnames and whether or not they'll be kept or discarded
    goodqname <- unique(bamdt[,list(qname, keep=FALSE, multiplicity=NA_integer_)])## , strand=bamdt$strand[1])])
    setkey(goodqname, qname)
    goodqname[bamdt[is.na(pos), qname], ':='(keep=TRUE,
                                             multiplicity=0L)] #save the unmapped reads because RSEM uses them to estimate error rates
    dtdup <- bamdt[!is.na(pos),]
    ## first matching alignment for each read
    dtdup <- dtdup[unique(dtdup$qname),,mult='first']
    if(debug>0)
    #statistics relating to edit distances and multiplicities
    multiplicity <- ed <- rep(NA, min(1e4, nrow(dtdup)))
    #extract UMI from qname
    dtdup[,umi:=stri_extract_last(qname, regex=umipattern)]
    #sort by rname, pos and umi
    setkey(dtdup, rname, pos, umi, squal)
    ## for each equivalence class of reads
    ## save first read
    ## test rest for similarity to first read
    ## delete matches
    ## repeat until no more deletions
    ## (do it this way to prevent quadratic complexity in equivalence class size)
    tic <- Sys.time()
    j <- i <- 1
    while(nrow(dtdup)>0){
        if(debug>2) message("Round ", j)
        if(debug>1 & nrow(dtdup)>1000) print(noquote(paste0(nrow(dtdup), ' records to process, ', sum(goodqname$keep), ' records kept.')))
        ##count number of records with same rname pos and umi
        dtdup[,umidup:=.N, key=list(rname, pos, umi)]
        ## only one record.  keep them.  (only needs to happen on first round)
        dt1 <- dtdup[umidup==1,.(qname, strand)]
        if(debug>2) message(nrow(dt1), " singletons kept")
        goodqname[dt1$qname, ':='(keep=TRUE,
                          multiplicity=1)]#, strand=dt1$strand)]
        ## current set of potential dups
        dtdup <- dtdup[umidup>1,]
	#distance between first sequence and all other members of equivalence class
        dtdup[,Dist:={
            adist(seq, seq[1], partial=FALSE)
        }, key=list(rname, pos, umi)]
        if(nrow(dtdup)==0) break
        dtdup[,':='(umidup=sum(Dist<max.edit.dist),
                    isDup=Dist<max.edit.dist)
                 , key=list(rname, pos, umi)]
        setorder(dtdup, rname, pos, umi, squal) #sort by quality
        ## pop off first read, which should be highest quality
        dt1 <- dtdup[,.(qname=qname[1], strand=strand[1], umidup=umidup[1]),keyby=list(rname, pos, umi)]
        if(debug>2) message(nrow(dt1), " first members of equivalence class kept")
        setkey(dt1, qname)
        goodqname[dt1, c('keep', 'strand', 'multiplicity'):=list(TRUE, strand, umidup)]
        ## index of edit distance matches
	##save a subsample of edit distances
        samp <- min(nrow(dtdup), 100)
        ed[i:(i+samp-1)] <- sample(dtdup[,Dist], samp)
        ##keep non-matches (delete matches)
        if(debug>2) message(sum(dtdup$isDup), " other members of equivalence class dropped")
        dtdup <- dtdup[isDup==FALSE,]
        i <- i + samp
        j <- j +1
    }
    toc <- Sys.time()
    if(debug>1)
    message('Spent ', round(toc-tic, 0), ' seconds deduplicating. ', sum(goodqname$keep), ' unique reads.')
    badqnames <- goodqname[keep==FALSE,qname]
    stopifnot(goodqname[keep==TRUE,sum(multiplicity)]==length(unique(bamdt[!is.na(pos), qname])))
    return(list(badqnames=badqnames, goodqnames=goodqname[keep==TRUE,.(qname, strand, multiplicity)], potentialDupED=ed, multiplicity=goodqname[keep==TRUE,multiplicity]))
}

