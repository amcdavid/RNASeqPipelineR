## deduplicate version .4
## 12/07/15
## Fixes an issue with strand (which wasn't needed)
## Uses the primary alignment for deduplication, which should speed things up
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

RSEMstripExtension <- function(x){
    stri_replace(x,  '', regex='(\\.transcript)?\\.bam$')
}

##' Deduplicate a folder of .bams by using UMIs
##'
##' Each read is compared to an \emph{equivalence class} of reads that have the same mapping coordinates and UMI.
##' If the sequences match within a certain edit distance (see \code{getUniqueQname}), then all except one read (the highest quality read) are removed.
##' Note that for paired end reads, only the left-most (3 prime) mate will typically be compared in this process
##' @param ncores \code{integer} number of parallel processes to use. Be careful not to exceed the available memory--children processes appear to burst up to 2 gigs ram per million lines of .BAM to be read in at a time.  (You can control this by setting chunksize to a smaller value).
##' @param bam.path If provided, \code{character} giving the folder of input .bams.  If missing, then use the \code{RSEM} folder.
##' @param destination.prefix If provided, a \code{character} giving the folder of output .bams.  If missing, then use '../DEDUPLICATED_BAM' (evaluated relative to bam.path)
##' @param ... additional arguments passed to \code{writeDeduplicatedBam}
##' @return data.table with statistics
##' @import ShortRead
##' @import stringi
##' @import data.table
##' @import Rsamtools
##' @export
deduplicateBam <- function(ncores=1, bam.path, destination.prefix='../BAM/', ...){
    if(missing(bam.path)){
        bam.path <- getConfig()[["subdirs"]][["RSEM"]]
    }
    bamin <- list.files(path=bam.path, pattern='*.(transcript.)?bam$', full.names=TRUE)

    isRelative <- !stri_detect(destination.prefix, regex='^/')
    dest <- if(isRelative) file.path(dirname(bamin[1]), destination.prefix) else destination.prefix
    if(!dir.exists(dest)) {
        message('Making directory ', dest)
        system(paste0('mkdir -p ', dest))
    }

    bamout <- list.files(dest, pattern='*.(transcript.)?bam$', full.names=TRUE)
    no <-  RSEMstripExtension(basename(bamin)) %in% RSEMstripExtension(basename(bamout))
    bamyes <- bamin[!no]
    if(length(bamyes)>0){
        message("Deduplicating ", length(bamyes), " files to ", dest)
        if(ncores>1){
            out <- parallel::mclapply(bamyes, function(x) writeDeduplicatedBam(x, destination.prefix=dest, ...), mc.cores=ncores)
        } else{
            out <- lapply(bamyes, function(x) writeDeduplicatedBam(x, destination.prefix=dest, ...))
        }
        err <- sapply(out, inherits, 'try-error')
        if(any(err)){
            warning('Problems deduplicating file(s) ', paste(bamyes[err], collapse=','))
            warning('Errors were:')
            warning(out[err])
        }
        existingStats <- list()
        try({
            statsRDS <- file.path(getConfig()[["subdirs"]][["STATS"]], 'deduplicateStats.rds')
            
            if(file.exists(statsRDS)) existingStats <- readRDS(statsRDS)
        })
        out <- c(existingStats, out)
        saveRDS(out, file=statsRDS)
        invisible(out)
    } else{
        message("Nothing to deduplicate")
    }
}

SEQ_IDX_START <- 5
SEQ_IDX_LEN <- 7
collapseLowQualityReads <- function(dt, minProb=.3, maxProb=.8){
    minQ <- floor(-10*log10(1-minProb)+.5)
    maxQ <- floor(-10*log10(1-maxProb)+.5)
    setkey(dt, qname, mapq)
    dt[,':='(nminQ=sum(mapq>minQ),
           nhiQ=sum(mapq>maxQ)), key=list(qname)]
    dt[!is.na(pos) & nminQ<1, rname:=substr(seq, SEQ_IDX_START, SEQ_IDX_START+SEQ_IDX_LEN-1)]
    dt[,rank:=seq_len(.N),keyby=list(qname, mapq)]
    dt[,keep:=FALSE]
    dt[rank==1, keep:=TRUE] #both hiQ>0 or nminQ < 1
    dt[rank<5 & nminQ>0 & nhiQ<1, keep:=TRUE]
    dt
}

##bamfilename: character giving filename to bam
##destination.prefix: string (maybe directory name) to prepend to output
## index=FALSE: is the bam filed indexed by its coordinate system?  Leave FALSE.
## chunksize: how many records to read in at once?
## trim.len: how far should the reads be trimmed for the purposes of finding duplicates
## return.stats: should statistics regarding the deduplication be returned
## write: should the deduplicated files be written (or just stats returned).
writeDeduplicatedBam <- function(bamfilename, destination.prefix='../DEDUPLICATED_BAM/', chunksize=2e6, trim.len=50, return.stats=TRUE,write=TRUE,debug=0, ...){
    what <- c("qname", "rname", "pos", "seq", "qual", "mapq")
    ## For the deduplication, it suffices to get the *primary* mapping of the forward read
    ## (Because we just use the first qname at the moment, anyways)
    flag <- scanBamFlag(isSecondaryAlignment=NA)
    bf <- open(BamFile(bamfilename, yieldSize=chunksize))
    uniq <- list()
    i <- 1
    repeat{
        nrbam <- 0
        if(debug>1) message('Processing chunk ', i, ' of ', bamfilename)
        bamdtlist <- list()
        j <- 1
        while(nrbam<chunksize){
            bam <- scanBam(bf, param=ScanBamParam(what=what, flag=flag))[[1]]
            if(length(bam$qname)==0) break
        ## squal is negative sum of qualities (so that sorting in ascending order gives highest quality reads first)
            bamdtlist[[j]] <- data.table(qname = bam$qname, rname=bam$rname, pos=bam$pos, idx=seq_along(bam$qname), seq=substr(as.character(bam$seq), 1, trim.len),  squal=-alphabetScore(bam$qual), mapq=bam$mapq)
            bamdtlist[[j]] <- collapseLowQualityReads(bamdtlist[[j]])
            nrbam <- nrow(bamdtlist[[j]]) + nrbam
            j <- j +1
        }
        bamdt <- rbindlist(bamdtlist)
        rm(bam)
        gc()
        uniq[[i]] <- getUniqueQname(bamdt, debug=debug, ...)
        rm(bamdt)
        i <- i+1
    }
    close(bf)
    ## in case we had to process the file in chunks
    badqname <- unique(do.call(c, lapply(uniq, '[[', 'badqnames')))
    goodqnametable <- unique(rbindlist(lapply(uniq, '[[', 'goodqnames')))
    goodqname <- sort(goodqnametable$qname)
    ## gqn <- as.list(rep(TRUE, length(goodqname)))
    ## names(gqn) <- goodqname
    ## setkey(goodqnametable, qname)
    message('In file ', bamfilename, ' ', length(goodqname)+length(badqname), 
' reads total')     
    message('Keeping ', nrow(goodqnametable), ' of which ', 
         nrow(goodqnametable[multiplicity>0]), ' are mapped.')
    filterFunc <- function(DF) !(DF$qname %in% badqname)
    FR <- IRanges::FilterRules(filterFunc)
    #FR <- IRanges::FilterRules(function(DF) !(DF$qname %in% badqname))
    destname <- file.path(destination.prefix, paste0(RSEMstripExtension(basename(bamfilename))))
    dest <- stats <- NULL
    if(write){
        tmpdest <- paste0(destname, '.incomplete') #so that if we error out we don't leave incomplete bam files sitting around
        dest <- filterBam(bamfilename, destination=tmpdest, index=character(0), filter=FR, params=ScanBamParam(what='qname'), indexDestination=FALSE, yieldSize=floor(chunksize/5))
        file.rename(tmpdest, paste0(destname, '.bam'))
    }
    if(return.stats){
        multiplicityUnmapped <- do.call(c, lapply(uniq, '[[', 'multiplicity'))      
        stats <- list(potentialDupED=do.call(c, lapply(uniq, '[[', 'potentialDupED')),
                      multiplicity=multiplicityUnmapped[multiplicityUnmapped>0])
    }
    
    list(ndiscard=length(badqname), nkeep=length(setdiff(goodqname, badqname)), dest=RSEMstripExtension(basename(bamfilename)), stats=stats)
}



##' Generate a candidate list of *duplicate* reads using UMI/Sequence
##' 
##' The mapping position (POS/rname) is used as heuristic to find reads with
##' high sequence similarities
##'
##' @param bamdt a data.table with qname, rname, pos and seq
##' @param umipattern regular expression matching the UMIs in the qnames
##' @param max.edit.dist distance below which we declare two candidate to be duplicates
##' @param debug higher numbers result in more verbose output
##' @return list with elements badqname (a character vector of qnames to discard)
## goodqname (data.table of qnames to keep)
## ed (sample of edit distances of potential duplicates, useful to tune max.edit.dist)
## multiplicity (number of times each duplicate appeared, of interest to estimate fragment dropout rate)
getUniqueQname <- function(bamdt,  umipattern='[ACGT]+$', max.edit.dist=8, debug=1){
    tic <- Sys.time()
    setkey(bamdt, qname, rname, pos, squal)
    
    ## list of qnames
    ## keep==TRUE when they'll be kept
    goodqname <- unique(bamdt[,list(qname, keep=FALSE, multiplicity=NA_integer_)])
    setkey(goodqname, qname)

    ## dtdup contains the list of unprocessed reads
    ## save the unmapped reads because RSEM uses them to estimate error rates
    dtdup <- bamdt[!is.na(pos),]
    goodqname[bamdt[is.na(pos), qname], ':='(keep=TRUE,
                                             multiplicity=0L)]

    ## first matching alignment for each read
    ##    (because we use alignments as heuristic to place similar reads into
    ##     equivalence classes )
    dtdup <- dtdup[unique(dtdup$qname),,mult='first']
    ## track statistics relating to edit distances and multiplicities
    multiplicity <- ed <- rep(NA, min(1e4, nrow(dtdup)))
    ## extract UMI from qname
    dtdup[,umi:=stri_extract_last(qname, regex=umipattern)]
    ## sort by rname, pos and umi
    setkey(dtdup, rname, pos, umi, squal)

    
    ## for each equivalence class of reads
    ## save first read
    ## test rest for similarity to first read
    ## delete matches
    ## repeat until no more deletions
    ## (do it this way to prevent quadratic complexity in equivalence class size)
    j <- i <- 1
    while(nrow(dtdup)>0){
        if(debug>2) message("Round ", j)
        if(debug>1 & nrow(dtdup)>1000) print(noquote(paste0(nrow(dtdup), ' records to process, ', sum(goodqname$keep), ' records kept.')))
        ## number of records with same rname, pos and umi
        ## --> number of records in equivalence class
        dtdup[,umidup:=.N, key=list(rname, pos, umi)]
        
        ## only one record.  keep them.
        ## Then pop these records from dtdup
        dt1 <- dtdup[umidup==1,.(qname)]
        if(debug>2) message(nrow(dt1), " singletons kept")
        goodqname[dt1$qname, ':='(keep=TRUE,
                          multiplicity=1)]
        dtdup <- dtdup[umidup>1,]
        if(nrow(dtdup)==0) break
        
	## distance between first sequence and all other members of equivalence class
        dtdup[,Dist:={
            adist(seq, seq[1], partial=FALSE)
        }, key=list(rname, pos, umi)]

        dtdup[,':='(umidup=sum(Dist<max.edit.dist),
                    isDup=Dist<max.edit.dist)
                 , key=list(rname, pos, umi)]
        setorder(dtdup, rname, pos, umi, squal) #sort by quality
        ## pop off first read, which should be highest quality
        dt1 <- dtdup[,.(qname=qname[1], umidup=umidup[1]),keyby=list(rname, pos, umi)]
        if(debug>2) message(nrow(dt1), " first members of equivalence class kept")
        setkey(dt1, qname)
        goodqname[dt1, c('keep', 'multiplicity'):=list(TRUE, umidup)]
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
    if(debug>1) message('Spent ', round(toc-tic, 1), ' seconds deduplicating. ', sum(goodqname$keep), ' unique reads.')
    badqnames <- goodqname[keep==FALSE,qname]
    
    ## sum of multiplicities should equal the number of mapped reads
    stopifnot(goodqname[keep==TRUE,sum(multiplicity)]==length(unique(bamdt[!is.na(pos), qname])))
    
    return(list(badqnames=badqnames,
                goodqnames=goodqname[keep==TRUE,.(qname, multiplicity)],
                potentialDupED=ed,
                multiplicity=goodqname[keep==TRUE,multiplicity],
                rname=goodqname[keep==TRUE,rname],
                umi=goodqname[keep==TRUE,umi]                
                ))
}

