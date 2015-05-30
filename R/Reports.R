##' Generate reports on partitioning/trimming/qc
##'
##' This copies 'adaptorStats.Rmd' to the current working directory
##' and runs \code{knit2html('adaptorStats.Rmd')} to generate a report on "partition_stats.rds".
##'
##' You can edit the .Rmd as you see fit to customize the analysis.
##' @return writes a .html
##' @export
partitionReport <- function(){
    file.copy(file.path(system.file(package='RNASeqPipelineR'), 'Rmd', 'adaptorStats.Rmd'), getwd(), overwrite=FALSE)
    knit2html('adaptorStats.Rmd')
}

##' Generate a report on mapping/deduplication
##'
##' @param ... named arguments giving directories in which RSEM was run.  By default, looks at \code{getConfig()[["RSEM"]]}
##' @describeIn partitionReport copies 'qc.Rmd' and generates a report on "deduplicateStats.RData" and on mapping statistics located in various RSEM directories 
##' @export
mappingReport <- function(...){
    rsemloc <- list(...)
    if(length(rsemloc)==0) rsemloc <- list(dup=getConfig()[['subdirs']][["RSEM"]])

    file.copy(file.path(system.file(package='RNASeqPipelineR'), 'Rmd', 'qc.Rmd'), getwd(), overwrite=FALSE)
    knit2html('qc.Rmd')
}

processDupStats <- function(ddb){
dest <- sapply(ddb, '[[', 'dest')
nkeep <- sapply(ddb, '[[', 'nkeep')
ndiscard <- sapply(ddb, '[[', 'ndiscard')
stats <- lapply(ddb, '[[', 'stats')

nkeepDT <- data.table(dest, nkeep, ndiscard, L1=seq_along(dest), key='L1')
mStats <- reshape2::melt(stats)
setDT(mStats)
setkey(mStats, L1)
mStats <- mStats[nkeepDT]
## why are there NAs?
mStats <- mStats[!is.na(value)]

multipl <- mStats[L2=='multiplicity',]
multipl <- multipl[,.SD[seq(1, min(500, .N)),],keyby=dest]
multipl[,':='(medianMult=median(value),
            pct75=quantile(value, .75),
              avgMult=mean(value)), keyby=dest]
setkey(multipl, dest)
ed <- mStats[L2=='potentialDupED',]
list(mQuantile=multipl, mtotal=mStats[L2=='multiplicity'], ed=ed)
}

makeSCA <- function(tpmfile, fdatfile, cdat){
    ee <- as.matrix(fread(tpmfile))
    rownames(ee) <- ee[,1]
    ee <- ee[,-1]
    fdat <- fread(fdatfile)
    fdat[,primerid:=str_replace(entrez_id, ',.*', '')]
    fdat[,primerid:=make.unique(ifelse(is.na(primerid), 'unknown', primerid))]
    ee <- ee[fdat$gene_id,]
    stopifnot(nrow(fdat) == nrow(ee))
    rownames(ee) <- fdat$primerid
    ee <- log2(t(ee)+1)
    cdat[,wellKey:=sampleID]
    setkey(cdat, wellKey)
    cdat <- cdat[rownames(ee)]
    fm <- MAST::FromMatrix('RNASeqAssay', ee, cdat, fdat)
    fm <- fm[,MAST::freq(fm)>.02]
    ngo <- rowSums(MAST::exprs(fm)>0)
    log10size <- rowSums(MAST::exprs(fm))
    fm <- MAST::combine(fm, data.frame(ngeneson=ngo, cngeneson=ngo-mean(ngo), log10size))
    fm
}
