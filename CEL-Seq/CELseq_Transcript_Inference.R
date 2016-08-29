
dumi = function(files, EM="global", mer.filter=0.01, aer.filter=0.01, verbose=TRUE) {

  bfl = BamFileList(files)
  
  ## Get summary stats for all files
  bam.align = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE, isNotPassingQualityControls=FALSE)))
  bam.unalign = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=TRUE)))
  bam.low.q = countBam(bfl, param=ScanBamParam(scanBamFlag(isUnmappedQuery=FALSE, isNotPassingQualityControls=TRUE)))
  
  ## Perform EM for each cell individually

  cell.mer = rep(NA, length(files))
  cell.mer.counts = rep(NA, length(files))
  cell.aer = rep(NA, length(files))  
  cell.aer.counts = rep(NA, length(files))
    
  for(i in 1:length(files)) {
    bem = barcodeEM(files(i), verbose=verbose)
    cell.aer[i] = bem$aer
    cell.mer[i] = bem$mer
    
    
  }
}





barcodeEM = function(filename, window=10, umi.flag=5, mismatch.error.rate = 0.05, alignment.error.rate = rep(1, (window*2)+1) / ((window*2)+1), max.iter=10, verbose=TRUE) {
  suppressPackageStartupMessages(require(stringr))
  suppressPackageStartupMessages(require(GenomicAlignments))
  suppressPackageStartupMessages(require(data.table))
  suppressPackageStartupMessages(require(stringdist))

  bamGA = readGAlignments(filename, param=ScanBamParam(what=c("qname", "flag")))
  
  umi = read.umi(mcols(bamGA)$qname, umi.flag)

  position = start(bamGA)
  strand = as.factor(strand(bamGA)) == "+"
  position[!strand] = end(bamGA)[!strand]

  info = data.table(data.frame(chr=seqnames(bamGA), position=position, strand=strand, umi=umi, inferred_position=position, inferred_umi=umi))
  info.counts = info[, `:=` (COUNT = .N, INITIAL_GRP=.GRP, IX=1:.N), by=list(chr, strand, inferred_position, inferred_umi)]
  info.counts = info.counts[,ID := 1:nrow(info.counts)]
  info.counts[IX > 1,COUNT := 0]
  info.counts.u = subset(info.counts, IX==1)
     
  ## Set up shorthand variables for use in rest of script
  aer = alignment.error.rate
  names(aer) = as.character(-(window):(window))
  umi.len = nchar(umi[1])
  mer = mismatch.error.rate
  w = window
  w.size = (2*w)+1
  
  ## Set up variables to save info after each iteration
  previous.umi = info.counts$inferred_umi
  previous.position = info.counts$inferred_position
  mer.list = list()
  aer.list = list()
  info.list = list(info.counts.u)

  
  iter = 1
  continue = TRUE
  while(continue == TRUE & iter <= max.iter) {
    
	next.umi = c()
	next.position = c()

    ## Set up mismatch and alignment error probability matrix
    probs = get.probability.matrix(umi.len, mer, aer)

	if(verbose) print(sprintf("Expectation ... Calculating most likely fragment memberships"))
	
	for(ix in 1:nrow(info.counts.u)) {
	  current.umi = as.character(info.counts.u[ix,inferred_umi])
	  current.position = info.counts.u[ix,inferred_position]
	  current.strand = info.counts.u[ix,strand]
	  current.chr = info.counts.u[ix,chr]
	
	  reads.ix = info.counts.u$chr == current.chr & info.counts.u$inferred_position > (current.position-w) & info.counts.u$inferred_position < (current.position+w) & info.counts.u$strand == current.strand

      ## Calculate UMI edit distance and genomic distance for each read within window
      reads.dist = stringdist(current.umi, info.counts.u[reads.ix,inferred_umi])
      if(current.strand == TRUE) {
	    reads.pos = info.counts.u$inferred_position[reads.ix] - current.position
      } else {
        reads.pos = -(info.counts.u$inferred_position[reads.ix] - current.position)
      }

      reads.prob = diag(probs[as.character(reads.dist), as.character(reads.pos), drop=FALSE])
      ll = log(reads.prob*info.counts.u$COUNT[reads.ix])

	  ties = ifelse(current.strand == TRUE, "left", "right")
	  ll.max.ix = get.max.ll(ll, ties=ties)

	  next.umi = c(next.umi, as.character(info.counts.u$inferred_umi[reads.ix][ll.max.ix]))
	  next.position = c(next.position, info.counts.u$inferred_position[reads.ix][ll.max.ix])
	}
 
	if(verbose) print(sprintf("Maximization ... Calculating new mismatch error rate"))
	mer.temp = calcMismatchErrorRate(info.counts.u$umi, next.umi, info.counts.u$COUNT)
	mer = mer.temp[1]
	mer.list = c(mer.list, list(mer.temp))
  
	if(verbose) print(sprintf("Maximization ... Calculating new alignment error rate"))
	aer.temp = calcAlignmentErrorRate(info.counts.u$position, next.position, info.counts.u$strand, info.counts.u$COUNT, window=w)
	aer = aer.temp[[1]]
	aer.list = c(aer.list, list(aer.temp))
  
    
    info.counts.u$inferred_umi = next.umi
    info.counts.u$inferred_position = next.position
    
    info.counts.u = info.counts.u[, `:=` (NEW_COUNT = sum(COUNT), NEW_IX=1:.N), by=list(chr, strand, inferred_position, inferred_umi)]
    info.counts.u[NEW_IX > 1,NEW_COUNT := 0]

    id = info.counts.u$ID
    info.counts$COUNT[id] = info.counts.u$NEW_COUNT
    info.counts$inferred_umi[id] = info.counts.u$inferred_umi
    info.counts$inferred_position[id] = info.counts.u$inferred_position
 
    info.counts.u = subset(info.counts.u, IX == 1)
    info.counts.u[,`:=`(IX=NEW_IX,COUNT=NEW_COUNT,NEW_COUNT=NULL,NEW_IX=NULL)]      
    
    # Update current UMI and position assignments
    info.list = c(info.list, list(info.counts.u))
  
	if(verbose) print(sprintf("Completed iteration %d", iter))
	iter = iter + 1
	
	total.diff = sum(!(info.counts$inferred_umi == previous.umi & info.counts$inferred_position == previous.position))
	if(verbose) print(sprintf("Total number of reads that changed fragment membership: %d", total.diff))
	if(total.diff == 0) {
	  continue = FALSE
	}
	previous.umi = info.counts$inferred_umi
	previous.position = info.counts$inferred_position
  }
  
  ## Make a final read matrix with duplicate infomation
  final.counts = info.counts
  final.counts[,`:=` (COUNT = NULL, IX = NULL)]
  same.fragment = final.counts$umi == final.counts$inferred_umi & final.counts$position == final.counts$position
  final.counts[,`:=` (same=same.fragment, duplicate = FALSE)]
  final.counts = final.counts[, `:=` (FINAL_GRP = .GRP, COUNT = .N, IX=1:.N), by=list(chr, strand, inferred_position, inferred_umi, same)]
  ind = final.counts$IX > 1 | !same.fragment
  final.counts$duplicate[ind] = TRUE 

  return(list(final.counts, list(info.list, info.counts, aer, mer, aer.list, mer.list, next.umi, next.position)))
}

calcMismatchErrorRate = function(actual.umi, inferred.umi, counts) {
  actual.umi = as.character(actual.umi)
  inferred.umi = as.character(inferred.umi)
  c = cbind(as.character(actual.umi), as.character(inferred.umi))
  l = unlist(lapply(1:nrow(c), function(i) stringDist(c(c[i,1], c[i,2]), method="hamming")[1]))
  
  mer.num = sum(l*counts)
  mer.denom = sum(counts*nchar(actual.umi[1]))
  new.mer =  mer.num / mer.denom
  
  return(c(new.mer, mer.num, mer.denom))
}


calcAlignmentErrorRate = function(actual.position, inferred.position, strand, counts, window=10, pseudo=1) {
  
  window.len = (2*window)+1
  new.aer = rep(0, window.len)
  names(new.aer) = (-window):(window)
  
  new.aer.forward = factor(inferred.position[strand] - actual.position[strand], levels=names(new.aer))
  new.aer.reverse = factor(-(inferred.position[!strand] - actual.position[!strand]), levels=names(new.aer))
  
  new.aer.forward.agg = aggregate(counts[strand], by=list(new.aer.forward), sum)
  new.aer.reverse.agg = aggregate(counts[!strand], by=list(new.aer.reverse), sum) 
  
  new.aer[new.aer.forward.agg[,1]] = new.aer[new.aer.forward.agg[,1]] + new.aer.forward.agg[,2]
  new.aer[new.aer.reverse.agg[,1]] = new.aer[new.aer.reverse.agg[,1]] + new.aer.reverse.agg[,2]
  
  new.aer.raw = new.aer
  new.aer = (new.aer + pseudo) / (sum(new.aer) + pseudo*window.len)
  
  return(list(new.aer, new.aer.raw))
}





get.max.ll = function(v, ties=c("left", "right")) {
  ties = match.arg(ties)

  v.max = max(v, na.rm=TRUE)
  v.max.ix = which(v == v.max)

  l = length(v.max.ix)+1
  ix = ifelse(ties == "left", floor(l/2), ceiling(l/2))
  return(v.max.ix[ix])
}



read.umi = function(qname, umi.flat) {
  if(length(umi.flag) == 1 & class(umi.flag) %in% c("numeric", "integer")) {
    umi = str_sub(qname, start=-umi.flag)
  } else if(length(umi.flag) == 2 & class(umi.flag) %in% c("numeric", "integer")) {  
    umi = str_sub(qname, start=umi.flag[1], end=umi.flag[2])
  } else if(length(umi.flag) == 1 & class(umi.flag) == "character") {
    umi = str_match(qname, pattern=umi.flag)[,2]
    if(sum(is.na(umi)) > 0) {
      stop(paste("UMIs not found for ", sum(is.na(umi)), " reads. Please check matching pattern", sep=""))
    }
  } else {
    stop("umi.flag not in a recognized format. See ?barcodeEM for details.")
  }
  return(umi)
}


get.probability.matrix = function(umi.len, mer, aer) {

  probs = do.call(cbind, lapply(1:length(aer), function(i) aer[i] * dbinom(x=0:umi.len, size=umi.len, prob=mer)))
  rownames(probs) = 0:umi.len
  colnames(probs) = names(aer)
  
  return(probs)
}

