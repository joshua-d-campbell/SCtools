scmap <- function (x, scale=TRUE, z.trim=c(-3,3), zero.val=-3.001, RowSideColors=NULL, ColSideColors=NULL, na.col="grey70",
                      row.clustering="unsupervised", col.clustering="unsupervised",
		      row.dendrogram=(row.clustering == "unsupervised"), col.dendrogram=(col.clustering=="unsupervised"),
		      hclustfun=hclust, distfun=dist, reorderfun = function(d, w) reorder(d, w),
		      margins = c(5, 5),
		      cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
		      labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
		      verbose = getOption("verbose"), col=c("grey90", "white", "red"), ...)
{

	# This is a replacement for the functions heatmap() and heatmap.plus().  It has a number of substantial improvements over the original versions, including:
	#      the ability to automatically resize row and column labels to fit the specified margins
	#      the use of split.screen() instead of layout(), which allows multiple heatmaps to be drawn on a single device
	#      parameters to suppress the plotting of row and column dendrograms
	#      a parameter for applying z-score cutoffs to the z-normalized heatmap matrix
	#
	# Written by Adam Gower, 2008-2009
	# 
	# INPUT
	# x               numeric matrix that will be plotted as heatmap
	# scale           one of either "row", "column", or "none"; specifies the direction in which z-normalization will be applied to x
	# z.trim          two-element vector that specifies the lower and upper z-score cutoffs to apply to the z-normalized matrix
	# RowSideColors   either a vector of length nrow(x) or a matrix with nrow(x) rows, of type "character"; specifies colors for sidebar labeling
	# ColSideColors   either a vector of length ncol(x) or a matrix with ncol(x) columns, of type "character"; specifies colors for sidebar labeling
	# row.clustering  one of either "supervised", "unsupervised", or "semisupervised"; specifies method for clustering rows of heatmap
	# col.clustering  one of either "supervised", "unsupervised", or "semisupervised"; specifies method for clustering columns of heatmap
	# row.dendrogram  logical value; specifies whether or not to draw the dendrogram of row clustering (suppressed if row.clustering is not "unsupervised")
	# col.dendrogram  logical value; specifies whether or not to draw the dendrogram of row clustering (suppressed if row.clustering is not "unsupervised")
	# hclustfun       specifies the hierarchical clustering function to be used
	# distfun         specifies the distance function to be used for hierarchical clustering
	# reorderfun      not used; only added for backwards compatibility
	# margins         two-element vector that specifies the margins for column and row names, respectively
	# cexRow          character expansion (cex) factor for row labels
	# cexCol          character expansion (cex) factor for column labels
	# labRow          character vector specifying row labels
	# labCol          character vector specifying column labels
	# main            character string specifying heatmap title
	# xlab            character string specifying x axis label
	# ylab            character string specifying y axis label
	# keep.dendro     logical value; specifies whether or not to return the information describing the dendrograms
	# verbose         logical value; specifies whether to print messages during plotting
	# col             character string specifying the color palette to use

	compute.cex <- function (label, width, height, units="inches", cex.lim=1) {
		# Adam Gower, 2009
		# Computes the cex factor needed to place a label (or the largest in a vector of labels) in a given plot region
	
		# label:  string of text to be drawn, or a character vector of text strings to be drawn, the largest of which will be used for computations
		# width:  width of the plot region in which to draw the text
		# height: height of the plot region in which to draw the text
		# units:  units in which width and height are reported
	
		# Determine the aspect ratio (W:H) of the plot region
		region.aspect <- width / height;
	
		# Get the width, height, and aspect ratio of the (largest) label
		label.width <- max(strwidth(label, units=units));
		label.height <- max(strheight(label, units=units));
		label.aspect <- label.width / label.height;
	
		# Compute the cex factor that will place the text completely within the plot region 
		cex.label <- ifelse(label.aspect > region.aspect, width/label.width, height/label.height);
		# If there is a maximum cex limit specified, apply it
		if (is.numeric(cex.lim)) cex.label <- min(cex.label, cex.lim);
	
		return(cex.label);
	}

	semisupervised.clustering <- function (x, col.class) {
		# Function to reorder the columns of a matrix using hierarchical clustering within two specified classes
		# Adam Gower, 2008
		#
		# x:         a matrix that has samples in columns and observations (e.g., probesets or genes) in rows
		# col.class: an integer vector, of length ncol(x), that denotes the classes of the samples;
	
		classes <- unique(col.class);
		
		# The dist function works on the rows, not on the columns, so we need to transpose x first
		tx <- t(x);
		# Use the functions hclust and dist to perform hierarchical clustering, once within each class, and concatenate the result
		col.order <- c();
		for (i in 1:length(classes)) {
			class.indices <- which(col.class == classes[i]); 
			col.order <- c(col.order, if (length(class.indices) > 1) class.indices[hclust(dist(tx[class.indices,]))$order] else class.indices);
		}
		return(col.order);
	}

	##### FIX-IT LIST
	# Automatically set margins according to maximum length of labRow, labCol (put after they are computed) based on some reasonable cex?
	# Also, automatically set title margin according to number of lines
	# i.e., length(unlist(strsplit("foo\nbar\nbaz","\n"))) = 3
	# And, automatically set cex.main according to length of maximum line, i.e., max(sapply(unlist(strsplit("foobar\nbarbarian\nbazilisk","\n")), nchar)) = 9
	#    and the size of the window

	# Check for correct margins argument
	if (!is.numeric(margins) || length(margins) != 2) stop("'margins' must be a numeric vector of length 2");

	# Get dimensions of matrix x and check for proper size
	nr <- nrow(x); nc <- ncol(x);
#	if (nr < 2 || nc < 2) stop("'x' must have at least 2 rows and 2 columns");

	# Check to make sure that RowSideColors and ColSideColors arguments are the correct size
	if (!is.null(RowSideColors) & (NROW(RowSideColors) != nr)) stop("'RowSideColors' must be a vector of length nrow(x) or a matrix with nrow(x) rows.");
	if (!is.null(ColSideColors) & (NROW(ColSideColors) != nc)) stop("'ColSideColors' must be a vector of length ncol(x) or a matrix with ncol(x) rows.");

	# Just in case, set the dendrograms accordingly if none are generated
	if (row.clustering != "unsupervised") row.dendrogram <- FALSE;
	if (col.clustering != "unsupervised") col.dendrogram <- FALSE;

	############################################## CREATE HEATMAP LAYOUT ##############################################
	# Create zeroed-out layout matrix with enough rows and columns to hold sidebars and dendrograms
	lrows <- 1+(!is.null(ColSideColors))+col.dendrogram+!is.null(main);
	lcols <- 1+(!is.null(RowSideColors))+row.dendrogram; 
	lmat <- matrix(0, nrow=lrows, ncol=lcols);
	# Plot heatmap of length 4 and width 4 first
	lmat[lrows, lcols] <- panel <- 1;
	lwid <- lhei <- c(4);
	# Add row side color bar of width 0.2 if specified
	if (!is.null(RowSideColors)) {
		lmat[lrows, lcols-1] <- panel <- panel + 1;
		lwid <- c(0.2, lwid);
	}
	# Add column side color bar of height 0.2 if specified
	if (!is.null(ColSideColors)) {
		lmat[lrows-1, lcols] <- panel <- panel + 1;
		lhei <- c(0.2, lhei);
	}
	# Add row dendrogram window of width 1 if specified
	if (row.dendrogram) {
		lmat[lrows, lcols-1-(!is.null(RowSideColors))] <- panel <- panel + 1;
		lwid <- c(1, lwid);
	}
	# Add column dendrogram window of height 1 if specified
	if (col.dendrogram) {
		lmat[lrows-1-(!is.null(ColSideColors)), lcols] <- panel <- panel + 1;
		lhei <- c(1, lhei);
	}
	# Add title window of height 0.5 if specified
	if (!is.null(main)) {
		lmat[1, lcols] <- panel <- panel + 1;
		lhei <- c(0.5, lhei);
	}
	if (verbose) {
		cat("layout: widths = ", lwid, ", heights = ", lhei, "\nlmat=\n");
		print(lmat);
	}

	heatmap.figs <- matrix(nrow=max(lmat), ncol=4);
	lwid <- lwid / sum(lwid);
	lhei <- lhei / sum(lhei);
	cumsum.lwid <- cumsum(c(0,lwid));
	cumsum.lhei <- 1 - cumsum(c(0,lhei));
	for (i in 1:lrows) {
		for (j in 1:lcols) {
			if (lmat[i,j])
				heatmap.figs[lmat[i,j],] <- c(left=cumsum.lwid[j], right=cumsum.lwid[j+1], bottom=cumsum.lhei[i+1], top=cumsum.lhei[i]);
		}
	}
	heatmap.screens <- split.screen(figs=heatmap.figs);

	# Normalize heatmap if requested, trimming at specified z-value cutoffs and setting zeros to specified values
	if (scale == TRUE) {
      x = t(apply(x, 1, sc.scale.cluster, cutoff=z.trim, zero.val=zero.val))
	} 

	# Get order of row indices based on clustering type specified
	rowInd <- switch(row.clustering, "semisupervised" = {
					 	if (NCOL(RowSideColors) == 1) {
							semisupervised.clustering(t(x), drop(RowSideColors))
						} else {
							# stop ("Cannot use semisupervised row clustering if RowSideColors is a matrix with > 1 column");
							semisupervised.clustering(t(x), RowSideColors[,1]);
						}
					 },
				     	 "supervised"     = 1:nr,
					 "unsupervised"   = {
						 ddr <- as.dendrogram(hclustfun(distfun(x)));
						 order.dendrogram(ddr);
					 });
	# Get order of column indices based on clustering type specified
	colInd <- switch(col.clustering, "semisupervised" = {
					 	if (NCOL(ColSideColors) == 1) {
							semisupervised.clustering(x, drop(ColSideColors))
						} else {
							# stop ("Cannot use semisupervised column clustering if ColSideColors is a matrix with > 1 column");
							semisupervised.clustering(x, ColSideColors[,1]);
						}
					 },
					 "supervised"     = 1:nc,
					 "unsupervised"   = {
						 ddc <- as.dendrogram(hclustfun(distfun(t(x))));
						 order.dendrogram(ddc);
					 });

	
	# Reorder the heatmap matrix according to the row and column indices
	x <- x[rowInd, colInd];
	
	# If row labels were not specified, create them from the row names of matrix x; otherwise, reorder the labels specified in labRow
	labRow <- if (is.null(labRow)) {
	      		if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x);
		  } else labRow[rowInd];
	# If column labels were not specified, create them from the column names of matrix x; otherwise, reorder the labels specified in labCol
        labCol <- if (is.null(labCol)) {
		  	if (is.null(colnames(x))) (1:nc)[colInd] else colnames(x);
                  } else labCol[colInd];

	############################################## PLOT HEATMAP ELEMENTS ##############################################

	# Set margins for heatmap, leaving bottom margin for column labels and right margin for row labels
	par(mar = c(margins[1], 0, 0, margins[2]));
	frame();
	# Plot main heatmap
	#image(1:nc, 1:nr, t(x), xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col=col, ...);
	temp.x = x
	temp.x[x == zero.val] = NA
	x.col = image.sc(temp.x, col, na.col=na.col)
	
	# Plot the column labels and any x axis label that was supplied
	if (all(!is.na(labCol))) {
		cexCol <- 0.95 * compute.cex(width = (margins[1] - 0.5) * par("cin")[1], height = par("fin")[1] / nc, label = labCol);
		axis(1, 1:nc, labels = rev(labCol), las = 2, line = -0.5, tick = 0, cex.axis = cexCol);
	}
	if (!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25);
	# Plot the row labels and any y axis label that was supplied 
	if (all(!is.na(labRow))) {
		cexRow <- 0.95 * compute.cex(width = (margins[2] - 0.5)* par("cin")[1], height = par("fin")[2] / nr, label = labRow);
		axis(4, 1:nr, labels = rev(labRow), las = 2, line = -0.5, tick = 0, cex.axis = cexRow);
	}
	if (!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25);
	close.screen(screen());

	# Plot row side colors if they were supplied
	if (!is.null(RowSideColors)) {
		# Coerce row side colors to matrix
		RowSideColors <- as.matrix(RowSideColors);
		# Get any track names that were specified
		row.tracks <- colnames(RowSideColors);
		# Convert row side colors to factor (in case color names were specified)
		RowSideColors <- as.factor(RowSideColors[rowInd[nr:1],]);
		# Set margins, leaving bottom margin for column labels and 10% left and right margins
		par(mai = c(margins[1]*par("csi"), 0.1*par("fin")[1], 0, 0.1*par("fin")[1]));
		# Plot row side colors; if color names were specified in the RowSideColors argument, they will be passed as 'col' argument
		image(1:(length(RowSideColors)/nr), 1:nr, t(matrix(as.numeric(RowSideColors), nrow=nr)), col=levels(RowSideColors), axes = FALSE, xlab=NA, ylab=NA);
		# If track names were specified, write them
		if (!is.null(row.tracks)) {
			par(cex.axis = 0.95 * compute.cex(width = (margins[1] - 0.5)* par("cin")[1], height = 0.8 * par("fin")[1] / length(row.tracks), label = row.tracks));
			axis(side=1, at=1:length(row.tracks), labels=row.tracks, line=-0.5, tick=FALSE, las=2);
		}
		close.screen(screen());
	}
	# Plot column side colors if they were supplied
	if (!is.null(ColSideColors)) {
		# Coerce column side colors to matrix
		ColSideColors <- as.matrix(ColSideColors);
		# Get any track names that were specified
		col.tracks <- colnames(ColSideColors);
		# Convert column side colors to factor (in case color names were specified)
		ColSideColors <- as.factor(ColSideColors[colInd[nc:1],]);
		# Set margins, leaving right margin for row labels and 10% bottom and top margins
		par(mai = c(0.1*par("fin")[2], 0, 0.1*par("fin")[2], margins[2]*par("csi")));
		# Plot column side colors; if color names were specified in the RowSideColors argument, they will be passed as 'col' argument
		image(1:nc, 1:(length(ColSideColors)/nc), matrix(as.numeric(ColSideColors), nrow=nc), col=levels(ColSideColors), axes = FALSE, xlab=NA, ylab=NA);
		# If track names were specified, write them
		if (!is.null(col.tracks)) {
			par(cex.axis = 0.95 * compute.cex(width = (margins[2] - 0.5) * par("cin")[1], height = 0.8 * par("fin")[2] / length(col.tracks), label = col.tracks));
			cat(sprintf("cex.axis for col.tracks is %g\n", par("cex.axis"))); 
			axis(side=4, at=1:length(col.tracks), labels=col.tracks, line=-0.5, tick=FALSE, las=2);
		}
		close.screen(screen());
	}

	# Plot row dendrogram if requested
	if (row.dendrogram) {
		# Set margins for row dendrogram panel, leaving bottom margin for column labels
		par(mar = c(margins[1], 0, 0, 0));
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none");
		close.screen(screen());
	}

	# Plot column dendrogram if requested
	if (col.dendrogram) {
		# Set margins for column dendrogram panel, leaving right margin for row labels
		par(mar = c(0, 0, 0, margins[2]));
		plot(ddc, horiz = FALSE, axes = FALSE, xaxs = "i", leaflab = "none");
		close.screen(screen());
	}
	
	# Add main title if requested
	if (!is.null(main)) {
		par(mar=c(0,0,0,margins[2]));
		frame();
		cex.main <- (1/1.04)*compute.cex(width = par("fin")[1] - (margins[2] * par("csi")), height = par("fin")[2], label = main);
		text(0.5, 0.5, main, font=1, cex=cex.main);
		close.screen(screen());
	}

	# Return list of row and column indices, and dendrograms, if requested
	invisible(list(matrix=x, matrix.col=x.col, rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && row.dendrogram) ddr, Colv = if (keep.dendro && col.dendrogram) ddc))
}




image.sc = function(m, col.str, na.col, na.val) {
  col = colorRampPalette(col.str)(100)
  
  nc = ncol(m)
  nr = nrow(m)
  
  m.min = min(m[!is.na(m)])
  m.max = max(m[!is.na(m)])

  m.col = matrix(NA, nrow=nr, ncol=nc)
  col.ix = seq(m.min, m.max, length=length(col))
  br = abs((abs(col.ix[2]) - abs(col.ix[1])) / 2)
  
  for(i in 1:length(col.ix)) {
    m.dist = abs(m-col.ix[i])
    m.col[m.dist < br] = col[i]
  }
  m.col[is.na(m.col)] = na.col
  
  plot(0, xlim=c(1,nc), ylim=c(1,nr), type="n", axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")
  rect(rep((nc:1)-0.5, nr), rep((nr:1)-0.5, each=nc), rep((nc:1)+0.5, nr), rep((nr:1)+0.5, each=nc), col=t(m.col), border=NA)
  invisible(return(m.col))
}


sc.scale.cluster = function(r, cutoff=c(-3,3), zero.val=-3.001) {
  k = r > 0
  r[k] = (r[k] - mean(r[k])) / sd(r[k])
  i = r[k] < cutoff[1]
  r[k][i] = cutoff[1]
  i = r[k] > cutoff[2]
  r[k][i] = cutoff[2]
  
  r[!k] = zero.val
  return(r)
}


