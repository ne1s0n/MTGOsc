# Copyright 2019 Nelson Nazzicari
# This file is part of MTGOsc
#
# MTGOsc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MTGOsc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MTGOsc If not, see <http://www.gnu.org/licenses/>.

#' Compute and save a coexpression matrix
#'
#' This function uses the passed coexpression function (default: Pearson's correlation \link[stats]{cor})
#' to compute a coexpression matrix between each pair of genes expression arrays. The matrix is
#' saved in the passed directory.
#'
#' @param geneExpression a Seurat dgCMatrix (from the @data field in a Seurat object), containing the expression data.
#' @param outfolder the data folder where to save the results
#' @param overwrite boolean. Should an existing coexpression file be overwritten on disk?
#' @param fun the coexpression function, it will be fed two genes' expression arrays. Defaults to \link[stats]{cor}
#' @param ... all extra parameters are passed to coexpression function
#'
#' @return the coexpression matrix, in long form
#' @export
write.coexpressionMatrix = function(geneExpression, outfolder, overwrite = FALSE, fun = cor, ...){
  fn = get.filenames(outfolder)

  if (file.exists(fn$coexpression.filename) & !overwrite){
    stop(paste('File', fn$coexpression.filename, 'already exists but *overwrite* is set to FALSE'))
  }

  #if we get here we can proceed
  mycols = c('gene1', 'gene2', 'coexpr')
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

  #computing coexpression as square, symmetric matrix
  args = list(
    as.matrix(geneExpression),
    ...
  )
  coexpr.square = do.call(fun, args)

  #if it's not square and symmetric it means we had problems
  if (!isSymmetric(coexpr.square)){
    stop('Prolems with selected coexpression function: returned matrix is not square and symmetrical')
  }

  #if we get here we can pass from wide to long form and
  #remove the duplicate values
  coexpr.square[upper.tri(coexpr.square, diag = TRUE)] = NA
  coexpr.long = reshape2::melt(coexpr.square, na.rm = TRUE)
  colnames(coexpr.long) = mycols

  #we can now save and return
  write.table(file = fn$coexpression.filename, x = coexpr.long, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

  return(coexpr.long)

  #BELOW: old implementation, to be removed as soon as the
  #current one is tested
  #
  # #room for returning the result
  # l = ngenes * (ngenes-1) / 2 - 1
  # g1 = array(dim=l)
  # g2 = array(dim=l)
  # coexpr = array(dim=l)
  # cnt = 0
  #
  # #computing for each possible pair of genes the correlation, and writing right
  # #away in the output file
  # for (gene1 in 1:(ngenes-1)){
  #   for(gene2 in (gene1+1):ngenes){
  #     #computing correlation between two genes
  #     args = list(
  #       geneExpression@data[gene1,],
  #       geneExpression@data[gene2,],
  #       ...
  #     )
  #     res.curr = do.call(fun, args)
  #
  #     #removing NAs, indicating those genes that never vary
  #     if(is.na(res.curr)){
  #       next
  #     }
  #
  #     #ready to save
  #     cnt = cnt + 1
  #     g1[cnt] = rownames(geneExpression@data)[gene1]
  #     g2[cnt] = rownames(geneExpression@data)[gene2]
  #     coexpr[cnt] = res.curr
  #   }
  # }
  #
  # #we can now save and return
  # coexpr.long = data.frame(
  #   gene1=g1[1:cnt],
  #   gene2=g2[1:cnt],
  #   coexpr=coexpr[1:cnt])
  #
  # #ready to save
  # write.table(file = fn$coexpression.filename, x = coexpr.long, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  #
  # return(coexpr.long)
}

#' Compute and save network edges
#'
#' This function takes a network (in the form of a coexpression matrix) and trims some of the
#' edges using the passed function \code{fun}. The resulting network is both returned and
#' saved on disk.
#'
#' @param coexpression a coexpression matrix, in long form, as the one returned by \link{write.coexpressionMatrix}
#' @param outfolder the data folder where to save the results
#' @param keep.weights boolean. If TRUE arc weights are kept. Defaults to FALSE.
#' @param overwrite boolean. Should an existing edge file be overwritten on disk?
#' @param fun edge filtering function, defaults to \link{abs_threshold}
#' @param ... extra arguments are passed to edge filtering function
#'
#' @return a trimmed network
#' @export
write.edges = function(coexpression, outfolder, keep.weights = TRUE, overwrite = TRUE, fun = abs_threshold, ...){
  fn = get.filenames(outfolder)

  if (file.exists(fn$edges.filename) & !overwrite){
    stop(paste('File', fn$edges.filename, 'already exists but *overwrite* is set to FALSE'))
  }

  #if the user passed a filename, we read the data in it
  if (is.character(coexpression)){
    coexpr.filename = coexpression
    coexpression = read.table(file = coexpression, header = FALSE, sep='\t')
    if(!all(colnames(coexpression) == c("gene1", "gene2", "coexpr"))){
      stop(paste('coexpression file', coexpr.filename, 'not in the right format'))
    }
  }

  #we take the absolute value of coexpression
  coexpression$coexpr = abs(coexpression$coexpr)

  #subsetting using the passed function, if any
  if (!is.null(fun)){
    args = list(coexpression, ...)
    coexpression = do.call(fun, args)
  }

  #should we keep the weights?
  if (!keep.weights){
    coexpression = coexpression[,c("gene1", "gene2")]
  }

  #ready to save edges
  write.table(file = fn$edges.filename, x = coexpression, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

  #return the edges anyway
  return(coexpression)
}

#' Writes a MTGO parameter file on disk
#'
#' This function creates and saves a MTGO parameter file on disk.
#'
#' @param outfolder the data folder where to save the param file
#' @param overwrite boolean. Should an existing parameter file be overwritten on disk?
#' @param MinSize MTGO MinSize parameter
#' @param MaxSize MTGO MaxSize parameter
#'
#' @return nothing
#' @export
write.paramFile = function(outfolder, overwrite = FALSE, MinSize = 3, MaxSize = 30){
  #file management
  fn = get.filenames(outfolder)
  if (file.exists(fn$param.filename) & !overwrite){
    stop(paste('File', fn$param.filename, 'already exists but *overwrite* is set to FALSE'))
  }
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

  #MTGO expexts a trailing slash when a folder is inputed
  outfolder.slash = paste(sep='', outfolder, .Platform$file.sep)

  #if we get here we can proceed
  fin = file(description = fn$param.filename, open = 'w')
  writeLines(con = fin, paste(sep='\t', 'PathEdge', fn$edges.filename))
  writeLines(con = fin, paste(sep='\t', 'PathGO', fn$GO.filename))
  writeLines(con = fin, paste(sep='\t', 'PathResult', outfolder.slash))
  writeLines(con = fin, paste(sep='\t', 'MinSize', MinSize))
  writeLines(con = fin, paste(sep='\t', 'MaxSize', MaxSize))
  close(fin)
}

#' Calls for MTGO execution
#'
#' This function invokes MTGO. You need to have created a MTGO parameter file (see \link{write.paramFile}),
#' an edge file (see \link{write.edges}) and a dictionary file (see \link{write.dictionary}). Files formats and
#' names are standard so it's better
#' to use the suggested functions to create them.
#'
#' @param outfolder the folder where to read data and write results
#' @param verbose if TRUE errors and warnings are printed after the run completes
#'
#' @return nothing
#' @export
call.MTGO = function(outfolder, verbose = TRUE){
  #checking for correct java version
  msg = verify.java()
  if (!is.null(msg)){
    stop(msg)
  }

  fn = get.filenames(outfolder)

  #executing MTGO, capturing output
  system2(
    stdout = fn$MTGO.res.out, stderr = fn$MTGO.res.err,
    command = 'java',
    args = c('-jar', shQuote(fn$MTGO.jar), shQuote(fn$param.filename))
  )

  if(verbose){
    writeLines('======> MTGO executed with the following OUTPUT:')
    writeLines(readLines(fn$MTGO.res.out))
    writeLines('\n======> MTGO executed with the following ERRORS:')
    writeLines(readLines(fn$MTGO.res.err))
  }
}


#' Create and save a dictionary file
#'
#' This function takes to sets (genes and corresponding terms) and saves a dictionary
#' file on disk.
#'
#' @param genes character array, list of genes
#' @param terms character array, a list of terms associated with the genes (same order, same size)
#' @param outfolder the data folder where to save the dictionary file
#' @param overwrite boolean. Should an existing dictionary file be overwritten on disk?
#'
#' @return a data frame with genes and terms as columns
#' @export
write.dictionary = function(genes, terms, outfolder, overwrite = FALSE){
  if (length(genes) != length(terms)){
    stop('genes array and terms array must be of the same length')
  }

  fn = get.filenames(outfolder)
  if (file.exists(fn$GO.filename) & !overwrite){
    stop(paste('File', fn$GO.filename, 'already exists but *overwrite* is set to FALSE'))
  }
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

  df = data.frame(g = genes, t = terms)
  write.table(df, file = fn$GO.filename, row.names = FALSE, col.names = FALSE, sep='\t', quote = FALSE)

  return(df)
}

#' List standard filenames
#'
#' All MTGOsc files have standard naming. This folders returns a list of them, starting from
#' a root folder.
#'
#' @param outfolder the root folder to be used.
#'
#' @return a list with all standard names listed as fields
#' @export
get.filenames = function(outfolder){
  res = list()
  res$MTGO.jar = file.path(find.package('MTGOsc'), 'java', 'MTGO.jar')
  res$coexpression.filename = file.path(outfolder, 'coexpression_full.tsv')
  res$edges.filename = file.path(outfolder, 'edges.tsv')
  res$gene_list.filename = file.path(outfolder, 'genes_list.txt')
  res$param.filename = file.path(outfolder, 'params.txt')
  res$GO.filename = file.path(outfolder, 'GO.txt')
  res$MTGO.res.out = file.path(outfolder, 'MTGO.out.txt')
  res$MTGO.res.err = file.path(outfolder, 'MTGO.err.txt')
  return(res)
}

#' Subset a coexpression network for thresholded absolute values
#'
#' This function removes from the passed coexpression network \code{x} all edges that,
#' in absolute value, are less of the passed \code{threshold}.
#'
#' @param x a coexpression matrix, as returned by \link{write.coexpressionMatrix}
#' @param threshold the value to be used for subsetting
#'
#' @return a smaller version of the original matrix
#' @export
abs_threshold = function(x, threshold = 0.5){
  sel = abs(x$coexpr) >= threshold
  return(x[sel,])
}

#' Subset a coexpression network to maximise free scale fit
#'
#' This function removes edges from the passed network so that the resulting subnetwork
#' is maximally scale free. First, for each value in \code{thresholds} argument, a subnetwork
#' is created. Each subnetwork is obtained from the original \code{x} removing all edges
#' smaller, in absolute value, than the corresponding threshold.
#' Then each subnetwork is fitted a power law, and the one whose gamma is closest
#' to \code{target.gamma} is selected.
#'
#' @param x a coexpression network, as returned by \link{write.coexpressionMatrix}
#' @param thresholds numeric array of thresholds to be considered
#' @param target.gamma for biological networks target gamma is somewhere between 2 and 3
#' @param verbose should extra informations be printed?
#'
#' @return a smaller version of the original matrix
#' @export
thinning_scale_free = function(x, thresholds = seq(from=0.1, to=0.9, by=0.1), target.gamma = 2, verbose=TRUE){
  best.gamma = NULL
  best.delta = Inf
  best.network = NULL
  best.threshold = NULL
  for (threshold in thresholds){
    #subsetting the network to current threshold
    x.curr = x[x$coexpr >= threshold,]

    if (nrow(x.curr) == 0){
      #we don't have any more data :(
      next()
    }

    #creating a graph from igraph package
    gra = igraph::graph_from_edgelist(as.matrix(x.curr[,c(1,2)]), directed = FALSE)
    d = igraph::degree(gra)
    fit = igraph::fit_power_law(d)
    gamma = fit$alpha
    delta = abs(gamma - target.gamma)

    #do we have a new best?
    if (delta < best.delta){
      best.delta = delta
      best.gamma = gamma
      best.network = x.curr
      best.threshold = threshold
    }
  }

  #should we tell the user about the choice of threshold?
  if(verbose){
    writeLines(paste('gamma:', best.gamma, paste(sep='', '(target:', target.gamma, ')')))
    writeLines(paste('threshold:', best.threshold))
    writeLines(paste('network edges:', nrow(best.network)))
  }

  return(best.network)
}

#' Subset a coexpression network to the desired percentile
#'
#' This function subsets the edges from the passed network and keeps only
#' those comprised in the desired top percentile.
#'
#' @param x a coexpression network, as returned by \link{write.coexpressionMatrix}
#' @param top_percentile the desired percentile (e.g. 0.1 means "keep only top 10%")
#'
#' @return a smaller version of the original matrix
#' @export
thinning_percentile = function(x, top_percentile = 0.1){
  threshold = quantile(x$coexpr, 1 - top_percentile)
  return(subset(x, coexpr >= threshold))
}

#' Verifies if the correct Java version is installed on the system
#'
#' This function verifies that the correct version of Java (TM) is installed on
#' the system. MTGO requires the Java (TM) Runtime Environment (version >= 1.8)
#' released by Oracle to run.
#'
#'
#' @return NULL if the right Java is present, or a character vector containing an error message if
#' something is wrong
#' @export
verify.java = function(){
  #executing MTGO, capturing output
  res = system2(
    stdout = TRUE, stderr = TRUE,
    command = 'java',
    args = c('-version')
  )

  all.good = TRUE

  #Java the original (Java(TM))
  all.good = all.good & grepl(pattern = 'Java\\(TM\\)', res[2])

  #at least version 1.8
  tmp = gsub(pattern = 'java version "', replacement = '', res[1])
  pieces = strsplit(tmp, split = '.', fixed = TRUE)[[1]]
  if (length(pieces) < 2){
    all.good = FALSE
  }else{
    version = as.numeric(pieces[1]) + as.numeric(pieces[2])/10
    all.good = all.good & version >= 1.8
  }

  if (all.good){
    return(NULL)
  }else{
    return('Java(TM) Runtime Environment version 1.8 or more recent is required to run MTGO.\n
Openjdk is not fine, you need to install the original one from Oracle.')
  }
}

#' Compute coexpression via propr::propr
#'
#' This function is a wrapper around \link{propr::propr} to compute coexpression starting from
#' a Seurat dgCMatrix (from the @data field in a Seurat object), containing the expression data.
#' This function returns only the @matrix field of the propr object.
#'
#' @param geneExpression a Seurat dgCMatrix (from the @data field in a Seurat object)
#' @param verbose if TRUE the messages from \link{propr::propr} are not suppressed
#' @param ... extra parameters are passed to \link{propr::propr}
#'
#' @return the @matrix field of a propr object
#' @export
coexpr_propr = function(geneExpression, verbose=FALSE, ...){
  if (verbose){
    #let's keep the messages
    res = propr::propr(t(as.matrix(geneExpression)), ...)
  }else{
    #suppressed messages for better logs
    suppressMessages({res = propr::propr(t(as.matrix(geneExpression)), ...)})
  }
  return(res@matrix)
}


