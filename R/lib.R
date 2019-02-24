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
#' @param SeuratObject the Seurat object containing the expression data
#' @param outfolder the data folder where to save the results
#' @param outfile.name optional filename for coexpression file
#' @param overwrite boolean. Should an existing coexpression file be overwritten on disk?
#' @param fun the coexpression function, it will be fed two genes' expression arrays. Defaults to \link[stats]{cor}
#' @param ... all extra parameters are passed to coexpression function
#'
#' @return the coexpression matrix, in long form
#' @export
write.coexpressionMatrix = function(SeuratObject, outfile.name = 'coexpression_full.tsv', overwrite = FALSE, fun = cor, ...){
  outfolder = dirname(outfile.name)
  fn = get.filenames(outfolder)

  if (!is.null(outfile.name)){
    fn$coexpression.filename = file.path(outfolder, outfile.name)
  }

  if (file.exists(fn$coexpression.filename) & !overwrite){
    stop(paste('File', fn$coexpression.filename, 'already exists but *overwrite* is set to FALSE'))
  }

  #if we get here we can proceed
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)
  fp = file(fn$coexpression.filename, open = 'w')
  writeLines(con=fp, paste(sep='\t', 'gene1', 'gene2', 'coexpr'))

  #number of genes
  ngenes = length(rownames(SeuratObject@data))

  #room for returning the result
  l = ngenes * (ngenes-1) / 2 - 1
  g1 = array(dim=l)
  g2 = array(dim=l)
  coexpr = array(dim=l)
  cnt = 0

  #computing for each possible pair of genes the correlation, and writing right
  #away in the output file
  for (gene1 in 1:(ngenes-1)){
    for(gene2 in (gene1+1):ngenes){
      #computing correlation between two genes
      args = list(
        SeuratObject@data[gene1,],
        SeuratObject@data[gene2,],
        ...
      )
      res.curr = do.call(fun, args)

      #removing NAs, indicating those genes that never vary
      if(is.na(res.curr)){
        next
      }

      #ready to save
      cnt = cnt + 1
      g1[cnt] = rownames(SeuratObject@data)[gene1]
      g2[cnt] = rownames(SeuratObject@data)[gene2]
      coexpr[cnt] = res.curr
    }
  }

  #we can now save and return
  coexpr.long = data.frame(
    gene1=g1[1:cnt],
    gene2=g2[1:cnt],
    coexpr=coexpr[1:cnt])

  #ready to save
  write.table(file = fn$coexpression.filename, x = coexpr.long, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)

  return(coexpr.long)
}


#' Compute and save network edges
#'
#' This function uses the passed function ()
#'
#' @param coexpression a coexpression matrix, in long form, as the one returned by \link{write.coexpressionMatrix}
#' @param outfolder the data folder where to save the results
#' @param outfile.name optional filename for coexpression file
#' @param keep.weights boolean. If TRUE arc weights are kept. Defaults to FALSE.
#' @param overwrite boolean. Should an existing edge file be overwritten on disk?
#' @param fun edge filtering function, defaults to \link{abs.threshold}
#' @param ... extra arguments are passed to edge filtering function
#'
#' @return f
#' @export
write.edges = function(coexpression, outfolder, keep.weights = TRUE, outfile.name = NULL, overwrite = TRUE, fun = abs.threshold, ...){
  fn = get.filenames(outfolder)

  if (!is.null(outfile.name)){
    fn$edges.filename = file.path(outfolder, outfile.name)
  }

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

write.paramFile = function(outfolder, overwrite = FALSE, MinSize = 3, MaxSize = 30){
  #file management
  fn = get.filenames(outfolder)
  if (file.exists(fn$param.filename) & !overwrite){
    stop(paste('File', fn$param.filename, 'already exists but *overwrite* is set to FALSE'))
  }
  dir.create(outfolder, showWarnings = FALSE, recursive = TRUE)

  #if we get here we can proceed
  fin = file(description = fn$param.filename, open = 'w')
  writeLines(con = fin, paste(sep='\t', 'PathEdge', fn$edges.filename))
  writeLines(con = fin, paste(sep='\t', 'PathGO', fn$GO.filename))
  writeLines(con = fin, paste(sep='\t', 'PathResult', outfolder))
  writeLines(con = fin, paste(sep='\t', 'MinSize', MinSize))
  writeLines(con = fin, paste(sep='\t', 'MaxSize', MaxSize))
  close(fin)
}

call.MTGO = function(outfolder, verbose = TRUE){
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

get.filenames = function(outfolder){
  res = list()
  res$Seurat2MTGO.folder = '/home/nelson/research/Seurat2MTGO' #this need to change to a more packet style approach
  res$MTGO.jar = file.path(res$Seurat2MTGO.folder, 'inst', 'java', 'MTGO.jar')
  res$coexpression.filename = file.path(outfolder, 'coexpression_full.tsv')
  res$edges.filename = file.path(outfolder, 'edges.tsv')
  res$gene_list.filename = file.path(outfolder, 'genes_list.txt')
  res$param.filename = file.path(outfolder, 'params.txt')
  res$GO.filename = file.path(outfolder, 'GO.txt')
  res$MTGO.res.out = file.path(outfolder, 'MTGO.out.txt')
  res$MTGO.res.err = file.path(outfolder, 'MTGO.err.txt')
  return(res)
}

cor.threshold = function(x, y, threshold = NA, use = "everything", method = c("pearson", "kendall", "spearman")){
  #if threshold is NA we just fall back to old correlation function
  #otherwise, we remove all values which are, for both arrays and in absolute value
  #less then passed threshold

  if (!is.na(threshold)){
    #in this case we remove some samples
    x.bad = abs(x) <= threshold
    y.bad = abs(y) <= threshold
    bad = x.bad & y.bad
    x = x[!bad]
    y = y[!bad]
  }

  #ready to run
  return(cor(x=x, y=y, use=use, method=method))
}

#' Subset data for thresholded absolute values
#'
#' This function removes from the passed data \code{x} all entried that, in absolute value, are
#' less of the passed \code{threshold}.
#'
#' @param x the data to be subset
#' @param threshold the value to be used for subsetting
#'
#' @return a smaller version of the original data
#' @export
abs.threshold = function(x, threshold = 0.5){
  sel = abs(x$coexpr) >= threshold
  return(x[sel,])
}

top.values = function(x, percentile = 0.5){
  sel = x$coexpr >= quantile(x$coexpr, percentile)
  return(x[sel,])
}

scale.free.threshold = function(x, thresholds = seq(from=0.1, to=0.9, by=0.1), target.gamma = 2, verbose=TRUE){
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
