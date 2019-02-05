#' Clean a string so that it's a valid file name
#' 
#' Removes all characters that are not: literals, numbers, underscore, dash, dot
#'
#' @param old_filename filename to be cleaned (without the full path!)
#'
#' @return cleaned filename
#' @keywords internal
#'
#' @examples
clean.filename = function(old_filename){
  newname = gsub('[^a-zA-Z0-9_\\.\\-]', '', old_filename)
  return(newname)
}

#' Export functional modules network
#' 
#' Functional modules for a single cluster, after MTGO run. You can either collapse the
#' nodes of the same functional module, or keep them separated.
#'
#' @param infolder the cluster folder
#' @param collapse.modules should all elements of a module be collapsed into a single node?
#'
#' @return a visNetwork network object
#' @export
export.network.modules = function(infolder, collapse.modules = TRUE){
  if (!dir.exists(infolder)){
    stop(paste('Passed folder does not exist:', infolder))
  }
  
  # DATA LOAD AND SETUP -----------------------------------------------------
  #where to save
  outfolder = file.path(infolder, 'Networks')
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  
  #nodes
  data.filename = file.path(infolder, 'Nodes_Best_QGO.txt')
  data=read.table(data.filename, header=TRUE, stringsAsFactors=FALSE, sep='\t')
  
  #edges
  edge.filename = file.path(infolder, 'edges.tsv')
  edge = read.csv(edge.filename, header=FALSE, stringsAsFactors=FALSE, sep='\t')
  colnames(edge) = c('from', 'to')
  
  #removing duplicated edges (a-b and b-a are synonyms)
  sel = edge$from > edge$to
  edge[sel, 'tmp'] = edge[sel, 'from']
  edge[sel, 'from'] = edge[sel, 'to']
  edge[sel, 'to'] = edge[sel, 'tmp']
  edge$tmp = NULL
  edge=unique(edge)
  
  #placeholders
  outfile.name = NULL
  net.title = NULL

  if (collapse.modules){
    # COLLAPSED MODULES -------------------------------------------------------
    net.title = 'Functional modules'
    
    #NODES
    nodes = data.frame(table(data$Nodes.Gene.Ontology), stringsAsFactors = FALSE)
    colnames(nodes) = c('id', 'value')
    
    #labels are shown to the user
    nodes$label = nodes$id
    
    #tooltip with the components of this functional group
    tmp = plyr::ddply(data, plyr::.(Nodes.Gene.Ontology), .fun = function(x){
      return(paste(x$Nodes, collapse = '<br>'))
    })
    nodes$title = plyr::mapvalues(nodes$id, from=tmp$Nodes.Gene.Ontology, to=tmp$V1, warn_missing = FALSE)
    
    #EDGES
    #so far edges connect proteins/genes, but we need to pass to
    #functional modules
    edge$from = plyr::mapvalues(x=edge$from, from = data$Nodes, to = data$Nodes.Gene.Ontology, warn_missing = FALSE)
    edge$to   = plyr::mapvalues(x=edge$to,   from = data$Nodes, to = data$Nodes.Gene.Ontology, warn_missing = FALSE)
    
    #ensuring symmetry
    sel = edge$from > edge$to
    edge[sel, 'tmp'] = edge[sel, 'from']
    edge[sel, 'from'] = edge[sel, 'to']
    edge[sel, 'to'] = edge[sel, 'tmp']
    edge$tmp = NULL
    
    #collapsing redundant edges, but keeping count
    edge = data.frame(dplyr::count(edge, from, to), stringsAsFactors = FALSE)
    colnames(edge) = c('from', 'to', 'value')
    
    #OUTFILE
    outfile.name = file.path(outfolder, "ClusterNetwork.html")
  }else{
    # FULL MODULES ------------------------------------------------------------
    net.title = 'Full network'
    
    #NODES
    nodes = data.frame(
      id = data$Nodes,
      group = data$Nodes.Gene.Ontology,
      label = data$Nodes,
      title = data$Nodes.Gene.Ontology,
      stringsAsFactors = FALSE)
    
    #EDGES
    #Edges are already correct
    
    #OUTFILE
    outfile.name = file.path(outfolder, "FullNetwork.html")
  }
  
  # SAVING NETWORK ---------------------------------------------------------
  network = visNetwork::visNetwork(nodes=nodes, edges=edge, height="800px",width="100%", main = net.title, submain = infolder)
  network = visNetwork::visExport(network)  
  network = visNetwork::visNodes(network, size = 10)
  #network = visNetwork::visPhysics(network, solver = "repulsion",repulsion = list(gravitationalConstant = 0,centralGravity= 0, springConstant= 0),stabilization = TRUE)
  #network = visNetwork::visPhysics(network, solver = "repulsion", stabilization = TRUE)
  #network = visNetwork::visIgraphLayout(network, layout = "layout_with_drl")
  network = visNetwork::visIgraphLayout(network, layout = "layout_in_circle")
  
  if (!collapse.modules){
    network = visNetwork::visOptions(network, selectedBy = "group")
    network = visNetwork::visLegend(network)
    network = visNetwork::visEdges(network, color = list(color = "steelblue", highlight = "red"))
  }
  
  visNetwork::visSave(network, file = outfile.name)
  
  #returning the created network
  return(network)
}

#' Export single module network(s)
#' 
#' Export single module network(s)
#'
#' @param infolder working dir after MTGO run
#' @param modules name of the module(s) to be exported, or NULL to export all of them
#'
#' @return a list of visNetwork networks, indexed by module name
#' @export
export.network.single.module = function(infolder, modules = NULL){
  if (!dir.exists(infolder)){
    stop(paste('Passed folder does not exist:', infolder))
  }
  
  # DATA LOAD AND SETUP -----------------------------------------------------
  #where to save
  outfolder = file.path(infolder, 'Networks')
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  
  #nodes
  data.filename = file.path(infolder, 'Nodes_Best_QGO.txt')
  data=read.table(data.filename, header=TRUE, stringsAsFactors=FALSE, sep='\t')
  
  #edges
  edge.filename = file.path(infolder, 'edges.tsv')
  edge = read.csv(edge.filename, header=FALSE, stringsAsFactors=FALSE, sep='\t')
  colnames(edge) = c('from', 'to')
  
  #removing duplicated edges (a-b and b-a are synonyms)
  sel = edge$from > edge$to
  edge[sel, 'tmp'] = edge[sel, 'from']
  edge[sel, 'from'] = edge[sel, 'to']
  edge[sel, 'to'] = edge[sel, 'tmp']
  edge$tmp = NULL
  edge=unique(edge)
  
  # CREATING NETWORK --------------------------------------------------------
  #do we have a single module or should we do all of them?
  if (is.null(modules)){
    modules = unique(data$Nodes.Gene.Ontology)  
  }
  
  #room for the final, returned, list of network
  result = list()
  
  #for each of the examined modules
  for(module.curr in modules){
    #subset to only nodes in the current module
    data.curr = subset(data, Nodes.Gene.Ontology == module.curr)
    
    #preparing nodes object for visNetwork
    nodes = data.frame(
      stringsAsFactors = FALSE,
      id = data.curr$Nodes,
      label = data.curr$Nodes
    )
    
    #subset edges to current vertexes only
    edge.curr = edge[(edge$from %in% nodes$id) & (edge$to %in% nodes$id),]
    
    
    #building the network
    network = visNetwork::visNetwork(nodes=nodes, edges=edge.curr, height="800px",width="100%", main = module.curr, submain = infolder)
    network = visNetwork::visExport(network)
    network = visNetwork::visPhysics(network, solver = "repulsion",repulsion = list(gravitationalConstant = 0,centralGravity= 0, springConstant= 0),stabilization = TRUE)
    network = visNetwork::visIgraphLayout(network, layout = "layout_with_drl")
    
    #saving
    outfile.name = file.path(outfolder, clean.filename(paste(sep='', 'module_', module.curr,'.html')))
    visNetwork::visSave(network, file = outfile.name)
    result[[module.curr]] = network
  }
  
  return(result)
}