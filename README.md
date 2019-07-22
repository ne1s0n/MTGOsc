<!-- README.md is generated from README.Rmd. Please edit that file -->
MTGO-SC
=======

Extracting gene modules from single cell RNA-seq cell clusters
--------------------------------------------------------------

MTGO-SC is an adaptation for single cell RNA-sequencing (scRNA-seq) of [MTGO](https://gitlab.com/d1vella/MTGO), a biological network module detection algorithm. MTGO-SC integrates external gene annotations, such as the Gene Ontology terms or Reactome pathways with the gene expression networks obtained from single-cell DGE matrices.

A post-processing step after cell clustering
--------------------------------------------

The typical scRNA-seq pipeline is designed to group cells into meaningful clusters, representing cells of similar (sub)type or stare. MTGO-SC provides the opportunity for a further step in the analysis by extracting the gene interaction network from each cluster, and detecting the gene functional modules. Each module is labeled with an annotation from the source provided by the user (for example, Reactome pathways). MTGO-SC is designed to be integrated with [Seurat](https://github.com/satijalab/seurat), a toolkit for single cell data analysis.

Example of MTGO-SC application
------------------------------

A practical example of application on [MCA](http://bis.zju.edu.cn/MCA/) cluster, along with enrichment and term literature search, is provided in the [Vignette](https://github.com/ne1s0n/MTGOsc/blob/master/doc/MTGOsc.vignette.pdf).
