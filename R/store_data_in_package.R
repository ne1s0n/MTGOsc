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

# This is an utility script to load a Seurat object from (my) computer
# into the MTGOsc package. It's here for handyness, repeatibility, and to show
# the process to the curious investigator. If you want to run the code do a
# big uncomment below this point.

# setwd("~/research/MTGOsc/")
#
# #loading bladder full dataset
# bladder.big = readRDS('~/research/Seurat2MTGO.notes/dataset5.2/data/Bladder.dataset5.RDS')
#
# #going to-lower to simplify match
# bladder.big@var.genes = tolower(bladder.big@var.genes)
# rownames(bladder.big@raw.data) = tolower(rownames(bladder.big@raw.data))
# rownames(bladder.big@data) = tolower(rownames(bladder.big@data))
# rownames(bladder.big@hvg.info) = tolower(rownames(bladder.big@hvg.info))
# bladder.big@calc.params$ScaleData$genes.use = tolower(bladder.big@calc.params$ScaleData$genes.use)
#
# #finding differentially expressed genes for each cluster
# markers = FindAllMarkers(bladder.big, logfc.threshold = 1)
# table(markers$cluster)
#
# #listing the available clusters, with cell count
# table(bladder.big@ident)
# #subsetting to just three clusters: a big one, a medium one, a small one
#
# #Stromal cell_Dpt high(Bladder) : 651 cells
# #Basal epithelial cell(Bladder) : 327 cells
# #Umbrella cell(Bladder) : 50 cells
# target = c('Stromal cell_Dpt high(Bladder)', 'Basal epithelial cell(Bladder)', 'Umbrella cell(Bladder)')
# bladder = SubsetData(bladder.big, ident.use = target, do.clean = TRUE, subset.raw = TRUE)
#
# #and subsetting also to differentially expressed genes only
# markers = subset(markers, cluster %in% target)
# bladder@data = bladder@data[unique(markers$gene),]
#
# #removing the heavy parts of the Seurat object, not used by MTGOsc
# bladder@raw.data = NULL
# bladder@scale.data = NULL
#
# #storing in MTGOsc
# usethis::use_data(bladder, markers, overwrite = TRUE)
#
# #loading the pathway dictionary
# mouse.pathways = read.table('~/research/Seurat2MTGO.notes/dataset5.2/data/Ensembl2Reactome_All_Levels.mus.musculus_mgi_small_caps.txt', stringsAsFactors = FALSE, sep='\t', quote='')
# colnames(mouse.pathways) = c('gene', 'pathway')
#
# #storing in MTGOsc
# usethis::use_data(mouse.pathways, overwrite = TRUE)
#
# #loading the ground truth, for method comparison
# load(file = '~/research/MTGOsc.notes/script/danila_AffinityScore/ground_truth.save')
#
# #storing in MTGOsc
# usethis::use_data(gGT, overwrite = TRUE)

