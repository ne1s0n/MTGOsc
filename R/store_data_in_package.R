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
# #finding differentially expressed genes for each cluster
# markers = FindAllMarkers(bladder, logfc.threshold = 1)
# table(markers$cluster)
#
# #listing the available clusters, with cell count
# table(bladder@ident)
#
# #subsetting to just three clusters: a big one, a medium one, a small one
#
# #Stromal cell_Dpt high(Bladder) : 651 cells
# #Urothelium(Bladder) : 277 cells
# #Umbrella cell(Bladder) : 50 cells
# target = c('Stromal cell_Dpt high(Bladder)', 'Urothelium(Bladder)', 'Umbrella cell(Bladder)')
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
# devtools::use_data(bladder, markers, overwrite = TRUE)
