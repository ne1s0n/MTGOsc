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

#' Dummy function for import management
#'
#' This function allows to have a single place where to place all the generic @import.
#' Usually only a couple are required, but for things used all over the places (like stats package)
#' there is no real single place where to put the import. It's cleaner in this way.
#'
#' @return nothing
#' @keywords internal
#' @import stats
#' @import utils
dummy.function = function(){}
