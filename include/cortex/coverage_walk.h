/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk)
 *
 * MetaCortex is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MetaCortex is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MetaCortex.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

/************************************************************************
 * coverage_walk.h
 *
 * Highest coverage walk used by MetaCortex.
 ************************************************************************/
#ifndef COVERAGE_WALK_H
#define COVERAGE_WALK_H

int coverage_walk_get_path(dBNode * node, Orientation orientation, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path, boolean check_bubbles);

void coverage_walk_get_path_forwards_and_backwards(dBNode * node, void (*node_action) (dBNode * node), dBGraph * db_graph, Path * path, boolean check_bubbles, int explore_length);

#endif //COVERAGE_WALK_H