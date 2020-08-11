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
 ************************************************************************
 *
 * This file is modified from source that was part of CORTEX. The
 * original license notice for that is given below.
 *
 ************************************************************************
 *
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 *
 ************************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 ************************************************************************/

/************************************************************************
 * branches.h
 ************************************************************************/
#ifndef BRACHES_H
#define BRANCHES_H

void branches_get_path(char * filename, int max_length, boolean with_coverages, dBGraph * db_graph) ;

PathArray * branches_get_all_paths_from(pathStep * first, int max_length, dBGraph * db_graph);

void branches_get_all_paths_from_with_callback(pathStep * first, int max_length, void (*path_action) (Path * path), dBGraph * db_graph);

Path * branches_get_path_between(pathStep * first, pathStep * last, int max_length,int max_depth,  Path * path,  dBGraph * db_graph);

Path * branches_get_path_between_paths(Path * first_path, Path * last_path, int max_length, int max_depth, Path * path,  dBGraph * db_graph);

#endif // BRANCHES_H
