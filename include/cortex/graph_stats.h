/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team:
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
 * **********************************************************************
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
 * **********************************************************************
 */

//#define MAX_SEEDS 200000000
//#define METACORTEX_QUEUE_SIZE 10000000
//#define MAX_EXPLORE_PATH_LENGTH 200000
//#define MAX_EXPLORE_NODES 200

//typedef struct {
//     dBNode* seed_node;
//     int graph_size;
//} SubGraphInfo;

typedef struct {
    int total_size;
    int branch_nodes;
    int end_nodes;
    int Y_degree_rev[4];
    int Y_degree_for[4];
    int X_degrees[8];
} GraphInfo;

void find_subgraph_stats(dBGraph* graph, char* consensus_contigs_filename);

void log_and_screen_print_stats(GraphInfo * nodes_in_graph);
