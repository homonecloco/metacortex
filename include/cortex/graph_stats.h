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
 * graph_stats.h
 ************************************************************************/
#ifndef GRAPH_STATS_H
#define GRAPH_STATS_H

#include <binary_kmer.h>
#include <node_queue.h>

typedef struct {
    int total_size;
    int largest_subgraph;
    int subgraph_dist[10];
    int num_subgraphs;
    int num_subgraphs_2k;
    int branch_nodes;
    int branch_nodes_total;
    int end_nodes;
    int node_degree[5][5];
    int best_coverage[5];
    int simple_bubbles;
    BinaryKmer kmer[5];
    BinaryKmer temp_kmer;
    BinaryKmer current_kmer;
    BinaryKmer highest_cov_in_subgraph;
    int highest_cov;
} GraphInfo;

int grow_graph_from_node_stats( dBNode* start_node, dBNode** best_node, 
                                dBGraph* graph, Queue* graph_queue, 
                                GraphInfo* nodes_in_graph, float delta);

void find_subgraph_stats(   dBGraph* graph, char* consensus_contigs_filename, 
                            int min_subgraph_kmers, int min_contig_size, 
                            int max_node_edges, float delta_coverage, 
                            int linked_list_max_size, int walk_paths, 
                            boolean gfa_fastg_output);

void print_degree_stats(GraphInfo * info, FILE* fp_degrees);

void initialise_GraphInfo(GraphInfo * info);

void new_GraphInfo(GraphInfo * info);

int explore_subgraphs(dBNode* start_node, dBGraph* graph, GraphInfo* nodes_in_graph);

#endif /* GRAPH_STATS_H */
