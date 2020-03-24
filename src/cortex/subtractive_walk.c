/************************************************************************
 *
 * This file is part of MetaCortex
 *
 * Authors:
 *     Richard M. Leggett (richard.leggett@earlham.ac.uk) and
 *     Martin Ayling (martin.ayling@earlham.ac.uk) and
 *     Samuel Martin (samuel.martin@earlham.ac.uk)
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <sys/stat.h>
#include <libgen.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include "cleaning.h"
#include "coverage_walk.h"
#include "dB_graph.h"
#include "element.h"
#include "graph_stats.h"
#include "logger.h"
#include "metacortex.h"
#include "metagraphs.h"
#include "subtractive_walk.h"

#define COVERAGE_BINS 10
#define COVERAGE_BIN_SIZE 1
#define MAX_BRANCHES 10
#define GRAPH_LOG10_LIMIT 10 // little hacky to do this here, because it needs to match size of subgraph_dist in graph_stats.h
#define NUM_BEST_NODES 5


// ----------------------------------------------------------------------
// Work through graph, count coverage, X, Y nodes
// ----------------------------------------------------------------------
void subtractive_walk(dBGraph * graph, char* consensus_contigs_filename,
                         int min_subgraph_kmers, int min_contig_size, int max_node_edges, float delta_coverage,
                         int linked_list_max_size, int walk_paths)
{
    FILE* fp_analysis;
    FILE* fp_report;
    FILE* fp_degrees;
    FILE* fp_contigs_fasta;
    long int Contig_Branches[MAX_BRANCHES];
    char* seq = calloc(256, 1);
    long int total_nodes = 0;
    int i;  int j;
    int counter= 0;
    int min_distance = 0; //10 * (graph->kmer_size);  // NOTE: needs to be a cmd_line option

    char cwd[1024];

    if (getcwd(cwd, sizeof(cwd)) != NULL){
        // do NOTHING
    }
    else{
        log_and_screen_printf("CWD command returned NULL\n");
    }

    char*  graph_wd = calloc(256, 1);

    Path *simple_path = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_fwd = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);
    Path *path_rev = path_new(MAX_EXPLORE_PATH_LENGTH, graph->kmer_size);

    GraphInfo* nodes_in_graph = calloc(1,sizeof(GraphInfo));
    new_GraphInfo(nodes_in_graph);

    // array to bin coverage 0-5, 5-10, 10-15..95-100
    long int Coverage_Dist[COVERAGE_BINS*COVERAGE_BIN_SIZE]; // will this work?
    char analysis_filename[256];
    char degrees_filename[256];

    Queue* graph_queue;
    for(i=0;i<MAX_BRANCHES;i++){
        Contig_Branches[i]=0;
    }
    // Initialise Coverage_Dist  int i;
    for(i=0;i<(COVERAGE_BINS*COVERAGE_BIN_SIZE);i++){
        Coverage_Dist[i]=0;
    }


    /* Open contigs file */
    fp_contigs_fasta = fopen(consensus_contigs_filename, "w");
    if (!fp_contigs_fasta) {
        log_and_screen_printf("ERROR: Can't open contig file.\n%s\n", consensus_contigs_filename);
        exit(-1);
    }

    remove_file_extension(consensus_contigs_filename);

    /* Open the analysis file */
    sprintf(analysis_filename, "%s.analysis", consensus_contigs_filename);
    fp_analysis = fopen(analysis_filename, "w");
    if (!fp_analysis) {
        log_and_screen_printf("ERROR: Can't open analysis file.\n");
        exit(-1);
    }

    /* Open the sugraph degree file */
    sprintf(degrees_filename, "%s.degrees", consensus_contigs_filename);
    fp_degrees = fopen(degrees_filename, "w");
    if (!fp_degrees) {
        log_and_screen_printf("ERROR: Can't open degrees file.\n");
        exit(-1);
    }

    // check for graphs dir existance
    if (basename(consensus_contigs_filename)==consensus_contigs_filename){
        log_and_screen_printf("(Relative path for contig output given, prefixing CWD)\n");
        // returns '.' which breaks other paths later on
        sprintf(graph_wd, "%s/graphs/", cwd);
    }
    else{
        // dirname modifies 'consensus_contigs_filename' on some platforms, shifted in here to avoid that
        sprintf(graph_wd, "%s/graphs/", dirname(consensus_contigs_filename));
        sprintf(analysis_filename, "%s%s.tex", graph_wd, basename(consensus_contigs_filename));
    }

    mkdir(graph_wd, 777);

    log_and_screen_printf("graphs dir\t%s\n", graph_wd);
    log_and_screen_printf("graphs\t%s\n", analysis_filename);

    /* Open the DIGEST file */
    fp_report = fopen(analysis_filename, "w");
    if (!fp_report) {
        log_and_screen_printf("ERROR: Can't open analysis (DIGEST) file.\n\t%s\n", analysis_filename);
        exit(-1);
    }

    // header line for degrees file
    for(i=0;i<5;i++){
        for(j=0;j<5;j++){
            fprintf(fp_degrees,"for[%d]rev[%d]\t", i, j);
        }
    }
    fprintf(fp_degrees,"total\n");

    db_graph_reset_flags(graph);

    // called by stats_traversal(), walks paths from a branch
    void find_path_length_with_first_edge_all_colours(Nucleotide n, dBNode * node, int * path_length, Orientation orientation) {
        if (db_node_edge_exist_any_colour(node, n, orientation)) {
            pathStep first_step;
            Path * new_path;
            new_path = path_new(MAX_EXPLORE_NODES, graph->kmer_size);
            first_step.node = node;
            first_step.orientation = orientation;
            first_step.label = n;
            first_step.flags = 0;

            db_graph_get_perfect_path_with_first_edge_all_colours(&first_step, &db_node_action_do_nothing, new_path, graph);
            new_path->length=1;
            * path_length += new_path->length;
            path_destroy(new_path);
        }
    }

    // Hash table iterator to label nodes
    void stats_traversal(dBNode * node) {
        //if (!db_node_check_flag_visited(node)) {
        int this_coverage = element_get_coverage_all_colours(node) - 1;
        int edges_forward= db_node_edges_count_all_colours(node, forward);
        int edges_reverse = db_node_edges_count_all_colours(node, reverse);
        int all_edges = edges_forward + edges_reverse;
        int local_distance = 1;
        int orientation;
        if (this_coverage<0) {
            log_and_screen_printf("Error: Coverage is <1 in the graph?\n");
            exit(-1);
        }


        // GRAPH DENSITY ESTIMATES
        // hash_table_traverse if edges_forward+edges reverse>2
        //   if(db_node_edge_exist_any_colour)(node, n, orientation)
        //   perfect path each edge
        //   count length
        //   average length > kmer?

        Contig_Branches[all_edges-1]++;

        if ((all_edges>2) && (all_edges<=max_node_edges)){

            if (linked_list_max_size){
                add_item(node, this_coverage, linked_list_max_size);
            }

            //log_and_screen_printf("\nWalking branch node...\n");
            // Look at all paths out from here

            orientation = forward;
            int i;
            for (i = 0; i < 4; i++) {
                find_path_length_with_first_edge_all_colours(i, node, &local_distance, orientation);
            }

            orientation = reverse;
            // nucleotide_iterator(&walk_if_exists);
            for (i = 0; i < 4; i++) {
                find_path_length_with_first_edge_all_colours(i, node, &local_distance, orientation);
            }

            local_distance = local_distance / all_edges;
        }
        else{
            // could check
            local_distance = min_distance;
        }

        // PARTITIONING - REMOVE EXTREMELY BRANCHED NODES
        if ((all_edges<=max_node_edges) && (local_distance>=min_distance))
        {
            if(this_coverage>COVERAGE_BINS*COVERAGE_BIN_SIZE-1)
            {
                this_coverage = COVERAGE_BINS*COVERAGE_BIN_SIZE-1;
            }

            Coverage_Dist[this_coverage]++;
            total_nodes++;

            // Look for Y shape branch forward orientation
            // The nodes at the top of the Y should contain different colours
            if (edges_forward > 1
                && edges_reverse == 1) {
                db_node_action_set_flag(node, BRANCH_NODE_FORWARD);
            }
            // Look for Y shape branch reverse orientation
            if (edges_reverse > 1
                && edges_forward == 1) {
                db_node_action_set_flag(node, BRANCH_NODE_REVERSE);
            }
            // Look for X-shaped branch
            if (edges_reverse > 1
                && edges_forward > 1) {
                db_node_action_set_flag(node, X_NODE);
            }
        }
        else{
            log_and_screen_printf("\nPruning node:\tedges\t%i\tdistance\t%i\n", all_edges, local_distance);
            cleaning_prune_db_node(node, graph);
        }
    } // stats_traversal()

    graph_queue = node_queue_new(METACORTEX_QUEUE_SIZE);
    if (!graph_queue) {
        log_and_screen_printf("Couldn't get memory for graph queue.\n");
        exit(-1);
    }
    /* Initialise temporaray path array buffers */
    path_array_initialise_buffers(graph->kmer_size);

    // Hash table iterator to initially walk all
    void explore_graph_size(dBNode * node) {
        if(db_node_check_for_any_flag(node, PRUNED | VISITED) == false){

            initialise_GraphInfo(nodes_in_graph);
            // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
            log_printf("\nGrowing graph from node");
            graph_queue->number_of_items = 0;

            log_printf("\n");

            // now with a subgraph, walk the graph counting degrees by graph
            explore_subgraphs(node, graph, nodes_in_graph);

            if (nodes_in_graph->total_size ==1) {
                // ignore; pruned node
            }
            else if (nodes_in_graph->total_size) {
                // print out the size of the current subgraph
                log_printf("graph size\t%i\n",nodes_in_graph->total_size);
                fprintf(fp_analysis, "%i\t%i\t",nodes_in_graph->branch_nodes,nodes_in_graph->total_size);
                binary_kmer_to_seq(&nodes_in_graph->highest_cov_in_subgraph, graph->kmer_size, seq);
                fprintf(fp_analysis, "%s\n", seq);

                // update graph wide stats
                print_degree_stats(nodes_in_graph, fp_degrees);
                if (nodes_in_graph->total_size>nodes_in_graph->largest_subgraph) {
                    nodes_in_graph->largest_subgraph=nodes_in_graph->total_size;
                }
                nodes_in_graph->branch_nodes_total=nodes_in_graph->branch_nodes_total+nodes_in_graph->branch_nodes;
                nodes_in_graph->num_subgraphs++;
                i=log10(nodes_in_graph->total_size);
                if(i>=GRAPH_LOG10_LIMIT){
                    i=GRAPH_LOG10_LIMIT-1;
                }
                nodes_in_graph->subgraph_dist[i]++;
                if(nodes_in_graph->total_size>min_subgraph_kmers){
                    nodes_in_graph->num_subgraphs_2k++;
                }
            }
            if (nodes_in_graph->branch_nodes>(MAX_BRANCHES-1)){
                nodes_in_graph->branch_nodes=MAX_BRANCHES-1;
            }
        }
    } // end of &explore_graph_size

    // Hash table iterator to walk graphs, produce paths
    void traversal_for_contigs(dBNode * node) {
        if(db_node_check_for_any_flag(node, PRUNED | VISITED) == false){

            dBNode* seed_node;
            initialise_GraphInfo(nodes_in_graph);

            // Grow graph from this node, returning the 'best' (highest coverage) node to store as seed point
            log_printf("\nGrowing graph from node");
            graph_queue->number_of_items = 0;

            log_printf("\n");

            // now with a subgraph, walk the graph counting degrees by graph
            // - this sets VISITED flag as true for many of the nodes in the graph.
            grow_graph_from_node_stats(node, &seed_node, graph, graph_queue, nodes_in_graph, delta_coverage);

            if (nodes_in_graph->total_size ==1) {
                // ignore; pruned node
                cleaning_prune_db_node(node, graph);
                db_node_action_set_flag(node, VISITED);
                log_printf("\t[Singleton pruned.]\n");
            }
            else if (seed_node == NULL) {
                printf("ERROR: Seed node is NULL, nodes in graph is %d\n", nodes_in_graph->total_size);
            } else if (nodes_in_graph->total_size) {
                int kmer_size = graph->kmer_size;
                char kmer_string[kmer_size + 1];
                binary_kmer_to_seq(&seed_node->kmer, kmer_size, kmer_string);
                log_printf("Seed node: %s\n", kmer_string);
                /* enough nodes to bother with? If so, get consensus contig */
                if (walk_paths && (nodes_in_graph->total_size >= min_subgraph_kmers)) {

                    // should be a perfect path? might be two paths though, if we started in the middle
                    // NOTE: unecessary coverage element but repeating the whole path finding without coverage
                    //  is more work than necessary I think. See what processing time it changes?
                    coverage_walk_get_path(seed_node, forward, NULL, graph, path_fwd);
                    coverage_walk_get_path(seed_node, reverse, NULL, graph, path_rev);

                    path_reverse(path_fwd, simple_path);
                    path_append(simple_path, path_rev);

                    simple_path->id = counter;
                    if (simple_path->length > (min_contig_size - graph->kmer_size)) {
                        log_printf("Write path of size %d\n", simple_path->length);
                        log_printf("graph size\t%i\n",nodes_in_graph->total_size);

                        // could save the path walking again here if needed, hold these figures in path structure
                        double average_coverage=0;
                        double standard_deviation = 0;
                        int min_coverage=0;
                        int max_coverage=0;
                        path_get_statistics(&average_coverage, &min_coverage, &max_coverage, simple_path);
                        path_get_coverage_standard_deviation(&standard_deviation, average_coverage, simple_path);
                        // NOTE: decision - minimum cov or average cov dictates confidence threshold met?
                        // Output for alternative formats
                        path_to_fasta(simple_path, fp_contigs_fasta);
                        
                        if(min_coverage < average_coverage - 2 * standard_deviation)
                        {
                            min_coverage = average_coverage - 2 * standard_deviation;
                        }
                        
                        for(i = 0; i < simple_path->length; i++) 
                        {
                            //TODO: make COLOUR-safe
                            simple_path->nodes[i]->coverage[0] -= min_coverage;
                            if(simple_path->nodes[i]->coverage[0] <= 0)
                            {
                                cleaning_prune_db_node(simple_path->nodes[i], graph);
                                db_node_action_set_flag(simple_path->nodes[i], VISITED);
                            }
                        }                            
                        
                        counter++;
                    } else {
                        log_printf("Didn't write path of size %d\n", simple_path->length);
                    }

            
                    /*	HERE - INSTEAD OF 'VISITING' AN X-NODE, REMOVE EDGES FROM BEFORE
                     AND AFTER IT ON CURRENT PATH */
                    
                    //unset VISITED flags from grow_graph_from_node_stats
                    dBNode* queue_node;
                    while (graph_queue->number_of_items > 0) {
                        queue_node = (dBNode*)queue_pop(graph_queue);
                        db_node_action_unset_flag(queue_node, VISITED);
                    }

                    /* Reset paths */
                    path_reset(simple_path);
                    path_reset(path_fwd);
                    path_reset(path_rev);

                } else  {
                    log_printf("  Number of nodes (%i) too small. Not outputting contig.\n", nodes_in_graph->total_size);
                } // end size graph size check


            } else {
                // catch graph size of zero? Not sure why this happens - grow-graph must be failing
                log_printf("graph size of zero?\n");
            }
            if (nodes_in_graph->branch_nodes>(MAX_BRANCHES-1)){
                nodes_in_graph->branch_nodes=MAX_BRANCHES-1;
            }
        }
    } // traversal_for_contigs

    // check each node in the graph, FLAG X&Y nodes (mark all nodes as visited)
    log_and_screen_printf("Stats traversal started...");
    hash_table_traverse(&stats_traversal, graph);
    log_and_screen_printf("DONE\n");

    // first line for stats output file
    fprintf(fp_analysis, "\n#Subgraph sizes\n");
    log_and_screen_printf("Graph size traversal started...");
    // first travesal - build subgraphs out, produce stats
    hash_table_traverse(&explore_graph_size, graph);
    log_and_screen_printf("DONE\n");

    // clear linked list
    log_and_screen_printf("Unique kmers before clearing:\t %lld\n", graph->unique_kmers);
    log_and_screen_printf("linked_list_max_size:\t %i\n", linked_list_max_size);
    if (linked_list_max_size){
        clear_list(graph);
    }
    log_and_screen_printf("Unique kmers after clearing:\t %lld\n", graph->unique_kmers);

    db_graph_reset_flags(graph);
    // second travesal - build subgraphs out, produce contigs
    log_and_screen_printf("Full traversal started...");
    hash_table_traverse(&traversal_for_contigs, graph);
    log_and_screen_printf("DONE\n");
       
    fclose(fp_analysis);
    fclose(fp_degrees);   
    db_graph_reset_flags(graph);
}