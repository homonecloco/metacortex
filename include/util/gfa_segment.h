/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   gfa_segment_array.h
 * Author: martins
 *
 * Created on 24 April 2019, 11:57
 */

#ifndef GFA_SEGMENT_H
#define GFA_SEGMENT_H

#include "global.h"

typedef struct 
{
    int m_segment_id;
    char* m_nucleotide_sequence;
    Orientation m_orientation;
} gfa_segment;

typedef struct
{
    FILE* m_file;
    int m_segment_count;
    long long m_path_id;
} gfa_file_wrapper;

typedef struct
{
    int m_length;
    int m_max_size;
    gfa_segment** m_array;
} gfa_segment_array;

gfa_segment* gfa_segment_new(int id, char* sequence, Orientation orientation);

void gfa_segment_destroy(gfa_segment* segment);

void write_gfa_segment(const gfa_segment* const segment, gfa_file_wrapper* gfa_file);

void write_gfa_edge(const gfa_segment* const first_segment, const gfa_segment* const second_segment, gfa_file_wrapper* gfa_file);

gfa_segment_array* gfa_segment_array_new(int inital_size);

void gfa_segment_array_append(gfa_segment_array* array, int id, char* sequence, Orientation orientation);

void gfa_segment_array_merge(gfa_segment_array* dst_array, const gfa_segment_array* src_array);

gfa_segment* gfa_segment_array_get(const gfa_segment_array* array, int i);

void gfa_segment_array_destroy(gfa_segment_array* array);

void write_gfa_segment_array(const gfa_segment_array* array, gfa_file_wrapper* gfa_file);

boolean gfa_segment_array_allocate_extra_space(gfa_segment_array* array, int extra_space);

#endif /* GFA_SEGMENT_H */

