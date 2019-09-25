/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include "gfa_segment.h"

/**
 * Allocates memory for a new gfa_segment, creates new gfa_segment from parameters. Each segments represents a node in the 
 * GFA graph. Also allocates memory for the nucleotide sequence, which is deep-copied.
 * @param id - The id of the segment. Currently managed by the caller, should probably be automatic.
 * @param sequence - The nucleotide sequence that this segment represents.
 * @param orientation - The orientation of the sequence, as to be used in the GFA file.
 * @return pointer to the segment.
 */
gfa_segment* gfa_segment_new(int id, char* sequence, Orientation orientation)
{
    gfa_segment* new_segment = malloc(sizeof(gfa_segment));
    int length = strlen(sequence);
    new_segment->m_nucleotide_sequence = malloc(sizeof(char) * (length + 1));
    new_segment->m_nucleotide_sequence[length] = '\0';
    strcpy(new_segment->m_nucleotide_sequence, sequence);
    new_segment->m_segment_id = id;
    new_segment->m_orientation = orientation;
    return new_segment;
}

/**
 * Frees memory allocated by and to this segment.
 * @param segment
 */
void gfa_segment_destroy(gfa_segment* segment)
{
    if(segment != NULL)
    {
        free(segment->m_nucleotide_sequence);
        free(segment);
        segment = NULL;
    }
}

/**
 * Writes a segment (S) line for this segment to the GFA file.
 * @param segment
 * @param gfa_file
 */
void write_gfa_segment(const gfa_segment* const segment, gfa_file_wrapper* gfa_file)
{
    fprintf(gfa_file->m_file, "S\tp%llds%i\t%lu\t%s\n", 
                        gfa_file->m_path_id,
                        segment->m_segment_id, 
                        strlen(segment->m_nucleotide_sequence), 
                        segment->m_nucleotide_sequence);  
}

/**
 * Writes an edge (E) line to the GFA file. This represents an edge from the first to the second segment.
 * Note that we assume the end of first segment is joined to the beginning of the second segment (i.e. there is no overlap).
 * @param first_segment
 * @param second_segment
 * @param gfa_file
 */
void write_gfa_edge(const gfa_segment* const first_segment, const gfa_segment* const second_segment, gfa_file_wrapper* gfa_file)
{
    assert(first_segment != NULL);
    assert(second_segment != NULL);
    long unsigned first_length = strlen(first_segment->m_nucleotide_sequence);
    char first_orientation = first_segment->m_orientation == forward ? '+' : '-';
    char second_orientation = second_segment->m_orientation == forward ? '+' : '-';
    fprintf(gfa_file->m_file, "E\tp%llds%i_p%llds%i\tp%llds%i%c\tp%llds%i%c\t%lu$\t%lu$\t0\t0\t*\n",
                        gfa_file->m_path_id, first_segment->m_segment_id, gfa_file->m_path_id, second_segment->m_segment_id,
                        gfa_file->m_path_id, first_segment->m_segment_id, first_orientation,
                        gfa_file->m_path_id, second_segment->m_segment_id, second_orientation,
                        first_length, first_length);                       
}

/**
 * Allocates memory for and creates a new gfa_segment_array. This is a wrapper for an array of gfa_segments that will grow as necessary 
 * (by reallocating memory if there is not enough). Memory is allocated for an array of gfa_segments of an inital size. It is expected 
 * that almost all gfa_segment_arrays will only need to contain one gfa_segment
 * 
 * Note that the actual array is an array of pointers. gfa_segment_arrays take ownership of the objects pointed to in this array.
 * gfa_segments should never be appended to the array manually, but only through the append function. Failure to do this will
 * almost certainly result in a memory access error (if only we were using an object orientated language...)
 * @param inital_size
 * @return Pointer to gfa_segment_array.
 */
gfa_segment_array* gfa_segment_array_new(int inital_size)
{
    gfa_segment_array* new_segment = calloc(1, sizeof(gfa_segment_array));
    new_segment->m_max_size = inital_size;
    new_segment->m_array = calloc(new_segment->m_max_size, sizeof(gfa_segment*));
    new_segment->m_length = 0;
    return new_segment;
}

/**
 * Frees memory allocated to and by this gfa_segment_array. Note that this includes all the gfa_segments
 * that are pointed to in the array, and the array itself.
 * @param array
 */
void gfa_segment_array_destroy(gfa_segment_array* array)
{
    if(array != NULL)
    {
        for(int i = 0; i < array->m_max_size; i++)
        {
            if(i < array->m_length)
            {
                free(array->m_array[i]);
                array->m_array[i] = NULL;
            }
            else
            {
                assert(array->m_array[i] == NULL);
            }       
        }
        free(array->m_array);
        array->m_array = NULL;
        free(array);
        array = NULL;
    }
}

/**
 * Access the ith element of the array. Probably unnecessary but here out of principle.
 * @param array
 * @param i
 * @return the ith pointer in the array. Returns NULL if i is out of range.
 */
gfa_segment* gfa_segment_array_get(const gfa_segment_array* array, int i)
{
    if(i < array->m_length)
    {
        return array->m_array[i];
    }
    else
    {
        printf("[gfa_segment_array_get] Error: attempting to access out of bounds index. Returning NULL");
        return NULL;
    }
}

/**
 * Append a new gfa_segment to the array. The user passes the values needed to create a new gfa_segment, instead of creating a new
 * gfa_segment themselves. This way, the gfa_segment_array creates the gfa_segment, and manages the memory allocated to it. Since
 * it's expected that almost every array will only hold one or two segments, only one extra space is allocated if the limit is reached.
 * @param array
 * @param id
 * @param sequence
 * @param orientation
 */
void gfa_segment_array_append(gfa_segment_array* array, int id, char* sequence, Orientation orientation)
{
    assert(array->m_length < array->m_max_size);
    if(array->m_length == array->m_max_size)
    {
        if(!gfa_segment_array_allocate_extra_space(array, 1))
        {
            printf("[gfa_segment_array_append] Error: Could not allocate enough extra space!");
            return;            
        }
    }
    gfa_segment* segment = gfa_segment_new(id, sequence, orientation);
    array->m_array[array->m_length++] = segment;
}

/**
 * Copy src_array and append it to dst_array. This performs a deep copy of all the gfa_segments in src_array, so it is safe to 
 * destroy src_array after calling this function.
 * Since this could be painfully slow if there are lots of elements in src_array, (only one space is reallocated at a time
 * in append) we reallocate exactly the amount of space that dst_array will need before appending.
 * @param dst_array
 * @param src_array
 */
void gfa_segment_array_merge(gfa_segment_array* dst_array, const gfa_segment_array* src_array)
{
    if(dst_array == NULL || dst_array->m_array == NULL)
    {
        printf("[gfa_segment_array_merge] Error: dst array not properly initialized. Cannot merge.");
        return;       
    }    
    if(src_array == NULL || src_array->m_array == NULL)
    {
        printf("[gfa_segment_array_merge] Error: src array not properly initialized. Cannot merge.");
        return;         
    }
    int available_space = dst_array->m_max_size - dst_array->m_length;
    int required_space = src_array->m_length;
    if(required_space > available_space)
    {
        if(!gfa_segment_array_allocate_extra_space(dst_array, required_space - available_space))
        {
            printf("[gfa_segment_array_merge] Error: Could not allocate enough extra space!");
            return;
        }
    }
    for(int i = 0; i < src_array->m_length; i++)
    {
        gfa_segment* segment = gfa_segment_array_get(src_array, i);
        gfa_segment_array_append(dst_array, segment->m_segment_id, segment->m_nucleotide_sequence, segment->m_orientation);
    }
}

/**
 * Write a segment (S) line for each gfa_segment in the array to gfa_file.
 * @param array
 * @param gfa_file
 */
void write_gfa_segment_array(const gfa_segment_array* array, gfa_file_wrapper* gfa_file)
{
    for(int i = 0; i < array->m_length; i++)
    {
        gfa_segment* segment = gfa_segment_array_get(array, i);
        write_gfa_segment(segment, gfa_file);
    }
}

/**
 * Private function, only to be called internally. Reallocates the array of gfa_segment pointers
 * to a new address with extra_space more slots available.
 * @param array
 * @param extra_space
 * @return Whether reallocation was successful.
 */
boolean gfa_segment_array_allocate_extra_space(gfa_segment_array* array, int extra_space)
{
    if(extra_space <=0)
    {
        printf("[gfa_segment_array_allocate_extra_space] Error: Extra space to allocate must be positive.");
        return false;
    }
    int new_max = array->m_max_size + extra_space;
    gfa_segment** new_array = realloc(array->m_array, sizeof(gfa_segment*) * new_max);
    if(new_array == NULL)
    {
        printf("[gfa_segment_array_allocate_extra_space] Error: Could not reallocate array, returning.");
        return false;
    }
    array->m_max_size = new_max;
    array->m_array = new_array;
    for(int i = array->m_length; i < array->m_max_size; i++)
    {
        array->m_array[i] = NULL;
    }
    return true;
}