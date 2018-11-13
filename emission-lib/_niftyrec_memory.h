/*
 *  _niftyrec_memory.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */

#ifndef _NIFTYREC_MEMORY_H
#define _NIFTYREC_MEMORY_H



#include <stdlib.h>
#include <stdio.h>
#include "nifti1_io.h"

#ifdef _USE_CUDA
    #include "_reg_cudaCommon.h"
#endif

#define ALLOCTYPE_HOST       0
#define ALLOCTYPE_CUDA       1
#define ALLOCTYPE_CUDA_ARRAY 2
#define ALLOCTYPE_NIFTI      3
#define RECORD_MAXELEMENTS   4096

typedef struct
{
    void *memory_section; 
    int platform; 
} alloc_record_element; 

typedef struct
{
    alloc_record_element *alloc_record_elements;
    int n_elements;
    int max_elements;
} alloc_record; 


alloc_record *alloc_record_create(int max_elements);
int alloc_record_destroy(alloc_record *record);
int free_record_element(alloc_record_element *element);
int alloc_record_add(alloc_record *record, void* memory_section, int platform);
int alloc_record_free_all(alloc_record *record);


#endif



