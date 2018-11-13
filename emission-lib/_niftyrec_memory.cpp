/*
 *  _niftyrec_memory.cpp
 *  
 *  NiftyRec
 *  Stefano Pedemonte, May 2012.
 *  CMIC - Centre for Medical Image Computing 
 *  UCL - University College London. 
 *  Released under BSD licence, see LICENSE.txt 
 */


#include "_niftyrec_memory.h"

alloc_record *alloc_record_create(int max_elements)
{
    alloc_record *record = (alloc_record *) malloc(sizeof(alloc_record)); 
    alloc_record_element *elements = (alloc_record_element *) malloc(max_elements*sizeof(alloc_record_element)); 
    record->alloc_record_elements = elements; 
    record->n_elements = 0; 
    record->max_elements = max_elements;
    for(int i=0;i<max_elements;i++)
    {
        record->alloc_record_elements[i].memory_section=NULL;
        record->alloc_record_elements[i].platform=0;
    }
    return record; 
}


int alloc_record_destroy(alloc_record *record)
{
    int status;
    //free all in case alloc_record_free_all() was not called. 
    status = alloc_record_free_all(record); 
    free(record->alloc_record_elements); 
    free(record); 
    return status; 
}


int free_record_element(alloc_record_element *element)
{
    if (element->platform == ALLOCTYPE_HOST)
        {
        //fprintf(stderr,"_niftyrec_memory: freeing GUEST memory %d\n",(int)element->memory_section);
        free(element->memory_section);
        element->memory_section=NULL;
        }
#ifdef _USE_CUDA
    else if (element->platform == ALLOCTYPE_CUDA)
        {
        //fprintf(stderr,"_niftyrec_memory: freeing CUDA memory %d\n",(int)element->memory_section);
        cudaFree(element->memory_section);
        element->memory_section=NULL;
        }
    else if (element->platform == ALLOCTYPE_CUDA_ARRAY)
        {
        //fprintf(stderr,"_niftyrec_memory: freeing CUDA ARRAY memory %d\n",(int)element->memory_section);
        cudaFreeArray((cudaArray*)element->memory_section);
        element->memory_section=NULL;
        }
#endif
    else if (element->platform == ALLOCTYPE_NIFTI)
        {
        //fprintf(stderr,"_niftyrec_memory: freeing NIFTI IMAGE memory %d\n",(int)element->memory_section);
        nifti_image_free((nifti_image*)element->memory_section);
        element->memory_section=NULL;
        }
    else return 1;
    return 0; 
}


int alloc_record_add(alloc_record *record, void* memory_section, int platform)
{
    if (record->n_elements>=record->max_elements)
        return 1;
    //fprintf(stderr,"_niftyrec_memory: adding record, type %d (%d records): %d \n",platform,record->n_elements+1,memory_section);
    record->alloc_record_elements[record->n_elements].memory_section = memory_section; 
    record->alloc_record_elements[record->n_elements].platform = platform; 
    record->n_elements+=1;  
    return 0; 
}

int alloc_record_free_all(alloc_record *record)
{
    int status = 0;
    for (int i=0;i<record->n_elements;i++)
        if (free_record_element(&(record->alloc_record_elements[i]))) {
        fprintf(stderr,"Unable to free memory at %lld.\n",(long long)record->alloc_record_elements[i].memory_section);
        status=1; 
        }
    record->n_elements=0;
    return status; 
}


