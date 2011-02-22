#include <stdlib.h>

#include "gmx_fatal.h"

#include "cudautils.h"
#include "cupmalloc.h"

/* page-locked alloc */
void pmalloc(void **h_ptr, size_t nbytes)
{
    cudaError_t stat;
    char        strbuf[50]; // FIXME what's the gmx macro for default small char buffers?
    int         flag = cudaHostAllocDefault; // TODO put here flag selection

    if (nbytes <= 0)
    {
        *h_ptr = NULL;
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaMallocHost(h_ptr, nbytes, flag);    
    sprintf(strbuf, "cudaMallocHost of size %d bytes failed", (int)nbytes);
    CU_RET_ERR(stat, strbuf);  
}

void pmalloc_wc(void **h_ptr, size_t nbytes)
{
    cudaError_t stat;
    char        strbuf[50]; // FIXME what's the gmx macro for default small char buffers?
    int         flag = cudaHostAllocDefault || cudaHostAllocWriteCombined; // TODO put here flag selection

    if (nbytes <= 0)
    {
        *h_ptr = NULL;
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaMallocHost(h_ptr, nbytes, flag);    
    sprintf(strbuf, "cudaMallocHost of size %d bytes failed", (int)nbytes);
    CU_RET_ERR(stat, strbuf);  
}

/* page locked free */
void pfree(void *h_ptr) 
{
    cudaError_t stat; 

    if (h_ptr == NULL)
    {        
        return;
    }

    CU_CHECK_PREV_ERR();

    stat = cudaFreeHost(h_ptr);
    CU_RET_ERR(stat, "cudaFreeHost failed");
}
