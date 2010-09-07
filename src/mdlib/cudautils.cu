#include <stdlib.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cuda.h"
#include "cudautils.h"

/*** General CUDA data operations ***/
/* TODO: create a cusmalloc module that implements similar things as smalloc */

int download_cudata(void * h_dest, void * d_src, size_t bytes)
{
    cudaError_t stat;
    
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    stat = cudaMemcpy(h_dest, d_src, bytes, cudaMemcpyDeviceToHost);
    CU_RET_ERR(stat, "DtoH cudaMemcpy failed");

    return 0;
}

int download_cudata_alloc(void ** h_dest, void * d_src, size_t bytes)
{ 
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    smalloc(*h_dest, bytes);

    return download_cudata(*h_dest, d_src, bytes);
}

int upload_cudata(void * d_dest, void * h_src, size_t bytes)
{   
    cudaError_t stat;

    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    stat = cudaMemcpy(d_dest, h_src, bytes, cudaMemcpyHostToDevice);
    CU_RET_ERR(stat, "HtoD cudaMemcpy failed");

    return 0;
}

int upload_cudata_alloc(void ** d_dest, void * h_src, size_t bytes)
{
    cudaError_t stat;
    
    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    stat = cudaMalloc(d_dest, bytes);
    CU_RET_ERR(stat, "cudaMalloc failed in upload_cudata_alloc");

    return upload_cudata(*d_dest, h_src, bytes);
}

/* pinned alloc */
void * pmalloc(size_t bytes, FILE *fplog)
{
    cudaError_t err;
    void *      ptr;

    if (bytes == 0)
	return NULL;

    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(fplog, "Just caught a previously occured CUDA error: %s, will try to continue.\n", 
                    cudaGetErrorString(err));
        fprintf(stderr, "Just caught a previously occured CUDA error: %s, will try to continue.\n", 
                    cudaGetErrorString(err));
    }

    cudaMallocHost(&ptr, bytes);
	
    if ((err = cudaGetLastError()) != cudaSuccess)
    {
    	gmx_fatal(FARGS, "palloc of size %d bytes failed: %s (%d)\n", 
                    bytes, cudaGetErrorString(err), err);
                    
	exit(1);
    }
    return ptr;
}

/* pinned free */
void pfree(void *h_ptr, FILE *fplog) 
{
    if (h_ptr)
	cudaFreeHost(h_ptr);
}

void calc_grid_block_conf(int *threds_per_block, int *nb_blocks, int dim,
			int min_blocks, int max_blocks,
			int min_threds, int default_threds, int max_threds)
{
    *threds_per_block = default_threds;
    *nb_blocks = (dim % default_threds == 0) ?
		    (dim / default_threds) :
		    (dim / default_threds) + 1;	

//	if (*nbr_ctas > max_ctas) 
//		*nbr_ctas = max_ctas;

    if (*nb_blocks > GRID_MAX_DIM) 
        *nb_blocks = GRID_MAX_DIM;

}
