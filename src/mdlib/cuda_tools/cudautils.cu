#include <stdlib.h>

#include "gmx_fatal.h"
#include "smalloc.h"

#include "cuda.h"
#include "cudautils.h"

/*** General CUDA data operations ***/
/* TODO: create a cusmalloc module that implements similar things as smalloc */

static int cu_copy_D2H_generic(void * h_dest, void * d_src, size_t bytes, 
                               gmx_bool async = FALSE, cudaStream_t stream = 0)
{
    cudaError_t stat;
    
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    if (async)
    {
        stat = cudaMemcpyAsync(h_dest, d_src, bytes, cudaMemcpyDeviceToHost, stream);
        CU_RET_ERR(stat, "DtoH cudaMemcpyAsync failed");

    }
    else
    {
        stat = cudaMemcpy(h_dest, d_src, bytes, cudaMemcpyDeviceToHost);
        CU_RET_ERR(stat, "DtoH cudaMemcpy failed");
    }

    return 0;
}

int cu_copy_D2H(void * h_dest, void * d_src, size_t bytes)
{
    return cu_copy_D2H_generic(h_dest, d_src, bytes, FALSE);
}

int cu_copy_D2H_async(void * h_dest, void * d_src, size_t bytes, cudaStream_t stream = 0)
{
    return cu_copy_D2H_generic(h_dest, d_src, bytes, TRUE, stream);
}

int cu_copy_D2H_alloc(void ** h_dest, void * d_src, size_t bytes)
{ 
    if (h_dest == 0 || d_src == 0 || bytes <= 0)
        return -1;

    smalloc(*h_dest, bytes);

    return cu_copy_D2H(*h_dest, d_src, bytes);
}


static int cu_copy_H2D_generic(void * d_dest, void * h_src, size_t bytes, 
                               gmx_bool async = FALSE, cudaStream_t stream = 0)
{
    cudaError_t stat;

    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    if (async)
    {
        stat = cudaMemcpyAsync(d_dest, h_src, bytes, cudaMemcpyHostToDevice, stream);
        CU_RET_ERR(stat, "HtoD cudaMemcpyAsync failed");
    }
    else
    {
        stat = cudaMemcpy(d_dest, h_src, bytes, cudaMemcpyHostToDevice);
        CU_RET_ERR(stat, "HtoD cudaMemcpy failed");
    }

    return 0;
}

int cu_copy_H2D(void * d_dest, void * h_src, size_t bytes)
{   
    return cu_copy_H2D_generic(d_dest, h_src, bytes, FALSE);
}

int cu_copy_H2D_async(void * d_dest, void * h_src, size_t bytes, cudaStream_t stream = 0)
{   
    return cu_copy_H2D_generic(d_dest, h_src, bytes, TRUE, stream);
}

int cu_copy_H2D_alloc(void ** d_dest, void * h_src, size_t bytes)
{
    cudaError_t stat;

    if (d_dest == 0 || h_src == 0 || bytes <= 0)
        return -1;

    stat = cudaMalloc(d_dest, bytes);
    CU_RET_ERR(stat, "cudaMalloc failed in cu_copy_H2D_alloc");

    return cu_copy_H2D(*d_dest, h_src, bytes);
}

int cu_wait_event(cudaEvent_t stop, cudaEvent_t start, float *time)
{
    cudaError_t s;

    s = cudaEventSynchronize(stop);
    CU_RET_ERR(s, "cudaEventSynchronize failed in cu_wait_event");

    s = cudaEventElapsedTime(time, start, stop);
    CU_RET_ERR(s, "cudaEventElapsedTime failed in cu_wait_event");

    return 0;
}

/* Binds texture with name tex_name to the GPU global memory (of size elements) 
   pointed by d_ptr.
   Returns the offset that needs to be used when fetching from the texture.
 */
template <typename T>
size_t cu_bind_texture(const char *tex_name, const T *d_ptr, int size)
{
    cudaError_t             stat;
    cudaChannelFormatDesc   cd;
    const textureReference  *tex;
    char                    str[100];

    size_t offset;

    stat = cudaGetTextureReference(&tex, tex_name);
    sprintf(str, "cudaGetTextureReference on %s failed", tex_name);
    CU_RET_ERR(stat, str);
    cd = cudaCreateChannelDesc<T>();

    stat = cudaBindTexture(&offset, tex, d_ptr, &cd, size*sizeof(*d_ptr));
    sprintf(str, "cudaBindTexture on %s failed ", tex_name);
    CU_RET_ERR(stat, str);

    return offset;
}

/* Instantiate cu_bind_texture with float */
template size_t cu_bind_texture<float>(const char *, const float *, int);

/*! Unbinds texture with name tex_name. */
void cu_unbind_texture(const char *tex_name)
{
    cudaError_t             stat;
    const textureReference  *tex;
    char                    str[100];

    stat = cudaGetTextureReference(&tex, tex_name);
    sprintf(str, "cudaGetTextureReference on %s failed", tex_name);
    CU_RET_ERR(stat, str);
    stat = cudaUnbindTexture(tex);
    sprintf(str, "cudaUnbindTexture on %s failed ", tex_name);
    CU_RET_ERR(stat, str);
}

/*! Caculates and returns the time difference between event start and stop. */
float cu_event_elapsed(cudaEvent_t start, cudaEvent_t stop)
{
    float t = 0.0;
    cudaError_t stat;

    stat = cudaEventElapsedTime(&t, start, stop);
    CU_RET_ERR(stat, "cudaEventElapsedTime failed in cu_event_elapsed");

    return t;
}


/**** Operation on buffered arrays (arrays with "over-allocation" in gmx wording) */
/*! Frees the device memory pointed by d_ptr and resets the associated 
 *  size and allocation size variables to -1.
 */
void cu_free_buffered(void *d_ptr, int *n, int *nalloc)
{
    cudaError_t stat;

    if (d_ptr)
    {
        stat = cudaFree(d_ptr);
        CU_RET_ERR(stat, "cudaFree failed");
    }

    if (n)
    {
        *n = -1;
    }

    if (nalloc)
    {
        *nalloc = -1;
    }
}

/*! Reallocates the device memory pointed by d_ptr and copies the data from the 
 * location pointed by h_src host-side pointer. Allocation is buffered and 
 * therefor freeing is only needed if the previously allocated space is not 
 * enough. 
 */
void cu_realloc_buffered(void **d_dest, void *h_src, size_t type_size,
                                    int *curr_size, int *curr_alloc_size,
                                    int req_size,
                                    cudaStream_t stream,
                                    gmx_bool doAsync)
{
    cudaError_t stat;

    if (d_dest == NULL || req_size < 0)
    {
        return;
    }

    /* reallocate only if the data does not fit = allocation size is smaller 
       than the current requested size */
    if (req_size > *curr_alloc_size)
    {
        /* only free if the array has already been initialized */
        if (*curr_alloc_size >= 0)
        {
            cu_free_buffered(*d_dest, curr_size, curr_alloc_size);
        }

        *curr_alloc_size = 1.2 * req_size + 100;

        stat = cudaMalloc(d_dest, *curr_alloc_size * type_size);
        CU_RET_ERR(stat, "cudaMalloc failed in cu_free_buffered");
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if (doAsync)
        {
            cu_copy_H2D_async(*d_dest, h_src, *curr_size * type_size, stream);
        }
        else
        {
            cu_copy_H2D(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}
