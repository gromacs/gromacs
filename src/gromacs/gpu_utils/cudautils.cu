/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#include "gmxpre.h"

#include "cudautils.cuh"

#include <cassert>
#include <cstdlib>

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/utility/smalloc.h"

/*** Generic CUDA data operation wrappers ***/

/*! Launches synchronous or asynchronous host to device memory copy.
 *
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
static int cu_copy_D2H_generic(void * h_dest, void * d_src, size_t bytes,
                               bool bAsync = false, cudaStream_t s = 0)
{
    cudaError_t stat;

    if (h_dest == NULL || d_src == NULL || bytes == 0)
    {
        return -1;
    }

    if (bAsync)
    {
        stat = cudaMemcpyAsync(h_dest, d_src, bytes, cudaMemcpyDeviceToHost, s);
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
    return cu_copy_D2H_generic(h_dest, d_src, bytes, false);
}

/*!
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_D2H_async(void * h_dest, void * d_src, size_t bytes, cudaStream_t s = 0)
{
    return cu_copy_D2H_generic(h_dest, d_src, bytes, true, s);
}

/*! Launches synchronous or asynchronous device to host memory copy.
 *
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
static int cu_copy_H2D_generic(void * d_dest, void * h_src, size_t bytes,
                               bool bAsync = false, cudaStream_t s = 0)
{
    cudaError_t stat;

    if (d_dest == NULL || h_src == NULL || bytes == 0)
    {
        return -1;
    }

    if (bAsync)
    {
        stat = cudaMemcpyAsync(d_dest, h_src, bytes, cudaMemcpyHostToDevice, s);
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
    return cu_copy_H2D_generic(d_dest, h_src, bytes, false);
}

/*!
 *  The copy is launched in stream s or if not specified, in stream 0.
 */
int cu_copy_H2D_async(void * d_dest, void * h_src, size_t bytes, cudaStream_t s = 0)
{
    return cu_copy_H2D_generic(d_dest, h_src, bytes, true, s);
}

/**** Operation on buffered arrays (arrays with "over-allocation" in gmx wording) *****/

/*!
 * If the pointers to the size variables are NULL no resetting happens.
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

/*!
 *  Reallocation of the memory pointed by d_ptr and copying of the data from
 *  the location pointed by h_src host-side pointer is done. Allocation is
 *  buffered and therefore freeing is only needed if the previously allocated
 *  space is not enough.
 *  The H2D copy is launched in stream s and can be done synchronously or
 *  asynchronously (the default is the latter).
 */
void cu_realloc_buffered(void **d_dest, void *h_src,
                         size_t type_size,
                         int *curr_size, int *curr_alloc_size,
                         int req_size,
                         cudaStream_t s,
                         bool bAsync = true)
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

        *curr_alloc_size = over_alloc_large(req_size);

        stat = cudaMalloc(d_dest, *curr_alloc_size * type_size);
        CU_RET_ERR(stat, "cudaMalloc failed in cu_free_buffered");
    }

    /* size could have changed without actual reallocation */
    *curr_size = req_size;

    /* upload to device */
    if (h_src)
    {
        if (bAsync)
        {
            cu_copy_H2D_async(*d_dest, h_src, *curr_size * type_size, s);
        }
        else
        {
            cu_copy_H2D(*d_dest, h_src,  *curr_size * type_size);
        }
    }
}

/*! \brief Return whether texture objects are used on this device.
 *
 * \param[in]   pointer to the GPU device info structure to inspect for texture objects support
 * \return      true if texture objects are used on this device
 */
static inline bool use_texobj(const gmx_device_info_t *dev_info)
{
    assert(!c_disableCudaTextures);
    /* Only device CC >= 3.0 (Kepler and later) support texture objects */
    return (dev_info->prop.major >= 3);
}

/*! \brief Set up texture object for an array of type T.
 *
 * Set up texture object for an array of type T and bind it to the device memory
 * \p d_ptr points to.
 *
 * \tparam[in] T        Raw data type
 * \param[out] texObj   texture object to initialize
 * \param[in]  d_ptr    pointer to device global memory to bind \p texObj to
 * \param[in]  sizeInBytes  size of memory area to bind \p texObj to
 */
template <typename T>
static void setup1DTexture(cudaTextureObject_t &texObj,
                           void                *d_ptr,
                           size_t               sizeInBytes)
{
    assert(!c_disableCudaTextures);

    cudaError_t      stat;
    cudaResourceDesc rd;
    cudaTextureDesc  td;

    memset(&rd, 0, sizeof(rd));
    rd.resType                = cudaResourceTypeLinear;
    rd.res.linear.devPtr      = d_ptr;
    rd.res.linear.desc        = cudaCreateChannelDesc<T>();
    rd.res.linear.sizeInBytes = sizeInBytes;

    memset(&td, 0, sizeof(td));
    td.readMode                 = cudaReadModeElementType;
    stat = cudaCreateTextureObject(&texObj, &rd, &td, NULL);
    CU_RET_ERR(stat, "cudaCreateTextureObject failed");
}

/*! \brief Set up texture reference for an array of type T.
 *
 * Set up texture object for an array of type T and bind it to the device memory
 * \p d_ptr points to.
 *
 * \tparam[in] T        Raw data type
 * \param[out] texObj   texture reference to initialize
 * \param[in]  d_ptr    pointer to device global memory to bind \p texObj to
 * \param[in]  sizeInBytes  size of memory area to bind \p texObj to
 */
template <typename T>
static void setup1DTexture(const struct texture<T, 1, cudaReadModeElementType> *texRef,
                           const void                                          *d_ptr,
                           size_t                                              sizeInBytes)
{
    assert(!c_disableCudaTextures);

    cudaError_t           stat;
    cudaChannelFormatDesc cd;

    cd   = cudaCreateChannelDesc<T>();
    stat = cudaBindTexture(nullptr, texRef, d_ptr, &cd, sizeInBytes);
    CU_RET_ERR(stat, "cudaBindTexture failed");
}

template <typename T>
void initParamLookupTable(T                        * &d_ptr,
                          cudaTextureObject_t       &texObj,
                          const struct texture<T, 1, cudaReadModeElementType> *texRef,
                          const T                   *h_ptr,
                          int                        numElem,
                          const gmx_device_info_t   *devInfo)
{
    const size_t sizeInBytes = numElem * sizeof(*d_ptr);
    cudaError_t  stat        = cudaMalloc((void **)&d_ptr, sizeInBytes);
    CU_RET_ERR(stat, "cudaMalloc failed in initParamLookupTable");
    cu_copy_H2D(d_ptr, (void *)h_ptr, sizeInBytes);

    if (!c_disableCudaTextures)
    {
        if (use_texobj(devInfo))
        {
            setup1DTexture<T>(texObj, d_ptr, sizeInBytes);
        }
        else
        {
            setup1DTexture<T>(texRef, d_ptr, sizeInBytes);
        }
    }
}

template <typename T>
void destroyParamLookupTable(T                       *d_ptr,
                             cudaTextureObject_t      texObj,
                             const struct texture<T, 1, cudaReadModeElementType> *texRef,
                             const gmx_device_info_t *devInfo)
{
    if (!c_disableCudaTextures)
    {
        if (use_texobj(devInfo))
        {
            CU_RET_ERR(cudaDestroyTextureObject(texObj), "cudaDestroyTextureObject on texObj failed");
        }
        else
        {
            CU_RET_ERR(cudaUnbindTexture(texRef), "cudaUnbindTexture on texRef failed");
        }
    }
    CU_RET_ERR(cudaFree(d_ptr), "cudaFree failed");
}

/*! \brief Add explicit instantiations of init/destroyParamLookupTable() here as needed.
 * One should also verify that the result of cudaCreateChannelDesc<T>() during texture setup
 * looks reasonable, when instantiating the templates for new types - just in case.
 */
template void initParamLookupTable<float>(float * &, cudaTextureObject_t &, const texture<float, 1, cudaReadModeElementType> *, const float *, int, const gmx_device_info_t *);
template void destroyParamLookupTable<float>(float *, cudaTextureObject_t, const texture<float, 1, cudaReadModeElementType> *, const gmx_device_info_t *);
template void initParamLookupTable<int>(int * &, cudaTextureObject_t &, const texture<int, 1, cudaReadModeElementType> *, const int *, int, const gmx_device_info_t *);
template void destroyParamLookupTable<int>(int *, cudaTextureObject_t, const texture<int, 1, cudaReadModeElementType> *, const gmx_device_info_t *);
