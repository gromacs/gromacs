/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#ifndef GMX_GPU_UTILS_TEXTURESUPPORT
#define GMX_GPU_UTILS_TEXTURESUPPORT

#include "gromacs/gpu_utils/cuda_arch_utils.cuh"
#include "gromacs/gpu_utils/cudautils.cuh"

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
struct TextureWrapperType
{
    T *d_pointer;
#if DISABLE_CUDA_TEXTURES
#elif GMX_PTX_ARCH >= 300
    cudaTextureObject_t                     object;
#else
    texture<T, 1, cudaReadModeElementType> *reference;
#endif
};

template <typename T>
void fillTextureWrapper(TextureWrapperType<T> *wrapper,
                        const T               *h_pointer,
                        int                    numElem)
{
    const size_t sizeInBytes = numElem * sizeof(T);
    cudaError_t  stat        = cudaMalloc((void **)&wrapper->d_pointer, sizeInBytes);
    CU_RET_ERR(stat, "cudaMalloc failed in fillTextureWrapper");
    cu_copy_H2D(wrapper->d_pointer, (void *)h_pointer, sizeInBytes);

    if (!c_disableCudaTextures)
    {
#if !DISABLE_CUDA_TEXTURES
#  if GMX_PTX_ARCH >= 300
        setup1DTexture<T>(wrapper->object, wrapper->d_pointer, sizeInBytes);
#  else
        setup1DTexture<T>(wrapper->reference, wrapper->d_pointer, sizeInBytes);
#  endif
#endif
    }
}

template <typename T>
void destroyTextureWrapper(TextureWrapperType<T> *wrapper)
{

    if (!c_disableCudaTextures)
    {
#if !DISABLE_CUDA_TEXTURES
#  if GMX_PTX_ARCH >= 300
        CU_RET_ERR(cudaDestroyTextureObject(wrapper->object), "cudaDestroyTextureObject on wrapper->object failed");
#  else
        CU_RET_ERR(cudaUnbindTexture(wrapper->reference), "cudaUnbindTexture on wrapper->reference failed");
#  endif
#endif
    }
    CU_RET_ERR(cudaFree(wrapper->d_pointer), "cudaFree failed");
}

#endif
