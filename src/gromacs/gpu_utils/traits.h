#ifndef TRAITS_H
#define TRAITS_H //FIXME

#include "config.h"

#if GMX_GPU == GMX_GPU_CUDA

template <typename ValueType>
using DeviceBuffer = ValueType *;

#endif

#if GMX_GPU == GMX_GPU_OPENCL

#define OPENCL_HOST_PASS !defined(__OPENCL_C_VERSION__)

#if OPENCL_HOST_PASS
/*!
 *\brief
 * This allows to use the familiar data types without the cl_ prefix in the host code.
 * Add more if needed.
 */
using float4 = cl_float4;

#endif

template <typename ValueType>
using DeviceBuffer = cl_mem;

#endif

#if GMX_GPU == GMX_GPU_NONE
#error "No reason to include this file in non-GPU builds"
#endif

#endif
