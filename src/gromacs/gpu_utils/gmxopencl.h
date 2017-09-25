/*! \libinternal \file
 * \brief
 * Wraps the complexity of including OpenCL in Gromacs.
 *
 * Because OpenCL 2.0 is not officially supported widely, \Gromacs
 * uses earlier interfaces. Some of those have been deprecated in 2.0,
 * and generate warnings, which we need to suppress.
 *
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GMXOPENCL_H
#define GMX_GPU_UTILS_GMXOPENCL_H

/*! \brief Declare to OpenCL SDKs that we intend to use OpenCL API
   features that were deprecated in 2.0, so that they don't warn about
   it. */
#  define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#  ifdef __APPLE__
#    include <OpenCL/opencl.h>
#  else
#    include <CL/opencl.h>
#  endif

#endif
