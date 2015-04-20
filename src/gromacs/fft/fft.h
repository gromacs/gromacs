/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2003 David van der Spoel, Erik Lindahl, University of Groningen.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \file
 *  \brief Fast Fourier Transforms.
 *
 *  This file provides an abstract Gromacs interface to Fourier transforms,
 *  including multi-dimensional and real-to-complex transforms.
 *
 *  Internally it is implemented as wrappers to external libraries such
 *  as FFTW or the Intel Math Kernel Library, but we also have a built-in
 *  version of FFTPACK in case the faster alternatives are unavailable.
 *
 *  We also provide our own multi-dimensional transform setups even when
 *  the underlying library does not support it directly.
 *
 * \inpublicapi
 * \ingroup module_fft
 */
#ifndef GMX_FFT_FFT_H
#define GMX_FFT_FFT_H

#include <stdio.h>

#include "gromacs/math/gmxcomplex.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif



/*! \brief Datatype for FFT setup
 *
 *  The gmx_fft_t type contains all the setup information, e.g. twiddle
 *  factors, necessary to perform an FFT. Internally it is mapped to
 *  whatever FFT library we are using, or the built-in FFTPACK if no fast
 *  external library is available.
 *
 *  Since some of the libraries (e.g. MKL) store work array data in their
 *  handles this datatype should only be used for one thread at a time, i.e.
 *  they should allocate one instance each when executing in parallel.
 */
typedef struct gmx_fft *
    gmx_fft_t;




/*! \brief Specifier for FFT direction.
 *
 *  The definition of the 1D forward transform from input x[] to output y[] is
 *  \f[
 *  y_{k} = \sum_{j=0}^{N-1} x_{j} \exp{-i 2 \pi j k /N}
 *  \f]
 *
 *  while the corresponding backward transform is
 *
 *  \f[
 *  y_{k} = \sum_{j=0}^{N-1} x_{j} \exp{i 2 \pi j k /N}
 *  \f]
 *
 *  A forward-backward transform pair will this result in data scaled by N.
 *
 *  For complex-to-complex transforms you can only use one of
 *  GMX_FFT_FORWARD or GMX_FFT_BACKWARD, and for real-complex transforms you
 *  can only use GMX_FFT_REAL_TO_COMPLEX or GMX_FFT_COMPLEX_TO_REAL.
 */
typedef enum gmx_fft_direction
{
    GMX_FFT_FORWARD,         /**< Forward complex-to-complex transform  */
    GMX_FFT_BACKWARD,        /**< Backward complex-to-complex transform */
    GMX_FFT_REAL_TO_COMPLEX, /**< Real-to-complex valued FFT            */
    GMX_FFT_COMPLEX_TO_REAL  /**< Complex-to-real valued FFT            */
} gmx_fft_direction;

/*! \brief Specifier for FFT flags.
 *
 *  Some FFT libraries (FFTW, in particular) can do timings and other
 *  tricks to try and optimize the FFT for the current architecture. However,
 *  this can also lead to results that differ between consecutive runs with
 *  identical input.
 *  To avoid this, the conservative flag will attempt to disable such
 *  optimization, but there are no guarantees since we cannot control what
 *  the FFT libraries do internally.
 */

typedef int gmx_fft_flag;
/** Macro to indicate no special flags for FFT routines. */
static const int GMX_FFT_FLAG_NONE = 0;
/** Flag to disable FFT optimizations based on timings, see ::gmx_fft_flag. */
static const int GMX_FFT_FLAG_CONSERVATIVE = (1<<0);

/*! \brief Setup a 1-dimensional complex-to-complex transform
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform
 *  \param flags  FFT options
 *
 *  \return status - 0 or a standard error message.
 *
 *  \note Since some of the libraries (e.g. MKL) store work array data in their
 *        handles this datatype should only be used for one thread at a time,
 *        i.e. you should create one copy per thread when executing in parallel.
 */
int
gmx_fft_init_1d        (gmx_fft_t *       fft,
                        int               nx,
                        gmx_fft_flag      flags);


/*! \brief Setup multiple 1-dimensional complex-to-complex transform
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform
 *  \param howmany Howmany 1D FFT
 *  \param flags  FFT options
 *
 *  \return status - 0 or a standard error message.
 *
 *  \note Since some of the libraries (e.g. MKL) store work array data in their
 *        handles this datatype should only be used for one thread at a time,
 *        i.e. you should create one copy per thread when executing in parallel.
 */
int
gmx_fft_init_many_1d        (gmx_fft_t *       fft,
                             int               nx,
                             int               howmany,
                             gmx_fft_flag      flags);


/*! \brief Setup a 1-dimensional real-to-complex transform
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform in real space
 *  \param flags  FFT options
 *
 *  \return status - 0 or a standard error message.
 *
 *  \note Since some of the libraries (e.g. MKL) store work array data in their
 *        handles this datatype should only be used for one thread at a time,
 *        i.e. you should create one copy per thread when executing in parallel.
 */
int
gmx_fft_init_1d_real        (gmx_fft_t *       fft,
                             int               nx,
                             gmx_fft_flag      flags);


/*! \brief Setup multiple 1-dimensional real-to-complex transform
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform in real space
 *  \param howmany Homany 1D FFTs
 *  \param flags  FFT options
 *
 *  \return status - 0 or a standard error message.
 *
 *  \note Since some of the libraries (e.g. MKL) store work array data in their
 *        handles this datatype should only be used for one thread at a time,
 *        i.e. you should create one copy per thread when executing in parallel.
 */
int
gmx_fft_init_many_1d_real        (gmx_fft_t *       fft,
                                  int               nx,
                                  int               howmany,
                                  gmx_fft_flag      flags);


/*! \brief Setup a 2-dimensional real-to-complex transform
 *
 *  \param fft    Pointer to opaque Gromacs FFT datatype
 *  \param nx     Length of transform in first dimension
 *  \param ny     Length of transform in second dimension
 *  \param flags  FFT options
 *
 *  The normal space is assumed to be real, while the values in
 *  frequency space are complex.
 *
 *  \return status - 0 or a standard error message.
 *
 *  \note Since some of the libraries (e.g. MKL) store work array data in their
 *        handles this datatype should only be used for one thread at a time,
 *        i.e. you should create one copy per thread when executing in parallel.
 */
int
gmx_fft_init_2d_real        (gmx_fft_t *         fft,
                             int                 nx,
                             int                 ny,
                             gmx_fft_flag        flags);


/*! \brief Perform a 1-dimensional complex-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from gmx_fft_init_1d()
 *  \param dir       Forward or Backward
 *  \param in_data   Input grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *  \param out_data  Output grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers
 *       depending on your grid type.
 */
int
gmx_fft_1d               (gmx_fft_t                  setup,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data);


/*! \brief Perform many 1-dimensional complex-to-complex transforms
 *
 *  Performs many instances of a transform previously initiated.
 *
 *  \param setup     Setup returned from gmx_fft_init_1d()
 *  \param dir       Forward or Backward
 *  \param in_data   Input grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *  \param out_data  Output grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers
 *       depending on your grid type.
 */
int
gmx_fft_many_1d          (gmx_fft_t                  setup,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data);


/*! \brief Perform a 1-dimensional real-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from gmx_fft_init_1d_real()
 *  \param dir       Real-to-complex or complex-to-real
 *  \param in_data   Input grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *  \param out_data  Output grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * If you are doing an in-place transform, the array must be padded up to
 * an even integer length so n/2 complex numbers can fit. Out-of-place arrays
 * should not be padded (although it doesn't matter in 1d).
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers
 *       depending on transform direction.
 */
int
gmx_fft_1d_real          (gmx_fft_t                  setup,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data);

/*! \brief Perform many 1-dimensional real-to-complex transforms
 *
 *  Performs many instances of a transform previously initiated.
 *
 *  \param setup     Setup returned from gmx_fft_init_1d_real()
 *  \param dir       Real-to-complex or complex-to-real
 *  \param in_data   Input grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *  \param out_data  Output grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * If you are doing an in-place transform, the array must be padded up to
 * an even integer length so n/2 complex numbers can fit. Out-of-place arrays
 * should not be padded (although it doesn't matter in 1d).
 *
 * \return 0 on success, or an error code.
 *
 * \note Data pointers are declared as void, to avoid casting pointers
 *       depending on transform direction.
 */
int
gmx_fft_many_1d_real     (gmx_fft_t                  setup,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data);

/*! \brief Perform a 2-dimensional real-to-complex transform
 *
 *  Performs an instance of a transform previously initiated.
 *
 *  \param setup     Setup returned from gmx_fft_init_1d_real()
 *  \param dir       Real-to-complex or complex-to-real
 *  \param in_data   Input grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *  \param out_data  Output grid data. This should be allocated with gmx_new()
 *                   to make it 16-byte aligned for better performance.
 *                   You can provide the same pointer for in_data and out_data
 *                   to perform an in-place transform.
 *
 * \return 0 on success, or an error code.
 *
 * \note If you are doing an in-place transform, the last dimension of the
 * array MUST be padded up to an even integer length so n/2 complex numbers can
 * fit. Thus, if the real grid e.g. has dimension 5*3, you must allocate it as
 * a 5*4 array, where the last element in the second dimension is padding.
 * The complex data will be written to the same array, but since that dimension
 * is 5*2 it will now fill the entire array. Reverse complex-to-real in-place
 * transformation will produce the same sort of padded array.
 *
 * The padding does NOT apply to out-of-place transformation. In that case the
 * input array will simply be 5*3 of real, while the output is 5*2 of complex.
 *
 * \note Data pointers are declared as void, to avoid casting pointers
 *       depending on transform direction.
 */
int
gmx_fft_2d_real          (gmx_fft_t                  setup,
                          enum gmx_fft_direction     dir,
                          void *                     in_data,
                          void *                     out_data);

/*! \brief Release an FFT setup structure
 *
 *  Destroy setup and release all allocated memory.
 *
 *  \param setup Setup returned from gmx_fft_init_1d(), or one
 *		 of the other initializers.
 *
 */
void
gmx_fft_destroy          (gmx_fft_t                 setup);

/*! \brief Release a many FFT setup structure
 *
 *  Destroy setup and release all allocated memory.
 *
 *  \param setup Setup returned from gmx_fft_init_1d(), or one
 *		 of the other initializers.
 *
 */
void
gmx_many_fft_destroy          (gmx_fft_t                 setup);


/*! \brief Transpose 2d complex matrix, in-place or out-of-place.
 *
 * This routines works when the matrix is non-square, i.e. nx!=ny too,
 * without allocating an entire matrix of work memory, which is important
 * for huge FFT grids.
 *
 * \param in_data    Input data, to be transposed
 * \param out_data   Output, transposed data. If this is identical to
 *                   in_data, an in-place transpose is performed.
 * \param nx         Number of rows before transpose
 * \param ny         Number of columns before transpose
 *
 * \return GMX_SUCCESS, or an error code from gmx_errno.h
 */
int
gmx_fft_transpose_2d   (t_complex *       in_data,
                        t_complex *       out_data,
                        int               nx,
                        int               ny);

/*! \brief Cleanup global data of FFT
 *
 *  Any plans are invalid after this function. Should be called
 *  after all plans have been destroyed.
 */
void gmx_fft_cleanup();

/*! \brief Return string describing the underlying FFT implementation.
 *
 * Used to print out information about the used FFT library where needed.
 */
const char *gmx_fft_get_version_info();

#ifdef __cplusplus
}
#endif

#endif
