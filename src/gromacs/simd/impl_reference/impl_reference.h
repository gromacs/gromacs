/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_REFERENCE_H
#define GMX_SIMD_IMPL_REFERENCE_H

/*! \libinternal \file
 *
 * \brief Reference SIMD implementation, including SIMD documentation.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */


#include <math.h>

#include "gromacs/utility/fatalerror.h"

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name SIMD implementation capability definitions
 *  \{
 */

/*! \brief Temporary indication of second-generation extended Gromacs SIMD.
 *
 * This define will be removed when the new intrinsics and routines are
 * implemented for all architectures. Used for testing.
 */
#define GMX_SIMD_EXTENDED

/*! \brief
 * Defined when SIMD float support is present.
 *
 * You should only use this to specifically check for single precision SIMD,
 * support, even when the rest of Gromacs uses double precision.
 * \sa GMX_SIMD_HAVE_REAL, GMX_SIMD_HAVE_DOUBLE
 */
#define GMX_SIMD_HAVE_FLOAT

/*! \brief Defined if SIMD double support is present. */
#define GMX_SIMD_HAVE_DOUBLE

/*! \brief Defined if SIMD is implemented with real hardware instructions. */
#define GMX_SIMD_HAVE_HARDWARE
#undef  GMX_SIMD_HAVE_HARDWARE /* Reference implementation setting - previous line was for Doxygen. */

/*! \brief Defined if the SIMD implementation supports unaligned loads. */
#define GMX_SIMD_HAVE_LOADU

/*! \brief Defined if the SIMD implementation supports unaligned stores. */
#define GMX_SIMD_HAVE_STOREU

/*! \brief Defined if SIMD implementation has logical operations on floating-point data. */
#define GMX_SIMD_HAVE_LOGICAL

/*! \brief Defined if SIMD fused multiply-add uses hardware instructions */
#define GMX_SIMD_HAVE_FMA
#undef  GMX_SIMD_HAVE_FMA  /* Reference implementation setting - previous line was for Doxygen. */

/*! \brief Defined if the SIMD fraction has a direct hardware instruction. */
#define GMX_SIMD_HAVE_FRACTION
#undef  GMX_SIMD_HAVE_FRACTION /* Reference implementation setting - previous line was for Doxygen. */

/*! \brief Support for extracting integers from \ref gmx_simd_fint32_t. */
#define GMX_SIMD_HAVE_FINT32_EXTRACT

/*! \brief Defined if SIMD logical operations are supported for \ref gmx_simd_fint32_t */
#define GMX_SIMD_HAVE_FINT32_LOGICAL

/*! \brief Defined if SIMD arithmetic operations are supported for \ref gmx_simd_fint32_t */
#define GMX_SIMD_HAVE_FINT32_ARITHMETICS

/*! \brief Support for extracting integer from \ref gmx_simd_dint32_t */
#define GMX_SIMD_HAVE_DINT32_EXTRACT

/*! \brief Defined if logical operations are supported for \ref gmx_simd_dint32_t */
#define GMX_SIMD_HAVE_DINT32_LOGICAL

/*! \brief Defined if SIMD arithmetic operations are supported for \ref gmx_simd_dint32_t */
#define GMX_SIMD_HAVE_DINT32_ARITHMETICS

/*! \brief Defined if implementation provides an efficient \ref gmx_simd_loadu_2_transpose_f */
#define GMX_SIMD_HAVE_LOADU_2_TRANSPOSE_FLOAT

/*! \brief Defined if implementation provides an efficient \ref gmx_simd_loadu_2_transpose_d */
#define GMX_SIMD_HAVE_LOADU_2_TRANSPOSE_DOUBLE

/*! \brief Defined if implementation provides float half-register load/store/reduce utils */
#define GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL

/*! \brief Defined if implementation provides double half-register load/store/reduce utils */
#define GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL

/*! \brief Defined if the implementation provides \ref gmx_simd4_float_t. */
#define GMX_SIMD4_HAVE_FLOAT

/*! \brief Defined if the implementation provides \ref gmx_simd4_double_t. */
#define GMX_SIMD4_HAVE_DOUBLE

#ifdef GMX_SIMD_REF_FLOAT_WIDTH
#    define GMX_SIMD_FLOAT_WIDTH             GMX_SIMD_REF_FLOAT_WIDTH
#else
/*! \brief Width of the \ref gmx_simd_float_t datatype. */
#    define GMX_SIMD_FLOAT_WIDTH             4
#endif

#ifdef GMX_SIMD_REF_DOUBLE_WIDTH
#    define GMX_SIMD_DOUBLE_WIDTH            GMX_SIMD_REF_DOUBLE_WIDTH
#else
/*! \brief Width of the \ref gmx_simd_double_t datatype. */
#    define GMX_SIMD_DOUBLE_WIDTH            4
#endif

/*! \brief Width of the \ref gmx_simd_fint32_t datatype. */
#define GMX_SIMD_FINT32_WIDTH            GMX_SIMD_FLOAT_WIDTH

/*! \brief Width of the \ref gmx_simd_dint32_t datatype. */
#define GMX_SIMD_DINT32_WIDTH            GMX_SIMD_DOUBLE_WIDTH

/*! \brief Accuracy of SIMD 1/sqrt(x) lookup. Used to determine number of iterations. */
#define GMX_SIMD_RSQRT_BITS             23

/*! \brief Accuracy of SIMD 1/x lookup. Used to determine number of iterations. */
#define GMX_SIMD_RCP_BITS               23

/*! \}
 *
 * \name SIMD implementation data types
 * \{
 */
/*! \libinternal \brief Float SIMD variable. Supported with GMX_SIMD_HAVE_FLOAT.
 */
typedef struct
{
    float r[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_float_t;

/*! \libinternal \brief Floating-point SIMD variable type in double precision.
 *
 * Supported with GMX_SIMD_HAVE_DOUBLE.
 */
typedef struct
{
    double r[GMX_SIMD_DOUBLE_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_double_t;

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from float.
 *
 * This is also the widest integer SIMD type. Available with GMX_SIMD_HAVE_FLOAT.
 *
 * \note The integer SIMD type will always be available, but on architectures
 * that do not have any real integer SIMD support it might be defined as the
 * floating-point type. This will work fine, since there are separate defines
 * for whether the implementation can actually do any operations on integer
 * SIMD types.
 */
typedef struct
{
    gmx_int32_t i[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fint32_t;

/*! \libinternal \brief Integer SIMD variable type to use for conversions to/from double.
 *
 * Available with GMX_SIMD_HAVE_DOUBLE.
 *
 * \note The integer SIMD type will always be available, but on architectures
 * that do not have any real integer SIMD support it might be defined as the
 * floating-point type. This will work fine, since there are separate defines
 * for whether the implementation can actually do any operations on integer
 * SIMD types.
 *
 * \note The Gromacs SIMD module works entirely with 32 bit integers, both
 * in single and double precision, since some platforms do not support 64 bit
 * SIMD integers at all. In particular, this means it is up to each
 * implementation to get this working even if the architectures internal
 * representation uses 64 bit integers when converting to/from double SIMD
 * variables. For now we will try HARD to use conversions, packing or shuffling
 * so the integer datatype has the same width as the floating-point type, i.e.
 * if you use double precision SIMD with a width of 8, we want the integers
 * we work with to also use a SIMD width of 8 to make it easy to load/store
 * indices from arrays. This refers entirely to the function calls
 * and how many integers we load/store in one call; the actual SIMD registers
 * might be wider for integers internally (e.g. on x86 gmx_simd_dint32_t will
 * only fill half the register), but this is none of the user's business.
 * While this works for all current architectures, and we think it will work
 * for future ones, we might have to alter this decision in the future. To
 * avoid rewriting every single instance that refers to the SIMD width we still
 * provide separate defines for the width of SIMD integer variables that you
 * should use.
 */
typedef struct
{
    gmx_int32_t i[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dint32_t;

/*! \libinternal \brief Boolean type for float SIMD data.
 *
 * You should likely use gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_FLOAT_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fbool_t;

/*! \libinternal \brief Boolean type for double precision SIMD data.
 *
 * Use the generic gmx_simd_bool_t
 * (for gmx_simd_real_t) instead, unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_DOUBLE_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dbool_t;

/*! \libinternal \brief Boolean type for integer datatypes corresponding to float SIMD. */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_FINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_fibool_t;

/*! \libinternal \brief Boolean type for integer datatypes corresponding to double SIMD.
 *
 * You should likely use gmx_simd_ibool_t (for gmx_simd_int32_t) instead,
 * unless you really know what you are doing.
 */
typedef struct
{
    gmx_int32_t b[GMX_SIMD_DINT32_WIDTH]; /**< Implementation dependent. Don't touch. */
}
gmx_simd_dibool_t;


/*! \}
 *
 * \name SIMD implementation load/store operations for single precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_FLOAT_WIDTH numbers from aligned memory.
 *
 * \param m Pointer to memory aligned to the SIMD width.
 * \return SIMD variable with data loaded.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_load_f(const float *m)
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Set all SIMD variable elements to float pointed to by m (unaligned).
 *
 * \param m Pointer to single value in memory.
 * \return SIMD variable with all elements set to *m.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_load1_f(const float *m)
{
    gmx_simd_float_t  a;
    int               i;
    float             f = *m;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = f;
    }
    return a;
}

/*! \brief Set all SIMD float variable elements to the value r.
 *
 *  \param r floating-point constant
 *  \return SIMD variable with all elements set to r.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_set1_f(float r)
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}

/*! \brief Set all SIMD float variable elements to 0.0f.
 *
 *  \return The value 0.0 in all elements of a SIMD variable.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_setzero_f()
{
    gmx_simd_float_t  a;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }
    return a;
}

/*! \brief Store the contents of the SIMD float variable pr to aligned memory m.
 *
 * \param[out] m Pointer to memory, aligned to SIMD width.
 * \param a SIMD variable to store
 */
static gmx_inline void
gmx_simd_store_f(float *m, gmx_simd_float_t a)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD float from unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirement.
 * \return SIMD variable with data loaded.
 */
#define gmx_simd_loadu_f gmx_simd_load_f

/*! \brief Store SIMD float to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a SIMD variable to store.
 */
#define gmx_simd_storeu_f gmx_simd_store_f

/*! \}
 *
 * \name SIMD implementation load/store operations for double precision floating point
 * \{
 */

/*! \brief Load \ref GMX_SIMD_DOUBLE_WIDTH numbers from aligned memory.
 *
 * \copydetails gmx_simd_load_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load_d(const double *m)
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = m[i];
    }
    return a;
}

/*! \brief Set all SIMD variable elements to double pointed to by m (unaligned).
 *
 * \copydetails gmx_simd_load1_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load1_d(const double *m)
{
    gmx_simd_double_t  a;
    int                i;
    double             d = *m;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = d;
    }
    return a;
}

/*! \brief Set all SIMD double variable elements to the value r.
 *
 * \copydetails gmx_simd_set1_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_set1_d(double r)
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = r;
    }
    return a;
}

/*! \brief Set all SIMD double variable elements to 0.0.
 *
 * \copydetails gmx_simd_setzero_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_setzero_d()
{
    gmx_simd_double_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        a.r[i] = 0.0;
    }
    return a;
}

/*! \brief Store the contents of the SIMD double variable pr to aligned memory m.
 *
 * \copydetails gmx_simd_store_f
 */
static gmx_inline void
gmx_simd_store_d(double *m, gmx_simd_double_t a)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        m[i] = a.r[i];
    }
}

/*! \brief Load SIMD double from unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \copydetails gmx_simd_loadu_f
 */
#define gmx_simd_loadu_d gmx_simd_load_d

/*! \brief Store SIMD double to unaligned memory.
 *
 * Available with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \copydetails gmx_simd_storeu_f
 */
#define gmx_simd_storeu_d gmx_simd_store_d

/*! \}
 *
 * \name SIMD implementation load/store operations for integers (corresponding to float)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_load_i.
 *
 * \param m Pointer to memory, aligned to integer SIMD width.
 * \return SIMD integer variable.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_load_fi(const gmx_int32_t * m)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_set1_i.
 *
 *  \param b integer value to set variable to.
 *  \return SIMD variable with all elements set to b.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_set1_fi(gmx_int32_t b)
{
    gmx_simd_fint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_setzero_i.
 *
 * \return SIMD integer variable with all bits set to zero.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_setzero_fi()
{
    gmx_simd_fint32_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_store_i.
 *
 * \param m Memory aligned to integer SIMD width.
 * \param a SIMD variable to store.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_store_fi(int * m, gmx_simd_fint32_t a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_loadu_i.
 *
 * Supported with \ref GMX_SIMD_HAVE_LOADU.
 *
 * \param m Pointer to memory, no alignment requirements.
 * \return SIMD integer variable.
 */
#define gmx_simd_loadu_fi  gmx_simd_load_fi

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx_simd_float_t.
 *
 * You should typically call the real-precision \ref gmx_simd_storeu_i.
 *
 * Supported with \ref GMX_SIMD_HAVE_STOREU.
 *
 * \param m Memory pointer, no alignment requirements.
 * \param a SIMD variable to store.
 */
#define gmx_simd_storeu_fi gmx_simd_store_fi

/*! \brief Extract element with index i from \ref gmx_simd_fint32_t.
 *
 * You should typically call the real-precision \ref gmx_simd_extract_i.
 *
 * Available with \ref GMX_SIMD_HAVE_FINT32_EXTRACT.
 *
 * \param a SIMD variable
 * \param index Position to extract integer from
 * \return Single integer from position index in SIMD variable.
 */
static gmx_inline gmx_int32_t
gmx_simd_extract_fi(gmx_simd_fint32_t a, int index)
{
    return a.i[index];
}

/*! \}
 *
 * \name SIMD implementation load/store operations for integers (corresponding to double)
 * \{
 */

/*! \brief Load aligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_load_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_load_di(const gmx_int32_t * m)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = m[i];
    }
    return a;
};

/*! \brief Set SIMD from integer, width corresponds to \ref gmx_simd_double_t.
 *
 *  \copydetails gmx_simd_set1_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_set1_di(gmx_int32_t b)
{
    gmx_simd_dint32_t  a;
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = b;
    }
    return a;
}

/*! \brief Set all SIMD variable elements to 0, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_setzero_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_setzero_di()
{
    gmx_simd_dint32_t  a;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        a.i[i] = 0;
    }
    return a;
}

/*! \brief Store aligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_store_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_store_di(gmx_int32_t * m, gmx_simd_dint32_t a)
{
    int                i;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        m[i] = a.i[i];
    }
    return a;
};

/*! \brief Load unaligned integer SIMD data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_loadu_fi
 */
#define gmx_simd_loadu_di  gmx_simd_load_di

/*! \brief Store unaligned SIMD integer data, width corresponds to \ref gmx_simd_double_t.
 *
 * \copydetails gmx_simd_storeu_fi
 */
#define gmx_simd_storeu_di gmx_simd_store_di

/*! \brief Extract element with index i from \ref gmx_simd_dint32_t.
 *
 * \copydetails gmx_simd_extract_fi
 */
static gmx_inline gmx_int32_t
gmx_simd_extract_di(gmx_simd_dint32_t a, int index)
{
    return a.i[index];
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point bitwise logical operations
 * \{
 */

/*! \brief Bitwise and for two SIMD float variables. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_and_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 & data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_and_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
#ifdef __cplusplus
    gmx_int32_t       val1, val2, res;
#else
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 & val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise andnot for SIMD float. c=(~a) & b. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_andnot_r.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_andnot_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
#ifdef __cplusplus
    gmx_int32_t       val1, val2, res;
#else
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = (~val1) & val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = (~conv1.i) & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise or for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_or_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 | data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_or_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
#ifdef __cplusplus
    gmx_int32_t       val1, val2, res;
#else
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 | val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i | conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise xor for SIMD float. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_r.
 *
 * \param a data1
 * \param b data2
 * \return data1 ^ data2
 */
static gmx_inline gmx_simd_float_t
gmx_simd_xor_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;
#ifdef __cplusplus
    gmx_int32_t       val1, val2, res;
#else
    union
    {
        float        r;
        gmx_int32_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<int &>(a.r[i]);
        val2   = reinterpret_cast<int &>(b.r[i]);
        res    = val1 ^ val2;
        c.r[i] = reinterpret_cast<float &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i ^ conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point arithmetics
 * \{
 */

/*! \brief Add two float SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_add_r.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static gmx_inline gmx_simd_float_t
gmx_simd_add_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/*! \brief Subtract two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_sub_r.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static gmx_inline gmx_simd_float_t
gmx_simd_sub_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/*! \brief Multiply two SIMD variables.
 *
 * You should typically call the real-precision \ref gmx_simd_mul_r.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_mul_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = a.r[i] * b.r[i];
    }
    return c;
}

/*! \brief Fused-multiply-add. Result is a*b+c.
 *
 * You should typically call the real-precision \ref gmx_simd_fmadd_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fmadd_f(a, b, c) gmx_simd_add_f(gmx_simd_mul_f(a, b), c)


/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * You should typically call the real-precision \ref gmx_simd_fmsub_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_mul_f(a, b), c)


/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * You should typically call the real-precision \ref gmx_simd_fnmadd_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b+c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fnmadd_f(a, b, c) gmx_simd_sub_f(c, gmx_simd_mul_f(a, b))


/*! \brief Fused-negated-multiply-sub. Result is -a*b-c.
 *
 * You should typically call the real-precision \ref gmx_simd_fnmsub_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \return -a*b-c
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
#define gmx_simd_fnmsub_f(a, b, c) gmx_simd_sub_f(gmx_simd_setzero_f(), gmx_simd_fmadd_f(a, b, c))

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * You should typically call the real-precision \ref gmx_simd_rsqrt_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_rsqrt_f(gmx_simd_float_t x)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = 1.0f / sqrtf(x.r[i]);
    }
    return b;
};

/*! \brief SIMD 1.0/x lookup.
 *
 * You should typically call the real-precision \ref gmx_simd_rcp_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x!=0
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_rcp_f(gmx_simd_float_t x)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = 1.0f / x.r[i];
    }
    return b;
}

/*! \brief Multiply two SIMD variables, masked version.
 *
 * You should typically call the real-precision \ref gmx_simd_mul_r.
 *
 * \param a factor1
 * \param b factor2
 * \param m mask
 * \return a*b where mask is true, 0.0 otherwise.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_mul_mask_f(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_fbool_t m)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = m.b[i] ? (a.r[i] * b.r[i]) : 0.0;
    }
    return c;
}

/*! \brief Fused-multiply-add. Result is a*b+c, masked version.
 *
 * You should typically call the real-precision \ref gmx_simd_fmadd_r.
 *
 *  If \ref GMX_SIMD_HAVE_FMA is defined this is a single hardware instruction.
 *
 * \param a value
 * \param b value
 * \param c value
 * \param m mask
 * \return a*b+c where mask is true, 0.0 otherwise.
 *
 * For some implementations you save an instruction if you assign the result
 * to c.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fmadd_mask_f(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_float_t c,
                      gmx_simd_fbool_t m)
{
    gmx_simd_float_t  d;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d.r[i] = m.b[i] ? (a.r[i] * b.r[i] + c.r[i]) : 0.0;
    }
    return d;
}

/*! \brief SIMD 1.0/sqrt(x) lookup, masked version.
 *
 * You should typically call the real-precision \ref gmx_simd_rsqrt_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0 for entries where mask is true.
 * \param m Mask
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 *         The result for masked-out entries will be 0.0.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_rsqrt_mask_f(gmx_simd_float_t x, gmx_simd_fbool_t m)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (m.b[i] != 0) ? 1.0f / sqrtf(x.r[i]) : 0.0f;
    }
    return b;
}

/*! \brief SIMD 1.0/x lookup, masked version.
 *
 * You should typically call the real-precision \ref gmx_simd_rcp_r.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the reciprocal in simd_math.h.
 *
 * \param x Argument, x>0 for entries where mask is true.
 * \param m Mask
 * \return Approximation of 1/x, accuracy is \ref GMX_SIMD_RCP_BITS.
 *         The result for masked-out entries will be 0.0.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_rcp_mask_f(gmx_simd_float_t x, gmx_simd_fbool_t m)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = (m.b[i] != 0) ? 1.0f / x.r[i] : 0.0f;
    }
    return b;
}

/*! \brief SIMD Floating-point fabs().
 *
 * You should typically call the real-precision \ref gmx_simd_fabs_r.
 *
 * \param a any floating point values
 * \return fabs(a) for each element.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fabs_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = fabsf(a.r[i]);
    }
    return c;
}

/*! \brief SIMD floating-point negate.
 *
 * You should typically call the real-precision \ref gmx_simd_fneg_r.
 *
 * \param a Any floating-point value
 * \return -a
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fneg_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * You should typically call the real-precision \ref gmx_simd_max_r.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_max_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * You should typically call the real-precision \ref gmx_simd_min_r.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return min(a,b) for each element.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_min_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_float_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Round to nearest integer value (in floating-point format).
 *
 * You should typically call the real-precision \ref gmx_simd_round_r.
 *
 * \param a Any floating-point value
 * \return The nearest integer, represented in floating-point format.
 *
 * \note The reference implementation rounds exact half-way cases
 * away from zero, whereas most SIMD intrinsics will round to nearest even.
 * This could be fixed by using rint/rintf, but the bigger problem is that
 * MSVC does not support full C99, and none of the round or rint
 * functions are defined. It's much easier to approximately implement
 * round() than rint(), so we do that and hope we never get bitten in
 * testing. (Thanks, Microsoft.)
 */
static gmx_inline gmx_simd_float_t
gmx_simd_round_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
#ifdef _MSC_VER
        int temp = (a.r[i] >= 0.0f) ? (a.r[i] + 0.5f) : (a.r[i] - 0.5f);
        b.r[i] = temp;
#else
        b.r[i] = roundf(a.r[i]);
#endif
    }
    return b;
}

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * You should typically call the real-precision \ref gmx_simd_trunc_r.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This is truncation towards zero, not floor(). The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but floor() frequently isn't.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_trunc_f(gmx_simd_float_t a)
{
    gmx_simd_float_t  b;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.r[i] = truncf(a.r[i]);
    }
    return b;
}


/*! \brief Fraction of the SIMD floating point number.
 *
 * You should typically call the real-precision \ref gmx_simd_fraction_r.
 *
 * \param a Any floating-point value
 * \return a-trunc(r)
 *
 * To maximize compatibility, we use the same definition of fractions as used
 * e.g. for the AMD64 hardware instructions. This relies on truncation towards
 * zero for the integer part, and the remaining fraction can thus be either
 * positive or negative. As an example, -1.42 would return the fraction -0.42.
 *
 * Hardware support with \ref GMX_SIMD_HAVE_FRACTION, otherwise emulated.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_fraction_f(gmx_simd_float_t a)
{
    return gmx_simd_sub_f(a, gmx_simd_trunc_f(a));
}

/*! \brief Extract (integer) exponent from single precision SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_get_exponent_r.
 *
 * \param a Any floating-point value
 * \return Exponent value, represented in floating-point format.
 *
 * The IEEE754 exponent field is selected, the bias removed, and it is converted to
 * a normal floating-point SIMD.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_get_exponent_f(gmx_simd_float_t a)
{
    /* Mask with ones for the exponent field of single precision fp */
    const gmx_int32_t  expmask = 0x7f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* Keep exponent, shift 23 right (float mantissa), remove bias (127) */
        b.r[i] = ((conv.i & expmask) >> 23) - 127;
    }
    return b;
}

/*! \brief Get SIMD mantissa.
 *
 * You should typically call the real-precision \ref gmx_simd_get_mantissa_r.
 *
 * \param a Any floating-point value
 * \return Mantissa, represented in floating-point format.
 *
 * The mantissa field is selected, and a new neutral exponent created.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_get_mantissa_f(gmx_simd_float_t a)
{
    const gmx_int32_t  mantmask = 0x007fffff;
    const gmx_int32_t  one      = 0x3f800000;
    gmx_simd_float_t   b;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        conv.f = a.r[i];
        /* remove current exponent, add a biased exponent for 1.0 (i.e., 2^0=1) */
        conv.i = (conv.i & (mantmask)) | one;
        b.r[i] = conv.f;
    }
    return b;
}

/*! \brief Set (integer) exponent from single precision floating-point SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_set_exponent_r.
 *
 * \param a A floating point value that will not overflow as 2^a.
 * \return 2^(round(a)).
 *
 * The input is \a rounded to the nearest integer, the exponent bias is added
 * to this integer, and the bits are shifted to the IEEE754 exponent part of the number.
 *
 * \note The argument will be \a rounded to nearest integer since that is what
 * we need for the exponential functions, and this integer x will be set as the
 * exponent so the new floating-point number will be 2^x.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_set_exponent_f(gmx_simd_float_t a)
{
    gmx_simd_float_t   b;
    gmx_int32_t        iexp;
    int                i;
    union
    {
        float        f;
        gmx_int32_t  i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        /* Critical to use same algorithm as for gmx_simd_round_f() */
#ifdef _MSC_VER
        iexp = (a.r[i] >= 0.0f) ? (a.r[i] + 0.5f) : (a.r[i] - 0.5f);
#else
        iexp = roundf(a.r[i]);
#endif
        /* Add bias (127), and shift 23 bits left (mantissa size) */
        conv.i = (iexp + 127) << 23;
        b.r[i] = conv.f;
    }
    return b;
}

/*! \}
 *
 * \name SIMD implementation single precision floating-point comparisons, boolean, selection.
 * \{
 */
/*! \brief SIMD a==b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpeq_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a==b.
 *
 * Beware that exact floating-point comparisons are difficult.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmpeq_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/*! \brief SIMD a!=0 for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpnz_r.
 *
 * \param a value
 * \return Each element of the boolean will be true if any bit in a is nonzero.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmpnz_f(gmx_simd_float_t a)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] != 0.0);
    }
    return c;
}

/*! \brief SIMD a<b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmplt_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<b.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmplt_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/*! \brief SIMD a<=b for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmple_r.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<=b.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cmple_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    gmx_simd_fbool_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}

/*! \brief Logical \a and on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx_simd_and_r.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a \& b are true.
 *
 * \note This is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa gmx_simd_and_ib
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_and_fb(gmx_simd_fbool_t a, gmx_simd_fbool_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical \a or on single precision SIMD booleans.
 *
 * You should typically call the real-precision \ref gmx_simd_or_r.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a or b is true.
 *
 * Note that this is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 *
 * \sa gmx_simd_or_ib
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_or_fb(gmx_simd_fbool_t a, gmx_simd_fbool_t b)
{
    gmx_simd_fbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx_simd_anytrue_b.
 *
 * \param a Logical variable.
 * \return non-zero if any element in a is true, otherwise 0.
 *
 * The actual return value for truth will depend on the architecture,
 * so any non-zero value is considered truth.
 */
static gmx_inline int
gmx_simd_anytrue_fb(gmx_simd_fbool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/*! \brief Select from single precision SIMD variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx_simd_blendzero_r.
 *
 * \param a Floating-point variable to select from
 * \param sel Boolean selector
 * \return  For each element, a is selected for true, 0 for false.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_blendzero_f(gmx_simd_float_t a, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? a.r[i] : 0.0;
    }
    return c;
}

/*! \brief Select from single precision SIMD variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx_simd_blendnotzero_r.
 *
 * \param a Floating-point variable to select from
 * \param sel Boolean selector
 * \return  For each element, a is selected for false, 0 for true (sic).
 */
static gmx_inline gmx_simd_float_t
gmx_simd_blendnotzero_f(gmx_simd_float_t a, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t   c;
    int                i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? 0.0 : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx_simd_blendv_r.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_blendv_f(gmx_simd_float_t a, gmx_simd_float_t b, gmx_simd_fbool_t sel)
{
    gmx_simd_float_t  d;
    int               i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}

/*! \brief Return sum of all elements in SIMD float variable.
 *
 * You should typically call the real-precision \ref gmx_simd_reduce_r.
 *
 * \param a SIMD variable to reduce/sum.
 * \return The sum of all elements in the argument variable.
 *
 */
static gmx_inline float
gmx_simd_reduce_f(gmx_simd_float_t a)
{
    float     sum = 0.0;
    int       i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        sum += a.r[i];
    }
    return sum;
}

/*! \}
 *
 * \name SIMD implementation double precision floating-point bitwise logical operations
 * \{
 */
/*! \brief Bitwise and for two SIMD double variables. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_and_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_and_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
#ifdef __cplusplus
    gmx_int64_t        val1, val2, res;
#else
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 & val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise andnot for SIMD double. c=(~a) & b. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_andnot_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_andnot_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
#ifdef __cplusplus
    gmx_int64_t        val1, val2, res;
#else
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = (~val1) & val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise or for SIMD double. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_or_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_or_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
#ifdef __cplusplus
    gmx_int64_t        val1, val2, res;
#else
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 | val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \brief Bitwise xor for SIMD double. Supported with \ref GMX_SIMD_HAVE_LOGICAL.
 *
 * \copydetails gmx_simd_xor_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_xor_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;
#ifdef __cplusplus
    gmx_int64_t        val1, val2, res;
#else
    union
    {
        double       r;
        gmx_int64_t  i;
    }
    conv1, conv2;
#endif

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef __cplusplus
        val1   = reinterpret_cast<gmx_int64_t &>(a.r[i]);
        val2   = reinterpret_cast<gmx_int64_t &>(b.r[i]);
        res    = val1 ^ val2;
        c.r[i] = reinterpret_cast<double &>(res);
#else
        conv1.r = a.r[i];
        conv2.r = b.r[i];
        conv1.i = conv1.i & conv2.i;
        c.r[i]  = conv1.r;
#endif
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation double precision floating-point arithmetics
 * \{
 */
/*! \brief Add two double SIMD variables.
 *
 * \copydetails gmx_simd_add_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_add_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] + b.r[i];
    }
    return c;
}

/*! \brief Add two float SIMD variables.
 *
 * \copydetails gmx_simd_sub_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_sub_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] - b.r[i];
    }
    return c;
}

/*! \brief Multiply two SIMD variables.
 *
 * \copydetails gmx_simd_mul_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_mul_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = a.r[i] * b.r[i];
    }
    return c;
}

/*! \brief Fused-multiply-add, double. Result is a*b+c.
 *
 * \copydetails gmx_simd_fmadd_f
 */
#define gmx_simd_fmadd_d(a, b, c) gmx_simd_add_d(gmx_simd_mul_d(a, b), c)

/*! \brief Fused-multiply-subtract. Result is a*b-c.
 *
 * \copydetails gmx_simd_fmsub_f
 */
#define gmx_simd_fmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_mul_d(a, b), c)

/*! \brief Fused-negated-multiply-add. Result is -a*b+c.
 *
 * \copydetails gmx_simd_fnmadd_f
 */
#define gmx_simd_fnmadd_d(a, b, c) gmx_simd_sub_d(c, gmx_simd_mul_d(a, b))

/*! \brief Fused-negated-multiply-add. Result is -a*b-c.
 *
 * \copydetails gmx_simd_fnmsub_f
 */
#define gmx_simd_fnmsub_d(a, b, c) gmx_simd_sub_d(gmx_simd_setzero_d(), gmx_simd_fmadd_d(a, b, c))

/*! \brief SIMD 1.0/sqrt(x) lookup.
 *
 * \copydetails gmx_simd_rsqrt_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_rsqrt_d(gmx_simd_double_t x)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Sic - we only need single precision for the reference lookup, since
         * we have defined GMX_SIMD_RSQRT_BITS to 23.
         */
        b.r[i] = 1.0 / sqrtf(x.r[i]);
    }
    return b;
};

/*! \brief 1.0/x lookup.
 *
 * \copydetails gmx_simd_rcp_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_rcp_d(gmx_simd_double_t x)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.r[i] = 1.0 / x.r[i];
    }
    return b;
};

/*! \brief Multiply two SIMD doubles, masked version.
 *
 * \copydetails gmx_simd_mul_mask_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_mul_mask_d(gmx_simd_double_t a, gmx_simd_double_t b, gmx_simd_dbool_t m)
{
    gmx_simd_double_t c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = m.b[i] ? (a.r[i] * b.r[i]) : 0.0;
    }
    return c;
}

/*! \brief Fused-multiply-add, double. Result is a*b+c, masked version.
 *
 * \copydetails gmx_simd_fmadd_mask_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fmadd_mask_d(gmx_simd_double_t a, gmx_simd_double_t b, gmx_simd_double_t c,
                      gmx_simd_dbool_t m)
{
    gmx_simd_double_t d;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = m.b[i] ? (a.r[i] * b.r[i] + c.r[i]) : 0.0;
    }
    return d;
}

/*! \brief SIMD 1.0/sqrt(x) lookup, masked version.
 *
 * \copydetails gmx_simd_rsqrt_mask_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_rsqrt_mask_d(gmx_simd_double_t x, gmx_simd_dbool_t m)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Sic - we only need single precision for the reference lookup, since
         * we have defined GMX_SIMD_RSQRT_BITS to 23.
         */
        b.r[i] = (m.b[i] != 0) ? 1.0 / sqrtf(x.r[i]) : 0.0;
    }
    return b;
}

/*! \brief 1.0/x lookup, masked version.
 *
 * \copydetails gmx_simd_rcp_mask_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_rcp_mask_d(gmx_simd_double_t x, gmx_simd_dbool_t m)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.r[i] = (m.b[i] != 0) ? 1.0 / x.r[i] : 0.0;
    }
    return b;
};

/*! \brief SIMD Floating-point fabs().
 *
 * \copydetails gmx_simd_fabs_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fabs_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = fabs(a.r[i]);
    }
    return c;
}

/*! \brief SIMD floating-point negate.
 *
 * \copydetails gmx_simd_fneg_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fneg_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = -a.r[i];
    }
    return c;
}

/*! \brief Set each SIMD element to the largest from two variables.
 *
 * \copydetails gmx_simd_max_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_max_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] >= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Set each SIMD element to the smallest from two variables.
 *
 * \copydetails gmx_simd_min_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_min_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = (a.r[i] <= b.r[i] ? a.r[i] : b.r[i]);
    }
    return c;
}

/*! \brief Round to nearest integer value (in double floating-point format).
 *
 * \copydetails gmx_simd_round_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_round_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
#ifdef _MSC_VER
        int temp = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
        b.r[i] = temp;
#else
        b.r[i] = round(a.r[i]);
#endif
    }
    return b;
}

/*! \brief Truncate SIMD, i.e. round towards zero - common hardware instruction.
 *
 * \copydetails gmx_simd_trunc_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_trunc_d(gmx_simd_double_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.r[i] = trunc(a.r[i]);
    }
    return b;
}

/*! \brief Fraction of the SIMD floating point number.
 *
 * \copydetails gmx_simd_fraction_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_fraction_d(gmx_simd_double_t a)
{
    return gmx_simd_sub_d(a, gmx_simd_trunc_d(a));
}


/*! \brief Extract (integer) exponent from double precision SIMD.
 *
 * \copydetails gmx_simd_get_exponent_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_get_exponent_d(gmx_simd_double_t a)
{
    /* Mask with ones for the exponent field of double precision fp */
    const gmx_int64_t      expmask = 0x7ff0000000000000LL;
    gmx_simd_double_t      b;
    int                    i;
    union
    {
        double             d;
        gmx_int64_t        i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv.d = a.r[i];
        /* Zero everything but exponent field (remove sign),
         * shift 23 bits right (mantissa size), and remove exponent bias (1023).
         */
        b.r[i] = ((conv.i & expmask) >> 52) - 1023;
    }
    return b;
}

/*! \brief Get SIMD doublemantissa.
 *
 * \copydetails gmx_simd_get_mantissa_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_get_mantissa_d(gmx_simd_double_t a)
{
    const gmx_int64_t      mantmask = 0x000fffffffffffffLL;
    const gmx_int64_t      one      = 0x3ff0000000000000LL;
    gmx_simd_double_t      b;
    int                    i;
    union
    {
        double          d;
        gmx_int64_t     i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        conv.d = a.r[i];
        conv.i = (conv.i & (mantmask)) | one;
        b.r[i] = conv.d;
    }
    return b;
}

/*! \brief Set (integer) exponent from single precision floating-point SIMD.
 *
 * \copydetails gmx_simd_set_exponent_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_set_exponent_d(gmx_simd_double_t a)
{
    gmx_simd_double_t      b;
    int                    i;
    gmx_int64_t            iexp;
    union
    {
        double          d;
        gmx_int64_t     i;
    }
    conv;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        /* Critical to use same algorithm as for gmx_simd_round_d() */
#ifdef _MSC_VER
        iexp = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        iexp = round(a.r[i]);
#endif
        /* Add bias (1023), and shift 52 bits left (mantissa size) */
        conv.i = (iexp + 1023) << 52;
        b.r[i] = conv.d;
    }
    return b;
}

/*! \}
 *
 * \name SIMD implementation double precision floating-point comparison, boolean, selection.
 * \{
 */
/*! \brief SIMD a==b for double SIMD.
 *
 * \copydetails gmx_simd_cmpeq_f
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmpeq_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] == b.r[i]);
    }
    return c;
}

/*! \brief SIMD a!=0 for single SIMD.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpnz_r.
 *
 * \param a value
 * \return Each element of the boolean will be true if any bit in a is nonzero.
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmpnz_d(gmx_simd_double_t a)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] != 0.0);
    }
    return c;
}

/*! \brief SIMD a<b for double SIMD.
 *
 * \copydetails gmx_simd_cmplt_f
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmplt_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] < b.r[i]);
    }
    return c;
}

/*! \brief SIMD a<=b for double SIMD.
 *
 * \copydetails gmx_simd_cmple_f
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cmple_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.r[i] <= b.r[i]);
    }
    return c;
}


/*! \brief Logical \a and on double precision SIMD booleans.
 *
 * \copydetails gmx_simd_and_fb
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_and_db(gmx_simd_dbool_t a, gmx_simd_dbool_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical \a or on double precision SIMD booleans.
 *
 * \copydetails gmx_simd_or_fb
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_or_db(gmx_simd_dbool_t a, gmx_simd_dbool_t b)
{
    gmx_simd_dbool_t  c;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}


/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * \copydetails gmx_simd_anytrue_fb
 */
static gmx_inline int
gmx_simd_anytrue_db(gmx_simd_dbool_t a)
{
    int         anytrue;
    int         i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}


/*! \brief Select from double SIMD variable where boolean is true.
 *
 * \copydetails gmx_simd_blendzero_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_blendzero_d(gmx_simd_double_t a, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? a.r[i] : 0.0;
    }
    return c;
}

/*! \brief Select from double SIMD variable where boolean is false.
 *
 * \copydetails gmx_simd_blendnotzero_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_blendnotzero_d(gmx_simd_double_t a, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        c.r[i] = sel.b[i] ? 0.0 : a.r[i];
    }
    return c;
}

/*! \brief Vector-blend double SIMD selection.
 *
 * \copydetails gmx_simd_blendv_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_blendv_d(gmx_simd_double_t a, gmx_simd_double_t b, gmx_simd_dbool_t sel)
{
    gmx_simd_double_t  d;
    int                i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = sel.b[i] ? b.r[i] : a.r[i];
    }
    return d;
}

/*! \brief Return sum of all elements in SIMD double variable.
 *
 * \copydetails gmx_simd_reduce_f
 *
 */
static gmx_inline double
gmx_simd_reduce_d(gmx_simd_double_t a)
{
    double    sum = 0.0;
    int       i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        sum += a.r[i];
    }
    return sum;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) bitwise logical operations
 * \{
 */

/*! \brief SIMD integer shift left logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx_simd_slli_i.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  left, while zeros are shifted in from the right. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is defined.
 *
 * \param a integer data to shift
 * \param n number of positions to shift left. n<=32.
 * \return shifted values
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_slli_fi(gmx_simd_fint32_t a, int n)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/*! \brief SIMD integer shift right logical, based on immediate value.
 *
 * You should typically call the real-precision \ref gmx_simd_srli_i.
 *
 *  Logical shift. Each element is shifted (independently) up to 32 positions
 *  right, while zeros are shifted in from the left. Only available if
 * \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single) or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
 *  (double) is defined.
 *
 * \param a integer data to shift
 * \param n number of positions to shift right. n<=32.
 * \return shifted values
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_srli_fi(gmx_simd_fint32_t a, int n)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/*! \brief Integer SIMD bitwise and.
 *
 * You should typically call the real-precision \ref gmx_simd_and_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is defined.
 *
 * \note You can \a not use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx_simd_blendzero_i instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \& b (bitwise and)
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_and_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise not-and.
 *
 * You should typically call the real-precision \ref gmx_simd_andnot_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is defined.
 *
 * Note that you can NOT use this operation directly to select based on a boolean
 * SIMD variable, since booleans are separate from integer SIMD. If that
 * is what you need, have a look at \ref gmx_simd_blendnotzero_i instead.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return (~a) \& b (bitwise andnot)
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_andnot_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise or.
 *
 * You should typically call the real-precision \ref gmx_simd_or_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is defined.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a \| b (bitwise or)
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_or_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/*! \brief Integer SIMD bitwise xor.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_LOGICAL (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_LOGICAL (double) is defined.
 *
 * \param a first integer SIMD
 * \param b second integer SIMD
 * \return a ^ b (bitwise xor)
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_xor_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] ^ b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) arithmetics
 * \{
 */
/*! \brief Add SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_add_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/*! \brief Subtract SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_sub_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/*! \brief Multiply SIMD integers.
 *
 * You should typically call the real-precision \ref gmx_simd_xor_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 *
 * \note Only the low 32 bits are retained, so this can overflow.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_mul_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] * b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to float) comparisons, boolean, selection
 * \{
 */

/*! \brief Equality comparison of two integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx_simd_cmpeq_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a==b
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cmpeq_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/*! \brief Less-than comparison of two SIMD integers corresponding to float values.
 *
 * You should typically call the real-precision \ref gmx_simd_cmplt_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD integer1
 * \param b SIMD integer2
 * \return SIMD integer boolean with true for elements where a<b
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cmplt_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/*! \brief Logical AND on gmx_simd_fibool_t.
 *
 * You should typically call the real-precision \ref gmx_simd_and_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_and_fib(gmx_simd_fibool_t a, gmx_simd_fibool_t b)
{
    gmx_simd_fibool_t c;
    int               i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical OR on gmx_simd_fibool_t.
 *
 * You should typically call the real-precision \ref gmx_simd_or_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD boolean 1
 * \param b SIMD boolean 2
 * \return True for elements where both a and b are true.
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_or_fib(gmx_simd_fibool_t a, gmx_simd_fibool_t b)
{
    gmx_simd_fibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns non-zero if any of the boolean in x is True, otherwise 0.
 *
 * You should typically call the real-precision \ref gmx_simd_anytrue_ib.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * The actual return value for "any true" will depend on the architecture.
 * Any non-zero value should be considered truth.
 *
 * \param a SIMD boolean
 * \return Nonzero integer if any of the elements in a is true, otherwise 0.
 */
static gmx_inline int
gmx_simd_anytrue_fib(gmx_simd_fibool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/*! \brief Select from \ref gmx_simd_fint32_t variable where boolean is true.
 *
 * You should typically call the real-precision \ref gmx_simd_blendzero_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD integer to select from
 * \param sel Boolean selector
 * \return Elements from a where sel is true, 0 otherwise.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendzero_fi(gmx_simd_fint32_t a, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? a.i[i] : 0.0;
    }
    return c;
}

/*! \brief Select from \ref gmx_simd_fint32_t variable where boolean is false.
 *
 * You should typically call the real-precision \ref gmx_simd_blendnotzero_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a SIMD integer to select from
 * \param sel Boolean selector
 * \return Elements from a where sel is false, 0 otherwise (sic).
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendnotzero_fi(gmx_simd_fint32_t a, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? 0.0 : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection.
 *
 * You should typically call the real-precision \ref gmx_simd_blendv_i.
 *
 * This routine is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS (single)
 *  or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS (double) is defined.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_blendv_fi(gmx_simd_fint32_t a, gmx_simd_fint32_t b, gmx_simd_fibool_t sel)
{
    gmx_simd_fint32_t d;
    int               i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        d.i[i] = sel.b[i] ? b.i[i] : a.i[i];
    }
    return d;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) bitwise logical operations
 * \{
 */

/*! \brief SIMD integer shift left, based on immediate value.
 *
 * \copydetails gmx_simd_slli_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_slli_di(gmx_simd_dint32_t a, int n)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] << n;
    }
    return c;
}

/*! \brief SIMD integer shift right, based on immediate value.
 *
 * \copydetails gmx_simd_srli_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_srli_di(gmx_simd_dint32_t a, int n)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] >> n;
    }
    return c;
}

/*! \brief Integer bitwise and for SIMD variables.
 *
 * \copydetails gmx_simd_and_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_and_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] & b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise not-and for SIMD variables.
 *
 * \copydetails gmx_simd_andnot_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_andnot_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = (~a.i[i]) & b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise or for SIMD variables.
 *
 * \copydetails gmx_simd_or_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_or_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] | b.i[i];
    }
    return c;
}

/*! \brief Integer bitwise xor for SIMD variables.
 *
 * \copydetails gmx_simd_xor_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_xor_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] ^ b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) arithmetics
 * \{
 */
/*! \brief Add SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_add_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_add_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] + b.i[i];
    }
    return c;
}

/*! \brief Subtract SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_sub_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_sub_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] - b.i[i];
    }
    return c;
}

/*! \brief Multiply SIMD integers, corresponding to double precision.
 *
 * \copydetails gmx_simd_mul_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_mul_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = a.i[i] * b.i[i];
    }
    return c;
}

/*! \}
 *
 * \name SIMD implementation integer (corresponding to double) comparisons, boolean selection
 * \{
 */

/*! \brief Equality comparison of two ints corresponding to double SIMD data.
 *
 * \copydetails gmx_simd_cmpeq_fi
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cmpeq_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] == b.i[i]);
    }
    return c;
}

/*! \brief Less-than comparison of two ints corresponding to double SIMD data.
 *
 * \copydetails gmx_simd_cmplt_fi
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cmplt_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.i[i] < b.i[i]);
    }
    return c;
}

/*! \brief Logical AND on gmx_simd_dibool_t.
 *
 * \copydetails gmx_simd_and_fib
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_and_dib(gmx_simd_dibool_t a, gmx_simd_dibool_t b)
{
    gmx_simd_dibool_t c;
    int               i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] && b.b[i]);
    }
    return c;
}

/*! \brief Logical OR on gmx_simd_dibool_t.
 *
 * \copydetails gmx_simd_or_fib
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_or_dib(gmx_simd_dibool_t a, gmx_simd_dibool_t b)
{
    gmx_simd_dibool_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.b[i] = (a.b[i] || b.b[i]);
    }
    return c;
}

/*! \brief Returns non-zero if any of the double-int SIMD booleans in x is True, otherwise 0.
 *
 * \copydetails gmx_simd_anytrue_fib
 */
static gmx_inline int
gmx_simd_anytrue_dib(gmx_simd_dibool_t a)
{
    int             anytrue;
    int             i;

    anytrue = 0;
    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        anytrue = anytrue || a.b[i];
    }
    return anytrue;
}

/*! \brief Select from SIMD ints (corresponding to double) where boolean is true.
 *
 * \copydetails gmx_simd_blendzero_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendzero_di(gmx_simd_dint32_t a, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? a.i[i] : 0.0;
    }
    return c;
}

/*! \brief Select from SIMD ints (corresponding to double) where boolean is false.
 *
 * \copydetails gmx_simd_blendnotzero_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendnotzero_di(gmx_simd_dint32_t a, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  c;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        c.i[i] = sel.b[i] ? 0.0 : a.i[i];
    }
    return c;
}

/*! \brief Vector-blend SIMD selection for double-int SIMD.
 *
 * \copydetails gmx_simd_blendv_fi
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_blendv_di(gmx_simd_dint32_t a, gmx_simd_dint32_t b, gmx_simd_dibool_t sel)
{
    gmx_simd_dint32_t  d;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        d.i[i] = sel.b[i] ? b.i[i] : a.i[i];
    }
    return d;
}

/*! \}
 *
 * \name SIMD implementation conversion operations
 * \{
 */

/*! \brief Round single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_r2i.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, rounded to nearest integer.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_cvt_f2i(gmx_simd_float_t a)
{
    gmx_simd_fint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
#ifdef _MSC_VER
        b.i[i] = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        b.i[i] = roundf(a.r[i]);
#endif
    }
    return b;
};

/*! \brief Truncate single precision floating point to integer.
 *
 * You should typically call the real-precision \ref gmx_simd_cvtt_r2i.
 *
 * \param a SIMD floating-point
 * \return SIMD integer, truncated towards zero.
 */
static gmx_inline gmx_simd_fint32_t
gmx_simd_cvtt_f2i(gmx_simd_float_t a)
{
    gmx_simd_fint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/*! \brief Convert integer to single precision floating-point.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_i2r.
 *
 * \param a SIMD integer
 * \return SIMD floating-pint
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_i2f(gmx_simd_fint32_t a)
{
    gmx_simd_float_t   b;
    int                i;

    for (i = 0; i < GMX_SIMD_FINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/*! \brief Round double precision floating point to integer.
 *
 * \copydetails gmx_simd_cvt_f2i
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_cvt_d2i(gmx_simd_double_t a)
{
    gmx_simd_dint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
#ifdef _MSC_VER
        b.i[i] = (a.r[i] >= 0.0) ? (a.r[i] + 0.5) : (a.r[i] - 0.5);
#else
        b.i[i] = round(a.r[i]);
#endif
    }
    return b;
};

/*! \brief Truncate double precision floating point to integer.
 *
 * \copydetails gmx_simd_cvtt_f2i
 */
static gmx_inline gmx_simd_dint32_t
gmx_simd_cvtt_d2i(gmx_simd_double_t a)
{
    gmx_simd_dint32_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.i[i] = a.r[i];
    }
    return b;
};

/*! \brief Convert integer to single precision floating-point.
 *
 * \copydetails gmx_simd_cvt_i2f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_cvt_i2d(gmx_simd_dint32_t a)
{
    gmx_simd_double_t  b;
    int                i;

    for (i = 0; i < GMX_SIMD_DINT32_WIDTH; i++)
    {
        b.r[i] = a.i[i];
    }
    return b;
};

/*! \brief Convert from float boolean to corresponding integer boolean.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_b2ib.
 *
 * \param a Boolean corresponding to SIMD floating-point
 * \return Boolean that can be applied to SIMD integer operations.
 */
static gmx_inline gmx_simd_fibool_t
gmx_simd_cvt_fb2fib(gmx_simd_fbool_t a)
{
    gmx_simd_fibool_t  b;
    int                i;

    /* Integer width >= float width */
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \brief Convert from integer boolean (corresponding to float) to float boolean.
 *
 * You should typically call the real-precision \ref gmx_simd_cvt_ib2b.
 *
 * \param a Boolean corresponding to SIMD integer
 * \return Boolean that can be applied to SIMD floating-point.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cvt_fib2fb(gmx_simd_fibool_t a)
{
    gmx_simd_fbool_t  b;
    int               i;

    /* Integer width >= float width */
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \brief Convert from double boolean to corresponding integer boolean.
 *
 * \copydetails gmx_simd_cvt_fb2fib
 */
static gmx_inline gmx_simd_dibool_t
gmx_simd_cvt_db2dib(gmx_simd_dbool_t a)
{
    gmx_simd_dibool_t  b;
    int                i;

    /* Integer width >= double width */
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \brief Convert from integer boolean (corresponding to double) to double boolean.
 *
 * \copydetails gmx_simd_cvt_fib2fb
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cvt_dib2db(gmx_simd_dibool_t a)
{
    gmx_simd_dbool_t  b;
    int               i;

    /* Integer width >= double width */
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        b.b[i] = a.b[i];
    }
    return b;
}

/*! \brief Convert SIMD float to double.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD variable
 * \return Double-precision SIMD variable of the same width
 */
static gmx_inline gmx_simd_double_t
gmx_simd_cvt_f2d(gmx_simd_float_t f)
{
    gmx_simd_double_t d;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int               i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.r[i] = f.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_f2d() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    d.r[0] = f.r[0];
#endif
    return d;
}

/*! \brief Convert SIMD double to float.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d Double-precision SIMD variable
 * \return Single-precision SIMD variable of the same width
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_d2f(gmx_simd_double_t d)
{
    gmx_simd_float_t f;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i] = d.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_d2f() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    f.r[0] = d.r[0];
#endif
    return f;
}

/*! \brief Convert SIMD float to double.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD variable
 * \param[out] d0 Double-precision SIMD variable, first half of values from f.
 * \param[out] d1 Double-precision SIMD variable, second half of values from f.
 */
static gmx_inline void
gmx_simd_cvt_f2dd(gmx_simd_float_t f, gmx_simd_double_t *d0, gmx_simd_double_t *d1)
{
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d0->r[i] = f.r[i];
        d1->r[i] = f.r[GMX_SIMD_DOUBLE_WIDTH + i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_f2dd() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments */
    d0->r[0] = f.r[0];
    d1->r[0] = f.r[0];
#endif
}

/*! \brief Convert SIMD double to float.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH.
 *
 * Float/double conversions are complex since the SIMD width could either
 * be different (e.g. on x86) or identical (e.g. IBM QPX). This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d0 Double-precision SIMD variable, first half of values to put in f.
 * \param d1 Double-precision SIMD variable, second half of values to put in f.
 * \return Single-precision SIMD variable with all values.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_cvt_dd2f(gmx_simd_double_t d0, gmx_simd_double_t d1)
{
    gmx_simd_float_t f;
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.r[i]                         = d0.r[i];
        f.r[GMX_SIMD_DOUBLE_WIDTH + i] = d1.r[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_dd2f() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments & uninitialized f */
    f.r[0] = d0.r[0] + d1.r[0];
#endif
    return f;
}

/*! \brief Convert SIMD float bool to double bool.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH.  Since this version returns one gmx_simd_dbool_t
 * we use a single "db" in the function name (repeating it indicates two).
 *
 * The float/double boolean conversions have similar challenges with the simd
 * width as the floating-point versions. This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD boolean variable
 * \return Double-precision SIMD boolean variable of the same width
 */
static gmx_inline gmx_simd_dbool_t
gmx_simd_cvt_fb2db(gmx_simd_fbool_t f)
{
    gmx_simd_dbool_t  d;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int               i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d.b[i] = f.b[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_fb2db() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    d.b[0] = f.b[0];
#endif
    return d;
}

/*! \brief Convert SIMD double bool to float bool.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is identical to
 * \ref GMX_SIMD_DOUBLE_WIDTH. Since this version only uses a gmx_simd_dbool_t
 * we use a single "db" in the function name (repeating it indicates two).
 *
 * The float/double boolean conversions have similar challenges with the simd
 * width as the floating-point versions. This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d Double-precision SIMD boolean variable
 * \return Single-precision SIMD boolean variable of the same width
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cvt_db2fb(gmx_simd_dbool_t d)
{
    gmx_simd_fbool_t f;
#if (GMX_SIMD_FLOAT_WIDTH == GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.b[i] = d.b[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_db2fb() requires GMX_SIMD_FLOAT_WIDTH==GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings */
    f.b[0] = d.b[0];
#endif
    return f;
}

/*! \brief Convert SIMD float bool to double bool.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH. Since this version uses two gmx_simd_dbool_t
 * arguments, we use "dbdb" in the function name.
 *
 * The float/double boolean conversions have similar challenges with the simd
 * width as the floating-point versions. This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param f Single-precision SIMD boolean variable
 * \param[out] d0 Double-precision SIMD boolean, first half of values from f.
 * \param[out] d1 Double-precision SIMD boolean, second half of values from f.
 */
static gmx_inline void
gmx_simd_cvt_fb2dbdb(gmx_simd_fbool_t f, gmx_simd_dbool_t *d0, gmx_simd_dbool_t *d1)
{
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        d0->b[i] = f.b[i];
        d1->b[i] = f.b[GMX_SIMD_DOUBLE_WIDTH + i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_fb2dbdb() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments */
    d0->b[0] = f.b[0];
    d1->b[0] = f.b[0];
#endif
}

/*! \brief Convert SIMD double bool to float bool.
 *
 * This version is available if \ref GMX_SIMD_FLOAT_WIDTH is twice as large
 * as \ref GMX_SIMD_DOUBLE_WIDTH. Since this version uses two gmx_simd_dbool_t
 * arguments, we use "dbdb" in the function name.
 *
 * The float/double boolean conversions have similar challenges with the simd
 * width as the floating-point versions. This means you will
 * need to check for the width in the code, and have different code paths.
 *
 * \param d0 Double-precision SIMD boolean, first half of values to put in f.
 * \param d1 Double-precision SIMD boolean, second half of values to put in f.
 * \return Single-precision SIMD boolean with all values.
 */
static gmx_inline gmx_simd_fbool_t
gmx_simd_cvt_dbdb2fb(gmx_simd_dbool_t d0, gmx_simd_dbool_t d1)
{
    gmx_simd_fbool_t f;
#if (GMX_SIMD_FLOAT_WIDTH == 2*GMX_SIMD_DOUBLE_WIDTH)
    int              i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        f.b[i]                       = d0.b[i];
        f.b[GMX_SIMD_DOUBLE_WIDTH+i] = d1.b[i];
    }
#else
    gmx_fatal(FARGS, "gmx_simd_cvt_dbdb2fb() requires GMX_SIMD_FLOAT_WIDTH==2*GMX_SIMD_DOUBLE_WIDTH");
    /* Avoid compiler warnings about unused arguments & uninitialized f */
    f.b[0] = d0.b[0] + d1.b[0];
#endif
    return f;
}

/*! \}
 *
 * \name Higher-level SIMD utility functions, single precision.
 *
 * These include generic functions to work with triplets of data, typically
 * coordinates, and a few utility functions to load and update data in the
 * nonbonded kernels.
 * These functions should be available on all implementations, although
 * some wide SIMD implementations (width>=8) also provide special optional
 * versions to work with half or quarter registers to improve the performance
 * in the nonbonded kernels.
 * \{
 */

/*! \brief Load float triplets from non-adjacent unaligned memory and transpose.
 *
 * \param      base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param[out] v0     First component of each triplet, base[offset[i]] for each i.
 * \param[out] v1     Second component of each triplet, base[offset[i] + 1] for each i.
 * \param[out] v2     Third component of each triplet, base[offset[i] + 2] for each i.
 *
 * The offset array should have the same length as the SIMD width.
 * The memory locations do not have to be aligned, and the routine will not
 * try to load other memory before or beyond the data (no padding required).
 *
 * \note You should NOT scale offsets by 3 before calling this routine; it is
 *       done internally, where some architectures can do it more efficiently.
 * \note On most architectures this will be implemented with slow multiple loads
 *       and shuffling operations, so try to work with simd-layout data instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
static gmx_inline void
gmx_simd_gather_load_3_transpose_f(const float * base, const gmx_int32_t offset[],
                                   gmx_simd_float_t * v0, gmx_simd_float_t * v1, gmx_simd_float_t * v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[3 * offset[i]];
        v1->r[i] = base[3 * offset[i] + 1];
        v2->r[i] = base[3 * offset[i] + 2];
    }
}

/*! \brief Transpose and store float triplets to non-adjacent unaligned memory.
 *
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     First component of triplets, written to base[offset[i]].
 * \param      v1     Second component of triplets, written to base[offset[i] + 1].
 * \param      v2     Third component of triplets, written to base[offset[i] + 2].
 *
 * The offset array should have the same length as the SIMD width.
 * The memory locations do not have to be aligned, and the routine will not
 * store to other memory before or beyond the data (no padding required).
 *
 * \note You should NOT scale offsets by 3 before calling this routine; it is
 *       done internally, where some architectures can do it more efficiently.
 * \note On most architectures this will be implemented with slow multiple loads
 *       and shuffling operations, so try to work with simd-layout data instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
static gmx_inline void
gmx_simd_transpose_scatter_store_3_f(float * base, const gmx_int32_t offset[],
                                     gmx_simd_float_t v0, gmx_simd_float_t v1, gmx_simd_float_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[3 * offset[i]]     = v0.r[i];
        base[3 * offset[i] + 1] = v1.r[i];
        base[3 * offset[i] + 2] = v2.r[i];
    }
}

/*! \brief Transpose and add to non-adjacent unaligned memory float triplets.
 *
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     First component of triplets, added to base[offset[i]].
 * \param      v1     Second component of triplets, added to base[offset[i] + 1].
 * \param      v2     Third component of triplets, added to base[offset[i] + 2].
 *
 * The offset array should have the same length as the SIMD width.
 * The memory locations do not have to be aligned, and the routine will not
 * store to other memory before or beyond the data (no padding required).
 *
 * \note You should NOT scale offsets by 3 before calling this routine; it is
 *       done internally, where some architectures can do it more efficiently.
 * \note On most architectures this will be implemented with slow multiple loads
 *       and shuffling operations, so try to work with simd-layout data instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
static gmx_inline void
gmx_simd_transpose_scatter_incr_3_f(float * base, const gmx_int32_t offset[],
                                    gmx_simd_float_t v0, gmx_simd_float_t v1, gmx_simd_float_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[3 * offset[i]]     += v0.r[i];
        base[3 * offset[i] + 1] += v1.r[i];
        base[3 * offset[i] + 2] += v2.r[i];
    }
}

/*! \brief Transpose and subtract from non-adjacent unaligned memory float triplets.
 *
 * \param[out] base   Pointer to the start of the memory area
 * \param      offset Aligned array with offsets to the start of each triplet.
 * \param      v0     First component of triplets, subtracted from base[offset[i]].
 * \param      v1     Second component of triplets, subtracted from base[offset[i] + 1].
 * \param      v2     Third component of triplets, subtracted from base[offset[i] + 2].
 *
 * The offset array should have the same length as the SIMD width.
 * The memory locations do not have to be aligned, and the routine will not
 * store to other memory before or beyond the data (no padding required).
 *
 * \note You should NOT scale offsets by 3 before calling this routine; it is
 *       done internally, where some architectures can do it more efficiently.
 * \note On most architectures this will be implemented with slow multiple loads
 *       and shuffling operations, so try to work with simd-layout data instead.
 * \note This routine uses a normal array for the offsets, since we typically
 *       load the data from memory. On the architectures we have tested this
 *       is faster even when a SIMD integer datatype is present.
 */
static gmx_inline void
gmx_simd_transpose_scatter_decr_3_f(float * base, const gmx_int32_t offset[],
                                    gmx_simd_float_t v0, gmx_simd_float_t v1, gmx_simd_float_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        base[3 * offset[i]]     -= v0.r[i];
        base[3 * offset[i] + 1] -= v1.r[i];
        base[3 * offset[i] + 2] -= v2.r[i];
    }
}

/*! \brief Expand floats into consecutive triplets in three outputs.
 *
 * \param      scalar    Floating-point input, e.g. [s0 s1 s2 s3] if width=4.
 * \param[out] triplets0 First output, e.g. [s0 s0 s0 s1] if width=4.
 * \param[out] triplets1 Second output, e.g. [s1 s1 s2 s2] if width=4.
 * \param[out] triplets2 Third output, e.g. [s2 s3 s3 s3] if width=4.
 *
 * This routine is meant to use for things like scalar-vector multiplication,
 * where the vectors are stored in a merged format like [x0 y0 z0 x1 y1 z1 ...],
 * while the scalars are stored as [s0 s1 s2...], and the data cannot easily
 * be changed to SIMD-friendly layout.
 *
 * In this case, load 3 full-width SIMD variables from the vector array (This
 * will always correspond to GMX_SIMD_FLOAT_WIDTH/GMX_SIMD_DOUBLE_WIDTH
 * triplets), load a single full-width variable from the scalar array, and
 * call this routine to expand the data. You can then simply multiply the
 * first, second and third pair of SIMD variables, and store the three
 * results back into a suitable vector-format array.
 */
static gmx_inline void
gmx_simd_expand_scalars_to_triplets_f(gmx_simd_float_t   scalar,
                                      gmx_simd_float_t * triplets0,
                                      gmx_simd_float_t * triplets1,
                                      gmx_simd_float_t * triplets2)
{
    int i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        triplets0->r[i] = scalar.r[i / 3];
        triplets1->r[i] = scalar.r[(i + GMX_SIMD_FLOAT_WIDTH) / 3];
        triplets2->r[i] = scalar.r[(i + 2 * GMX_SIMD_FLOAT_WIDTH) / 3];
    }
}

/*! \brief Alignment to use for specially aligned float pair data.
 *
 *  The nonbonded kernels depend on being able to load a pair of parameters
 *  (for Lennard-Jones interactions) quickly. Some SIMD implementations provide
 *  support for doing this efficiently, even from unaligned memory, while others
 *  might require padding, which wastes cache. To allow each architecture to
 *  use the most optimal form, we use a constant that code outside the SIMD
 *  module should use to store things properly. It must be at least 2. For
 *  example, a value of 2 means the two parameters A & B are stored as
 *  [A0 B0 A1 B1] while stride-4 means [A0 B0 - - A1 B1 - -].
 *
 *  This alignment depends on the efficiency of partial-register load/store
 *  operations, and will depend on the architecture.
 */
static const int gmx_simd_pairs_storage_alignment_f = 4;


/*! \brief Load pairs that are custom-aligned in memory and transpose.
 *
 * \param      base   Pointer to start of aligned memory
 * \param      offset Offset to the start of each pair
 * \param[out] v0     First element in each pair
 * \param[out] v1     Second element in each pair
 *
 * The offset array should have the same length as the SIMD width.
 *
 * \note All pairs should have alignment specified by
 * \ref gmx_simd_pairs_storage_alignment_f (see that description for
 * an explanation of the format). This also means the offset
 * is scaled by this number internally, so you should NOT scale any offsets
 * it manually too.
 */
static gmx_inline void
gmx_simd_load_pairs_transpose_f(const float * base,    gmx_int32_t offset[],
                                gmx_simd_float_t * v0, gmx_simd_float_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[gmx_simd_pairs_storage_alignment_f * offset[i]];
        v1->r[i] = base[gmx_simd_pairs_storage_alignment_f * offset[i] + 1];
    }
}


/*! \brief Load two unaligned float elements and transpose.
 *
 * Only available when \ref GMX_SIMD_HAVE_LOADU_2_TRANSPOSE_FLOAT is defined for
 * single, or \ref GMX_SIMD_HAVE_LOADU_2_TRANSPOSE_DOUBLE for double.
 *
 * \param      base   Pointer to start of pair data
 * \param      offset Offset to the start of each pair
 * \param[out] v0     First element in each pair, base[offset[i]] for each i.
 * \param[out] v1     Second element in each pair, base[offset[i] + 1] for each i.
 *
 * The offset array should have the same length as the SIMD width.
 *
 * \note This routine is intended for reading (linear) table data, which
 * means the offset is not scaled - the two values in memory do not form a single
 * pair, but represent adjacent (packed) table points. This is also the
 * reason why the data cannot be aligned.
 */
static gmx_inline void
gmx_simd_loadu_2_transpose_f(const float * base,    gmx_int32_t offset[],
                             gmx_simd_float_t * v0, gmx_simd_float_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[offset[i]];
        v1->r[i] = base[offset[i] + 1];
    }
}


/*! \brief Load float quartets from non-adjacent memory aligned-to-4 and transpose.
 *
 * \param      base   Pointer to the start of the memory area
 * \param      offset SIMD integer type with offsets to the start of each triplet.
 * \param[out] v0     First component, base[4 * offset[i]] for each i.
 * \param[out] v1     Second component, base[4 * offset[i] + 1] for each i.
 * \param[out] v2     Third component, base[4 * offset[i] + 2] for each i.
 * \param[out] v3     Fourth component, base[4 * offset[i] + 3] for each i.
 *
 * The memory locations must be aligned, but only to four elements even if the
 * SIMD implementation is larger width.
 *
 * \note that you should NOT scale the offsets by 4 before calling this routine;
 * it is done internally, where some architectures can do it more efficiently.
 *
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible.
 */
static gmx_inline void
gmx_simd_load_4_transpose_f(const float * base,    gmx_simd_fint32_t offset,
                            gmx_simd_float_t * v0, gmx_simd_float_t * v1,
                            gmx_simd_float_t * v2, gmx_simd_float_t * v3)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[4 * offset.i[i]];
        v1->r[i] = base[4 * offset.i[i] + 1];
        v2->r[i] = base[4 * offset.i[i] + 2];
        v3->r[i] = base[4 * offset.i[i] + 3];
    }
}

/*! \brief Load two floats from non-adjacent memory aligned-to-4 and transpose.
 *
 * \param      base   Pointer to the start of the memory area
 * \param      offset Offsets to the start of each triplet.
 * \param[out] v0     First component, base[4 * offset[i]] for each i.
 * \param[out] v1     Second component, base[4 * offset[i] + 1] for each i.
 *
 * The memory locations must be aligned, but only to four elements even if the
 * SIMD implementation is larger width.
 *
 * \note that you should NOT scale the offsets by 4 before calling this routine;
 * it is done internally, where some architectures can do it more efficiently.
 *
 * \note This is a special routine primarily intended for loading Gromacs
 *       table data as efficiently as possible.
 */
static gmx_inline void
gmx_simd_load_4_transpose_f(const float * base,    gmx_simd_fint32_t offset,
                            gmx_simd_float_t * v0, gmx_simd_float_t * v1,
                            gmx_simd_float_t * v2, gmx_simd_float_t * v3)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH; i++)
    {
        v0->r[i] = base[4 * offset.i[i]];
        v1->r[i] = base[4 * offset.i[i] + 1];
        v2->r[i] = base[4 * offset.i[i] + 2];
        v3->r[i] = base[4 * offset.i[i] + 3];
    }
}

/*! \brief Reduce four SIMD floats, increment adjacent floats in mem, return sum.
 *
 * \param m   Pointer to memory where four floats should be incremented
 * \param v0  SIMD variable whose sum should be added to m[0]
 * \param v1  SIMD variable whose sum should be added to m[1]
 * \param v2  SIMD variable whose sum should be added to m[2]
 * \param v3  SIMD variable whose sum should be added to m[3]
 *
 * \return Sum of all elements in the four SIMD variables.
 *
 * \note This is a special routine intended for the Gromacs nonbonded kernels.
 * It is used in the epilogue of the outer loop, where the variables will
 * contain unrolled forces for one outer-loop-particle each, corresponding to
 * a single coordinate (i.e, say, four x-coordinate force variables). These
 * should be summed and added to the force array in memory. Since we always work
 * with contiguous SIMD-layout , we can use efficient aligned loads/stores.
 * When calculating the virial, we also need the total sum of all forces for
 * each coordinate. This is provided as the return value. For routines that
 * do not need these, this extra code will be optimized away completely if you
 * just ignore the return value (Checked with gcc-4.9.1 and clang-3.6 for AVX).
 */
static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_f(float * m,
                                    gmx_simd_float_t v0, gmx_simd_float_t v1,
                                    gmx_simd_float_t v2, gmx_simd_float_t v3)
{
    /* Note that the 4 here corresponds to the 4 elements, not any SIMD width */
    float sum[4];

    sum[0] = gmx_simd_reduce_f(v0);
    sum[1] = gmx_simd_reduce_f(v1);
    sum[2] = gmx_simd_reduce_f(v2);
    sum[3] = gmx_simd_reduce_f(v3);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}


/*! \}
 *
 * \name Higher-level SIMD utility functions, double precision.
 *
 * These include generic functions to work with triplets of data, typically
 * coordinates, and a few utility functions to load and update data in the
 * nonbonded kernels.
 * These functions should be available on all implementations.
 * \{
 */

/*! \brief Load double triplets from non-adjacent unaligned memory and transpose.
 *
 * \copydetails gmx_simd_gather_load_3_transpose_f
 */
static gmx_inline void
gmx_simd_gather_load_3_transpose_d(const double * base, const gmx_int32_t offset[],
                                   gmx_simd_double_t * v0, gmx_simd_double_t * v1, gmx_simd_double_t * v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0->r[i] = base[3 * offset[i]];
        v1->r[i] = base[3 * offset[i] + 1];
        v2->r[i] = base[3 * offset[i] + 2];
    }
}

/*! \brief Transpose and store double triplets to non-adjacent unaligned memory.
 *
 * \copydetails gmx_simd_transpose_scatter_store_3_f
 */
static gmx_inline void
gmx_simd_transpose_scatter_store_3_d(double * base, const gmx_int32_t offset[],
                                     gmx_simd_double_t v0, gmx_simd_double_t v1, gmx_simd_double_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[3 * offset[i]]     = v0.r[i];
        base[3 * offset[i] + 1] = v1.r[i];
        base[3 * offset[i] + 2] = v2.r[i];
    }
}

/*! \brief Transpose and add to non-adjacent unaligned memory double triplets.
 *
 * \copydetails gmx_simd_transpose_scatter_incr_3_f
 */
static gmx_inline void
gmx_simd_transpose_scatter_incr_3_d(double * base, const gmx_int32_t offset[],
                                    gmx_simd_double_t v0, gmx_simd_double_t v1, gmx_simd_double_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[3 * offset[i]]     += v0.r[i];
        base[3 * offset[i] + 1] += v1.r[i];
        base[3 * offset[i] + 2] += v2.r[i];
    }
}

/*! \brief Transpose and subtract from non-adjacent unaligned memory double triplets.
 *
 * \\copydetails gmx_simd_transpose_scatter_decr_3_f
 */
static gmx_inline void
gmx_simd_transpose_scatter_decr_3_d(double * base, const gmx_int32_t offset[],
                                    gmx_simd_double_t v0, gmx_simd_double_t v1, gmx_simd_double_t v2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        base[3 * offset[i]]     -= v0.r[i];
        base[3 * offset[i] + 1] -= v1.r[i];
        base[3 * offset[i] + 2] -= v2.r[i];
    }
}

/*! \brief Expand floats into consecutive triplets in three outputs.
 *
 * \copydetails gmx_simd_expand_scalars_to_triplets_f
 */
static gmx_inline void
gmx_simd_expand_scalars_to_triplets_d(gmx_simd_double_t   scalar,
                                      gmx_simd_double_t * triplets0,
                                      gmx_simd_double_t * triplets1,
                                      gmx_simd_double_t * triplets2)
{
    int i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        triplets0->r[i] = scalar.r[i / 3];
        triplets1->r[i] = scalar.r[(i + GMX_SIMD_DOUBLE_WIDTH) / 3];
        triplets2->r[i] = scalar.r[(i + 2 * GMX_SIMD_DOUBLE_WIDTH) / 3];
    }
}

/*! \brief Alignment to use for specially aligned float pair data.
 *
 *  \copydetails gmx_simd_pairs_storage_alignment_f
 */
static const int gmx_simd_pairs_storage_alignment_d = 2;


/*! \brief Load pairs that are custom-aligned in memory and transpose.
 *
 * \copydetails gmx_simd_load_pairs_transpose_f
 */
static gmx_inline void
gmx_simd_load_pairs_transpose_d(const double * base,    gmx_int32_t offset[],
                                gmx_simd_double_t * v0, gmx_simd_double_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0->r[i] = base[gmx_simd_pairs_storage_alignment_d * offset[i]];
        v1->r[i] = base[gmx_simd_pairs_storage_alignment_d * offset[i] + 1];
    }
}

/*! \brief Load two unaligned double elements and transpose.
 *
 * \copydetails gmx_simd_loadu_2_transpose_f
 */
static gmx_inline void
gmx_simd_loadu_2_transpose_d(const double * base,    gmx_int32_t offset[],
                             gmx_simd_double_t * v0, gmx_simd_double_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0->r[i] = base[offset[i]];
        v1->r[i] = base[offset[i] + 1];
    }
}

/*! \brief Load float quartets from non-adjacent memory aligned-to-4 and transpose.
 *
 * \copydetails gmx_simd_load_4_transpose_f
 */
static gmx_inline void
gmx_simd_load_4_transpose_d(const double * base,    gmx_simd_dint32_t offset,
                            gmx_simd_double_t * v0, gmx_simd_double_t * v1,
                            gmx_simd_double_t * v2, gmx_simd_double_t * v3)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0->r[i] = base[4 * offset.i[i]];
        v1->r[i] = base[4 * offset.i[i] + 1];
        v2->r[i] = base[4 * offset.i[i] + 2];
        v3->r[i] = base[4 * offset.i[i] + 3];
    }
}

/*! \brief Load two floats from non-adjacent memory aligned-to-4 and transpose.
 *
 * \copydetails gmx_simd_load_2of4_transpose_f
 */
static gmx_inline void
gmx_simd_load_2of4_transpose_d(const double * base,    gmx_simd_dint32_t offset,
                               gmx_simd_double_t * v0, gmx_simd_double_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH; i++)
    {
        v0->r[i] = base[4 * offset.i[i]];
        v1->r[i] = base[4 * offset.i[i] + 1];
    }
}

/*! \brief Reduce four SIMD floats, increment adjacent floats in mem, return sum.
 *
 * \copydetails gmx_simd_reduce_incr_4_return_sum_f
 */
static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_d(double * m,
                                    gmx_simd_double_t v0, gmx_simd_double_t v1,
                                    gmx_simd_double_t v2, gmx_simd_double_t v3)
{
    /* The 4 here corresponds to the 4 elements, not any SIMD width */
    double sum[4];

    sum[0] = gmx_simd_reduce_d(v0);
    sum[1] = gmx_simd_reduce_d(v1);
    sum[2] = gmx_simd_reduce_d(v2);
    sum[3] = gmx_simd_reduce_d(v3);

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}

/*! \}
 *
 * \name Higher-level SIMD utilities accessing partial (half-width) SIMD floats.
 *
 * These functions are optional. The are only useful for SIMD implementation
 * where the width is 8 or larger, and where it would be inefficient
 * to process 4*8, 8*8, or more, interactions in parallel.
 *
 * Currently, only Intel provides very wide SIMD implementations, but these
 * also come with excellent support for loading, storing, accessing and
 * shuffling parts of the register in so-called 'lanes' of 4 bytes each.
 * We can use this to load separate parts into the low/high halves of the
 * register in the inner loop of the nonbonded kernel, which e.g. makes it
 * possible to process 4*4 nonbonded interactions as a pattern of 2*8. We
 * can also use implementations with width 16 or greater.
 *
 * To make this more generic, when \ref GMX_SIMD_HAVE_HSIMD_UTIL is
 * defined, the SIMD implementation provides seven special routines that:
 *
 * - Load the low/high parts of a SIMD variable from different pointers
 * - Load half the SIMD width from one pointer, and duplicate in low/high parts
 * - Load two reals, put 1st one in all low elements, and 2nd in all high ones.
 * - Store the low/high parts of a SIMD variable to different pointers
 * - Subtract both SIMD halves from a single half-SIMD-width memory location.
 * - Load aligned pairs (LJ parameters) from two base pointers, with a common
 *   offset list, and put these in the low/high SIMD halves.
 * - Reduce each half of two SIMD registers (i.e., 4 parts in total), increment
 *   four adjacent memory positions, and return the total sum.
 *
 * But remember, this is ONLY used when the native SIMD width is large. You will
 * just waste time if you implement it for normal 16-byte SIMD architectures.
 * \{
 */

/*! \brief Load low & high parts of SIMD float from different locations
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 *
 * \return SIMD variable with low part loaded from m0, high from m1.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_load_dual_hsimd_f(const float * m0, const float * m1)
{
    gmx_simd_float_t a;
    int              i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m0[i];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = m1[i];
    }
    return a;
}

/*! \brief Load half-SIMD-width float data, spread to both halves.
 *
 * \param m Pointer to memory aligned to half SIMD width.
 *
 * \return SIMD variable with both halves loaded from m..
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_loaddup_hsimd_f(const float * m)
{
    gmx_simd_float_t a;
    int              i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m[i];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = a.r[i];
    }
    return a;
}

/*! \brief Load two floats, spread 1st in low half, 2nd in high half.
 *
 * \param m Pointer to two adjacent floating-point values.
 *
 * \return SIMD variable where all elements in the low half have been set
 *         to m[0], and all elements in high half to m[1].
 *
 * \note This routine always loads two values and sets the halves separately.
 *       If you want to set all elements to the same value, simply use
 *       the standard \ref gmx_simd_load1_r.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline gmx_simd_float_t
gmx_simd_load1_dual_hsimd_f(const float * m)
{
    gmx_simd_float_t a;
    int              i;

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        a.r[i]                          = m[0];
        a.r[GMX_SIMD_FLOAT_WIDTH/2 + i] = m[1];
    }
    return a;
}


/*! \brief Load low & high parts of SIMD float from different locations
 *
 * \param m0 Pointer to memory aligned to half SIMD width.
 * \param m1 Pointer to memory aligned to half SIMD width.
 * \param a  SIMD variable. Low half should be stored to m0, high to m1.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline void
gmx_simd_store_dual_hsimd_f(float * m0, float * m1, gmx_simd_float_t a)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        m0[i] = a.r[i];
        m1[i] = a.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }
}

/*! \brief Subtract both halves from half-SIMD-width memory float data.
 *
 * \param m  half-width aligned memory, from which sum of the halves will be subtracted.
 * \param a  SIMD variable. Upper & lower halves will first be added.
 *
 * If the SIMD width is 8 and contains [a b c d e f g h], the
 * memory will be modified to [m[0]-(a+e) m[1]-(b+f) m[2]-(c+g) m[3]-(d+h)].
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline void
gmx_simd_decr_hsimd_f(float * m, gmx_simd_float_t a)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        m[i] -= a.r[i] + a.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }
}


/*! \brief Load pairs custom-aligned from two memory bases, transpose to simd
 *
 * \param      base0  Pointer to base of first aligned memory
 * \param      base1  Pointer to base of second aligned memory
 * \param      offset Offset to the start of each pair
 * \param[out] v0     1st element in each pair, base0 in low and base1 in high half.
 * \param[out] v1     2nd element in each pair, base0 in low and base1 in high half.
 *
 * The offset array should be of half the SIMD width length, so it corresponds
 * to the half-SIMD-registers.
 *
 * This routine is primarily designed to load nonbonded parameters in the
 * kernels. It is the equivalent of the full-width routine
 * \ref gmx_simd_load_pairs_transpose_r, but just
 * as the other hsimd routines it will pick half-SIMD-width data from base0
 * and put in the lower half, while the upper half comes from base1.
 *
 * For an example, assume the SIMD width is 8, pair alignment is 2, that
 * base0 is [A0 A1 B0 B1 C0 C1 D0 D1 ...], and base1 [E0 E1 F0 F1 G0 G1 H0 H1...].
 *
 * Then we will get v0 as [A0 B0 C0 D0 E0 F0 G0 H0] and v1 as [A1 B1 C1 D1 E1 F1 G1 H1].
 *
 * \note All pairs should have alignment specified by
 * \ref gmx_simd_pairs_storage_alignment_f (see that description for
 * an explanation of the format). This also means the offset
 * is scaled by this number internally, so you should NOT scale the offsets
 * manually too.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline void
gmx_simd_load_pairs_transpose_hsimd_f(const float * base0, const float * base1,
                                      gmx_int32_t offset[],
                                      gmx_simd_float_t * v0, gmx_simd_float_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        v0->r[i] = base0[gmx_simd_pairs_storage_alignment_f * offset[i]];
        v1->r[i] = base0[gmx_simd_pairs_storage_alignment_f * offset[i] + 1];
        v0->r[GMX_SIMD_FLOAT_WIDTH/2 + i] = base1[gmx_simd_pairs_storage_alignment_f * offset[i]];
        v1->r[GMX_SIMD_FLOAT_WIDTH/2 + i] = base1[gmx_simd_pairs_storage_alignment_f * offset[i] + 1];
    }
}


/*! \brief Reduce the 4 half-SIMD floats in 2 variables, increment mem, return sum.
 *
 * \param m    Pointer to memory where the four values should be incremented
 * \param v0   Variable whose half-SIMD sums should be added to m[0]/m[1], respectively.
 * \param v1   Variable whose half-SIMD sums should be added to m[2]/m[3], respectively.
 *
 * \return Sum of all elements in the four SIMD variables.
 *
 * \note This is the half-SIMD-width version of
 * \ref gmx_simd_reduce_incr_4_return_sum_r. The only difference is that the
 *      four half-SIMD inputs needed are present in the low/high halves of the
 *      two SIMD arguments.
 *
 * Available when \ref GMX_SIMD_HAVE_HSIMD_FLOAT_UTIL is defined for single,
 * or \ref GMX_SIMD_HAVE_HSIMD_DOUBLE_UTIL for double.
 */
static gmx_inline float
gmx_simd_reduce_incr_4_return_sum_hsimd_f(float * m, gmx_simd_float_t v0, gmx_simd_float_t v1)
{
    /* The 4 here corresponds to the 4 elements, not any SIMD width */
    float sum[4];
    int   i;

    for (i = 0; i < 4; i++)
    {
        sum[i] = 0;
    }

    for (i = 0; i < GMX_SIMD_FLOAT_WIDTH/2; i++)
    {
        sum[0] += v0.r[i];
        sum[1] += v0.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
        sum[2] += v1.r[i];
        sum[3] += v1.r[GMX_SIMD_FLOAT_WIDTH/2 + i];
    }

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}


/*! \}
 *
 * \name Higher-level SIMD utilities accessing partial (half-width) SIMD doubles.
 *
 * See the single-precision versions for documentation. Since double precision
 * is typically half the width of single, this double version is likely only
 * useful with 512-bit and larger implementations.
 *
 * \{
 */

/*! \brief Load low & high parts of SIMD float from different locations
 *
 * \copydetails gmx_simd_load_dual_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load_dual_hsimd_d(const double * m0, const double * m1)
{
    gmx_simd_double_t a;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m0[i];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = m1[i];
    }
    return a;
}

/*! \brief Load half-SIMD-width float data, spread to both halves.
 *
 * \copydetails gmx_simd_loaddup_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_loaddup_hsimd_d(const double * m)
{
    gmx_simd_double_t a;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m[i];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = a.r[i];
    }
    return a;
}

/*! \brief Load two floats, spread 1st in low half, 2nd in high half.
 *
 * \copydetails gmx_simd_load1_dual_hsimd_f
 */
static gmx_inline gmx_simd_double_t
gmx_simd_load1_dual_hsimd_d(const double * m)
{
    gmx_simd_double_t a;
    int               i;

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        a.r[i]                           = m[0];
        a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = m[1];
    }
    return a;
}


/*! \brief Load low & high parts of SIMD float from different locations
 *
 * \copydetails gmx_simd_store_dual_hsimd_f
 */
static gmx_inline void
gmx_simd_store_dual_hsimd_d(double * m0, double * m1, gmx_simd_double_t a)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        m0[i] = a.r[i];
        m1[i] = a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }
}

/*! \brief Subtract both halves from half-SIMD-width memory float data.
 *
 * \copydetails gmx_simd_decr_hsimd_f
 */
static gmx_inline void
gmx_simd_decr_hsimd_d(double * m, gmx_simd_double_t a)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        m[i] -= a.r[i] + a.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }
}


/*! \brief Load pairs custom-aligned from two memory bases, transpose to simd
 *
 * \copydetails gmx_simd_load_pairs_transpose_hsimd_f
 */
static gmx_inline void
gmx_simd_load_pairs_transpose_hsimd_d(const double * base0, const double * base1,
                                      gmx_int32_t offset[],
                                      gmx_simd_double_t * v0, gmx_simd_double_t * v1)
{
    int i;
    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        v0->r[i] = base0[gmx_simd_pairs_storage_alignment_d * offset[i]];
        v1->r[i] = base0[gmx_simd_pairs_storage_alignment_d * offset[i] + 1];
        v0->r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = base1[gmx_simd_pairs_storage_alignment_d * offset[i]];
        v1->r[GMX_SIMD_DOUBLE_WIDTH/2 + i] = base1[gmx_simd_pairs_storage_alignment_d * offset[i] + 1];
    }
}


/*! \brief Reduce the 4 half-SIMD floats in 2 variables, increment mem, return sum.
 *
 * \copydetails gmx_simd_reduce_incr_4_return_sum_hsimd_f
 */
static gmx_inline double
gmx_simd_reduce_incr_4_return_sum_hsimd_d(double * m, gmx_simd_double_t v0, gmx_simd_double_t v1)
{
    /* The 4 here corresponds to the 4 elements, not any SIMD width */
    double sum[4];
    int    i;

    for (i = 0; i < 4; i++)
    {
        sum[i] = 0;
    }

    for (i = 0; i < GMX_SIMD_DOUBLE_WIDTH/2; i++)
    {
        sum[0] += v0.r[i];
        sum[1] += v0.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
        sum[2] += v1.r[i];
        sum[3] += v1.r[GMX_SIMD_DOUBLE_WIDTH/2 + i];
    }

    m[0] += sum[0];
    m[1] += sum[1];
    m[2] += sum[2];
    m[3] += sum[3];

    return sum[0] + sum[1] + sum[2] + sum[3];
}





/*! \} */

/*! \name SIMD4. Constant width-4 SIMD types and instructions
 * \{
 */

#if (GMX_SIMD_FLOAT_WIDTH == 4) || (defined DOXYGEN)


/*! \brief SIMD4 float type. Available with \ref GMX_SIMD4_HAVE_FLOAT.
 *
 * Unless you specifically want a single-precision type you should check
 * \ref gmx_simd4_real_t instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define gmx_simd4_float_t    gmx_simd_float_t

/*! \brief Load SIMD4 float from aligned memory.
 *  \copydetails gmx_simd_load_f
 */
#    define gmx_simd4_load_f     gmx_simd_load_f

/*! \brief Set all elements of SIMD4 float from single pointer.
 *  \copydetails gmx_simd_load1_f
 */
#    define gmx_simd4_load1_f    gmx_simd_load1_f

/*! \brief Set all SIMD4 float elements to the value r.
 *  \copydetails gmx_simd_set1_f
 */
#    define gmx_simd4_set1_f     gmx_simd_set1_f

/*! \brief Store the contents of SIMD4 float pr to aligned memory m.
 *  \copydetails gmx_simd_store_f
 */
#    define gmx_simd4_store_f    gmx_simd_store_f

/*! \brief Load SIMD4 float from unaligned memory.
 * \copydetails gmx_simd_loadu_f
 */
#    define gmx_simd4_loadu_f    gmx_simd_loadu_f

/*! \brief Store SIMD4 float to unaligned memory.
 * \copydetails gmx_simd_storeu_f
 */
#    define gmx_simd4_storeu_f   gmx_simd_storeu_f

/*! \brief Set all SIMD4 float elements to 0.
 * \copydetails gmx_simd_setzero_f
 */
#    define gmx_simd4_setzero_f  gmx_simd_setzero_f

/*! \brief Bitwise and for two SIMD4 float variables.
 * \copydetails gmx_simd_and_f
 */
#    define gmx_simd4_and_f      gmx_simd_and_f

/*! \brief Bitwise andnot for two SIMD4 float variables. c=(~a) & b.
 * \copydetails gmx_simd_andnot_f
 */
#    define gmx_simd4_andnot_f   gmx_simd_andnot_f

/*! \brief Bitwise or for two SIMD4 float variables.
 * \copydetails gmx_simd_or_f
 */
#    define gmx_simd4_or_f       gmx_simd_or_f

/*! \brief Bitwise xor for two SIMD4 float variables.
 * \copydetails gmx_simd_xor_f
 */
#    define gmx_simd4_xor_f      gmx_simd_xor_f

/*! \brief Add two SIMD4 float variables.
 * \copydetails gmx_simd_add_f
 */
#    define gmx_simd4_add_f      gmx_simd_add_f

/*! \brief Subtract two SIMD4 float variables.
 * \copydetails gmx_simd_sub_f
 */
#    define gmx_simd4_sub_f      gmx_simd_sub_f

/*! \brief Multiply two SIMD4 float variables.
 * \copydetails gmx_simd_mul_f
 */
#    define gmx_simd4_mul_f      gmx_simd_mul_f

/*! \brief Fused-multiply-add for SIMD4 float. Result is a*b+c.
 * \copydetails gmx_simd_fmadd_f
 */
#    define gmx_simd4_fmadd_f    gmx_simd_fmadd_f

/*! \brief Fused-multiply-subtract for SIMD4 float. Result is a*b-c.
 * \copydetails gmx_simd_fmsub_f
 */
#    define gmx_simd4_fmsub_f    gmx_simd_fmsub_f

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b+c.
 * \copydetails gmx_simd_fnmadd_f
 */
#    define gmx_simd4_fnmadd_f   gmx_simd_fnmadd_f

/*! \brief Fused-negated-multiply-add for SIMD4 float. Result is -a*b-c.
 * \copydetails gmx_simd_fnmsub_f
 */
#    define gmx_simd4_fnmsub_f   gmx_simd_fnmsub_f

/*! \brief Lookup of approximate 1/sqrt(x) for SIMD4 float.
 * \copydetails gmx_simd_rsqrt_f
 */
#    define gmx_simd4_rsqrt_f    gmx_simd_rsqrt_f

/*! \brief Floating-point absolute value for SIMD4 float.
 * \copydetails gmx_simd_fabs_f
 */
#    define gmx_simd4_fabs_f     gmx_simd_fabs_f

/*! \brief Floating-point negate for SIMD4 float.
 * \copydetails gmx_simd_fneg_f
 */
#    define gmx_simd4_fneg_f     gmx_simd_fneg_f

/*! \brief Set each SIMD4 float element to the largest from two variables.
 * \copydetails gmx_simd_max_f
 */
#    define gmx_simd4_max_f      gmx_simd_max_f

/*! \brief Set each SIMD4 float element to the smallest from two variables.
 * \copydetails gmx_simd_min_f
 */
#    define gmx_simd4_min_f      gmx_simd_min_f

/*! \brief Round to nearest integer value for SIMD4 float.
 * \copydetails gmx_simd_round_f
 */
#    define gmx_simd4_round_f    gmx_simd_round_f

/*! \brief Round to largest integral value for SIMD4 float.
 * \copydetails gmx_simd_trunc_f
 */
#    define gmx_simd4_trunc_f    gmx_simd_trunc_f

/*! \brief Return dot product of two single precision SIMD4 variables.
 *
 * The dot product is calculated between the first three elements in the two
 * vectors, while the fourth is ignored. The result is returned as a scalar.
 *
 * \param a vector1
 * \param b vector2
 * \result a[0]*b[0]+a[1]*b[1]+a[2]*b[2], returned as scalar. Last element is ignored.
 */
static gmx_inline float
gmx_simd4_dotproduct3_f(gmx_simd_float_t a, gmx_simd_float_t b)
{
    return a.r[0] * b.r[0] + a.r[1] * b.r[1] + a.r[2] * b.r[2];
}

/*! \brief SIMD4 float transpose
 *
 * \param[in,out] v0  Row 0 on input, column 0 on output
 * \param[in,out] v1  Row 1 on input, column 1 on output
 * \param[in,out] v2  Row 2 on input, column 2 on output
 * \param[in,out] v3  Row 3 on input, column 3 on output
 */
static gmx_inline void
gmx_simd4_transpose_f(gmx_simd4_float_t * v0, gmx_simd4_float_t * v1,
                      gmx_simd4_float_t * v2, gmx_simd4_float_t * v3)
{
    gmx_simd4_float_t t0, t1, t2, t3;
    t0       = *v0;
    t1       = *v1;
    t2       = *v2;
    t3       = *v3;
    v0->r[0] = t0.r[0];
    v0->r[1] = t1.r[0];
    v0->r[2] = t2.r[0];
    v0->r[3] = t3.r[0];
    v1->r[0] = t0.r[1];
    v1->r[1] = t1.r[1];
    v1->r[2] = t2.r[1];
    v1->r[3] = t3.r[1];
    v2->r[0] = t0.r[2];
    v2->r[1] = t1.r[2];
    v2->r[2] = t2.r[2];
    v2->r[3] = t3.r[2];
    v3->r[0] = t0.r[3];
    v3->r[1] = t1.r[3];
    v3->r[2] = t2.r[3];
    v3->r[3] = t3.r[3];
}

/*! \brief SIMD4 variable type to use for logical comparisons on floats.
 * \copydetails gmx_simd_fbool_t
 */
#    define gmx_simd4_fbool_t   gmx_simd_fbool_t

/*! \brief Equality comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmpeq_f
 */
#    define gmx_simd4_cmpeq_f   gmx_simd_cmpeq_f

/*! \brief Return true for nonzero SIMD4 elements (i.e., with bits set)
 * \copydetails gmx_simd_cmpnz_f
 */
#    define gmx_simd4_cmpnz_f   gmx_simd_cmpnz_f

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmplt_f
 */
#    define gmx_simd4_cmplt_f   gmx_simd_cmplt_f

/*! \brief Less-than comparison of two single precision SIMD4.
 * \copydetails gmx_simd_cmple_f
 */
#    define gmx_simd4_cmple_f   gmx_simd_cmple_f

/*! \brief Logical AND on float SIMD4 booleans.
 * \copydetails gmx_simd_and_fb
 */
#    define gmx_simd4_and_fb gmx_simd_and_fb

/*! \brief Logical OR on float SIMD4 booleans.
 * \copydetails gmx_simd_or_fb
 */
#    define gmx_simd4_or_fb gmx_simd_or_fb

/*! \brief Returns non-zero if any of the SIMD4 boolean in x is True.
 * \copydetails gmx_simd_anytrue_fb
 */
#    define gmx_simd4_anytrue_fb gmx_simd_anytrue_fb

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 * \copydetails gmx_simd_blendzero_f
 */
#    define gmx_simd4_blendzero_f gmx_simd_blendzero_f

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 * \copydetails gmx_simd_blendnotzero_f
 */
#    define gmx_simd4_blendnotzero_f gmx_simd_blendnotzero_f

/*! \brief Vector-blend instruction form SIMD4 float.
 * \copydetails gmx_simd_blendv_f
 */
#    define gmx_simd4_blendv_f  gmx_simd_blendv_f

/*! \brief Return sum of all elements in SIMD4 float.
 * \copydetails gmx_simd_reduce_f
 */
#    define gmx_simd4_reduce_f  gmx_simd_reduce_f

#else /* GMX_SIMD_FLOAT_WIDTH!=4 */
#    undef GMX_SIMD4_HAVE_FLOAT
#endif


#if (GMX_SIMD_DOUBLE_WIDTH == 4) || (defined DOXYGEN)

/*! \brief SIMD4 double type. Available with \ref GMX_SIMD4_HAVE_DOUBLE.
 *
 * Unless you specifically want a double-precision type you should check
 * \ref gmx_simd4_real_t instead.
 *
 * While the SIMD4 datatype is identical to the normal SIMD type in the
 * reference implementation, this will often not be the case for
 * other architectures.
 */
#    define gmx_simd4_double_t   gmx_simd_double_t

/*! \brief Double precision SIMD4 load aligned.
 * \copydetails gmx_simd_load_d
 */
#    define gmx_simd4_load_d     gmx_simd_load_d

/*! \brief Double precision SIMD4 load single value to all elements.
 * \copydetails gmx_simd_load1_d
 */
#    define gmx_simd4_load1_d    gmx_simd_load1_d

/*! \brief Double precision SIMD4 set all elements from value.
 * \copydetails gmx_simd_set1_d
 */
#    define gmx_simd4_set1_d     gmx_simd_set1_d

/*! \brief Double precision SIMD4 store to aligned memory.
 * \copydetails gmx_simd_store_d
 */
#    define gmx_simd4_store_d   gmx_simd_store_d

/*! \brief Load unaligned SIMD4 double.
 * \copydetails gmx_simd_loadu_d
 */
#    define gmx_simd4_loadu_d   gmx_simd_loadu_d

/*! \brief Store unaligned SIMD4 double.
 * \copydetails gmx_simd_storeu_d
 */
#    define gmx_simd4_storeu_d  gmx_simd_storeu_d

/*! \brief Set all elements in SIMD4 double to 0.0.
 * \copydetails gmx_simd_setzero_d
 */
#    define gmx_simd4_setzero_d gmx_simd_setzero_d

/*! \brief Bitwise and for two SIMD4 double variables.
 * \copydetails gmx_simd_and_d
 */
#    define gmx_simd4_and_d     gmx_simd_and_d

/*! \brief Bitwise andnot for SIMD4 double. c=(~a) & b.
 * \copydetails gmx_simd_andnot_d
 */
#    define gmx_simd4_andnot_d  gmx_simd_andnot_d

/*! \brief Bitwise or for SIMD4 double.
 * \copydetails gmx_simd_or_d
 */
#    define gmx_simd4_or_d      gmx_simd_or_d

/*! \brief Bitwise xor for SIMD4 double.
 * \copydetails gmx_simd_xor_d
 */
#    define gmx_simd4_xor_d     gmx_simd_xor_d

/*! \brief Add two SIMD4 double values.
 * \copydetails gmx_simd_add_d
 */
#    define gmx_simd4_add_d     gmx_simd_add_d

/*! \brief Subtract two SIMD4 double values.
 * \copydetails gmx_simd_sub_d
 */
#    define gmx_simd4_sub_d     gmx_simd_sub_d

/*! \brief Multiply two SIMD4 double values.
 * \copydetails gmx_simd_mul_d
 */
#    define gmx_simd4_mul_d     gmx_simd_mul_d

/*! \brief Fused-multiply-add for SIMD4 double. Result is a*b+c.
 * \copydetails gmx_simd_fmadd_d
 */
#    define gmx_simd4_fmadd_d   gmx_simd_fmadd_d

/*! \brief Fused-multiply-subtract for SIMD4 double. Result is a*b-c.
 * \copydetails gmx_simd_fmsub_d
 */
#    define gmx_simd4_fmsub_d   gmx_simd_fmsub_d

/*! \brief Fused-negated-multiply-add for SIMD4 double. Result is -a*b+c.
 * \copydetails gmx_simd_fnmadd_d
 */
#    define gmx_simd4_fnmadd_d  gmx_simd_fnmadd_d

/*! \brief Fused-negated-multiply-sub for SIMD4 double. Result is -a*b-c.
 * \copydetails gmx_simd_fnmsub_d
 */
#    define gmx_simd4_fnmsub_d  gmx_simd_fnmsub_d

/*! \brief SIMD4 double 1.0/sqrt(x) lookup.
 * \copydetails gmx_simd_rsqrt_d
 */
#    define gmx_simd4_rsqrt_d   gmx_simd_rsqrt_d

/*! \brief SIMD4 double Floating-point fabs().
 * \copydetails gmx_simd_fabs_d
 */
#    define gmx_simd4_fabs_d    gmx_simd_fabs_d

/*! \brief SIMD4 double floating-point negate.
 * \copydetails gmx_simd_fneg_d
 */
#    define gmx_simd4_fneg_d    gmx_simd_fneg_d

/*! \brief Set each SIMD4 element to the largest from two variables.
 * \copydetails gmx_simd_max_d
 */
#    define gmx_simd4_max_d     gmx_simd_max_d

/*! \brief Set each SIMD4 element to the smallest from two variables.
 * \copydetails gmx_simd_min_d
 */
#    define gmx_simd4_min_d     gmx_simd_min_d

/*!  \brief Round SIMD4 double to nearest integer value (in floating-point format).
 * \copydetails gmx_simd_round_d
 */
#    define gmx_simd4_round_d   gmx_simd_round_d

/*! \brief Truncate SIMD4 double, i.e. round towards zero.
 * \copydetails gmx_simd_trunc_d
 */
#    define gmx_simd4_trunc_d   gmx_simd_trunc_d

/*! \brief Return dot product of two double precision SIMD4 variables.
 * \copydetails gmx_simd_setzero_f
 */
static gmx_inline double
gmx_simd4_dotproduct3_d(gmx_simd_double_t a, gmx_simd_double_t b)
{
    return a.r[0] * b.r[0] + a.r[1] * b.r[1] + a.r[2] * b.r[2];
}

/*! \brief SIMD4 double transpose
 *
 * \param[in,out] v0  Row 0 on input, column 0 on output
 * \param[in,out] v1  Row 1 on input, column 1 on output
 * \param[in,out] v2  Row 2 on input, column 2 on output
 * \param[in,out] v3  Row 3 on input, column 3 on output
 */
static gmx_inline void
gmx_simd4_transpose_d(gmx_simd4_double_t * v0, gmx_simd4_double_t * v1,
                      gmx_simd4_double_t * v2, gmx_simd4_double_t * v3)
{
    gmx_simd4_double_t t0, t1, t2, t3;
    t0       = *v0;
    t1       = *v1;
    t2       = *v2;
    t3       = *v3;
    v0->r[0] = t0.r[0];
    v0->r[1] = t1.r[0];
    v0->r[2] = t2.r[0];
    v0->r[3] = t3.r[0];
    v1->r[0] = t0.r[1];
    v1->r[1] = t1.r[1];
    v1->r[2] = t2.r[1];
    v1->r[3] = t3.r[1];
    v2->r[0] = t0.r[2];
    v2->r[1] = t1.r[2];
    v2->r[2] = t2.r[2];
    v2->r[3] = t3.r[2];
    v3->r[0] = t0.r[3];
    v3->r[1] = t1.r[3];
    v3->r[2] = t2.r[3];
    v3->r[3] = t3.r[3];
}

/*! \brief SIMD4 variable type to use for logical comparisons on doubles.
 * \copydetails gmx_simd_dbool_t
 */
#    define gmx_simd4_dbool_t   gmx_simd_dbool_t

/*! \brief Equality comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmpeq_d
 */
#    define gmx_simd4_cmpeq_d   gmx_simd_cmpeq_d

/*! \brief Return true for nonzero SIMD4 elements (i.e., with bits set)
 * \copydetails gmx_simd_cmpnz_d
 */
#    define gmx_simd4_cmpnz_d   gmx_simd_cmpnz_d

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmplt_d
 */
#    define gmx_simd4_cmplt_d   gmx_simd_cmplt_d

/*! \brief Less-than comparison of two double precision SIMD4 values.
 * \copydetails gmx_simd_cmple_d
 */
#    define gmx_simd4_cmple_d   gmx_simd_cmple_d

/*! \brief Logical AND on double SIMD4 booleans.
 * \copydetails gmx_simd_and_db
 */
#    define gmx_simd4_and_db gmx_simd_and_db

/*! \brief Logical OR on double SIMD4 booleans.
 * \copydetails gmx_simd_or_db
 */
#    define gmx_simd4_or_db gmx_simd_or_db

/*! \brief Returns non-zero if any of the SIMD4 booleans in x is True.
 * \copydetails gmx_simd_anytrue_db
 */
#    define gmx_simd4_anytrue_db gmx_simd_anytrue_db

/*! \brief Select from double precision SIMD4 variable where boolean is true.
 * \copydetails gmx_simd_blendzero_d
 */
#    define gmx_simd4_blendzero_d gmx_simd_blendzero_d

/*! \brief Select from double precision SIMD4 variable where boolean is false.
 * \copydetails gmx_simd_blendnotzero_d
 */
#    define gmx_simd4_blendnotzero_d gmx_simd_blendnotzero_d

/*! \brief Vector-blend instruction for SIMD4 double.
 * \copydetails gmx_simd_blendv_d
 */
#    define gmx_simd4_blendv_d  gmx_simd_blendv_d

/*! \brief Return sum of all elements in SIMD4 double.
 * \copydetails gmx_simd_reduce_d
 */
#    define gmx_simd4_reduce_d  gmx_simd_reduce_d

#else /* GMX_SIMD4_DOUBLE_WIDTH!=4 */
#    undef GMX_SIMD4_HAVE_DOUBLE
#endif

/*! \} */


/*! \brief Return 1 if SIMD floating-point ops have overflowed, and reset check.

 * This function to check whether SIMD operations have resulted in overflow,
 * and returns 1 if it occured, 0 otherwise.
 * For now, this is unfortunately a dummy for all architectures except x86.
 */
static int
gmx_simd_check_and_reset_overflow(void)
{
    return 0;
}

/*! \} */
/*! \endcond */

#endif /* GMX_SIMD_IMPL_REFERENCE_H */
