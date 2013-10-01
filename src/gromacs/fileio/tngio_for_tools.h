/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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

#ifndef GMX_FILEIO_TNGIO_FOR_TOOLS_H
#define GMX_FILEIO_TNGIO_FOR_TOOLS_H

#if !defined(GMX_FILEIO_TNGIO_H) && !defined(GMX_FILEIO_TRXIO_H)
typedef struct tng_trajectory *tng_trajectory_t;
#endif

#define TNG_TRAJ_POSITIONS              0x0000000010000001LL

#ifdef USE_STD_INTTYPES_H
#include <inttypes.h>
#else
/* Visual Studio does not contain inttypes.h and stdint.h. Some defines and
 * typedefs are used from the GNU C Library */
#ifdef _MSC_VER

typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
typedef __int64 int64_t;
typedef unsigned __int64 uint64_t;

#else
#include <stdint.h>
#endif

/* This is from inttypes.h  (GNU C Library) */
/* The ISO C99 standard specifies that these macros must only be
   defined if explicitly requested.  */
#if !defined __cplusplus || defined __STDC_FORMAT_MACROS

# if __WORDSIZE == 64
#  define __PRI64_PREFIX        "l"
#  define __PRIPTR_PREFIX       "l"
# else
#  define __PRI64_PREFIX        "ll"
#  define __PRIPTR_PREFIX
# endif

/* From stdint.h (GNU C Library) */
/* Macros for printing format specifiers. */
/* Decimal notation.  */
# define PRId8          "d"
# define PRId16         "d"
# define PRId32         "d"
# define PRId64         __PRI64_PREFIX "d"
#endif

#endif

#include "gromacs/legacyheaders/typedefs.h"
#include "trxio.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void prepare_tng_writing(const char              *filename,
                         char                     mode,
                         tng_trajectory_t        *in,
                         tng_trajectory_t        *out,
                         int                      nAtoms);
/* TODO */

void open_tng_for_reading(const char              *fn,
                          tng_trajectory_t        *input);
/* TODO */

void write_tng_from_trxframe(tng_trajectory_t        output,
                             t_trxframe             *frame);
/* TODO */

void tng_tools_close(tng_trajectory_t *tng);
/* TODO */

gmx_bool read_next_tng_frame(tng_trajectory_t        input,
                             t_trxframe             *fr,
                             int64_t                *requestedIds);
/* Read the first/next TNG frame. */

void tng_set_compression_precision(tng_trajectory_t tng,
                                   real prec);
/* Set the precision of TNG compression */

void print_tng_molecule_system(tng_trajectory_t input,
                               FILE *stream);
/* Print the molecule system to stream */

gmx_bool get_tng_data_block_types_of_next_frame(tng_trajectory_t input,
                                                int frame,
                                                int nRequestedIds,
                                                int64_t *requestedIds,
                                                int64_t *nextFrame,
                                                int64_t *nBlocks,
                                                int64_t **blockIds);
/* Get a list of block IDs present in the next frame with data. */

gmx_bool get_tng_data_next_frame_of_block_type(tng_trajectory_t input,
                                               int64_t blockId,
                                               real **values,
                                               int64_t *frameNumber,
                                               double *frameTime,
                                               int64_t *nValuesPerFrame,
                                               int64_t *nAtoms,
                                               real *prec,
                                               char *name,
                                               int maxLen,
                                               gmx_bool *bOK);
/* Get data of the next frame with data from the data block with the specified
 * block ID */

#ifdef __cplusplus
}
#endif

#endif
