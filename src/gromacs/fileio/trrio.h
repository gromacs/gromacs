/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_FILEIO_TRRIO_H
#define GMX_FILEIO_TRRIO_H

#include <cstdint>

#include <filesystem>
#include <limits>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/**************************************************************
 *
 * These routines handle trr (trajectory) I/O, they read and
 * write trr files. The routines should be able to read single
 * and double precision files without the user noting it.
 * The files are backward compatible, therefore the header holds
 * some unused variables.
 *
 * The routines in the corresponding c-file trrio.cpp
 * are based on the lower level routines in gmxfio.cpp
 * The file handle returned from gmx_trr_open()
 * can also be used with the routines in gmxfio.h
 *
 * Note that TRR was designed to represent a step number as a default
 * integer, which depends on the implementation, but is typically and
 * 32 bit. We didn't design the format to be extensible, so we can't
 * fix the fact that after 2^31 frames, step numbers will wrap to be
 * negative. Fortunately, this tends not to cause serious problems,
 * and we've fixed it in TNG. Meanwhile, the implementation pretends
 * to the rest of GROMACS that it functions with int64_t like all
 * other step numbers, but the actual range in practice depends on the
 * defaults of the implementation in use now (or when the file was
 * written).
 *
 **************************************************************/

struct t_fileio;

static constexpr int64_t sc_trrMaxAtomCount = std::numeric_limits<unsigned int>::max() / 3;

/* This struct describes the order and the  */
/* sizes of the structs in a trr file, sizes are given in bytes. */
typedef struct gmx_trr_header_t
{
    gmx_bool     bDouble;   /* Double precision?                   */
    int          ir_size;   /* Backward compatibility              */
    int          e_size;    /* Backward compatibility              */
    int          box_size;  /* Non zero if a box is present        */
    int          vir_size;  /* Backward compatibility              */
    int          pres_size; /* Backward compatibility              */
    int          top_size;  /* Backward compatibility              */
    int          sym_size;  /* Backward compatibility              */
    unsigned int x_size;    /* Non zero if coordinates are present */
    unsigned int v_size;    /* Non zero if velocities are present  */
    unsigned int f_size;    /* Non zero if forces are present      */

    int     natoms;    /* The total number of atoms           */
    int64_t step;      /* Current step number                 */
    int     nre;       /* Backward compatibility              */
    real    t;         /* Current time                        */
    real    lambda;    /* Current value of lambda             */
    int     fep_state; /* Current value of alchemical state   */
} gmx_trr_header_t;

struct t_fileio* gmx_trr_open(const std::filesystem::path& fn, const char* mode);
/* Open a trr file */

void gmx_trr_close(struct t_fileio* fio);
/* Close it */

gmx_bool gmx_trr_read_frame_header(struct t_fileio* fio, gmx_trr_header_t* header, gmx_bool* bOK);
/* Read the header of a trr file. Return FALSE if there is no frame.
 * bOK will be FALSE when the header is incomplete.
 */

gmx_bool gmx_trr_read_frame_data(struct t_fileio* fio, gmx_trr_header_t* sh, rvec* box, rvec* x, rvec* v, rvec* f);
/* Extern read a frame except the header (that should be pre-read,
 * using routine gmx_trr_read_frame_header(), see above) from a trr file.
 * Return FALSE on error
 */

gmx_bool gmx_trr_read_frame(struct t_fileio* fio,
                            int64_t*         step,
                            real*            t,
                            real*            lambda,
                            rvec*            box,
                            int*             natoms,
                            rvec*            x,
                            rvec*            v,
                            rvec*            f);
/* Read a trr frame, including the header from fp. box, x, v, f may
 * be NULL, in which case the data will be skipped over.
 * return FALSE on error
 */

void gmx_trr_write_frame(struct t_fileio* fio,
                         int64_t          step,
                         real             t,
                         real             lambda,
                         const rvec*      box,
                         int              natoms,
                         const rvec*      x,
                         const rvec*      v,
                         const rvec*      f);
/* Write a trr frame to file fp, box, x, v, f may be NULL */

void gmx_trr_read_single_header(const std::filesystem::path& fn, gmx_trr_header_t* header);
/* Read the header of a trr file from fn, and close the file afterwards.
 */

void gmx_trr_read_single_frame(const std::filesystem::path& fn,
                               int64_t*                     step,
                               real*                        t,
                               real*                        lambda,
                               rvec*                        box,
                               int*                         natoms,
                               rvec*                        x,
                               rvec*                        v,
                               rvec*                        f);
/* Read a single trr frame from file fn, which is closed afterwards
 */

void gmx_trr_write_single_frame(const std::filesystem::path& fn,
                                int64_t                      step,
                                real                         t,
                                real                         lambda,
                                const rvec*                  box,
                                int                          natoms,
                                const rvec*                  x,
                                const rvec*                  v,
                                const rvec*                  f);
/* Write a single trr frame to file fn, which is closed afterwards */


#endif
