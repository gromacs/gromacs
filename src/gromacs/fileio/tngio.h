/*
 * This file is part of the GROMACS molecular simulation package.
 *
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

#ifndef GMX_FILEIO_TNGIO_H
#define GMX_FILEIO_TNGIO_H

#include "tng/tng_io_fwd.h"

#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

struct gmx_mtop_t;

/*! \brief Open a TNG trajectory file
 *
 * \param filename   Name of file to open
 * \param mode       Can be set to 'r', 'w' or 'a' for reading, writing or appending respectively.
 * \param tng_data_p Pointer to an allocated tng_trajectory_t into which a handle to a TNG trajectory will be stored.
 *
 * Handles all I/O errors internally via fatal error
 */
void gmx_tng_open(const char       *filename,
                  char              mode,
                  tng_trajectory_t *tng_data_p);

/*! \brief Finish writing a TNG trajectory file */
void gmx_tng_close(tng_trajectory_t *tng);

/*!\brief Add molecular topology information to TNG output (if
 * available)
 *
 * \param tng   Valid handle to a TNG trajectory
 * \param mtop  Pointer to a topology (can be NULL)
 */
void gmx_tng_add_mtop(tng_trajectory_t         tng,
                      const struct gmx_mtop_t *mtop);

/*! \brief Do all TNG preparation for full-precision whole-system
 * trajectory writing during MD simulations.
 *
 * \param tng   Valid handle to a TNG trajectory
 * \param mtop  Global topology
 * \param ir    Input settings (for writing frequencies)
 */
void gmx_tng_prepare_md_writing(tng_trajectory_t         tng,
                                const struct gmx_mtop_t *mtop,
                                const t_inputrec        *ir);

/*! \brief Set the default compression precision for TNG writing
 *
 * \param tng   Valid handle to a TNG trajectory
 * \param prec  GROMACS-style precision setting (i.e. 1000 for 3 digits of precision) */
void gmx_tng_set_compression_precision(tng_trajectory_t tng,
                                       real             prec);

/*! \brief Do all TNG preparation for low-precision selection-based
 * trajectory writing during MD simulations.
 *
 * \param tng   Valid handle to a TNG trajectory
 * \param mtop  Global topology
 * \param ir    Input settings (for writing frequencies)
 */
void gmx_tng_prepare_low_prec_writing(tng_trajectory_t         tng,
                                      const struct gmx_mtop_t *mtop,
                                      const t_inputrec        *ir);

/*! \brief Write a frame to a TNG file
 *
 * \param tng                  Valid handle to a TNG trajectory
 * \param bUseLossyCompression Whether to use lossy compression
 * \param step                 MD step number
 * \param elapsedPicoSeconds   Elapsed MD time
 * \param lambda               Free-energy lambda value
 * \param box                  Simulation box
 * \param nAtoms               Number of atoms (i.e. vector lengths)
 * \param x                    Vector of position coordinates
 * \param v                    Vector of elocities
 * \param f                    Vector of forces
 *
 * The pointers tng, x, v, f may be NULL, which triggers not writing
 * (that component). box can only be NULL if x is also NULL. */
void gmx_fwrite_tng(tng_trajectory_t tng,
                    const gmx_bool   bUseLossyCompression,
                    int              step,
                    real             elapsedPicoSeconds,
                    real             lambda,
                    const rvec      *box,
                    int              nAtoms,
                    const rvec      *x,
                    const rvec      *v,
                    const rvec      *f);

/*! \brief Write the current frame set to disk. Perform compression
 * etc.
 *
 * \param tng Valid handle to a TNG trajectory
 */
void fflush_tng(tng_trajectory_t tng);

/*! \brief Get the time (in picoseconds) of the final frame in the
 * trajectory.
 *
 * \param tng Valid handle to a TNG trajectory
 */
float gmx_tng_get_time_of_final_frame(tng_trajectory_t tng);

#ifdef __cplusplus
}
#endif

#endif /* GMX_FILEIO_TNGIO_H */
