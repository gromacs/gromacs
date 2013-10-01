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

#ifndef GMX_FILEIO_TNGIO_H
#define GMX_FILEIO_TNGIO_H

#include "tng_io.h"

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

void tng_open(const char       *filename,
              char              mode,
              tng_trajectory_t *tng_data_p);
/* Open a TNG trajectory file */

void tng_close(tng_trajectory_t *tng);
/* Finish writing a TNG trajectory file */

void tng_prepare_md_writing(tng_trajectory_t  tng,
                            const gmx_mtop_t *mtop,
                            const t_inputrec *ir);
/* Do all TNG preparation for MD trajectory writing */

void tng_prepare_low_prec_writing(tng_trajectory_t  tng,
                                  const gmx_mtop_t *mtop,
                                  const t_inputrec *ir);
/* Do all TNG preparation for low precision MD trajectory
 * writing */

void tng_add_selection_groups(tng_trajectory_t  tng,
                              const gmx_mtop_t *mtop);
/* Create a TNG molecule representing the selection groups
 * to write */

void tng_add_mtop(tng_trajectory_t  tng,
                  const gmx_mtop_t *mtop);
/* Converts the current topology to TNG molecular data */

void tng_set_frames_per_frame_set(tng_trajectory_t  tng,
                                  const gmx_bool    bUseLossyCompression,
                                  const t_inputrec *ir);
/* Set the number of frames per frame set according to output
 * intervals. The default is that 100 frames are written of
 * the data that is written most often. */

void tng_set_writing_intervals(tng_trajectory_t  tng,
                               const gmx_bool    bUseLossyCompression,
                               const t_inputrec *ir);
/* Set output frequency and number of frames per frame set */

void fwrite_tng(tng_trajectory_t tng,
                const gmx_bool   bUseLossyCompression,
                int              step,
                real             elapsedPicoSeconds,
                real             lambda,
                const rvec      *box,
                int              natoms,
                const rvec      *x,
                const rvec      *v,
                const rvec      *f);
/* Write a frame to a TNG file fp, box, x, v, f may be NULL */

#ifdef __cplusplus
}
#endif

#endif /* GMX_FILEIO_TNGIO_H */
