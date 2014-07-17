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
#ifndef GMX_FILEIO_TNGIO_FOR_TOOLS_H
#define GMX_FILEIO_TNGIO_FOR_TOOLS_H

#include <stdio.h>

#include "tng/tng_io_fwd.h"

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

struct gmx_mtop_t;
struct t_trxframe;

/*! \brief Prepare to write TNG output from trajectory conversion tools */
void gmx_prepare_tng_writing(const char              *filename,
                             char                     mode,
                             tng_trajectory_t        *in,
                             tng_trajectory_t        *out,
                             int                      nAtoms,
                             const struct gmx_mtop_t *mtop,
                             const atom_id           *index,
                             const char              *indexGroupName);

/*! \brief Write a trxframe to a TNG file
 *
 * \param output Trajectory to write to
 * \param frame  Frame data to write
 * \param natoms Number of atoms to actually write
 *
 * The natoms field in frame is the number of atoms in the system. The
 * parameter natoms supports writing an index-group subset of the
 * atoms.
 */
void gmx_write_tng_from_trxframe(tng_trajectory_t        output,
                                 struct t_trxframe      *frame,
                                 int                     natoms);

/*! \brief Creates a molecule containing only the indexed atoms and sets
 * the number of all other molecules to 0. Works similar to a
 * selection group. */
void gmx_tng_setup_atom_subgroup(tng_trajectory_t tng,
                                 const int        nind,
                                 const atom_id   *ind,
                                 const char      *name);

/*! \brief Read the first/next TNG frame. */
gmx_bool gmx_read_next_tng_frame(tng_trajectory_t            input,
                                 struct t_trxframe          *fr,
                                 gmx_int64_t                *requestedIds,
                                 int                         numRequestedIds);

/*! \brief Print the molecule system to stream */
void gmx_print_tng_molecule_system(tng_trajectory_t input,
                                   FILE            *stream);

/*! \brief Get a list of block IDs present in the next frame with data. */
gmx_bool gmx_get_tng_data_block_types_of_next_frame(tng_trajectory_t     input,
                                                    int                  frame,
                                                    int                  nRequestedIds,
                                                    gmx_int64_t         *requestedIds,
                                                    gmx_int64_t         *nextFrame,
                                                    gmx_int64_t         *nBlocks,
                                                    gmx_int64_t        **blockIds);

/*! \brief Get data of the next frame with data from the data block
 * with the specified block ID. */
gmx_bool gmx_get_tng_data_next_frame_of_block_type(tng_trajectory_t     input,
                                                   gmx_int64_t          blockId,
                                                   real               **values,
                                                   gmx_int64_t         *frameNumber,
                                                   double              *frameTime,
                                                   gmx_int64_t         *nValuesPerFrame,
                                                   gmx_int64_t         *nAtoms,
                                                   real                *prec,
                                                   char                *name,
                                                   int                  maxLen,
                                                   gmx_bool            *bOK);

#ifdef __cplusplus
}
#endif

#endif
