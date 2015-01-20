/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#ifndef _checkpoint_h
#define _checkpoint_h


#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the name of the environment variable to disable fsync failure checks with */
#define GMX_IGNORE_FSYNC_FAILURE_ENV "GMX_IGNORE_FSYNC_FAILURE"

/* Write a checkpoint to <fn>.cpt
 * Appends the _step<step>.cpt with bNumberAndKeep,
 * otherwise moves the previous <fn>.cpt to <fn>_prev.cpt
 */
void write_checkpoint(const char *fn, gmx_bool bNumberAndKeep,
                      FILE *fplog, t_commrec *cr,
                      int eIntegrator, int simulation_part,
                      gmx_bool bExpanded, int elamstats,
                      gmx_int64_t step, double t,
                      t_state *state);

/* Loads a checkpoint from fn for run continuation.
 * Generates a fatal error on system size mismatch.
 * The master node reads the file
 * and communicates all the modified number of steps and the parallel setup,
 * but not the state itself.
 * When bAppend is set, lock the log file and truncate the existing output
 * files so they can be appended.
 * With bAppend and bForceAppend: truncate anyhow if the system does not
 * support file locking.
 */
void load_checkpoint(const char *fn, FILE **fplog,
                     t_commrec *cr, ivec dd_nc,
                     t_inputrec *ir, t_state *state,
                     gmx_bool *bReadEkin,
                     gmx_bool bAppend, gmx_bool bForceAppend);

/* Read the state from checkpoint file.
 * Arrays in state that are NULL are allocated.
 * If bReadRNG=TRUE a RNG state compatible with the current
 * number of nodes was read.
 */
void read_checkpoint_state(const char *fn, int *simulation_part,
                           gmx_int64_t *step, double *t, t_state *state);

/* Read everything that can be stored in t_trxframe from a checkpoint file */
void read_checkpoint_trxframe(t_fileio *fp, t_trxframe *fr);

/* Print the complete contents of checkpoint file fn to out */
void list_checkpoint(const char *fn, FILE *out);

/* ! \brief Read simulation step and part from a checkpoint file
 *
 * Used by tune_pme to handle tuning with a checkpoint file as part of the input.
 *
 * \param[in]  filename         Name of checkpoint file
 * \param[out] simulation_part  The part of the simulation that wrote the checkpoint
 * \param[out] step             The final step number of the simulation that wrote the checkpoint
 *
 * The output variables will both contain 0 if filename is NULL, the file
 * does not exist, or is not readable. */
void read_checkpoint_part_and_step(const char  *filename,
                                   int         *simulation_part,
                                   gmx_int64_t *step);

/* ! \brief Read simulation part and output filenames from a checkpoint file
 *
 * Used by mdrun to handle restarts
 *
 * \param[in]  fp               Handle to open checkpoint file
 * \param[out] simulation_part  The part of the simulation that wrote the checkpoint
 * \param[out] nfiles           Number of output files from the previous run
 * \param[out] outputfiles      Pointer to array of output file names from the previous run. Pointer is allocated in this function.
 */
void
read_checkpoint_simulation_part_and_filenames(t_fileio             *fp,
                                              int                  *simulation_part,
                                              int                  *nfiles,
                                              gmx_file_position_t **outputfiles);

#ifdef __cplusplus
}
#endif

#endif
