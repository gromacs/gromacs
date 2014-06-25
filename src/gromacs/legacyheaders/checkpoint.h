/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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

#ifndef _checkpoint_h
#define _checkpoint_h


#include "typedefs.h"
#include "../fileio/gmxfio.h"
#include "../fileio/filenm.h"

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

/* Read just the simulation 'generation' and with bAppendReq check files.
 * This is necessary already at the beginning of mdrun,
 * to be able to rename the logfile correctly.
 * When file appending is requested, checks which output files are present:
 * all present: return TRUE,
 * none present: return FALSE,
 * part present: fatal error.
 * When TRUE is returned, bAddPart will tell whether the simulation part
 * needs to be added to the output file name.
 */
gmx_bool read_checkpoint_simulation_part(const char *filename, int *simulation_part,
                                         gmx_int64_t *step, t_commrec *cr,
                                         gmx_bool bAppendReq,
                                         int nfile, const t_filenm fnm[],
                                         const char *part_suffix, gmx_bool *bAddPart);

#ifdef __cplusplus
}
#endif

#endif
