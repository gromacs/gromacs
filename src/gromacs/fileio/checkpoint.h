/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

class energyhistory_t;
struct gmx_file_position_t;
struct ObservablesHistory;
struct t_commrec;
struct t_fileio;
struct t_inputrec;
class t_state;
struct t_trxframe;

/* the name of the environment variable to disable fsync failure checks with */
#define GMX_IGNORE_FSYNC_FAILURE_ENV "GMX_IGNORE_FSYNC_FAILURE"

// TODO Replace this mechanism with std::array<char, 1024> or similar.
#define CPTSTRLEN 1024

// TODO Expand this into being a container of all data for
// serialization of a checkpoint, which can be stored by the caller
// (e.g. so that mdrun doesn't have to open the checkpoint twice).
// This will separate issues of allocation from those of
// serialization, help separate comparison from reading, and have
// better defined transformation functions to/from trajectory frame
// data structures.
//
// Several fields were once written to checkpoint file headers, but
// have been removed. So that old files can continue to be read,
// the names of such fields contain the string "_UNUSED" so that it
// is clear they should not be used.
struct CheckpointHeaderContents
{
    int         file_version;
    char        version[CPTSTRLEN];
    char        btime_UNUSED[CPTSTRLEN];
    char        buser_UNUSED[CPTSTRLEN];
    char        bhost_UNUSED[CPTSTRLEN];
    int         double_prec;
    char        fprog[CPTSTRLEN];
    char        ftime[CPTSTRLEN];
    int         eIntegrator;
    int         simulation_part;
    int64_t     step;
    double      t;
    int         nnodes;
    ivec        dd_nc;
    int         npme;
    int         natoms;
    int         ngtc;
    int         nnhpres;
    int         nhchainlength;
    int         nlambda;
    int         flags_state;
    int         flags_eks;
    int         flags_enh;
    int         flagsPullHistory;
    int         flags_dfh;
    int         flags_awhh;
    int         nED;
    int         eSwapCoords;
};

/* Write a checkpoint to <fn>.cpt
 * Appends the _step<step>.cpt with bNumberAndKeep,
 * otherwise moves the previous <fn>.cpt to <fn>_prev.cpt
 */
void write_checkpoint(const char *fn, gmx_bool bNumberAndKeep,
                      FILE *fplog, const t_commrec *cr,
                      ivec domdecCells, int nppnodes,
                      int eIntegrator, int simulation_part,
                      gmx_bool bExpanded, int elamstats,
                      int64_t step, double t,
                      t_state *state, ObservablesHistory *observablesHistory);

/* Loads a checkpoint from fn for run continuation.
 * Generates a fatal error on system size mismatch.
 * The master node reads the file
 * and communicates all the modified number of steps,
 * but not the state itself.
 * With reproducibilityRequested warns about version, build, #ranks differences.
 */
void load_checkpoint(const char *fn, t_fileio *logfio,
                     const t_commrec *cr, const ivec dd_nc,
                     t_inputrec *ir, t_state *state,
                     ObservablesHistory *observablesHistory,
                     gmx_bool reproducibilityRequested);

/* Read everything that can be stored in t_trxframe from a checkpoint file */
void read_checkpoint_trxframe(struct t_fileio *fp, t_trxframe *fr);

/* Print the complete contents of checkpoint file fn to out */
void list_checkpoint(const char *fn, FILE *out);

/*!\brief Read simulation step and part from a checkpoint file
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
                                   int64_t     *step);

/*!\brief Return header information from an open checkpoint file.
 *
 * Used by mdrun to handle restarts
 *
 * \param[in]  fp               Handle to open checkpoint file
 * \param[out] outputfiles      Container of output file names from the previous run. */
CheckpointHeaderContents
read_checkpoint_simulation_part_and_filenames(t_fileio                         *fp,
                                              std::vector<gmx_file_position_t> *outputfiles);

#endif
