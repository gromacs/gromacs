/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_MDOUTF_H
#define GMX_MDLIB_MDOUTF_H

#include <stdio.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"

class energyhistory_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct MdrunOptions;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;
struct t_inputrec;

namespace gmx
{
class IMDOutputProvider;
}

typedef struct gmx_mdoutf *gmx_mdoutf_t;

/*! \brief Allocate and initialize object to manager trajectory writing output
 *
 * Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */
gmx_mdoutf_t init_mdoutf(FILE                   *fplog,
                         int                     nfile,
                         const t_filenm          fnm[],
                         const MdrunOptions     &mdrunOptions,
                         const t_commrec        *cr,
                         gmx::IMDOutputProvider *outputProvider,
                         const t_inputrec       *ir,
                         gmx_mtop_t             *mtop,
                         const gmx_output_env_t *oenv,
                         gmx_wallcycle_t         wcycle);

/*! \brief Getter for file pointer */
ener_file_t mdoutf_get_fp_ene(gmx_mdoutf_t of);

/*! \brief Getter for file pointer */
FILE *mdoutf_get_fp_dhdl(gmx_mdoutf_t of);

/*! \brief Getter for wallcycle timer */
gmx_wallcycle_t mdoutf_get_wcycle(gmx_mdoutf_t of);

/*! \brief Close TNG files if they are open.
 *
 * This also measures the time it takes to close the TNG
 * files.
 */
void mdoutf_tng_close(gmx_mdoutf_t of);

/*! \brief Close all open output files and free the of pointer */
void done_mdoutf(gmx_mdoutf_t of);

/*! \brief Routine that writes trajectory-like frames.
 *
 * Writes data to trn, xtc and/or checkpoint. What is written is
 * determined by the mdof_flags defined below. Data is collected to
 * the master node only when necessary. Without domain decomposition
 * only data from state_local is used and state_global is ignored.
 */
void mdoutf_write_to_trajectory_files(FILE *fplog, const t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      gmx_mtop_t *top_global,
                                      int64_t step, double t,
                                      t_state *state_local, t_state *state_global,
                                      ObservablesHistory *observablesHistory,
                                      gmx::ArrayRef<gmx::RVec> f_local);

/*! \brief Get the output interval of box size of uncompressed TNG output.
 * Returns 0 if no uncompressed TNG file is open.
 */
int mdoutf_get_tng_box_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of lambda of uncompressed TNG output.
 * Returns 0 if no uncompressed TNG file is open.
 */
int mdoutf_get_tng_lambda_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of box size of compressed TNG output.
 * Returns 0 if no compressed TNG file is open.
 */
int mdoutf_get_tng_compressed_box_output_interval(gmx_mdoutf_t of);

/*! \brief Get the output interval of lambda of compressed TNG output.
 * Returns 0 if no compressed TNG file is open.
 */
int mdoutf_get_tng_compressed_lambda_output_interval(gmx_mdoutf_t of);

#define MDOF_X                 (1<<0)
#define MDOF_V                 (1<<1)
#define MDOF_F                 (1<<2)
#define MDOF_X_COMPRESSED      (1<<3)
#define MDOF_CPT               (1<<4)
#define MDOF_IMD               (1<<5)
#define MDOF_BOX               (1<<6)
#define MDOF_LAMBDA            (1<<7)
#define MDOF_BOX_COMPRESSED    (1<<8)
#define MDOF_LAMBDA_COMPRESSED (1<<9)

#endif
