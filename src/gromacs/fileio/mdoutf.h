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

#ifndef GMX_FILEIO_MDOUTF_H
#define GMX_FILEIO_MDOUTF_H

#include <stdio.h>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/oenv.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_mtop_t;

typedef struct gmx_mdoutf *gmx_mdoutf_t;

/*! \brief Allocate and initialize object to manager trajectory writing output
 *
 * Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */
gmx_mdoutf_t init_mdoutf(FILE              *fplog,
                         int                nfile,
                         const t_filenm     fnm[],
                         int                mdrun_flags,
                         const t_commrec   *cr,
                         const t_inputrec  *ir,
                         struct gmx_mtop_t *mtop,
                         const output_env_t oenv,
                         gmx_wallcycle_t    wcycle);

/*! \brief Getter for file pointer */
FILE *mdoutf_get_fp_field(gmx_mdoutf_t of);

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
 * the master node only when necessary.
 */
void mdoutf_write_to_trajectory_files(FILE *fplog, t_commrec *cr,
                                      gmx_mdoutf_t of,
                                      int mdof_flags,
                                      struct gmx_mtop_t *top_global,
                                      gmx_int64_t step, double t,
                                      t_state *state_local, t_state *state_global,
                                      rvec *f_local, rvec *f_global);

#define MDOF_X            (1<<0)
#define MDOF_V            (1<<1)
#define MDOF_F            (1<<2)
#define MDOF_X_COMPRESSED (1<<3)
#define MDOF_CPT          (1<<4)
#define MDOF_IMD          (1<<5)

#ifdef __cplusplus
}
#endif

#endif /* GMX_FILEIO_MDOUTF_H */
