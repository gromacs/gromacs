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

#include "filenm.h"
#include <stdio.h>
#include "../legacyheaders/types/simple.h"
#include "enxio.h"
#include "gmxfio.h"

typedef struct {
    t_fileio     *fp_trn;
    t_fileio     *fp_xtc;
    int           xtc_prec;
    ener_file_t   fp_ene;
    const char   *fn_cpt;
    gmx_bool      bKeepAndNumCPT;
    int           eIntegrator;
    gmx_bool      bExpanded;
    int           elamstats;
    int           simulation_part;
    FILE         *fp_dhdl;
    FILE         *fp_field;
    int           natoms_global;
    int           natoms_xtc;
    gmx_groups_t *groups; /* for XTC writing */
} gmx_mdoutf_t;

gmx_mdoutf_t *init_mdoutf(int nfile, const t_filenm fnm[],
                          int mdrun_flags,
                          const t_commrec *cr, const t_inputrec *ir,
                          gmx_mtop_t *top_global,
                          const output_env_t oenv);
/* Returns a pointer to a data structure with all output file pointers
 * and names required by mdrun.
 */

void done_mdoutf(gmx_mdoutf_t *of);
/* Close all open output files and free the of pointer */

#define MDOF_X   (1<<0)
#define MDOF_V   (1<<1)
#define MDOF_F   (1<<2)
#define MDOF_XTC (1<<3)
#define MDOF_CPT (1<<4)

#endif /* GMX_FILEIO_MDOUTF_H */
