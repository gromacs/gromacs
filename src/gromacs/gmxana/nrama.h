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

#ifndef GMX_GMXANA_NRAMA_H
#define GMX_GMXANA_NRAMA_H

#include "gromacs/fileio/trxio.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    gmx_bool bShow;
    char    *label;
    int      iphi, ipsi; /* point in the dih array of xr... */
} t_phipsi;

typedef struct {
    atom_id ai[4];
    int     mult;
    real    phi0;
    real    ang;
} t_dih;

typedef struct {
    int          ndih;
    t_dih       *dih;
    int          npp;
    t_phipsi    *pp;
    t_trxstatus *traj;
    int          natoms;
    int          amin, amax;
    real         t;
    rvec        *x;
    matrix       box;
    t_idef      *idef;
    int          ePBC;
    output_env_t oenv;
} t_xrama;

t_topology *init_rama(const output_env_t oenv, const char *infile,
                      const char *topfile, t_xrama *xr, int mult);

gmx_bool new_data(t_xrama *xr);

#ifdef __cplusplus
}
#endif

#endif  /* GMX_GMXANA_NRAMA_H */
