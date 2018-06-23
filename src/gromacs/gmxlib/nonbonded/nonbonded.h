/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2018, by the GROMACS development team, led by
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
#ifndef GMX_GMXLIB_NONBONDED_NONBONDED_H
#define GMX_GMXLIB_NONBONDED_NONBONDED_H

#include <stdio.h>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_grppairener_t;
struct t_forcerec;

void
gmx_nonbonded_setup(t_forcerec *fr,
                    gmx_bool    bGenericKernelOnly);





void
gmx_nonbonded_set_kernel_pointers(FILE *       fplog,
                                  t_nblist *   nl,
                                  gmx_bool     bElecAndVdwSwitchDiffers);



#define GMX_NONBONDED_DO_FORCE          (1<<1)
#define GMX_NONBONDED_DO_SHIFTFORCE     (1<<2)
#define GMX_NONBONDED_DO_FOREIGNLAMBDA  (1<<3)
#define GMX_NONBONDED_DO_POTENTIAL      (1<<4)
#define GMX_NONBONDED_DO_SR             (1<<5)

void
do_nonbonded(const t_forcerec  *fr,
             rvec               x[],
             rvec               f_shortrange[],
             const t_mdatoms   *md,
             const t_blocka    *excl,
             gmx_grppairener_t *grppener,
             t_nrnb            *nrnb,
             real              *lambda,
             real               dvdlambda[],
             int                nls,
             int                eNL,
             int                flags);

#endif
