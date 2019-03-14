/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 *
 * \brief This file contains declarations necessary for low-level
 * functions for computing energies and forces for position
 * restraints.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_listed-forces
 */

#ifndef GMX_LISTED_FORCES_POSITION_RESTRAINTS_H
#define GMX_LISTED_FORCES_POSITION_RESTRAINTS_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;
struct gmx_wallcycle;
struct t_forcerec;
struct t_idef;
struct t_lambda;
struct t_nrnb;
struct t_pbc;

namespace gmx
{
class ForceWithVirial;
}

/*! \brief Helper function that wraps calls to posres */
void
posres_wrapper(t_nrnb               *nrnb,
               const t_idef         *idef,
               const struct t_pbc   *pbc,
               const rvec           *x,
               gmx_enerdata_t       *enerd,
               const real           *lambda,
               const t_forcerec     *fr,
               gmx::ForceWithVirial *forceWithVirial);

/*! \brief Helper function that wraps calls to posres for free-energy
    pertubation */
void
posres_wrapper_lambda(struct gmx_wallcycle *wcycle,
                      const t_lambda       *fepvals,
                      const t_idef         *idef,
                      const struct t_pbc   *pbc,
                      const rvec            x[],
                      gmx_enerdata_t       *enerd,
                      const real           *lambda,
                      const t_forcerec     *fr);

/*! \brief Helper function that wraps calls to fbposres for
    free-energy perturbation */
void fbposres_wrapper(t_nrnb               *nrnb,
                      const t_idef         *idef,
                      const struct t_pbc   *pbc,
                      const rvec           *x,
                      gmx_enerdata_t       *enerd,
                      const t_forcerec     *fr,
                      gmx::ForceWithVirial *forceWithVirial);

#endif
