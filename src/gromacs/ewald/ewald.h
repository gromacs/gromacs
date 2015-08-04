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
/*! \libinternal \defgroup module_ewald Ewald-family treatments of long-ranged forces
 * \ingroup group_mdrun
 *
 * \brief Computes energies and forces for long-ranged interactions
 * using the Ewald decomposition. Includes plain Ewald, PME, P3M for
 * Coulomb, PME for Lennard-Jones, load-balancing for PME, and
 * supporting code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Erik Lindahl <erik@kth.se>
 * \author Roland Schulz <roland@rschulz.eu>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Christan Wennberg <cwennberg@kth.se>
 */

/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the plain-Ewald long-ranged part,
 * and the correction for overall system charge for all Ewald-family
 * methods.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_EWALD_H
#define GMX_EWALD_EWALD_H

#include <stdio.h>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

/* Forward declaration of type for managing Ewald tables */
struct gmx_ewald_tab_t;

/*! \brief Initialize the tables used in the Ewald long-ranged part */
void
init_ewald_tab(struct gmx_ewald_tab_t **et, const t_inputrec *ir,
               FILE *fp);

/*! \brief Do the long-ranged part of an Ewald calculation */
real
do_ewald(t_inputrec *ir,
         rvec x[],        rvec f[],
         real chargeA[],  real chargeB[],
         rvec box,
         t_commrec *cr,  int natoms,
         matrix lrvir,   real ewaldcoeff,
         real lambda,    real *dvdlambda,
         struct gmx_ewald_tab_t *et);

/*! \brief Calculate the correction to the Ewald sum, due to a net system
 * charge.
 *
 * Should only be called on one thread. */
real
ewald_charge_correction(t_commrec *cr, t_forcerec *fr, real lambda, matrix box,
                        real *dvdlambda, tensor vir);

#endif
