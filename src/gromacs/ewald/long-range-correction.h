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
/*! \libinternal \file
 *
 * \brief This file contains function declarations necessary for
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Erik Lindahl <erik@kth.se>
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_LONG_RANGE_CORRECTION_H
#define GMX_EWALD_LONG_RANGE_CORRECTION_H

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/forcerec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

/*! \brief Calculate long-range Ewald correction terms.
 *
 * For the group cutoff scheme (only), calculates the correction to
 * the Ewald sums (electrostatic and/or LJ) due to pairs excluded from
 * the long-ranged part.
 *
 * For both cutoff schemes, but only for Coulomb interactions,
 * calculates correction for surface dipole terms. */
void
ewald_LRcorrection(int start, int end,
                   t_commrec *cr, int thread, t_forcerec *fr,
                   real *chargeA, real *chargeB,
                   real *C6A, real *C6B,
                   real *sigmaA, real *sigmaB,
                   real *sigma3A, real *sigma3B,
                   gmx_bool bHaveChargeOrTypePerturbed,
                   gmx_bool calc_excl_corr,
                   t_blocka *excl, rvec x[],
                   matrix box, rvec mu_tot[],
                   int ewald_geometry, real epsilon_surface,
                   rvec *f, tensor vir_q, tensor vir_lj,
                   real *Vcorr_q, real *Vcorr_lj,
                   real lambda_q, real lambda_lj,
                   real *dvdlambda_q, real *dvdlambda_lj);

#endif
