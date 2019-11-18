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

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct t_commrec;
struct t_forcerec;
struct t_inputrec;

/*! \brief Calculate long-range Ewald correction terms.
 *
 * Calculate correction for electrostatic surface dipole terms.
 */
void ewald_LRcorrection(int               numAtomsLocal,
                        const t_commrec*  cr,
                        int               numThreads,
                        int               thread,
                        const t_forcerec& fr,
                        const t_inputrec& ir,
                        const real*       chargeA,
                        const real*       chargeB,
                        gmx_bool          bHaveChargePerturbed,
                        const rvec        x[],
                        const matrix      box,
                        const rvec        mu_tot[],
                        rvec*             f,
                        real*             Vcorr_q,
                        real              lambda_q,
                        real*             dvdlambda_q);

#endif
