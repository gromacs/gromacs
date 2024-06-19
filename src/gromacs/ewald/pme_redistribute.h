/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief Declares functions for redistributing atoms over the PME domains
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_REDISTRIBUTE_H
#define GMX_EWALD_PME_REDISTRIBUTE_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "pme_internal.h"

namespace gmx
{
template<typename T>
class ArrayRef;
} // namespace gmx
struct t_commrec;

//! Redistributes forces along the dimension gives by \p atc
void dd_pmeredist_f(struct gmx_pme_t* pme, PmeAtomComm* atc, gmx::ArrayRef<gmx::RVec> f, gmx_bool bAddF);

//! Redistributes coefficients and when \p bFirst=true coordinates over MPI ranks
void do_redist_pos_coeffs(struct gmx_pme_t*              pme,
                          const t_commrec*               cr,
                          gmx_bool                       bFirst,
                          gmx::ArrayRef<const gmx::RVec> x,
                          gmx::ArrayRef<const real>      data);

#endif
