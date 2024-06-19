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
/*! \libinternal \file
 *
 * \brief This file declares functions for "pair" interactions
 * (i.e. listed non-bonded interactions, e.g. 1-4 interactions)
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_PAIRS_H
#define GMX_LISTED_FORCES_PAIRS_H

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/real.h"

struct gmx_grppairener_t;
struct t_forcerec;
struct t_pbc;
union t_iparams;

namespace gmx
{
class StepWorkload;
template<typename>
class ArrayRef;
} // namespace gmx

/*! \brief Calculate VdW/charge listed pair interactions (usually 1-4
 * interactions).
 *
 * global_atom_index is only passed for printing error messages.
 */
void do_pairs(int                                 ftype,
              int                                 nbonds,
              const t_iatom                       iatoms[],
              const t_iparams                     iparams[],
              const rvec                          x[],
              rvec4                               f[],
              rvec                                fshift[],
              const struct t_pbc*                 pbc,
              const real*                         lambda,
              real*                               dvdl,
              gmx::ArrayRef<const real>           chargeA,
              gmx::ArrayRef<const real>           chargeB,
              gmx::ArrayRef<const bool>           atomIsPerturbed,
              gmx::ArrayRef<const unsigned short> cENER,
              int                                 numEnergyGroups,
              const t_forcerec*                   fr,
              bool                                havePerturbedPairs,
              const gmx::StepWorkload&            stepWork,
              gmx_grppairener_t*                  grppener,
              int*                                global_atom_index);

#endif
