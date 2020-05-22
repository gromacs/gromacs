/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Declares interface to LINCS code.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_LINCS_H
#define GMX_MDLIB_LINCS_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct gmx_multisim_t;
class InteractionDefinitions;
struct t_commrec;
struct t_inputrec;
struct t_nrnb;
struct t_pbc;

namespace gmx
{

template<typename>
class ArrayRefWithPadding;
enum class ConstraintVariable : int;
class Lincs;
template<typename>
class ListOfLists;


/*! \brief Return the data for determining constraint RMS relative deviations. */
ArrayRef<real> lincs_rmsdData(Lincs* lincsd);

/*! \brief Return the RMSD of the constraint. */
real lincs_rmsd(const Lincs* lincsd);

/*! \brief Initializes and returns the lincs data struct. */
Lincs* init_lincs(FILE*                            fplog,
                  const gmx_mtop_t&                mtop,
                  int                              nflexcon_global,
                  ArrayRef<const ListOfLists<int>> atomsToConstraintsPerMolType,
                  bool                             bPLINCS,
                  int                              nIter,
                  int                              nProjOrder);

/*! \brief Destructs the lincs object when it is not nullptr. */
void done_lincs(Lincs* li);

/*! \brief Initialize lincs stuff */
void set_lincs(const InteractionDefinitions& idef,
               int                           numAtoms,
               const real*                   invmass,
               real                          lambda,
               bool                          bDynamics,
               const t_commrec*              cr,
               Lincs*                        li);

/*! \brief Applies LINCS constraints.
 *
 * \returns true if the constraining succeeded. */
bool constrain_lincs(bool                            computeRmsd,
                     const t_inputrec&               ir,
                     int64_t                         step,
                     Lincs*                          lincsd,
                     const real*                     invmass,
                     const t_commrec*                cr,
                     const gmx_multisim_t*           ms,
                     ArrayRefWithPadding<const RVec> x,
                     ArrayRefWithPadding<RVec>       xprime,
                     ArrayRef<RVec>                  min_proj,
                     const matrix                    box,
                     t_pbc*                          pbc,
                     bool                            hasMassPerturbed,
                     real                            lambda,
                     real*                           dvdlambda,
                     real                            invdt,
                     ArrayRef<RVec>                  v,
                     bool                            bCalcVir,
                     tensor                          vir_r_m_dr,
                     ConstraintVariable              econq,
                     t_nrnb*                         nrnb,
                     int                             maxwarn,
                     int*                            warncount);

} // namespace gmx

#endif
