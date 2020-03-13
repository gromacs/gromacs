/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
/*! \defgroup module_listed_forces Interactions between lists of particles
 * \ingroup group_mdrun
 *
 * \brief Handles computing energies and forces for listed
 * interactions.
 *
 * Located here is the code for
 * - computing energies and forces for interactions between a small
     number of particles, e.g bonds, position restraints and listed
     non-bonded interactions (e.g. 1-4).
 * - high-level functions used by mdrun for computing a set of such
     quantities
 * - managing thread-wise decomposition, thread-local buffer output,
     and reduction of output data across threads.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 */
/*! \libinternal \file
 *
 * \brief This file contains declarations of high-level functions used
 * by mdrun to compute energies and forces for listed interactions.
 *
 * Clients of libgromacs that want to evaluate listed interactions
 * should call functions declared here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_listed_forces
 */
#ifndef GMX_LISTED_FORCES_LISTED_FORCES_H
#define GMX_LISTED_FORCES_LISTED_FORCES_H

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_enerdata_t;
struct gmx_grppairener_t;
struct gmx_multisim_t;
class history_t;
class InteractionDefinitions;
struct t_commrec;
struct t_fcdata;
struct t_forcerec;
struct t_lambda;
struct t_mdatoms;
struct t_nrnb;
class t_state;

namespace gmx
{
class ForceOutputs;
class StepWorkload;
template<typename>
class ArrayRef;
} // namespace gmx

//! Type of CPU function to compute a bonded interaction.
using BondedFunction = real (*)(int              nbonds,
                                const t_iatom    iatoms[],
                                const t_iparams  iparams[],
                                const rvec       x[],
                                rvec4            f[],
                                rvec             fshift[],
                                const t_pbc*     pbc,
                                real             lambda,
                                real*            dvdlambda,
                                const t_mdatoms* md,
                                t_fcdata*        fcd,
                                int*             ddgatindex);

//! Getter for finding a callable CPU function to compute an \c ftype interaction.
BondedFunction bondedFunction(int ftype);

/*! \brief Do all aspects of energy and force calculations for mdrun
 * on the set of listed interactions
 *
 * xWholeMolecules only needs to contain whole molecules when orientation
 * restraints need to be computed and can be empty otherwise.
 */
void do_force_listed(struct gmx_wallcycle*          wcycle,
                     const matrix                   box,
                     const t_lambda*                fepvals,
                     const t_commrec*               cr,
                     const gmx_multisim_t*          ms,
                     const InteractionDefinitions&  idef,
                     const rvec                     x[],
                     gmx::ArrayRef<const gmx::RVec> xWholeMolecules,
                     history_t*                     hist,
                     gmx::ForceOutputs*             forceOutputs,
                     const t_forcerec*              fr,
                     const struct t_pbc*            pbc,
                     gmx_enerdata_t*                enerd,
                     t_nrnb*                        nrnb,
                     const real*                    lambda,
                     const t_mdatoms*               md,
                     struct t_fcdata*               fcd,
                     int*                           global_atom_index,
                     const gmx::StepWorkload&       stepWork);

/*! \brief Returns true if there are position, distance or orientation restraints. */
bool haveRestraints(const InteractionDefinitions& idef, const t_fcdata& fcd);

/*! \brief Returns true if there are CPU (i.e. not GPU-offloaded) bonded interactions to compute. */
bool haveCpuBondeds(const t_forcerec& fr);

/*! \brief Returns true if there are listed interactions to compute.
 *
 * NOTE: the current implementation returns true if there are position restraints
 * or any bonded interactions computed on the CPU.
 */
bool haveCpuListedForces(const t_forcerec& fr, const InteractionDefinitions& idef, const t_fcdata& fcd);

#endif
