/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares interface to SHAKE code.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 * \inlibraryapi
 */

#ifndef GMX_MDLIB_SHAKE_H
#define GMX_MDLIB_SHAKE_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/real.h"

struct InteractionList;
class InteractionDefinitions;
struct t_inputrec;
struct t_nrnb;
struct t_pbc;

namespace gmx
{
template<typename T>
class ArrayRef;

enum class ConstraintVariable : int;

/*! \libinternal
 * \brief Working data for the SHAKE algorithm
 */
struct shakedata
{
    //! Returns the number of SHAKE blocks */
    int numShakeBlocks() const { return sblock.size() - 1; }

    //! The reference constraint vectors
    std::vector<RVec> rij;
    //! The reduced mass of the two atoms in each constraint times 0.5
    std::vector<real> half_of_reduced_mass;
    //! Multiplicative tolerance on the difference in the square of the constrained distance
    std::vector<real> distance_squared_tolerance;
    //! The reference constraint distances squared
    std::vector<real> constraint_distance_squared;
    /* SOR stuff */
    //! SOR delta
    real delta = 0.1;
    //! SOR omega
    real omega = 1.0;
    //! SOR gamma
    real gamma = 1000000;
    //! The SHAKE blocks, block i contains constraints sblock[i]/3 to sblock[i+1]/3 */
    std::vector<int> sblock = { 0 };
    /*! \brief Scaled Lagrange multiplier for each constraint.
     *
     * Value is -2 * eta from p. 336 of the paper, divided by the
     * constraint distance. */
    std::vector<real> scaled_lagrange_multiplier;
};

//! Make SHAKE blocks when not using DD.
void make_shake_sblock_serial(shakedata* shaked, InteractionDefinitions* idef, int numAtoms);

//! Make SHAKE blocks when using DD.
void make_shake_sblock_dd(shakedata* shaked, const InteractionList& ilcon);

/*! \brief Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return TRUE when OK, FALSE when shake-error
 */
bool constrain_shake(FILE*                         log,     /* Log file			*/
                     shakedata*                    shaked,  /* Total number of atoms	*/
                     gmx::ArrayRef<const real>     invmass, /* Atomic masses		*/
                     const InteractionDefinitions& idef,    /* The interaction def		*/
                     const t_inputrec&             ir,      /* Input record		        */
                     ArrayRef<const RVec>          x_s,     /* Coords before update		*/
                     ArrayRef<RVec>                xprime,  /* Output coords when constraining x */
                     ArrayRef<RVec>                vprime,  /* Output coords when constraining v */
                     const t_pbc*                  pbc,     /* PBC information              */
                     t_nrnb*                       nrnb,    /* Performance measure          */
                     real                          lambda,  /* FEP lambda                   */
                     real*                         dvdlambda,  /* FEP force                    */
                     real                          invdt,      /* 1/delta_t                    */
                     ArrayRef<RVec>                v,          /* Also constrain v if not empty  */
                     bool                          bCalcVir,   /* Calculate r x m delta_r      */
                     tensor                        vir_r_m_dr, /* sum r x m delta_r            */
                     bool                          bDumpOnError, /* Dump debugging stuff on error*/
                     ConstraintVariable econq); /* which type of constraint is occurring */

/*! \brief Regular iterative shake */
void cshake(const int            iatom[],
            int                  ncon,
            int*                 nnit,
            int                  maxnit,
            ArrayRef<const real> dist2,
            ArrayRef<RVec>       xp,
            const t_pbc*         pbc,
            ArrayRef<const RVec> rij,
            ArrayRef<const real> half_of_reduced_mass,
            real                 omega,
            ArrayRef<const real> invmass,
            ArrayRef<const real> distance_squared_tolerance,
            ArrayRef<real>       scaled_lagrange_multiplier,
            int*                 nerror);

} // namespace gmx

#endif
