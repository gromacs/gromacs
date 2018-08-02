/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

/*! \libinternal
 * \ingroup group_mdrun
 * \brief
 * Implements a part of a hybrid Monte Carlo integrator for md-vv.
 *
 * \author Sebastian Wingbermuehle
 */

/*! \libinternal \file
 *
 * \brief
 * Declares the HybridMCMDVelocities class and an interface to a Metropolis step.
 *
 * This class draws the initial velocities of the short NVE simulations from a
 * Boltzmann distribution using the Andersen-massive routines. Moreoever, it
 * calculates the corresponding kinetic energy (after removing the centre-of-mass
 * motion if desired) and reports it to the Metropolis step.
 *
 * \author Sebastian Wingbermuehle
 * \inlibraryapi
 * \ingroup module_hybridMCMD
 */

#ifndef GMX_HYBRID_MC_MD_VELOCITIES_H
#define GMX_HYBRID_MC_MD_VELOCITIES_H

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct gmx_global_stat;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
class t_state;

namespace gmx
{
class Constraints;
class SimulationSignaller;

class IMetropolisStepVelocities
{
    public:
        //! \brief Set the initial kinetic energy to be used in the Metropolis criterion
        virtual void setInitialKineticEnergy(const double initialKineticEnergy) = 0;

        //! \brief Destructor
        virtual ~IMetropolisStepVelocities() {};
};

class HybridMCMDVelocities
{
    public:
        //! \brief Constructor
        HybridMCMDVelocities(IMetropolisStepVelocities *metropolisStep);
        void draw(t_inputrec                *ir,
                  gmx_int64_t                step,
                  t_commrec                 *cr,
                  t_mdatoms                 *mdatoms,
                  t_state                   *state,
                  Constraints               *constr,
                  tensor                     tmp_vir,
                  gmx_wallcycle             *wcycle,
                  gmx_bool                   bCalcVir,
                  bool                       do_log,
                  bool                       do_ene,
                  FILE                      *fplog,
                  gmx_global_stat           *gstat,
                  t_forcerec                *fr,
                  gmx_ekindata_t            *ekind,
                  t_nrnb                    *nrnb,
                  t_vcm                     *vcm,
                  gmx_enerdata_t            *enerd,
                  rvec                       mu_tot,
                  SimulationSignaller       *nullSignaller,
                  gmx_bool                  *bSumEkinhOld);

    private:
        IMetropolisStepVelocities *metropolisStep_;
};

} // namespace

#endif
