/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief
 * Implements the HybridMCMDVelocities class and an interface to a Metropolis step.
 *
 * \author Sebastian Wingbermuehle
 * \ingroup module_hybridMCMD
 */

#include "gmxpre.h"

#include "hybridmcmdvelocities.h"
#include "metropolisinterfaces.h"

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/timing/wallcycle.h"

namespace gmx
{

HybridMCMDVelocities::HybridMCMDVelocities(IMetropolisStepVelocities &metropolisStep) : metropolisStep_(metropolisStep)
{
}

void HybridMCMDVelocities::draw(t_inputrec                *ir,
                                int64_t                    step,
                                t_commrec                 *cr,
                                t_mdatoms                 *mdatoms,
                                t_state                   *state,
                                Constraints               *constr,
                                tensor                     tmp_vir,
                                gmx_wallcycle             *wcycle,
                                gmx_bool                   bCalcVir,
                                gmx_bool                   do_log,
                                gmx_bool                   do_ene,
                                FILE                      *fplog,
                                gmx_global_stat           *gstat,
                                gmx_ekindata_t            *ekind,
                                t_nrnb                    *nrnb,
                                t_vcm                     *vcm,
                                gmx_enerdata_t            *enerd,
                                SimulationSignaller       *nullSignaller)
{
    if (do_per_step(step, ir->hybridMCMDParams->nstMetropolis))
    {

        /* use functions for Andersen-massive thermostat to generate random initial velocities
         * we don't need to pass:
         * - real rate (TODO pass parameter rate to allow for partial randomization of velocities needed for a generalized hybrid Monte Carlo)
         * - const gmx_bool *randomize (there is no thermostat with hybrid MC/MD)
         * - const real *boltzfac (there is no thermostat with hybrid MC/MD)
         */
        std::vector<bool> nullVectorBool;
        std::vector<real> nullVectorReal;
        andersen_tcoupl(ir, step, cr, mdatoms, state->v, 0, nullVectorBool, nullVectorReal);

        /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
        if (constr)
        {
            constrain_velocities(step, nullptr,
                                 state,
                                 tmp_vir,
                                 constr,
                                 bCalcVir, do_log, do_ene);
        }

        /* For md-vv, we just need the full-step kinetic energy here. There is no need to interfere with the logic of bSumEkinhOld. */
        bool bSumEkinhOldDummy = false;

        /* remove centre-of-mass motion if required */
        if (ir->comm_mode != ecmNO)
        {
            /* enerd is not changed in this call, but global_stat always copies some of its entries to a buffer (that may not be used later)
             * pass CGLO_GSTAT to ensure communication (flags allow to reduce (MPI) only the desired quantities)
             */
            compute_globals(fplog, gstat, cr, ir, nullptr, nullptr, state, mdatoms, nrnb, vcm,
                            wcycle, enerd, nullptr, nullptr, nullptr, nullptr, nullptr,
                            nullptr, nullSignaller, nullptr, nullptr, &bSumEkinhOldDummy,
                            CGLO_GSTAT | CGLO_STOPCM);
        }
        /* calculate kinetic energy */
        compute_globals(nullptr, gstat, cr, ir, nullptr, ekind, state, mdatoms, nrnb, nullptr,
                        wcycle, enerd, nullptr, nullptr, nullptr, nullptr, nullptr,
                        nullptr, nullSignaller, nullptr, nullptr, &bSumEkinhOldDummy,
                        CGLO_GSTAT | CGLO_TEMPERATURE);

        /* pass kinetic energy to Metropolis criterion (will also be printed to the log file at the end of the step) */
        metropolisStep_.setInitialKineticEnergy(enerd->term[F_EKIN]);
    }
}

} // namespace
