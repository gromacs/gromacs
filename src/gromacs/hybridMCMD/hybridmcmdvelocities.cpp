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

/*! \internal \file
 * \brief
 * Implements the HybridMCMDVelocities class and an interface to a Metropolis step.
 *
 * \author Sebastian Wingbermuehle
 * \ingroup module_hybridMCMD
 */

#include "gmxpre.h"

#include "hybridmcmdvelocities.h"

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

HybridMCMDVelocities::HybridMCMDVelocities(IMetropolisStepVelocities *metropolisStep) : metropolisStep_(metropolisStep)
{
}

//! \brief Draws initial velocities from a Maxwell-Boltzmann distribution for the short NVE simulation between two Metropolis steps.
void HybridMCMDVelocities::draw(t_inputrec                *ir,
                                int64_t                    step,
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
                                gmx_bool                  *bSumEkinhOld)
{
    /* use Andersen-massive thermostat (do_per_step should activate it in step 0)
     * as there is no thermostat with hybrid MC/MD, we don't need to pass gmx_update_t
     */
    bool bIfRandomize;
    bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state, nullptr, constr);
    /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
    if (bIfRandomize)
    {
        bool bStopCM = (ir->comm_mode != ecmNO);
        if (constr)
        {
            constrain_velocities(step, nullptr,
                                 state,
                                 tmp_vir,
                                 wcycle, constr,
                                 bCalcVir, do_log, do_ene);
        }
        /* calculate kinetic energy (as done for RE in md.cpp)
         * remove centre-of-mass motion if required
         * if centre-of-mass motion has to be removed, call compute_globals twice (first remove centre-of-mass motion, then calculate Ekin)
         * always pass CGLO_GSTAT to ensure communication
         */
        // Is CGLO_SCALEEKIN needed here?
        compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                        wcycle, enerd, nullptr, nullptr, nullptr, nullptr, mu_tot,
                        constr, nullSignaller, state->box,
                        nullptr, bSumEkinhOld,
                        CGLO_GSTAT | (!bStopCM ? CGLO_TEMPERATURE : 0) | (bStopCM ? CGLO_STOPCM : 0) | CGLO_SCALEEKIN);
        if (bStopCM)
        {
            compute_globals(fplog, gstat, cr, ir, fr, ekind, state, mdatoms, nrnb, vcm,
                            wcycle, enerd, nullptr, nullptr, nullptr, nullptr, mu_tot,
                            constr, nullSignaller, state->box,
                            nullptr, bSumEkinhOld,
                            CGLO_GSTAT | CGLO_TEMPERATURE | CGLO_SCALEEKIN);
        }

        // pass kinetic energy to Metropolis criterion (will also be printed to the log file at the end of the step)
        metropolisStep_->setInitialKineticEnergy(enerd->term[F_EKIN]);
    }
}

} // namespace
