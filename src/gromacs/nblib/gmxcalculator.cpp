/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Implements a force calculator based on GROMACS data structures.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include "gmxcalculator.h"

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/nblib/simulationstate.h"
#include "gromacs/nbnxm/nbnxm.h"

#include "nbkerneloptions.h"

namespace nblib
{

static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

interaction_const_t setupInteractionConst(const std::shared_ptr<NBKernelOptions> options)
{
    interaction_const_t interactionConst;
    interactionConst.vdwtype      = evdwCUT;
    interactionConst.vdw_modifier = eintmodPOTSHIFT;
    interactionConst.rvdw         = options->pairlistCutoff;

    switch (options->coulombType)
    {
        case BenchMarkCoulomb::Pme: interactionConst.eeltype = eelPME; break;
        case BenchMarkCoulomb::Cutoff: interactionConst.eeltype = eelCUT; break;
        case BenchMarkCoulomb::ReactionField: interactionConst.eeltype = eelRF; break;
        case BenchMarkCoulomb::Count:
            GMX_THROW(gmx::InvalidInputError("Unsupported electrostatic interaction"));
    }
    interactionConst.coulomb_modifier = eintmodPOTSHIFT;
    interactionConst.rcoulomb         = options->pairlistCutoff;
    // Note: values correspond to ic.coulomb_modifier = eintmodPOTSHIFT
    interactionConst.dispersion_shift.cpot = -1.0 / gmx::power6(interactionConst.rvdw);
    interactionConst.repulsion_shift.cpot  = -1.0 / gmx::power12(interactionConst.rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    interactionConst.epsilon_r  = 1.0;
    interactionConst.epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (interactionConst.epsilon_r != 0)
    {
        interactionConst.epsfac = ONE_4PI_EPS0 / interactionConst.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        interactionConst.epsfac = 0;
    }

    calc_rffac(nullptr, interactionConst.epsilon_r, interactionConst.epsilon_rf,
               interactionConst.rcoulomb, &interactionConst.k_rf, &interactionConst.c_rf);

    if (EEL_PME_EWALD(interactionConst.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        interactionConst.ewaldcoeff_q = ewaldCoeff(1e-5, options->pairlistCutoff);
        GMX_RELEASE_ASSERT(interactionConst.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        interactionConst.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &interactionConst);
    }
    return interactionConst;
}

gmx::StepWorkload setupStepWorkload(const std::shared_ptr<NBKernelOptions> options)
{
    gmx::StepWorkload stepWork;
    stepWork.computeForces          = true;
    stepWork.computeNonbondedForces = true;

    if (options->computeVirialAndEnergy)
    {
        stepWork.computeVirial = true;
        stepWork.computeEnergy = true;
    }
    return stepWork;
}

GmxForceCalculator::GmxForceCalculator(SimulationState                        simState,
                                       const std::shared_ptr<NBKernelOptions> options) :
    interactionConst_(setupInteractionConst(options)),
    stepWork_(setupStepWorkload(options)),
    enerd_(gmx_enerdata_t(1, 0))
{
    gmx::fillLegacyMatrix(simState.box().matrix(), box_);
}

gmx::PaddedHostVector<gmx::RVec> GmxForceCalculator::compute()
{

    nbv_->dispatchNonbondedKernel(gmx::InteractionLocality::Local, interactionConst_, stepWork_,
                                  enbvClearFNo, forcerec_, &enerd_, &nrnb_);

    nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::All, verletForces_);

    return verletForces_;
}

} // namespace nblib
