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

void GmxForceCalculator::setupInteractionConst(const std::shared_ptr<NBKernelOptions> options)
{
    interactionConst_.vdwtype      = evdwCUT;
    interactionConst_.vdw_modifier = eintmodPOTSHIFT;
    interactionConst_.rvdw         = options->pairlistCutoff;

    switch (options->coulombType)
    {
        case BenchMarkCoulomb::Pme: interactionConst_.eeltype = eelPME; break;
        case BenchMarkCoulomb::Cutoff: interactionConst_.eeltype = eelCUT; break;
        case BenchMarkCoulomb::ReactionField: interactionConst_.eeltype = eelRF; break;
        case BenchMarkCoulomb::Count:
            GMX_THROW(gmx::InvalidInputError("Unsupported electrostatic interaction"));
    }
    interactionConst_.coulomb_modifier = eintmodPOTSHIFT;
    interactionConst_.rcoulomb         = options->pairlistCutoff;
    //! Note: values correspond to ic.coulomb_modifier = eintmodPOTSHIFT
    interactionConst_.dispersion_shift.cpot = -1.0 / gmx::power6(interactionConst_.rvdw);
    interactionConst_.repulsion_shift.cpot  = -1.0 / gmx::power12(interactionConst_.rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    interactionConst_.epsilon_r  = 1.0;
    interactionConst_.epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (interactionConst_.epsilon_r != 0)
    {
        interactionConst_.epsfac = ONE_4PI_EPS0 / interactionConst_.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        interactionConst_.epsfac = 0;
    }

    calc_rffac(nullptr, interactionConst_.epsilon_r, interactionConst_.epsilon_rf,
               interactionConst_.rcoulomb, &interactionConst_.k_rf, &interactionConst_.c_rf);

    if (EEL_PME_EWALD(interactionConst_.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        interactionConst_.ewaldcoeff_q = ewaldCoeff(1e-5, options->pairlistCutoff);
        GMX_RELEASE_ASSERT(interactionConst_.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        interactionConst_.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &interactionConst_);
    }
}

void GmxForceCalculator::setupStepWorkload(const std::shared_ptr<NBKernelOptions> options)
{
    stepWork_.computeForces          = true;
    stepWork_.computeNonbondedForces = true;

    if (options->computeVirialAndEnergy)
    {
        stepWork_.computeVirial = true;
        stepWork_.computeEnergy = true;
    }
}

GmxForceCalculator::GmxForceCalculator(const std::shared_ptr<SimulationState> system,
                                       const std::shared_ptr<NBKernelOptions> options) :
    verletForces_({}),
    enerd_(1, 0)
{
    setupInteractionConst(options);

    gmx::fillLegacyMatrix(system->box().matrix(), box_);

    setupStepWorkload(options);
}

gmx::PaddedHostVector<gmx::RVec> GmxForceCalculator::compute()
{

    nbv_->dispatchNonbondedKernel(gmx::InteractionLocality::Local, interactionConst_, stepWork_,
                                  enbvClearFNo, forcerec_, &enerd_, &nrnb_);

    nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::All, verletForces_);

    return verletForces_;
}

} // namespace nblib
