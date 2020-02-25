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

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/rf_util.h"
#include "gromacs/gmxlib/nrnb.h"

#include "gmxcalculator.h"
#include "nbkerneloptions.h"
#include "simulationstate.h"

namespace nblib
{

static real ewaldCoeff(const real ewald_rtol, const real pairlistCutoff)
{
    return calc_ewaldcoeff_q(pairlistCutoff, ewald_rtol);
}

//! Return an interaction constants struct with members used in the benchmark set appropriately
static interaction_const_t setupInteractionConst(const std::shared_ptr<NBKernelOptions> options)
{
    interaction_const_t ic;

    ic.vdwtype      = evdwCUT;
    ic.vdw_modifier = eintmodPOTSHIFT;
    ic.rvdw         = options->pairlistCutoff;

    switch (options->coulombType)
    {
        case BenchMarkCoulomb::Pme: ic.eeltype = eelPME; break;
        case BenchMarkCoulomb::Cutoff: ic.eeltype = eelCUT; break;
        case BenchMarkCoulomb::ReactionField: ic.eeltype = eelRF; break;
        case BenchMarkCoulomb::Count:
            GMX_THROW(gmx::InvalidInputError("Unsupported electrostatic interaction"));
    }
    ic.coulomb_modifier = eintmodPOTSHIFT;
    ic.rcoulomb         = options->pairlistCutoff;
    //! Note: values correspond to ic.coulomb_modifier = eintmodPOTSHIFT
    ic.dispersion_shift.cpot = -1.0 / gmx::power6(ic.rvdw);
    ic.repulsion_shift.cpot  = -1.0 / gmx::power12(ic.rvdw);

    // These are the initialized values but we leave them here so that later
    // these can become options.
    ic.epsilon_r  = 1.0;
    ic.epsilon_rf = 1.0;

    /* Set the Coulomb energy conversion factor */
    if (ic.epsilon_r != 0)
    {
        ic.epsfac = ONE_4PI_EPS0 / ic.epsilon_r;
    }
    else
    {
        /* eps = 0 is infinite dieletric: no Coulomb interactions */
        ic.epsfac = 0;
    }

    calc_rffac(nullptr, ic.epsilon_r, ic.epsilon_rf, ic.rcoulomb, &ic.k_rf, &ic.c_rf);

    if (EEL_PME_EWALD(ic.eeltype))
    {
        // Ewald coefficients, we ignore the potential shift
        ic.ewaldcoeff_q = ewaldCoeff(1e-5, options->pairlistCutoff);
        GMX_RELEASE_ASSERT(ic.ewaldcoeff_q > 0, "Ewald coefficient should be > 0");
        ic.coulombEwaldTables = std::make_unique<EwaldCorrectionTables>();
        init_interaction_const_tables(nullptr, &ic);
    }

    return ic;
}

GmxForceCalculator::GmxForceCalculator(const std::shared_ptr<SimulationState> system,
                                       const std::shared_ptr<NBKernelOptions> options) :
    enerd_(1, 0),
    verletForces_({})
{
    interactionConst_ = setupInteractionConst(options);

    gmx::fillLegacyMatrix(system->box().matrix(), box_);

    stepWork_.computeForces          = true;
    stepWork_.computeNonbondedForces = true;

    if (options->computeVirialAndEnergy)
    {
        stepWork_.computeVirial = true;
        stepWork_.computeEnergy = true;
    }

    forcerec_.ntype = system->topology().numParticles();
}

gmx::PaddedHostVector<gmx::RVec> GmxForceCalculator::compute()
{
    t_nrnb nrnb = { 0 };

    nbv_->dispatchNonbondedKernel(gmx::InteractionLocality::Local, interactionConst_, stepWork_,
                                  enbvClearFNo, forcerec_, &enerd_, &nrnb);

    nbv_->atomdata_add_nbat_f_to_f(gmx::AtomLocality::All, verletForces_);

    return verletForces_;
}

} // namespace nblib