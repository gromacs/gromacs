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
 * \brief
 * Tests for functionality of the density fitting module.
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/densityfitting.h"

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

TEST(DensityFittingTest, Options)
{
    auto densityFittingModule(gmx::createDensityFittingModule());

    // Prepare MDP inputs
    gmx::KeyValueTreeBuilder mdpValueBuilder;
    mdpValueBuilder.rootObject().addValue("density-guided-simulation-active", "yes");
    KeyValueTreeObject       densityFittingMdpValues = mdpValueBuilder.build();

    // set up options
    gmx::Options densityFittingModuleOptions;
    densityFittingModule->mdpOptionProvider()->initMdpOptions(&densityFittingModuleOptions);

    // Add rules to transform mdp inputs to densityFittingModule data
    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
    densityFittingModule->mdpOptionProvider()->initMdpTransform(transform.rules());

    // Execute the transform on the mdpValues
    auto transformedMdpValues = transform.transform(densityFittingMdpValues, nullptr);
    gmx::assignOptionsFromKeyValueTree(&densityFittingModuleOptions, transformedMdpValues.object(), nullptr);

    // Build the force provider, once all input data is gathered
    ForceProviders densityFittingForces;
    densityFittingModule->initForceProviders(&densityFittingForces);

    // Build a minimal simulation system.
    t_mdatoms               mdAtoms;
    mdAtoms.homenr = 1;
    PaddedVector<gmx::RVec> x             = {{0, 0, 0}};
    matrix                  simulationBox = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    const double            t             = 0.0;
    t_commrec              *cr            = init_commrec();
    gmx::ForceProviderInput forceProviderInput(x, mdAtoms, t, simulationBox, *cr);

    // The forces that the force-provider is to update
    PaddedVector<gmx::RVec>  f = {{0, 0, 0}};
    gmx::ForceWithVirial     forceWithVirial(f, false);

    gmx_enerdata_t           energyData(1, 0);
    gmx::ForceProviderOutput forceProviderOutput(&forceWithVirial, &energyData);

    // update the forces
    densityFittingForces.calculateForces(forceProviderInput, &forceProviderOutput);

    // check that calculated forces match expected forces
    std::vector<RVec> expectedForce = {{0, 0, 0}};
    EXPECT_THAT(expectedForce, Pointwise(test::RVecEq(test::defaultFloatTolerance()),
                                         forceProviderOutput.forceWithVirial_.force_));

    // clean up C-style commrec, so this test leaks no memory
    // \todo remove once commrec manages its own memory
    done_commrec(cr);
}

} // namespace

} // namespace gmx
