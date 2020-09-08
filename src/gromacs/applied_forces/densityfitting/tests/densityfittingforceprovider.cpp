/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * Tests amplitude lookup for density fitting
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/densityfitting/densityfittingforceprovider.h"

#include <gtest/gtest.h>

namespace gmx
{

TEST(DensityFittingForceProviderState, RoundTripSaving)
{
    DensityFittingForceProviderState state;
    // set-up state
    state.adaptiveForceConstantScale_                   = 1.0;
    state.stepsSinceLastCalculation_                    = 0;
    state.exponentialMovingAverageState_.increasing_    = false;
    state.exponentialMovingAverageState_.weightedCount_ = 0;
    state.exponentialMovingAverageState_.weightedSum_   = 0;

    KeyValueTreeBuilder kvtBuilder;
    const std::string   identifier = "test-module";
    state.writeState(kvtBuilder.rootObject(), identifier);
    KeyValueTreeObject stateStoredInKvt = kvtBuilder.build();

    // invalidate state
    state.adaptiveForceConstantScale_                   = -1;
    state.stepsSinceLastCalculation_                    = -1;
    state.exponentialMovingAverageState_.increasing_    = true;
    state.exponentialMovingAverageState_.weightedCount_ = -1;
    state.exponentialMovingAverageState_.weightedSum_   = -1;

    // read back the original state
    state.readState(stateStoredInKvt, identifier);

    EXPECT_EQ(state.adaptiveForceConstantScale_, 1.0);
    EXPECT_EQ(state.stepsSinceLastCalculation_, 0);

    EXPECT_EQ(state.exponentialMovingAverageState_.increasing_, false);
    EXPECT_EQ(state.exponentialMovingAverageState_.weightedCount_, 0);
    EXPECT_EQ(state.exponentialMovingAverageState_.weightedSum_, 0);
}


} // namespace gmx
