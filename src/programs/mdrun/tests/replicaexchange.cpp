/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
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
 * Tests for the mdrun replica-exchange functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <gtest/gtest.h>

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

#include "multisimtest.h"

namespace gmx
{
namespace test
{

//! Convenience typedef
typedef MultiSimTest ReplicaExchangeEnsembleTest;

TEST_P(ReplicaExchangeEnsembleTest, ExitsNormally)
{
    mdrunCaller_->addOption("-replex", 1);
    runExitsNormallyTest();
}

/* Note, not all preprocessor implementations nest macro expansions
   the same way / at all, if we would try to duplicate less code. */
#if GMX_LIB_MPI
INSTANTIATE_TEST_CASE_P(WithDifferentControlVariables,
                        ReplicaExchangeEnsembleTest,
                        ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_WithDifferentControlVariables,
                        ReplicaExchangeEnsembleTest,
                        ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#endif

//! Convenience typedef
typedef MultiSimTest ReplicaExchangeTerminationTest;

TEST_F(ReplicaExchangeTerminationTest, WritesCheckpointAfterMaxhTerminationAndThenRestarts)
{
    mdrunCaller_->addOption("-replex", 1);
    runMaxhTest();
}

} // namespace test
} // namespace gmx
