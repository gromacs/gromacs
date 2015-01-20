/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include "multisimtest.h"

namespace gmx
{
namespace test
{

//! Convenience typedef
typedef MultiSimTest ReplicaExchangeTest;

/* This test ensures mdrun can run NVT REMD under the supported
 * conditions. It runs one replica per MPI rank.
 *
 * See also comments about MultiSimTest */
TEST_P(ReplicaExchangeTest, ExitsNormally)
{
    if (size_ <= 1)
    {
        /* Can't test replica exchange without multiple ranks. */
        return;
    }

    const char *pcoupl = GetParam();
    organizeMdpFile(pcoupl);
    /* Call grompp on every rank - the standard callGrompp() only runs
       grompp on rank 0. */
    EXPECT_EQ(0, runner_.callGromppOnThisRank());

    // mdrun names the files without the rank suffix
    runner_.tprFileName_ = mdrunTprFileName_;
    mdrunCaller_->addOption("-replex", 1);
    ASSERT_EQ(0, runner_.callMdrun(*mdrunCaller_));
}

/* Note, not all preprocessor implementations nest macro expansions
   the same way / at all, if we would try to duplicate less code. */
#ifdef GMX_LIB_MPI
INSTANTIATE_TEST_CASE_P(WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#else
INSTANTIATE_TEST_CASE_P(DISABLED_WithDifferentControlVariables, ReplicaExchangeTest,
                            ::testing::Values("pcoupl = no", "pcoupl = Berendsen"));
#endif

} // namespace
} // namespace
