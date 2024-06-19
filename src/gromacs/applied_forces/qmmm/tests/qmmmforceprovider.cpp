/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Tests for QMMMForceProvider class for QMMM MDModule
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/qmmm/qmmmforceprovider.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/qmmm/qmmmtypes.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

enum class PbcType : int;

namespace gmx
{

class QMMMForceProviderTest : public ::testing::Test
{
public:
    void setDefaultParameters()
    {
        parameters_.active_                = true;
        std::vector<gmx::Index> qmIndicies = { 0, 1, 2 };
        std::vector<gmx::Index> mmIndicies = { 3, 4, 5 };
        LocalAtomSet            set1 = atomSetManager_.add(ArrayRef<const gmx::Index>(qmIndicies));
        LocalAtomSet            set2 = atomSetManager_.add(ArrayRef<const gmx::Index>(mmIndicies));
        qmAtomSet_                   = std::make_unique<LocalAtomSet>(set1);
        mmAtomSet_                   = std::make_unique<LocalAtomSet>(set2);
    }

protected:
    QMMMParameters                parameters_;
    LocalAtomSetManager           atomSetManager_;
    std::unique_ptr<LocalAtomSet> qmAtomSet_;
    std::unique_ptr<LocalAtomSet> mmAtomSet_;
    PbcType                       pbcType_;
    MDLogger                      logger_;
};

TEST_F(QMMMForceProviderTest, CanConstructOrNot)
{
    setDefaultParameters();

    // GMX_CP2K is defined in CMakeList.txt trough set_source_files_properties()
    if (GMX_CP2K)
    {
        // if libcp2k linked then we do not expect throws
        EXPECT_NO_THROW(QMMMForceProvider forceProvider(
                parameters_, *qmAtomSet_, *mmAtomSet_, pbcType_, logger_));
    }
    else
    {
        // if libcp2k not linked then we expect throw from constructor
        EXPECT_ANY_THROW(QMMMForceProvider forceProvider(
                parameters_, *qmAtomSet_, *mmAtomSet_, pbcType_, logger_));
    }
}

} // namespace gmx
