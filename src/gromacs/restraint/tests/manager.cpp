/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

#include "gmxpre.h"

#include "gromacs/restraint/manager.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmxapi
{
class SessionResources;
} // namespace gmxapi

namespace gmx
{
namespace test
{
namespace
{

class DummyRestraint : public gmx::IRestraintPotential
{
public:
    ~DummyRestraint() override = default;

    gmx::PotentialPointData evaluate(gmx::Vector gmx_unused r1,
                                     gmx::Vector gmx_unused r2,
                                     double gmx_unused      t) override
    {
        return {};
    }

    void update(gmx::Vector gmx_unused v, gmx::Vector gmx_unused v0, double gmx_unused t) override
    {
    }

    std::vector<int> sites() const override { return std::vector<int>(); }

    void bindSession(gmxapi::SessionResources* session) override { (void)session; }
};

TEST(RestraintManager, restraintList)
{
    auto managerInstance = gmx::RestraintManager();
    managerInstance.addToSpec(std::make_shared<DummyRestraint>(), "a");
    managerInstance.addToSpec(std::make_shared<DummyRestraint>(), "b");
    EXPECT_EQ(managerInstance.countRestraints(), 2);
    managerInstance.clear();
    EXPECT_EQ(managerInstance.countRestraints(), 0);
    managerInstance.addToSpec(std::make_shared<DummyRestraint>(), "c");
    managerInstance.addToSpec(std::make_shared<DummyRestraint>(), "d");
    EXPECT_EQ(managerInstance.countRestraints(), 2);
}

} // end namespace
} // namespace test
} // namespace gmx
