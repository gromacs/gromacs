/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Helper class declarations for testing gmx_ana tools.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 */

#include "gmxpre.h"

#include <initializer_list>
#include <memory>

#include "testutils/cmdlinetest.h"

namespace gmx
{

template<typename>
class ArrayRef;

namespace test
{

//! \libinternal \brief Test fixture for tests that call a gmx analysis tool with input/output files.
class GmxAnaTestBase : public CommandLineTestBase
{
public:
    GmxAnaTestBase();
    ~GmxAnaTestBase() override;

    //! \brief Run command and check the resulting output.
    void runAndCheckResults();

    //! \brief Select groups when prompted for when running the tool.
    void selectGroups(const std::initializer_list<const char*> groups);

    //! \brief Select groups when prompted for when running the tool.
    void selectGroups(const ArrayRef<const std::string> groups);

private:
    //! \brief Function for tool to test, must be set by derived class.
    virtual int gmxTool(int, char**) const = 0;

    class Impl;

    //! Handle to implementation object
    std::unique_ptr<Impl> impl_;
};

} // namespace test
} // namespace gmx
