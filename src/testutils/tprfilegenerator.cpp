/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Implements helper for generating reusuable TPR files for tests within the same test binary.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_testutils
 */

#include "gmxpre.h"

#include "testutils/tprfilegenerator.h"

#include <filesystem>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

TprAndFileManager::TprAndFileManager(const std::string& baseName, const std::string& mdpContent) :
    tprFileName_{ fileManager_.getTemporaryFilePath(baseName + ftp2ext_with_dot(efTPR)) },
    mdpInputFileName_{ fileManager_.getTemporaryFilePath(baseName + ftp2ext_with_dot(efMDP)) },
    inputConf_{ baseName + ftp2ext_with_dot(efPDB) },
    topology_{ baseName + ftp2ext_with_dot(efTOP) }
{
    TextWriter::writeFileFromString(mdpInputFileName_, mdpContent);
}

void TprAndFileManager::setInputConf(const std::filesystem::path& name)
{
    if (hasGeneratedTpr_)
    {
        GMX_THROW(APIError("Cannot set file name after a TPR file has been generated"));
    }
    inputConf_ = name;
}

void TprAndFileManager::setPositionRestraintConf(const std::filesystem::path& name)
{
    if (hasGeneratedTpr_)
    {
        GMX_THROW(APIError("Cannot set file name after a TPR file has been generated"));
    }
    positionRestraints_ = name;
}

void TprAndFileManager::setTopology(const std::filesystem::path& name)
{
    if (hasGeneratedTpr_)
    {
        GMX_THROW(APIError("Cannot set file name after a TPR file has been generated"));
    }
    topology_ = name;
}

const std::filesystem::path& TprAndFileManager::tprName() const
{
    if (!hasGeneratedTpr_)
    {
        CommandLine caller;
        caller.append("grompp");

        caller.addOption("-f", mdpInputFileName_);
        caller.addOption("-p", fileManager_.getInputFilePath(topology_));
        caller.addOption("-c", fileManager_.getInputFilePath(inputConf_));
        if (!positionRestraints_.empty())
        {
            caller.addOption("-r", fileManager_.getInputFilePath(positionRestraints_));
        }
        caller.addOption("-o", tprFileName_);
        // Ensure mdout.mdp file is cleaned up
        caller.addOption("-po", fileManager_.getTemporaryFilePath("mdout.mdp"));

        EXPECT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        hasGeneratedTpr_ = true;
    }

    return tprFileName_;
}

} // namespace test
} // namespace gmx
