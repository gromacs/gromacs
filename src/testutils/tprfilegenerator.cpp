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
 * Implements helper for generating reusuable TPR files for tests within the same test binary.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_testutils
 */

#include "gmxpre.h"

#include "tprfilegenerator.h"

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"

namespace gmx
{
namespace test
{

TprAndFileManager::TprAndFileManager(const std::string& name)
{
    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath(name + ".mdp");
    gmx::TextWriter::writeFileFromString(mdpInputFileName, "");
    tprFileName_ = fileManager_.getTemporaryFilePath(name + ".tpr");
    {
        CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", TestFileManager::getInputFilePath(name + ".top"));
        caller.addOption("-c", TestFileManager::getInputFilePath(name + ".pdb"));
        caller.addOption("-o", tprFileName_);
        EXPECT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }
}

} // namespace test
} // namespace gmx
