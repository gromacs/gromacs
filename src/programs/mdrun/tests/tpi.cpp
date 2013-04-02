/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Tests for the mdrun test-particle insertion functionality
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include <gtest/gtest.h>
#include "../legacyheaders/macros.h"
#include "moduletest.h"

#ifdef __cplusplus
extern "C" {
#endif

int grompp_cmain(int argc, char *argv[]);
int mdrun_cmain(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

namespace
{

//! Test fixture for tpi integrator (which uses the re-run machinery)
typedef gmx::test::MdrunTestFixture TpiTest;

TEST_F(TpiTest, TpiExitsNormally)
{
    // We need named string variables so that later c_str() can return
    // allocated memory.

    // Input files
    std::string mdpFileName = fileManager.getInputFilePath("simple.mdp");
    std::string topFileName = fileManager.getInputFilePath("spc2.top");
    std::string groFileName = fileManager.getInputFilePath("spc2.gro");
    std::string trrFileName = fileManager.getInputFilePath("spc2.trr");

    // Output files
    std::string mdpOutputFileName = fileManager.getTemporaryFilePath(".mdp");
    std::string tprFileName = fileManager.getTemporaryFilePath(".tpr");
    std::string logFileName = fileManager.getTemporaryFilePath(".log");
    std::string edrFileName = fileManager.getTemporaryFilePath(".edr");
    std::string numThreadsStr = std::to_string(numThreads);
    std::string numOpenMPThreadsStr = std::to_string(numOpenMPThreads);

    ASSERT_EQ(numThreads, 1) << "The toy system is so small only one thread can be used";
    ASSERT_EQ(numOpenMPThreads, 1) << "The toy system is so small only one OpenMP thread can be used";

    // Don't fake argv[0] if you need at run time to use the machinery
    // that finds share/top in the source tree.
    static char *grompp_argv[] = {
        (char*) programName.c_str(),
        (char*) "-f",  (char*) mdpFileName.c_str(),
        (char*) "-p",  (char*) topFileName.c_str(),
        (char*) "-c",  (char*) groFileName.c_str(),
        (char*) "-o",  (char*) tprFileName.c_str(),
        (char*) "-po", (char*) mdpOutputFileName.c_str(),
    };

    // This could be dressed up in the test fixture?
    EXPECT_EQ(0,grompp_cmain(asize(grompp_argv),grompp_argv));

    static char *mdrun_argv[] = {
        (char*) programName.c_str(),
        (char*) "-s",     (char*) tprFileName.c_str(),
        (char*) "-g",     (char*) logFileName.c_str(),
        (char*) "-e",     (char*) edrFileName.c_str(),
        (char*) "-rerun", (char*) trrFileName.c_str(),
        (char*) "-nt",    (char*) numThreadsStr.c_str(),
        (char*) "-ntomp", (char*) numOpenMPThreadsStr.c_str(),
    };

    ASSERT_EQ(0, mdrun_cmain(asize(mdrun_argv), mdrun_argv));
}

} // namespace
