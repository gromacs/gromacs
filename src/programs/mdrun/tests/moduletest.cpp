/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Implements classes in moduletest.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include "moduletest.h"

#include "testutils/integrationtests.h"
#include "testutils/integrationtests_impl.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "testutils/testoptions.h"
#include "gromacs/utility/exceptions.h"

#ifdef __cplusplus
extern "C" {
#endif

int grompp_cmain(int argc, char *argv[]);
int mdrun_cmain(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

namespace gmx
{
namespace test
{
 
/********************************************************************
 * MdrunTestFixture
 */

MdrunTestFixture::MdrunTestFixture()
    : numThreads(0),
      numOpenMPThreads(0),
      // TODO: link default behaviours in mdrun test fixtures to those
      // for production mdrun
      mdpFileName(fileManager.getTemporaryFilePath("input.mdp")),
      topFileName(),
      groFileName(),
      trrFileName(),
      mdpOutputFileName(fileManager.getTemporaryFilePath("output.mdp")),
      tprFileName(fileManager.getTemporaryFilePath(".tpr")),
      logFileName(fileManager.getTemporaryFilePath(".log")),
      edrFileName(fileManager.getTemporaryFilePath(".edr"))
{
    // Set up the handling of the command-line options to mdrun-test
    // (these are not the options for the GROMACS tools it will call!)
    gmx::Options options("mdrun-test", "mdrun test standard options");
    options.addOption(IntegerOption("nt").description("Number of tMPI threads to use").store(&numThreads));
    options.addOption(IntegerOption("ntomp").description("Number of OpenMP threads to use").store(&numOpenMPThreads));
    parseTestOptions(&options);

    configureProgramWithDefaultOptions("grompp");
    configureProgramWithDefaultOptions("mdrun");
}

MdrunTestFixture::~MdrunTestFixture()
{
}

void
MdrunTestFixture::configureProgramWithDefaultOptions(std::string const& programName)
{
    if (0 == programName.compare("grompp"))
    {
        impl_->createProgramCaller(programName, grompp_cmain);

        impl_->defineProgramOption<FileNameOption>(programName, "f", "Input .mdp file");
        impl_->defineProgramOption<FileNameOption>(programName, "p", "Input .top file");
        impl_->defineProgramOption<FileNameOption>(programName, "c", "Input structure file");
        impl_->defineProgramOption<FileNameOption>(programName, "po", "Output post-processed .mdp file");
        impl_->defineProgramOption<FileNameOption>(programName, "o", "Output .tpr file");
    }
    else if (0 == programName.compare("mdrun"))
    {
        impl_->createProgramCaller(programName, mdrun_cmain);
        impl_->defineProgramOption<FileNameOption>(programName, "s", "Input .tpr file");
        impl_->defineProgramOption<IntegerOption>(programName, "nt", "Number of tMPI threads to use");
        impl_->defineProgramOption<IntegerOption>(programName, "ntomp", "Number of OpenMP threads to use");
        impl_->defineProgramOption<FileNameOption>(programName, "rerun", "Input .trr file for rerun");
        impl_->defineProgramOption<FileNameOption>(programName, "g", "Output .log file");
        impl_->defineProgramOption<FileNameOption>(programName, "e", "Output .edr file");
        impl_->defineProgramOption<FileNameOption>(programName, "tpi", "Output .xvg file for TPI");
        impl_->defineProgramOption<FileNameOption>(programName, "tpid", "Output .xvg file for TPI distances");
        impl_->defineProgramOption<FileNameOption>(programName, "t", "Output .trr file");
    }
    else
    {
        GMX_THROW(APIError("In MdrunTestFixture::configureProgramWithDefaultOptions() an unknown program name was specified"));
    }
}

int
MdrunTestFixture::callGrompp()
{
    std::string programName("grompp");

    setProgramOption(programName, "f", mdpFileName);
    setProgramOption(programName, "p", topFileName);
    setProgramOption(programName, "c", groFileName);
    setProgramOption(programName, "po", mdpOutputFileName);
    setProgramOption(programName, "o", tprFileName);

    return runProgram(programName);
}

int
MdrunTestFixture::callMdrun()
{
    std::string programName("mdrun");

    setProgramOption(programName, "s", tprFileName);
    setProgramOption(programName, "g", logFileName);
    setProgramOption(programName, "e", edrFileName);
    setProgramOption(programName, "nt", numThreads);
    setProgramOption(programName, "ntomp", numOpenMPThreads);

    return runProgram(programName);
}

} // namespace test
} // namespace gmx
