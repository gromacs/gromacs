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
#include "testutils/testoptions.h"
#include "testutils/argsbuilder.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/init.h"
#include "programs/gmx/legacymodules.h"
#include "programs/mdrun/mdrun_main.h"

#include <sstream>

extern "C"
{

// TODO should there be a header file for these?
// src/programs/gmx/legacymodules.cpp needs them also
/*! \cond */
int gmx_gmxcheck(int argc, char *argv[]);
int gmx_gmxdump(int argc, char *argv[]);
int gmx_grompp(int argc, char *argv[]);
int gmx_pdb2gmx(int argc, char *argv[]);
int gmx_protonate(int argc, char *argv[]);
int gmx_tpbconv(int argc, char *argv[]);
int gmx_x2top(int argc, char *argv[]);
/*! \endcond */

}

namespace gmx
{
namespace test
{

/********************************************************************
 * MdrunTestFixture
 */

MdrunTestFixture::MdrunTestFixture()
    : numThreads(1),       //TODO consider how varying this should be managed
      numOpenMPThreads(1), //TODO consider how varying this should be managed
      // TODO: link default behaviours in mdrun test fixtures to those
      // for production mdrun
      topFileName(),
      groFileName(),
      trrFileName(),
      mdpInputFileName(fileManager_.getTemporaryFilePath("input.mdp")),
      mdpOutputFileName(fileManager_.getTemporaryFilePath("output.mdp")),
      tprFileName(fileManager_.getTemporaryFilePath(".tpr")),
      logFileName(fileManager_.getTemporaryFilePath(".log")),
      edrFileName(fileManager_.getTemporaryFilePath(".edr"))
{
}

MdrunTestFixture::~MdrunTestFixture()
{
}

void
MdrunTestFixture::useEmptyMdpFile()
{
    useStringAsMdpFile("");
}

void
MdrunTestFixture::useStringAsMdpFile(const char *mdpString)
{
    gmx::File::writeFileFromString(mdpInputFileName, mdpString);
}

void
MdrunTestFixture::useTopAndGroFromDatabase(const char *name)
{
    topFileName = fileManager_.getInputFilePath((std::string(name) + ".top").c_str());
    groFileName = fileManager_.getInputFilePath((std::string(name) + ".gro").c_str());
    if (0 == strcmp(name, "spc2"))
    {
        ASSERT_EQ(numThreads, 1) << "The toy system is so small only one thread can be used";
        ASSERT_EQ(numOpenMPThreads, 1) << "The toy system is so small only one OpenMP thread can be used";
    }
}

int
MdrunTestFixture::callGrompp()
{
    gmx::test::ProgramCaller caller("grompp");
    caller.addOption(FileNameOption("f").description("Input .mdp file").store(&mdpInputFileName));
    caller.addOption(FileNameOption("p").description("Input .top file").store(&topFileName));
    caller.addOption(FileNameOption("c").description("Input structure file").store(&groFileName));

    caller.addOption(FileNameOption("po").description("Output post-processed .mdp file").store(&mdpOutputFileName));
    caller.addOption(FileNameOption("o").description("Output .tpr file").store(&tprFileName));

    bool bQuiet = true;
    caller.addOption(BooleanOption("quiet").store(&bQuiet));

    ArgsBuilder args(caller.getProgramName());
    args.visitSubSection(caller.getOptions());
    int         rc = gmx::CommandLineModuleManager::runAsMainCMain(args.getArgc(), args.getArgv(), &gmx_grompp, false);

    return rc;
}

int
MdrunTestFixture::callMdrun()
{
    gmx::test::ProgramCaller caller(std::string("mdrun"));

    caller.addOption(FileNameOption("s").description("Input .tpr file").store(&tprFileName));
    caller.addOption(FileNameOption("rerun").description("Input trajectory file for rerun").store(&rerunFileName));

    caller.addOption(FileNameOption("g").description("Output .log file").store(&logFileName));
    caller.addOption(FileNameOption("e").description("Output .edr file").store(&edrFileName));
    caller.addOption(FileNameOption("o").description("Output trajectory file").store(&trrFileName));
    caller.addOption(FileNameOption("x").description("Output .xtc file").store(&xtcFileName));

#ifdef GMX_THREAD_MPI
    caller.addOption(IntegerOption("nt").description("Number of tMPI threads to use").store(&numThreads));
#endif
#ifdef GMX_LIB_MPI
    caller.addOption(IntegerOption("np").description("Number of MPI ranks to use").store(&numThreads));
#endif
#ifdef GMX_OPENMP
    caller.addOption(IntegerOption("ntomp").description("Number of OpenMP threads to use").store(&numOpenMPThreads));
#endif
    bool bQuiet = true;
    caller.addOption(BooleanOption("quiet").store(&bQuiet));

    ArgsBuilder args(caller.getProgramName());
    args.visitSubSection(caller.getOptions());
    int         rc = gmx::CommandLineModuleManager::runAsMainCMain(args.getArgc(), args.getArgv(), &gmx_mdrun, false);

    return rc;
}

} // namespace test
} // namespace gmx
