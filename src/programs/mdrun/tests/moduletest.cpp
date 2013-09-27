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

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include "testutils/integrationtests.h"
#include "testutils/testoptions.h"
#include "testutils/argsbuilder.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/init.h"
#include "gromacs/legacyheaders/network.h"
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

int MdrunTestFixture::numRanks = 1;
int MdrunTestFixture::numThreads = 1;
int MdrunTestFixture::numOpenMPThreads = 1;

GMX_TEST_OPTIONS(MdrunTestOptions, options)
{
#ifdef GMX_THREAD_MPI
    options->addOption(IntegerOption("nt").store(&MdrunTestFixture::numThreads)
                       .description("Number of thread-MPI threads/ranks for child mdrun call"));
#endif
#ifdef GMX_OPENMP
    options->addOption(IntegerOption("nt_omp").store(&MdrunTestFixture::numOpenMPThreads)
                       .description("Number of OpenMP threads for child mdrun call"));
#endif
}

MdrunTestFixture::MdrunTestFixture() :
    topFileName(),
    groFileName(),
    trrFileName(),
    mdpInputFileName(fileManager_.getTemporaryFilePath("input.mdp")),
    mdpOutputFileName(fileManager_.getTemporaryFilePath("output.mdp")),
    tprFileName(fileManager_.getTemporaryFilePath(".tpr")),
    logFileName(fileManager_.getTemporaryFilePath(".log")),
    edrFileName(fileManager_.getTemporaryFilePath(".edr"))
{
#ifdef GMX_LIB_MPI
    assert(gmx_mpi_initialized());
#endif
    numRanks = gmx_node_num();
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
}

int
MdrunTestFixture::callGrompp()
{
#ifdef GMX_LIB_MPI
    // When compiled with external MPI, only call one instance of the grompp function
    if (0 != gmx_node_rank())
    {
        return 0;
    }
#endif

    gmx::test::ProgramCaller caller("grompp");
    caller.addOption(FileNameOption("f").description("Input .mdp file").store(&mdpInputFileName));
    caller.addOption(FileNameOption("p").description("Input .top file").store(&topFileName));
    caller.addOption(FileNameOption("c").description("Input structure file").store(&groFileName));

    caller.addOption(FileNameOption("po").description("Output post-processed .mdp file").store(&mdpOutputFileName));
    caller.addOption(FileNameOption("o").description("Output .tpr file").store(&tprFileName));

    ArgsBuilder args(caller.getProgramName());
    args.visitSubSection(caller.getOptions());
    int         rc = gmx_grompp(args.getArgc(), args.getArgv());

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
#ifdef GMX_OPENMP
    caller.addOption(IntegerOption("ntomp").description("Number of OpenMP threads to use").store(&numOpenMPThreads));
#endif

    ArgsBuilder args(caller.getProgramName());
    args.visitSubSection(caller.getOptions());
    int         rc = gmx_mdrun(args.getArgc(), args.getArgv());

    return rc;
}

} // namespace test
} // namespace gmx
