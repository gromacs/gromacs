/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * Tests for PlumedOptions and PlumedOptionProvider classes of Plumed MDModule.
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/plumed/plumedOptions.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
namespace gmx
{
namespace
{
TEST(PlumedOptionsTest, defaultConstructor)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    // PlumedOptions should be not active if no file are passed to it
    PlumedOptionProvider options_;
    checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
}

TEST(PlumedOptionsTest, setTimeStep)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    PlumedOptionProvider options_;

    const real myTimeStep = 1.0;
    options_.setSimulationTimeStep(myTimeStep);
    checker.checkValue(options_.options().simulationTimeStep_, "Setting the timestep to 1.0");
    checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
}

TEST(PlumedOptionsTest, setStartingBehavior)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    PlumedOptionProvider            options_;

    const StartingBehavior myStart = StartingBehavior::NewSimulation;
    options_.setStartingBehavior(myStart);
    checker.checkValue(options_.options().startingBehavior_ == StartingBehavior::NewSimulation,
                       "Setting the starting behavior to NewSimulation");
    checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
}


TEST(PlumedOptionsTest, setPlumedFile)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    PlumedOptionProvider            options_;

    const std::optional<std::string> myFile = "myFile";
    options_.setPlumedFile(myFile);

    checker.checkValue(options_.options().plumedFile_, "Setting the plumed file to \"myFile\"");
    // this is the ONLY setter with the extra side effect of setting the active flag
    checker.checkBoolean(options_.active(), "Plumed should be activated when the filename is set");
}

TEST(PlumedOptionsTest, setPlumedFileNotSet)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    PlumedOptionProvider            options_;

    const std::optional<std::string> myFile{};

    options_.setPlumedFile(myFile);

    // checker.checkValue(options_.options().plumedFile_, "Setting the plumed file to \"myFile\"");
    // this is the ONLY setter with the extra side effect of setting the active flag
    checker.checkBoolean(options_.active(),
                         "Plumed should not be activated when the filename is not set");
}

TEST(PlumedOptionsTest, setEnsembleTemperature_data)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    PlumedOptionProvider            options_;
    real                            temp = 300.0;
    t_inputrec                      ir;
    ir.ensembleTemperatureSetting = EnsembleTemperatureSetting::Constant;
    ir.ensembleTemperature        = temp;

    EnsembleTemperature ensembleT(ir);

    options_.setEnsembleTemperature(ensembleT);
    checker.checkValue(options_.options().ensembleTemperature_.value(),
                       "Setting the ensemble temperature to 300.0");
    checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
}


TEST(PlumedOptionsTest, setEnsembleTemperature_notConstant)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    {
        SCOPED_TRACE("using variable temperature");
        PlumedOptionProvider options_;
        real                 temp = 300.0;
        t_inputrec           ir;
        ir.ensembleTemperatureSetting = EnsembleTemperatureSetting::Variable;
        ir.ensembleTemperature        = temp;

        EnsembleTemperature ensembleT(ir);

        options_.setEnsembleTemperature(ensembleT);
        checker.checkValue(
                (options_.options().ensembleTemperature_.has_value()),
                "If the temperature is not constant, its optional value should not be set");
        checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
    }
    {
        SCOPED_TRACE("using not available temperature");
        PlumedOptionProvider options_;
        t_inputrec           ir;
        ir.ensembleTemperatureSetting = EnsembleTemperatureSetting::NotAvailable;

        EnsembleTemperature ensembleT(ir);

        options_.setEnsembleTemperature(ensembleT);

        checker.checkValue(
                (options_.options().ensembleTemperature_.has_value()),
                "If the temperature is not constant, its optional value should not be set");
        checker.checkBoolean(options_.active(), "Plumed should be deactivated by default");
    }
}


void getTopologyFromSimulationDatabase(gmx_mtop_t* const mtop, const std::string& fileName)
{
    gmx::test::TestFileManager  fileManager_;
    const std::filesystem::path simDatabase =
            gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();

    // Generate empty mdp file
    const std::string mdpInputFileName = fileManager_.getTemporaryFilePath(fileName + ".mdp").string();
    gmx::TextWriter::writeFileFromString(mdpInputFileName, "");

    // Generate tpr file
    const std::string tprName = fileManager_.getTemporaryFilePath(fileName + ".tpr").string();
    {
        gmx::test::CommandLine caller;
        caller.append("grompp");
        caller.addOption("-f", mdpInputFileName);
        caller.addOption("-p", (simDatabase / fileName).replace_extension(".top").string());
        caller.addOption("-c", (simDatabase / fileName).replace_extension(".gro").string());
        caller.addOption("-o", tprName);
        ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
    }

    // Load topology
    bool    fullTopology;
    matrix  box;
    PbcType pbcType;
    readConfAndTopology(tprName.c_str(), &fullTopology, mtop, &pbcType, nullptr, nullptr, box);
}

TEST(PlumedOptionsTest, setTopology)
{
    gmx::test::TestReferenceData    data{};
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    for (const auto* const fnm :
         { "4water", "angles1", "argon5832", "argon4", "dipoles", "spc_and_methane" })
    {
        std::string fname(fnm);
        gmx_mtop_t  myTopology{};

        PlumedOptionProvider options_;
        getTopologyFromSimulationDatabase(&myTopology, fname);
        options_.setTopology(myTopology);

        checker.checkBoolean(options_.options().natoms_ == myTopology.natoms,
                             ("Correct number of atoms (" + fname + ")").c_str());
        checker.checkBoolean(options_.active(),
                             ("Plumed should be deactivated by default(" + fname + ")").c_str());
    }
}

} // namespace
} // namespace gmx
