/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Tests for ColvarsOptions class of Colvars MDModule.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/colvars/colvarsoptions.h"

#include "config.h"

#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/colvars/colvarsMDModule.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

enum class PbcType : int;

namespace gmx
{

static const std::string colvarsConfig = "colvars_sample.dat";

//! Convenience method to avoid specifying the template parameter repetitively
static ColvarsOptions fillOptionsFromMdpValues(const KeyValueTreeObject& moduleMdpValues)
{
    return test::fillOptionsFromMdpValuesTemplate<ColvarsOptions>(moduleMdpValues);
}

class ColvarsOptionsTest : public ::testing::Test
{
public:
    static KeyValueTreeObject ColvarsBuildDefaulMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(ColvarsModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        return mdpValueBuilder.build();
    }

    static KeyValueTreeObject ColvarsBuildInputMdpValues()
    {
        // Prepare MDP inputs
        KeyValueTreeBuilder mdpValueBuilder;
        mdpValueBuilder.rootObject().addValue(std::string(ColvarsModuleInfo::sc_name) + "-active",
                                              std::string("true"));
        mdpValueBuilder.rootObject().addValue(
                std::string(ColvarsModuleInfo::sc_name) + "-configfile", colvarsConfig);
        mdpValueBuilder.rootObject().addValue(std::string(ColvarsModuleInfo::sc_name) + "-seed",
                                              std::string("12345789"));
        return mdpValueBuilder.build();
    }

#if GMX_HAVE_COLVARS

    void prepareInputColvarsPreProcessor(const std::string& fileName, ColvarsOptions* colvarsOptions)
    {

        // Path to the sample colvars input file
        std::string colvarsConfigFile =
                gmx::test::TestFileManager::getInputFilePath("colvars_sample.dat").string();

        gmx::test::TestFileManager  fileManager_;
        const std::filesystem::path simData =
                gmx::test::TestFileManager::getTestSimulationDatabaseDirectory();

        // Generate empty mdp file
        const std::string mdpInputFileName =
                fileManager_.getTemporaryFilePath(fileName + ".mdp").string();
        gmx::TextWriter::writeFileFromString(mdpInputFileName, "");

        // Generate tpr file
        const std::string tprName = fileManager_.getTemporaryFilePath(fileName + ".tpr").string();
        {
            gmx::test::CommandLine caller;
            caller.append("grompp");
            caller.addOption("-f", mdpInputFileName);
            caller.addOption("-p", (simData / fileName).replace_extension(".top").string());
            caller.addOption("-c", (simData / fileName).replace_extension(".gro").string());
            caller.addOption("-o", tprName);
            ASSERT_EQ(0, gmx_grompp(caller.argc(), caller.argv()));
        }


        bool       fullTopology;
        PbcType    pbcType;
        matrix     box;
        gmx_mtop_t mtop;

        // Load topology
        readConfAndTopology(tprName.c_str(), &fullTopology, &mtop, &pbcType, &coords, nullptr, box);

        atoms = gmx_mtop_global_atoms(mtop);
        ArrayRef<const RVec> x =
                gmx::constArrayRefFromArray(reinterpret_cast<gmx::RVec*>(coords), atoms.nr);

        // Populate attributes outside the use of the defined callbacks.
        colvarsOptions->setParameters(colvarsConfigFile, atoms, x, pbcType, box, 300);
    }

    void deleteInputColvarsPreProcessor()
    {
        sfree(coords);
        done_atom(&atoms);
    }

#endif // GMX_HAVE_COLVARS

protected:
    rvec*   coords;
    t_atoms atoms;
};


// The following tests should work either with or without Colvars compiled
// because GROMACS without colvars should still be able to read a tpr file with colvars mdp keywords.
// (For example, gmx dump should be able to output colvars values even without Colvars compiled.)
TEST_F(ColvarsOptionsTest, OutputNoDefaultValuesWhenInactive)
{
    // Test buildMdpOutput()
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    ColvarsOptions colvarsOptions;
    colvarsOptions.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(ColvarsOptionsTest, OutputDefaultValuesWhenActive)
{

    // Activate colvars
    ColvarsOptions colvarsOptions = fillOptionsFromMdpValues(ColvarsBuildDefaulMdpValues());

    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    colvarsOptions.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(ColvarsOptionsTest, OutputValuesWhenActive)
{

    // Activate colvars
    ColvarsOptions colvarsOptions = fillOptionsFromMdpValues(ColvarsBuildInputMdpValues());

    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    colvarsOptions.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

#if GMX_HAVE_COLVARS

TEST_F(ColvarsOptionsTest, OptionSetsActive)
{
    ColvarsOptions colvarsOptions;
    EXPECT_FALSE(colvarsOptions.isActive());
    colvarsOptions = fillOptionsFromMdpValues(ColvarsBuildDefaulMdpValues());
    EXPECT_TRUE(colvarsOptions.isActive());
}

TEST_F(ColvarsOptionsTest, InternalsToKvtAndBack)
{

    // Activate colvars
    ColvarsOptions colvarsOptions = fillOptionsFromMdpValues(ColvarsBuildInputMdpValues());
    // Set up parameters with a test system
    prepareInputColvarsPreProcessor("4water", &colvarsOptions);

    // Write parameters to the KVT
    KeyValueTreeBuilder builder;
    MDLogger            logger;
    colvarsOptions.setLogger(logger);
    colvarsOptions.writeInternalParametersToKvt(builder.rootObject());
    const auto inputTree = builder.build();

    // Copy internal parameters
    auto refColvarsInputContent = colvarsOptions.colvarsConfigContent();
    auto refColvarsCoordinates  = colvarsOptions.colvarsAtomCoords();
    auto refTemperature         = colvarsOptions.colvarsEnsTemp();
    auto refSeed                = colvarsOptions.colvarsSeed();

    // Retrieve paramaters from the KVT
    colvarsOptions.readInternalParametersFromKvt(inputTree);

    // Check parameters taken back from KVT
    EXPECT_EQ(refColvarsInputContent, colvarsOptions.colvarsConfigContent());

    auto actualColvarsCoordinates = colvarsOptions.colvarsAtomCoords();
    EXPECT_REAL_EQ(refColvarsCoordinates[0][XX], actualColvarsCoordinates[0][XX]);
    EXPECT_REAL_EQ(refColvarsCoordinates[0][YY], actualColvarsCoordinates[0][YY]);
    EXPECT_REAL_EQ(refColvarsCoordinates[0][ZZ], actualColvarsCoordinates[0][ZZ]);
    EXPECT_REAL_EQ(refColvarsCoordinates[1][XX], actualColvarsCoordinates[1][XX]);
    EXPECT_REAL_EQ(refColvarsCoordinates[1][YY], actualColvarsCoordinates[1][YY]);
    EXPECT_REAL_EQ(refColvarsCoordinates[1][ZZ], actualColvarsCoordinates[1][ZZ]);

    EXPECT_EQ(refTemperature, colvarsOptions.colvarsEnsTemp());
    EXPECT_EQ(refSeed, colvarsOptions.colvarsSeed());

    deleteInputColvarsPreProcessor();
}


TEST_F(ColvarsOptionsTest, RetrieveEdrFilename)
{
    // Activate colvars
    ColvarsOptions colvarsOptions = fillOptionsFromMdpValues(ColvarsBuildInputMdpValues());

    std::string refEdrFilename = "output/ener.edr";
    colvarsOptions.processEdrFilename(EdrOutputFilename{ refEdrFilename });
    const std::string ref = std::filesystem::path("output/ener").make_preferred().string();
    EXPECT_EQ(ref, colvarsOptions.colvarsOutputPrefix());

    refEdrFilename = "sim.part1.edr";
    colvarsOptions.processEdrFilename(EdrOutputFilename{ refEdrFilename });
    EXPECT_EQ("sim.part1", colvarsOptions.colvarsOutputPrefix());
}
#endif // GMX_HAVE_COLVARS

} // namespace gmx
