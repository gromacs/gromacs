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
#include "gmxpre.h"

#include "gromacs/fmm/fmm_mdpoptions.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fmm/fmm_mdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreemdpwriter.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

/**
 * @brief Convenience method to create FmmMdpOptions from transformed MDP key-value input.
 *
 * Mimics grompp-style option processing for testing.
 *
 * @param moduleMdpValues Key-value tree representing MDP input values.
 * @return Filled FmmMdpOptions object.
 */
static FmmMdpOptions fillOptionsFromMdpValues(const KeyValueTreeObject& moduleMdpValues)
{
    return test::fillOptionsFromMdpValuesTemplate<FmmMdpOptions>(moduleMdpValues);
}

/*! \brief Returns the string name for a given FMM backend.
 *
 * \param backend The FMM backend enum value.
 * \return The corresponding backend name as a std::string.
 */
static std::string fmmBackendName(ActiveFmmBackend backend)
{
    return std::string(c_activeFmmBackendNames[backend]);
}

class FmmMdpOptionsTest : public ::testing::Test
{
protected:
    FmmMdpOptions fmmMdpOptions_;
};

TEST_F(FmmMdpOptionsTest, ThrowsOnExaFmmOptionsAccessWhenInactive)
{
    EXPECT_THROW(fmmMdpOptions_.exaFmmOptions(), gmx::InternalError);
}

TEST_F(FmmMdpOptionsTest, ThrowsOnFMSolvrOptionsAccessWhenInactive)
{
    EXPECT_THROW(fmmMdpOptions_.fmSolvrOptions(), gmx::InternalError);
}

TEST_F(FmmMdpOptionsTest, ExaFmmOptionsMatchDefaultsWhenActive)
{
    // Activate ExaFMM backend via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    const ExaFmmOptions  defaultExaFmmOptionValues;
    const ExaFmmOptions& opts = fmmMdpOptions_.exaFmmOptions();

    EXPECT_EQ(opts.order, defaultExaFmmOptionValues.order);
    EXPECT_EQ(opts.directRange, defaultExaFmmOptionValues.directRange);
    EXPECT_EQ(opts.directProvider, defaultExaFmmOptionValues.directProvider);
}

TEST_F(FmmMdpOptionsTest, FMSolvrOptionsMatchDefaultsWhenActive)
{
    // Activate FMSolvr backend via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::FMSolvr));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    const FMSolvrOptions  defaultFMSolvrOptionValues;
    const FMSolvrOptions& opts = fmmMdpOptions_.fmSolvrOptions();

    EXPECT_EQ(opts.order, defaultFMSolvrOptionValues.order);
    EXPECT_EQ(opts.directRange, defaultFMSolvrOptionValues.directRange);
    EXPECT_EQ(opts.directProvider, defaultFMSolvrOptionValues.directProvider);
    EXPECT_EQ(opts.treeDepth, defaultFMSolvrOptionValues.treeDepth);
    EXPECT_EQ(opts.dipoleCompensation, defaultFMSolvrOptionValues.dipoleCompensation);
    EXPECT_EQ(opts.sparse, defaultFMSolvrOptionValues.sparse);
}

TEST_F(FmmMdpOptionsTest, ParsesNonDefaultExaFmmOrder)
{
    // Activate ExaFMM backend with non-default order via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string activeBackendKey;
    activeBackendKey.append(FmmModuleInfo::sc_name);
    activeBackendKey.append("-");
    activeBackendKey.append(c_fmmActiveOptionName);

    std::string orderKey;
    orderKey.append(FmmModuleInfo::sc_name);
    orderKey.append("-");
    orderKey.append(c_fmmExaFmmOrderOptionName);

    mdpValueBuilder.rootObject().addValue(activeBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    mdpValueBuilder.rootObject().addValue(orderKey, std::string("10"));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    const ExaFmmOptions& opts = fmmMdpOptions_.exaFmmOptions();
    EXPECT_EQ(opts.order, 10); // Non-default value
}

TEST_F(FmmMdpOptionsTest, ParsesNonDefaultFMSolvrOrder)
{
    // Activate FMSolvr backend with non-default order via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string activeBackendKey;
    activeBackendKey.append(FmmModuleInfo::sc_name);
    activeBackendKey.append("-");
    activeBackendKey.append(c_fmmActiveOptionName);

    std::string orderKey;
    orderKey.append(FmmModuleInfo::sc_name);
    orderKey.append("-");
    orderKey.append(c_fmmFMSolvrOrderOptionName);

    mdpValueBuilder.rootObject().addValue(activeBackendKey, fmmBackendName(ActiveFmmBackend::FMSolvr));
    mdpValueBuilder.rootObject().addValue(orderKey, std::string("10"));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    const FMSolvrOptions& opts = fmmMdpOptions_.fmSolvrOptions();
    EXPECT_EQ(opts.order, 10); // Non-default value
}


TEST_F(FmmMdpOptionsTest, OutputNoDefaultValuesWhenInactive)
{
    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    fmmMdpOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(FmmMdpOptionsTest, OutputDefaultExaFmmValuesWhenActive)
{
    // Activate ExaFMM backend via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string activeBackendKey;
    activeBackendKey.append(FmmModuleInfo::sc_name);
    activeBackendKey.append("-");
    activeBackendKey.append(c_fmmActiveOptionName);

    mdpValueBuilder.rootObject().addValue(activeBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    fmmMdpOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

TEST_F(FmmMdpOptionsTest, OutputDefaultFMSolvrValuesWhenActive)
{
    // Activate FMSolvr backend via MDP configuration
    KeyValueTreeBuilder mdpValueBuilder;

    std::string activeBackendKey;
    activeBackendKey.append(FmmModuleInfo::sc_name);
    activeBackendKey.append("-");
    activeBackendKey.append(c_fmmActiveOptionName);

    mdpValueBuilder.rootObject().addValue(activeBackendKey, fmmBackendName(ActiveFmmBackend::FMSolvr));

    fmmMdpOptions_ = fillOptionsFromMdpValues(mdpValueBuilder.build());

    // Transform module data into a flat key-value tree for output.
    StringOutputStream        stream;
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder builderObject = builder.rootObject();

    fmmMdpOptions_.buildMdpOutput(&builderObject);
    {
        TextWriter writer(&stream);
        writeKeyValueTreeAsMdp(&writer, builder.build());
    }
    stream.close();

    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());

    checker.checkString(stream.toString(), "Mdp output");
}

} // namespace test
} // namespace gmx
