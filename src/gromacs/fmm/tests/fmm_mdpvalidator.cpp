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
 * \brief Tests for FmmMdpValidator.
 *
 * \author Muhammad Umair Sadiq <mumairsadiq1@gmail.com>
 */

#include "gmxpre.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/warninp.h"
#include "gromacs/fmm/fmm_mdmodule.h"
#include "gromacs/fmm/fmm_mdpoptions.h"
#include "gromacs/fmm/fmmoptions.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "testutils/testasserts.h"

namespace gmx
{

namespace test
{

class FmmMdpValidatorTest : public ::testing::Test
{
public:
    FmmMdpValidatorTest() : warningHandler_(false, 0) {}

    void notifyFmmModulePreProcessingSetup()
    {
        fmmModule_->subscribeToPreProcessingNotifications(&notifiers_);
        notifiers_.preProcessingNotifier_.notify(&warningHandler_);
        notifiers_.preProcessingNotifier_.notify(mdCoulombType_);
    }

protected:
    std::unique_ptr<IMDModule> fmmModule_;
    MDModulesNotifiers         notifiers_;

    // default settings
    MdModulesCoulombTypeInfo mdCoulombType_;
    WarningHandler           warningHandler_;
};

TEST_F(FmmMdpValidatorTest, ReportsErrorWhenExaFmmExpansionOrderIsZero)
{
    mdCoulombType_ = { CoulombInteractionType::Fmm };

    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    std::string orderKey;
    orderKey.append(FmmModuleInfo::sc_name);
    orderKey.append("-");
    orderKey.append(c_fmmExaFmmOrderOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    mdpValueBuilder.rootObject().addValue(orderKey, std::string("0"));

    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModulePreProcessingSetup();

    EXPECT_EQ(warningHandler_.errorCount(), 1);
}

TEST_F(FmmMdpValidatorTest, ReportsErrorWithFMSolvrNegativeTreeDepth)
{
    mdCoulombType_ = { CoulombInteractionType::Fmm };

    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    std::string treeDepthKey;
    treeDepthKey.append(FmmModuleInfo::sc_name);
    treeDepthKey.append("-");
    treeDepthKey.append(c_fmmFMSolvrTreeDepthOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::FMSolvr));
    mdpValueBuilder.rootObject().addValue(treeDepthKey, std::string("-1"));

    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModulePreProcessingSetup();

    EXPECT_EQ(warningHandler_.errorCount(), 1);
}

TEST_F(FmmMdpValidatorTest, ReportsErrorOnInvalidExaFmmDirectRangeForGmx)
{
    mdCoulombType_ = { CoulombInteractionType::Fmm };

    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    std::string directProviderKey;
    directProviderKey.append(FmmModuleInfo::sc_name);
    directProviderKey.append("-");
    directProviderKey.append(c_fmmExaFmmDirectProviderOptionName);

    std::string directRangeKey;
    directRangeKey.append(FmmModuleInfo::sc_name);
    directRangeKey.append("-");
    directRangeKey.append(c_fmmExaFmmDirectRangeOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    mdpValueBuilder.rootObject().addValue(directProviderKey,
                                          fmmDirectProviderName(FmmDirectProvider::Gromacs));
    mdpValueBuilder.rootObject().addValue(directRangeKey, std::string("1"));

    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModulePreProcessingSetup();

    // If DirectProvider is GROMACS, DirectRange must be 2 for ExaFMM
    EXPECT_EQ(warningHandler_.errorCount(), 1);
}

TEST_F(FmmMdpValidatorTest, ReportsNoErrorsOnValidMdpSettings)
{
    mdCoulombType_ = { CoulombInteractionType::Fmm };

    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    std::string directProviderKey;
    directProviderKey.append(FmmModuleInfo::sc_name);
    directProviderKey.append("-");
    directProviderKey.append(c_fmmExaFmmDirectProviderOptionName);

    std::string directRangeKey;
    directRangeKey.append(FmmModuleInfo::sc_name);
    directRangeKey.append("-");
    directRangeKey.append(c_fmmExaFmmDirectRangeOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));
    mdpValueBuilder.rootObject().addValue(directProviderKey,
                                          fmmDirectProviderName(FmmDirectProvider::Gromacs));
    mdpValueBuilder.rootObject().addValue(directRangeKey, std::string("2"));

    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModulePreProcessingSetup();

    EXPECT_EQ(warningHandler_.errorCount(), 0);
}

TEST_F(FmmMdpValidatorTest, ReportsErrorWhenColumbTypeIsNotFmm)
{
    mdCoulombType_ = { CoulombInteractionType::Pme };

    KeyValueTreeBuilder mdpValueBuilder;
    std::string         fmmActiveBackendKey;
    fmmActiveBackendKey.append(FmmModuleInfo::sc_name);
    fmmActiveBackendKey.append("-");
    fmmActiveBackendKey.append(c_fmmActiveOptionName);

    mdpValueBuilder.rootObject().addValue(fmmActiveBackendKey, fmmBackendName(ActiveFmmBackend::ExaFmm));

    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    fmmModule_ = FmmModuleInfo::create();
    fillOptionsFromMdpValues(mdpOptionsTree, fmmModule_->mdpOptionProvider());
    notifyFmmModulePreProcessingSetup();

    EXPECT_EQ(warningHandler_.errorCount(), 1);
}

} // namespace test

} // namespace gmx
