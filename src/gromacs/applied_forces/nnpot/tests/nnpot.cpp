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
 * Tests for functionality of the NNPot MDModule
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/nnpot/nnpot.h"

#include <gtest/gtest.h>

#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdrunutility/plainpairlistranges.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/stringcompare.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace test
{

class NNPotTest : public ::testing::Test
{
public:
    void addMdpOptionNNPotActive()
    {
        mdpValueBuilder_.rootObject().addValue("nnpot-active", std::string("true"));
    }

    //! Request the given model inputs and a pair cutoff through the mdp options
    void addMdpOptionModelInputs(const std::vector<std::string>& inputs, const std::string& pairCutoff)
    {
        for (size_t i = 0; i < inputs.size(); ++i)
        {
            const std::string key = "nnpot-model-input" + std::to_string(i + 1);
            mdpValueBuilder_.rootObject().addValue(key, inputs[i]);
        }
        mdpValueBuilder_.rootObject().addValue("nnpot-pair-cutoff", pairCutoff);
    }

    //! Collect the pair-list range the module requests during simulation setup
    std::vector<real> requestedPairlistRanges()
    {
        MDModulesNotifiers notifiers;
        NNPotModule_->subscribeToSimulationSetupNotifications(&notifiers);
        const gmx_mtop_t    mtop;
        const t_inputrec    inputrec;
        PlainPairlistRanges ranges(mtop, inputrec);
        notifiers.simulationSetupNotifier_.notify(&ranges);
        return { ranges.ranges().begin(), ranges.ranges().end() };
    }

    //! Whether the module subscribes to the pair-list construction signal during the run
    bool subscribesToPairlistConstruction()
    {
        MDModulesNotifiers notifiers;
        NNPotModule_->subscribeToSimulationRunNotifications(&notifiers);
        using Signal = const MDModulesPairlistConstructedSignal&;
        return notifiers.simulationRunNotifier_.haveSubscribers<Signal>();
    }

    //! build an mdp options tree that sets the options for the NNPot module
    void makeNNPotModuleWithSetOptions()
    {
        KeyValueTreeObject mdpOptionsTree = mdpValueBuilder_.build();

        NNPotModule_ = NNPotModuleInfo::create();

        test::fillOptionsFromMdpValues(mdpOptionsTree, NNPotModule_->mdpOptionProvider());
    }

    void initializeForceProviders()
    {
        MDModulesNotifiers notifiers;
        NNPotModule_->subscribeToSimulationSetupNotifications(&notifiers);
        notifiers.simulationSetupNotifier_.notify(logger_);
        notifiers.simulationSetupNotifier_.notify(mpiComm_);
        NNPotModule_->initForceProviders(&NNPotForces_);
    }

protected:
    KeyValueTreeBuilder        mdpValueBuilder_;
    ForceProviders             NNPotForces_;
    std::unique_ptr<IMDModule> NNPotModule_;
    MDLogger                   logger_;
    MpiComm                    mpiComm_{ MpiComm::SingleRank{} };
};

//! A model that takes the pair list requests a pair-list range equal to the cutoff
//! and is told when the pair list is built
TEST_F(NNPotTest, AtomPairsInputRequestsPairlist)
{
    addMdpOptionNNPotActive();
    addMdpOptionModelInputs({ "atom-positions", "atom-numbers", "atom-pairs" }, "0.5");
    makeNNPotModuleWithSetOptions();

    const std::vector<real> ranges = requestedPairlistRanges();
    ASSERT_EQ(ranges.size(), 1U);
    EXPECT_FLOAT_EQ(ranges[0], 0.5);
    EXPECT_TRUE(subscribesToPairlistConstruction());
}

//! A model that takes only the periodic pair shifts needs the pair list just as much
TEST_F(NNPotTest, PairShiftsInputRequestsPairlist)
{
    addMdpOptionNNPotActive();
    addMdpOptionModelInputs({ "atom-positions", "atom-numbers", "pair-shifts" }, "0.5");
    makeNNPotModuleWithSetOptions();

    const std::vector<real> ranges = requestedPairlistRanges();
    ASSERT_EQ(ranges.size(), 1U);
    EXPECT_FLOAT_EQ(ranges[0], 0.5);
    EXPECT_TRUE(subscribesToPairlistConstruction());
}

//! Without pair-list inputs the module neither requests a range nor asks to be told
//! about the pair list
TEST_F(NNPotTest, NoPairlistInputRequestsNoPairlist)
{
    addMdpOptionNNPotActive();
    addMdpOptionModelInputs({ "atom-positions", "atom-numbers" }, "0.5");
    makeNNPotModuleWithSetOptions();

    EXPECT_TRUE(requestedPairlistRanges().empty());
    EXPECT_FALSE(subscribesToPairlistConstruction());
}

TEST_F(NNPotTest, ForceProviderLackingInputThrows)
{
    // Prepare MDP inputs
    addMdpOptionNNPotActive();

    // Initialize with default options
    makeNNPotModuleWithSetOptions();

    // Build the force provider
    EXPECT_ANY_THROW(initializeForceProviders());
}

} // namespace test

} // namespace gmx
