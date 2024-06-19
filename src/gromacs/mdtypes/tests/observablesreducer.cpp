/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * Tests for ObservablesReducer.
 *
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "gromacs/mdtypes/observablesreducer.h"

#include <cstddef>

#include <algorithm>
#include <functional>
#include <numeric>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx::test
{
namespace
{

TEST(ObservablesReducerTest, CanMoveAssign)
{
    ObservablesReducerBuilder builder;
    ObservablesReducer        observablesReducer = builder.build();
    EXPECT_FALSE(observablesReducer.isReductionRequired())
            << "no reduction required when no subscribers requested reduction";
    EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
            << "no buffer available when no subscribers requested reduction";
}

TEST(ObservablesReducerTest, CanMoveConstruct)
{
    ObservablesReducerBuilder builder;
    ObservablesReducer        observablesReducerOriginal = builder.build();
    ObservablesReducer        observablesReducer(std::move(observablesReducerOriginal));
    EXPECT_FALSE(observablesReducer.isReductionRequired())
            << "no reduction required when no subscribers requested reduction";
    EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
            << "no buffer available when no subscribers requested reduction";
    observablesReducer.markAsReadyToReduce();
}

TEST(ObservablesReducerTest, CanBuildAndUseWithNoSubscribers)
{
    ObservablesReducerBuilder builder;

    ObservablesReducer observablesReducer = builder.build();
    EXPECT_FALSE(observablesReducer.isReductionRequired())
            << "no reduction required when no subscribers requested reduction";
    EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
            << "no buffer available when no subscribers requested reduction";
    observablesReducer.reductionComplete(0);

    EXPECT_FALSE(observablesReducer.isReductionRequired())
            << "no reduction required when no subscribers requested reduction";
    EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
            << "no buffer available after reductionComplete()";
    observablesReducer.markAsReadyToReduce();
}

TEST(ObservablesReducerTest, CanBuildAndUseWithOneSubscriber)
{
    ObservablesReducerBuilder builder;

    // This test implements the caller, the builder and the
    // ObservablesReducer all in the one scope, which likely does not
    // resemble any actual use case. More realistic test cases are
    // found below.

    std::optional<int>                                    stepUponWhichReductionOccured;
    ObservablesReducerBuilder::CallbackToRequireReduction callbackToRequireReduction;
    ArrayRef<double>                                      bufferView;
    ObservablesReducerBuilder::CallbackFromBuilder        callbackFromBuilder =
            [&](ObservablesReducerBuilder::CallbackToRequireReduction&& c, ArrayRef<double> b) {
                callbackToRequireReduction = std::move(c);
                bufferView                 = b;
            };

    ObservablesReducerBuilder::CallbackAfterReduction callbackAfterReduction =
            [&stepUponWhichReductionOccured](Step step) { stepUponWhichReductionOccured = step; };
    const int requiredBufferSize = 2;
    builder.addSubscriber(
            requiredBufferSize, std::move(callbackFromBuilder), std::move(callbackAfterReduction));

    ObservablesReducer observablesReducer = builder.build();
    EXPECT_FALSE(observablesReducer.isReductionRequired())
            << "no reduction required when no subscribers requested reduction";
    EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
            << "no buffer available when no subscribers requested reduction";
    ASSERT_EQ(requiredBufferSize, bufferView.size());
    ASSERT_NE(callbackToRequireReduction, nullptr)
            << "must have valid callback supplied by the builder";
    EXPECT_FALSE(stepUponWhichReductionOccured.has_value())
            << "no callbacks until reductionComplete() is called";

    // Fill some dummy data, so we can check the zeroing later
    bufferView[0] = 3.0;
    bufferView[1] = 4.0;

    {
        SCOPED_TRACE("Test that ReductionRequirement::Eventually doesn't trigger behavior");

        EXPECT_EQ(callbackToRequireReduction(ReductionRequirement::Eventually),
                  ObservablesReducerStatus::ReadyToReduce);
        EXPECT_FALSE(observablesReducer.isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(observablesReducer.communicationBuffer(false).empty())
                << "no buffer available when the only subscribers requested reduction eventually";
        EXPECT_FALSE(stepUponWhichReductionOccured.has_value())
                << "no callbacks until reductionComplete() is called";

        // Note that there's nothing else to check here, because the
        // empty buffer means that no reduction should take place.
        observablesReducer.markAsReadyToReduce();
    }
    {
        SCOPED_TRACE("Test that ReductionRequirement::Soon does trigger behavior");

        EXPECT_EQ(callbackToRequireReduction(ReductionRequirement::Soon),
                  ObservablesReducerStatus::ReadyToReduce);
        EXPECT_EQ(observablesReducer.communicationBuffer(false).size(), requiredBufferSize)
                << "buffer available when a subscriber requested reduction soon";
        EXPECT_FALSE(stepUponWhichReductionOccured.has_value())
                << "no callbacks until reductionComplete() is called";

        // In the intended use case, some external component must do
        // the actual reduction across ranks using the buffer at this
        // point. Here, we just pretend it happened.

        int step = 2;
        observablesReducer.reductionComplete(step);
        ASSERT_TRUE(stepUponWhichReductionOccured.has_value()) << "reduction callbacks took place";
        EXPECT_EQ(stepUponWhichReductionOccured.value(), step)
                << "reduction step is passed through correctly";
        EXPECT_THAT(bufferView, testing::AllOf(testing::SizeIs(requiredBufferSize), testing::Each(0.0)))
                << "buffer is zeroed after reduction";
        observablesReducer.markAsReadyToReduce();
    }
}

// Integration tests of ObservablesReducer, builder, and fake
// subscriber(s). These will model multiple ranks each with multiple
// subscribers. Building tests that used actual MPI would be extra
// complexity that is not needed at this time.

//! Helper class that models an MD module that needs to make a subscription to \c ObservablesReducer
class Subscriber
{
public:
    //! Ensure that each subscriber sends an interesting amount of data
    static constexpr int s_subscriberBufferMinimumSize = 3;
    /*! \brief Base value used to ensure the data reduced by each
     * subscriber is distinct, to help diagnose bugs.
     *
     * Also contributes to ensuring that the reduced total is
     * never zero.
     *
     * Note that in a real use case, the subscribers will generally be
     * located in multiple modules. */
    static constexpr double s_subscriberOffset = 1000;
    //! Constructor
    Subscriber(int subscriberIndex, int numRanks) :
        // Ensure each subscriber sends a different amount of data, to expose bugs
        sizeRequired_(s_subscriberBufferMinimumSize + subscriberIndex),
        // Ensure each subscriber sends a distinct range of data, to expose bugs
        valueOffset_(s_subscriberOffset * (subscriberIndex + 1)),
        numRanks_(numRanks),
        subscriberIndex_(subscriberIndex)
    {
    }

    //! Make the subscription via the \c observablesReducerBuilder
    void makeSubscription(ObservablesReducerBuilder* observablesReducerBuilder)
    {
        observablesReducerBuilder->addSubscriber(
                sizeRequired_,
                [this](ObservablesReducerBuilder::CallbackToRequireReduction callback,
                       ArrayRef<double>                                      bufferView) {
                    this->callbackWhenBufferAvailable(std::move(callback), bufferView);
                },
                [this](Step step) { this->callbackAfterReduction(step); });
    }

    //! Callback to receive data from the builder
    void callbackWhenBufferAvailable(ObservablesReducerBuilder::CallbackToRequireReduction&& callbackToRequireReduction,
                                     ArrayRef<double> bufferView)
    {
        SCOPED_TRACE("In callback from builder");

        callbackToRequireReduction_ = std::move(callbackToRequireReduction);
        communicationBuffer_        = bufferView;
        EXPECT_THAT(communicationBuffer_, testing::AllOf(testing::SizeIs(sizeRequired_), testing::Each(0.0)))
                << "size of buffer did not match request";
    }

    //! Pretend to do some simulation work characteristic of \c step
    void doSimulationWork(Step step, ReductionRequirement reductionRequirement) const
    {
        // In a real case, MD simulation work for this PP rank and
        // step would go here.
        // ...
        // Then we put values that model its intermediate output into
        // the communication buffer. Those values vary with the step,
        // so that we can test for correctness over multiple reduction
        // events.
        std::iota(communicationBuffer_.begin(), communicationBuffer_.end(), valueOffset_ + double(step));
        // Then we require reduction.
        EXPECT_EQ(callbackToRequireReduction_(reductionRequirement), ObservablesReducerStatus::ReadyToReduce);
    }

    //! After the reduction, check the values for this subscriber are as expected
    void callbackAfterReduction(Step step)
    {
        SCOPED_TRACE(formatString("In callback after reduction for subscriber %d", subscriberIndex_));

        // Expected values are different for each subscriber, and
        // vary with step and number of ranks.
        std::vector<double> expectedResult(communicationBuffer_.size());
        std::iota(expectedResult.begin(), expectedResult.end(), valueOffset_ + double(step));
        std::for_each(expectedResult.begin(), expectedResult.end(), [this](auto& v) {
            v *= this->numRanks_;
        });
        EXPECT_THAT(communicationBuffer_, testing::Pointwise(testing::Eq(), expectedResult))
                << "wrong values were reduced";
        // Ensuring that zero is never computed by a reduction helps
        // test that the zeroing of the communication buffer is
        // working correctly, as we will only observe zero after
        // zeroing and no subsequent activity.
        EXPECT_THAT(communicationBuffer_, testing::Not(testing::Each(0)))
                << "zero may not be the result of an reduction during testing";
    }

    //! The number of doubles required to reduce
    int sizeRequired_;
    //! The callback used to require reduction
    ObservablesReducerBuilder::CallbackToRequireReduction callbackToRequireReduction_;
    //! The buffer used for communication, supplied by an \c ObservablesReducer
    ArrayRef<double> communicationBuffer_;
    //! Offset that differentiates the values reduced by each subscriber
    double valueOffset_;
    //! Number of ranks, used in constructing test expectations
    int numRanks_;
    //! Index within the group of subscribers
    int subscriberIndex_;
};

//! Test fixture class
class ObservablesReducerIntegrationTest : public testing::TestWithParam<std::tuple<int, int>>
{
public:
    //! Helper struct to model data on a single MPI rank
    struct RankData
    {
        /*! \brief Builder of \c observablesReducer for this "rank,"
         * valid until after its build() method has been called. */
        std::optional<ObservablesReducerBuilder> builder = ObservablesReducerBuilder{};
        //! Subscribers to \c observablesReducer
        std::vector<Subscriber> subscribers;
        /*! \brief Manages reduction of observables on behalf of this
         * "rank", valid only after the ObserbalesReducerBuilder
         * builds it. */
        std::optional<ObservablesReducer> observablesReducer;
    };

    //! Constructor
    ObservablesReducerIntegrationTest() : numSubscribers_(std::get<0>(GetParam()))
    {
        int numRanks(std::get<1>(GetParam()));

        rankData_.resize(numRanks);
        for (auto& rankData : rankData_)
        {
            for (int i = 0; i < numSubscribers_; ++i)
            {
                // Ensure each subscriber sends a different (but small) amount of data
                rankData.subscribers.emplace_back(i, numRanks);
            }
            // Now that the addresses of the subscribers are
            // stable, set up the build-time callback.
            for (auto& subscriber : rankData.subscribers)
            {
                subscriber.makeSubscription(&rankData.builder.value());
            }
        }
    }

    /*! \brief Performs the equivalent of MPI_Allreduce on the
     * communication buffer over \c rankData_ */
    void fakeMpiAllReduce(const bool reductionRequiredExternally)
    {
        std::vector<double> reducedValues(rankData_[0]
                                                  .observablesReducer.value()
                                                  .communicationBuffer(reductionRequiredExternally)
                                                  .size(),
                                          0.0);
        // Reduce the values across "ranks"
        for (auto& rankData : rankData_)
        {
            for (size_t i = 0; i != reducedValues.size(); ++i)
            {
                reducedValues[i] += rankData.observablesReducer.value().communicationBuffer(
                        reductionRequiredExternally)[i];
            }
        }
        // Copy the reduced values to all "ranks"
        for (auto& rankData : rankData_)
        {
            auto buffer = rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally);
            std::copy(reducedValues.begin(), reducedValues.end(), buffer.begin());
        }
    }

    //! The number of subscribers
    int numSubscribers_;
    //! Models data distributed over MPI ranks
    std::vector<RankData> rankData_;
};

TEST_P(ObservablesReducerIntegrationTest, CanBuildAndUseSimply)
{
    const bool reductionRequiredExternally = false;
    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer = rankData.builder.value().build();
        rankData.builder.reset();
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available when no subscribers requested reduction";
    }

    Step step = 0;
    for (auto& rankData : rankData_)
    {
        for (auto& subscriber : rankData.subscribers)
        {
            subscriber.doSimulationWork(step, ReductionRequirement::Soon);
        }
        EXPECT_NE(numSubscribers_ == 0, rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_EQ(
                numSubscribers_ == 0,
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "buffer should be available only when there are active subscribers";
    }

    // This does reduction work, and calls the callbacks that check
    // the buffer contents.
    fakeMpiAllReduce(reductionRequiredExternally);

    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer.value().reductionComplete(step);
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available after reductionComplete()";
        rankData.observablesReducer.value().markAsReadyToReduce();
    }
}

TEST_P(ObservablesReducerIntegrationTest, CanBuildAndUseOverMultipleSteps)
{
    const bool reductionRequiredExternally = false;
    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer = rankData.builder.value().build();
        rankData.builder.reset();
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available when no subscribers requested reduction";
    }

    for (Step step = 0; step < 20; step += 10)
    {
        for (auto& rankData : rankData_)
        {
            for (auto& subscriber : rankData.subscribers)
            {
                subscriber.doSimulationWork(step, ReductionRequirement::Soon);
            }
            EXPECT_NE(numSubscribers_ == 0, rankData.observablesReducer.value().isReductionRequired())
                    << "no reduction required when no subscribers requested reduction";
            EXPECT_EQ(numSubscribers_ == 0,
                      rankData.observablesReducer.value()
                              .communicationBuffer(reductionRequiredExternally)
                              .empty())
                    << "buffer should be available only when there are subscribers";
        }

        // This does reduction work, and calls the callbacks that
        // check the buffer contents.
        fakeMpiAllReduce(reductionRequiredExternally);

        for (auto& rankData : rankData_)
        {
            rankData.observablesReducer.value().reductionComplete(step);
            EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                    << "no reduction required after reductionComplete()";
            EXPECT_TRUE(rankData.observablesReducer.value()
                                .communicationBuffer(reductionRequiredExternally)
                                .empty())
                    << "no buffer available after reductionComplete()";
            rankData.observablesReducer.value().markAsReadyToReduce();
        }
    }
}

TEST_P(ObservablesReducerIntegrationTest, CanBuildAndUseWithoutAllNeedingReduction)
{
    if (numSubscribers_ == 0)
    {
        // Test is meaningless with no subscribers
        return;
    }

    const bool reductionRequiredExternally = false;
    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer = rankData.builder.value().build();
        rankData.builder.reset();
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available when no subscribers requested reduction";
    }

    // Only one subscriber does work leading to reduction
    size_t subscriberNeedingReduction = 0;
    Step   step                       = 0;
    for (auto& rankData : rankData_)
    {
        auto& subscriber = rankData.subscribers[subscriberNeedingReduction];
        subscriber.doSimulationWork(step, ReductionRequirement::Soon);
        EXPECT_TRUE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_FALSE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "buffer should be available when there is an active subscriber";
    }

    // This does reduction work, and calls the callbacks that check
    // the buffer contents.
    fakeMpiAllReduce(reductionRequiredExternally);

    // Check that other subscribers didn't reduce anything
    for (auto& rankData : rankData_)
    {
        for (size_t r = 0; r != rankData.subscribers.size(); ++r)
        {
            if (r == subscriberNeedingReduction)
            {
                continue;
            }
            EXPECT_THAT(rankData.subscribers[r].communicationBuffer_, testing::Each(0.0))
                    << "buffer for non-subscribers should be zero";
        }
    }

    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer.value().reductionComplete(step);
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required after reductionComplete()";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available after reductionComplete()";
        rankData.observablesReducer.value().markAsReadyToReduce();
    }
}

TEST_P(ObservablesReducerIntegrationTest, CanBuildAndUseWhenASubscriberUsesEventually)
{
    if (numSubscribers_ < 2)
    {
        // Test is meaningful only with multiple subscribers
        return;
    }

    const bool reductionRequiredExternally = false;
    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer = rankData.builder.value().build();
        rankData.builder.reset();
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available when no subscribers requested reduction";
    }

    // Only one subscriber does work leading to reduction
    size_t subscriberUsingEventually = 1;
    Step   step                      = 1;
    for (auto& rankData : rankData_)
    {
        auto& subscriber = rankData.subscribers[subscriberUsingEventually];
        subscriber.doSimulationWork(step, ReductionRequirement::Eventually);
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "reduction should not be required when the only active subscriber used "
                   "ReductionRequirement::Eventually";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "buffer should not be available when the only active subscriber used "
                   "ReductionRequirement::Eventually";
    }

    // This will do nothing, as all the communication buffers are
    // empty, but we can't directly test that nothing
    // occurred. Instead, we will later do some
    // ReductionRequirement::Soon work and observe that result is
    // consistent with exactly one reduction.
    fakeMpiAllReduce(reductionRequiredExternally);

    for (auto& rankData : rankData_)
    {
        for (size_t i = 0; i != rankData.subscribers.size(); ++i)
        {
            if (i == subscriberUsingEventually)
            {
                continue;
            }
            rankData.subscribers[i].doSimulationWork(step, ReductionRequirement::Soon);
        }
        EXPECT_TRUE(rankData.observablesReducer.value().isReductionRequired())
                << "reduction should be required since there are subscribers";
        EXPECT_FALSE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "buffer should be available since there are subscribers";
    }

    // This does reduction work, and calls the callbacks that check
    // the buffer contents.
    fakeMpiAllReduce(reductionRequiredExternally);

    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer.value().reductionComplete(step);
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required after reductionComplete()";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available after reductionComplete()";
        rankData.observablesReducer.value().markAsReadyToReduce();
    }
}

TEST_P(ObservablesReducerIntegrationTest, CanBuildAndUseWhenAllSubscribersUseEventually)
{
    if (numSubscribers_ < 2)
    {
        // Test is meaningful only with multiple subscribers
        return;
    }

    const bool reductionRequiredExternally = true;
    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer = rankData.builder.value().build();
        rankData.builder.reset();
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required when no subscribers requested reduction";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available when no subscribers requested reduction";
    }

    Step step = 1;
    // All subscribers do work leading to reduction eventually
    for (auto& rankData : rankData_)
    {
        for (size_t i = 0; i != rankData.subscribers.size(); ++i)
        {
            rankData.subscribers[i].doSimulationWork(step, ReductionRequirement::Eventually);
        }
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "reduction should not be required since there are no subscribers using "
                   "ReductionRequirement::Soon";
        EXPECT_TRUE(rankData.observablesReducer.value().communicationBuffer(false).empty())
                << "buffer should not be available unless reduction is required externally";
        EXPECT_FALSE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "buffer should be available since there are subscribers and reduction is "
                   "required externally";
    }

    // This does reduction work, and calls the callbacks that check
    // the buffer contents.
    fakeMpiAllReduce(reductionRequiredExternally);

    for (auto& rankData : rankData_)
    {
        rankData.observablesReducer.value().reductionComplete(step);
        EXPECT_FALSE(rankData.observablesReducer.value().isReductionRequired())
                << "no reduction required after reductionComplete()";
        EXPECT_TRUE(
                rankData.observablesReducer.value().communicationBuffer(reductionRequiredExternally).empty())
                << "no buffer available after reductionComplete() even when reduction required "
                   "externally";
        rankData.observablesReducer.value().markAsReadyToReduce();
    }
}

//! Help GoogleTest name our test cases
std::string namesOfTests(const testing::TestParamInfo<ObservablesReducerIntegrationTest::ParamType>& info)
{
    // NB alphanumeric characters only
    return formatString("numSubscribers%dnumRanks%d", std::get<0>(info.param), std::get<1>(info.param));
}
INSTANTIATE_TEST_SUITE_P(WithVariousSubscriberCounts,
                         ObservablesReducerIntegrationTest,
                         testing::Combine(testing::Values(0, 1, 2, 3), // subscriber counts
                                          testing::Values(1, 2, 3)),   // rank counts
                         namesOfTests);

} // namespace
} // namespace gmx::test
