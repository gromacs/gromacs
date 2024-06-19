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
 * Defines gmx::ObservablesReducer and its builder
 *
 * These are defined in the same translation unit so that the Impl
 * object of the ObservablesReducer can be built by the builder.
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */
#include "gmxpre.h"

#include "observablesreducer.h"

#include <cstddef>

#include <algorithm>
#include <numeric>
#include <utility>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

//! Impl class for ObservablesReducer
class ObservablesReducer::Impl
{
public:
    //! Constructor
    Impl(std::vector<double>&&                                            communicationBuffer,
         std::vector<ObservablesReducerBuilder::CallbackAfterReduction>&& registeredCallbacks) :
        communicationBuffer_(std::move(communicationBuffer)),
        registeredCallbacks_(std::move(registeredCallbacks))
    {
    }

    /*! \brief May be called by any subscribed module, via the callback
     * supplied by ObservablesReducerBuilder::build().
     *
     * If this method has been called on this rank since the last call
     * to reduceComplete() with a \c requirement to communicate soon,
     * then it will make the communication buffer available via
     * communicationBuffer() to \c compute_globals(), so it can copy
     * it to the buffer it uses for MPI communication.
     *
     * It is the subscribers' responsibility to coordinate so that
     * all subscribers on all ranks agree on the need to
     * communicate, e.g. by orchestating communication based on
     * the current step number or a previous message.
     *
     * Does not check that the callback corresponds to a module that
     * subscribed to the builder().
     *
     * Returns the status of the ObservablesReducer about whether
     * reduction has already been called on this step.
     */
    ObservablesReducerStatus requireReduction(int callbackIndex, ReductionRequirement requirement)
    {
        if (requirement == ReductionRequirement::Soon)
        {
            reduceSoon_ = true;
        }
        callbacksAfterReduction_.push_back(callbackIndex);
        return status_;
    }

    /*! \brief Storage for communication buffer
     *
     * Must never be resized, because that would potentially
     * invalidate views of it held by the subscribers. */
    std::vector<double> communicationBuffer_;
    //! Registered callbacks that might be used after communication
    std::vector<ObservablesReducerBuilder::CallbackAfterReduction> registeredCallbacks_;
    //! Indices into registeredCallbacks_ of callbacks to use after the next reduction
    std::vector<int> callbacksAfterReduction_;
    /*! \brief Whether the reduction will occur soon because a
     * module required it
     *
     * "Soon" means this step or next, depending when during the step
     * it was required, as there is only one point during a normal
     * simulator step where observables reduction might occur. */
    bool reduceSoon_ = false;
    //! Whether reduction has taken place this step
    ObservablesReducerStatus status_ = ObservablesReducerStatus::ReadyToReduce;
};

ObservablesReducer::ObservablesReducer(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}

ObservablesReducer::ObservablesReducer(ObservablesReducer&&) noexcept = default;

ObservablesReducer::~ObservablesReducer() = default;

ObservablesReducer& ObservablesReducer::operator=(ObservablesReducer&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

bool ObservablesReducer::isReductionRequired() const
{
    return impl_->reduceSoon_;
}

ArrayRef<double> ObservablesReducer::communicationBuffer(const bool reductionRequiredExternally)
{
    // If there's reduction to do, and some module here or externally
    // required reduction, then return a view of the reduction buffer.
    if (!impl_->callbacksAfterReduction_.empty() && (impl_->reduceSoon_ || reductionRequiredExternally))
    {
        return impl_->communicationBuffer_;
    }
    // Nothing to reduce
    return {};
}

void ObservablesReducer::reductionComplete(Step step)
{
    impl_->status_ = ObservablesReducerStatus::AlreadyReducedThisStep;
    for (auto& callbackIndex : impl_->callbacksAfterReduction_)
    {
        impl_->registeredCallbacks_[callbackIndex](step);
    }
    // Prepare for the next reduction
    std::fill(impl_->communicationBuffer_.begin(), impl_->communicationBuffer_.end(), 0.0);
    impl_->callbacksAfterReduction_.clear();
    impl_->reduceSoon_ = false;
}

void ObservablesReducer::markAsReadyToReduce()
{
    impl_->status_ = ObservablesReducerStatus::ReadyToReduce;
}

//! Impl class for ObservablesReducerBuilder
class ObservablesReducerBuilder::Impl
{
public:
    //! Data required to set up a subscription
    struct Subscription
    {
        //! Number of doubles required to reduce
        int sizeRequired;
        //! The callback to notify of the view and future callback to require reduction
        ObservablesReducerBuilder::CallbackFromBuilder callbackFromBuilder;
        //! The callback later used by ObservablesReducer after reduction events
        ObservablesReducerBuilder::CallbackAfterReduction callbackAfterReduction;
    };

    //! Contains all subscriptions received
    std::vector<Subscription> subscriptions_;
    //! Whether build() has already been called on the owner object
    bool buildHasBeenCalled_ = false;
};

ObservablesReducerBuilder::ObservablesReducerBuilder() : impl_(std::make_unique<Impl>()) {}

ObservablesReducerBuilder::ObservablesReducerBuilder(ObservablesReducerBuilder&&) noexcept = default;

ObservablesReducerBuilder& ObservablesReducerBuilder::operator=(ObservablesReducerBuilder&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

ObservablesReducerBuilder::~ObservablesReducerBuilder() = default;

void ObservablesReducerBuilder::addSubscriber(const int                sizeRequired,
                                              CallbackFromBuilder&&    callbackFromBuilder,
                                              CallbackAfterReduction&& callbackAfterReduction)
{
    GMX_RELEASE_ASSERT(!impl_->buildHasBeenCalled_,
                       "Cannot add subscribers to a builder once build() has been called");
    impl_->subscriptions_.emplace_back(Impl::Subscription{
            sizeRequired, std::move(callbackFromBuilder), std::move(callbackAfterReduction) });
}

ObservablesReducer ObservablesReducerBuilder::build()
{
    GMX_RELEASE_ASSERT(!impl_->buildHasBeenCalled_,
                       "Cannot build ObservablesReducer again from the same builder");

    // Prepare the communication buffer
    const int           totalSizeRequired = std::accumulate(impl_->subscriptions_.begin(),
                                                  impl_->subscriptions_.end(),
                                                  0.0,
                                                  [](int subtotal, const auto& subscription) {
                                                      return subtotal + subscription.sizeRequired;
                                                  });
    std::vector<double> communicationBuffer(totalSizeRequired);
    // Set up a view of the communication buffer that we can use after
    // ownership has been transferred to the impl object.
    const ArrayRef<double> bufferView(communicationBuffer);

    std::vector<ObservablesReducerBuilder::CallbackAfterReduction> registeredCallbacksAfterReduction;
    registeredCallbacksAfterReduction.reserve(impl_->subscriptions_.size());
    for (const Impl::Subscription& subscription : impl_->subscriptions_)
    {
        registeredCallbacksAfterReduction.emplace_back(subscription.callbackAfterReduction);
    }

    // Make the impl object so we can set up callbacks to it. This is
    // safe because the impl object is allocated on the heap, so we
    // have a stable value to share with the subscribers.
    auto implPtr = std::make_unique<ObservablesReducer::Impl>(
            std::move(communicationBuffer), std::move(registeredCallbacksAfterReduction));
    auto* impl = implPtr.get();

    // Then use the impl object to make the real one. Note that there
    // is no need for an ObservablesReducer to keep track of all the
    // modules that might notify it, so it doesn't do that.
    ObservablesReducer observablesReducer(std::move(implPtr));

    // Now let the subscribers know how to require reduction in
    // future, and which memory they should use for input and output.
    size_t start                         = 0;
    int    indexToCallbackAfterReduction = 0;
    for (const Impl::Subscription& subscription : impl_->subscriptions_)
    {
        // Construct the callback that will hereafter be owned by the
        // subscriber.
        CallbackToRequireReduction callbackToRequireReduction =
                [impl, indexToCallbackAfterReduction](ReductionRequirement requirement) {
                    return impl->requireReduction(indexToCallbackAfterReduction, requirement);
                };
        subscription.callbackFromBuilder(std::move(callbackToRequireReduction),
                                         bufferView.subArray(start, subscription.sizeRequired));
        start += subscription.sizeRequired;
        ++indexToCallbackAfterReduction;
    }

    impl_->buildHasBeenCalled_ = true;

    return observablesReducer;
}

} // namespace gmx
