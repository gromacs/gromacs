/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

/*! \file
 * \brief Implementation details for MD signalling support.
 *
 * \ingroup gmxapi_md
 */

#include "gmxapi/md/mdsignals.h"

#include <algorithm>
#include <atomic>
#include <memory>
#include <utility>

#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdrun/runner.h"
#include "gromacs/utility/gmxassert.h"

#include "gmxapi/exceptions.h"
#include "gmxapi/session.h"

#include "mdsignals.h"
#include "sessionresources.h"

namespace gmxapi
{

//! \cond
Signal::Signal(Signal&&) noexcept = default;
Signal& Signal::operator=(Signal&&) noexcept = default;

Signal::Signal(std::unique_ptr<SignalImpl> impl) : impl_{ std::move(impl) } {}

Signal::~Signal() = default;

void Signal::operator()()
{
    impl_->call();
}
//! \endcond

void SignalManager::addSignaller(const std::string& name)
{
    called_[name].store(false);
}

/*!
 * Implement the SignalImpl interface to provide a logical AND for managed MD signals.
 *
 * The class is a signal issuer and a signal receiver, but
 * \todo signals received by this operation and received by Mdrunner do not yet have a common interface.
 *
 * Tracks whether each registered input has issued a signal to this operation. When the
 * final registered input issues `call()`, the LogicalAND issues `call()` on the output
 * signal path.
 *
 * State is managed by the parent SignalManager. Client code should get a short-lived handle
 * to a Signal wrapping this implementation object by calling SignalManager::getSignal()
 * with the unique workflow operation name for the block of client code and a gmxapi::md::signals::STOP
 * signal argument.
 *
 * Currently explicitly supports the MD stop signal only.
 *
 * Version gmxapi 0.0.6:  Also, all registered restraints
 * are automatically in the set of ANDed inputs.
 *
 * \ingroup gmxapi_md
 */
class SignalManager::LogicalAND : public Signal::SignalImpl
{
public:
    /*!
     * \brief Create short-lived signal issuer implementation.
     *
     * \param manager
     * \param name
     *
     * Caller is responsible for ensuring that the object pointed to by
     * manager remains valid for the life time of a LogicalAND instance.
     */
    LogicalAND(SignalManager* manager, std::string name) : name_(std::move(name)), manager_(manager)
    {
    }

    //! \cond
    ~LogicalAND() override = default;
    //! \endcond

    /*!
     * \brief Sets the stop condition when the last issuer issues.
     *
     * Once all participating signal issuers have called for a stop signal,
     * the stop condition state is updated to stopAtNextNSStep.
     */
    void call() override
    {
        auto& callCounter = manager_->called_.at(name_);
        callCounter.store(true);
        using pairType = typename decltype(manager_->called_)::value_type;
        if (std::all_of(manager_->called_.cbegin(), manager_->called_.cend(), [](const pairType& p) {
                return p.second.load();
            }))
        {
            *manager_->state_ = gmx::StopSignal::stopAtNextNSStep;
        }
    }

private:
    //! Named signal issuer for the current operation.
    const std::string name_;

    //! The manager that generated this function object.
    SignalManager* manager_;
};

Signal SignalManager::getSignal(const std::string& name, md::signals signal)
{
    if (called_.find(name) == called_.end())
    {
        std::string message = name + " is not registered for this signal.";
        throw gmxapi::ProtocolError(std::move(message));
    }

    if (signal != md::signals::STOP)
    {
        throw gmxapi::MissingImplementationError("This signaller only handles stop signals.");
    }

    auto signalImpl = std::make_unique<LogicalAND>(this, name);
    auto functor    = Signal(std::move(signalImpl));
    return functor;
}

Signal getMdrunnerSignal(SessionResources* resources, md::signals signal)
{
    // while there is only one choice...
    if (signal != md::signals::STOP)
    {
        throw gmxapi::MissingImplementationError("This signaller only handles stop signals.");
    }

    if (resources == nullptr)
    {
        throw gmxapi::UsageError(
                "Caller must provide a valid SessionResources to getMdrunnerSignal.");
    }

    auto signaller = resources->getMdrunnerSignal(signal);

    return signaller;
}

} // end namespace gmxapi
