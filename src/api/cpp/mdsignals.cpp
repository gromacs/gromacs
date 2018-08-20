/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
//
// Created by Eric Irrgang on 5/18/18.
//

/*! \file
 * \brief Implementation details for MD signalling support.
 *
 * \ingroup gmxapi_md
 */

#include <algorithm>

#include <atomic>

#include "mdsignals-impl.h"
#include "sessionresources-impl.h"
#include "gmxapi/exceptions.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdsignals.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdrun/runner.h"

namespace gmxapi
{


Signal::Signal(Signal &&) noexcept            = default;
Signal &Signal::operator=(Signal &&) noexcept = default;

Signal::Signal(std::unique_ptr<gmxapi::Signal::SignalImpl> &&impl) :
    impl_ {std::move(impl)}
{
}

Signal::~Signal() = default;

void Signal::operator()()
{
    impl_->call();
}

StopSignal::StopSignal(gmx::Mdrunner *runner) : runner_ {runner}
{}

void StopSignal::call()
{
    auto signals = runner_->signals();
    // sig > 0 stops at next NS step. sig < 0 stops at next step.
    signals->at(eglsSTOPCOND).sig = -1;
}

void SignalManager::addSignaller(std::string name)
{
    called_[name].store(false);
}

/*!
 * Implement the SignalImpl interface to provide a logical AND for managed MD signals.
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
        LogicalAND(SignalManager* manager, std::string name) :
            name_ {std::move(name)},
        manager_ {manager}
        {}

        void call() override
        {
            auto &callCounter = manager_->called_.at(name_);
            callCounter.store(true);
            using pairType = typename decltype(manager_->called_) ::value_type;
            if (std::all_of(manager_->called_.cbegin(),
                            manager_->called_.cend(),
                            [](const pairType &p){ return p.second.load(); }))
            {
                StopSignal(manager_->runner_).call();
            }
        }

    private:
        const std::string name_;
        SignalManager   * manager_;
};

Signal SignalManager::getSignal(std::string name,
                                md::signals signal)
{
    if (called_.find(name) == called_.end())
    {
        std::string message = name + " is not registered for this signal.";
        throw gmxapi::ProtocolError(std::move(message));
    }

    if (signal != md::signals::STOP)
    {
        throw gmxapi::NotImplementedError("This signaller only handles stop signals.");
    }

    auto   signalImpl = gmx::compat::make_unique<LogicalAND>(this, name);
    Signal functor {
        std::move(signalImpl)
    };
    return functor;
}

Signal getMdrunnerSignal(SessionResources *resources,
                         md::signals       signal)
{
//// while there is only one choice...
//    if (signal == md::signals::STOP)
//    {
    if (signal != md::signals::STOP)
    {
        throw gmxapi::NotImplementedError("This signaller only handles stop signals.");
    }
    ;

    assert(resources);

    auto signaller = resources->getMdrunnerSignal(signal);

    return signaller;
}

} // end namespace gmxapi
