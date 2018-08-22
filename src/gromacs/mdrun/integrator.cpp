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
/*! \internal
 * \brief Defines the dispatch function for the .mdp integrator field.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */

// Including this header before "integrator.h" triggers a bunch of linter errors
// that seem like they might come from conflicting standards-related compiler
// macros (they all seem related to signatures that have std::unique_ptr),
// except that including gmxpre.h first in integrator.h does not solve the linter
// errors. Anyway, it doesn't seem to affect building and testing, but seems like
// something to look into.
#include "gmxpre.h"

#include "integrator.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"

#include "md.h"

namespace gmx
{

//! \brief Run the correct integrator function.
void IntegratorDispatcher::run()
{
    unsigned int ei {
        method_
    };
    switch (ei)
    {
        case eiSteep:
            do_steep();
            break;
        case eiCG:
            do_cg();
            break;
        case eiNM:
            do_nm();
            break;
        case eiLBFGS:
            do_lbfgs();
            break;
        case eiTPI:
        case eiTPIC:
            if (!EI_TPI(ei))
            {
                GMX_THROW(APIError("do_tpi integrator would be called for a non-TPI integrator"));
            }
            do_tpi();
            break;
        case eiSD2_REMOVED:
            GMX_THROW(NotImplementedError("SD2 integrator has been removed"));
        default:
            GMX_THROW(APIError("Non existing integrator selected"));
    }
}

IIntegrator::~IIntegrator() = default;

SimulationMethod IntegratorDispatcher::getMethod() const
{
    return method_;
}

IntegratorBuilder::Base &IntegratorBuilder::Base::addContext(const md::Context &context)
{
    (void)context;
    // default behavior is to ignore the argument. Maybe a usage error should be
    // thrown instead?
    return *this;
}

IntegratorBuilder::DataSentry IntegratorBuilder::Base::setAggregateAdapter
    (std::unique_ptr<IntegratorAggregateAdapter> context)
{
    (void)context;
    // default no-op. There should probably be some sanity checking here.
    return {};
}

IntegratorBuilder::Base::~Base() = default;

class IntegratorDispatcherBuilder : public IntegratorBuilder::Base
{
    public:
        /*!
         * \brief To avoid ambiguity, this type can only be instantiated with
         * its main parameter initialized.
         */
        IntegratorDispatcherBuilder() = delete;
        IntegratorDispatcherBuilder(const IntegratorDispatcherBuilder &)                = delete;
        IntegratorDispatcherBuilder &operator=(const IntegratorDispatcherBuilder &)     = delete;
        IntegratorDispatcherBuilder(IntegratorDispatcherBuilder &&) noexcept            = default;
        IntegratorDispatcherBuilder &operator=(IntegratorDispatcherBuilder &&) noexcept = default;
        /*!
         * \brief Initialize a new builder.
         *
         * \param method Simulation method that our dispatcher will dispatch to.
         */
        explicit IntegratorDispatcherBuilder(SimulationMethod method) :
            method_ {method}
        {
            // For some reason at least one linter thinks this constructor is never used,
            // but it is used by IntegratorBuilder::create
        };

        /*!
         * \brief Create the dispatcher with the catch-all parameter pack.
         */
        IntegratorBuilder::DataSentry setAggregateAdapter(std::unique_ptr<IntegratorAggregateAdapter> container)
        override
        {
            if (paramsContainer_)
            {
                GMX_THROW(APIError("setParams has already been called on this builder."));
            }
            else
            {
                paramsContainer_ = std::move(container);
            }
            IntegratorBuilder::DataSentry sentry;
            // This is probably where we should bind the integrator and data sentry when we have an API to do so.
            return sentry;
        }

        std::unique_ptr<IIntegrator> build() override
        {
            // Indicates a usage error, not a violation of an established invariant.
            if (!paramsContainer_)
            {
                GMX_THROW(APIError("This builder requires setParams() to be called before build()."));
            }

            auto integrator = gmx::compat::make_unique<IntegratorDispatcher>(method_, *paramsContainer_);
            paramsContainer_ = nullptr;
            return integrator;
        }


    private:
        /*!
         * \brief Method that our dispatcher will be dispatching for.
         */
        const SimulationMethod method_;
        /*!
         * \brief Dispatcher instance is not created until setParams is called.
         */
        std::unique_ptr<IntegratorAggregateAdapter> paramsContainer_ {nullptr};
};

IntegratorBuilder::IntegratorBuilder(std::unique_ptr<IntegratorBuilder::Base> impl) :
    impl_ {std::move(impl)}
{
    if (!impl_)
    {
        GMX_THROW(APIError("Builder should not be created with an uninitialized implementation object."));
    }
}

IntegratorBuilder::~IntegratorBuilder() = default;

IntegratorBuilder::IntegratorBuilder(IntegratorBuilder &&) noexcept = default;

IntegratorBuilder &IntegratorBuilder::operator=(IntegratorBuilder &&) noexcept = default;

std::unique_ptr<IIntegrator> IntegratorBuilder::build()
{
    return impl_->build();
}

IntegratorBuilder IntegratorBuilder::create(const SimulationMethod &integratorType)
{
    // Dispatch creation of an appropriate builder for the integration method.
    // Initially, the factory is not extensible.
    std::unique_ptr<IntegratorBuilder::Base> builderImpl;
    auto method = static_cast<decltype(eiMD)>(integratorType.method_);
    switch (method)
    {
        case eiMD:
        case eiBD:
        case eiSD1:
        case eiVV:
        case eiVVAK:
            if (!EI_DYNAMICS(method))
            {
                GMX_THROW(APIError("do_md integrator would be called for a non-dynamical integrator"));
            }
            builderImpl = gmx::compat::make_unique<gmx::MDIntegrator::Builder>();
            break;
        default:
            builderImpl = gmx::compat::make_unique<IntegratorDispatcherBuilder>(integratorType);
    }

    IntegratorBuilder builder {
        std::move(builderImpl)
    };

    assert(builder.impl_);
    return builder;
}

IntegratorBuilder &IntegratorBuilder::addContext(const md::Context &context)
{
    impl_->addContext(context);
    return *this;
}

IntegratorBuilder::DataSentry::~DataSentry() = default;
IntegratorBuilder::DataSentry::DataSentry()  = default;
IntegratorBuilder::DataSentry::DataSentry(IntegratorBuilder::DataSentry &&) noexcept = default;
IntegratorBuilder::DataSentry &
IntegratorBuilder::DataSentry::operator=(IntegratorBuilder::DataSentry &&) noexcept = default;

}  // namespace gmx
