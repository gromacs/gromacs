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
#include <memory>

#include "testingconfiguration.h"
#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"
#include "gmxapi/md/mdmodule.h"
#include <gtest/gtest.h>

#include "gromacs/compat/make_unique.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"

#include "testingconfiguration.h"

namespace
{

const auto          filename = gmxapi::testing::sample_tprfilename;

class DummyMDModule final : public gmx::IMDModule
{
    private:
        class OptionProvider : public   gmx::IMdpOptionProvider
        {
            public:
                void initMdpTransform(gmx::IKeyValueTreeTransformRules*) override
                {
                }

                void initMdpOptions(gmx::IOptionsContainerWithSections*) override
                {
                }

                void buildMdpOutput(gmx::KeyValueTreeObjectBuilder*) const override
                {
                }
        };
        std::shared_ptr<OptionProvider> optionprovider {std::make_shared<OptionProvider>()};

        class OutputProvider : public   gmx::IMDOutputProvider
        {
            public:
                void initOutput(FILE                  *,
                                int,
                                const t_filenm        *,
                                bool,
                                const gmx_output_env_t*) override
                {
                }

                void finishOutput() override
                {
                }
        };
        std::shared_ptr<OutputProvider> outputprovider {std::make_shared<OutputProvider>()};

        class ForceProvider : public    gmx::IForceProvider
        {
            public:

                void calculateForces(const gmx::ForceProviderInput &,
                                     gmx::ForceProviderOutput      *) override
                {
                    force_called++;
                }

                unsigned int force_called {0};
        };
        std::shared_ptr<ForceProvider> forceprovider {std::make_shared<ForceProvider>()};

        gmx::IForceProvider* getForceProvider()
        {
            return forceprovider.get();
        };
    public:
        gmx::IMdpOptionProvider *mdpOptionProvider() override
        {
            return optionprovider.get();
        }

        gmx::IMDOutputProvider *outputProvider() override
        {
            return outputprovider.get();
        }

        void initForceProviders(ForceProviders *forceProviders) override
        {
            forceProviders->addForceProvider(getForceProvider());
        }

        unsigned int force_called() { return forceprovider->force_called; };
};

TEST(ApiRunner, BasicMD)
{

    auto system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
        ASSERT_TRUE(context != nullptr);
        ASSERT_TRUE(system != nullptr);
        gmxapi::MDArgs args = gmxapi::testing::mdArgs;
        args.emplace_back("-nsteps");
        args.emplace_back("10");
        context->setMDArgs(args);
        auto           session = system->launch(context);
        ASSERT_TRUE(session != nullptr);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
//        ASSERT_NO_THROW(session->run(1000));
        ASSERT_TRUE(status.success());
        status = session->close();
        ASSERT_TRUE(status.success());
    }
}

/*!
 * \brief Test our ability to reinitialize the libgromacs environment between simulations.
 */
TEST(ApiRunner, Reinitialize)
{
    std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
    gmxapi::MDArgs                   args    = gmxapi::testing::mdArgs;
    args.emplace_back("-nsteps");
    args.emplace_back("20");

    {
        context->setMDArgs(args);
        auto system  = gmxapi::fromTprFile(filename);
        auto session = system->launch(context);

        // Try to simulate an interrupt signal to catch.
        gmx_set_stop_condition(gmx_stop_cond_next_ns);

        session->run();

        // If this assertion fails, it is not an error, but it indicates expected behavior has
        // changed and we need to consider the impact of whatever changes caused this.
        ASSERT_NE(gmx_get_stop_condition(), gmx_stop_cond_none);

        session->close();
    }   // allow system and session to be destroyed.

    {
        context->setMDArgs(args);
        auto system = gmxapi::fromTprFile(filename);

        // If this assertion fails, it is not an error, but it indicates expected behavior has
        // changed and we need to consider the impact of whatever changes caused this.
        // We are expecting that the libgromacs state has retained the stop condition from the
        // previously issued SIGINT
        ASSERT_NE(gmx_get_stop_condition(), gmx_stop_cond_none);

        auto session = system->launch(context);

        // Launching a session should clear the stop condition
        ASSERT_EQ(gmx_get_stop_condition(), gmx_stop_cond_none);

        session->run();

        // Stop condition should still be clear.
        ASSERT_EQ(gmx_get_stop_condition(), gmx_stop_cond_none);

        session->close();
    }

}

TEST(ApiRunner, ContinuedMD)
{
    // Run a simulation, then extend the target number of steps and continue the simulation
    auto system = gmxapi::fromTprFile(filename);

    {
        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();

        {
            ASSERT_TRUE(context != nullptr);
            ASSERT_TRUE(system != nullptr);
            gmxapi::MDArgs args = gmxapi::testing::mdArgs;
            args.emplace_back("-nsteps");
            args.emplace_back("20");
            context->setMDArgs(args);
            auto           session = system->launch(context);
            ASSERT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            ASSERT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            ASSERT_TRUE(status.success());
        }

        // Reuse the context. Add MD parameters. Run a new session extending the previous trajectory.
        {
            gmxapi::MDArgs args = gmxapi::testing::mdArgs;
            args.emplace_back("-nsteps");
            args.emplace_back("20");
            context->setMDArgs(args);
            auto           session = system->launch(context);
            ASSERT_TRUE(session != nullptr);
            gmxapi::Status status;
            ASSERT_NO_THROW(status = session->run());
            ASSERT_TRUE(status.success());
            ASSERT_NO_THROW(status = session->close());
            ASSERT_TRUE(status.success());
        }
    }
}

} // end anonymous namespace
