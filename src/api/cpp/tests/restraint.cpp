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

#include "gromacs/math/vectypes.h"
#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

namespace
{

//! Input file for testing is built by CMake script and filename is compiled into in testingconfiguration binary.
const auto &filename = gmxapi::testing::sample_tprfilename;

/*!
 * \brief Restraint that does nothing other than note whether it has been used.
 */
class NullRestraint : public gmx::IRestraintPotential
{
    public:
        /*! \cond Implement IRestraintPotential */
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double      t) override
        {
            (void)r1;
            (void)r2;
            (void)t;
            hasBeenCalled_ = true;
            return {{0., 0., 0.}, 0.};
        }

        std::vector<int> sites() const override
        {
            return {{0, 1}};
        }

        void bindSession(gmxapi::SessionResources* resources) override
        {
            (void)resources;
        }
        //! \endcond

        /*!
         * \brief Check whether the restraint has been called as an IForceProvider.
         *
         * \return false until restraint is used in MD loop, then true.
         */
        bool hasBeenCalled() const
        {
            return hasBeenCalled_;
        }

    private:
        //! State: false until after the first call to evaluate()
        bool hasBeenCalled_ = false;
};

/*!
 * \brief Wrap a NullRestraint for testing purposes.
 */
class SimpleApiModule : public gmxapi::MDModule
{
    public:
        /*! \cond
         * Implement gmxapi::MDModule interface.
         */
        SimpleApiModule() :
            restraint_(std::make_shared<NullRestraint>())
        {}

        const char *name() const override
        {
            return "SimpleApiModule";
        }

        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            return restraint_;
        }
        //! \endcond

        /*!
         * \brief Check whether the restraint has been called as an IForceProvider.
         *
         * \return false until restraint is used in MD loop, then true.
         */
        bool hasBeenCalled() const
        {
            return restraint_->hasBeenCalled();
        }

    private:
        //! restraint to provide to client or MD simulator, as well as to use to implement hasBeenCalled.
        std::shared_ptr<NullRestraint> restraint_;
};

/*!
 * \brief Check that we can attach a restraint and have it called.
 */
TEST(ApiRunner, RestrainedMD)
{

    auto system = gmxapi::fromTprFile(filename);

    {
        auto           context = std::make_shared<gmxapi::Context>();
        gmxapi::MDArgs args    = gmxapi::testing::mdArgs;
        args.emplace_back("-nsteps");
        args.emplace_back("2");
        // Work around unclean working directory.
        args.emplace_back("-noappend");
        context->setMDArgs(args);

        auto           restraint = std::make_shared<SimpleApiModule>();

        auto           session = system.launch(context);
        EXPECT_TRUE(session != nullptr);
        EXPECT_EQ(restraint->hasBeenCalled(), false);

        gmxapi::addSessionRestraint(session.get(), restraint);
        gmxapi::Status status;
        ASSERT_NO_THROW(status = session->run());
        EXPECT_TRUE(status.success());
        EXPECT_EQ(restraint->hasBeenCalled(), true);

        status = session->close();
        EXPECT_TRUE(status.success());
        EXPECT_EQ(restraint->hasBeenCalled(), true);
    }
}

} // end anonymous namespace
