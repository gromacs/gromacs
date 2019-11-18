//
// Created by Eric Irrgang on 11/9/17.
//

#include "testingconfiguration.h"

#include <memory>
#include <vector>

#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"

#include <gtest/gtest.h>

namespace
{

const auto filename = plugin::testing::sample_tprfilename;

class NullRestraint : public gmx::IRestraintPotential
{
    public:
        gmx::PotentialPointData evaluate(gmx::Vector r1,
                                         gmx::Vector r2,
                                         double t) override
        {
            return {};
        }

        std::vector<int> sites() const override
        {
            return {0,0};
        }
};

class SimpleApiModule : public gmxapi::MDModule
{
    public:
        const char *name() const override
        {
            return "NullApiModule";
        }

        std::shared_ptr<gmx::IRestraintPotential> getRestraint() override
        {
            auto restraint = std::make_shared<NullRestraint>();
            return restraint;
        }
};

} // end anonymous namespace
