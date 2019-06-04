//
// Created by Eric Irrgang on 10/13/17.
//

#include "testingconfiguration.h"

#include <iostream>
#include <memory>

#include "harmonicpotential.h"

#include "gmxapi/context.h"
#include "gmxapi/md.h"
#include "gmxapi/session.h"
#include "gmxapi/status.h"
#include "gmxapi/system.h"

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/arrayref.h"

#include <gtest/gtest.h>

namespace {

std::ostream& operator<<(std::ostream& stream, const ::gmx::Vector& vec)
{
    stream << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
    return stream;
}

const auto filename = plugin::testing::sample_tprfilename;

TEST(HarmonicPotentialPlugin, Build)
{
    ASSERT_TRUE(true);
    ASSERT_FALSE(false);

    plugin::Harmonic puller;
}

TEST(HarmonicPotentialPlugin, ForceCalc)
{
    const ::gmx::Vector zerovec = {0, 0, 0};
    // define some unit vectors
    const ::gmx::Vector e1{real(1), real(0), real(0)};
    const ::gmx::Vector e2{real(0), real(1), real(0)};
    const ::gmx::Vector e3{real(0), real(0), real(1)};

    const real R0{1.0};
    const real k{1.0};

    // store temporary values long enough for inspection
    ::gmx::Vector force{};

    plugin::Harmonic puller{R0, k};

    // When input vectors are equal, output vector is meaningless and magnitude is set to zero.
    auto calculateForce = [&puller](const ::gmx::Vector& a, const ::gmx::Vector& b) { return puller.calculate(a,b,0).force; };
    EXPECT_FLOAT_EQ(0., norm(calculateForce(e1, e1)));

    // Default equilibrium distance is 1.0, so force should be zero when norm(r12) == 1.0.
    force = calculateForce(zerovec, e1);
    EXPECT_FLOAT_EQ(0., norm(zerovec - force)) << " where force is (" << force[0] << ", " << force[1] << ", " <<
    force[2] <<
    ")\n";

    force = calculateForce(e1, zerovec);
    EXPECT_FLOAT_EQ(0., norm(zerovec - force)) << " where force is (" << force[0] << ", " << force[1] << ", " <<
    force[2] << ")"
                                                                                                                   "\n";

    force = calculateForce(e1,
                           static_cast<real>(2)*e1);
    EXPECT_FLOAT_EQ(0., norm(zerovec - force)) << " where force is (" << force[0] << ", " << force[1] << ", " <<
    force[2] << ")\n";

    // -kx should give vector (1, 0, 0) when vector r1 == r2 - (2, 0, 0)
    force = calculateForce(static_cast<real>(-2)*e1, zerovec);
    EXPECT_FLOAT_EQ(1., force[0]);
    force = calculateForce(static_cast<real>(-2)*e1, zerovec);
    EXPECT_FLOAT_EQ(0., norm(e1 - force)) << " where force is (" << force[0] << ", " << force[1] << ", " << force[2] <<
    ")\n";

    // -kx should give vector (-2, 0, 0) when vector r1 == r2 + (2, 0, 0)
    force = calculateForce(static_cast<real>(2)*e1,
                           static_cast<real>(-1)*e1);
    EXPECT_FLOAT_EQ(0., norm(static_cast<real>(-2)*e1 - force)) << " where force is (" << force[0] << ", " << force[1]
    <<
    ", " <<
    force[2] << ")\n";
}

TEST(HarmonicPotentialPlugin, EnergyCalc)
{
    const ::gmx::Vector zerovec = {0, 0, 0};
    // define some unit vectors
    const ::gmx::Vector e1{real(1), real(0), real(0)};
    const ::gmx::Vector e2{real(0), real(1), real(0)};
    const ::gmx::Vector e3{real(0), real(0), real(1)};

    const real R0{1.0};
    const real k{1.0};

    // store temporary values long enough for inspection
    real energy{0};

    plugin::Harmonic puller{R0, k};

    // When input vectors are equal, potential energy is still calculable.
    auto calculateEnergy = [&puller](const ::gmx::Vector& a, const ::gmx::Vector& b) { return puller.calculate(a,b,0).energy; };
    EXPECT_FLOAT_EQ(real(0.5*k*R0*R0), calculateEnergy(e1, e1));

    // Default equilibrium distance is 1.0, so energy should be zero when norm(r12) == 1.0.
    energy = calculateEnergy(zerovec, e1);
    EXPECT_FLOAT_EQ(0, energy) << " where energy is " << energy << "\n";

    energy = calculateEnergy(e1, zerovec);
    EXPECT_FLOAT_EQ(0, energy) << " where energy is " << energy << "\n";

    energy = calculateEnergy(e1,
                             static_cast<real>(2)*e1);
    EXPECT_FLOAT_EQ(0, energy) << " where energy is " << energy << "\n";

    // -kx should give vector (1, 0, 0) when vector r1 == r2 - (2, 0, 0)
    energy = calculateEnergy(static_cast<real>(-2)*e1, zerovec);
    EXPECT_FLOAT_EQ(real(0.5*k*R0*R0), energy) << " where energy is " << energy << "\n";

    // -kx should give vector (-2, 0, 0) when vector r1 == r2 + (2, 0, 0)
    energy = calculateEnergy(static_cast<real>(2)*e1,
                             static_cast<real>(-1)*e1);
    EXPECT_FLOAT_EQ(real(0.5*k*4*R0*R0), energy) << " where energy is " << energy << "\n";
}

// This should be part of a validation test, not a unit test.
//TEST(HarmonicPotentialPlugin, Bind)
//{
//
//    {
//        std::string waterfile = "water.tpr";
//        auto system = gmxapi::fromTprFile(waterfile);
//        std::shared_ptr<gmxapi::Context> context = gmxapi::defaultContext();
//
//        auto module = std::make_shared<plugin::HarmonicModule>();
//        module->setParams(1, 4, 2.0, 100.0);
//        system->setRestraint(module);
//
//        auto session = system->launch(context);
//
//        gmxapi::Status status;
//        ASSERT_NO_THROW(status = session->run());
////        ASSERT_TRUE(module->force_called() > 0);
////        ASSERT_NO_THROW(session->run(1000));
//        ASSERT_TRUE(status.success());
//    }
//
//}

} // end anonymous namespace

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
