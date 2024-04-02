//
// Created by Eric Irrgang on 3/24/18.
//

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "ensemblepotential.h"
#include "testingconfiguration.h"

namespace
{

using ::gmx::Vector;

std::ostream& operator<<(std::ostream& stream, const gmx::Vector& vec)
{
    stream << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
    return stream;
}

const auto filename = plugin::testing::sample_tprfilename;

TEST(EnsembleBoundingPotentialPlugin, ForceCalc)
{
    const Vector zerovec{ 0, 0, 0 };
    // define some unit vectors
    const Vector e1{ real(1), real(0), real(0) };
    const Vector e2{ real(0), real(1), real(0) };
    const Vector e3{ real(0), real(0), real(1) };

    const real R0{ 1.0 };
    const real k{ 1.0 };

    // store temporary values long enough for inspection
    Vector force{};

    // Get a dummy EnsembleResources. We aren't testing that here.
    auto dummyFunc = [](const plugin::Matrix<double>&, plugin::Matrix<double>*) { return; };
    auto resource  = std::make_shared<plugin::Resources>(dummyFunc);

    // Define a reference distribution with a triangular peak at the 1.0 bin.
    const std::vector<double> experimental{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };


    plugin::EnsemblePotential restraint{
        10,           // nbins
        1.0,          // binWidth
        5.0,          // minDist
        5.0,          // maxDist
        experimental, // experimental reference histogram
        1,            // nSamples
        0.001,        // samplePeriod
        1,            // nWindows
        100.,         // k
        1.0           // sigma
    };

    auto calculateForce = [&restraint](const Vector& a, const Vector& b, double t) {
        return restraint.calculate(a, b, t).force;
    };

    // Atoms should be driven towards each other when above maxDist and and away under minDist.
    force = calculateForce(e1, static_cast<real>(3) * e1, 0.001);
    ASSERT_LT(force[0], 0.) << " where force is (" << force[0] << ", " << force[1] << ", "
                            << force[2] << ")\n";
    force = calculateForce(e1, static_cast<real>(7) * e1, 0.001);
    ASSERT_GT(force[0], 0.) << " where force is (" << force[0] << ", " << force[1] << ", "
                            << force[2] << ")\n";
}

} // end anonymous namespace
