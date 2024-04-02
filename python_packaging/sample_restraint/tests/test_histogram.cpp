//
// Created by Eric Irrgang on 3/24/18.
//

#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "ensemblepotential.h"
#include "sessionresources.h"
#include "testingconfiguration.h"

using ::gmx::Vector;

namespace
{

std::ostream& operator<<(std::ostream& stream, const Vector& vec)
{
    stream << "(" << vec[0] << "," << vec[1] << "," << vec[2] << ")";
    return stream;
}

TEST(EnsembleHistogramPotentialPlugin, ForceCalc)
{
    const Vector zerovec = { 0, 0, 0 };
    // define some unit vectors
    const Vector e1{ real(1), real(0), real(0) };
    const Vector e2{ real(0), real(1), real(0) };
    const Vector e3{ real(0), real(0), real(1) };

    const real R0{ 1.0 };
    const real k{ 1.0 };

    // store temporary values long enough for inspection
    Vector force{};

    /*! We need to be able to mock up a Session or otherwise get a dummy SessionResources.
    // Get a dummy EnsembleResources. We aren't testing that here.
    auto dummyFunc = [](const plugin::Matrix<double>&, plugin::Matrix<double>*){
        return;};
    auto resource = std::make_shared<plugin::EnsembleResources>(dummyFunc);
    */

    // Define a reference distribution with a triangular peak at the 1.0 bin.
    const std::vector<double> experimental{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };


    plugin::EnsemblePotential restraint{
        10,           // nbins
        1.0,          // binWidth
        0.0,          // minDist
        10.0,         // maxDist
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

    // With the initial histogram (all zeros) the force should be zero no matter where the particles are.
    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, e1, 0.)));
    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, e2, 0.)));
    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, static_cast<real>(-1) * e1, 0.)));

    /* In 0.0.7, we cannot have a SessionResource without a Session. We can't
     * really do this mock-up test in the 0.0.7 infrastructure without some
     * colossal kludges, so we'll remove this until the new SessionResources
     * protocol is further along.
     * See https://github.com/kassonlab/gmxapi/issues/77
     * See https://github.com/kassonlab/gmxapi/issues/186

    // Establish a history of the atoms being 2.0 apart.
    restraint.callback(e1, static_cast<real>(3)*e1, 0.001, *resource);

    // Atoms should now be driven towards each other where the difference in experimental and historic distributions is greater.
    force = calculateForce(e1, static_cast<real>(3)*e1, 0.001);
    ASSERT_GT(force[0], 0.) << " where force is (" << force[0] << ", " << force[1] << ", " << force[2] << ")\n";
    force = calculateForce(static_cast<real>(3)*e1, e1, 0.001);
    ASSERT_LT(force[0], 0.) << " where force is (" << force[0] << ", " << force[1] << ", " << force[2] << ")\n";

    // When input vectors are equal, output vector is meaningless and magnitude is set to zero.
    ASSERT_EQ(static_cast<real>(0.0), norm(calculateForce(e1, e1, 0.001)));
    */
}

} // end anonymous namespace
