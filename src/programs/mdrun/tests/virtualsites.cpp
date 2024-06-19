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
 * \brief Sanity checks for virtual sites
 *
 * The tests in this file test the virtual site implementation in two ways.
 * 1) An artificial test system containing all virtual site types is run
 *    end-to-end. The virtual sites are recalculated using a reference
 *    implementation, and compared to the trajectory values. This ensures
 *    that no relevant (real or virtual) coordinates were changed between
 *    virtual site computation, and that the mdrun implementation agrees
 *    with the reference implementation. The latter has the advantage to
 *    be written closer to the analytical expressions, hence easier to
 *    check by eye. Unlike the mdrun implementation, it can also be unit-
 *    tested (see 2)). Since this is an end-to-end test, it also ensures
 *    that virtual sites can be processed by grompp and run by mdrun.
 * 2) The reference implementation is tested to have corresponding positions
 *    and velocities. This is achieved by comparing virtual site positions
 *    calculated from propagated real positions to virtual site positions
 *    calculated by propagating virtual sites using the virtual velocities.
 *
 * Note, this only ensures that the position and velocity implementation
 * match, not that they are actually correct. Some regression test systems
 * include virtual sites, so there is some testing that no bugs are introduced.
 * It would be good to have unit tests, though. This can either be achieved by
 * refactoring the mdrun implementation, or adding them to the reference
 * implementation here. See also #3911.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/mpitest.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"
#include "testutils/trajectoryreader.h"

#include "moduletest.h"
#include "simulatorcomparison.h"

namespace gmx::test
{
namespace
{
using VirtualSiteTestParams = std::tuple<std::string, std::string, std::string>;
class VirtualSiteTest : public MdrunTestFixture, public ::testing::WithParamInterface<VirtualSiteTestParams>
{
public:
    struct VirtualSite;

    /*! \brief Check against reference implementation
     *
     * Check that the positions and velocities of virtual sites are equal to the reference
     * implementation, within a given tolerance. This loops over the trajectory, calculates
     * the virtual site positions and velocities using a reference implementation, and
     * compares it with the reported virtual site positions and velocities.
     *
     * The expectation is that the trajectory and the reference calculations are identical.
     * This ensures that neither the real atoms nor the virtual sites are changed between
     * virtual site calculation and trajectory printing. The reference implementation is
     * also closer to the analytically derived equations than the mdrun implementation,
     * making it easier to verify. Finally, the reference implementation is tested (within
     * this file), which is a reasonable work-around for the fact that the actual implementation
     * doesn't have unit tests.
     *
     * Note that the reference implementation does not take into account PBC. It's intended
     * to be used with simple test cases in which molecules are ensured not to be broken
     * across periodic boundaries.
     */
    static void checkVirtualSitesAgainstReferenceImplementation(const std::string& trajectoryName,
                                                                ArrayRef<const VirtualSite> virtualSites,
                                                                FloatingPointTolerance tolerance)
    {
        SCOPED_TRACE(
                "Checking virtual site positions and velocities against reference implementation.");
        TrajectoryFrameReader trajectoryFrameReader(trajectoryName);
        while (trajectoryFrameReader.readNextFrame())
        {
            const auto frame = trajectoryFrameReader.frame();
            SCOPED_TRACE(formatString("Checking frame %s", frame.frameName().c_str()));
            std::vector<RVec> positions;
            std::vector<RVec> velocities;
            std::vector<RVec> refPositions;
            std::vector<RVec> refVelocities;
            for (const auto& vSite : virtualSites)
            {
                auto [refPosition, refVelocity] = vSite.calculate(frame.x(), frame.v());
                refPositions.emplace_back(refPosition);
                refVelocities.emplace_back(refVelocity);
                positions.emplace_back(frame.x().at(vSite.atomIdx));
                velocities.emplace_back(frame.v().at(vSite.atomIdx));
            }
            EXPECT_THAT(refPositions, Pointwise(RVecEq(tolerance), positions));
            EXPECT_THAT(refVelocities, Pointwise(RVecEq(tolerance), velocities));
        }
    }

    /*! \brief Check the reference implementation
     *
     * This tests that the reference implementation position and velocities correspond.
     * This is done by
     *   a) generating real atom starting positions (x(t)) and half-step velocities (v(t+dt/2))
     *   b) propagating the real atoms positions by one time step x(t+dt) = x(t) + dt*v(t+dt/2)
     *   c) calculating the half-step positions x(t+dt/2) = 0.5 * (x(t) + x(t+dt))
     *   d) calculating the virtual position xv(t) := xv(x(t)) and xv1(t+dt) := xv(x(t+dt))
     *      using the reference implementation
     *   e) calculating the virtual velocity vv(t+dt/2) := vv(x(t+dt/2), v(t+dt/2))
     *      using the reference implementation
     *   f) calculating the virtual position xv2(t+dt) = xv(t) + dt*xv(t+dt/2)
     *   g) comparing xv1(t+dt) and xv2(t+dt)
     * If the calculation of the virtual positions and velocities correspond, xv1 and xv2 will
     * be identical up to some integration error.
     *
     * Maybe unused because this test runs only in double precision.
     */
    [[maybe_unused]] static void checkReferenceImplementation()
    {
        SCOPED_TRACE("Checking virtual site reference implementation.");
        // Randomly generated real atom positions and velocities
        std::vector<RVec> startPositions     = { { 2.641321, 2.076298, 2.138602 },
                                             { 3.776765, 3.154901, 1.556379 },
                                             { 2.376669, 1.166706, 2.457044 },
                                             { 3.242320, 2.142465, 2.023578 } };
        std::vector<RVec> halfStepVelocities = { { 0.154667, 0.319010, 0.458749 },
                                                 { -0.010590, -0.191858, -0.096820 },
                                                 { -0.008609, 0.004656, 0.448852 },
                                                 { 0.411874, -0.038205, -0.151459 } };
        // Virtual site definitions with randomly generated parameters
        std::vector<VirtualSite> virtualSites = {
            { F_VSITE1, 6, { 0 }, {} },
            { F_VSITE2, 6, { 0, 1 }, { 0.710573 } },
            { F_VSITE2FD, 6, { 0, 1 }, { 0.292430 } },
            { F_VSITE3, 6, { 0, 1, 2 }, { 0.060990, 0.543636 } },
            { F_VSITE3FD, 6, { 0, 1, 2 }, { 0.125024, 0.444587 } },
            { F_VSITE3FAD, 6, { 0, 1, 2 }, { 0.414850, 0.349767 } },
            { F_VSITE3OUT, 6, { 0, 1, 2 }, { 0.779323, 0.093773, 0.743164 } },
            { F_VSITE4FDN, 6, { 0, 1, 2, 3 }, { 0.975111, 0.952180, 0.757594 } }
        };

        // Make integration step
        const real        timeStep = 1e-5;
        std::vector<RVec> endPositions;
        std::vector<RVec> halfStepPositions;
        GMX_RELEASE_ASSERT(startPositions.size() == halfStepVelocities.size(),
                           "Need positions and velocities for every real atom.");
        const auto numRealAtoms = startPositions.size();
        for (auto idx = decltype(numRealAtoms){ 0 }; idx < numRealAtoms; idx++)
        {
            endPositions.emplace_back(startPositions[idx] + timeStep * halfStepVelocities[idx]);
            halfStepPositions.emplace_back(startPositions[idx]
                                           + real(0.5) * timeStep * halfStepVelocities[idx]);
        }

        // Check that displacement equals the calculated velocities
        for (const auto& vSite : virtualSites)
        {
            SCOPED_TRACE(formatString("Checking %s", interaction_function[vSite.type].longname));

            /* Calculate start and end virtual position
             *
             * The reference implementation always calculates the virtual velocity, but
             * since we don't have real velocities at full steps, these virtual velocities
             * are unprecise. We don't need them anyway, so we'll just ignore them.
             */
            const auto [startVPosition, vVelocityUnused1] =
                    vSite.calculate(startPositions, halfStepVelocities);
            const auto [endVPosition1, vVelocityUnused2] =
                    vSite.calculate(endPositions, halfStepVelocities);

            /* Calculate virtual velocity at half step using reference implementation
             *
             * The virtual positions are exact, but we don't need them.
             */
            const auto [halfStepVPositionUnused, halfStepVVelocity] =
                    vSite.calculate(halfStepPositions, halfStepVelocities);

            // We can now integrate the virtual positions using the half step velocity
            const auto endVPosition2 = startVPosition + timeStep * halfStepVVelocity;

            /* We can now calculate the displacement of the virtual site in two ways:
             *   (1) Using endVPosition1, calculated by advancing the real positions
             *       and calculating the virtual position from them.
             *   (2) Using endVPosition2, calculated by advancing the virtual position
             *       by the virtual velocity.
             * Comparing the difference of the displacement with relative tolerance makes
             * this at least somewhat independent of the choice of time step and coordinates -
             * assuming that the displacement doesn't go to zero, of course!
             */
            const auto displacement1 = endVPosition1 - startVPosition;
            const auto displacement2 = endVPosition2 - startVPosition;
            ASSERT_GT(std::abs(displacement1[XX]), 0)
                    << "adjust choice of coordinates or time step.";
            ASSERT_GT(std::abs(displacement1[YY]), 0)
                    << "adjust choice of coordinates or time step.";
            ASSERT_GT(std::abs(displacement1[ZZ]), 0)
                    << "adjust choice of coordinates or time step.";

            EXPECT_REAL_EQ_TOL(displacement1[XX],
                               displacement2[XX],
                               relativeToleranceAsFloatingPoint(displacement1[XX], 1e-9));
            EXPECT_REAL_EQ_TOL(displacement1[YY],
                               displacement2[YY],
                               relativeToleranceAsFloatingPoint(displacement1[YY], 1e-9));
            EXPECT_REAL_EQ_TOL(displacement1[ZZ],
                               displacement2[ZZ],
                               relativeToleranceAsFloatingPoint(displacement1[ZZ], 1e-9));
        }
    }

    //! Holds parameters of virtual site and allows calculation
    struct VirtualSite
    {
        //! Type of virtual site
        int type;
        //! Index of virtual site
        int atomIdx;
        //! Indices of constructing atoms
        std::vector<int> constructingAtomIdx;
        //! Construction parameters
        std::vector<real> parameters;

        //! Dispatch function to compute position and velocity of virtual site from reference implementation based on \p type
        [[nodiscard]] std::tuple<RVec, RVec> calculate(ArrayRef<const RVec> constructingPositions,
                                                       ArrayRef<const RVec> constructingVelocities) const;

    private:
        //! Templated reference implementation of virtual site position and velocity calculation
        template<int vsiteType>
        [[nodiscard]] std::tuple<RVec, RVec> calculateVSite(ArrayRef<const RVec> positions,
                                                            ArrayRef<const RVec> velocities) const;
    };

    /*! \brief Helper function returning a list of virtual sites from the topology
     *
     * This also prints the indices of the virtual sites. If any tests fail, this
     * can be used to understand which type is failing.
     */
    static std::vector<VirtualSite> vSiteList(const TopologyInformation& topologyInformation)
    {
        std::vector<VirtualSite> virtualSites;
        const auto&              localTopology = *topologyInformation.expandedTopology();
        printf("Reading virtual site types...\n");
        for (int vsiteType = F_VSITE1; vsiteType <= F_VSITEN; vsiteType++)
        {
            const auto& interactionList = localTopology.idef.il.at(vsiteType);
            if (vsiteType == F_VSITE4FD || interactionList.empty())
            {
                // 4FD is deprecated. Interaction list empty means system doesn't contain this type.
                continue;
            }
            const int   numConstructingAtoms = interaction_function[vsiteType].nratoms - 1;
            const int   defaultIncrement     = numConstructingAtoms + 2;
            std::string indexString;

            for (int i = 0; i < interactionList.size();)
            {
                const int parameterIdx   = interactionList.iatoms[i];
                const int virtualSiteIdx = interactionList.iatoms[i + 1];
                if (!indexString.empty())
                {
                    indexString += ", ";
                }
                indexString += toString(virtualSiteIdx);

                if (vsiteType == F_VSITEN)
                {
                    const int vSiteNConstructingAtoms =
                            localTopology.idef.iparams[parameterIdx].vsiten.n;
                    VirtualSite vSite{ vsiteType, virtualSiteIdx, { interactionList.iatoms[i + 2] }, {} };
                    for (int j = 3; j < 3 * vSiteNConstructingAtoms; j += 3)
                    {
                        vSite.constructingAtomIdx.push_back(interactionList.iatoms[j + 2]);
                        vSite.parameters.push_back(
                                localTopology.idef.iparams[interactionList.iatoms[j]].vsiten.a);
                    }
                    virtualSites.push_back(std::move(vSite));
                    i += 3 * vSiteNConstructingAtoms;
                }
                else
                {
                    virtualSites.emplace_back(VirtualSite{
                            vsiteType,
                            virtualSiteIdx,
                            { interactionList.iatoms.data() + i + 2,
                              interactionList.iatoms.data() + i + 2 + numConstructingAtoms },
                            { localTopology.idef.iparams[parameterIdx].generic.buf,
                              localTopology.idef.iparams[parameterIdx].generic.buf + MAXFORCEPARAM } });
                    i += defaultIncrement;
                }
            }
        }
        return virtualSites;
    }
};

// check-source gets confused by these
//! \cond
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE1>(ArrayRef<const RVec> positions,
                                                       ArrayRef<const RVec> velocities) const
{
    return { positions[constructingAtomIdx.at(0)], velocities[constructingAtomIdx.at(0)] };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE2>(ArrayRef<const RVec> positions,
                                                       ArrayRef<const RVec> velocities) const
{
    const auto& a = parameters[0];
    return { (1 - a) * positions[constructingAtomIdx.at(0)] + a * positions[constructingAtomIdx.at(1)],
             (1 - a) * velocities[constructingAtomIdx.at(0)] + a * velocities[constructingAtomIdx.at(1)] };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE2FD>(ArrayRef<const RVec> positions,
                                                         ArrayRef<const RVec> velocities) const
{
    const auto& a   = parameters[0];
    const auto& ri  = positions[constructingAtomIdx.at(0)];
    const auto& rj  = positions[constructingAtomIdx.at(1)];
    const auto& vi  = velocities[constructingAtomIdx.at(0)];
    const auto& vj  = velocities[constructingAtomIdx.at(1)];
    const auto  rij = rj - ri;
    const auto  vij = vj - vi;

    return { ri + (a / rij.norm()) * rij,
             vi + (a / rij.norm()) * (vij - rij * (vij.dot(rij) / rij.norm2())) };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE3>(ArrayRef<const RVec> positions,
                                                       ArrayRef<const RVec> velocities) const
{
    const auto& a  = parameters[0];
    const auto& b  = parameters[1];
    const auto& ri = positions[constructingAtomIdx.at(0)];
    const auto& rj = positions[constructingAtomIdx.at(1)];
    const auto& rk = positions[constructingAtomIdx.at(2)];
    const auto& vi = velocities[constructingAtomIdx.at(0)];
    const auto& vj = velocities[constructingAtomIdx.at(1)];
    const auto& vk = velocities[constructingAtomIdx.at(2)];

    return { (1 - a - b) * ri + a * rj + b * rk, (1 - a - b) * vi + a * vj + b * vk };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE3FD>(ArrayRef<const RVec> positions,
                                                         ArrayRef<const RVec> velocities) const
{
    const auto& a   = parameters[0];
    const auto& b   = parameters[1];
    const auto& ri  = positions[constructingAtomIdx.at(0)];
    const auto& rj  = positions[constructingAtomIdx.at(1)];
    const auto& rk  = positions[constructingAtomIdx.at(2)];
    const auto& vi  = velocities[constructingAtomIdx.at(0)];
    const auto& vj  = velocities[constructingAtomIdx.at(1)];
    const auto& vk  = velocities[constructingAtomIdx.at(2)];
    const auto  rij = rj - ri;
    const auto  rjk = rk - rj;
    const auto  vij = vj - vi;
    const auto  vjk = vk - vj;

    // TODO: Should be uncommented after resolution of #3909
    const auto rijk = /*(1 - a) **/ rij + a * rjk;
    const auto vijk = /*(1 - a) **/ vij + a * vjk;

    return { ri + (b / rijk.norm()) * rijk,
             vi + (b / rijk.norm()) * (vijk - rijk * (vijk.dot(rijk) / rijk.norm2())) };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE3FAD>(ArrayRef<const RVec> positions,
                                                          ArrayRef<const RVec> velocities) const
{
    // Note: a = d * cos(theta)
    //       b = d * sin(theta)
    const auto& a   = parameters[0];
    const auto& b   = parameters[1];
    const auto& ri  = positions[constructingAtomIdx.at(0)];
    const auto& rj  = positions[constructingAtomIdx.at(1)];
    const auto& rk  = positions[constructingAtomIdx.at(2)];
    const auto& vi  = velocities[constructingAtomIdx.at(0)];
    const auto& vj  = velocities[constructingAtomIdx.at(1)];
    const auto& vk  = velocities[constructingAtomIdx.at(2)];
    const auto  rij = rj - ri;
    const auto  rjk = rk - rj;
    const auto  vij = vj - vi;
    const auto  vjk = vk - vj;

    const auto rPerp        = rjk - rij * (rij.dot(rjk) / rij.norm2());
    const auto dtRijNormRij = (1 / rij.norm()) * (vij - rij * (vij.dot(rij) / rij.norm2()));
    const auto vPerp        = vjk
                       - rij
                                 * ((vij.dot(rjk) + rij.dot(vjk)) / rij.norm2()
                                    - 2 * rij.dot(rjk) * rij.dot(vij) / rij.norm2() / rij.norm2())
                       - vij * (rij.dot(rjk) / rij.norm2());
    const auto dtRPerpNormRPerp =
            (1 / rPerp.norm()) * (vPerp - rPerp * (vPerp.dot(rPerp) / rPerp.norm2()));

    return { ri + (a / rij.norm()) * rij + (b / rPerp.norm()) * rPerp,
             vi + a * dtRijNormRij + b * dtRPerpNormRPerp };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE3OUT>(ArrayRef<const RVec> positions,
                                                          ArrayRef<const RVec> velocities) const
{
    const auto& a   = parameters[0];
    const auto& b   = parameters[1];
    const auto& c   = parameters[2];
    const auto& ri  = positions[constructingAtomIdx.at(0)];
    const auto& rj  = positions[constructingAtomIdx.at(1)];
    const auto& rk  = positions[constructingAtomIdx.at(2)];
    const auto& vi  = velocities[constructingAtomIdx.at(0)];
    const auto& vj  = velocities[constructingAtomIdx.at(1)];
    const auto& vk  = velocities[constructingAtomIdx.at(2)];
    const auto  rij = rj - ri;
    const auto  rik = rk - ri;
    const auto  vij = vj - vi;
    const auto  vik = vk - vi;

    return { ri + a * rij + b * rik + c * rij.cross(rik),
             vi + a * vij + b * vik + c * (vij.cross(rik) + rij.cross(vik)) };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITE4FDN>(ArrayRef<const RVec> positions,
                                                          ArrayRef<const RVec> velocities) const
{
    const auto& a   = parameters[0];
    const auto& b   = parameters[1];
    const auto& c   = parameters[2];
    const auto& ri  = positions[constructingAtomIdx.at(0)];
    const auto& rj  = positions[constructingAtomIdx.at(1)];
    const auto& rk  = positions[constructingAtomIdx.at(2)];
    const auto& rl  = positions[constructingAtomIdx.at(3)];
    const auto& vi  = velocities[constructingAtomIdx.at(0)];
    const auto& vj  = velocities[constructingAtomIdx.at(1)];
    const auto& vk  = velocities[constructingAtomIdx.at(2)];
    const auto& vl  = velocities[constructingAtomIdx.at(3)];
    const auto  rij = rj - ri;
    const auto  rik = rk - ri;
    const auto  ril = rl - ri;
    const auto  vij = vj - vi;
    const auto  vik = vk - vi;
    const auto  vil = vl - vi;

    const auto rja = a * rik - rij;
    const auto rjb = b * ril - rij;
    const auto rm  = rja.cross(rjb);

    const auto vja = a * vik - vij;
    const auto vjb = b * vil - vij;
    const auto vm  = vja.cross(rjb) + rja.cross(vjb);

    return { ri + (c / rm.norm()) * rm, vi + (c / rm.norm()) * (vm - rm * (vm.dot(rm) / rm.norm2())) };
}

template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<F_VSITEN>(ArrayRef<const RVec> positions,
                                                       ArrayRef<const RVec> velocities) const
{
    const auto& ri = positions[constructingAtomIdx.at(0)];
    const auto& vi = velocities[constructingAtomIdx.at(0)];
    GMX_RELEASE_ASSERT(constructingAtomIdx.size() == parameters.size() - 1,
                       "VSITEN atom / parameters mismatch.");

    RVec rSum(0, 0, 0);
    RVec vSum(0, 0, 0);

    const auto parameterSize = parameters.size();
    for (auto idx = decltype(parameterSize){ 0 }; idx < parameterSize; idx++)
    {
        rSum += parameters[idx] * (positions[constructingAtomIdx[idx + 1]] - ri);
        vSum += parameters[idx] * (velocities[constructingAtomIdx[idx + 1]] - vi);
    }

    return { ri + rSum, vi + vSum };
}

[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculate(ArrayRef<const RVec> constructingPositions,
                                        ArrayRef<const RVec> constructingVelocities) const
{
    switch (type)
    {
        case F_VSITE1:
            return calculateVSite<F_VSITE1>(constructingPositions, constructingVelocities);
        case F_VSITE2:
            return calculateVSite<F_VSITE2>(constructingPositions, constructingVelocities);
        case F_VSITE2FD:
            return calculateVSite<F_VSITE2FD>(constructingPositions, constructingVelocities);
        case F_VSITE3:
            return calculateVSite<F_VSITE3>(constructingPositions, constructingVelocities);
        case F_VSITE3FD:
            return calculateVSite<F_VSITE3FD>(constructingPositions, constructingVelocities);
        case F_VSITE3FAD:
            return calculateVSite<F_VSITE3FAD>(constructingPositions, constructingVelocities);
        case F_VSITE3OUT:
            return calculateVSite<F_VSITE3OUT>(constructingPositions, constructingVelocities);
        case F_VSITE4FDN:
            return calculateVSite<F_VSITE4FDN>(constructingPositions, constructingVelocities);
        case F_VSITEN:
            return calculateVSite<F_VSITEN>(constructingPositions, constructingVelocities);
        default: throw NotImplementedError("Unknown virtual site type");
    }
}
//! \endcond

TEST(VirtualSiteVelocityTest, ReferenceIsCorrect)
{
    // Test is too sensitive to run in single precision
    if constexpr (GMX_DOUBLE)
    {
        VirtualSiteTest::checkReferenceImplementation();
    }
}

TEST_P(VirtualSiteTest, WithinToleranceOfReference)
{
    const auto& params         = GetParam();
    const auto& integrator     = std::get<0>(params);
    const auto& tcoupling      = std::get<1>(params);
    const auto& pcoupling      = std::get<2>(params);
    const real  timeStep       = 0.001;
    const auto& simulationName = "vsite_test";

    if (integrator == "md-vv" && pcoupling == "parrinello-rahman")
    {
        // Parrinello-Rahman is not implemented in md-vv
        return;
    }

    if ((integrator == "sd" || integrator == "bd") && tcoupling != "no")
    {
        // bd and sd handle temperature coupling implicitly and would set tcoupling to "no" anyway
        return;
    }

    // Prepare mdp input
    auto mdpFieldValues = prepareMdpFieldValues(simulationName, integrator, tcoupling, pcoupling);
    mdpFieldValues["nsteps"]      = "8";
    mdpFieldValues["nstxout"]     = "4";
    mdpFieldValues["nstvout"]     = "4";
    mdpFieldValues["dt"]          = toString(timeStep);
    mdpFieldValues["constraints"] = "none";
    if (tcoupling != "no" || integrator == "sd" || integrator == "bd")
    {
        mdpFieldValues["tc-grps"] = "system";
        mdpFieldValues["ref-t"]   = "298";
        mdpFieldValues["tau-t"]   = "1";
    }
    if (pcoupling == "parrinello-rahman")
    {
        mdpFieldValues["tau-p"] = "2";
    }

    if (pcoupling == "c-rescale" && tcoupling == "no" && integrator != "sd" && integrator != "bd")
    {
        mdpFieldValues["ensemble-temperature-setting"] = "constant";
        mdpFieldValues["ensemble-temperature"]         = "298";
    }


    // Run grompp
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
    runGrompp(&runner_);
    // Run mdrun
    runMdrun(&runner_);

    TopologyInformation topologyInformation;
    topologyInformation.fillFromInputFile(runner_.tprFileName_);
    const auto virtualSites = vSiteList(topologyInformation);

    // This is in line with other tests (e.g. exact continuation, rerun), which
    // never reach the same reproducibility for BD as for the other integrators.
    const auto tolerance =
            (integrator == "bd") ? relativeToleranceAsUlp(1.0, 100) : defaultRealTolerance();

    checkVirtualSitesAgainstReferenceImplementation(
            runner_.fullPrecisionTrajectoryFileName_, virtualSites, tolerance);
}

INSTANTIATE_TEST_SUITE_P(
        VelocitiesConformToExpectations,
        VirtualSiteTest,
        ::testing::Combine(::testing::Values("md", "md-vv", "sd", "bd"),
                           ::testing::Values("no", "v-rescale", "nose-hoover"),
                           ::testing::Values("no", "c-rescale", "parrinello-rahman")));

} // namespace
} // namespace gmx::test
