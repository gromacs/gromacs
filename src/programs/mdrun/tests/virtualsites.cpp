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
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/hardware/hw_info.h"
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
#include "gromacs/utility/vectypes.h"

#include "testutils/mpitest.h"
#include "testutils/naming.h"
#include "testutils/setenv.h"
#include "testutils/simulationdatabase.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"
#include "testutils/trajectoryreader.h"

#include "moduletest.h"
#include "simulatorcomparison.h"
#include "trajectorycomparison.h"

namespace gmx
{
namespace test
{
namespace
{
/*! \brief Simulation name, bond type, whether to test N-type virtual sites,
 * integrator, temperature coupling and pressure coupling
 *
 * The bond type is either InteractionFunction::GROMOS96Bonds or
 * InteractionFunction::Bonds (harmonic). GROMOS96 bonds can not be
 * computed on GPUs, so they are computed on the host also when the other
 * bonded interactions are offloaded, while harmonic bonds can be offloaded.
 *
 * When N-type virtual sites are not tested, the test system constructs
 * them with other types instead, if it has any.
 */
using VirtualSiteTestParams =
        std::tuple<std::string, InteractionFunction, bool, std::string, std::string, std::string>;

//! Returns the name of the bond type, for use in test and file names
std::string bondTypeName(const InteractionFunction bondType)
{
    return interaction_function[bondType].name;
}

//! Returns a name for whether N-type virtual sites are tested, for use in test and file names
std::string virtualSitesNName(const bool testVirtualSitesN)
{
    return testVirtualSitesN ? "VSITEN" : "NO_VSITEN";
}

//! Names the tests from their parameters
const NameOfTestFromTuple<VirtualSiteTestParams> sc_testNamer{
    std::make_tuple(useString, bondTypeName, virtualSitesNName, useString, useString, useString)
};

class VirtualSiteTest : public MdrunTestFixture, public ::testing::WithParamInterface<VirtualSiteTestParams>
{
public:
    struct VirtualSite;

    /*! \brief Check against reference implementation
     *
     * Check that the positions and velocities of virtual sites are equal to the reference
     * implementation, within a given tolerance. This loops over the trajectory, calculates
     * the virtual site positions and velocities using a reference implementation, and
     * compares it with the reported virtual site positions and velocities. It also checks that
     * forces on virtual sites are cleared after spreading.
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
            // Output intervals differ, so a frame can contain any subset of x, v and f.
            // Positions can be checked in all frames with x, velocities only in frames
            // which also have x, since the reference velocities depend on the positions.
            const bool haveX = !frame.x().empty();
            const bool haveV = !frame.v().empty();
            const bool haveF = !frame.f().empty();
            if (haveX)
            {
                const std::vector<RVec> zeroVelocities(haveV ? 0 : frame.x().size(), { 0, 0, 0 });
                const ArrayRef<const RVec> frameVelocities =
                        haveV ? frame.v() : ArrayRef<const RVec>(zeroVelocities);
                std::vector<RVec> positions;
                std::vector<RVec> velocities;
                std::vector<RVec> refPositions;
                std::vector<RVec> refVelocities;
                for (const auto& vSite : virtualSites)
                {
                    auto [refPosition, refVelocity] = vSite.calculate(frame.x(), frameVelocities);
                    refPositions.emplace_back(refPosition);
                    refVelocities.emplace_back(refVelocity);
                    positions.emplace_back(frame.x().at(vSite.atomIdx));
                    velocities.emplace_back(frameVelocities[vSite.atomIdx]);
                }
                EXPECT_THAT(refPositions, Pointwise(RVecEq(tolerance), positions));
                if (haveV)
                {
                    EXPECT_THAT(refVelocities, Pointwise(RVecEq(tolerance), velocities));
                }
            }
            if (haveF)
            {
                // The forces on the virtual sites are spread to the constructing atoms
                // and then cleared, so the output must have exactly zero force on them.
                std::vector<RVec> forces;
                for (const auto& vSite : virtualSites)
                {
                    forces.emplace_back(frame.f().at(vSite.atomIdx));
                }
                const std::vector<RVec> zeroForces(forces.size(), { 0, 0, 0 });
                EXPECT_THAT(zeroForces, Pointwise(RVecEq(absoluteTolerance(0)), forces));
            }
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
            { InteractionFunction::VirtualSite1, 6, { 0 }, {} },
            { InteractionFunction::VirtualSite2, 6, { 0, 1 }, { 0.710573 } },
            { InteractionFunction::VirtualSite2FixedDistance, 6, { 0, 1 }, { 0.292430 } },
            { InteractionFunction::VirtualSite3, 6, { 0, 1, 2 }, { 0.060990, 0.543636 } },
            { InteractionFunction::VirtualSite3FixedDistance, 6, { 0, 1, 2 }, { 0.125024, 0.444587 } },
            { InteractionFunction::VirtualSite3FixedAngleDistance, 6, { 0, 1, 2 }, { 0.414850, 0.349767 } },
            { InteractionFunction::VirtualSite3Outside, 6, { 0, 1, 2 }, { 0.779323, 0.093773, 0.743164 } },
            { InteractionFunction::VirtualSite4FixedDistanceNormalization,
              6,
              { 0, 1, 2, 3 },
              { 0.975111, 0.952180, 0.757594 } },
            // The weight of the first constructing atom is implied by the weights summing to 1
            { InteractionFunction::VirtualSiteN, 6, { 0, 1, 2, 3 }, { 0.211335, 0.197166, 0.394332 } }
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

        /* The displacements are differences of positions of magnitude up to ~4 nm (see vsite_test.gro),
         * so they cannot be more accurate than a few ULP of those positions. When a displacement
         * component is small (e.g. when the weighted velocities nearly cancel), a purely
         * relative tolerance would fall below that rounding floor. We therefore accept
         * differences within either an absolute rounding floor or a relative tolerance.
         * This is still far tighter than any error from a wrong velocity expression,
         * which would be of the order of the displacement itself.
         */
        constexpr double c_maxPositionMagnitude = 4.0;
        constexpr double c_absoluteTolerance    = 8 * GMX_DOUBLE_EPS * c_maxPositionMagnitude;
        const FloatingPointTolerance displacementTolerance(
                0, c_absoluteTolerance, 0, 1e-9, UINT64_MAX, UINT64_MAX, false);

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

            EXPECT_REAL_EQ_TOL(displacement1[XX], displacement2[XX], displacementTolerance);
            EXPECT_REAL_EQ_TOL(displacement1[YY], displacement2[YY], displacementTolerance);
            EXPECT_REAL_EQ_TOL(displacement1[ZZ], displacement2[ZZ], displacementTolerance);
        }
    }

    //! Holds parameters of virtual site and allows calculation
    struct VirtualSite
    {
        //! Type of virtual site
        InteractionFunction type;
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
        template<InteractionFunction vsiteType>
        [[nodiscard]] std::tuple<RVec, RVec> calculateVSite(ArrayRef<const RVec> positions,
                                                            ArrayRef<const RVec> velocities) const;
    };

    /*! \brief Run a simulation and check its virtual sites against the reference implementation
     *
     * Prepares the tpr file for the given simulation and settings, runs mdrun (passing
     * \p mdrunOptions on to mdrun), and compares the virtual sites in the resulting
     * trajectory to the reference implementation.
     */
    void runTest(const std::vector<SimulationOptionTuple>& mdrunOptions)
    {
        prepareTpr();
        if (HasFatalFailure())
        {
            return;
        }
        runMdrunAndCheckVirtualSites(std::get<3>(GetParam()), mdrunOptions);
    }

    //! Runs grompp to make the tpr file of \c runner_ for the settings of this test
    void prepareTpr()
    {
        const auto& [simulationName, bondType, testVirtualSitesN, integrator, tcoupling, pcoupling] =
                GetParam();

        const real timeStep = 0.001;

        // Prepare mdp input
        auto mdpFieldValues = prepareMdpFieldValues(simulationName, integrator, tcoupling, pcoupling);
        // The output and coupling intervals are chosen to be coprime, so that steps
        // with every combination of x, v and f output occur, as well as every
        // combination of virial computation with any two kinds of output.
        // (Every combination of all four would need 140 steps.)
        // Choosing them as multiples of each other would hide errors in quantities
        // which are only computed or transferred for one kind of output.
        mdpFieldValues["nsteps"]        = "36";
        mdpFieldValues["nstxout"]       = "3";
        mdpFieldValues["nstvout"]       = "4";
        mdpFieldValues["nstfout"]       = "5";
        mdpFieldValues["nstcalcenergy"] = "7";
        mdpFieldValues["nstenergy"]     = "7";
        mdpFieldValues["nstcomm"]       = "7";
        mdpFieldValues["nsttcouple"]    = "7";
        mdpFieldValues["nstpcouple"]    = "7";
        mdpFieldValues["dt"]            = toString(timeStep);
        mdpFieldValues["constraints"]   = "none";
        GMX_RELEASE_ASSERT(bondType == InteractionFunction::GROMOS96Bonds
                                   || bondType == InteractionFunction::Bonds,
                           "The test systems only have GROMOS96 and harmonic bond variants");
        std::string defines;
        if (bondType == InteractionFunction::Bonds)
        {
            defines += " -DHARMONIC_BONDS";
        }
        if (!testVirtualSitesN && simulationName == "vsite_test")
        {
            defines += " -DREPLACE_VIRTUAL_SITES_N";
        }
        if (!defines.empty())
        {
            mdpFieldValues["define"] = defines;
        }
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

        if (pcoupling == "c-rescale" && tcoupling == "no" && integrator != "sd"
            && integrator != "bd")
        {
            mdpFieldValues["ensemble-temperature-setting"] = "constant";
            mdpFieldValues["ensemble-temperature"]         = "298";
        }

        runner_.useTopGroAndNdxFromDatabase(simulationName);
        runner_.useStringAsMdpFile(prepareMdpFileContents(mdpFieldValues));
        ASSERT_EQ(0, runner_.callGrompp());
    }

    /*! \brief Run mdrun on the prepared tpr file and check the virtual sites in the trajectory
     *
     * Passes \p mdrunOptions on to mdrun, and compares the virtual sites in the resulting
     * trajectory to the reference implementation.
     */
    void runMdrunAndCheckVirtualSites(const std::string&                        integrator,
                                      const std::vector<SimulationOptionTuple>& mdrunOptions)
    {
        runMdrun(&runner_, mdrunOptions);

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
        for (InteractionFunction vsiteType : { InteractionFunction::VirtualSite1,
                                               InteractionFunction::VirtualSite2,
                                               InteractionFunction::VirtualSite2FixedDistance,
                                               InteractionFunction::VirtualSite3,
                                               InteractionFunction::VirtualSite3FixedDistance,
                                               InteractionFunction::VirtualSite3FixedAngleDistance,
                                               InteractionFunction::VirtualSite3Outside,
                                               InteractionFunction::VirtualSite4FixedDistance,
                                               InteractionFunction::VirtualSite4FixedDistanceNormalization,
                                               InteractionFunction::VirtualSiteN })
        {
            const auto& interactionList = localTopology.idef.il[vsiteType];
            if (vsiteType == InteractionFunction::VirtualSite4FixedDistance || interactionList.empty())
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

                if (vsiteType == InteractionFunction::VirtualSiteN)
                {
                    const int vSiteNConstructingAtoms =
                            localTopology.idef.iparams[parameterIdx].vsiten.n;
                    VirtualSite vSite{ vsiteType, virtualSiteIdx, { interactionList.iatoms[i + 2] }, {} };
                    for (int j = 3; j < 3 * vSiteNConstructingAtoms; j += 3)
                    {
                        vSite.constructingAtomIdx.push_back(interactionList.iatoms[i + j + 2]);
                        vSite.parameters.push_back(
                                localTopology.idef.iparams[interactionList.iatoms[i + j]].vsiten.a);
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite1>(ArrayRef<const RVec> positions,
                                                                                ArrayRef<const RVec> velocities) const
{
    return { positions[constructingAtomIdx.at(0)], velocities[constructingAtomIdx.at(0)] };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite2>(ArrayRef<const RVec> positions,
                                                                                ArrayRef<const RVec> velocities) const
{
    const auto& a = parameters[0];
    return { (1 - a) * positions[constructingAtomIdx.at(0)] + a * positions[constructingAtomIdx.at(1)],
             (1 - a) * velocities[constructingAtomIdx.at(0)] + a * velocities[constructingAtomIdx.at(1)] };
}
template<>
[[nodiscard]] std::tuple<RVec, RVec>
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite2FixedDistance>(
        ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite3>(ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite3FixedDistance>(
        ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite3FixedAngleDistance>(
        ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite3Outside>(
        ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSite4FixedDistanceNormalization>(
        ArrayRef<const RVec> positions,
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
VirtualSiteTest::VirtualSite::calculateVSite<InteractionFunction::VirtualSiteN>(ArrayRef<const RVec> positions,
                                                                                ArrayRef<const RVec> velocities) const
{
    const auto& ri = positions[constructingAtomIdx.at(0)];
    const auto& vi = velocities[constructingAtomIdx.at(0)];
    GMX_RELEASE_ASSERT(parameters.size() == constructingAtomIdx.size() - 1,
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
        case InteractionFunction::VirtualSite1:
            return calculateVSite<InteractionFunction::VirtualSite1>(constructingPositions,
                                                                     constructingVelocities);
        case InteractionFunction::VirtualSite2:
            return calculateVSite<InteractionFunction::VirtualSite2>(constructingPositions,
                                                                     constructingVelocities);
        case InteractionFunction::VirtualSite2FixedDistance:
            return calculateVSite<InteractionFunction::VirtualSite2FixedDistance>(
                    constructingPositions, constructingVelocities);
        case InteractionFunction::VirtualSite3:
            return calculateVSite<InteractionFunction::VirtualSite3>(constructingPositions,
                                                                     constructingVelocities);
        case InteractionFunction::VirtualSite3FixedDistance:
            return calculateVSite<InteractionFunction::VirtualSite3FixedDistance>(
                    constructingPositions, constructingVelocities);
        case InteractionFunction::VirtualSite3FixedAngleDistance:
            return calculateVSite<InteractionFunction::VirtualSite3FixedAngleDistance>(
                    constructingPositions, constructingVelocities);
        case InteractionFunction::VirtualSite3Outside:
            return calculateVSite<InteractionFunction::VirtualSite3Outside>(constructingPositions,
                                                                            constructingVelocities);
        case InteractionFunction::VirtualSite4FixedDistanceNormalization:
            return calculateVSite<InteractionFunction::VirtualSite4FixedDistanceNormalization>(
                    constructingPositions, constructingVelocities);
        case InteractionFunction::VirtualSiteN:
            return calculateVSite<InteractionFunction::VirtualSiteN>(constructingPositions,
                                                                     constructingVelocities);
        default: throw NotImplementedError("Unknown virtual site type");
    }
}
//! \endcond

TEST(VirtualSiteVelocityTest, ReferenceIsCorrect)
{
    if (!GMX_DOUBLE)
    {
        GTEST_SKIP() << "Test is too sensitive to run in single precision";
    }
    VirtualSiteTest::checkReferenceImplementation();
}

TEST_P(VirtualSiteTest, WithinToleranceOfReference)
{
    const auto& [simulationName, bondType, testVirtualSitesN, integrator, tcoupling, pcoupling] =
            GetParam();

    if (integrator == "md-vv" && pcoupling == "parrinello-rahman")
    {
        GTEST_SKIP() << "Parrinello-Rahman is not implemented in md-vv";
    }

    if ((integrator == "sd" || integrator == "bd") && tcoupling != "no")
    {
        GTEST_SKIP() << "bd and sd handle temperature coupling implicitly and would set tcoupling "
                        "to \"no\" anyway";
    }

    // The test binary is registered with CTest both with one rank and with multiple ranks.
    // The single-rank registration is needed for VirtualSiteGpuTest, but it also runs these
    // cases again. When the GROMACS CI system runs the tests under compute-sanitizer
    // racecheck, that brings the single-rank run close to its timeout, so we skip these
    // CPU-only cases under racecheck.
    const char* const computeSanitizerTool = std::getenv("GMX_CI_COMPUTE_SANITIZER");
    if (getNumberOfTestMpiRanks() == 1 && computeSanitizerTool != nullptr
        && std::string(computeSanitizerTool) == "racecheck")
    {
        GTEST_SKIP() << "Under racecheck in the GROMACS CI system, the CPU-update cases only "
                        "run with multiple ranks";
    }

    runTest({});
}

INSTANTIATE_TEST_SUITE_P(
        VelocitiesConformToExpectations,
        VirtualSiteTest,
        // The bond type only matters for which forces are computed on the GPU, so
        // the CPU tests use only the default.
        ::testing::Combine(::testing::Values("vsite_test", "vsite_test_linear"),
                           ::testing::Values(InteractionFunction::GROMOS96Bonds),
                           ::testing::Values(true),
                           ::testing::Values("md", "md-vv", "sd", "bd"),
                           ::testing::Values("no", "v-rescale", "nose-hoover"),
                           ::testing::Values("no", "c-rescale", "parrinello-rahman")),
        sc_testNamer);

/*! \brief Test fixture for virtual sites constructed on the GPU
 *
 * Runs the same checks as VirtualSiteTest, but with the update and constraints
 * (and hence the virtual site construction) offloaded to the GPU.
 */
class VirtualSiteGpuTest : public VirtualSiteTest
{
};

TEST_P(VirtualSiteGpuTest, WithinToleranceOfReference)
{
    if (!GpuConfigurationCapabilities::UpdateVsite)
    {
        GTEST_SKIP() << "GPU update does not support virtual sites with this GPU backend.";
    }
    if (getCompatibleDevices(s_hwinfo->deviceInfoList).empty())
    {
        GTEST_SKIP() << "No compatible GPUs to test on.";
    }
    if (getNumberOfTestMpiRanks() > 1)
    {
        // With domain decomposition, the GPU update requires all virtual sites to be
        // within update groups. The test systems are not designed to allow update
        // groups to be used, so this can only be tested with a single rank. The test
        // binary is registered with CTest both with a single rank and with multiple
        // ranks, so this is still covered.
        GTEST_SKIP() << "Virtual sites with GPU update require a single rank for this system.";
    }

    const auto& [simulationName, bondType, testVirtualSitesN, integrator, tcoupling, pcoupling] =
            GetParam();
    if (integrator == "sd" && (tcoupling != "no" || pcoupling == "parrinello-rahman"))
    {
        GTEST_SKIP() << "sd handles temperature coupling implicitly, and the GPU update does not "
                        "support it with Parrinello-Rahman pressure coupling.";
    }

    // Constructing virtual sites on the GPU is experimental, and needs to be enabled explicitly.
    GmxEnvGuard gpuVsiteEnvGuard("GMX_GPU_VSITES", "1");

    prepareTpr();
    if (HasFatalFailure())
    {
        return;
    }

    // The test system uses plain cut-off electrostatics, so the non-bonded
    // interactions have to be offloaded for the GPU update to be possible.
    // The bonded interactions are required to run on the GPU too. With harmonic
    // bonds, all forces are then computed on the GPU. GROMOS96 bonds can not run
    // on GPUs, so they are then computed on the host, which tests that the host
    // forces are also spread from the virtual sites.
    runMdrunAndCheckVirtualSites(integrator,
                                 { { "-update", "gpu" }, { "-nb", "gpu" }, { "-bonded", "gpu" } });
    const std::string gpuUpdateTrajectoryName = runner_.fullPrecisionTrajectoryFileName_;

    // The above only checks frames individually, so we also compare the whole trajectory to
    // one using the CPU update. The non-bonded and bonded interactions stay on the GPU,
    // so that the main part of any difference stems from the update and virtual sites.
    runner_.fullPrecisionTrajectoryFileName_ =
            fileManager_.getTemporaryFilePath("cpu-update.trr").string();
    runMdrunAndCheckVirtualSites(integrator,
                                 { { "-update", "cpu" }, { "-nb", "gpu" }, { "-bonded", "gpu" } });

    TrajectoryFrameMatchSettings matchSettings;
    matchSettings.mustCompareBox        = true;
    matchSettings.coordinatesComparison = ComparisonConditions::CompareIfReferenceFound;
    matchSettings.velocitiesComparison  = ComparisonConditions::CompareIfReferenceFound;
    matchSettings.forcesComparison      = ComparisonConditions::CompareIfReferenceFound;
    // Set the tolerances for the trajectory comparison
    // Loosen the tolerance for the forces and esp. the SD integrator, which has larger rounding
    // differences between CPU and GPU implementations
    TrajectoryTolerances tolerances = TrajectoryComparison::s_defaultTrajectoryTolerances;
    tolerances.coordinates          = relativeToleranceAsFloatingPoint(1.0, 1.0e-5);
    tolerances.velocities = relativeToleranceAsFloatingPoint(1.0, integrator == "sd" ? 1.0e-3 : 5.0e-4);
    tolerances.forces = relativeToleranceAsFloatingPoint(100.0, integrator == "sd" ? 1.0e-2 : 5.0e-3);
    compareTrajectories(runner_.fullPrecisionTrajectoryFileName_,
                        gpuUpdateTrajectoryName,
                        TrajectoryComparison(matchSettings, tolerances));
}

INSTANTIATE_TEST_SUITE_P(
        VelocitiesConformToExpectations,
        VirtualSiteGpuTest,
        // GPU update only supports the linear virtual site types so far.
        // N-type virtual sites are not tested, as GPU update will support them later
        // than the other types; this can be changed once it does.
        ::testing::Combine(::testing::Values("vsite_test_linear"),
                           ::testing::Values(InteractionFunction::GROMOS96Bonds, InteractionFunction::Bonds),
                           ::testing::Values(false),
                           ::testing::Values("md", "sd"),
                           ::testing::Values("no", "v-rescale"),
                           ::testing::Values("no", "c-rescale", "parrinello-rahman")),
        sc_testNamer);

} // namespace
} // namespace test
} // namespace gmx
