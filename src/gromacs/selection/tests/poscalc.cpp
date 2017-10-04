/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests the position mapping engine.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/selection/poscalc.h"

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"

#include "toputils.h"

namespace
{

/********************************************************************
 * PositionCalculationTest
 */

class PositionCalculationTest : public ::testing::Test
{
    public:
        PositionCalculationTest();
        ~PositionCalculationTest();

        void generateCoordinates();

        gmx_ana_poscalc_t *createCalculation(e_poscalc_t type, int flags);
        void setMaximumGroup(gmx_ana_poscalc_t *pc, const gmx::ArrayRef<const int> &atoms);
        gmx_ana_pos_t *initPositions(gmx_ana_poscalc_t *pc, const char *name);

        void checkInitialized();
        void updateAndCheck(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p,
                            const gmx::ArrayRef<const int> &atoms,
                            gmx::test::TestReferenceChecker *checker,
                            const char *name);

        void testSingleStatic(e_poscalc_t type, int flags, bool bExpectTop,
                              const gmx::ArrayRef<const int> &atoms,
                              const gmx::ArrayRef<const int> &index = gmx::EmptyArrayRef());
        void testSingleDynamic(e_poscalc_t type, int flags, bool bExpectTop,
                               const gmx::ArrayRef<const int> &initAtoms,
                               const gmx::ArrayRef<const int> &evalAtoms,
                               const gmx::ArrayRef<const int> &index = gmx::EmptyArrayRef());

        gmx::test::TestReferenceData        data_;
        gmx::test::TestReferenceChecker     checker_;
        gmx::test::TopologyManager          topManager_;
        gmx::PositionCalculationCollection  pcc_;

    private:
        typedef std::unique_ptr<gmx_ana_pos_t> PositionPointer;

        struct PositionTest
        {
            PositionTest(PositionPointer pos, gmx_ana_poscalc_t *pc,
                         const char *name)
                : pos(std::move(pos)), pc(pc), name(name)
            {
            }

            // Default move constructor and assignment. Only needed for old compilers.
            PositionTest(PositionTest &&o)
                : pos(std::move(o.pos)), pc(o.pc), name(o.name)
            {
            }

            PositionTest &operator= (PositionTest &&o)
            {
                pos  = std::move(o.pos);
                pc   = o.pc;
                name = o.name;
                return *this;
            }

            PositionPointer                 pos;
            gmx_ana_poscalc_t              *pc;
            const char                     *name;
        };

        typedef std::vector<PositionTest> PositionTestList;

        void setTopologyIfRequired();
        void checkPositions(gmx::test::TestReferenceChecker *checker,
                            const char *name, gmx_ana_pos_t *p,
                            bool bCoordinates);

        std::vector<gmx_ana_poscalc_t *>    pcList_;
        PositionTestList                    posList_;
        bool                                bTopSet_;
};

PositionCalculationTest::PositionCalculationTest()
    : checker_(data_.rootChecker()), bTopSet_(false)
{
    topManager_.requestFrame();
}

PositionCalculationTest::~PositionCalculationTest()
{
    std::vector<gmx_ana_poscalc_t *>::reverse_iterator pci;
    for (pci = pcList_.rbegin(); pci != pcList_.rend(); ++pci)
    {
        gmx_ana_poscalc_free(*pci);
    }
}

void PositionCalculationTest::generateCoordinates()
{
    t_atoms    &atoms = topManager_.atoms();
    t_trxframe *frame = topManager_.frame();
    for (int i = 0; i < atoms.nr; ++i)
    {
        frame->x[i][XX] = i;
        frame->x[i][YY] = atoms.atom[i].resind;
        frame->x[i][ZZ] = 0.0;
        if (frame->bV)
        {
            copy_rvec(frame->x[i], frame->v[i]);
            frame->v[i][ZZ] = 1.0;
        }
        if (frame->bF)
        {
            copy_rvec(frame->x[i], frame->f[i]);
            frame->f[i][ZZ] = -1.0;
        }
    }
}

gmx_ana_poscalc_t *
PositionCalculationTest::createCalculation(e_poscalc_t type, int flags)
{
    pcList_.reserve(pcList_.size() + 1);
    pcList_.push_back(pcc_.createCalculation(type, flags));
    return pcList_.back();
}

void PositionCalculationTest::setMaximumGroup(gmx_ana_poscalc_t              *pc,
                                              const gmx::ArrayRef<const int> &atoms)
{
    setTopologyIfRequired();
    gmx_ana_index_t g;
    g.isize = atoms.size();
    g.index = const_cast<int *>(atoms.data());
    gmx_ana_poscalc_set_maxindex(pc, &g);
}

gmx_ana_pos_t *
PositionCalculationTest::initPositions(gmx_ana_poscalc_t *pc, const char *name)
{
    posList_.reserve(posList_.size() + 1);
    PositionPointer p(new gmx_ana_pos_t());
    gmx_ana_pos_t  *result = p.get();
    posList_.emplace_back(std::move(p), pc, name);
    gmx_ana_poscalc_init_pos(pc, result);
    return result;
}

void PositionCalculationTest::checkInitialized()
{
    gmx::test::TestReferenceChecker  compound(
            checker_.checkCompound("InitializedPositions", nullptr));
    PositionTestList::const_iterator pi;
    for (pi = posList_.begin(); pi != posList_.end(); ++pi)
    {
        checkPositions(&compound, pi->name, pi->pos.get(), false);
    }
}

void PositionCalculationTest::updateAndCheck(
        gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p, const gmx::ArrayRef<const int> &atoms,
        gmx::test::TestReferenceChecker *checker, const char *name)
{
    gmx_ana_index_t g;
    g.isize = atoms.size();
    g.index = const_cast<int *>(atoms.data());
    gmx_ana_poscalc_update(pc, p, &g, topManager_.frame(), nullptr);
    checkPositions(checker, name, p, true);
}

void PositionCalculationTest::testSingleStatic(
        e_poscalc_t type, int flags, bool bExpectTop,
        const gmx::ArrayRef<const int> &atoms,
        const gmx::ArrayRef<const int> &index)
{
    t_trxframe *frame = topManager_.frame();
    if (frame->bV)
    {
        flags |= POS_VELOCITIES;
    }
    if (frame->bF)
    {
        flags |= POS_FORCES;
    }
    gmx_ana_poscalc_t *pc = createCalculation(type, flags);
    const bool         requiresTopology
        = gmx_ana_poscalc_required_topology_info(pc)
            != gmx::PositionCalculationCollection::RequiredTopologyInfo::None;
    EXPECT_EQ(bExpectTop, requiresTopology);
    setMaximumGroup(pc, atoms);
    gmx_ana_pos_t *p = initPositions(pc, nullptr);
    checkInitialized();
    {
        generateCoordinates();
        if (!index.empty())
        {
            topManager_.initFrameIndices(index);
        }
        pcc_.initEvaluation();
        pcc_.initFrame(frame);
        gmx::test::TestReferenceChecker frameCompound(
                checker_.checkCompound("EvaluatedPositions", "Frame0"));
        updateAndCheck(pc, p, atoms, &frameCompound, nullptr);
    }
}

void PositionCalculationTest::testSingleDynamic(
        e_poscalc_t type, int flags, bool bExpectTop,
        const gmx::ArrayRef<const int> &initAtoms,
        const gmx::ArrayRef<const int> &evalAtoms,
        const gmx::ArrayRef<const int> &index)
{
    gmx_ana_poscalc_t *pc = createCalculation(type, flags | POS_DYNAMIC);
    const bool         requiresTopology
        = gmx_ana_poscalc_required_topology_info(pc)
            != gmx::PositionCalculationCollection::RequiredTopologyInfo::None;
    EXPECT_EQ(bExpectTop, requiresTopology);
    setMaximumGroup(pc, initAtoms);
    gmx_ana_pos_t *p = initPositions(pc, nullptr);
    checkInitialized();
    {
        generateCoordinates();
        if (!index.empty())
        {
            topManager_.initFrameIndices(index);
        }
        pcc_.initEvaluation();
        pcc_.initFrame(topManager_.frame());
        gmx::test::TestReferenceChecker frameCompound(
                checker_.checkCompound("EvaluatedPositions", "Frame0"));
        updateAndCheck(pc, p, evalAtoms, &frameCompound, nullptr);
    }
}

void PositionCalculationTest::setTopologyIfRequired()
{
    if (bTopSet_)
    {
        return;
    }
    std::vector<gmx_ana_poscalc_t *>::const_iterator pci;
    for (pci = pcList_.begin(); pci != pcList_.end(); ++pci)
    {
        const bool         requiresTopology
            = gmx_ana_poscalc_required_topology_info(*pci)
                != gmx::PositionCalculationCollection::RequiredTopologyInfo::None;
        if (requiresTopology)
        {
            bTopSet_ = true;
            pcc_.setTopology(topManager_.topology());
            return;
        }
    }
}

void PositionCalculationTest::checkPositions(
        gmx::test::TestReferenceChecker *checker,
        const char *name, gmx_ana_pos_t *p, bool bCoordinates)
{
    gmx::test::TestReferenceChecker compound(
            checker->checkCompound("Positions", name));
    compound.checkInteger(p->count(), "Count");
    const char *type = "???";
    switch (p->m.type)
    {
        case INDEX_UNKNOWN: type = "unknown";   break;
        case INDEX_ATOM:    type = "atoms";     break;
        case INDEX_RES:     type = "residues";  break;
        case INDEX_MOL:     type = "molecules"; break;
        case INDEX_ALL:     type = "single";    break;
    }
    compound.checkString(type, "Type");
    compound.checkSequenceArray(p->count() + 1, p->m.mapb.index, "Block");
    for (int i = 0; i < p->count(); ++i)
    {
        gmx::test::TestReferenceChecker posCompound(
                compound.checkCompound("Position", nullptr));
        posCompound.checkSequence(&p->m.mapb.a[p->m.mapb.index[i]],
                                  &p->m.mapb.a[p->m.mapb.index[i+1]],
                                  "Atoms");
        posCompound.checkInteger(p->m.refid[i], "RefId");
        if (bCoordinates)
        {
            posCompound.checkVector(p->x[i], "Coordinates");
        }
        if (bCoordinates && p->v != nullptr)
        {
            posCompound.checkVector(p->v[i], "Velocity");
        }
        if (bCoordinates && p->f != nullptr)
        {
            posCompound.checkVector(p->f[i], "Force");
        }
        int originalIdIndex = (p->m.refid[i] != -1 ? p->m.refid[i] : i);
        EXPECT_EQ(p->m.orgid[originalIdIndex], p->m.mapid[i]);
    }
}

/********************************************************************
 * Actual tests
 */

TEST_F(PositionCalculationTest, ComputesAtomPositions)
{
    const int group[] = { 0, 1, 2, 3 };
    topManager_.requestVelocities();
    topManager_.requestForces();
    topManager_.initAtoms(4);
    testSingleStatic(POS_ATOM, 0, false, group);
}

TEST_F(PositionCalculationTest, ComputesResidueCOGPositions)
{
    const int group[] = { 0, 1, 2, 3, 4, 8 };
    topManager_.requestVelocities();
    topManager_.requestForces();
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);
    testSingleStatic(POS_RES, 0, true, group);
}

TEST_F(PositionCalculationTest, ComputesResidueCOMPositions)
{
    const int group[] = { 0, 1, 2, 3, 4, 8 };
    topManager_.requestVelocities();
    topManager_.requestForces();
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);
    testSingleStatic(POS_RES, POS_MASS, true, group);
}

TEST_F(PositionCalculationTest, ComputesGroupCOGPositions)
{
    const int group[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    topManager_.requestVelocities();
    topManager_.requestForces();
    topManager_.initAtoms(9);
    // Topology (masses) is requires for computing the force
    testSingleStatic(POS_ALL, 0, true, group);
}

TEST_F(PositionCalculationTest, ComputesGroupCOMPositions)
{
    const int group[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    topManager_.requestVelocities();
    topManager_.requestForces();
    topManager_.initAtoms(9);
    testSingleStatic(POS_ALL, POS_MASS, true, group);
}

TEST_F(PositionCalculationTest, ComputesPositionsWithCompleteWhole)
{
    const int group[] = { 0, 1, 2, 3, 4, 8 };
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);
    testSingleStatic(POS_RES, POS_COMPLWHOLE, true, group);
}

TEST_F(PositionCalculationTest, ComputesPositionsWithCompleteMax)
{
    const int maxGroup[]  = { 0, 1, 4, 5, 6, 8 };
    const int evalGroup[] = { 0, 1, 5, 6 };
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);
    testSingleDynamic(POS_RES, POS_COMPLMAX, true, maxGroup, evalGroup);
}

TEST_F(PositionCalculationTest, ComputesPositionMask)
{
    const int maxGroup[]  = { 0, 1, 2, 3, 4, 5 };
    const int evalGroup[] = { 1, 2, 4 };
    topManager_.initAtoms(6);
    testSingleDynamic(POS_ATOM, POS_MASKONLY, false, maxGroup, evalGroup);
}

// TODO: Check for POS_ALL_PBC

TEST_F(PositionCalculationTest, HandlesFramesWithLessAtoms)
{
    const int group[] = { 2, 3, 5, 6 };
    topManager_.initAtoms(10);
    testSingleStatic(POS_ATOM, 0, false, group, group);
}

TEST_F(PositionCalculationTest, HandlesFramesWithLessAtoms2)
{
    const int group[] = { 2, 3, 6, 7 };
    const int index[] = { 1, 2, 3, 4, 6, 7 };
    topManager_.initAtoms(10);
    testSingleStatic(POS_ATOM, 0, false, group, index);
}

TEST_F(PositionCalculationTest, HandlesIdenticalStaticCalculations)
{
    const int group[] = { 0, 1, 4, 5, 6, 7 };
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);

    gmx_ana_poscalc_t *pc1 = createCalculation(POS_RES, 0);
    gmx_ana_poscalc_t *pc2 = createCalculation(POS_RES, 0);
    gmx_ana_poscalc_t *pc3 = createCalculation(POS_RES, 0);
    setMaximumGroup(pc1, group);
    setMaximumGroup(pc2, group);
    setMaximumGroup(pc3, group);
    gmx_ana_pos_t *p1 = initPositions(pc1, "Positions");
    gmx_ana_pos_t *p2 = initPositions(pc2, "Positions");
    gmx_ana_pos_t *p3 = initPositions(pc3, "Positions");
    checkInitialized();
    {
        generateCoordinates();
        pcc_.initEvaluation();
        pcc_.initFrame(topManager_.frame());
        gmx::test::TestReferenceChecker frameCompound(
                checker_.checkCompound("EvaluatedPositions", "Frame0"));
        updateAndCheck(pc1, p1, group, &frameCompound, "Positions");
        updateAndCheck(pc2, p2, group, &frameCompound, "Positions");
        updateAndCheck(pc3, p3, group, &frameCompound, "Positions");
    }
}

TEST_F(PositionCalculationTest, HandlesOverlappingStaticCalculations)
{
    const int group1[] = { 0, 1, 4, 5 };
    const int group2[] = { 4, 5, 7, 8 };
    topManager_.initAtoms(9);
    topManager_.initUniformResidues(3);

    gmx_ana_poscalc_t *pc1 = createCalculation(POS_RES, 0);
    gmx_ana_poscalc_t *pc2 = createCalculation(POS_RES, 0);
    setMaximumGroup(pc1, group1);
    setMaximumGroup(pc2, group2);
    gmx_ana_pos_t *p1 = initPositions(pc1, "P1");
    gmx_ana_pos_t *p2 = initPositions(pc2, "P2");
    checkInitialized();
    {
        generateCoordinates();
        pcc_.initEvaluation();
        pcc_.initFrame(topManager_.frame());
        gmx::test::TestReferenceChecker frameCompound(
                checker_.checkCompound("EvaluatedPositions", "Frame0"));
        updateAndCheck(pc1, p1, group1, &frameCompound, "P1");
        updateAndCheck(pc2, p2, group2, &frameCompound, "P2");
    }
}

// TODO: Check for handling of more multiple calculation cases

} // namespace
