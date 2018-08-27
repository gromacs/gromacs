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
/*! \internal \file
 * \brief
 * Tests whether the class AcceptOrRewind is working as it should.
 *
 * \author Sebastian Wingbermuehle
 */

#include "gmxpre.h"

#include "gromacs/hybridMCMD/acceptorrewind.h"

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"

namespace gmx
{

namespace test
{

//! Database of 51 water atom input positions (DIM reals per atom, taken from spc216.gro) for use as test inputs.
const real g_positions[] = {
    .130, -.041, -.291,
    .120, -.056, -.192,
    .044, -.005, -.327,
    -.854, -.406, .477,
    -.900, -.334, .425,
    -.858, -.386, .575,
    .351, -.061, .853,
    .401, -.147, .859,
    .416, .016, .850,
    -.067, -.796, .873,
    -.129, -.811, .797,
    -.119, -.785, .958,
    -.635, -.312, -.356,
    -.629, -.389, -.292,
    -.687, -.338, -.436,
    .321, -.919, .242,
    .403, -.880, .200,
    .294, -1.001, .193,
    -.404, .735, .728,
    -.409, .670, .803,
    -.324, .794, .741,
    .461, -.596, -.135,
    .411, -.595, -.221,
    .398, -.614, -.059,
    -.751, -.086, .237,
    -.811, -.148, .287,
    -.720, -.130, .152,
    .202, .285, -.364,
    .122, .345, -.377,
    .192, .236, -.278,
    -.230, -.485, .081,
    -.262, -.391, .071,
    -.306, -.548, .069,
    .464, -.119, .323,
    .497, -.080, .409,
    .540, -.126, .258,
    -.462, .107, .426,
    -.486, .070, .336,
    -.363, .123, .430,
    .249, -.077, -.621,
    .306, -.142, -.571,
    .233, -.110, -.714,
    -.922, -.164, .904,
    -.842, -.221, .925,
    -.971, -.204, .827,
    .382, .700, .480,
    .427, .610, .477,
    .288, .689, .513,
    .781, .264, -.113,
    .848, .203, -.070,
    .708, .283, -.048
};

//! Dummy Metropolis step that returns the boolean required by AcceptOrRewind
class MetropolisStepDummy : public IMetropolisStep
{
    public:
        //! Constructor
        MetropolisStepDummy()
        {
        }

        bool accept(const int64_t step, const gmx_enerdata_t *enerd) override
        {
            if (step == 2)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        //! Destructor
        virtual ~MetropolisStepDummy() override {}
};

class AcceptOrRewindTest : public ::testing::Test
{
    public:
        //! Constructor
        AcceptOrRewindTest() : acceptOrRewind_(&metropolisStepDummy_, &localState_, &commrec_)
        {
            // read in water atoms into localState_.x
            int                    i = 0;
            gmx::BasicVector<real> r;
            while (i < 255)
            {
                r[i % DIM] = g_positions[i];
                if (i % DIM == DIM - 1)
                {
                    localState_.x.push_back(r);
                }
                i++;
            }
        }

        void updatePositions()
        {
            // Perturb the atom positions, to appear like an "update"
            unsigned int i, j;
            for (i = 0; i < localState_.x.size(); i++)
            {
                for (j = 0; j < DIM; j++)
                {
                    localState_.x[i][j] += 0.5;
                }
            }
        }

        bool checkVectors()
        {
            bool         isEqual   = true;
            real         tolerance = 0.001;
            unsigned int i, j;
            for (i = 0; i < localState_.x.size(); i++)
            {
                if (!isEqual)
                {
                    break;
                }
                for (j = 0; j < DIM; j++)
                {
                    if ((localState_.x[i][j] < xRef_[i][j] - tolerance) || (localState_.x[i][j] > xRef_[i][j] + tolerance))
                    {
                        isEqual = false;
                        break;
                    }
                }
            }
            return isEqual;
        }

        void runTest(const int64_t step)
        {
            // step 0 is needed to back-up coordinates
            if (step > 0)
            {
                updatePositions();
            }
            // all steps but step 2 are accepted => update reference vector accordingly
            if (step < 2)
            {
                xRef_ = localState_.x;
            }
            accepted_ = acceptOrRewind_.updateTState(step, &localState_, nullptr, nullptr, &commrec_, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
            // compare vectors
            EXPECT_TRUE(checkVectors());
        }

    private:
        MetropolisStepDummy         metropolisStepDummy_;
        t_state                     localState_;
        t_commrec                   commrec_;
        AcceptOrRewind              acceptOrRewind_;
        bool                        accepted_;
        gmx::HostVector<gmx::RVec>  xRef_;
};

TEST_F(AcceptOrRewindTest, AcceptOrRewindWithoutDDWorks)
{
    runTest(0);
    runTest(1);
    runTest(2);
}

} // namespace
} // namespace
