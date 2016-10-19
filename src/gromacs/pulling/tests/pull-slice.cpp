/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements test of some pulling routines
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_pulling
 */
#include "gmxpre.h"

#include <cmath>

#include <algorithm>

#include <gtest/gtest.h>

#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace
{

using gmx::test::defaultRealTolerance;

class SlicedPullTest : public ::testing::Test
{
    protected:
        SlicedPullTest() {}

        void test(rvec x[], int nat, int ind[], int nat_grp, int nat_slice)
        {
            GMX_ASSERT(nat >= nat_grp, "The number of atoms inside the group has to be smaller than the total number of atoms");

            {
                pull_t pull;
                pull.ngroup = 2;
                snew(pull.group, pull.ngroup);

                pull_group_work_t *pgrp_sliced = &pull.group[0];
                pgrp_sliced->params.bSliced     = true;
                pgrp_sliced->params.slice_x_min = 0.1;
                pgrp_sliced->params.slice_x_max = 0.2;
                pgrp_sliced->nat_loc            = nat_grp;
                snew(pgrp_sliced->ind_loc, pgrp_sliced->nat_loc);
                for (int i = 0; i < nat_grp; i++)
                {
                    pgrp_sliced->ind_loc[i] = ind[i];
                }
                pgrp_sliced->weight_loc = NULL;

                pgrp_sliced->nat_loc_sliced = pgrp_sliced->nat_loc;
                snew(pgrp_sliced->ind_loc_sliced, pgrp_sliced->nat_loc_sliced);
                memcpy(pgrp_sliced->ind_loc_sliced, pgrp_sliced->ind_loc, pgrp_sliced->nat_loc_sliced);
                snew(pgrp_sliced->weight_loc_sliced, pgrp_sliced->nat_loc_sliced);

                pull_group_work_t *pgrp_nonsliced = &pull.group[1];
                pgrp_nonsliced->params.bSliced = false;
                pgrp_nonsliced->nat_loc        = nat_grp;
                snew(pgrp_nonsliced->ind_loc, pgrp_nonsliced->nat_loc);
                for (int i = 0; i < nat_grp; i++)
                {
                    pgrp_nonsliced->ind_loc[i] = ind[i];
                }
                pgrp_nonsliced->weight_loc = NULL;

                int nat_nonslice = pgrp_nonsliced->nat_loc_sliced;

                update_sliced_pull_groups(NULL, &pull, x);
                EXPECT_EQ(pull.group[0].nat_loc_sliced, nat_slice);
                EXPECT_EQ(pull.group[1].nat_loc_sliced, nat_nonslice);

                int   nat_grp_sliced, nalloc_grp_sliced, nat_grp_nonsliced, nalloc_grp_nonsliced;
                int  *ind_grp_sliced, *ind_grp_nonsliced;
                real *weight_grp_sliced, *weight_grp_nonsliced;

                pull_get_indices_weights(pgrp_sliced, nat_grp_sliced, nalloc_grp_sliced, ind_grp_sliced, weight_grp_sliced);
                EXPECT_EQ(pgrp_sliced->nat_loc_sliced, nat_grp_sliced);
                EXPECT_EQ(pgrp_sliced->nalloc_loc_sliced, nalloc_grp_sliced);
                EXPECT_EQ(pgrp_sliced->ind_loc_sliced, ind_grp_sliced);
                EXPECT_EQ(pgrp_sliced->weight_loc_sliced, weight_grp_sliced);

                pull_get_indices_weights(pgrp_nonsliced, nat_grp_nonsliced, nalloc_grp_nonsliced, ind_grp_nonsliced, weight_grp_nonsliced);
                EXPECT_EQ(pgrp_nonsliced->nat_loc, nat_grp_nonsliced);
                EXPECT_EQ(pgrp_nonsliced->nalloc_loc, nalloc_grp_nonsliced);
                EXPECT_EQ(pgrp_nonsliced->ind_loc, ind_grp_nonsliced);
                EXPECT_EQ(pgrp_nonsliced->weight_loc, weight_grp_nonsliced);
            }
        }
};

TEST_F (SlicedPullTest, 2of3AtomsInSlice)
{
    int  nat = 5;
    rvec x[] = { { 0.15, 0.09, 0.68 },
                 { 0.17, 0.85, 0.37 },
                 { 0.25, 0.19, 0.10 },
                 { 0.35, 0.15, 0.04 },
                 { 0.16, 0.14, 0.26 } };
    int  nat_grp   = 3;
    int  ind[]     = { 1, 2, 4 };
    int  nat_slice = 2;

    test(x, nat, ind, nat_grp, nat_slice);
}


TEST_F (SlicedPullTest, 4of4AtomsInSlice)
{
    int  nat = 5;
    rvec x[] = { { 0.11, 0.84, 0.78 },
                 { 0.16, 0.53, 0.68 },
                 { 0.19, 0.80, 0.63 },
                 { 0.18, 0.56, 0.11 },
                 { 0.12, 0.60, 0.35 } };
    int  nat_grp   = 4;
    int  ind[]     = { 1, 2, 3, 4 };
    int  nat_slice = 4;

    test(x, nat, ind, nat_grp, nat_slice);
}

}

}
