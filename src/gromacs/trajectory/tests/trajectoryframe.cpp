/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for t_trxframe and supporting functions.
 *
 * \ingroup module_trajectory
 */

#include "gmxpre.h"

#include "gromacs/trajectory/trajectoryframe.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

// How do I set the working directory or input path appropriately?
/*! \internal \brief input file path
 *
 * \{
 */
const std::string filepath = std::string(SOURCEDIR) + "/simple.gro";
const char* const trjfile  = filepath.c_str();
/*! \}
 */

/// trxframe test case
TEST(TrxFrameTest, Construct)
{
    t_trxframe f;
    clear_trxframe(&f, true);
    ASSERT_EQ(f.x, nullptr);
}

/// read a frame for reference
TEST(TrxFrameTest, Read)
{
    t_trxstatus *status {
        nullptr
    };
    t_trxframe *fr;
    snew(fr, 1);

    time_unit_t time_unit {};
    //    = static_cast<time_unit_t>(settings_.timeUnit() + 1);

    gmx_output_env_t *oenv {
        nullptr
    };

    output_env_init(&oenv, gmx::getProgramContext(), time_unit, FALSE, exvgNONE, 0);

    int frflags = TRX_NEED_X;

    // Check succesfful read
    ASSERT_TRUE(read_first_frame(oenv, &status, trjfile, fr, frflags));
    ASSERT_NE(fr->x, nullptr);

    close_trj(status);
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}

/// make a deep copy of a frame
TEST(TrxFrameTest, Copy)
{
    t_trxstatus *status {
        nullptr
    };
    t_trxframe *fr;
    snew(fr, 1);

    time_unit_t time_unit {};
    //    = static_cast<time_unit_t>(settings_.timeUnit() + 1);

    gmx_output_env_t *oenv {
        nullptr
    };

    output_env_init(&oenv, gmx::getProgramContext(), time_unit, FALSE, exvgNONE, 0);

    int frflags = TRX_NEED_X;

    // Check succesfful read
    ASSERT_TRUE(read_first_frame(oenv, &status, trjfile, fr, frflags));
    ASSERT_NE(fr->x, nullptr);

    auto fr_copy = gmx::trajectory::trxframe_copy(*fr);
    ASSERT_NE(fr_copy->x, nullptr);
    ASSERT_NE(fr_copy->x, fr->x);
    ASSERT_NE(fr_copy->x[0][0], 0);
    ASSERT_EQ(fr_copy->x[0][0], fr->x[0][0]);

    close_trj(status);
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}

/// test deleter
/// not really a good way to test...
TEST(TrxFrameTest, Release)
{
    t_trxstatus *status {
        nullptr
    };
    t_trxframe *fr;
    snew(fr, 1);

    time_unit_t time_unit {};
    //    = static_cast<time_unit_t>(settings_.timeUnit() + 1);

    gmx_output_env_t *oenv {
        nullptr
    };

    output_env_init(&oenv, gmx::getProgramContext(), time_unit, FALSE, exvgNONE, 0);

    int frflags = TRX_NEED_X;

    // Check succesfful read
    ASSERT_TRUE(read_first_frame(oenv, &status, trjfile, fr, frflags));
    ASSERT_NE(fr->x, nullptr);

    std::weak_ptr<t_trxframe> fr_dangler;
    {
        // This copies to fr_tmp and transfers ownership to fr_copy, but that
        // doesn't prove much other than that such operations compile and run.
        auto fr_tmp = gmx::trajectory::trxframe_copy(*fr);
        std::shared_ptr<t_trxframe> fr_copy {
            std::move(fr_tmp)
        };
        // If trxframe_copy returns unique_ptr, fr_tmp should be empty now.
        ASSERT_EQ(fr_tmp.get(), nullptr);
        ASSERT_NE(fr_copy->x, nullptr);
        ASSERT_NE(fr_copy->x, fr->x);
        ASSERT_NE(fr_copy->x[0][0], 0);
        ASSERT_EQ(fr_copy->x[0][0], fr->x[0][0]);
        fr_dangler = fr_copy;
    }
    ASSERT_EQ(fr_dangler.use_count(), 0);

    close_trj(status);
    sfree(fr->x);
    sfree(fr->v);
    sfree(fr->f);
    sfree(fr->index);
    sfree(fr);
}
