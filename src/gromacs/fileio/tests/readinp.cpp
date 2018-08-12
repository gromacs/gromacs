/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 * Tests utilities for routines that parse fields e.g. from grompp input
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/fileio/readinp.h"

#include <gtest/gtest.h>

#include "gromacs/fileio/warninp.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{
namespace testing
{

class ReadTest : public ::testing::Test
{
    public:
        ReadTest() : inputField_ {{(t_inpfile(0, 0, false, false, false, "test", ""))}},
        wi_()

        {
            wi_ = init_warning(FALSE, 0);
            wiGuard_.reset(wi_);
        }

        std::vector<t_inpfile>                         inputField_;
        warninp_t                                      wi_;
        gmx::unique_cptr<struct warninp, free_warning> wiGuard_;
};

TEST_F(ReadTest, get_eint_ReadsInteger)
{
    inputField_.front().value_.assign("1");
    ASSERT_EQ(1, get_eint(&inputField_, "test", 2, wi_));
    ASSERT_FALSE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_eint_WarnsAboutFloat)
{
    inputField_.front().value_.assign("0.8");
    get_eint(&inputField_, "test", 2, wi_);
    ASSERT_TRUE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_eint_WarnsAboutString)
{
    inputField_.front().value_.assign("hello");
    get_eint(&inputField_, "test", 2, wi_);
    ASSERT_TRUE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_eint64_ReadsInteger)
{
    inputField_.front().value_.assign("1");
    ASSERT_EQ(1, get_eint64(&inputField_, "test", 2, wi_));
    ASSERT_FALSE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_eint64_WarnsAboutFloat)
{
    inputField_.front().value_.assign("0.8");
    get_eint64(&inputField_, "test", 2, wi_);
    ASSERT_TRUE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_eint64_WarnsAboutString)
{
    inputField_.front().value_.assign("hello");
    get_eint64(&inputField_, "test", 2, wi_);
    ASSERT_TRUE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_ereal_ReadsInteger)
{
    inputField_.front().value_.assign("1");
    ASSERT_EQ(1, get_ereal(&inputField_, "test", 2, wi_));
    ASSERT_FALSE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_ereal_ReadsFloat)
{
    inputField_.front().value_.assign("0.8");
    ASSERT_EQ(0.8, get_ereal(&inputField_, "test", 2, wi_));
    ASSERT_FALSE(warning_errors_exist(wi_));
}

TEST_F(ReadTest, get_ereal_WarnsAboutString)
{
    inputField_.front().value_.assign("hello");
    get_ereal(&inputField_, "test", 2, wi_);
    ASSERT_TRUE(warning_errors_exist(wi_));
}

}  // namespace testing
}  // namespace gmx
