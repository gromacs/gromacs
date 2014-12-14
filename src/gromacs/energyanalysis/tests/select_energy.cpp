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
 * Tests for analysis of energy files
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gromacs/energyanalysis/select_energy.h"

#include <cstdlib>

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/energyanalysis/analysismodule.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringstream.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/textblockmatchers.h"
#include "testutils/xvgtest.h"

namespace gmx
{

namespace energyanalysis
{

namespace
{

class SelectTest : public gmx::test::CommandLineTestBase
{
    public:
        SelectTest()
        {
            eNU[0] = { "Dog",  "Ear"   };
            eNU[1] = { "Cat",  "Tail"  };
            eNU[2] = { "Parrot",  "Beak"  };
            eNU[3] = { "  Pheasant",  "Paw"   };
        }
        //! Input strings
        std::array<EnergyNameUnit, 4> eNU;
};

TEST_F(SelectTest, Strings)
{
    std::array<const char *, 2> terms;
    terms[0] = { "Dog" };
    terms[1] = { "Cat" };
    std::vector<int>  set;
    StringInputStream input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(2, set.size());
    if (2 == set.size())
    {
        EXPECT_EQ(0, set[0]);
        EXPECT_EQ(1, set[1]);
    }
}

TEST_F(SelectTest, Numbers)
{
    std::array<const char *, 2> terms;
    terms[0] = { "2" };
    terms[1] = { "1" };
    std::vector<int>  set;
    StringInputStream input(terms);

    select_energies(eNU, false, &input, set);
    // Numbers that come back are indices in the array,
    // so one less than the input.
    EXPECT_EQ(2, set.size());
    if (2 == set.size())
    {
        EXPECT_EQ(1, set[0]);
        EXPECT_EQ(0, set[1]);
    }
}

TEST_F(SelectTest, NumberOutOfRange)
{
    std::array<const char *, 1> terms;
    terms[0] = { "7" };
    std::vector<int>            set;
    StringInputStream           input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(0, set.size());
}

TEST_F(SelectTest, NegativeNumber)
{
    std::array<const char *, 1> terms;
    terms[0] = { "-3" };
    std::vector<int>            set;
    StringInputStream           input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(0, set.size());
}

TEST_F(SelectTest, IncorrectString)
{
    std::array<const char *, 1> terms;
    terms[0] = { "koko" };
    std::vector<int>            set;
    StringInputStream           input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(0, set.size());
}

TEST_F(SelectTest, CorrectStringLeadingSpacesTarget)
{
    std::array<const char *, 1> terms;
    terms[0] = { "phea" };
    std::vector<int>            set;
    StringInputStream           input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(1, set.size());
}

TEST_F(SelectTest, CorrectStringLeadingSpacesKey)
{
    std::array<const char *, 1> terms;
    terms[0] = { " parr" };
    std::vector<int>            set;
    StringInputStream           input(terms);

    select_energies(eNU, false, &input, set);
    EXPECT_EQ(1, set.size());
}

}

} // namespace energyanalysis

} // namespace gmx
