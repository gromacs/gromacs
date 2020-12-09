/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 *
 * \brief Tests routines in checkpoint.h .
 *
 * \author Christian Blau <blau@kth.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/checkpoint.h"

#include <gtest/gtest.h>

#include "gromacs/utility/exceptions.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

TEST(Checkpoint, ReadingThrowsWhenValueNotPresent)
{
    KeyValueTreeObject kvtObject;
    std::int64_t       readValue = 0;
    EXPECT_THROW_GMX(readKvtCheckpointValue(compat::make_not_null(&readValue), "non", "sense", kvtObject),
                     gmx::InternalError);
}

TEST(Checkpoint, ReadingDoesNotThrowWhenValuePresent)
{
    int64_t             value      = 37;
    std::string         name       = "checkpointedInteger";
    std::string         identifier = "testingmodule";
    KeyValueTreeBuilder kvtBuilder;
    writeKvtCheckpointValue(value, name, identifier, kvtBuilder.rootObject());
    const auto kvtObject = kvtBuilder.build();
    int64_t    readValue = 0;
    EXPECT_NO_THROW_GMX(
            readKvtCheckpointValue(compat::make_not_null(&readValue), name, identifier, kvtObject));
}

TEST(Checkpoint, KvtRoundTripInt64)
{
    int64_t             value      = INT64_MAX;
    std::string         name       = "checkpointedInteger";
    std::string         identifier = "testingmodule";
    KeyValueTreeBuilder kvtBuilder;
    writeKvtCheckpointValue(value, name, identifier, kvtBuilder.rootObject());
    const auto kvtObject = kvtBuilder.build();
    int64_t    readValue = 0;
    readKvtCheckpointValue(compat::make_not_null(&readValue), name, identifier, kvtObject);
    EXPECT_EQ(value, readValue);
}

TEST(Checkpoint, KvtRoundTripReal)
{
    real                value      = 926.7;
    std::string         name       = "checkpointedReal";
    std::string         identifier = "testingmodule";
    KeyValueTreeBuilder kvtBuilder;
    writeKvtCheckpointValue(value, name, identifier, kvtBuilder.rootObject());
    const auto kvtObject = kvtBuilder.build();
    real       readValue = 0;
    readKvtCheckpointValue(compat::make_not_null(&readValue), name, identifier, kvtObject);
    EXPECT_EQ(value, readValue);
}


} // namespace
} // namespace test
} // namespace gmx
