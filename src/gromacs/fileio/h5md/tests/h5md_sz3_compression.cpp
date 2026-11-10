/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Tests for SZ3-specific compression routines.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include <hdf5.h>

#if GMX_HAVE_SZ3
#    include <H5Z_SZ3.hpp>
#endif

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/exceptions.h"
#include "gromacs/fileio/h5md/h5md_compression.h"
#include "gromacs/fileio/h5md/h5md_guard.h"

#if !GMX_HAVE_SZ3
#    include "gromacs/utility/exceptions.h"

#    include "testutils/testasserts.h"
#endif

namespace gmx
{
namespace test
{
namespace
{

#if GMX_HAVE_SZ3

TEST(H5mdSz3CompressionTest, FilterIsApplied)
{
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));

    setSz3Compression(creationPropertyList, 5e-4);
    EXPECT_EQ(H5Pget_nfilters(creationPropertyList), 1);
    // H5Pget_filter_by_id2 returns: >=0 if the filter (second argument) is set, else <0
    EXPECT_GE(H5Pget_filter_by_id2(
                      creationPropertyList, H5Z_FILTER_SZ3, nullptr, nullptr, nullptr, 0, nullptr, nullptr),
              0);
}

TEST(H5mdSz3CompressionTest, AbsoluteErrorIsSet)
{
    constexpr double absoluteError = 5e-4;

    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));
    setSz3Compression(creationPropertyList, absoluteError);

    SZ3::Config conf;
    get_SZ3_conf_from_H5(creationPropertyList, conf);
    EXPECT_DOUBLE_EQ(absoluteError, conf.absErrorBound);
}

TEST(H5mdSz3CompressionTest, ThrowsForNonPositiveAbsoluteError)
{
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));

    EXPECT_THROW(setSz3Compression(creationPropertyList, -5e-4), H5mdError);
    EXPECT_THROW(setSz3Compression(creationPropertyList, 0.0), H5mdError);
}

TEST(H5mdSz3CompressionTest, ThrowsForInvalidHandles)
{
    const auto [accessPropertyList, accessPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_ACCESS));
    EXPECT_THROW(setSz3Compression(accessPropertyList, 5e-4), H5mdError)
            << "Must throw for handle which is not a H5P_DATASET_CREATE property list";
    EXPECT_THROW(setSz3Compression(H5I_INVALID_HID, 5e-4), H5mdError)
            << "Must throw for invalid handle";
}

#else

TEST(H5mdSz3DisabledTest, SetSz3CompressionThrows)
{
    const auto [creationPropertyList, creationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_DATASET_CREATE));
    EXPECT_THROW_GMX(setSz3Compression(creationPropertyList, 5e-4), NotImplementedError);
}

#endif // GMX_HAVE_SZ3

} // namespace
} // namespace test
} // namespace gmx
