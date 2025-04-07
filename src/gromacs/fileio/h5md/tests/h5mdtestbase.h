/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief Test helper class for creating HDF5 files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_H5MD_TEST_BASE_H
#define GMX_FILEIO_H5MD_TEST_BASE_H

#include <hdf5.h>

#include <memory>

#include <gtest/gtest.h>

namespace gmx
{
namespace test
{

//! \brief Base class for creating an H5md file in write mode for a test suite.
class H5mdTestBase : public ::testing::Test
{
public:
    H5mdTestBase();
    ~H5mdTestBase() override;

    //! \brief Return a handle to the open H5md file.
    hid_t fileid();

private:
    class Impl;

    //! Handle to implementation object.
    std::unique_ptr<Impl> impl_;
};

} // namespace test
} // namespace gmx

#endif // GMX_FILEIO_H5MD_TEST_BASE_H
