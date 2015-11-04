/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
   \brief
   Tests for class shared_string and helper classes

   From the stringutil tests as template.
   For development, the tests can be run with a '-stdout' command-line option
   to print out the help to stdout instead of using the XML reference
   framework.

   \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>

   \ingroup module_utility
 */
#include "gmxpre.h"

#include <fstream>
#include <iomanip>  // ::std::setw, ::std::setfill
#include <ios>      // ::std::left, ::std::rightusing debug;
#include <iostream> // ::debug
#include <ostream>  // ::std::ostream (<<)
#include <stdexcept>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/path.h"

#include "testutils/refdata.h"
#include "testutils/stringtest.h"
#include "testutils/testfilemanager.h"
#include "testutils/testoptions.h"

namespace gmx
{

namespace data_struct_test
{


//! flag that enables (true) or disables (false) debug output
//! can be set to via commandline option flag -[no]debug
bool g_debugBool = false;
//! \cond
GMX_TEST_OPTIONS(MyTestOptions, options)
{
    options->addOption(gmx::BooleanOption("debug").store(&g_debugBool)
                           .description("Print additional debug output if available."));
}
//! \endcond


} // namespace data_struct_test

} // namespace gmx
