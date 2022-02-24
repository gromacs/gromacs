/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \libinternal \file
 * \brief
 * Helper for generating reusuable TPR files for tests within the same test binary.
 *
 * \ingroup module_testutils
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */
#ifndef GMX_TESTUTILS_TPRFILEGENERATOR_H
#define GMX_TESTUTILS_TPRFILEGENERATOR_H

#include <memory>
#include <string>

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

class TestFileManager;

/*! \libinternal \brief
 * Helper to bundle generated TPR and the file manager to clean it up.
 */
class TprAndFileManager
{
public:
    /*! \brief
     * Generates the file when needed.
     *
     * \param[in] name The basename of the input files and the generated TPR.
     */
    TprAndFileManager(const std::string& name);
    //! Access to the string.
    const std::string& tprName() const { return tprFileName_; }

private:
    //! Tpr file name.
    std::string tprFileName_;
    //! Filemanager, needed to clean up the file later.
    TestFileManager fileManager_;
};

} // namespace test
} // namespace gmx

#endif
