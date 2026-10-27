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

#include <filesystem>
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
     * \param[in] baseName The basename of the input files and the generated TPR.
     * \param[in] mdpContent Optionally, content of the MDP file used to generate the tpr
     */
    TprAndFileManager(const std::string& baseName, const std::string& mdpContent = "");

    //! Set the path to the input configuration file (-c for grompp)
    void setInputConf(const std::filesystem::path& name);

    //! Set the path to the position restraints file (-r for grompp)
    void setPositionRestraintConf(const std::filesystem::path& name);

    //! Set the path to the topology file (-p for grompp)
    void setTopology(const std::filesystem::path& name);

    /*! \brief Construct a .tpr file and return a path to it.
     *
     * \note Subsequent calls do not generate a new .tpr file, but only
     * return the same handle.
     */
    const std::filesystem::path& tprName() const;

private:
    //! Filemanager, needed to clean up the file later.
    TestFileManager fileManager_;

    /*! \brief Whether grompp has been called and a .tpr file generated
     *
     * Mutable so the lazy construction of the tpr
     * file does not inhibit use cases that are
     * const correct from the caller's perspective.
     */
    mutable bool hasGeneratedTpr_ = false;

    //! Output .tpr file (-o)
    const std::filesystem::path tprFileName_;

    //! Input .mdp file (-f)
    const std::filesystem::path mdpInputFileName_;

    //! Input configuration (-c)
    std::filesystem::path inputConf_;

    //! Position restraints file (-r), only used if not-empty
    std::filesystem::path positionRestraints_;

    //! Topology file (-p)
    std::filesystem::path topology_;
};

} // namespace test
} // namespace gmx

#endif
