/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * Declares gmx::LocalAtomSetManager
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_LOCALATOMSETMANAGER_H
#define GMX_DOMDEC_LOCALATOMSETMANAGER_H

#include <memory>
#include <type_traits>

#include "gromacs/utility/basedefinitions.h"

class gmx_ga2la_t;

namespace gmx
{
template<typename>
class ArrayRef;
class LocalAtomSet;

/*! \libinternal \brief
 * Hands out handles to local atom set indices and triggers index recalculation
 * for all sets upon domain decomposition if run in parallel.
 *
 * \inlibraryapi
 * \ingroup module_domdec
 */
class LocalAtomSetManager
{
public:
    LocalAtomSetManager();
    ~LocalAtomSetManager();
#ifndef DOXYGEN
    /*! \brief Add a new atom set to be managed and give back a handle.
     *
     * \todo remove this routine once all indices are represented as
     *       gmx::Index instead of int.
     *
     * \note Not created if the internal int type does match index
     *
     * \tparam T template parameter to use SFINAE for conditional function
     *           activation
     * \tparam U template parameter for conditional function activation
     *
     * \param[in] globalAtomIndex Indices of the atoms to be managed
     * \returns Handle to LocalAtomSet.
     */
    template<typename T = void, typename U = std::enable_if_t<!std::is_same_v<int, Index>, T>>
    LocalAtomSet add(ArrayRef<const int> globalAtomIndex);
#endif
    /*! \brief Add a new atom set to be managed and give back a handle.
     *
     * \param[in] globalAtomIndex Indices of the atoms to be managed
     * \returns Handle to LocalAtomSet.
     */
    LocalAtomSet add(ArrayRef<const Index> globalAtomIndex);

    /*! \brief Recalculate local and collective indices from ga2la.
     * Uses global atom to local atom lookup structure to
     * update atom indices.
     */
    void setIndicesInDomainDecomposition(const gmx_ga2la_t& ga2la);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
