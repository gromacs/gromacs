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

/*! \libinternal \file
 * \brief Declarations of utility functions to read and write index groups in H5MD files.
 *
 * Each index group is an HDF5 group named after it, holding at least one of two data sets,
 * whose union is the index group:
 *  - `indices`: [N] plain particle indices.
 *  - `ranges`:  [N, 3] atom ranges, each row [start, stop, step]. The interval is half-open
 *    (Python-style): start is inclusive, stop is exclusive.
 *
 * Examples of atom ranges:
 *  - IndexGroup [0, 1, 2, 3]       -> atom range [0, 4, 1]
 *  - IndexGroup [2, 4, 6, 8, 10]   -> atom range [2, 12, 2]
 *
 * Indices are stored as 0-based integers, matching the in-memory \c IndexGroup representation.
 * They may be in any order and may repeat, but a group is stored as a set of particles: it
 * reads back sorted and without repeats.
 *
 * A group with no particles is stored as an `indices` data set of length zero. Holding neither
 * data set is reserved for selections that are not index-based, so such a group is not an
 * index group.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#ifndef GMX_FILEIO_H5MD_INDEXGROUP_H
#define GMX_FILEIO_H5MD_INDEXGROUP_H

#include <hdf5.h>

#include <string>
#include <vector>

struct IndexGroup;

namespace gmx
{

template<typename>
class ArrayRef;

/*! \brief Get the names of all index groups in the base HDF5 container.
 *
 * \param[in] baseContainer The base HDF5 container for group information.
 * \returns A vector of index group names.
 *
 * \throws H5mdError If the links of \p baseContainer cannot be iterated.
 */
std::vector<std::string> getNamesOfIndexGroups(hid_t baseContainer);

/*! \brief Write an array of index groups to the base HDF5 container.
 *
 * The groups are stored in the format described in the file documentation.
 *
 * \param[in] baseContainer The base HDF5 container for group information.
 * \param[in] indexGroups   The vector of \c IndexGroup objects to write.
 *
 * \throws H5mdError If \p indexGroups is empty, if any group name is empty or contains the
 *                   HDF5 path separator, if a group of that name already exists in
 *                   \p baseContainer or is repeated within \p indexGroups, if any index of
 *                   a group is negative, or if a group or one of its datasets cannot be
 *                   created.
 */
void writeIndexGroups(hid_t baseContainer, ArrayRef<const IndexGroup> indexGroups);

/*! \brief Write a single index group to the base HDF5 container.
 *
 * \note
 * This function is a wrapper of \ref writeIndexGroups(), provided for convenience when
 * only a single group is available. Writing multiple groups in one call is preferred,
 * since the sanity check then sees the whole batch at once.
 *
 * \param[in] baseContainer The base HDF5 container for group information.
 * \param[in] indexGroup    The \c IndexGroup object to write.
 *
 * \throws H5mdError If \p indexGroup has an empty name, a name containing the HDF5 path
 *                   separator, a name that already exists in \p baseContainer, a negative
 *                   index, or if the group or one of its datasets cannot be created.
 */
void writeIndexGroup(hid_t baseContainer, const IndexGroup& indexGroup);

/*! \brief Read a single index group by its name.
 *
 * The group is read from the format described in the file documentation.
 *
 * \param[in]  baseContainer The base HDF5 container for group information.
 * \param[in]  groupName     The name of the index group to read.
 * \returns An \c IndexGroup object of the target group.
 *
 * \throws H5mdError If \p groupName does not exist in \p baseContainer or cannot be opened,
 *                   if it holds neither an `indices` nor a `ranges` dataset, or if its
 *                   `ranges` dataset is not [N, 3] or holds a row with a non-positive step
 *                   or a negative start.
 */
IndexGroup readIndexGroup(hid_t baseContainer, const char* groupName);

/*! \brief Read all registered index groups from the base container.
 *
 * \note
 * This function returns a `std::vector<IndexGroup>` to match the `init_index` function
 * when reading index groups from an .ndx file. The group names are enumerated from the
 * links of \p baseContainer via \ref getNamesOfIndexGroups(), so the groups are returned in
 * alphabetical order rather than in the order they were written. This function then wraps
 * \ref readIndexGroup() to read each group.
 *
 * \param[in] baseContainer The base HDF5 container for group information.
 * \returns A vector of \c IndexGroup objects.
 *
 * \throws H5mdError If the links of \p baseContainer cannot be listed, or if any of the
 *                   index groups it holds cannot be read.
 */
std::vector<IndexGroup> readIndexGroups(hid_t baseContainer);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_INDEXGROUP_H
