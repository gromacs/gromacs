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

/*! \brief Declares the topology utility functions for H5MD.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#ifndef GMX_FILEIO_H5MD_TOPOLOGYUTILS_H
#define GMX_FILEIO_H5MD_TOPOLOGYUTILS_H

#include <hdf5.h>

#include <string>
#include <unordered_map>
#include <vector>

#include "gromacs/topology/index.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

//! \brief Map the global indices of the selection to the local indices within the selection
using IndexMap = std::unordered_map<int32_t, int32_t>;

/*! \brief Map the selected atoms to local internal indices
 *
 * \param[in] selectedIndices An array of selected atom indices.
 * \return A index map with the original index as key and internal index as value.
 */
IndexMap mapSelectionToInternalIndices(const ArrayRef<const int32_t>& selectedIndices);

/*! \brief Write atomic properties to the H5MD file.
 *
 * The result atomic properties are directly stored in the \p baseContainer.
 * Atomic information includes:
 * - attribute:nr_particles: The number of particles in the group
 * - dataset:particle_id: The unique identifier for the particle
 * - dataset:particle_name: The name of the particle, stored as an index to a lookup table
 * - dataset:particle_name_table (string[numUniqueNames]): The lookup table for particle names
 * - dataset:particle_species: The atomic number of the particle
 * - dataset:charge: The partial charge of the particle
 * - dataset:mass: The mass of the particle
 *
 * \param[in] atomRange Atom range for iterating all atoms and accessing atomic properties.
 * \param[in] baseContainer The HDF5 container to write the atomic properties to.
 * \param[in] selectedAtomsIndexMap The index map of the selected atoms, not provided means considering all atoms. Providing an empty map throws an internal error.
 */
void writeAtomicProperties(AtomRange&                     atomRange,
                           const hid_t                    baseContainer,
                           const std::optional<IndexMap>& selectedAtomsIndexMap = std::nullopt);

/*! \brief Write residue and sequence information to the H5MD file
 *
 * The result residue-related information is stored in the group `baseContainer/<groupName>`.
 * Including:
 * - residue_id: The unique identifier for the residue
 * - sequence: The sequence of residues, stored as indices to the lookup table `residue_name_table`
 * - residue_name: The name of the residue, stored as indices to the lookup table `residue_name_table`
 * - residue_name_table (string[numUniqueResidueNames]): The lookup table for residue names
 *
 * \param[in] atomRange Atom range for iterating all atoms and accessing atomic properties.
 * \param[in] baseContainer The HDF5 container to write the residue information to.
 * \param[in] selectedAtomsIndexMap The index map of the selected atoms, not provided means considering all atoms. Providing an empty map throws an internal error.
 */
void writeResidueInfo(AtomRange&                     atomRange,
                      const hid_t                    baseContainer,
                      const std::optional<IndexMap>& selectedAtomsIndexMap = std::nullopt);

} // namespace gmx

#endif // GMX_FILEIO_H5MD_TOPOLOGYUTILS_H
