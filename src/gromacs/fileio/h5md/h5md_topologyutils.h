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

#include <unordered_map>
#include <vector>

#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

namespace detail
{

/*! \brief Create a topology from a molecule type.
 *
 * This function is designed for obtaining the AtomRange of the desired molecule type
 * for writing atomic properties and residue information.
 *
 * \param[out] mtop The pointer to the topology to create.
 * \param[in] sourceMolType The source molecule type to create the topology from.
 */
void mtopFromMolType(gmx_mtop_t* mtop, const gmx_moltype_t& sourceMolType);

} // namespace detail

//! \brief Map the global indices of the selection to the local indices within the selection
using IndexMap = std::unordered_map<int32_t, int32_t>;

//! \brief List of bond representations as pairs of atom indices
using BondPairs = std::vector<std::pair<int64_t, int64_t>>;

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

/*! \brief Write the covalent connectivity to the H5MD file.
 *
 * The result bond connectivity is stored in the group `<baseContainer>/bonds`.
 * The connectivity information is deduced from the topology and written as an [numBonds x 2] dataset of int64_t.
 * If a selection is applied, only the bonds where both bonded particles are within the selection are considered.
 *
 * \param[in] topology The topology of the simulation system.
 * \param[in] baseContainer The HDF5 container to write the connectivity to.
 * \param[in] selectedAtomsIndexMap The index map of the selected atoms, not provided means considering all atoms. Providing an empty map throws an internal error.
 */
void writeBonds(const gmx_mtop_t&              topology,
                const hid_t                    baseContainer,
                const std::optional<IndexMap>& selectedAtomsIndexMap = std::nullopt);

/*! \brief Write the additional annotation of disulfide bonds to the H5MD file.
 *
 * The result bond connectivity is stored in the group `<baseContainer>/disulfide_bonds`.
 *
 * \param[in] topology The topology of the simulation system.
 * \param[in] baseContainer The HDF5 container to write the connectivity to.
 * \param[in] selectedAtomsIndexMap The index map of the selected atoms, not provided means considering all atoms. Providing an empty map throws an internal error.
 */
void writeDisulfideBonds(const gmx_mtop_t&              topology,
                         const hid_t                    baseContainer,
                         const std::optional<IndexMap>& selectedAtomsIndexMap = std::nullopt);

/*! \brief Label the version of the internal topology module.
 *
 * \param[in] baseContainer The HDF5 container to write to.
 */
void labelInternalTopologyVersion(const hid_t baseContainer);

/*! \brief Label the name of the simulation system.
 *
 * \param[in] baseContainer The HDF5 container to write to.
 * \param[in] topName The name of the topology.
 */
void labelTopologyName(const hid_t baseContainer, const char* topName);

/*! \brief Write a set of molecule type blocks in \p moltypes to the HDF5 container \p baseContainer.
 *
 * The hierarchy is as follows:
 *
 * /h5md/modules/gromacs_topology  (group (baseContainer) for H5md internal topology module)
 * \++ molecule_names              (attribute for the names of molecules in the system)
 * \-- molecule1                   (group for molecule 1, same group name as molecule_names[0])
 *     \++ nr_particles            (attribute for the number of particles in molecule 1)
 *     \++ nr_residues             (attribute for the number of residues in molecule 1)
 *     \-- id                      (dataset for the atomic identifier in molecule 1)
 *     \-- mass                    (dataset for the atomic masses in molecule 1)
 *     \-- charge                  (dataset for the atomic charges in molecule 1)
 *     \-- species                 (dataset for the atomic species in molecule 1)
 *     \-- particle_name           (dataset for the atomic names (indices into particle name table) in molecule 1)
 *     \-- particle_name_table     (dataset for the atomic name lookup table in molecule 1)
 *     \-- residue_id              (dataset for the residue identifier in molecule 1)
 *     \-- sequence                (dataset for the residue sequence (indices into residue name table) in molecule 1)
 *     \-- residue_name            (dataset for the residue names (indices into residue name table) in molecule 1)
 *     \-- residue_name_table      (dataset for the residue name lookup table in molecule 1)
 *
 * \note: The molecule_names attribute is needed when writing the molecule blocks. Hence, firstly run
 *        writeMoleculeTypes() and then run writeMoleculeBlocks().
 *
 * \param[in] baseContainer The HDF5 container to write to
 * \param[in] moltypes The molecule types to write
 */
void writeMoleculeTypes(const hid_t baseContainer, const ArrayRef<const gmx_moltype_t> moltypes);

/*! \brief Write the molecule block information to the HDF5 container of GROMACS internal topology.
 *
 * \note Write molecule type information by calling writeMoleculeTypes(), before writing molecule blocks.
 *       Without the list of molecule type names, the function will throw a file IO error.
 *       The result molecule block information will be put into the corresponding
 *       `<baseContainer>/<moleculeTypeName>` group.
 *
 * \param[in] baseContainer The HDF5 container to write the topology information to.
 * \param[in] molBlocks The molecule blocks to write.
 */
void writeMoleculeBlocks(const hid_t baseContainer, const ArrayRef<const gmx_molblock_t> molBlocks);


} // namespace gmx

#endif // GMX_FILEIO_H5MD_TOPOLOGYUTILS_H
