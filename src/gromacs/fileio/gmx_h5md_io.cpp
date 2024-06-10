/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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

#include "gmxpre.h"

#include "h5md_io.h"

#include "config.h"

#include <strings.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <sys/_types/_int64_t.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

#include "gmx_h5md_io.h"
#include "h5md_datablock.h"
#include "h5md_io.h"
#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "external/SZ3-bio/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif

namespace
{

void setupSystemParticleProperties(gmx::h5mdio::GmxH5mdIo*  file,
                                   const t_atoms&           atoms,
                                   gmx::ArrayRef<const int> selectionIndices,
                                   std::string              selectionName)
{
    /* Vectors are used to keep the values in a continuous memory block. */
    std::vector<real> atomCharges;
    std::vector<real> atomMasses;
    std::vector<int>  atomSpecies;
    /* Since the system block contains all atoms it is not necessary to record the ID,
     * but we do that in order to allow changing the mapping or "remove" particles,
     * in order to enable grand canonical simulations. */
    std::vector<int> atomIds;

    const size_t numSelectedParticles = selectionIndices.size() > 0 ? selectionIndices.size() : atoms.nr;

    atomCharges.reserve(numSelectedParticles);
    atomMasses.reserve(numSelectedParticles);
    atomSpecies.reserve(numSelectedParticles);
    atomIds.reserve(numSelectedParticles);

    /* FIXME: Should use int64_t. Needs changes in atoms. */
    for (size_t i = 0; i < numSelectedParticles; i++)
    {
        size_t iParticle = selectionIndices.size() > 0 ? selectionIndices[i] : i;
        atomCharges.push_back(atoms.atom[iParticle].q);
        atomMasses.push_back(atoms.atom[iParticle].m);
        atomSpecies.push_back(atoms.atom[iParticle].type);
        atomIds.push_back(iParticle);
    }

    file->setNumericProperty("/particles/" + selectionName, "charge", atomCharges, "", false);
    file->setNumericProperty("/particles/" + selectionName, "mass", atomMasses, "amu", false);
    file->setNumericProperty("/particles/" + selectionName, "species", atomSpecies, "", false);
    file->setNumericProperty("/particles/" + selectionName, "id", atomIds, "", false);
}

/*! \brief Add atom type entries (species) for all different atom types in \p atoms.
 *
 * \param[in]     file           The H5MD file manager to use.
 * \param[in]     atoms          The GROMACS atoms to iterate through to add their corresponding
 * atom types (species)
 * \param[in,out] atomTypesAdded Keeps track of which atom types have been added already.
 * \throws FileIOError If there was en error adding the atom type entries.
 */
void addAtomTypesOfAtoms(gmx::h5mdio::GmxH5mdIo* file, const t_atoms& atoms, std::vector<bool>& atomTypesAdded)
{
    hid_t atomTypesGroup =
            file->createGroup(gmx::h5mdio::s_gromacsTopologyGroupName + "/atom_species");
    hid_t   dataType     = H5Tcopy(H5T_NATIVE_INT);
    hsize_t chunkDims[1] = { atomTypesAdded.size() };
    hid_t   atomTypeAtomicNumberDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(atomTypesGroup,
                                                "atomic_number",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::LosslessNoShuffle,
                                                0);
    for (int i = 0; i < atoms.nr; i++)
    {
        t_atom* atom = &atoms.atom[i];
        if (!atomTypesAdded[atom->type])
        {
            gmx::h5mdio::writeData<1, false>(atomTypeAtomicNumberDataSet, &atom->atomnumber, atom->type);
            atomTypesAdded[atom->type] = true;
        }
    }
}


/*! \brief Get the number of atoms of the molecule type specified by \p molTypeName.
 * \param[in] file               The H5MD file manager to use.
 * \param[in] molTypeName        The name of the molecule type.
 * \returns the number of atoms in the molecule type or -1 if the molecule type could not be found.
 * \throws FileIOError If there was an error reading the molecule type information.
 */
int64_t getNumberOfAtomsOfMoleculeTypeByName(gmx::h5mdio::GmxH5mdIo* file, std::string molTypeName)
{
#if GMX_USE_HDF5
    std::string moleculeTypesGroupName = gmx::h5mdio::s_gromacsTopologyGroupName + "/molecule_types";
    std::string moleculeTypeName       = moleculeTypesGroupName + "/" + molTypeName;
    hid_t       molelculeTypeGroup     = file->getGroupId(moleculeTypeName);

    if (molelculeTypeGroup < 0)
    {
        return -1;
    }
    std::int64_t numAtoms;
    gmx::h5mdio::getAttribute(molelculeTypeGroup, "number_of_atoms", &numAtoms);

    return numAtoms;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

/*! \brief Add a block consisting of a number of copies of a molecule type to the GROMACS topology section in the file.
 * \param[in] file          The H5MD file manager to use.
 * \param[in] moleculeTypeName The name of the molecule type of this molecule block.
 * \param[in] molBlockIndex The index of the molecule block.
 * \param[in] numMol The number of molecules of this type (in this molecule block).
 * \param[in] molBlockIndices The GROMACS data structure containing the indices of the molecule block.
 * \throws FileIOError If there was an error adding the molecule type information.
 */
void addBlockOfMoleculeType(gmx::h5mdio::GmxH5mdIo*     file,
                            const std::string&          moleculeTypeName,
                            size_t                      molBlockIndex,
                            size_t                      numMol,
                            const MoleculeBlockIndices& molBlockIndices)
{
#if GMX_USE_HDF5
    std::string moleculeBlocksName  = gmx::h5mdio::s_gromacsTopologyGroupName + "/molecule_blocks";
    hid_t       moleculeBlocksGroup = file->createGroup(moleculeBlocksName);

    hid_t stringDataType = H5Tcopy(H5T_C_S1);
    H5Tset_cset(stringDataType, H5T_CSET_UTF8);
    size_t maxNameStringLength = gmx::h5mdio::c_moleculeTypeStringLen;
    H5Tset_size(stringDataType, maxNameStringLength);
    hid_t   dataType     = H5Tcopy(H5T_NATIVE_INT64);
    hsize_t chunkDims[1] = { 1 };

    hid_t moleculeTypeNameDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "molecule_type",
                                                nullptr,
                                                stringDataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    gmx::h5mdio::writeData<1, false>(moleculeTypeNameDataSet, moleculeTypeName.c_str(), molBlockIndex);

    hid_t        numMolDataSet = gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                              "number_of_molecules",
                                                              nullptr,
                                                              dataType,
                                                              chunkDims,
                                                              gmx::h5mdio::CompressionAlgorithm::None,
                                                              0);
    std::int64_t tmpValue      = numMol;
    gmx::h5mdio::writeData<1, false>(numMolDataSet, &tmpValue, molBlockIndex);

    hid_t numAtomsPerMoleculeDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "num_atoms_per_molecule",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.numAtomsPerMolecule;
    gmx::h5mdio::writeData<1, false>(numAtomsPerMoleculeDataSet, &tmpValue, molBlockIndex);

    hid_t globalAtomsStartDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "global_atom_start",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.globalAtomStart;
    gmx::h5mdio::writeData<1, false>(globalAtomsStartDataSet, &tmpValue, molBlockIndex);

    hid_t globalAtomsEndDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "global_atom_end",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.globalAtomEnd;
    gmx::h5mdio::writeData<1, false>(globalAtomsEndDataSet, &tmpValue, molBlockIndex);

    hid_t globalResidueStartDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "global_residue_start",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.globalResidueStart;
    gmx::h5mdio::writeData<1, false>(globalResidueStartDataSet, &tmpValue, molBlockIndex);

    hid_t residueNumberStartDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "residue_number_start",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.residueNumberStart;
    gmx::h5mdio::writeData<1, false>(residueNumberStartDataSet, &tmpValue, molBlockIndex);

    hid_t moleculeIndexStartDataSet =
            gmx::h5mdio::openOrCreateDataSet<1>(moleculeBlocksGroup,
                                                "molecule_index_start",
                                                nullptr,
                                                dataType,
                                                chunkDims,
                                                gmx::h5mdio::CompressionAlgorithm::None,
                                                0);
    tmpValue = molBlockIndices.moleculeIndexStart;
    gmx::h5mdio::writeData<1, false>(moleculeIndexStartDataSet, &tmpValue, molBlockIndex);


#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}
/*! \brief Add a molecule type to the GROMACS topology section in the file.
 * \param[in] file The H5MD file manager to use.
 * \param[in] molType The molecule type to add.
 * \returns the H5MD ID of the molecule type group
 * \throws FileIOError If there was an error adding the molecule type information.
 */
hid_t addMoleculeType(gmx::h5mdio::GmxH5mdIo* file, const gmx_moltype_t& molType)
{
#if GMX_USE_HDF5
    std::string moleculeTypesGroupName = gmx::h5mdio::s_gromacsTopologyGroupName + "/molecule_types";
    file->createGroup(moleculeTypesGroupName);
    std::string moleculeTypeName  = moleculeTypesGroupName + "/" + (*molType.name);
    hid_t       moleculeTypeGroup = file->createGroup(moleculeTypeName);

    gmx::h5mdio::setAttribute(
            moleculeTypeGroup, "number_of_atoms", static_cast<int64_t>(molType.atoms.nr), H5T_NATIVE_INT64);

    std::vector<std::string>  atomNames(molType.atoms.nr);
    std::vector<int>          atomTypes(molType.atoms.nr);
    std::vector<int>          atomTypesB(molType.atoms.nr);
    std::vector<std::string>  residueNames(molType.atoms.nr);
    std::vector<std::int64_t> residueNumbers(molType.atoms.nr);
    std::vector<std::string>  chainIds(molType.atoms.nr);

    for (ssize_t i = 0; i < molType.atoms.nr; i++)
    {
        int residueIndex = molType.atoms.atom[i].resind;

        atomNames[i]      = *molType.atoms.atomname[i];
        atomTypes[i]      = molType.atoms.atom[i].type;
        atomTypesB[i]     = molType.atoms.atom[i].typeB;
        residueNames[i]   = *molType.atoms.resinfo[residueIndex].name;
        residueNumbers[i] = molType.atoms.resinfo[residueIndex].nr;
        chainIds[i]       = molType.atoms.resinfo[residueIndex].chainid;
    }

    file->setStringProperty(
            moleculeTypeName, "atom_name", atomNames, false, gmx::h5mdio::c_atomResidueStringLen);
    file->setNumericProperty(moleculeTypeName, "atom_species", atomTypes, "", false);
    file->setNumericProperty(moleculeTypeName, "atom_species_state_b", atomTypesB, "", false);
    file->setStringProperty(
            moleculeTypeName, "residue_name", residueNames, false, gmx::h5mdio::c_atomResidueStringLen);
    file->setNumericProperty(moleculeTypeName, "residue_number", residueNumbers, "", false);
    file->setStringProperty(moleculeTypeName, "chain_id", chainIds, false, 1);

    return moleculeTypeGroup;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

/*! \brief Adds chemical bonds (including constraints and settle) in a "connectivity" dataset
 * of the group \p molTypeGroup
 *
 * \param[in] file         The H5MD file manager to use.
 * \param[in] molTypeGroup The H5MD group of the molecule type.
 * \param[in] molType      The molecule type, of which to add bonds.
 * \param[in] numMols      Number of molecules of this type (in this molecule block).
 * \param[in] index        Indices of selected atoms.
 * \param[in] selectionName The name of the selection.
 * \param[in] allSystemBonds Pointer to a vector that will contain a full connectivity table of the system.
 * Can be a nullptr, in which case the bonds will not be added.
 * \param[in] selectionBonds Pointer to a vector that will contain a full connectivity table for a selection
 * of the system. Can be a nullptr, in which case the bonds will not be added.
 * \throws FileIOError If there was an error adding the connectivity records.
 */
void addMoleculeTypeBondsToTopology(gmx::h5mdio::GmxH5mdIo* gmx_unused        file,
                                    hid_t                                     molTypeGroup,
                                    const gmx_moltype_t&                      molType,
                                    int64_t                                   numMols,
                                    gmx::ArrayRef<const int>                  index,
                                    std::vector<std::pair<int64_t, int64_t>>* systemBonds,
                                    std::vector<std::pair<int64_t, int64_t>>* selectionBonds)
{
    std::vector<std::pair<int64_t, int64_t>> bonds;
    /* Bonds have to be deduced from interactions (constraints etc). Different
     * interactions have different sets of parameters. */
    for (int i = 0; i < F_NRE; i++)
    {
        if (IS_CHEMBOND(i))
        {
            const InteractionList& ilist         = molType.ilist[i];
            int                    fromAtomIndex = 1;
            while (fromAtomIndex < ilist.size())
            {
                // int64_t atoms[2] = {ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]};
                // gmx::h5mdio::writeData<2, false>(moleculeConnectivityDataSet, atoms, bondCount++);
                bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]);
                fromAtomIndex += 3;
            }
        }
    }
    /* Settle is described using three atoms */
    const InteractionList& ilist         = molType.ilist[F_SETTLE];
    int                    fromAtomIndex = 1;
    while (fromAtomIndex < ilist.size())
    {
        // int64_t atoms[2] = {ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]};
        // gmx::h5mdio::writeData<2, false>(moleculeConnectivityDataSet, atoms, bondCount++);
        // atoms[1] = ilist.iatoms[fromAtomIndex + 2];
        // gmx::h5mdio::writeData<2, false>(moleculeConnectivityDataSet, atoms, bondCount++);
        bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]);
        bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 2]);
        fromAtomIndex += 4;
    }

    if (systemBonds != nullptr || selectionBonds != nullptr)
    {
        for (int64_t molIterator = 0; molIterator < numMols; molIterator++)
        {
            for (size_t bondIterator = 0; bondIterator < bonds.size(); bondIterator++)
            {
                int64_t offset = molIterator * molType.atoms.nr;
                if (systemBonds != nullptr)
                {
                    systemBonds->emplace_back(bonds[bondIterator].first + offset,
                                              bonds[bondIterator].second + offset);
                }
                if (selectionBonds != nullptr)
                {
                    if (std::find(index.begin(), index.end(), bonds[bondIterator].first + offset)
                                != index.end()
                        && std::find(index.begin(), index.end(), bonds[bondIterator].second + offset)
                                   != index.end())
                    {
                        selectionBonds->emplace_back(bonds[bondIterator].first + offset,
                                                     bonds[bondIterator].second + offset);
                    }
                }
            }
        }
    }

    char molTypeGroupPath[gmx::h5mdio::c_maxFullNameLength];
    H5Iget_name(molTypeGroup, molTypeGroupPath, gmx::h5mdio::c_maxFullNameLength - 1);

    file->setNumericProperty(molTypeGroupPath, "connectivity", bonds, "", false);
}

/*! \brief Check whether there is a separate selection for output.
 * TODO: Investigate more flexible ways of specifying selection groups.
 *
 * \param[in] index    The selected atoms to include. If empty, check whether there is a separate
 * compression group.
 */
bool hasSeparateSelection(/*const gmx_mtop_t& topology,*/ gmx::ArrayRef<const int> index)
{
    /* We only need to create a separate selection group entry if not all atoms are part of it. */
    /* If a selection of atoms is explicitly provided then use that instead of the CompressedPositionOutput */
    bool separateSelection = false;
    if (index.ssize() > 0)
    {
        separateSelection = true;
    }
    /* TODO: Check if the below code is necessary.
    // else
    // {
    //     /* FIXME: Should use int64_t. Needs changes in topology. */
    //     for (int i = 0; i < topology.natoms; i++)
    //     {
    //         if (getGroupType(topology.groups, SimulationAtomGroupType::CompressedPositionOutput, i) != 0)
    //         {
    //             separateSelection = true;
    //             break;
    //         }
    //     }
    // }
    return separateSelection;
}


} // namespace

namespace gmx
{

void setH5mdAuthorAndCreator(h5mdio::GmxH5mdIo* file)
{
    char tmpUserName[gmx::h5mdio::c_maxFullNameLength];
    if (!gmx_getusername(tmpUserName, gmx::h5mdio::c_maxFullNameLength))
    {
        file->setAuthor(tmpUserName);
    }

    std::string precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif
    std::string programInfo = gmx::getProgramContext().displayName() + precisionString;
    file->setCreatorProgramName(programInfo);

    const std::string gmxVersion = gmx_version();
    file->setCreatorProgramVersion(gmxVersion);
}

void setupMolecularSystemParticleData(h5mdio::GmxH5mdIo*       file,
                                      const gmx_mtop_t&        topology,
                                      gmx::ArrayRef<const int> index,
                                      std::string              selectionName)
{
#if GMX_USE_HDF5
    t_atoms atoms = gmx_mtop_global_atoms(topology);

    if (atoms.nr == 0)
    {
        return;
    }

    setupSystemParticleProperties(file, atoms, gmx::ArrayRef<const int>(), "system");

    /* We only need to create a separate selection group entry if not all atoms are part of it. */
    /* If a selection of atoms is explicitly provided then use that instead of the CompressedPositionOutput */
    bool separateSelection = hasSeparateSelection(index);
    if (separateSelection)
    {
        std::string systemOutputName;
        if (index.ssize() > 0 && selectionName != "")
        {
            systemOutputName = selectionName;
        }
        /* If no name was specified fall back to using the selection group name of compressed output, if any. */
        else if (topology.groups.numberOfGroupNumbers(SimulationAtomGroupType::CompressedPositionOutput) != 0)
        {
            int nameIndex = topology.groups.groups[SimulationAtomGroupType::CompressedPositionOutput][0];
            systemOutputName = *topology.groups.groupNames[nameIndex];
        }
        setupSystemParticleProperties(file, atoms, index, systemOutputName);
    }

    done_atom(&atoms);
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

MoleculeBlockIndices getMoleculeBlockIndicesByIndex(h5mdio::GmxH5mdIo* file, size_t molBlockIndex)
{
#if GMX_USE_HDF5
    std::string moleculeBlocksName  = h5mdio::s_gromacsTopologyGroupName + "/molecule_blocks";
    hid_t       moleculeBlocksGroup = file->getGroupId(moleculeBlocksName);

    MoleculeBlockIndices molBlockIndices;
    if (moleculeBlocksGroup < 0)
    {
        return molBlockIndices;
    }

    std::vector<std::string> moleculeTypeNames =
            file->readStringProperty(moleculeBlocksName, "molecule_type");
    std::string moleculeTypeName = moleculeTypeNames[molBlockIndex];

    molBlockIndices.numAtomsPerMolecule = getNumberOfAtomsOfMoleculeTypeByName(file, moleculeTypeName);

    void* buffer                     = nullptr;
    hid_t numAtomsPerMoleculeDataSet = H5Dopen(moleculeBlocksGroup, "num_atoms_per_molecule", H5P_DEFAULT);
    h5mdio::readData<1>(numAtomsPerMoleculeDataSet, molBlockIndex, &buffer);
    molBlockIndices.numAtomsPerMolecule = *(static_cast<int64_t*>(buffer));

    hid_t globalAtomStartDataSet = H5Dopen(moleculeBlocksGroup, "global_atom_start", H5P_DEFAULT);
    h5mdio::readData<1>(globalAtomStartDataSet, molBlockIndex, &buffer);
    molBlockIndices.globalAtomStart = *(static_cast<int64_t*>(buffer));

    hid_t globalAtomEndDataSet = H5Dopen(moleculeBlocksGroup, "global_atom_end", H5P_DEFAULT);
    h5mdio::readData<1>(globalAtomEndDataSet, molBlockIndex, &buffer);
    molBlockIndices.globalAtomEnd = *(static_cast<int64_t*>(buffer));

    hid_t globalResidueStartDataSet = H5Dopen(moleculeBlocksGroup, "global_residue_start", H5P_DEFAULT);
    h5mdio::readData<1>(globalResidueStartDataSet, molBlockIndex, &buffer);
    molBlockIndices.globalResidueStart = *(static_cast<int64_t*>(buffer));

    hid_t residueNumberStartDataSet = H5Dopen(moleculeBlocksGroup, "residue_number_start", H5P_DEFAULT);
    h5mdio::readData<1>(residueNumberStartDataSet, molBlockIndex, &buffer);
    molBlockIndices.residueNumberStart = *(static_cast<int64_t*>(buffer));

    hid_t moleculeIndexStartDataSet = H5Dopen(moleculeBlocksGroup, "molecule_index_start", H5P_DEFAULT);
    h5mdio::readData<1>(moleculeIndexStartDataSet, molBlockIndex, &buffer);
    molBlockIndices.moleculeIndexStart = *(static_cast<int64_t*>(buffer));

    free(buffer);

    return molBlockIndices;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void setupMolecularSystemTopology(h5mdio::GmxH5mdIo*       file,
                                  const gmx_mtop_t&        topology,
                                  gmx::ArrayRef<const int> index,
                                  bool                     abortIfPresent,
                                  bool                     writeVmdStructureData)
{
#if GMX_USE_HDF5
    if (file == nullptr || !file->isFileOpen())
    {
        throw gmx::FileIOError("No file open for writing.");
    }

    size_t            numMolBlocks       = topology.molblock.size();
    size_t gmx_unused numMolBlockIndices = topology.moleculeBlockIndices.size();

    GMX_ASSERT(numMolBlocks == numMolBlockIndices,
               "The number of molecule blocks and molecule block indices do not match.");

    hid_t topologyGroup = file->getGroupId(h5mdio::s_gromacsTopologyGroupName);
    if (topologyGroup >= 0 && abortIfPresent)
    {
        return;
    }

    if (topologyGroup < 0)
    {
        topologyGroup = file->createGroup(h5mdio::s_gromacsTopologyGroupName);
    }

    h5mdio::setVersionAttribute(topologyGroup,
                                h5mdio::c_gmxH5mdParametersGroupMajorVersion,
                                h5mdio::c_gmxH5mdParametersGroupMinorVersion);

    std::vector<bool>                        atomTypesAdded(topology.ffparams.atnr, false);
    std::vector<std::pair<int64_t, int64_t>> systemBonds, selectionBonds;
    for (size_t i = 0; i < numMolBlocks; i++)
    {
        const gmx_molblock_t&       molBlock        = topology.molblock[i];
        const MoleculeBlockIndices& molBlockIndices = topology.moleculeBlockIndices[i];
        const gmx_moltype_t&        molType         = topology.moltype[molBlock.type];
        const std::string           molName         = *molType.name;
        const size_t                numMol          = molBlock.nmol;
        hid_t                       molTypeGroup    = addMoleculeType(file, molType);
        if (molTypeGroup < 0)
        {
            throw gmx::FileIOError("Cannot write molecule type group.");
        }
        addMoleculeTypeBondsToTopology(
                file, molTypeGroup, molType, numMol, index, &systemBonds, &selectionBonds);
        addAtomTypesOfAtoms(file, molType.atoms, atomTypesAdded);
        addBlockOfMoleculeType(file, molName, i, numMol, molBlockIndices);
    }
    if (!systemBonds.empty())
    {
        file->setNumericProperty("/connectivity", "system", systemBonds, "", false);
    }
    if (writeVmdStructureData)
    {
        std::vector<int64_t> firstAtomsInPairs, secondAtomsInPairs;
        if (selectionBonds.empty() || selectionBonds.size() == systemBonds.size())
        {
            firstAtomsInPairs.reserve(systemBonds.size());
            secondAtomsInPairs.reserve(systemBonds.size());
            std::transform(systemBonds.begin(),
                           systemBonds.end(),
                           std::back_inserter(firstAtomsInPairs),
                           [](auto const& pair) { return pair.first; });
            std::transform(systemBonds.begin(),
                           systemBonds.end(),
                           std::back_inserter(secondAtomsInPairs),
                           [](auto const& pair) { return pair.second; });
        }
        else if (!selectionBonds.empty())
        {
            firstAtomsInPairs.reserve(selectionBonds.size());
            secondAtomsInPairs.reserve(selectionBonds.size());
            std::transform(selectionBonds.begin(),
                           selectionBonds.end(),
                           std::back_inserter(firstAtomsInPairs),
                           [](auto const& pair) { return pair.first; });
            std::transform(selectionBonds.begin(),
                           selectionBonds.end(),
                           std::back_inserter(secondAtomsInPairs),
                           [](auto const& pair) { return pair.second; });
        }
        file->setNumericProperty("/parameters/vmd_structure", "bond_from", firstAtomsInPairs, "", false);
        file->setNumericProperty("/parameters/vmd_structure", "bond_to", secondAtomsInPairs, "", false);
    }

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void writeFrameToStandardDataBlocks(h5mdio::GmxH5mdIo* file,
                                    int64_t            step,
                                    real               time,
                                    real               lambda,
                                    const rvec*        box,
                                    const int64_t      numParticles,
                                    const rvec*        x,
                                    const rvec*        v,
                                    const rvec*        f,
                                    const double       xCompressionError,
                                    const std::string  selectionName)
{
#if GMX_USE_HDF5
    if (numParticles <= 0)
    {
        throw gmx::FileIOError("There must be particles/atoms when writing trajectory frames.");
    }
    if (file == nullptr || !file->isFileOpen())
    {
        throw gmx::FileIOError("No file open for writing.");
    }

    /* There is so little lambda data per frame that it is best to write multiple per chunk. */
    hsize_t     numFramesPerChunk = 20;
    std::string wantedName        = "/observables/lambda";
    file->writeDataFrame(
            step, time, wantedName, 1, 1, &lambda, "", numFramesPerChunk, h5mdio::CompressionAlgorithm::LosslessNoShuffle);

    if (x != nullptr)
    {
        wantedName = "/particles/" + selectionName + "/position";
        h5mdio::CompressionAlgorithm compressionAlgorithm =
                h5mdio::CompressionAlgorithm::LosslessWithShuffle;
        if (xCompressionError != 0)
        {
            /* Use no more than 20 frames per chunk (compression unit). Use fewer frames per chunk if there are many atoms. */
            numFramesPerChunk    = std::min(20, int(std::ceil(5e6f / numParticles)));
            compressionAlgorithm = h5mdio::CompressionAlgorithm::LossySz3;

            /* Register the SZ3 filter. This is not necessary when creating a dataset with the filter,
             * but must be done to append to an existing file (e.g. when restarting from checkpoint). */
            h5mdio::registerSz3FilterImplicitly();
        }
        file->writeDataFrame(step,
                             time,
                             wantedName,
                             numParticles,
                             DIM,
                             static_cast<const real*>(x[0]),
                             "nm",
                             numFramesPerChunk,
                             compressionAlgorithm,
                             xCompressionError);
    }

    if (box != nullptr)
    {
        /* There is so little box data per frame that it is best to write multiple per chunk. */
        numFramesPerChunk = 20;
        wantedName        = "/particles/" + selectionName + "/box/edges";
        file->writeDataFrame(step,
                             time,
                             wantedName,
                             DIM,
                             DIM,
                             static_cast<const real*>(box[0]),
                             "nm",
                             numFramesPerChunk,
                             h5mdio::CompressionAlgorithm::LosslessNoShuffle);
    }

    /* There is no temporal compression of velocities and forces. */
    numFramesPerChunk = 1;
    if (v != nullptr)
    {
        wantedName = "/particles/" + selectionName + "/velocity";
        file->writeDataFrame(step,
                             time,
                             wantedName,
                             numParticles,
                             DIM,
                             static_cast<const real*>(v[0]),
                             "nm ps-1",
                             numFramesPerChunk,
                             h5mdio::CompressionAlgorithm::LosslessWithShuffle);
    }
    if (f != nullptr)
    {
        wantedName = "/particles/" + selectionName + "/force";
        file->writeDataFrame(step,
                             time,
                             wantedName,
                             numParticles,
                             DIM,
                             static_cast<const real*>(f[0]),
                             "kJ mol-1 nm-1",
                             numFramesPerChunk,
                             h5mdio::CompressionAlgorithm::LosslessWithShuffle);
    }
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

bool readNextFrameOfStandardDataBlocks(h5mdio::GmxH5mdIo* file,
                                       int64_t*           step,
                                       real*              time,
                                       real*              lambda,
                                       rvec*              box,
                                       rvec*              x,
                                       rvec*              v,
                                       rvec*              f,
                                       real*              xCompressionError,
                                       bool*              readLambda,
                                       bool*              readBox,
                                       bool*              readX,
                                       bool*              readV,
                                       bool*              readF,
                                       const std::string  selectionName)
{
#if GMX_USE_HDF5
    if (file == nullptr || !file->isFileOpen())
    {
        throw gmx::FileIOError("No file open for reading.");
    }

    std::string particlesNameStem = "/particles/" + selectionName;
    *readLambda = *readBox = *readX = *readV = *readF = false;

    std::tuple<int64_t, real> temporaryStepTime = file->getNextStepAndTimeToRead();
    *step                                       = std::get<0>(temporaryStepTime);
    *time                                       = std::get<1>(temporaryStepTime);

    bool didReadFrame  = false;
    *xCompressionError = -1;

    if (lambda != nullptr)
    {
        if (file->readNextFrameOfDataBlock("/observables/lambda", lambda, *step))
        {
            *readLambda  = true;
            didReadFrame = true;
        }
    }
    if (box != nullptr)
    {
        std::string boxDataName = particlesNameStem + "/box/edges";
        if (file->readNextFrameOfDataBlock(boxDataName.c_str(), static_cast<real*>(box[0]), *step))
        {
            *readBox     = true;
            didReadFrame = true;
        }
    }
    if (x != nullptr)
    {
        std::string xDataName = particlesNameStem + "/position";
        if (file->readNextFrameOfDataBlock(xDataName.c_str(), static_cast<real*>(x[0]), *step))
        {
            *readX             = true;
            didReadFrame       = true;
            *xCompressionError = file->getLossyCompressionErrorOfDataBlock(xDataName.c_str());
        }
    }
    if (v != nullptr)
    {
        std::string boxDataName = particlesNameStem + "/velocity";
        if (file->readNextFrameOfDataBlock(boxDataName.c_str(), static_cast<real*>(v[0]), *step))
        {
            *readV       = true;
            didReadFrame = true;
        }
    }
    if (f != nullptr)
    {
        std::string boxDataName = particlesNameStem + "/force";
        if (file->readNextFrameOfDataBlock(boxDataName.c_str(), static_cast<real*>(f[0]), *step))
        {
            *readF       = true;
            didReadFrame = true;
        }
    }
    return didReadFrame;
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

bool copyProvenanceRecords(h5mdio::GmxH5mdIo* srcFile, h5mdio::GmxH5mdIo* destFile)
{
    hid_t srcModulesGroup = srcFile->getGroupId("/modules");
    if (srcModulesGroup < 0)
    {
        return false;
    }
    hid_t destModulesGroup = destFile->createGroup("/modules");

    if (H5Ocopy(srcModulesGroup, "provenance", destModulesGroup, "provenance", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        return false;
    }
    return true;
}

} // namespace gmx
