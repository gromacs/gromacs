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

/*! \brief Definition of the topology utility functions for H5MD.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gromacs/fileio/h5md/h5md_topologyutils.h"

#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

// NOTE: Since H5MD is an experiment feature, we start the versioning from 0.1 and
//       increment to 1.0 when the metadata definition is stable and published.
//! \brief The major version of the GROMACS internal topology.
constexpr int c_h5mdInternalTopologyVersionMajor = 0;
//! \brief The minor version of the GROMACS internal topology.
constexpr int c_h5mdInternalTopologyVersionMinor = 1;

//! \brief Unit for mass data in GROMACS.
constexpr char c_massUnit[] = "u";
//! \brief Unit for electric charge data in GROMACS.
constexpr char c_elementaryChargeUnit[] = "e";
//! \brief Name of dataset for the identifier for the particles.
constexpr char c_particleIdName[] = "id";
//! \brief Name of dataset for the mass of the particles.
constexpr char c_particleMassName[] = "mass";
//! \brief Name of dataset for the charge of the particles.
constexpr char c_particleChargeName[] = "charge";
//! \brief Name of dataset for the species (atomic number) of the particles.
constexpr char c_particleSpeciesName[] = "species";
//! \brief Name of dataset for the names of the particles.
constexpr char c_particleNameIndicesName[] = "particle_name";
//! \brief Name of dataset for the particle name lookup table.
constexpr char c_particleNameTableName[] = "particle_name_table";
//! \brief Name of attribute for the number of the particles.
constexpr char c_numParticlesAttributeKey[] = "nr_particles";
//! \brief Name of dataset for the residue identifier of the particles.
constexpr char c_residueIdName[] = "residue_id";
//! \brief Name of dataset for the residue names of the particles.
constexpr char c_residueNameIndicesName[] = "residue_name";
//! \brief Name of dataset for the residue name lookup table.
constexpr char c_residueNameTableName[] = "residue_name_table";
//! \brief Name of dataset for the residue sequence.
constexpr char c_residueSequenceName[] = "sequence";
//! \brief Name of attribute for the number of the residues.
constexpr char c_numResiduesAttributeKey[] = "nr_residues";
//! \brief Name of bond dataset inside connectivity group
constexpr char c_bondsConnectivityName[] = "bonds";
//! \brief Name of the attribute storing number of bonds
constexpr char c_numberOfBondsAttributeKey[] = "nr_bonds";
//! \brief Name of disulfide bonds within the special bonds group
constexpr char c_disulfideBondsName[] = "disulfide_bonds";
//! \brief Name of attribute for the version of internal topology.
constexpr char c_h5mdInternalTopologyVersionAttributeKey[] = "version";
//! \brief Name of attribute for topology name.
constexpr char c_h5mdInternalTopologyNameAttributeKey[] = "system_name";
//! \brief Name of attribute for the names of the molecule types.
constexpr char c_moleculeNamesAttributeKey[] = "molecule_names";

namespace
{
//! \brief Atomic number for sulfur element
constexpr int c_sulfurAtomicNumber = 16;

// Deduced bonds from interactions list
BondPairs deduceBondsFromMolecule(const gmx_moltype_t& moltype)
{
    BondPairs bonds;
    for (InteractionFunction itype : EnumerationWrapper<InteractionFunction>{})
    {
        if (IS_CHEMBOND(itype))
        {
            const InteractionList& ilist         = moltype.ilist[itype];
            int                    fromAtomIndex = 1;
            while (fromAtomIndex < ilist.size())
            {
                bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]);
                fromAtomIndex += 3;
            }
        }
        else if (itype == InteractionFunction::SETTLE)
        {
            // Specially treat the SETTLE interaction for water molecules
            const InteractionList& ilist         = moltype.ilist[itype];
            int                    fromAtomIndex = 1;
            while (fromAtomIndex < ilist.size())
            {
                bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 1]);
                bonds.emplace_back(ilist.iatoms[fromAtomIndex], ilist.iatoms[fromAtomIndex + 2]);
                fromAtomIndex += 4;
            }
        }
    }
    return bonds;
}

BondPairs deduceDisulfideBondsFromMolecule(const gmx_moltype_t& molType)
{
    BondPairs disulfideBonds;
    for (InteractionFunction itype : EnumerationWrapper<InteractionFunction>{})
    {
        if (IS_CHEMBOND(itype))
        {
            const InteractionList& ilist         = molType.ilist[itype];
            int                    fromAtomIndex = 1;
            while (fromAtomIndex < ilist.size())
            {
                int p1 = ilist.iatoms[fromAtomIndex];
                int p2 = ilist.iatoms[fromAtomIndex + 1];
                // Check if both atoms are sulfur
                if (molType.atoms.atom[p1].atomnumber == c_sulfurAtomicNumber
                    && molType.atoms.atom[p2].atomnumber == c_sulfurAtomicNumber)
                {
                    disulfideBonds.emplace_back(p1, p2);
                }
                fromAtomIndex += 3;
            }
        }
    }
    return disulfideBonds;
}

/*! \brief Cache bonds for \p numMols of molecule blocks into \p bondsBuffer
 *
 * This is an internal utility function when iterating over molecule types
 *
 * \param[in] bonds The bond pairs to cache.
 * \param[out] bondsBuffer The buffer to store the cached bonds.
 * \param[in] numMols The number of molecule blocks.
 * \param[in] atomsPerMol The number of atoms per molecule.
 * \param[in] selectedAtomsIndexMap The index map of the selected atoms, not provided means considering all atoms. Providing an empty map throws an internal error.
 * \param[in] globalOffset The global index offset for the current molecule type.
 */
void cacheBondForNBlocks(const BondPairs&               bonds,
                         std::vector<int64_t>*          bondsBuffer,
                         const int32_t                  numMols,
                         const int32_t                  atomsPerMol,
                         const std::optional<IndexMap>& selectedAtomsIndexMap,
                         const int64_t                  globalOffset)
{
    for (int32_t n = 0; n < numMols; ++n)
    {
        for (const auto& bond : bonds)
        {
            // Compute the global index of the first molecule in this block
            auto offset = globalOffset + n * atomsPerMol;
            // Applied a selection on atom indices if any
            if (selectedAtomsIndexMap.has_value())
            {
                const auto firstIt = selectedAtomsIndexMap->find(bond.first + offset);
                if (firstIt != selectedAtomsIndexMap->end())
                {
                    const auto secondIt = selectedAtomsIndexMap->find(bond.second + offset);
                    if (secondIt != selectedAtomsIndexMap->end())
                    {
                        bondsBuffer->push_back(firstIt->second);
                        bondsBuffer->push_back(secondIt->second);
                    }
                }
            }
            else
            {
                bondsBuffer->push_back(bond.first + offset);
                bondsBuffer->push_back(bond.second + offset);
            }
        }
    }
}

} // namespace

namespace detail
{

void mtopFromMolType(gmx_mtop_t* mtop, const gmx_moltype_t& sourceMolType)
{
    mtop->moltype.resize(1);
    copy_moltype(&sourceMolType, &mtop->moltype[0]);

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;
    mtop->natoms           = sourceMolType.atoms.nr * mtop->molblock[0].nmol;
    mtop->finalize();
}

} // namespace detail

IndexMap mapSelectionToInternalIndices(const ArrayRef<const int32_t>& selectedIndices)
{
    IndexMap selectedAtomsIndexMap;
    for (size_t i = 0; i < selectedIndices.size(); ++i)
    {
        selectedAtomsIndexMap[selectedIndices[i]] = i;
    }
    return selectedAtomsIndexMap;
}

void writeAtomicProperties(AtomRange&                     atomRange,
                           const hid_t                    baseContainer,
                           const std::optional<IndexMap>& selectedAtomsIndexMap)
{
    if (selectedAtomsIndexMap.has_value() && selectedAtomsIndexMap->size() == 0)
    {
        throw InternalError(
                "The index map is empty when writing atomic properties with selection.");
    }
    const int         totalAtoms = selectedAtomsIndexMap.has_value() ? selectedAtomsIndexMap->size()
                                                                     : atomRange.end()->globalAtomNumber();
    std::vector<real> atomCharges;
    atomCharges.reserve(totalAtoms);
    std::vector<real> atomMasses;
    atomMasses.reserve(totalAtoms);
    std::vector<int32_t> atomSpecies;
    atomSpecies.reserve(totalAtoms);
    std::vector<int64_t> atomIds;
    atomIds.reserve(totalAtoms);
    std::vector<int32_t> atomNameTable;
    atomNameTable.reserve(totalAtoms);

    // The lookup table for atom names
    std::vector<std::string> existingAtomNames;
    existingAtomNames.reserve(4096); // Guess a reasonable size for the lookup table
    int       maxStrLength   = 0;
    const int maxGlobalIndex = selectedAtomsIndexMap.has_value()
                                       ? std::max_element(
                                                 selectedAtomsIndexMap->begin(),
                                                 selectedAtomsIndexMap->end(),
                                                 [](const auto& a, const auto& b) {
                                                     return a.first < b.first;
                                                 })->first
                                                 + 1
                                       : atomRange.end()->globalAtomNumber();
    for (auto it = atomRange.begin(); it != atomRange.end(); ++it)
    {
        if (it->globalAtomNumber() > maxGlobalIndex)
        {
            break;
        }
        if (selectedAtomsIndexMap.has_value())
        {
            if (selectedAtomsIndexMap->find(it->globalAtomNumber()) == selectedAtomsIndexMap->end())
            {
                continue;
            }
            // Use 1-based indexing for atom IDs
            atomIds.push_back(selectedAtomsIndexMap->at(it->globalAtomNumber()) + 1);
        }
        else
        {
            // Use 1-based indexing for atom IDs
            atomIds.push_back(it->globalAtomNumber() + 1);
        }

        // Push unique atom names into the string lookup table
        if (std::find(existingAtomNames.begin(), existingAtomNames.end(), it->atomName())
            == existingAtomNames.end())
        {
            existingAtomNames.push_back(it->atomName());
            maxStrLength = std::max(maxStrLength, static_cast<int>(std::strlen(it->atomName())));
        }

        atomCharges.push_back(it->atom().q);
        atomMasses.push_back(it->atom().m);
        atomSpecies.push_back(it->atom().atomnumber);
        atomNameTable.push_back(std::find(existingAtomNames.begin(), existingAtomNames.end(), it->atomName())
                                - existingAtomNames.begin());
    }

    // Create the datasets for atomic properties
    H5mdFixedDataSet<int64_t> datasetAtomId = H5mdDataSetBuilder<int64_t>(baseContainer, c_particleIdName)
                                                      .withDimension({ atomIds.size() })
                                                      .build();
    H5mdFixedDataSet<int32_t> datasetAtomSpecies =
            H5mdDataSetBuilder<int32_t>(baseContainer, c_particleSpeciesName)
                    .withDimension({ atomSpecies.size() })
                    .build();
    H5mdFixedDataSet<real> datasetCharge = H5mdDataSetBuilder<real>(baseContainer, c_particleChargeName)
                                                   .withDimension({ atomCharges.size() })
                                                   .withUnit(c_elementaryChargeUnit)
                                                   .build();
    H5mdFixedDataSet<real> datasetMass = H5mdDataSetBuilder<real>(baseContainer, c_particleMassName)
                                                 .withDimension({ atomMasses.size() })
                                                 .withUnit(c_massUnit)
                                                 .build();
    H5mdFixedDataSet<int32_t> datasetAtomNameInLookupTable =
            H5mdDataSetBuilder<int32_t>(baseContainer, c_particleNameIndicesName)
                    .withDimension({ atomNameTable.size() })
                    .build();
    H5mdFixedDataSet<std::string> datasetAtomNameTable =
            H5mdDataSetBuilder<std::string>(baseContainer, c_particleNameTableName)
                    .withDimension({ existingAtomNames.size() })
                    .withMaxStringLength(maxStrLength + 1)
                    .build();

    // Write cached data to the datasets
    datasetAtomId.writeData(atomIds);
    datasetAtomSpecies.writeData(atomSpecies);
    datasetCharge.writeData(atomCharges);
    datasetMass.writeData(atomMasses);
    datasetAtomNameTable.writeData(existingAtomNames);
    datasetAtomNameInLookupTable.writeData(atomNameTable);

    setAttribute<int64_t>(baseContainer, c_numParticlesAttributeKey, atomIds.size());
}

void writeResidueInfo(AtomRange& atomRange, const hid_t baseContainer, const std::optional<IndexMap>& selectedAtomsIndexMap)
{
    if (selectedAtomsIndexMap.has_value() && selectedAtomsIndexMap->size() == 0)
    {
        throw InternalError(
                "The index map is empty when writing atomic properties with selection.");
    }
    const int totalAtoms = selectedAtomsIndexMap.has_value() ? selectedAtomsIndexMap->size()
                                                             : atomRange.end()->globalAtomNumber();
    std::vector<int32_t> residueIds;
    residueIds.reserve(totalAtoms);
    std::vector<int32_t> residueNames;
    residueNames.reserve(totalAtoms);
    std::vector<int32_t> sequence;
    sequence.reserve(totalAtoms / 3); // A rough estimation of the number of residues (3 atoms per residue)
    std::vector<int32_t> cachedResIDs;
    cachedResIDs.reserve(totalAtoms / 3);

    // Build the residue index map to re-order the residue indices if `selectedAtomsIndexMap` is given
    std::vector<std::string> existingResidueNames;
    existingResidueNames.reserve(1024); // Guess a reasonable size for the lookup table
    std::vector<int32_t> residueIndices(atomRange.end()->globalAtomNumber(), -1);
    int                  residueCount = 0;
    int                  maxStrLength = 0;
    for (auto it = atomRange.begin(); it != atomRange.end(); ++it)
    {
        const int resind = it->residueNumber();

        if (selectedAtomsIndexMap.has_value())
        {
            if (selectedAtomsIndexMap->find(it->globalAtomNumber()) == selectedAtomsIndexMap->end())
            {
                continue;
            }

            // This avoids the expensive find operation for every atom
            if (residueIndices[resind] == -1)
            {
                residueIndices[resind] = residueCount;
                residueCount += 1;
            }

            // Switch the indices to 1-based indexing
            residueIds.push_back(residueIndices[resind] + 1);
        }
        else
        {
            // Switch the indices to 1-based indexing
            residueIds.push_back(resind + 1);
        }

        // Put unique residue names into the string table
        if (std::find(existingResidueNames.begin(), existingResidueNames.end(), it->residueName())
            == existingResidueNames.end())
        {
            existingResidueNames.push_back(it->residueName());
            maxStrLength = std::max(maxStrLength, static_cast<int32_t>(std::strlen(it->residueName())));
        }

        // Get the index of the atom name and put to residue name list and sequence
        const int indexOfResidue =
                std::find(existingResidueNames.begin(), existingResidueNames.end(), it->residueName())
                - existingResidueNames.begin();

        residueNames.push_back(indexOfResidue);
        if (std::find(cachedResIDs.begin(), cachedResIDs.end(), resind) == cachedResIDs.end())
        {
            sequence.push_back(indexOfResidue);
            cachedResIDs.push_back(resind);
        }
    }

    // Create the datasets for residue information
    H5mdFixedDataSet<int32_t> datasetResidueId = H5mdDataSetBuilder<int32_t>(baseContainer, c_residueIdName)
                                                         .withDimension({ residueIds.size() })
                                                         .build();
    H5mdFixedDataSet<int32_t> datasetResidueName =
            H5mdDataSetBuilder<int32_t>(baseContainer, c_residueNameIndicesName)
                    .withDimension({ residueNames.size() })
                    .build();
    H5mdFixedDataSet<int32_t> datasetSequence =
            H5mdDataSetBuilder<int32_t>(baseContainer, c_residueSequenceName)
                    .withDimension({ sequence.size() })
                    .build();
    H5mdFixedDataSet<std::string> datasetResidueNameTable =
            H5mdDataSetBuilder<std::string>(baseContainer, c_residueNameTableName)
                    .withDimension({ existingResidueNames.size() })
                    .withMaxStringLength(maxStrLength + 1)
                    .build();

    // Write cached data to the datasets
    datasetResidueId.writeData(residueIds);
    datasetResidueName.writeData(residueNames);
    datasetSequence.writeData(sequence);
    datasetResidueNameTable.writeData(existingResidueNames);

    setAttribute<int32_t>(baseContainer, c_numResiduesAttributeKey, sequence.size());
}

void writeBonds(const gmx_mtop_t& topology, const hid_t baseContainer, const std::optional<IndexMap>& selectedAtomsIndexMap)
{
    if (selectedAtomsIndexMap.has_value() && selectedAtomsIndexMap->size() == 0)
    {
        throw InternalError(
                "The index map is empty when writing atomic properties with selection.");
    }
    // NOTE: Initialize the bond container and reserve enough space, assuming the
    //       number of bonds is smaller than the number of atoms
    std::vector<int64_t> systemBonds;
    if (selectedAtomsIndexMap.has_value())
    {
        systemBonds.reserve(selectedAtomsIndexMap->size() * 2);
    }
    else
    {
        systemBonds.reserve(topology.natoms * 2);
    }
    int64_t globalOffset = 0;
    for (size_t i = 0; i < topology.molblock.size(); i++)
    {
        const gmx_molblock_t& molBlock = topology.molblock[i];
        const gmx_moltype_t&  molType  = topology.moltype[molBlock.type];
        const size_t          numMols  = molBlock.nmol;
        const size_t          numAtoms = molType.atoms.nr;

        // Deduce bonds from the interaction list of the molecule type
        const BondPairs bonds = deduceBondsFromMolecule(molType);
        if (selectedAtomsIndexMap.has_value())
        {
            cacheBondForNBlocks(
                    bonds, &systemBonds, numMols, numAtoms, selectedAtomsIndexMap.value(), globalOffset);
        }
        else
        {
            cacheBondForNBlocks(bonds, &systemBonds, numMols, numAtoms, std::nullopt, globalOffset);
        }

        // Global index offset for the next molecule type
        globalOffset += numMols * molType.atoms.nr;
    }

    // Register the Root group explicitly for connectivity
    const size_t bCount = systemBonds.size() / 2;
    H5mdFixedDataSet<int64_t> dBonds = H5mdDataSetBuilder<int64_t>(baseContainer, c_bondsConnectivityName)
                                               .withDimension({ bCount, 2 })
                                               .build();
    dBonds.writeData(systemBonds);
    setAttribute<int64_t>(baseContainer, c_numberOfBondsAttributeKey, bCount);
}

void writeDisulfideBonds(const gmx_mtop_t&              topology,
                         const hid_t                    baseContainer,
                         const std::optional<IndexMap>& selectedAtomsIndexMap)
{
    if (selectedAtomsIndexMap.has_value() && selectedAtomsIndexMap->size() == 0)
    {
        throw InternalError(
                "The index map is empty when writing atomic properties with selection.");
    }
    // Obtain the disulfide bonds from the topology
    std::vector<int64_t> systemBonds;
    int64_t              globalOffset = 0;
    for (size_t i = 0; i < topology.molblock.size(); i++)
    {
        const gmx_molblock_t& molBlock = topology.molblock[i];
        const gmx_moltype_t&  molType  = topology.moltype[molBlock.type];
        const size_t          numMols  = molBlock.nmol;
        const size_t          numAtoms = molType.atoms.nr;

        const BondPairs disulfideBonds = deduceDisulfideBondsFromMolecule(molType);
        if (selectedAtomsIndexMap.has_value())
        {
            cacheBondForNBlocks(
                    disulfideBonds, &systemBonds, numMols, numAtoms, selectedAtomsIndexMap.value(), globalOffset);
        }
        else
        {
            cacheBondForNBlocks(disulfideBonds, &systemBonds, numMols, numAtoms, std::nullopt, globalOffset);
        }

        // Set the global index offset for the next molecule type
        globalOffset += numMols * molType.atoms.nr;
    }

    // Write the information of disulfide bonds
    if (systemBonds.size() > 0)
    {
        const size_t bondCount = systemBonds.size() / 2;
        H5mdFixedDataSet<int64_t> dBonds = H5mdDataSetBuilder<int64_t>(baseContainer, c_disulfideBondsName)
                                                   .withDimension({ bondCount, 2 })
                                                   .build();
        dBonds.writeData(systemBonds);
        setAttribute<int64_t>(baseContainer, c_numberOfBondsAttributeKey, bondCount);
    }
}

void labelInternalTopologyVersion(const hid_t baseContainer)
{
    // Similar to H5MD standard, the version is a simple array of two integers.
    setAttributeVector<int>(
            baseContainer,
            c_h5mdInternalTopologyVersionAttributeKey,
            std::vector<int>{ c_h5mdInternalTopologyVersionMajor, c_h5mdInternalTopologyVersionMinor });
}

void labelTopologyName(const hid_t baseContainer, const char* topName)
{
    setAttribute(baseContainer, c_h5mdInternalTopologyNameAttributeKey, topName);
}

void writeMoleculeTypes(const hid_t baseContainer, const ArrayRef<const gmx_moltype_t> moltypes)
{
    const int numMolTypes = moltypes.size();
    // Nothing to write
    if (moltypes.empty())
    {
        return;
    }

    std::vector<std::string> moltypeNames(numMolTypes);
    for (int index = 0; index < numMolTypes; ++index)
    {
        const gmx_moltype_t& moltype     = moltypes[index];
        const std::string&   moltypeName = *(moltype.name);
        moltypeNames[index]              = moltypeName;
        const auto [molContainer, molGuard] =
                makeH5mdGroupGuard(createGroup(baseContainer, moltypeName.c_str()));

        // Create an AtomRange to iterate over the molecule type
        gmx_mtop_t tempMtop;
        detail::mtopFromMolType(&tempMtop, moltype);
        AtomRange atomRange(tempMtop);

        // Write atomic properties information
        writeAtomicProperties(atomRange, molContainer);

        // Write sequence/residue information
        writeResidueInfo(atomRange, molContainer);

        // TODO: Write interaction list and exclusions list - MR !5524/!5523
    }
    setAttributeVector(baseContainer, c_moleculeNamesAttributeKey, moltypeNames);
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
