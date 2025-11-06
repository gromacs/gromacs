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
#include "gromacs/fileio/h5md/h5md_datasetbase.h"
#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/real.h"

// HDF5 constants use old style casts.
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")

namespace gmx
{

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
    sequence.reserve(totalAtoms / 3); // A rought estimation of the number of residues (3 atoms per residue)
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

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
