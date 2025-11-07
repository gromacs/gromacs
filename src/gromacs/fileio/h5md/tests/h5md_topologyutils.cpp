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

/*! \internal \file
 * \brief
 * Tests for topology-related utility functions.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_topologyutils.h"

#include <hdf5.h>

#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/tprfilegenerator.h"

namespace gmx
{
namespace test
{
namespace
{

inline std::vector<std::string> readFixedStringDataset(hid_t baseContainer, const char* datasetName)
{
    // Create the memory dataspace for string reading
    H5mdDataSetBase<std::string> dataset(baseContainer, datasetName);
    auto [dataSpace, dataSpaceGuard]     = makeH5mdDataSpaceGuard(H5Dget_space(dataset.id()));
    const auto [memSpace, memSpaceGuard] = makeH5mdDataSpaceGuard(
            H5Screate_simple(dataset.dims().size(), dataset.dims().data(), nullptr));

    // Read the data from the dataset
    int               nrStrings = static_cast<int>(dataset.dims()[0]);
    int               strLength = H5Tget_size(dataset.dataType());
    std::vector<char> buffer(nrStrings * strLength, '\0');
    if (H5Dread(dataset.id(), dataset.dataType(), memSpace, dataSpace, H5P_DEFAULT, buffer.data()) < 0)
    {
        throw FileIOError("Failed to read string data from dataset.");
    }

    // Unpack the string buffer into vector of strings
    std::vector<std::string> retStringArray(nrStrings);
    for (int i = 0; i < nrStrings; ++i)
    {
        retStringArray[i] = std::string(&buffer[i * strLength]);
    }
    return retStringArray;
}

void createDisulfideBondTopology(gmx_mtop_t*                topology,
                                 const ArrayRef<char*>&     atomNames,
                                 const ArrayRef<const int>& resIndices,
                                 const BondPairs&           disulfideBonds,
                                 const int                  numMols)
{
    topology->moltype.resize(1);
    topology->molblock.resize(1);

    gmx_moltype_t& molType = topology->moltype[0];
    molType.atoms.nr       = static_cast<int>(atomNames.size());
    snew(molType.atoms.atom, molType.atoms.nr);
    snew(molType.atoms.atomname, molType.atoms.nr);

    auto atomNameIterator = atomNames.begin();
    for (int i = 0; i < molType.atoms.nr; ++i, ++atomNameIterator)
    {
        molType.atoms.atom[i].resind = resIndices[i];
        molType.atoms.atomname[i]    = atomNameIterator.data();
        if (std::strcmp(*atomNameIterator, "SG") == 0)
        {
            molType.atoms.atom[i].atomnumber = 16;
        }
        else
        {
            molType.atoms.atom[i].atomnumber = 6;
        }
    }

    topology->molblock[0].type = 0;
    topology->molblock[0].nmol = numMols;

    for (const auto& bond : disulfideBonds)
    {
        molType.ilist[InteractionFunction::Bonds].push_back(
                0, std::array<int, 2>({ static_cast<int>(bond.first), static_cast<int>(bond.second) }));
    }
}

using H5mdTopologyUtilTest = H5mdTestBase;

TEST_F(H5mdTopologyUtilTest, WriteFullTopologyAtomProp)
{
    // Prepare the topology to write
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    topology.molblock[1].nmol = 10;
    int natoms                = 0;
    for (size_t i = 0; i < topology.moltype.size(); ++i)
    {
        natoms += topology.moltype[i].atoms.nr * topology.molblock[i].nmol;
    }
    topology.natoms = natoms;
    AtomRange atomRange(topology);

    {
        SCOPED_TRACE("Apply no selection and by default writes all atoms");

        const auto [baseContainer, baseContainerGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/systemNoSelection"));
        writeAtomicProperties(atomRange, baseContainer);
    }

    {
        SCOPED_TRACE("Select all atoms and write their properties");

        // Prepare the index group to select all atoms
        std::vector<int> indices(topology.natoms, 0);
        std::iota(indices.begin(), indices.end(), 0);
        IndexGroup selectedGroup = { "systemSelectAll", indices };
        IndexMap selectedAtomsIndexMap = mapSelectionToInternalIndices(selectedGroup.particleIndices);

        const auto [baseContainer, baseContainerGuard] = makeH5mdGroupGuard(createGroup(
                fileid(), formatString("/particles/%s", selectedGroup.name.c_str()).c_str()));
        writeAtomicProperties(atomRange, baseContainer, selectedAtomsIndexMap);
    }

    {
        SCOPED_TRACE("Check the correctness of the dumped atomic properties");

        // Prepare the reference data checker
        TestReferenceData    data;
        TestReferenceChecker checker(data.rootChecker());
        checker.setDefaultTolerance(test::absoluteTolerance(0.1));

        for (auto& name : { "systemNoSelection", "systemSelectAll" })
        {
            const auto [baseContainer, baseContainerGuard] =
                    makeH5mdGroupGuard(openGroup(fileid(), formatString("/particles/%s", name).c_str()));
            // Check the correctness of the dumped atomic properties
            auto retNrAtoms = getAttribute<int64_t>(baseContainer, "nr_particles");
            ASSERT_TRUE(retNrAtoms.has_value());
            // 1 protein (29) + 10 water (30) = 59
            ASSERT_EQ(retNrAtoms.value(), 59);

            // Match the reordered atom indices
            std::vector<int64_t> expectedID(retNrAtoms.value());
            std::iota(expectedID.begin(), expectedID.end(), 1);
            H5mdFixedDataSet<int64_t> dataAtomID = H5mdFixedDataSet<int64_t>(baseContainer, "id");
            std::vector<int64_t>      retAtomID(dataAtomID.numValues());
            dataAtomID.readData(retAtomID);
            EXPECT_EQ(retAtomID, expectedID);

            std::vector<std::string> retNameStringTable =
                    readFixedStringDataset(baseContainer, "particle_name_table");
            checker.checkSequence(makeConstArrayRef(retNameStringTable), "ParticleNameTable");

            H5mdFixedDataSet<int32_t> dataAtomicNumber =
                    H5mdFixedDataSet<int32_t>(baseContainer, "species");
            std::vector<int32_t> retAtomicNumber(dataAtomicNumber.numValues());
            dataAtomicNumber.readData(retAtomicNumber);
            checker.checkSequence(makeConstArrayRef(retAtomicNumber), "ParticleSpecies");

            H5mdFixedDataSet<real> dataCharges = H5mdFixedDataSet<real>(baseContainer, "charge");
            std::vector<real>      retCharges(dataCharges.numValues());
            dataCharges.readData(retCharges);
            checker.checkSequence(makeConstArrayRef(retCharges), "ParticleCharge");

            H5mdFixedDataSet<real> dataMasses = H5mdFixedDataSet<real>(baseContainer, "mass");
            std::vector<real>      retMasses(dataMasses.numValues());
            dataMasses.readData(retMasses);
            checker.checkSequence(makeConstArrayRef(retMasses), "ParticleMass");
        }
    }
}

TEST_F(H5mdTopologyUtilTest, WriteFullTopologyResInfo)
{
    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    // Only keep 10 water molecules to reduce the size of refdata
    topology.molblock[1].nmol = 10;
    int natoms                = 0;
    for (size_t i = 0; i < topology.moltype.size(); ++i)
    {
        natoms += topology.moltype[i].atoms.nr * topology.molblock[i].nmol;
    }
    topology.natoms = natoms;
    AtomRange atomRange(topology);

    {
        SCOPED_TRACE("Apply no selection and by default writes all sequence and residues");

        const auto [baseContainer, baseContainerGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/particles/systemNoSelection"));
        writeResidueInfo(atomRange, baseContainer);
    }

    {
        SCOPED_TRACE("Select all atoms and write their properties");

        // Prepare the index group to select all atoms
        std::vector<int> indices(topology.natoms, 0);
        std::iota(indices.begin(), indices.end(), 0);
        IndexGroup selectedGroup = { "systemSelectAll", indices };
        IndexMap selectedAtomsIndexMap = mapSelectionToInternalIndices(selectedGroup.particleIndices);

        const auto [baseContainer, baseContainerGuard] = makeH5mdGroupGuard(createGroup(
                fileid(), formatString("/particles/%s", selectedGroup.name.c_str()).c_str()));
        writeResidueInfo(atomRange, baseContainer, selectedAtomsIndexMap);
    }

    {
        SCOPED_TRACE("Checking the correctness of the dumped residue information");

        // Prepare the reference data checker
        TestReferenceData    data;
        TestReferenceChecker checker(data.rootChecker());
        checker.setDefaultTolerance(test::absoluteTolerance(0.1));

        for (auto& name : { "systemNoSelection", "systemSelectAll" })
        {
            const auto [baseContainer, baseContainerGuard] =
                    makeH5mdGroupGuard(openGroup(fileid(), formatString("/particles/%s", name).c_str()));

            auto retNrResidues = getAttribute<int32_t>(baseContainer, "nr_residues");
            ASSERT_TRUE(retNrResidues.has_value());
            // 1 protein (2) + 10 water (10)
            ASSERT_EQ(retNrResidues.value(), 12);

            std::vector<std::string> retResidueNameTable =
                    readFixedStringDataset(baseContainer, "residue_name_table");
            checker.checkSequence(makeConstArrayRef(retResidueNameTable), "ResidueNameTable");

            H5mdFixedDataSet<int32_t> dataSequence =
                    H5mdFixedDataSet<int32_t>(baseContainer, "sequence");
            std::vector<int32_t> retSequence(dataSequence.numValues());
            dataSequence.readData(retSequence);
            checker.checkSequence(makeConstArrayRef(retSequence), "Sequence");
        }
    }
}

TEST_F(H5mdTopologyUtilTest, WriteThreeSelectedWater)
{
    // Prepare the reference data checker
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(test::absoluteTolerance(0.1));

    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);
    AtomRange atomRange(topology);

    std::vector<int> indices       = { 29,  30,  31,  218, 219, 220, 314, 315,
                                       316, 566, 567, 568, 860, 861, 862 };
    IndexGroup       selectedGroup = { "random_water", indices };
    IndexMap selectedAtomsIndexMap = mapSelectionToInternalIndices(selectedGroup.particleIndices);

    const auto [baseContainer, baseContainerGuard] = makeH5mdGroupGuard(createGroup(
            fileid(), formatString("/particles/%s", selectedGroup.name.c_str()).c_str()));

    {
        writeResidueInfo(atomRange, baseContainer, selectedAtomsIndexMap);

        auto retNrResidues = getAttribute<int32_t>(baseContainer, "nr_residues");
        ASSERT_TRUE(retNrResidues.has_value());
        EXPECT_EQ(retNrResidues.value(), 5);

        H5mdFixedDataSet<int32_t> datasetResidueIDs =
                H5mdFixedDataSet<int32_t>(baseContainer, "residue_id");
        std::vector<int32_t> retResidueIDs(datasetResidueIDs.numValues());
        datasetResidueIDs.readData(retResidueIDs);
        checker.checkSequence(makeConstArrayRef(retResidueIDs), "ResidueIDs");

        // All atom name should 0 and direct to SOL
        H5mdFixedDataSet<int32_t> dataResidueNames =
                H5mdFixedDataSet<int32_t>(baseContainer, "residue_name");
        std::vector<int32_t> retResidueNames(dataResidueNames.numValues());
        dataResidueNames.readData(retResidueNames);
        for (size_t i = 0; i < retResidueNames.size(); ++i)
        {
            EXPECT_EQ(retResidueNames[i], 0);
        }

        std::vector<std::string> retResidueNameTable =
                readFixedStringDataset(baseContainer, "residue_name_table");
        checker.checkSequence(makeConstArrayRef(retResidueNameTable), "ResidueNameTable");

        // Increment residue ID every 3 atoms
        H5mdFixedDataSet<int32_t> dataResidueID =
                H5mdFixedDataSet<int32_t>(baseContainer, "residue_id");
        std::vector<int32_t> retResidueID(dataResidueID.numValues());
        dataResidueID.readData(retResidueID);
        for (size_t i = 0; i < retResidueID.size(); ++i)
        {
            EXPECT_EQ(retResidueID[i], (i / 3) + 1);
        }
    }


    {
        SCOPED_TRACE("Test writing atomic properties");

        writeAtomicProperties(atomRange, baseContainer, selectedAtomsIndexMap);

        // Match the number of particles
        auto retNrAtoms = getAttribute<int64_t>(baseContainer, "nr_particles");
        ASSERT_TRUE(retNrAtoms.has_value());
        EXPECT_EQ(retNrAtoms.value(), 15);

        // Match the reordered atom indices
        std::vector<int64_t> expectedID(15);
        std::iota(expectedID.begin(), expectedID.end(), 1);

        H5mdFixedDataSet<int64_t> dataAtomID = H5mdFixedDataSet<int64_t>(baseContainer, "id");
        std::vector<int64_t>      retAtomID(dataAtomID.numValues());
        dataAtomID.readData(retAtomID);
        EXPECT_EQ(retAtomID, expectedID);

        // Match the number of possible unique atom names in water molecule
        std::vector<std::string> retNameStringTable =
                readFixedStringDataset(baseContainer, "particle_name_table");
        checker.checkSequence(makeConstArrayRef(retNameStringTable), "ParticleNameTable");

        H5mdFixedDataSet<int32_t> dataAtomicNumber =
                H5mdFixedDataSet<int32_t>(baseContainer, "species");
        std::vector<int32_t> retAtomicNumber(dataAtomicNumber.numValues());
        dataAtomicNumber.readData(retAtomicNumber);
        checker.checkSequence(makeConstArrayRef(retAtomicNumber), "ParticleSpecies");

        H5mdFixedDataSet<real> dataCharges = H5mdFixedDataSet<real>(baseContainer, "charge");
        std::vector<real>      retCharges(dataCharges.numValues());
        dataCharges.readData(retCharges);
        checker.checkSequence(makeConstArrayRef(retCharges), "ParticleCharge");

        H5mdFixedDataSet<real> dataMasses = H5mdFixedDataSet<real>(baseContainer, "mass");
        std::vector<real>      retMasses(dataMasses.numValues());
        dataMasses.readData(retMasses);
        checker.checkSequence(makeConstArrayRef(retMasses), "ParticleMass");
    }
}


TEST_F(H5mdTopologyUtilTest, WriteProteinTopology)
{
    // Prepare the reference data checker
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(test::absoluteTolerance(0.1));

    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);
    AtomRange atomRange(topology);

    const int        protNumAtoms = 29;
    std::vector<int> indices(protNumAtoms);
    std::iota(indices.begin(), indices.end(), 0);
    IndexGroup selectedGroup = { "protein", indices };

    IndexMap selectedAtomsIndexMap = mapSelectionToInternalIndices(selectedGroup.particleIndices);
    const auto [baseContainer, baseContainerGuard] = makeH5mdGroupGuard(createGroup(
            fileid(), formatString("/particles/%s", selectedGroup.name.c_str()).c_str()));
    const auto [connectivity, connectivityGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/connectivity/bond"));

    {
        SCOPED_TRACE("Test writing residue information and sequence");

        writeResidueInfo(atomRange, baseContainer, selectedAtomsIndexMap);

        auto retNrResidues = getAttribute<int32_t>(baseContainer, "nr_residues");
        ASSERT_TRUE(retNrResidues.has_value());
        EXPECT_EQ(retNrResidues.value(), 2);

        std::vector<std::string> retResidueNameTable =
                readFixedStringDataset(baseContainer, "residue_name_table");
        checker.checkSequence(makeConstArrayRef(retResidueNameTable), "ResidueNameTable");

        H5mdFixedDataSet<int32_t> dataSequence =
                H5mdFixedDataSet<int32_t>(baseContainer, "sequence");
        std::vector<int32_t> retSequence(dataSequence.numValues());
        dataSequence.readData(retSequence);
        checker.checkSequence(makeConstArrayRef(retSequence), "Sequence");
    }

    {
        SCOPED_TRACE("Test writing atomic properties");

        writeAtomicProperties(atomRange, baseContainer, selectedAtomsIndexMap);

        auto retNrAtoms = getAttribute<int64_t>(baseContainer, "nr_particles");
        ASSERT_TRUE(retNrAtoms.has_value());
        EXPECT_EQ(retNrAtoms.value(), protNumAtoms);

        // Match the reordered atom indices
        std::vector<int64_t> expectedID(retNrAtoms.value());
        std::iota(expectedID.begin(), expectedID.end(), 1);
        H5mdFixedDataSet<int64_t> dataAtomID = H5mdFixedDataSet<int64_t>(baseContainer, "id");
        std::vector<int64_t>      retAtomID(dataAtomID.numValues());
        dataAtomID.readData(retAtomID);
        EXPECT_EQ(retAtomID, expectedID);

        std::vector<std::string> retNameStringTable =
                readFixedStringDataset(baseContainer, "particle_name_table");
        checker.checkSequence(makeConstArrayRef(retNameStringTable), "ParticleNameTable");

        H5mdFixedDataSet<int32_t> dataAtomicNumber =
                H5mdFixedDataSet<int32_t>(baseContainer, "species");
        std::vector<int32_t> retAtomicNumber(dataAtomicNumber.numValues());
        dataAtomicNumber.readData(retAtomicNumber);
        checker.checkSequence(makeConstArrayRef(retAtomicNumber), "ParticleSpecies");

        H5mdFixedDataSet<real> dataCharges = H5mdFixedDataSet<real>(baseContainer, "charge");
        std::vector<real>      retCharges(dataCharges.numValues());
        dataCharges.readData(retCharges);
        checker.checkSequence(makeConstArrayRef(retCharges), "ParticleCharge");

        H5mdFixedDataSet<real> dataMasses = H5mdFixedDataSet<real>(baseContainer, "mass");
        std::vector<real>      retMasses(dataMasses.numValues());
        dataMasses.readData(retMasses);
        checker.checkSequence(makeConstArrayRef(retMasses), "ParticleMass");
    }
}

TEST_F(H5mdTopologyUtilTest, ThrowUponEmptyIndexMap)
{
    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);
    AtomRange atomRange(topology);
    const auto [baseContainer, baseContainerGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/particles/empty_selection"));

    EXPECT_THROW(writeAtomicProperties(atomRange, baseContainer, IndexMap({})), InternalError);
    EXPECT_THROW(writeResidueInfo(atomRange, baseContainer, IndexMap({})), InternalError);
}

TEST_F(H5mdTopologyUtilTest, WriteConnectivityProteinPart)
{
    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    const int        protNrAtoms = 29;
    std::vector<int> indices(protNrAtoms);
    std::iota(indices.begin(), indices.end(), 0);

    // Prepare the reference data checker
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(test::absoluteTolerance(0.1));

    const auto [connectivity, connectivityGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/connectivity/bond"));

    writeBonds(topology, connectivity, mapSelectionToInternalIndices(indices));
    const auto nrBonds = getAttribute<int64_t>(connectivity, "nr_bonds");
    ASSERT_TRUE(nrBonds.has_value());
    EXPECT_EQ(nrBonds.value(), 22);

    H5mdFixedDataSet<int64_t> datasetBonds = H5mdFixedDataSet<int64_t>(connectivity, "bonds");
    std::vector<int64_t>      bondData(datasetBonds.numValues());
    datasetBonds.readData(bondData);
    EXPECT_EQ(bondData.size(), 2 * nrBonds.value());
    checker.checkSequence(makeConstArrayRef(bondData), "Bonds");
}


TEST_F(H5mdTopologyUtilTest, WriteConnectivityRandomSelectedWater)
{
    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    const std::vector<int> indices = { 29,  30,  31,  218, 219, 220, 314, 315,
                                       316, 566, 567, 568, 860, 861, 862 };

    const auto [connectivity, connectivityGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/connectivity/bond"));

    // Prepare the reference data checker
    TestReferenceData    data;
    TestReferenceChecker checker(data.rootChecker());
    checker.setDefaultTolerance(test::absoluteTolerance(0.1));

    writeBonds(topology, connectivity, mapSelectionToInternalIndices(indices));

    const auto nrBonds = getAttribute<int64_t>(connectivity, "nr_bonds");
    ASSERT_TRUE(nrBonds.has_value());
    EXPECT_EQ(nrBonds.value(), 10);

    H5mdFixedDataSet<int64_t> datasetBonds = H5mdFixedDataSet<int64_t>(connectivity, "bonds");
    std::vector<int64_t>      bondData(datasetBonds.numValues());
    datasetBonds.readData(bondData);
    EXPECT_EQ(bondData.size(), 2 * nrBonds.value());
    checker.checkSequence(makeConstArrayRef(bondData), "Bonds");
}

TEST_F(H5mdTopologyUtilTest, WriteDisulfideBonds)
{
    // Prepare the topology for disulfide bond writing
    const std::vector<std::string> atomNames = {
        "AA", "BB",       // Dummy atoms
        "SG", "CB", "CA", // CYS 1-4
        "AA", "BB",       // Dummy atoms
        "SG", "CB", "CA", // CYS 2-3
        "AA", "BB",       // Dummy atoms
        "SG", "CB", "CA", // 2-3
        "SG", "CB", "CA", // 1-4
        "AA", "BB",       // Dummy atoms
        "SG", "CB", "CA", // Not forming a disulfide bond here
        "SG", "CB", "CA", // CYS 5-6
        "SG", "CB", "CA", // 5-6
        "AA", "BB",       // Dummy atoms
    };
    std::vector<char*> charAtomNames;
    for (const auto& atomName : atomNames)
    {
        charAtomNames.push_back(const_cast<char*>(atomName.data()));
    }

    const std::vector<int> resinds = { 1, 1, 2, 2, 2, 3, 3, 4,  4,  4,  5,  5,  6,  6,  6, 7,
                                       7, 7, 8, 8, 9, 9, 9, 10, 10, 10, 11, 11, 11, 12, 12 };
    std::vector<std::pair<int64_t, int64_t>> expectedBonds = {
        { 2, 15 },  // Disulfide CYS 1-4
        { 7, 12 },  // Disulfide CYS 2-3
        { 23, 26 }, // Disulfide CYS 5-6
    };

    {
        SCOPED_TRACE("Create the disulfide bonds and check connectivity");

        gmx_mtop_t topology;
        createDisulfideBondTopology(&topology, charAtomNames, resinds, expectedBonds, 1);

        const auto [connectivity, connectivityGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/connectivity/special_bonds"));
        writeDisulfideBonds(topology, connectivity);

        H5mdFixedDataSet<int64_t> datasetBond =
                H5mdFixedDataSet<int64_t>(connectivity, "disulfide_bonds");
        std::vector<int64_t> retBonds(datasetBond.numValues());
        datasetBond.readData(retBonds);
        EXPECT_EQ(retBonds.size(), expectedBonds.size() * 2);
        for (size_t i = 0; i < expectedBonds.size(); ++i)
        {
            EXPECT_EQ(retBonds[2 * i], expectedBonds[i].first);
            EXPECT_EQ(retBonds[2 * i + 1], expectedBonds[i].second);
        }
    }

    {
        SCOPED_TRACE("Write multiple replication of the block that contains disulfide bonds");

        const int  blockReplication = 4;
        gmx_mtop_t topology;
        createDisulfideBondTopology(&topology, charAtomNames, resinds, expectedBonds, blockReplication);

        const auto [connectivity, connectivityGuard] =
                makeH5mdGroupGuard(createGroup(fileid(), "/connectivity_rep_blocks"));
        writeDisulfideBonds(topology, connectivity);

        H5mdFixedDataSet<int64_t> datasetBond =
                H5mdFixedDataSet<int64_t>(connectivity, "disulfide_bonds");
        std::vector<int64_t> retBonds(datasetBond.numValues());
        datasetBond.readData(retBonds);

        EXPECT_EQ(retBonds.size(), expectedBonds.size() * 2 * blockReplication);
        for (size_t i = 0; i < expectedBonds.size(); ++i)
        {
            for (int b = 0; b < blockReplication; ++b)
            {
                const size_t offset = b * topology.moltype[0].atoms.nr;
                const size_t index  = b * expectedBonds.size() * 2 + i * 2;
                EXPECT_EQ(retBonds[index], expectedBonds[i].first + offset);
                EXPECT_EQ(retBonds[index + 1], expectedBonds[i].second + offset);
            }
        }
    }
}

TEST_F(H5mdTopologyUtilTest, CreateMTopFromMolType)
{
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    for (const auto& moltype : topology.moltype)
    {
        gmx_mtop_t newMtop;
        gmx::detail::mtopFromMolType(&newMtop, moltype);

        ASSERT_EQ(moltype.atoms.nr, newMtop.moltype[0].atoms.nr);
        for (int i = 0; i < moltype.atoms.nr; ++i)
        {
            EXPECT_EQ(moltype.atoms.atom[i].type, newMtop.moltype[0].atoms.atom[i].type);
            EXPECT_EQ(moltype.atoms.atom[i].m, newMtop.moltype[0].atoms.atom[i].m);
            EXPECT_EQ(moltype.atoms.atom[i].q, newMtop.moltype[0].atoms.atom[i].q);
            EXPECT_EQ(moltype.atoms.atom[i].resind, newMtop.moltype[0].atoms.atom[i].resind);
            EXPECT_EQ(moltype.atoms.atom[i].atomnumber, newMtop.moltype[0].atoms.atom[i].atomnumber);
            EXPECT_EQ(std::string(moltype.atoms.atom[i].elem),
                      std::string(newMtop.moltype[0].atoms.atom[i].elem));
        }
    }
}

TEST_F(H5mdTopologyUtilTest, LabelVersionH5MDMTop)
{
    const auto [topologyContainer, topologyGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/h5md/modules/gromacs_topology"));
    const std::vector<int> expectedVersion = { 0, 1 };

    // Set headers for the internal topology
    labelInternalTopologyVersion(topologyContainer);

    // Read back the version and compare with the internal version constant
    const auto version = getAttributeVector<int32_t>(topologyContainer, "version");
    ASSERT_TRUE(version.has_value());
    EXPECT_EQ(version.value(), expectedVersion);
}

TEST_F(H5mdTopologyUtilTest, LabelSystemName)
{
    const auto [topologyContainer, topologyGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/h5md/modules/gromacs_topology"));
    const std::string systemName = "TestSystem";

    // Set headers for the internal topology
    labelTopologyName(topologyContainer, systemName.c_str());

    // Read back the version and compare with the internal version constant
    const auto version = getAttribute<std::string>(topologyContainer, "system_name");
    ASSERT_TRUE(version.has_value());
    EXPECT_EQ(version.value(), systemName);
}

TEST_F(H5mdTopologyUtilTest, WritesMoleculeTypes)
{
    // Prepare the topology of the molecule
    gmx_mtop_t        topology;
    bool              haveTopology;
    TprAndFileManager tprFileHandle("alanine_vsite_solvated");
    readConfAndTopology(
            tprFileHandle.tprName(), &haveTopology, &topology, nullptr, nullptr, nullptr, nullptr);

    const auto [gmxMol, gmxMolGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/h5md/modules/gromacs_topology"));

    // Write the molecule type and block information
    writeMoleculeTypes(gmxMol, makeConstArrayRef(topology.moltype));

    {
        // Prepare the reference data checker
        TestReferenceData    data;
        TestReferenceChecker checker(data.rootChecker());
        checker.setDefaultTolerance(test::absoluteTolerance(0.1));

        // Read back and check the molecule type names
        const std::vector<std::string> expectedMolNames = { "Alanine_dipeptide", "SOL" };

        const auto retMoleculeNames = getAttributeVector<std::string>(gmxMol, "molecule_names");
        ASSERT_TRUE(retMoleculeNames.has_value());
        ASSERT_EQ(retMoleculeNames.value(), expectedMolNames);

        for (const std::string& molName : expectedMolNames)
        {
            // charge - real
            const std::string      chargeName = (molName + "/charge");
            H5mdFixedDataSet<real> datasetCharge(gmxMol, chargeName.c_str());
            std::vector<real>      dataCharge(datasetCharge.numValues());
            datasetCharge.readData(dataCharge);
            checker.checkSequence(makeConstArrayRef(dataCharge), chargeName.c_str());

            // mass - real
            const std::string      massName = (molName + "/mass");
            H5mdFixedDataSet<real> datasetMass(gmxMol, massName.c_str());
            std::vector<real>      dataMass(datasetMass.numValues());
            datasetMass.readData(dataMass);
            checker.checkSequence(makeConstArrayRef(dataMass), massName.c_str());

            // particle_name - int
            const std::string     particleName = (molName + "/particle_name");
            H5mdFixedDataSet<int> datasetParticleName(gmxMol, particleName.c_str());
            std::vector<int>      dataParticleName(datasetParticleName.numValues());
            datasetParticleName.readData(dataParticleName);
            checker.checkSequence(makeConstArrayRef(dataParticleName), particleName.c_str());

            // particle_name_table - std::string
            const std::string        particleNameTable = (molName + "/particle_name_table");
            std::vector<std::string> retNameStringTable =
                    readFixedStringDataset(gmxMol, particleNameTable.c_str());
            checker.checkSequence(makeConstArrayRef(retNameStringTable), particleNameTable.c_str());

            // residue_name - int
            const std::string     residueName = (molName + "/residue_name");
            H5mdFixedDataSet<int> datasetResidueName(gmxMol, residueName.c_str());
            std::vector<int>      dataResidueName(datasetResidueName.numValues());
            datasetResidueName.readData(dataResidueName);
            checker.checkSequence(makeConstArrayRef(dataResidueName), residueName.c_str());

            // residue_name_table - std::string
            const std::string        residueNameTable = (molName + "/residue_name_table");
            std::vector<std::string> retResidueNameTable =
                    readFixedStringDataset(gmxMol, residueNameTable.c_str());
            checker.checkSequence(makeConstArrayRef(retResidueNameTable), residueNameTable.c_str());
        }
        // TODO: further check the detailed data when the actual writing functions are merged
    }
}

TEST_F(H5mdTopologyUtilTest, WriteEmptyMoleculeTypes)
{
    // Prepare the topology of the molecule
    gmx_mtop_t topology;
    const auto [gmxMol, gmxMolGuard] =
            makeH5mdGroupGuard(createGroup(fileid(), "/h5md/modules/gromacs_topology"));

    // Write the molecule type and block information
    writeMoleculeTypes(gmxMol, makeConstArrayRef(topology.moltype));

    EXPECT_EQ(getAttributeVector<std::string>(gmxMol, "molecule_names"), std::nullopt);
}

} // namespace
} // namespace test
} // namespace gmx
