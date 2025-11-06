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

/*! \brief Declares the GROMACS/H5MD i/o interface.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_topologyutils.h"

#include <hdf5.h>

#include <filesystem>
#include <memory>
#include <numeric>

#include <gtest/gtest.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/h5md/h5md_attribute.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/include/gromacs/topology/mtop_util.h"
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

// Fixture for tests that does not need a tpr file
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

    const int        protNrAtoms = 29;
    std::vector<int> indices(protNrAtoms);
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
        EXPECT_EQ(retNrAtoms.value(), protNrAtoms);

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

} // namespace
} // namespace test
} // namespace gmx
