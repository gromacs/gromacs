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

/*! \internal \file
 * \brief Tests for writing group information to H5MD files and reading it back.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_indexgroup.h"

#include <hdf5.h>

#include <climits>

#include <algorithm>
#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/h5md/h5md_datasetbuilder.h"
#include "gromacs/fileio/h5md/h5md_fixeddataset.h"
#include "gromacs/fileio/h5md/h5md_group.h"
#include "gromacs/fileio/h5md/h5md_guard.h"
#include "gromacs/fileio/h5md/h5md_util.h"
#include "gromacs/fileio/h5md/tests/h5mdtestbase.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

using H5mdIndexGroupTest = H5mdTestBase;

//! \brief An int32 data set exactly as it is stored in the file.
struct RawDataSet
{
    //! Dimensions of the data set.
    DataSetDims dims;
    //! All values of the data set, in row-major order.
    std::vector<int> values;
};

/*! \brief Read the int32 data set \p name in \p container without any interpretation.
 *
 * \param[in] container HDF5 container holding the data set.
 * \param[in] name      Name of the data set to read.
 * \return The dimensions and raw values of the data set.
 */
RawDataSet readRawDataSet(const hid_t container, const char* name)
{
    const H5mdFixedDataSet<int> dataSet(container, name);
    RawDataSet                  raw;
    raw.dims = dataSet.dims();
    raw.values.resize(dataSet.numValues());
    dataSet.readData(raw.values);
    return raw;
}

/*! \brief Create an index group holding raw `ranges` and `indices` data sets.
 *
 * The writer decomposes a group into disjoint parts and validates as it goes, so data sets
 * that the reader must nevertheless cope with can only be produced by writing them directly.
 *
 * \param[in] container HDF5 container to create the index group in.
 * \param[in] name      Name of the index group to create.
 * \param[in] rangeDims Dimensions of the `ranges` data set, unused when \p ranges is empty.
 * \param[in] ranges    Values of the `ranges` data set, or empty to omit the data set.
 * \param[in] indices   Values of the `indices` data set, or empty to omit the data set.
 */
void writeRawIndexGroup(const hid_t             container,
                        const char*             name,
                        const DataSetDims&      rangeDims,
                        const std::vector<int>& ranges,
                        const std::vector<int>& indices)
{
    const auto [child, childGuard] = makeH5mdGroupGuard(createGroup(container, name));
    if (!ranges.empty())
    {
        const H5mdFixedDataSet<int> rangeSet =
                H5mdDataSetBuilder<int>(child, "ranges").withDimension(rangeDims).build();
        rangeSet.writeData(ranges);
    }
    if (!indices.empty())
    {
        const H5mdFixedDataSet<int> indexSet =
                H5mdDataSetBuilder<int>(child, "indices").withDimension({ indices.size() }).build();
        indexSet.writeData(indices);
    }
}

/*! \brief Test reading and writing of multiple index groups at a time
 */
TEST_F(H5mdIndexGroupTest, ReadWriteMultipleIndexGroups)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    // Read index groups from a test file
    const std::string ndxfile = gmx::test::TestFileManager::getInputFilePath("alanine_vacuo.ndx");
    const std::vector<IndexGroup> indexGroups = init_index(ndxfile);

    {
        SCOPED_TRACE("Testing H5MD writing of index groups.");

        writeIndexGroups(group, indexGroups);

        // All written groups must be enumerated; getNamesOfIndexGroups returns alphabetical order
        std::vector<std::string> refGroupNames;
        refGroupNames.reserve(indexGroups.size());
        for (const auto& g : indexGroups)
        {
            refGroupNames.push_back(g.name);
        }
        std::sort(refGroupNames.begin(), refGroupNames.end());
        EXPECT_EQ(getNamesOfIndexGroups(group), refGroupNames);
    }

    {
        SCOPED_TRACE("Testing reading group index in batch (normal indices array).");

        const std::vector<IndexGroup> retIndexGroups = readIndexGroups(group);
        EXPECT_EQ(retIndexGroups.size(), indexGroups.size());
        for (const auto& g : indexGroups)
        {
            // Groups come back in alphabetical order, so match them up by name
            const auto it = std::find_if(retIndexGroups.begin(),
                                         retIndexGroups.end(),
                                         [&g](const IndexGroup& r) { return r.name == g.name; });
            ASSERT_NE(it, retIndexGroups.end()) << "Group not found: " << g.name;
            EXPECT_EQ(g.particleIndices.size(), it->particleIndices.size());
            EXPECT_EQ(g.particleIndices, it->particleIndices);
        }
    }

    {
        SCOPED_TRACE("Testing reading group index by its name (normal indices array).");

        // Read each group by its name and compare with the original one
        for (const auto& g : indexGroups)
        {
            IndexGroup retIndexGroup = readIndexGroup(group, g.name.c_str());
            EXPECT_EQ(retIndexGroup.name, g.name);
            EXPECT_EQ(retIndexGroup.particleIndices.size(), g.particleIndices.size());
            EXPECT_EQ(retIndexGroup.particleIndices, g.particleIndices);
        }
    }
}

/*! \brief One index group, and the data sets the writer is expected to store it in.
 *
 * The writer stores runs of three or more particles as rows of `ranges` and everything
 * else verbatim in `indices`, so which data sets exist is part of the contract and is
 * checked here alongside the round trip.
 */
struct EncodingCase
{
    //! Name of the index group, also used as the name of the test instance.
    const char* name;
    //! Particle indices of the group, in ascending order.
    std::vector<int> particleIndices;
    //! Whether the group is expected to hold an `indices` data set.
    bool hasIndices;
    //! Whether the group is expected to hold a `ranges` data set.
    bool hasRanges;
};

//! \brief Helper function for GTest to print encoding parameters.
void PrintTo(const EncodingCase& info, std::ostream* os)
{
    const auto toString = [](const int value) -> std::string { return std::to_string(value); };

    *os << info.name << " {" << formatAndJoin(info.particleIndices, ", ", toString) << "} -> "
        << (info.hasIndices ? "indices" : "no indices") << ", "
        << (info.hasRanges ? "ranges" : "no ranges");
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfEncodingCase(const ::testing::TestParamInfo<EncodingCase>& info)
{
    return info.param.name;
}

const EncodingCase g_encodingCases[] = {
    // A group with no particles still gets an indices data set, of length zero
    { "EmptyGroup", {}, true, false },
    // Runs shorter than three particles are not worth a range
    { "SingleElement", { 3 }, true, false },
    { "TwoElements", { 10, 11 }, true, false },
    // Runs of three or more particles are stored as ranges
    { "ThreeElements", { 4, 5, 6 }, false, true },
    { "NElements", { 4, 5, 6, 7, 8, 9, 10 }, false, true },
    // Particles that form no run at all are stored verbatim
    { "DiscreteElements", { 1, 21, 26, 37, 47 }, true, false },
    // A run may have a stride greater than one
    { "NonZeroStep", { 4, 8, 12, 16, 20, 24 }, false, true },
    // A group holding both runs and lone particles needs both data sets
    { "MixedSteps",
      { 1, 2, 3, 10, 20, 21, 22, 23, 44, 49, 54, 59, 77, 81, 82, 88, 90, 92, 95, 98, 101, 104 },
      true,
      true },
    { "MixedStep2", { 2, 5, 8, 9, 10, 15, 33, 44, 55, 105, 106, 107, 300, 400, 500 }, true, true },
};

//! \brief Parametrized test fixture for writing single index groups of different shapes.
class IndexGroupEncoding : public H5mdTestBase, public ::testing::WithParamInterface<EncodingCase>
{
};

TEST_P(IndexGroupEncoding, StoresExpectedDataSets)
{
    const EncodingCase& params     = GetParam();
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeIndexGroup(group, IndexGroup{ params.name, params.particleIndices });

    EXPECT_EQ(objectExists(group, formatString("%s/indices", params.name).c_str()), params.hasIndices);
    EXPECT_EQ(objectExists(group, formatString("%s/ranges", params.name).c_str()), params.hasRanges);
}

TEST_P(IndexGroupEncoding, RoundTripsUnchanged)
{
    const EncodingCase& params     = GetParam();
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeIndexGroup(group, IndexGroup{ params.name, params.particleIndices });

    const IndexGroup readBack = readIndexGroup(group, params.name);
    EXPECT_EQ(readBack.name, params.name);
    EXPECT_EQ(readBack.particleIndices, params.particleIndices);
}

TEST_P(IndexGroupEncoding, IsEnumerated)
{
    const EncodingCase& params     = GetParam();
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeIndexGroup(group, IndexGroup{ params.name, params.particleIndices });

    EXPECT_EQ(getNamesOfIndexGroups(group), std::vector<std::string>{ params.name })
            << "Every written group must be found when the container is enumerated.";
}

INSTANTIATE_TEST_SUITE_P(H5mdIndexGroupTest, IndexGroupEncoding, ::testing::ValuesIn(g_encodingCases), nameOfEncodingCase);

/*! \brief An index group in any order, as GROMACS allows, and the set it is stored as.
 *
 * The H5MD selection is the union of its indices and ranges, so order and repeats are not kept.
 */
struct NonAscendingCase
{
    //! Name of the index group, also used as the name of the test instance.
    const char* name;
    //! Particle indices of the group, as the caller gives them.
    std::vector<int> particleIndices;
    //! The indices read back: sorted, without repeats.
    std::vector<int> expectedSet;
};

//! \brief Helper function for GTest to print non-ascending group parameters.
void PrintTo(const NonAscendingCase& info, std::ostream* os)
{
    const auto toString = [](const int value) -> std::string { return std::to_string(value); };

    *os << info.name << " {" << formatAndJoin(info.particleIndices, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfNonAscendingCase(const ::testing::TestParamInfo<NonAscendingCase>& info)
{
    return info.param.name;
}

const NonAscendingCase g_nonAscendingCases[] = {
    { "DescendingPair", { 5, 4 }, { 4, 5 } },
    { "OneIndexOutOfPlace", { 3, 1, 2, 4, 5 }, { 1, 2, 3, 4, 5 } },
    { "ReversedRun", { 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 } },
    { "RunsInSwappedOrder", { 10, 11, 12, 13, 0, 1, 2, 3 }, { 0, 1, 2, 3, 10, 11, 12, 13 } },
    { "RepeatedIndex", { 3, 3, 3, 4, 5, 6 }, { 3, 4, 5, 6 } },
    { "RepeatedRun", { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 }, { 0, 1, 2, 3, 4 } },
    { "InterleavedStrides", { 0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11 }, { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 } },
    { "RandomOrder", { 42, 7, 19, 3, 100, 8, 9, 10, 55 }, { 3, 7, 8, 9, 10, 19, 42, 55, 100 } },
};

//! \brief Parametrized test fixture for writing index groups that are not strictly ascending.
class NonAscendingIndexGroup : public H5mdTestBase, public ::testing::WithParamInterface<NonAscendingCase>
{
};

TEST_P(NonAscendingIndexGroup, ReadsBackAsSortedSet)
{
    const NonAscendingCase& params = GetParam();
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeIndexGroup(group, IndexGroup{ params.name, params.particleIndices });

    EXPECT_EQ(readIndexGroup(group, params.name).particleIndices, params.expectedSet);
}

INSTANTIATE_TEST_SUITE_P(H5mdIndexGroupTest,
                         NonAscendingIndexGroup,
                         ::testing::ValuesIn(g_nonAscendingCases),
                         nameOfNonAscendingCase);

/*! \brief An index group holding no particles is stored as a zero-length `indices` data set.
 *
 * An empty group is a selection the user made, so it has to stay resolvable from the file.
 * Writing nothing at all would lose it, so the data set exists and holds no elements.
 */
TEST_F(H5mdIndexGroupTest, EmptyIndexGroupIsStoredAsZeroLengthIndices)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeIndexGroup(group, IndexGroup{ "empty_group", {} });

    ASSERT_TRUE(objectExists(group, "empty_group/indices"));
    EXPECT_EQ(readRawDataSet(group, "empty_group/indices").dims, DataSetDims({ 0 }));
    EXPECT_FALSE(objectExists(group, "empty_group/ranges"))
            << "A group with no particles has no runs to store as ranges.";
}

TEST_F(H5mdIndexGroupTest, SuccessAdditionWithoutDuplication)
{
    std::vector<IndexGroup> indexGroups;
    indexGroups.push_back(IndexGroup{ "group1", { 1, 2, 3 } });
    indexGroups.push_back(IndexGroup{ "group2", { 4, 5, 6, 7, 8 } });
    indexGroups.push_back(IndexGroup{ "group3", { 9, 10, 11, 12 } });

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));
    writeIndexGroups(group, indexGroups);

    std::vector<std::string> refGroupNames = { "group1", "group2", "group3" };
    EXPECT_EQ(getNamesOfIndexGroups(group), refGroupNames);

    {
        SCOPED_TRACE("Succeed upon writing index groups without duplication.");

        std::vector<IndexGroup> indexGroups2;
        indexGroups2.push_back(IndexGroup{ "group4", { 4, 5, 6 } });
        indexGroups2.push_back(IndexGroup{ "group5", { 7, 8, 9 } });
        writeIndexGroups(group, indexGroups2);

        // All groups written so far must be enumerated.
        const std::vector<std::string> expectedGroupNames = {
            "group1", "group2", "group3", "group4", "group5"
        };
        std::vector<std::string> sortedExpected = expectedGroupNames;
        std::sort(sortedExpected.begin(), sortedExpected.end());
        EXPECT_EQ(getNamesOfIndexGroups(group), sortedExpected);
    }
}

/*! \brief Test failure upon writing duplicated index groups
 */
TEST_F(H5mdIndexGroupTest, FailAdditionUponDuplicatedGroups)
{
    std::vector<IndexGroup> indexGroups;
    indexGroups.push_back(IndexGroup{ "group1", { 1, 2, 3 } });
    indexGroups.push_back(IndexGroup{ "group2", { 4, 5, 6, 7, 8 } });
    indexGroups.push_back(IndexGroup{ "group3", { 9, 10, 11, 12 } });

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));
    writeIndexGroups(group, indexGroups);

    std::vector<std::string> refGroupNames = { "group1", "group2", "group3" };
    EXPECT_EQ(getNamesOfIndexGroups(group), refGroupNames);

    {
        SCOPED_TRACE("Fail upon writing index groups with duplication.");

        std::vector<IndexGroup> indexGroupsWithDuplication;
        indexGroupsWithDuplication.push_back(IndexGroup{ "group1", { 1, 101, 1001, 10001, 100001 } });
        indexGroupsWithDuplication.push_back(IndexGroup{ "group2", { 3, 4, 5 } });

        // All should fail upon sanity check of duplicated group names.
        EXPECT_THROW(writeIndexGroup(group, indexGroupsWithDuplication[0]), FileIOError);
        EXPECT_THROW(writeIndexGroups(group, indexGroupsWithDuplication), FileIOError);

        // Failed writing should not add any group to the container.
        const std::vector<std::string> expectedGroupNames = { "group1", "group2", "group3" };
        std::vector<std::string>       sortedExpected     = expectedGroupNames;
        std::sort(sortedExpected.begin(), sortedExpected.end());
        EXPECT_EQ(getNamesOfIndexGroups(group), sortedExpected);
    }
}

//! \brief A name repeated within the input is rejected before any group is written.
TEST_F(H5mdIndexGroupTest, FailAdditionUponNameRepeatedInInput)
{
    const std::vector<IndexGroup> indexGroups = { IndexGroup{ "first", { 1, 2, 3 } },
                                                  IndexGroup{ "repeated", { 4, 5 } },
                                                  IndexGroup{ "repeated", { 6, 7 } } };

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    EXPECT_THROW(writeIndexGroups(group, indexGroups), FileIOError);
    EXPECT_TRUE(getNamesOfIndexGroups(group).empty())
            << "A rejected input must not leave its earlier groups behind in the container.";
}

/*! \brief One entry point of this module, called with the container handle it is given.
 *
 * Every entry point takes an HDF5 container as its first argument and has to reject an
 * invalid one rather than act on it. Listing them here makes that a checklist of the
 * module's public surface.
 */
struct EntryPoint
{
    //! Name of the entry point, also used as the name of the test instance.
    const char* name;
    //! Call the entry point with \p container, using arguments that are otherwise valid.
    std::function<void(hid_t container)> call;
};

//! \brief Helper function for GTest to print entry point parameters.
void PrintTo(const EntryPoint& info, std::ostream* os)
{
    *os << info.name;
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfEntryPoint(const ::testing::TestParamInfo<EntryPoint>& info)
{
    return info.param.name;
}

const EntryPoint g_entryPoints[] = {
    { "getNamesOfIndexGroups", [](const hid_t container) { getNamesOfIndexGroups(container); } },
    { "readIndexGroup", [](const hid_t container) { readIndexGroup(container, "group1"); } },
    { "readIndexGroups", [](const hid_t container) { readIndexGroups(container); } },
    { "writeIndexGroup",
      [](const hid_t container) { writeIndexGroup(container, IndexGroup{ "group1", { 1, 2, 3 } }); } },
    { "writeIndexGroups",
      [](const hid_t container)
      {
          const std::vector<IndexGroup> indexGroups     = { IndexGroup{ "group1", { 1, 2, 3 } } };
          writeIndexGroups(container, indexGroups);
      } },
};

//! \brief Parametrized test fixture for calling entry points with an invalid container.
class InvalidContainer : public H5mdTestBase, public ::testing::WithParamInterface<EntryPoint>
{
};

TEST_P(InvalidContainer, IsRejected)
{
    EXPECT_THROW(GetParam().call(H5I_INVALID_HID), FileIOError)
            << "An entry point must reject an invalid HDF5 handle rather than act on it.";
}

INSTANTIATE_TEST_SUITE_P(H5mdIndexGroupTest, InvalidContainer, ::testing::ValuesIn(g_entryPoints), nameOfEntryPoint);

/*! \brief An index group definition that the writer has to reject.
 *
 * The writer is the only guard against a file that the reader, or another tool reading the
 * specification, cannot make sense of, so a group it cannot represent faithfully must be
 * refused rather than written.
 */
struct RejectedIndexGroupCase
{
    //! Name of the case, used as the name of the test instance.
    const char* name;
    //! Why this group has to be rejected, reported when the expectation fails.
    const char* reason;
    //! The group to attempt to write.
    IndexGroup indexGroup;
};

//! \brief Helper function for GTest to print rejected group parameters.
void PrintTo(const RejectedIndexGroupCase& info, std::ostream* os)
{
    const auto toString = [](const int value) -> std::string { return std::to_string(value); };

    *os << info.name << " '" << info.indexGroup.name << "' {"
        << formatAndJoin(info.indexGroup.particleIndices, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfRejectedIndexGroupCase(const ::testing::TestParamInfo<RejectedIndexGroupCase>& info)
{
    return info.param.name;
}

const RejectedIndexGroupCase g_rejectedIndexGroupCases[] = {
    // A name that HDF5 reads as a path would silently create nested groups that the reader
    // cannot enumerate, and a group with no name cannot be referred to at all
    { "EmptyName", "An empty group name is not writable.", IndexGroup{ "", { 1, 2, 3 } } },
    { "NameWithSlash",
      "A group name containing the HDF5 path separator is not writable.",
      IndexGroup{ "with/slash", { 1, 2, 3 } } },
    // A particle index is an offset into the particle arrays, matching the reader's
    // rejection of a negative start in a range
    { "NegativeIndices",
      "A negative particle index is not a valid offset.",
      IndexGroup{ "negative", { -5, -3, 1 } } },
    { "NegativeRun",
      "A negative particle index is not a valid offset, including within a run.",
      IndexGroup{ "negative_run", { -3, -2, -1, 0 } } },
    { "NegativeAfterPositive",
      "A negative particle index is not a valid offset, wherever it appears.",
      IndexGroup{ "negative_last", { 5, 6, -1 } } },
};

//! \brief Parametrized test fixture for writing index groups that have to be rejected.
class RejectedIndexGroup : public H5mdTestBase, public ::testing::WithParamInterface<RejectedIndexGroupCase>
{
};

TEST_P(RejectedIndexGroup, IsRejected)
{
    const RejectedIndexGroupCase& params = GetParam();
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    EXPECT_THROW(writeIndexGroup(group, params.indexGroup), FileIOError) << params.reason;
    EXPECT_TRUE(getNamesOfIndexGroups(group).empty())
            << "A rejected group must not be left behind in the container.";
}

INSTANTIATE_TEST_SUITE_P(H5mdIndexGroupTest,
                         RejectedIndexGroup,
                         ::testing::ValuesIn(g_rejectedIndexGroupCases),
                         nameOfRejectedIndexGroupCase);

/*! \brief Writing nothing is a caller mistake rather than a no-op.
 *
 * Unlike the cases above this concerns the array, not any group in it, so there is nothing
 * to parametrize over.
 */
TEST_F(H5mdIndexGroupTest, WriteRejectsEmptyIndexGroupArray)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    const std::vector<IndexGroup> emptyIndexGroups;
    EXPECT_THROW(writeIndexGroups(group, emptyIndexGroups), FileIOError);
}

//! \brief A group that is not in the container cannot be read.
TEST_F(H5mdIndexGroupTest, ReadRejectsNonExistingIndexGroup)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    EXPECT_THROW(readIndexGroup(group, "non_existing_group"), FileIOError);
}

/*! \brief A group holding neither data set is not an index group, e.g. a string-based selection.
 *
 * The writer always gives a group one data set, so a raw group is the only way to cover this.
 */
TEST_F(H5mdIndexGroupTest, ReadRejectsGroupWithoutIndexData)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeRawIndexGroup(group, "string_selection", {}, {}, {});

    EXPECT_THROW(readIndexGroup(group, "string_selection"), FileIOError);
}

/*! \brief A negative index in the file is rejected, as a negative range start is.
 *
 * The writer refuses to produce this, so a raw data set is the only way to cover the reader.
 */
TEST_F(H5mdIndexGroupTest, ReadRejectsNegativeIndicesDataSet)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeRawIndexGroup(group, "negative_indices", {}, {}, { -1, 0, 1 });

    EXPECT_THROW(readIndexGroup(group, "negative_indices"), FileIOError);
}

/*! \brief Verify the on-disk encoding of index groups against the H5MD specification.
 *
 * Every other value check in this file is a write->read round trip, which passes for any
 * self-consistent encoding: the conversion applied on write is undone on read.
 */
TEST_F(H5mdIndexGroupTest, StoredEncodingFollowsSpecification)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    {
        SCOPED_TRACE("A consecutive run starting at particle 0 is stored 0-based.");
        writeIndexGroup(group, IndexGroup{ "from_zero", { 0, 1, 2 } });
        const auto [child, childGuard] = makeH5mdGroupGuard(openGroup(group, "from_zero"));
        const RawDataSet ranges        = readRawDataSet(child, "ranges");
        // 1-based storage would give [1, 4, 1]; an inclusive stop would give [0, 2, 1]
        EXPECT_EQ(ranges.dims, DataSetDims({ 1, 3 }));
        EXPECT_EQ(ranges.values, (std::vector<int>{ 0, 3, 1 }));
        EXPECT_FALSE(objectExists(child, "indices"))
                << "A purely consecutive group needs no plain indices data set.";
    }

    {
        SCOPED_TRACE("A strided run stores an exclusive stop one step past the last particle.");
        writeIndexGroup(group, IndexGroup{ "strided", { 2, 4, 6, 8, 10 } });
        const auto [child, childGuard] = makeH5mdGroupGuard(openGroup(group, "strided"));
        const RawDataSet ranges        = readRawDataSet(child, "ranges");
        // An inclusive stop would give [2, 10, 2]
        EXPECT_EQ(ranges.dims, DataSetDims({ 1, 3 }));
        EXPECT_EQ(ranges.values, (std::vector<int>{ 2, 12, 2 }));
    }

    {
        SCOPED_TRACE("Irregular indices are stored verbatim and 0-based.");
        writeIndexGroup(group, IndexGroup{ "irregular", { 1, 21, 26 } });
        const auto [child, childGuard] = makeH5mdGroupGuard(openGroup(group, "irregular"));
        const RawDataSet indices       = readRawDataSet(child, "indices");
        // 1-based storage would give [2, 22, 27]
        EXPECT_EQ(indices.dims, DataSetDims({ 3 }));
        EXPECT_EQ(indices.values, (std::vector<int>{ 1, 21, 26 }));
        EXPECT_FALSE(objectExists(child, "ranges"))
                << "A group with no runs of three or more needs no ranges data set.";
    }

    {
        SCOPED_TRACE("Several runs are stored as separate rows of the [N, 3] ranges data set.");
        writeIndexGroup(group, IndexGroup{ "two_runs", { 0, 1, 2, 10, 20, 30 } });
        const auto [child, childGuard] = makeH5mdGroupGuard(openGroup(group, "two_runs"));
        const RawDataSet ranges        = readRawDataSet(child, "ranges");
        EXPECT_EQ(ranges.dims, DataSetDims({ 2, 3 }));
        EXPECT_EQ(ranges.values, (std::vector<int>{ 0, 3, 1, 10, 40, 10 }));
    }

    {
        SCOPED_TRACE("A mixed group splits across both data sets without losing particles.");
        writeIndexGroup(group, IndexGroup{ "mixed", { 0, 1, 2, 10 } });
        const auto [child, childGuard] = makeH5mdGroupGuard(openGroup(group, "mixed"));
        const RawDataSet ranges        = readRawDataSet(child, "ranges");
        const RawDataSet indices       = readRawDataSet(child, "indices");
        EXPECT_EQ(ranges.dims, DataSetDims({ 1, 3 }));
        EXPECT_EQ(ranges.values, (std::vector<int>{ 0, 3, 1 }));
        EXPECT_EQ(indices.dims, DataSetDims({ 1 }));
        EXPECT_EQ(indices.values, (std::vector<int>{ 10 }));
    }
}

/*! \brief A `ranges` data set whose layout or values cannot be trusted.
 *
 * The reader parses files this process did not write, so a `ranges` data set that is not
 * [N, 3], or whose rows describe an impossible range, must be rejected rather than indexed
 * row-wise. The writer cannot produce any of these, so each is written as a raw data set.
 */
struct MalformedRangesCase
{
    //! Name of the index group, also used as the name of the test instance.
    const char* name;
    //! Why this data set has to be rejected, reported when the expectation fails.
    const char* reason;
    //! Dimensions of the raw `ranges` data set.
    DataSetDims dims;
    //! Values of the raw `ranges` data set.
    std::vector<int> values;
};

//! \brief Helper function for GTest to print malformed range parameters.
void PrintTo(const MalformedRangesCase& info, std::ostream* os)
{
    const auto toString = [](const int value) -> std::string { return std::to_string(value); };

    *os << info.name << " {" << formatAndJoin(info.values, ", ", toString) << "}";
}

//! \brief Helper function for GTest to construct test names.
std::string nameOfMalformedRangesCase(const ::testing::TestParamInfo<MalformedRangesCase>& info)
{
    return info.param.name;
}

const MalformedRangesCase g_malformedRangesCases[] = {
    { "RankOne", "A rank-1 ranges data set is not [N, 3].", { 3 }, { 0, 3, 1 } },
    { "WrongRowWidth", "A ranges data set with rows of four is not [N, 3].", { 1, 4 }, { 0, 3, 1, 0 } },
    { "ZeroStride", "A zero stride would never terminate the expansion.", { 1, 3 }, { 0, 3, 0 } },
    { "NegativeStride", "A negative stride would never terminate the expansion.", { 1, 3 }, { 0, 3, -1 } },
    { "NegativeStart", "A negative start is not a valid particle index.", { 1, 3 }, { -1, 3, 1 } },
    // span + step - 1 leaves the range of an int
    { "StopOverflowsCount",
      "A stop too close to the largest index would overflow while counting.",
      { 1, 3 },
      { 0, INT_MAX, 2 } },
    // (stop - 1) + step leaves the range of an int, wrapping below stop so that the
    // expansion would never reach its bound
    { "StopOverflowsExpansion",
      "A stop too close to the largest index would overflow while expanding.",
      { 1, 3 },
      { 2, INT_MAX, INT_MAX } },
};

//! \brief Parametrized test fixture for reading `ranges` data sets that must be rejected.
class MalformedRanges : public H5mdTestBase, public ::testing::WithParamInterface<MalformedRangesCase>
{
};

TEST_P(MalformedRanges, IsRejected)
{
    const MalformedRangesCase& params = GetParam();
    const auto [group, groupGuard]    = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeRawIndexGroup(group, params.name, params.dims, params.values, {});

    EXPECT_THROW(readIndexGroup(group, params.name), FileIOError) << params.reason;
}

INSTANTIATE_TEST_SUITE_P(H5mdIndexGroupTest,
                         MalformedRanges,
                         ::testing::ValuesIn(g_malformedRangesCases),
                         nameOfMalformedRangesCase);

/*! \brief A range whose stop does not lie above its start describes no particles.
 *
 * This is well-formed rather than malformed, so it is read rather than rejected. Read as
 * an unsigned bound a negative stop would instead wrap and expand for effectively ever.
 */
TEST_F(H5mdIndexGroupTest, RangeWithNonPositiveSpanYieldsNoParticles)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    writeRawIndexGroup(group, "negative_stop", { 1, 3 }, { 0, -1, 1 }, {});

    EXPECT_TRUE(readIndexGroup(group, "negative_stop").particleIndices.empty());
}

/*! \brief An index group is the union of its `ranges` and `indices` data sets.
 *
 * The specification defines the group as the union of both data sets, so an index the file
 * lists more than once must be reported once. The writer decomposes into disjoint parts and
 * so cannot produce an overlap, which makes raw data sets the only way to cover this.
 */
TEST_F(H5mdIndexGroupTest, OverlappingRangesAndIndicesAreMergedAsUnion)
{
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(fileid(), "index_groups"));

    {
        SCOPED_TRACE("An index listed in both data sets is reported once.");
        // ranges [0, 5, 1) -> {0, 1, 2, 3, 4}, indices -> {2, 4, 7}
        writeRawIndexGroup(group, "both_overlap", { 1, 3 }, { 0, 5, 1 }, { 2, 4, 7 });
        EXPECT_EQ(readIndexGroup(group, "both_overlap").particleIndices,
                  (std::vector<int>{ 0, 1, 2, 3, 4, 7 }));
    }

    {
        SCOPED_TRACE("Overlapping ranges are merged into a single ascending run.");
        // [0, 4, 1) -> {0, 1, 2, 3} and [2, 6, 1) -> {2, 3, 4, 5}
        writeRawIndexGroup(group, "overlapping_ranges", { 2, 3 }, { 0, 4, 1, 2, 6, 1 }, {});
        EXPECT_EQ(readIndexGroup(group, "overlapping_ranges").particleIndices,
                  (std::vector<int>{ 0, 1, 2, 3, 4, 5 }));
    }

    {
        SCOPED_TRACE("A repeated plain index is reported once.");
        writeRawIndexGroup(group, "repeated_indices", {}, {}, { 3, 3, 8 });
        EXPECT_EQ(readIndexGroup(group, "repeated_indices").particleIndices, (std::vector<int>{ 3, 8 }));
    }
}

} // namespace
} // namespace test
} // namespace gmx
