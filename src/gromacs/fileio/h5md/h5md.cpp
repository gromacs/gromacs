/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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

/*! \brief I/o interface to H5MD HDF5 files.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#include "gmxpre.h"

#include "h5md.h"

#include "config.h"

#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "h5md_attribute.h"
#    include "h5md_error.h"
#    include "h5md_framedatasetbuilder.h"
#    include "h5md_group.h"
#    include "h5md_guard.h"
#    include "h5md_particleblock.h"
#    include "h5md_timedatablock.h"
#    include "h5md_util.h"
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")
#else
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")
#endif // GMX_USE_HDF5

namespace gmx
{

// Only declare these variables if they will be used (i.e. if we're using HDF5)
#if GMX_USE_HDF5
//! \brief Path to particle block group from the file root.
constexpr char c_particlesGroupPath[] = "/particles";
//! \brief Name of position group inside the particle block in the H5md specification.
constexpr char c_positionGroupName[] = "position";
//! \brief Name of velocity group inside the particle block in the H5md specification.
constexpr char c_velocityGroupName[] = "velocity";
//! \brief Name of force group inside the particle block in the H5md specification.
constexpr char c_forceGroupName[] = "force";
//! \brief Name of box group inside the particle block in the H5md specification.
constexpr char c_boxGroupName[] = "box";
//! \brief Group name for all atoms in system.
constexpr char c_fullSystemGroupName[] = "system";
//! \brief Name of group inside the simulation box groups which contains the box size data.
constexpr char c_boxSizeName[] = "edges";
//! \brief Attribute name for number of dimensions of simulation box.
constexpr char c_boxDimensionAttributeKey[] = "dimension";
//! \brief Attribute name for periodic boundary definition of simulation box.
constexpr char c_boxBoundaryAttributeKey[] = "boundary";
//! \brief Name of the value data set in H5mdTimeDataBlock groups.
constexpr char c_valueName[] = "value";
//! \brief Name of the step data set in H5mdTimeDataBlock groups.
constexpr char c_stepName[] = "step";
//! \brief Name of the time data set in H5mdTimeDataBlock groups.
constexpr char c_timeName[] = "time";
//! \brief Version of H5MD specification used for this file ({ majorVersion, minorVersion }).
const std::vector<int> c_h5mdSpecificationVersion = { 1, 1 };
//! \brief Path to H5MD metadata group from the file root.
constexpr char c_h5mdMetaDataGroupName[] = "/h5md";
//! \brief Name of author group in the H5MD metadata group.
constexpr char c_h5mdAuthorGroupName[] = "author";
//! \brief Name of creator (program) group in the H5MD metadata group.
constexpr char c_h5mdCreatorGroupName[] = "creator";
//! \brief Attribute name for the author or creator (program) name.
constexpr char c_h5mdNameAttributeKey[] = "name";
//! \brief Attribute name for a version specification.
constexpr char c_h5mdVersionAttributeKey[] = "version";
//! \brief Maximum length of author names (used to allocate memory).
constexpr int c_maxUserNameLength = 4096;
#endif

H5md::H5md(const std::filesystem::path& fileName, const H5mdFileMode mode)
{
#if GMX_USE_HDF5
    /* Disable automatic HDF5 error output, e.g. when items are not found. Explicit H5EPrint2() will
     * still print error messages. */
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

    switch (mode)
    {
        case H5mdFileMode::Write:
            file_ = H5Fcreate(
                    fileName.string().c_str(), H5F_ACC_TRUNC, H5Pcreate(H5P_FILE_CREATE), H5P_DEFAULT);
            break;
        case H5mdFileMode::Read:
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
        default: throw NotImplementedError("Appending to H5MD is not implemented yet.");
    }
    gmx::throwUponInvalidHid(file_, "Cannot open H5MD file.");

    filemode_ = mode;

#else
    GMX_UNUSED_VALUE(fileName);
    GMX_UNUSED_VALUE(mode);
    throw NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

H5md::~H5md()
{
#if GMX_USE_HDF5
    if (handleIsValid(file_))
    {
        // Do not throw exceptions when flushing from the destructor.
        flush(false);
        H5Fclose(file_);
    }

    /* Do not throw, if GMX_USE_HDF5 is false, in the destructor. */

#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
hid_t H5md::fileid() const
{
#if GMX_USE_HDF5
    return file_;
#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::flush(bool throwExceptionUponError)
{
#if GMX_USE_HDF5
    if (throwExceptionUponError)
    {
        GMX_ASSERT(handleIsValid(file_), "Cannot flush an invalid H5MD file.");
        gmx::throwUponH5mdError(H5Fflush(file_, H5F_SCOPE_LOCAL) < 0, "Error flushing H5MD.");
    }
    else
    {
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }

#else
    GMX_UNUSED_VALUE(throwExceptionUponError);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setAuthor(const std::string& authorName)
{
#if GMX_USE_HDF5
    const auto [authorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/author"));
    setAttribute(authorGroup, "name", authorName);
#else
    GMX_UNUSED_VALUE(authorName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::author()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/author"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/author"));
        return getAttribute<std::string>(group, "name");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setCreatorProgramName(const std::string& creatorName)
{
#if GMX_USE_HDF5
    const auto [creatorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/creator"));
    setAttribute(creatorGroup, "name", creatorName);
#else
    GMX_UNUSED_VALUE(creatorName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::creatorProgramName()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/creator"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/creator"));
        return getAttribute<std::string>(group, "name");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setCreatorProgramVersion(const std::string& version)
{
#if GMX_USE_HDF5
    const auto [creatorGroup, groupGuard] =
            makeH5mdGroupGuard(openOrCreateGroup(file_, "h5md/creator"));
    setAttribute(creatorGroup, "version", version);
#else
    GMX_UNUSED_VALUE(version);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
std::optional<std::string> H5md::creatorProgramVersion()
{
#if GMX_USE_HDF5
    if (objectExists(file_, "h5md/creator"))
    {
        const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, "h5md/creator"));
        return getAttribute<std::string>(group, "version");
    }
    else
    {
        return std::nullopt;
    }

#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

#if GMX_USE_HDF5
/*! \brief Create and link the simulation box data sets inside \p selectionGroup.
 *
 * Per the H5md specification, if position data is stored we also store the box size
 * of the system. This data is stored in a subgroup similar to a `H5mdTimeDataBlock`
 * but with its step and time data sets hard linked to those of the position data set.
 *
 * This function creates a new data set for storing a box matrix (type: real[DIM][DIM])
 * and assigns it to the simulation box container in \p particleBlock. It then sets up
 * links to the step and time data sets of the position data set in \p selectionGroup.
 *
 * \warning Must be called after the position group has been set up.
 *
 * \param[in]  selectionGroup Handle to selection group inside /particles.
 * \param[out] blockBuilder   Builder to assign the simulation box data set to.
 *
 * \throws gmx::FileIOError if there was an error setting up the links to the position
 * step or time data sets.
 */
static void setupSimulationBoxDataSet(const hid_t selectionGroup, H5mdParticleBlockBuilder& blockBuilder)
{
    // Per the H5md spec we should write the box only once if it is constant.
    // At the time of writing we don't have access to the box values when the H5md
    // object is setup, so we'll return to this later and for now only set up
    // to write one box per frame.
    // const bool boxIsConstant = inputRecord.pressureCouplingOptions.epc == PressureCoupling::No;

    // Additionally, per the H5md spec we should write cubic boxes as RVecs
    // and non-cubic as matrix. Can we infer this at this point? For now we always
    // write the full matrix.
    // const bool boxIsCubic = (box[XX][YY] == 0.0) && (box[XX][ZZ] == 0.0) && (box[YY][XX] == 0.0)
    //                         && (box[YY][ZZ] == 0.0) && (box[ZZ][XX] == 0.0) && (box[ZZ][YY] == 0.0);
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(openGroup(selectionGroup, c_boxGroupName));
    const auto [edgesGroup, edgesGroupGuard] = makeH5mdGroupGuard(createGroup(boxGroup, c_boxSizeName));

    // Matrices are stored as real[DIM][DIM], which is the data set frame dimension
    blockBuilder.setBox(
            H5mdFrameDataSetBuilder<real>(edgesGroup, c_valueName).withFrameDimension({ DIM, DIM }).build());

    const auto [positionGroup, positionGroupGuard] =
            makeH5mdGroupGuard(openGroup(selectionGroup, c_positionGroupName));
    throwUponH5mdError(
            H5Lcreate_hard(positionGroup, c_stepName, edgesGroup, c_stepName, H5P_DEFAULT, H5P_DEFAULT) < 0,
            "Could not create hard link from position/step to box/edges/step");
    throwUponH5mdError(
            H5Lcreate_hard(positionGroup, c_timeName, edgesGroup, c_timeName, H5P_DEFAULT, H5P_DEFAULT) < 0,
            "Could not create hard link from position/time to box/edges/time");
}

/*! \brief Create the simulation box group inside \p selectionGroup.
 *
 * The simulation box attributes "dimension" and "boundary" are written to the box group.
 *
 * \warning Must be called after the position group has been set up.
 *
 * \param[in] selectionGroup Handle to group in which to create the box group.
 * \param[in] inputRecord    Simulation input record.
 */
static void setupSimulationBoxGroup(const hid_t selectionGroup, const t_inputrec& inputRecord)
{
    const auto [boxGroup, boxGroupGuard] =
            makeH5mdGroupGuard(createGroup(selectionGroup, c_boxGroupName));

    const std::vector<std::string> boundary = [&]() -> std::vector<std::string>
    {
        constexpr char periodicValue[] = "periodic";
        constexpr char noneValue[]     = "none";
        switch (inputRecord.pbcType)
        {
            case PbcType::Xyz: return { periodicValue, periodicValue, periodicValue };
            case PbcType::XY: return { periodicValue, periodicValue, noneValue };
            default: return { noneValue, noneValue, noneValue };
        }
    }();

    std::vector<char>                 outputBuffer;
    const ArrayRef<const std::string> boundaryRef = makeConstArrayRef(boundary);
    setAttributeStringVector(boxGroup,
                             c_boxBoundaryAttributeKey,
                             std::move(outputBuffer),
                             boundaryRef.begin(),
                             boundaryRef.end());
    setAttribute(boxGroup, c_boxDimensionAttributeKey, static_cast<int32_t>(boundary.size()));
}
#endif

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupParticleBlockForGroup(const gmx_mtop_t&   topology,
                                      ArrayRef<const int> selectionIndices,
                                      const std::string&  selectionName,
                                      const t_inputrec&   inputRecord)
{
#if GMX_USE_HDF5
    const auto [particlesGroup, particlesGroupGuard] =
            makeH5mdGroupGuard(createGroup(file_, c_particlesGroupPath));
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(createGroup(particlesGroup, selectionName.c_str()));
    setupSimulationBoxGroup(selectionGroup, inputRecord);

    const hsize_t numAtoms = selectionIndices.empty() ? topology.natoms : selectionIndices.size();
    throwUponH5mdError(
            numAtoms == 0,
            formatString("Cannot setup particle group '%s': no atoms in group", selectionName.c_str()));

    const DataSetDims frameDims = numAtoms > 0 ? DataSetDims{ numAtoms } : DataSetDims{};

    H5mdParticleBlockBuilder blockBuilder;
    if (inputRecord.nstxout > 0)
    {
        blockBuilder.setPosition(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_positionGroupName)
                                         .withFrameDimension(frameDims)
                                         .build());
        setupSimulationBoxDataSet(selectionGroup, blockBuilder);
    }
    if (inputRecord.nstvout > 0)
    {
        blockBuilder.setVelocity(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_velocityGroupName)
                                         .withFrameDimension(frameDims)
                                         .build());
    }
    if (inputRecord.nstfout > 0)
    {
        blockBuilder.setForce(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_forceGroupName)
                                      .withFrameDimension(frameDims)
                                      .build());
    }
    particleBlocks_.insert({ selectionName, TrajectoryReadCursor(blockBuilder.build()) });
#else
    GMX_UNUSED_VALUE(topology);
    GMX_UNUSED_VALUE(selectionIndices);
    GMX_UNUSED_VALUE(selectionName);
    GMX_UNUSED_VALUE(inputRecord);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupMetadataGroup()
{
#if GMX_USE_HDF5
    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(file_, c_h5mdMetaDataGroupName));
    setAttributeVector<int>(group, c_h5mdVersionAttributeKey, c_h5mdSpecificationVersion);

    const auto [authorGroup, authorGroupGuard] =
            makeH5mdGroupGuard(createGroup(group, c_h5mdAuthorGroupName));
    std::string username(c_maxUserNameLength, '\0');
    if (gmx_getusername(username.data(), c_maxUserNameLength) != 0)
    {
        username = "<unknown author>";
    }
    setAttribute(authorGroup, c_h5mdNameAttributeKey, username);

    const auto [creatorGroup, creatorGroupGuard] =
            makeH5mdGroupGuard(createGroup(group, c_h5mdCreatorGroupName));
    setAttribute(creatorGroup, c_h5mdNameAttributeKey, "GROMACS");
    setAttribute(creatorGroup, c_h5mdVersionAttributeKey, gmx_version());
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFileFromInput(const gmx_mtop_t& topology, const t_inputrec& inputRecord)
{
#if GMX_USE_HDF5
    setupMetadataGroup();
    setupParticleBlockForGroup(topology, {}, c_fullSystemGroupName, inputRecord);
#else
    GMX_UNUSED_VALUE(topology);
    GMX_UNUSED_VALUE(inputRecord);
#endif
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
