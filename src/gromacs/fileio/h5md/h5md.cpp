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
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 */

#include "gmxpre.h"

#include "h5md.h"

#include "config.h"

#include <filesystem>
#include <optional>
#include <string>

#include "gromacs/fileio/h5md/h5md_wrapper.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/selection/selection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
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
#    include "h5md_topologyutils.h"
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
//! \brief Name of modules group in the H5MD metadata group.
constexpr char c_h5mdModulesGroupName[] = "modules";
//! \brief Attribute name for the author or creator (program) name.
constexpr char c_h5mdNameAttributeKey[] = "name";
//! \brief Attribute name for a version specification.
constexpr char c_h5mdVersionAttributeKey[] = "version";
//! \brief Path to H5MD connectivity group from the file root.
constexpr char c_h5mdConnectivityGroupPath[] = "/connectivity";
//! \brief Maximum length of author names (used to allocate memory).
constexpr int c_maxUserNameLength = 4096;
//! \brief Unit for position data for simulations produced by mdrun.
constexpr char c_positionUnit[] = "nm";
//! \brief Unit for velocity data for simulations produced by mdrun.
constexpr char c_velocityUnit[] = "nm ps-1";
//! \brief Unit for force data for simulations produced by mdrun.
constexpr char c_forceUnit[] = "kJ mol-1 nm-1";

namespace
{
/*! \brief Name of attribute which stores a URL for a module.
 *
 * Storing a URL value is not a part of the H5MD module specification
 * but we can use this to add extra information for users by linking
 * them to specifications or just our own website.
 *
 * This does not break compatibility since setting this field
 * is fully independent of other data stored in a module.
 */
constexpr char c_urlAttributeKey[] = "url";

/*! \brief Set up the GROMACS module with metadata in \p modulesGroup.
 *
 * Creates a module group with name "gromacs" and writes version
 * metadata and a project URL to it.
 */
static void setupGromacsModule(const hid_t modulesGroup)
{
    constexpr char name[] = "gromacs";
    // The version attribute is required by the H5MD specification to follow
    // semantic versioning and to store two values: { majorVersion, minorVersion }
    //
    // To follow semantic versioning any change that breaks backwards compatibility
    // for data stored inside or defined by this module must also bump the major
    // version number, similar to how the TPR versioning works.
    //
    // The exception to the above is for major version = 0, which is treated as
    // experimental by semantic versioning rules. This module is marked as experimental
    // until further notice.
    constexpr int  version[] = { 0, 1 };
    constexpr char url[]     = "https://www.gromacs.org";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(modulesGroup, name));
    setAttributeVector(group, c_h5mdVersionAttributeKey, constArrayRefFromArray(version, 2));
    setAttribute(group, c_urlAttributeKey, url);
}

/*! \brief Set up the units module with metadata in \p modulesGroup.
 *
 * Creates a module group with name "units" in and writes required
 * metadata and a URL to its specification to it. Per the specification
 * a unit system attribute is also written to the created group.
 */
static void setupUnitsModule(const hid_t modulesGroup, const char* unitSystem = "SI")
{
    // Module specification per the official H5MD specification
    constexpr char name[]    = "units";
    constexpr int  version[] = { 1, 0 };
    constexpr char url[]     = "https://h5md.nongnu.org/modules/units.html";

    const auto [group, groupGuard] = makeH5mdGroupGuard(createGroup(modulesGroup, name));
    setAttributeVector(group, c_h5mdVersionAttributeKey, constArrayRefFromArray(version, 2));
    setAttribute(group, "system", unitSystem);
    setAttribute(group, c_urlAttributeKey, url);
}
} // namespace
#endif // GMX_USE_HDF5

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
    GMX_H5MD_THROW_UPON_INVALID_HID(file_, "Cannot open H5MD file.");

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
        GMX_H5MD_THROW_UPON_ERROR(H5Fflush(file_, H5F_SCOPE_LOCAL) < 0, "Error flushing H5MD.");
    }
    else
    {
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }

#else
    GMX_UNUSED_VALUE(throwExceptionUponError);
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
    // TODO: Per the H5md spec we should write the box only once if it is constant.
    // At the time of writing we don't have access to the box values when the H5md
    // object is setup, so we'll return to this later and for now only set up
    // to write one box per frame.
    // const bool boxIsConstant = inputRecord.pressureCouplingOptions.epc == PressureCoupling::No;

    // TODO: Additionally, per the H5md spec we should write rectangular boxes as RVecs
    // and triclinic boxes as matrix. Can we infer this at this point? For now we always
    // write the full matrix.
    // const bool boxIsCubic = (box[XX][YY] == 0.0) && (box[XX][ZZ] == 0.0) && (box[YY][XX] == 0.0)
    //                         && (box[YY][ZZ] == 0.0) && (box[ZZ][XX] == 0.0) && (box[ZZ][YY] == 0.0);
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(openGroup(selectionGroup, c_boxGroupName));
    const auto [edgesGroup, edgesGroupGuard] = makeH5mdGroupGuard(createGroup(boxGroup, c_boxSizeName));

    // Matrices are stored as real[DIM][DIM], which is the data set frame dimension
    blockBuilder.setBox(H5mdFrameDataSetBuilder<real>(edgesGroup, c_valueName)
                                .withFrameDimension({ DIM, DIM })
                                .withUnit(c_positionUnit)
                                .build());

    const auto [positionGroup, positionGroupGuard] =
            makeH5mdGroupGuard(openGroup(selectionGroup, c_positionGroupName));
    GMX_H5MD_THROW_UPON_ERROR(
            H5Lcreate_hard(positionGroup, c_stepName, edgesGroup, c_stepName, H5P_DEFAULT, H5P_DEFAULT) < 0,
            "Could not create hard link from position/step to box/edges/step");
    GMX_H5MD_THROW_UPON_ERROR(
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

    setAttributeVector(boxGroup, c_boxBoundaryAttributeKey, boundary);
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
    GMX_H5MD_THROW_UPON_ERROR(
            numAtoms == 0,
            formatString("Cannot setup particle group '%s': no atoms in group", selectionName.c_str()));

    const DataSetDims frameDims = numAtoms > 0 ? DataSetDims{ numAtoms } : DataSetDims{};

    H5mdParticleBlockBuilder blockBuilder;
    if (inputRecord.nstxout > 0)
    {
        blockBuilder.setPosition(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_positionGroupName)
                                         .withFrameDimension(frameDims)
                                         .withUnit(c_positionUnit)
                                         .build());
        setupSimulationBoxDataSet(selectionGroup, blockBuilder);
    }
    if (inputRecord.nstvout > 0)
    {
        blockBuilder.setVelocity(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_velocityGroupName)
                                         .withFrameDimension(frameDims)
                                         .withUnit(c_velocityUnit)
                                         .build());
    }
    if (inputRecord.nstfout > 0)
    {
        blockBuilder.setForce(H5mdTimeDataBlockBuilder<RVec>(selectionGroup, c_forceGroupName)
                                      .withFrameDimension(frameDims)
                                      .withUnit(c_forceUnit)
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
void H5md::setupGromacsTopology(const gmx_mtop_t& topology)
{
#if GMX_USE_HDF5
    if (topology.moltype.empty())
    {
        return;
    }
    // Create the HDF5 group for storing the GROMACS topology
    constexpr char    c_moduleNameGromacsTopology[] = "gromacs_topology";
    const std::string pathToGromacsTopology         = joinStrings(
            { c_h5mdMetaDataGroupName, c_h5mdModulesGroupName, c_moduleNameGromacsTopology }, "/");
    const auto [gmxTop, gmxTopGuard] =
            makeH5mdGroupGuard(createGroup(file_, pathToGromacsTopology.c_str()));

    writeMoleculeTypes(gmxTop, makeConstArrayRef(topology.moltype));
    writeMoleculeBlocks(gmxTop, makeConstArrayRef(topology.molblock));
    // The h5md_topologyutils module takes care of the versioning of the GROMACS topology
    labelInternalTopologyVersion(gmxTop);
    labelTopologyName(gmxTop, *(topology.name));
#else
    GMX_UNUSED_VALUE(topology);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupBondConnectivity(const gmx_mtop_t& topology)
{
#if GMX_USE_HDF5
    if (topology.moltype.empty())
    {
        return;
    }
    const auto [connectivityGroup, connectivityGroupGuard] =
            makeH5mdGroupGuard(createGroup(file_, c_h5mdConnectivityGroupPath));

    writeBonds(topology, connectivityGroup);
    // TODO: need to determine if there is a disulfide bond in the system
    // writeDisulfideBonds(topology, connectivityGroup);
#else
    GMX_UNUSED_VALUE(topology);
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

    const auto [modulesGroup, modulesGroupGuard] =
            makeH5mdGroupGuard(createGroup(group, c_h5mdModulesGroupName));
    setupGromacsModule(modulesGroup);
    setupUnitsModule(modulesGroup);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFileFromInput(const gmx_mtop_t& topology, const t_inputrec& inputRecord)
{
#if GMX_USE_HDF5
    setupMetadataGroup();
    setupParticleBlockForGroup(topology, {}, c_fullSystemGroupName, inputRecord);
    setupGromacsTopology(topology);
    setupBondConnectivity(topology);
#else
    GMX_UNUSED_VALUE(topology);
    GMX_UNUSED_VALUE(inputRecord);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupParticleBlockForGroupFromExistingFile(const std::string& selectionName)
{
#if GMX_USE_HDF5
    const auto [particlesGroup, particlesGroupGuard] =
            makeH5mdGroupGuard(openGroup(file_, c_particlesGroupPath));
    GMX_H5MD_THROW_UPON_INVALID_HID(particlesGroup, "H5md trajectory file has no /particles group");
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(openGroup(particlesGroup, selectionName.c_str()));
    GMX_H5MD_THROW_UPON_INVALID_HID(selectionGroup, "No trajectory data found for system group");

    H5mdParticleBlockBuilder blockBuilder;

    // We don't have access to t_inputrec from the trajectory file opening framework,
    // so open all available data sets (opening handles only is cheap).
    if (objectExists(selectionGroup, c_positionGroupName))
    {
        blockBuilder.setPosition(H5mdTimeDataBlock<RVec>(selectionGroup, c_positionGroupName));

        // Right now we always write the box size for every position frame
        // and thus must always find a corresponding data set here.
        // TODO: Revisit as needed when we've decided how to treat this when writing.

        // Construct the box size group name relative to the current selection group,
        // then open the box size value data set inside
        const std::string boxSizeGroupPath = joinStrings({ c_boxGroupName, c_boxSizeName }, "/");
        const auto [boxSizeGroup, boxSizeGroupGuard] =
                makeH5mdGroupGuard(openGroup(selectionGroup, boxSizeGroupPath.c_str()));
        blockBuilder.setBox(H5mdFrameDataSet<real>(boxSizeGroup, c_valueName));
    }
    if (objectExists(selectionGroup, c_velocityGroupName))
    {
        blockBuilder.setVelocity(H5mdTimeDataBlock<RVec>(selectionGroup, c_velocityGroupName));
    }
    if (objectExists(selectionGroup, c_forceGroupName))
    {
        blockBuilder.setForce(H5mdTimeDataBlock<RVec>(selectionGroup, c_forceGroupName));
    }
    particleBlocks_.insert({ selectionName, TrajectoryReadCursor(blockBuilder.build()) });
#else
    GMX_UNUSED_VALUE(selectionName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFromExistingFile()
{
#if GMX_USE_HDF5
    setupParticleBlockForGroupFromExistingFile(c_fullSystemGroupName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool H5md::readNextFrame(t_trxframe* frame, const std::string& selectionName)
{
#if GMX_USE_HDF5
    bool frameWasRead = false;

    // TODO: We should check this in the file (and use appropriate conversion when reading)
    //
    // The current data set opening framework does not support opening double-precision `real`
    // data as single-precision (this results in a throw during setup). So for now bDouble
    // always matches the build precision.
    //
    // Since reading trajectory data from any-precision builds is wanted we need to expand
    // this and then set bDouble from the actual precision stored in the data sets.
    //
    // See issue #5474
#    if GMX_DOUBLE
    frame->bDouble = true;
#    else
    frame->bDouble = false;
#    endif
    frame->bLambda = false;
    // TODO: This should be read from the file
    frame->bPrec = false;
    frame->prec  = 0.0;

    TrajectoryReadCursor& readCursor = particleBlocks_.at(selectionName);
    if (readCursor.nextFrameContents(
                &frame->bX, &frame->bV, &frame->bF, &frame->bBox, &frame->bStep, &frame->bTime))
    {
        ArrayRef<RVec> positions{};
        ArrayRef<RVec> velocities{};
        ArrayRef<RVec> forces{};

        frame->natoms = readCursor.block().numParticles();
        if (frame->bX)
        {
            srenew(frame->x, frame->natoms);
            positions = arrayRefFromArray(reinterpret_cast<RVec*>(frame->x), frame->natoms);
        }
        if (frame->bV)
        {
            srenew(frame->v, frame->natoms);
            velocities = arrayRefFromArray(reinterpret_cast<RVec*>(frame->v), frame->natoms);
        }
        if (frame->bF)
        {
            srenew(frame->f, frame->natoms);
            forces = arrayRefFromArray(reinterpret_cast<RVec*>(frame->f), frame->natoms);
        }

        double timeAsDouble;
        frameWasRead = readCursor.readNextFrame(
                positions, velocities, forces, frame->box, &frame->step, &timeAsDouble);

        if (frame->bTime)
        {
            frame->time = timeAsDouble;
        }
    }

    return frameWasRead;
#else
    throw gmx::NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
    GMX_UNUSED_VALUE(frame);
    GMX_UNUSED_VALUE(selectionName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::writeNextFrame(ArrayRef<const RVec> positions,
                          ArrayRef<const RVec> velocities,
                          ArrayRef<const RVec> forces,
                          const matrix         box,
                          const int64_t        step,
                          const double         time)
{
#if GMX_USE_HDF5
    // Helper for error message creation when writing each array to the data sets
    const auto blockNotFoundError = [&](const char* blockType)
    {
        return gmx::formatString(
                "Cannot write %s for group '%s': no data set exists", blockType, c_fullSystemGroupName);
    };

    // We don't want to do a try-catch block or expensive search for every write,
    // since we should always find the desired block after the H5md file has been set up.
    // For extra debugging we leave this assert as a bread crumb.
    GMX_ASSERT(particleBlocks_.find(c_fullSystemGroupName) != particleBlocks_.cend(),
               "Could not find particle block when writing trajectory data");
    H5mdParticleBlock& particleBlock = particleBlocks_.at(c_fullSystemGroupName).block();
    if (!positions.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!particleBlock.hasPosition(), blockNotFoundError("positions"));
        particleBlock.position()->writeNextFrame(positions, step, time);
        if (particleBlock.hasBox())
        {
            particleBlock.box()->writeNextFrame(
                    constArrayRefFromArray<real>(reinterpret_cast<const real*>(box), DIM * DIM));
        }
    }
    if (!velocities.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!particleBlock.hasVelocity(), blockNotFoundError("velocities"));
        particleBlock.velocity()->writeNextFrame(velocities, step, time);
    }
    if (!forces.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!particleBlock.hasForce(), blockNotFoundError("forces"));
        particleBlock.force()->writeNextFrame(forces, step, time);
    }
#else
    GMX_UNUSED_VALUE(positions);
    GMX_UNUSED_VALUE(velocities);
    GMX_UNUSED_VALUE(forces);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(time);
#endif
}

#if GMX_USE_HDF5
// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool H5md::TrajectoryReadCursor::nextFrameContents(bool* hasPosition,
                                                   bool* hasVelocity,
                                                   bool* hasForce,
                                                   bool* hasBox,
                                                   bool* hasStep,
                                                   bool* hasTime)
{
    const std::optional<int64_t> positionStep =
            block_.hasPosition() ? block_.position()->readStepAtIndex(nextPositionFrameToRead_)
                                 : std::nullopt;
    const std::optional<int64_t> velocityStep =
            block_.hasVelocity() ? block_.velocity()->readStepAtIndex(nextVelocityFrameToRead_)
                                 : std::nullopt;
    const std::optional<int64_t> forceStep =
            block_.hasForce() ? block_.force()->readStepAtIndex(nextForceFrameToRead_) : std::nullopt;

    // Positions, velocities and forces may be sampled at different frequencies.
    // To determine the next frame to read we need to calculate the minimum simulation
    // step of the next frame for all trajectory data. All data sets whose next simulation
    // step matches this minimum are those which we can read data from for the next frame.
    std::optional<int64_t> minimumStep = std::nullopt;
    for (const auto& blockStep : { positionStep, velocityStep, forceStep })
    {
        // The order in this condition is important: we first check if we *do not* have
        // a minimum step, in which case we always take the current block value.
        // Only after we have found a step from a block do we compare it to the new block.
        if (!minimumStep.has_value() || (blockStep.has_value() && blockStep.value() < minimumStep.value()))
        {
            minimumStep = blockStep;
        }
    }

    if (minimumStep.has_value())
    {
        *hasPosition = positionStep.has_value() && (positionStep.value() == minimumStep.value());
        *hasVelocity = velocityStep.has_value() && (velocityStep.value() == minimumStep.value());
        *hasForce    = forceStep.has_value() && (forceStep.value() == minimumStep.value());
        *hasBox      = *hasPosition;
        *hasStep     = true;
        *hasTime     = (*hasPosition && block_.position()->hasTime())
                   || (*hasVelocity && block_.velocity()->hasTime())
                   || (*hasForce && block_.force()->hasTime());
        return true;
    }
    else
    {
        return false;
    }
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool H5md::TrajectoryReadCursor::readNextFrame(ArrayRef<RVec> positions,
                                               ArrayRef<RVec> velocities,
                                               ArrayRef<RVec> forces,
                                               matrix         box,
                                               int64_t*       step,
                                               double*        time)
{
    // Helper for error message creation when writing each array to the data sets
    const auto blockNotFoundError = [&](const char* blockType)
    {
        return gmx::formatString(
                "Cannot read %s for group '%s': no data set exists", blockType, c_fullSystemGroupName);
    };

    bool                   frameWasRead    = false;
    std::optional<int64_t> stepThatWasRead = std::nullopt;
    if (!positions.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!block_.hasPosition(), blockNotFoundError("positions"));
        frameWasRead = block_.position()->readFrame(nextPositionFrameToRead_, positions, step, time)
                       || frameWasRead;
        // TODO: For a constant box we must also read it!
        GMX_H5MD_THROW_UPON_ERROR(!block_.hasBox(), blockNotFoundError("box"));
        block_.box()->readFrame(nextPositionFrameToRead_,
                                arrayRefFromArray<real>(reinterpret_cast<real*>(box), DIM * DIM));

        GMX_H5MD_THROW_UPON_ERROR(
                stepThatWasRead.has_value() && *step != stepThatWasRead.value(),
                "Tried to read trajectory data for different simulation steps as the next frame");
        stepThatWasRead = *step;
        ++nextPositionFrameToRead_;
    }
    if (!velocities.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!block_.hasVelocity(), blockNotFoundError("velocities"));
        frameWasRead = block_.velocity()->readFrame(nextVelocityFrameToRead_, velocities, step, time)
                       || frameWasRead;

        GMX_H5MD_THROW_UPON_ERROR(
                stepThatWasRead.has_value() && *step != stepThatWasRead.value(),
                "Tried to read trajectory data for different simulation steps as the next frame");
        stepThatWasRead = *step;
        ++nextVelocityFrameToRead_;
    }
    if (!forces.empty())
    {
        GMX_H5MD_THROW_UPON_ERROR(!block_.hasForce(), blockNotFoundError("forces"));
        frameWasRead = block_.force()->readFrame(nextForceFrameToRead_, forces, step, time) || frameWasRead;

        GMX_H5MD_THROW_UPON_ERROR(
                stepThatWasRead.has_value() && *step != stepThatWasRead.value(),
                "Tried to read trajectory data for different simulation steps as the next frame");
        stepThatWasRead = *step;
        ++nextForceFrameToRead_;
    }
    return frameWasRead;
}

#endif

H5md* makeH5md(const std::filesystem::path& fileName, H5mdFileMode mode)
{
    return new gmx::H5md(fileName, mode);
}

void setupFileFromInput(H5md* h5md, const gmx_mtop_t& topology, const t_inputrec& inputRecord)
{
    h5md->setupFileFromInput(topology, inputRecord);
}

void writeNextFrame(H5md*                h5md,
                    ArrayRef<const RVec> positions,
                    ArrayRef<const RVec> velocities,
                    ArrayRef<const RVec> forces,
                    const matrix         box,
                    const int64_t        step,
                    const double         time)
{
    h5md->writeNextFrame(positions, velocities, forces, box, step, time);
}

//! Deallocate \c h5md
void destroyH5md(H5md* h5md)
{
    delete h5md;
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
