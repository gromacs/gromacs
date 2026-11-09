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
#include "gromacs/utility/futil.h"
#include "gromacs/utility/stringutil.h"

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "exceptions.h"
#    include "h5md_attribute.h"
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
//! \brief Attribute name for the author email address.
constexpr char c_h5mdEmailAttributeKey[] = "email";
//! \brief Attribute name for a version specification.
constexpr char c_h5mdVersionAttributeKey[] = "version";
//! \brief Path to H5MD connectivity group from the file root.
constexpr char c_h5mdConnectivityGroupPath[] = "/connectivity";
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

/*! \brief Return the number of particles for the block at \c /particles/selectionName.
 *
 * \note Does not check consistency of trajectory data sets in the block.
 *
 * \param[in] fileid        Handle to H5md file.
 * \param[in] selectionName Name of particle group.
 *
 * \returns The number of particles of trajectory data inside the block,
 * or std::nullopt if no trajectory data exists for the block.
 *
 * \throws H5mdError if there was an error reading the number of particles
 * from an existing data set.
 */
static std::optional<int64_t> calcNumParticlesForSelection(const hid_t fileid, const char* selectionName)
{
    const auto [particlesGroup, particlesGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid, c_particlesGroupPath));
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(openGroup(particlesGroup, selectionName));

    for (const char* dataGroupName : { c_positionGroupName, c_velocityGroupName, c_forceGroupName })
    {
        const std::string valueDataSetPath = joinStrings({ dataGroupName, c_valueName }, "/");
        if (objectExists(selectionGroup, valueDataSetPath.c_str()))
        {
            const H5mdFrameDataSet<RVec> valueDataSet(selectionGroup, valueDataSetPath.c_str());
            GMX_H5MD_THROW_UPON_ERROR(
                    valueDataSet.frameDims().empty(),
                    formatString(
                            "Incorrect dimensions for data set at '%s': must not be 0-dimensional",
                            joinStrings({ c_particlesGroupPath, selectionName, dataGroupName, c_valueName },
                                        "/")
                                    .c_str()));
            return valueDataSet.frameDims()[0];
        }
    }
    return std::nullopt;
}

/*! \brief Return the PbcType for the block at \c /particles/selectionName/box.
 *
 * \param[in] fileid        Handle to H5md file.
 * \param[in] selectionName Name of particle group.
 *
 * \returns The PbcType or std::nullopt if it does not exist.
 *
 * \throws H5mdError if the simulation box group at \c /particles/selectionName/box
 * exists but the PbcType could not be determined from it.
 */
static std::optional<PbcType> getPbcForSelection(const hid_t fileid, const char* selectionName)
{
    const auto [particlesGroup, particlesGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid, c_particlesGroupPath));
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(openGroup(particlesGroup, selectionName));

    if (!objectExists(selectionGroup, c_boxGroupName))
    {
        return std::nullopt;
    }
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(openGroup(selectionGroup, c_boxGroupName));

    const std::optional<std::vector<std::string>> boundary =
            getAttributeVector<std::string>(boxGroup, c_boxBoundaryAttributeKey);
    GMX_H5MD_THROW_UPON_ERROR(!boundary.has_value(),
                              formatString("Cannot read periodicity for group '%s': required "
                                           "'boundary' attribute does not exist",
                                           getHandlePath(boxGroup).value_or("<unknown_path>").c_str()));

    const std::vector<std::string>& pbc = boundary.value();
    GMX_H5MD_THROW_UPON_ERROR(
            pbc.size() != DIM,
            formatString("Cannot read periodicity for group '%s': required 'boundary' "
                         "attribute must have 3 values for periodicity",
                         getHandlePath(boxGroup).value_or("<unknown_path>").c_str()));
    constexpr char periodicValue[] = "periodic";
    constexpr char noneValue[]     = "none";

    if (pbc[XX] == periodicValue && pbc[YY] == periodicValue && pbc[ZZ] == periodicValue)
    {
        return PbcType::Xyz;
    }
    else if (pbc[XX] == periodicValue && pbc[YY] == periodicValue && pbc[ZZ] == noneValue)
    {
        return PbcType::XY;
    }
    else if (pbc[XX] == noneValue && pbc[YY] == noneValue && pbc[ZZ] == noneValue)
    {
        return PbcType::No;
    }
    else
    {
        GMX_THROW(
                H5mdError(formatString("Cannot read periodicity for group '%s': boundary attribute "
                                       "does not map to GROMACS PBC types",
                                       getHandlePath(boxGroup).value_or("<unknown path>").c_str())));
    }
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
        {
            const auto [createPropertyList, createPropertyListGuard] =
                    makeH5mdPropertyListGuard(H5Pcreate(H5P_FILE_CREATE));
            make_backup(fileName);
            file_ = H5Fcreate(fileName.string().c_str(), H5F_ACC_TRUNC, createPropertyList, H5P_DEFAULT);

            /* Set up the H5MD metadata according to the H5md specification:
             *
             * /h5md               (H5md metadata group)
             * \++ version         (H5md specification version attribute)
             * \-- author          (author group)
             *     \++ name        (author name, or a dummy name if we cannot get it)
             *     \++ email       (author email, if available)
             * \-- creator         (creator group)
             *     \++ name        (GROMACS)
             *     \++ version     (GROMACS version)
             *
             * These groups and attributes are required by the H5MD specification to exist
             * in the file. After setting these up we can proceed with the rest of the file.
             */
            const auto [group, groupGuard] =
                    makeH5mdGroupGuard(createGroup(file_, c_h5mdMetaDataGroupName));
            setAttributeVector<int>(group, c_h5mdVersionAttributeKey, c_h5mdSpecificationVersion);

            const auto [authorGroup, authorGroupGuard] =
                    makeH5mdGroupGuard(createGroup(group, c_h5mdAuthorGroupName));
            // Author name and email are optionally set from environment variables,
            // write if available. The H5MD spec requires a name to be written, so
            // use a fallback if it is not set.
            const char* authorName = std::getenv("GMX_AUTHOR_NAME");
            setAttribute(authorGroup, c_h5mdNameAttributeKey, authorName != nullptr ? authorName : "N/A");
            if (const char* authorEmail = std::getenv("GMX_AUTHOR_EMAIL"); authorEmail != nullptr)
            {
                setAttribute(authorGroup, c_h5mdEmailAttributeKey, authorEmail);
            }

            const auto [creatorGroup, creatorGroupGuard] =
                    makeH5mdGroupGuard(createGroup(group, c_h5mdCreatorGroupName));
            setAttribute(creatorGroup, c_h5mdNameAttributeKey, "GROMACS");
            setAttribute(creatorGroup, c_h5mdVersionAttributeKey, gmx_version());
        }
        break;
        case H5mdFileMode::Read:
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
            break;
        case H5mdFileMode::Append:
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            break;
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
        flush<false>();
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
    throw NotImplementedError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template<bool throwExceptionUponError>
void H5md::flush()
{
#if GMX_USE_HDF5
    if constexpr (throwExceptionUponError)
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
 * \param[in]  fileid        Handle to H5md file.
 * \param[in]  groupName     Name of particle group.
 * \param[out] particleBlock Particle block to create box data set for.
 *
 * \throws H5mdError if there was an error setting up the links to the position
 * step or time data sets.
 */
static void setupSimulationBoxDataSet(const hid_t fileid, const char* groupName, H5mdParticleBlock* particleBlock)
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
    const std::string selectionGroupPath = joinStrings({ c_particlesGroupPath, groupName }, "/");
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(openGroup(fileid, selectionGroupPath.c_str()));
    const auto [boxGroup, boxGroupGuard] = makeH5mdGroupGuard(openGroup(selectionGroup, c_boxGroupName));
    const auto [edgesGroup, edgesGroupGuard] = makeH5mdGroupGuard(createGroup(boxGroup, c_boxSizeName));

    // Matrices are stored as real[DIM][DIM], which is the data set frame dimension
    particleBlock->addBox(H5mdFrameDataSetBuilder<real>(edgesGroup, c_valueName)
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

    const int64_t numParticles = selectionIndices.empty() ? topology.natoms : selectionIndices.size();
    GMX_H5MD_THROW_UPON_ERROR(
            numParticles == 0,
            formatString("Cannot setup particle group '%s': no atoms in group", selectionName.c_str()));

    particleBlocks_.insert({ selectionName, TrajectoryReadCursor(H5mdParticleBlock(numParticles)) });
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
    writeMoleculeBlocks(gmxTop, makeConstArrayRef(topology.molblock), makeConstArrayRef(topology.moltype));
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
void H5md::setupModules()
{
#if GMX_USE_HDF5
    const auto [group, groupGuard] = makeH5mdGroupGuard(openGroup(file_, c_h5mdMetaDataGroupName));
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
    setupFileFromInput(topology, inputRecord, {}, c_fullSystemGroupName);
#else
    GMX_UNUSED_VALUE(topology);
    GMX_UNUSED_VALUE(inputRecord);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFileFromInput(const gmx_mtop_t&   topology,
                              const t_inputrec&   inputRecord,
                              ArrayRef<const int> selectionInds,
                              const char*         selectionName)
{
#if GMX_USE_HDF5
    setupModules();
    setupParticleBlockForGroup(topology, selectionInds, selectionName, inputRecord);
    setupGromacsTopology(topology);
    setupBondConnectivity(topology);
#else
    GMX_UNUSED_VALUE(topology);
    GMX_UNUSED_VALUE(inputRecord);
    GMX_UNUSED_VALUE(selectionInds);
    GMX_UNUSED_VALUE(selectionName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFromOtherFile(const H5md& other, const ArrayRef<const int> selectionInds, const char* selectionName)
{
#if GMX_USE_HDF5
    const std::optional<int64_t> numParticles =
            calcNumParticlesForSelection(other.fileid(), selectionName);
    GMX_H5MD_THROW_UPON_ERROR(!numParticles.has_value(),
                              formatString("Cannot read number of particles from input H5MD file"));
    gmx_mtop_t topology;
    topology.natoms = numParticles.value();

    const std::optional<PbcType> pbcType = getPbcForSelection(other.fileid(), selectionName);
    GMX_H5MD_THROW_UPON_ERROR(!pbcType.has_value(), "Cannot read PbcType from input H5MD file");
    t_inputrec inputRecord;
    inputRecord.pbcType = pbcType.value();

    // NOTE: Only copy over the topology data if the new file contains data
    // for the same number of particles. Until we figure out how to deal with
    // sub group selections we should ensure that the topology represents all
    // trajectory data.
    const bool copyTopology = selectionInds.empty() || numParticles.value() == gmx::ssize(selectionInds);
    copyMetaDataFromFile(other, copyTopology);
    setupParticleBlockForGroup(topology, selectionInds, selectionName, inputRecord);
#else
    GMX_UNUSED_VALUE(other);
    GMX_UNUSED_VALUE(selectionInds);
    GMX_UNUSED_VALUE(selectionName);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupParticleBlockForGroupFromExistingFile(const std::string& selectionName,
                                                      const std::optional<int64_t> restartingFromStep,
                                                      const int64_t numParticles)
{
#if GMX_USE_HDF5
    const auto [particlesGroup, particlesGroupGuard] =
            makeH5mdGroupGuard(openGroup(file_, c_particlesGroupPath));
    const auto [selectionGroup, selectionGroupGuard] =
            makeH5mdGroupGuard(openGroup(particlesGroup, selectionName.c_str()));

    H5mdParticleBlock particleBlock(numParticles);
    if (objectExists(selectionGroup, c_positionGroupName))
    {
        particleBlock.addPosition(H5mdTimeDataBlock<RVec>(selectionGroup, c_positionGroupName));

        // Right now we always write the box size for every position frame
        // and thus must always find a corresponding data set here.
        // TODO: Revisit as needed when we've decided how to treat this when writing.

        // Construct the box size group name relative to the current selection group,
        // then open the box size value data set inside
        const std::string boxSizeGroupPath = joinStrings({ c_boxGroupName, c_boxSizeName }, "/");
        const auto [boxSizeGroup, boxSizeGroupGuard] =
                makeH5mdGroupGuard(openGroup(selectionGroup, boxSizeGroupPath.c_str()));
        particleBlock.addBox(H5mdFrameDataSet<real>(boxSizeGroup, c_valueName));
    }
    if (objectExists(selectionGroup, c_velocityGroupName))
    {
        particleBlock.addVelocity(H5mdTimeDataBlock<RVec>(selectionGroup, c_velocityGroupName));
    }
    if (objectExists(selectionGroup, c_forceGroupName))
    {
        particleBlock.addForce(H5mdTimeDataBlock<RVec>(selectionGroup, c_forceGroupName));
    }

    if (restartingFromStep.has_value())
    {
        particleBlock.trimDataSetsFromStep(restartingFromStep.value());
    }

    particleBlocks_.insert({ selectionName, TrajectoryReadCursor(std::move(particleBlock)) });

#else
    GMX_UNUSED_VALUE(selectionName);
    GMX_UNUSED_VALUE(numParticles);
    GMX_UNUSED_VALUE(restartingFromStep);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFromExistingFile()
{
#if GMX_USE_HDF5
    const std::optional<int64_t> numParticles =
            calcNumParticlesForSelection(fileid(), c_fullSystemGroupName);
    GMX_H5MD_THROW_UPON_ERROR(
            !numParticles.has_value(),
            formatString(
                    "Cannot set up particle block for group '%s': no trajectory data is available",
                    c_fullSystemGroupName));

    setupParticleBlockForGroupFromExistingFile(c_fullSystemGroupName, std::nullopt, numParticles.value());
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::setupFromExistingFileForAppending(const int64_t restartingFromStep, const int64_t numParticles)
{
#if GMX_USE_HDF5
    setupParticleBlockForGroupFromExistingFile(c_fullSystemGroupName, restartingFromStep, numParticles);
#else
    GMX_UNUSED_VALUE(restartingFromStep);
    GMX_UNUSED_VALUE(numParticles);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
void H5md::copyMetaDataFromFile(const H5md& other, const bool copyTopology)
{
#if GMX_USE_HDF5
    // Helper to create full paths to a H5MD module from its name,
    // i.e. from '<name>' to '/h5md/modules/<name>'
    auto pathToModule = [&](const char* moduleName) -> std::string
    { return joinStrings({ c_h5mdMetaDataGroupName, c_h5mdModulesGroupName, moduleName }, "/"); };

    std::vector<std::string> objectsToCopy = {
        pathToModule("gromacs"),
        pathToModule("units"),
    };

    if (copyTopology)
    {
        objectsToCopy.push_back(pathToModule("gromacs_topology"));
        objectsToCopy.push_back(c_h5mdConnectivityGroupPath);
    }

    // To copy nested paths we need to set the option to create intermediate
    // groups to the link creation property list (lcpl_id in the HDF5 docs).
    const auto [linkCreationPropertyList, linkCreationPropertyListGuard] =
            makeH5mdPropertyListGuard(H5Pcreate(H5P_LINK_CREATE));
    H5Pset_create_intermediate_group(linkCreationPropertyList, 1);

    for (const std::string& path : objectsToCopy)
    {
        if (objectExists(other.fileid(), path.c_str()))
        {
            GMX_H5MD_THROW_UPON_ERROR(
                    H5Ocopy(other.fileid(), path.c_str(), fileid(), path.c_str(), H5P_DEFAULT, linkCreationPropertyList)
                            < 0,
                    formatString("Encountered error when copying object at '%s'.", path.c_str()));
        }
    }

#else
    GMX_UNUSED_VALUE(other);
    GMX_UNUSED_VALUE(copyTopology);
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
bool H5md::readNextFrame(t_trxframe* frame, const std::string& selectionName)
{
#if GMX_USE_HDF5
    bool frameWasRead = false;

    frame->bLambda = false;
    // TODO: This should be read from the file
    frame->bPrec = false;
    frame->prec  = 0.0;
    // TODO: Due to current limitations, time is currently always read from double
    // precision data sets, but after this is addressed we need to read this from
    // the open file. See #5474.
    frame->timeIsDouble = true;

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

        frameWasRead = readCursor.readNextFrame(
                positions, velocities, forces, frame->box, &frame->step, &frame->time);
    }

    return frameWasRead;
#else
    throw NotImplementedError(
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
    // We don't want to do a try-catch block or expensive search for every write,
    // since we should always find the desired block after the H5md file has been set up.
    // For extra debugging we leave this assert as a bread crumb.
    GMX_ASSERT(particleBlocks_.find(c_fullSystemGroupName) != particleBlocks_.cend(),
               "Could not find particle block when writing trajectory data");

    const auto createTimeDataBlock = [&](const char* name, const char* unit, const ArrayRef<const RVec> values)
    {
        const std::string selectionGroupPath =
                joinStrings({ c_particlesGroupPath, c_fullSystemGroupName }, "/");
        const auto [selectionGroup, selectionGroupGuard] =
                makeH5mdGroupGuard(openGroup(fileid(), selectionGroupPath.c_str()));

        return H5mdTimeDataBlockBuilder<RVec>(selectionGroup, name)
                .withFrameDimension({ values.size() })
                .withUnit(unit)
                .build();
    };

    H5mdParticleBlock& particleBlock = particleBlocks_.at(c_fullSystemGroupName).block();
    if (!positions.empty())
    {
        if (!particleBlock.hasPosition())
        {
            // Check this before creating the time data block to not fill
            // the file with garbage in case the numbers do not match (also below)
            GMX_H5MD_THROW_UPON_ERROR(
                    gmx::ssize(positions) != particleBlock.numParticles(),
                    "Number of positions does not match number of particles in block");
            particleBlock.addPosition(createTimeDataBlock(c_positionGroupName, c_positionUnit, positions));
            setupSimulationBoxDataSet(fileid(), c_fullSystemGroupName, &particleBlock);
        }

        particleBlock.position()->writeNextFrame(positions, step, time);
        if (particleBlock.hasBox())
        {
            particleBlock.box()->writeNextFrame(
                    constArrayRefFromArray<real>(reinterpret_cast<const real*>(box), DIM * DIM));
        }
    }
    if (!velocities.empty())
    {
        if (!particleBlock.hasVelocity())
        {
            GMX_H5MD_THROW_UPON_ERROR(
                    gmx::ssize(velocities) != particleBlock.numParticles(),
                    "Number of velocities does not match number of particles in block");
            particleBlock.addVelocity(createTimeDataBlock(c_velocityGroupName, c_velocityUnit, velocities));
        }

        particleBlock.velocity()->writeNextFrame(velocities, step, time);
    }
    if (!forces.empty())
    {
        if (!particleBlock.hasForce())
        {
            GMX_H5MD_THROW_UPON_ERROR(
                    gmx::ssize(forces) != particleBlock.numParticles(),
                    "Number of forces does not match number of particles in block");
            particleBlock.addForce(createTimeDataBlock(c_forceGroupName, c_forceUnit, forces));
        }

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
        return formatString(
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
    return new H5md(fileName, mode);
}

void setupFileFromInput(H5md* h5md, const gmx_mtop_t& topology, const t_inputrec& inputRecord)
{
    h5md->setupFileFromInput(topology, inputRecord);
}

void setupFromExistingFileForAppending(H5md* h5md, const int64_t restartingFromStep, const int64_t numParticles)
{
    h5md->setupFromExistingFileForAppending(restartingFromStep, numParticles);
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

//! Flush \c h5md file.
void flushH5md(H5md* h5md)
{
    h5md->flush();
}

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
