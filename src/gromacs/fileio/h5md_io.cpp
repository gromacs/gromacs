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
/* This file was inspired by ch5md by Pierre de Buyl (BSD license). */

#include "gmxpre.h"

#include "h5md_io.h"

#include "config.h"

#include <cmath>
#include <cstring>

#include <algorithm>
#include <functional>
#include <limits>
#include <string>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

#include "h5md_datablock.h"
#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "external/SZ3-bio/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif

/*! \brief Iterates through groups with contents matching time dependent particles data blocks,
 * i.e., "step", "time" and "value". Then it creates corresponding H5MD data blocks.
 * Inspired by https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5G/h5ex_g_traverse.c
 */
static herr_t iterativeSetupTimeDataBlocks(hid_t            locationId,
                                           const char*      name,
                                           const H5L_info_t gmx_unused* info,
                                           void*                        operatorData)
{
    /*
     * Get type of the object. The name of the object is passed to this function by
     * the Library.
     */
    H5O_info_t infoBuffer;
    H5Oget_info_by_name(locationId, name, &infoBuffer, H5P_DEFAULT);
    herr_t            returnVal        = 0;
    const std::string stepDataSetName  = std::string(name) + std::string("/step");
    const std::string timeDataSetName  = std::string(name) + std::string("/time");
    const std::string valueDataSetName = std::string(name) + std::string("/value");
    switch (infoBuffer.type)
    {
        case H5O_TYPE_GROUP:
            if (objectExists(locationId, stepDataSetName.c_str())
                && objectExists(locationId, timeDataSetName.c_str())
                && objectExists(locationId, valueDataSetName.c_str()))
            {
                char containerFullName[c_maxFullNameLength];
                H5Iget_name(locationId, containerFullName, c_maxFullNameLength);
                GmxH5mdTimeDataBlock             dataBlock(locationId, name);
                std::list<GmxH5mdTimeDataBlock>* dataBlocks =
                        static_cast<std::list<GmxH5mdTimeDataBlock>*>(operatorData);

                dataBlock.updateNumWrittenFrames();
                dataBlocks->emplace_back(dataBlock);

                returnVal = 0;
            }
            else
            {
                returnVal = H5Literate_by_name(locationId,
                                               name,
                                               H5_INDEX_NAME,
                                               H5_ITER_NATIVE,
                                               nullptr,
                                               iterativeSetupTimeDataBlocks,
                                               operatorData,
                                               H5P_DEFAULT);
            }
            break;
        default: /* Ignore other contents */ break;
    }
    return returnVal;
}


GmxH5mdIo::GmxH5mdIo(const std::string fileName, const char mode)
{
    file_ = -1;
    if (fileName.length() > 0)
    {
        openFile(fileName.c_str(), mode);
    }
}

GmxH5mdIo::~GmxH5mdIo()
{
    if (file_ != -1)
    {
        closeFile();
    }
}

void GmxH5mdIo::openFile(const std::string fileName, const char mode)
{
#if GMX_USE_HDF5
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error output, e.g. when items are not found.

    closeFile();

    systemOutputName_ = "system";
    dataBlocks_.clear();

    if (debug)
    {
        fprintf(debug, "Opening H5MD file %s with mode %c\n", fileName.c_str(), mode);
    }
    if (mode == 'w' || mode == 'a')
    {
        bool fileExists = gmx_fexist(fileName);
        if (!fileExists || mode == 'w')
        {
            make_backup(fileName.c_str());
            hid_t createPropertyList = H5Pcreate(H5P_FILE_CREATE);
            file_ = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC, createPropertyList, H5P_DEFAULT);
            if (file_ < 0)
            {
                throw gmx::FileIOError("Cannot create H5MD file.");
            }
            setAuthorAndCreator();
        }
        else
        {
            file_ = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        /* Create H5MD groups. They should already be there if appending to a valid H5MD file, but it's better to be on the safe side. */
        openOrCreateGroup(file_, "h5md");
        openOrCreateGroup(file_, "particles");
        openOrCreateGroup(file_, "particles/system");
    }
    else
    {
        file_ = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    if (file_ < 0)
    {
        throw gmx::FileIOError("Cannot open H5MD file.");
    }
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::closeFile()
{
#if GMX_USE_HDF5
    if (file_ >= 0)
    {
        if (H5Fflush(file_, H5F_SCOPE_LOCAL) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error flushing H5MD file when closing.");
        }
        if (debug)
        {
            fprintf(debug, "Closing H5MD file.\n");
        }
        for (auto dataBlock : dataBlocks_)
        {
            dataBlock.closeAllDataSets();
        }
        H5Fclose(file_);
        file_ = -1;
    }
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::flush()
{
#if GMX_USE_HDF5
    if (file_ >= 0)
    {
        if (debug)
        {
            fprintf(debug, "Flushing H5MD file.\n");
        }
        if (H5Fflush(file_, H5F_SCOPE_LOCAL) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error flushing H5MD file when closing.");
        }
    }
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::initParticleDataBlocksFromFile()
{
    hid_t particlesGroup = H5Gopen(file_, "particles", H5P_DEFAULT);
    if (particlesGroup < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError(
                "Cannot find particles group when initializing particles data blocks. Invalid H5MD "
                "file?");
    }
    if (H5Literate(particlesGroup,
                   H5_INDEX_NAME,
                   H5_ITER_NATIVE,
                   nullptr,
                   iterativeSetupTimeDataBlocks,
                   static_cast<void*>(&dataBlocks_))
        < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Error iterating over particles data blocks.");
    }
}

void GmxH5mdIo::setUpParticlesDataBlocks(int     writeCoordinatesSteps,
                                         int     writeVelocitiesSteps,
                                         int     writeForcesSteps,
                                         int64_t numParticles,
                                         PbcType pbcType,
                                         double  xCompressionError)
{
    if (numParticles <= 0)
    {
        throw gmx::FileIOError("There must be particles/atoms when writing trajectory frames.");
    }
#if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif

    std::string name        = "particles/" + systemOutputName_;
    hid_t       systemGroup = openOrCreateGroup(file_, name.c_str());

    hsize_t              numFramesPerChunk    = 1;
    CompressionAlgorithm compressionAlgorithm = CompressionAlgorithm::LosslessWithShuffle;
    if (writeCoordinatesSteps > 0)
    {
        if (xCompressionError > 0)
        {
            /* Use no more than 10 frames per chunk (compression unit). Use fewer frames per chunk if there are many atoms. */
            numFramesPerChunk    = std::min(10, int(std::ceil(2e6f / numParticles)));
            compressionAlgorithm = CompressionAlgorithm::LossySz3;

            /* Register the SZ3 filter. This is not necessary when creating a dataset with the filter,
             * but must be done to append to an existing file (e.g. when restarting from checkpoint). */
            registerSz3FilterImplicitly();
        }


        hid_t boxGroup = openOrCreateGroup(systemGroup, "box");
        setBoxGroupAttributes(boxGroup, pbcType);
        GmxH5mdTimeDataBlock box(boxGroup,
                                 "edges",
                                 "nm",
                                 writeCoordinatesSteps,
                                 numFramesPerChunk,
                                 DIM,
                                 DIM,
                                 datatype,
                                 CompressionAlgorithm::LosslessNoShuffle, // Never compress box output
                                 0);
        dataBlocks_.emplace_back(box);

        GmxH5mdTimeDataBlock position(systemGroup,
                                      "position",
                                      "nm",
                                      writeCoordinatesSteps,
                                      numFramesPerChunk,
                                      numParticles,
                                      DIM,
                                      datatype,
                                      compressionAlgorithm,
                                      xCompressionError);
        dataBlocks_.emplace_back(position);
    }
    numFramesPerChunk    = 1;
    compressionAlgorithm = CompressionAlgorithm::LosslessWithShuffle;
    if (writeForcesSteps > 0)
    {
        GmxH5mdTimeDataBlock force(systemGroup,
                                   "force",
                                   "kJ mol-1 nm-1",
                                   writeForcesSteps,
                                   numFramesPerChunk,
                                   numParticles,
                                   DIM,
                                   datatype,
                                   compressionAlgorithm,
                                   0);
        dataBlocks_.emplace_back(force);
    }
    if (writeVelocitiesSteps > 0)
    {
        GmxH5mdTimeDataBlock velocity(systemGroup,
                                      "velocity",
                                      "nm ps-1",
                                      writeVelocitiesSteps,
                                      numFramesPerChunk,
                                      numParticles,
                                      DIM,
                                      datatype,
                                      compressionAlgorithm,
                                      0);
        dataBlocks_.emplace_back(velocity);
    }
}

void GmxH5mdIo::setAuthorAndCreator()
{
    std::string precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif

    std::string programInfo  = gmx::getProgramContext().displayName() + precisionString;
    hid_t       creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "name", programInfo.c_str());
    const std::string gmxVersion = gmx_version();
    setAttribute(creatorGroup, "version", gmxVersion.c_str());
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    char  username[c_maxFullNameLength];
    if (!gmx_getusername(username, c_maxFullNameLength))
    {
        setAttribute(authorGroup, "name", username);
    }
}

void GmxH5mdIo::setupMolecularSystem(const gmx_mtop_t&        topology,
                                     gmx::ArrayRef<const int> index,
                                     const std::string        index_group_name)
{
#if GMX_USE_HDF5
    if (file_ < 0)
    {
        throw gmx::FileIOError("No file open for writing");
    }

    /* Vectors are used to keep the values in a continuous memory block. */
    std::vector<real>        atomCharges;
    std::vector<real>        atomMasses;
    std::vector<std::string> atomNames;

    t_atoms atoms = gmx_mtop_global_atoms(topology);

    atomCharges.reserve(atoms.nr);
    atomMasses.reserve(atoms.nr);
    atomNames.reserve(atoms.nr);

    /* FIXME: The names could be copied directly to a char array instead. */
    for (int atomCounter = 0; atomCounter < atoms.nr; atomCounter++)
    {
        atomCharges.push_back(atoms.atom[atomCounter].q);
        atomMasses.push_back(atoms.atom[atomCounter].m);
        atomNames.push_back(*(atoms.atomname[atomCounter]));
    }

    hsize_t atomPropertiesChunkDims[1];

    /* Don't replace data that already exists and cannot (currently) change during the simulation */
    /* FIXME: Currently atom names cannot change during the simulation. */
    if (!H5Lexists(file_, "/particles/system/atomname", H5P_DEFAULT))
    {
        /* Is there a more convenient way to do this? std::string is nice above, but cannot be used for writing in HDF5. */
        /* Hard-code the atom name lengths to max 17 (max 4 char 4-byte UTF8). Flexible strings make a lot of unaccounted space,
         * which is wasted. For strings that are numerous, such as atom names, it is better to use fixed-length. */
        constexpr size_t atomNameLen = 17;
        char*            atomNamesChars;
        snew(atomNamesChars, atomNames.size() * atomNameLen);
        for (size_t i = 0; i < atomNames.size(); i++)
        {
            strncpy(&atomNamesChars[i * atomNameLen], atomNames[i].c_str(), atomNameLen);
        }

        hid_t stringDataType = H5Tcopy(H5T_C_S1);
        H5Tset_size(stringDataType, atomNameLen);
        // H5Tset_strpad(stringDataType, H5T_STR_NULLTERM);
        H5Tset_cset(stringDataType, H5T_CSET_UTF8);

        // hsize_t atomPropertiesChunkDims[2];
        // atomPropertiesChunkDims[0] = 1;
        // atomPropertiesChunkDims[1] = topology.natoms;
        atomPropertiesChunkDims[0] = topology.natoms;

        hid_t atomName = openOrCreateDataSet<1>(file_,
                                                "/particles/system/atomname",
                                                "",
                                                stringDataType,
                                                atomPropertiesChunkDims,
                                                CompressionAlgorithm::LosslessNoShuffle,
                                                0);
        // writeData<1, true>(atomName, atomNamesChars.data(), 0);
        // printf("atomNamesChars.data()[0]: %s %s %s %s\n", atomNamesChars.data()[0], atomNamesChars.data()[1], atomNamesChars.data()[2], atomNamesChars.data()[3]);
        writeData<1, true>(atomName, atomNamesChars, 0);
        H5Dclose(atomName);
        sfree(atomNamesChars);
    }

#    if GMX_DOUBLE
    const hid_t floatDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#    else
    const hid_t floatDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#    endif
    /* FIXME: Currently charges and masses cannot change during the simulation. For time dependent data use GmxH5mdDataBlock */
    if (!H5Lexists(file_, "/particles/system/charge", H5P_DEFAULT))
    {
        hid_t charge = openOrCreateDataSet<1>(file_,
                                              "/particles/system/charge",
                                              "",
                                              floatDatatype,
                                              atomPropertiesChunkDims,
                                              CompressionAlgorithm::LosslessNoShuffle,
                                              0);
        writeData<1, true>(charge, atomCharges.data(), 0);
        H5Dclose(charge);
    }

    if (!H5Lexists(file_, "/particles/system/mass", H5P_DEFAULT))
    {
        hid_t mass = openOrCreateDataSet<1>(file_,
                                            "/particles/system/mass",
                                            "",
                                            floatDatatype,
                                            atomPropertiesChunkDims,
                                            CompressionAlgorithm::LosslessNoShuffle,
                                            0);
        writeData<1, true>(mass, atomMasses.data(), 0);
        H5Dclose(mass);
    }

    /* We only need to create a separate selection group entry if not all atoms are part of it. */
    /* TODO: Write atom name, charge and mass for the selection group as well. */
    /* If a selection of atoms is explicitly provided then use that instead of the CompressedPositionOutput */
    bool all_atoms_selected = true;
    if (index.ssize() > 0 && index.ssize() != topology.natoms)
    {
        all_atoms_selected = false;
    }
    else
    {
        for (int i = 0; (i < topology.natoms); i++)
        {
            if (getGroupType(topology.groups, SimulationAtomGroupType::CompressedPositionOutput, i) != 0)
            {
                all_atoms_selected = false;
                break;
            }
        }
    }
    bool setupSeparateOutputGroup = false;
    if (!all_atoms_selected)
    {
        if (index.ssize() > 0 && index_group_name != "")
        {
            setupSeparateOutputGroup = true;
            systemOutputName_        = index_group_name;
        }
        /* If no name was specified fall back to using the selection group name of compressed output, if any. */
        else if (topology.groups.numberOfGroupNumbers(SimulationAtomGroupType::CompressedPositionOutput) != 0)
        {
            setupSeparateOutputGroup = true;
            int nameIndex = topology.groups.groups[SimulationAtomGroupType::CompressedPositionOutput][0];
            systemOutputName_ = *topology.groups.groupNames[nameIndex];
        }
    }
    if (!setupSeparateOutputGroup)
    {
        systemOutputName_ = "system";
    }
    std::string name = "particles/" + systemOutputName_;
    openOrCreateGroup(file_, name.c_str());

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::writeFrame(int64_t      step,
                           real         time,
                           real         lambda,
                           const rvec*  box,
                           const int    numParticles,
                           const rvec*  x,
                           const rvec*  v,
                           const rvec*  f,
                           const double xCompressionError)
{
#if GMX_USE_HDF5
    if (numParticles <= 0)
    {
        throw gmx::FileIOError("There must be particles/atoms when writing trajectory frames.");
    }
    if (file_ < 0)
    {
        throw gmx::FileIOError("No file open for writing");
    }


#    if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#    else
    const hid_t datatype      = H5Tcopy(H5T_NATIVE_FLOAT);
#    endif

    std::string name        = "particles/" + systemOutputName_;
    hid_t       systemGroup = openOrCreateGroup(file_, name.c_str());

    CompressionAlgorithm compressionAlgorithm = CompressionAlgorithm::LosslessWithShuffle;
    hsize_t              numFramesPerChunk    = 1;
    if (x != nullptr)
    {
        std::string wantedName = "/particles/" + systemOutputName_ + "/position";
        auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
        if (foundDataBlock == dataBlocks_.end())
        {
            if (xCompressionError != 0)
            {
                /* Use no more than 10 frames per chunk (compression unit). Use fewer frames per chunk if there are many atoms. */
                numFramesPerChunk    = std::min(10, int(std::ceil(2e6f / numParticles)));
                compressionAlgorithm = CompressionAlgorithm::LossySz3;

                /* Register the SZ3 filter. This is not necessary when creating a dataset with the filter,
                * but must be done to append to an existing file (e.g. when restarting from checkpoint). */
                registerSz3FilterImplicitly();
            }
            GmxH5mdTimeDataBlock position(systemGroup,
                                          "position",
                                          "nm",
                                          -1,
                                          numFramesPerChunk,
                                          numParticles,
                                          DIM,
                                          datatype,
                                          compressionAlgorithm,
                                          xCompressionError);
            dataBlocks_.emplace_back(position);
            foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
            if (foundDataBlock == dataBlocks_.end())
            {
                throw gmx::FileIOError("Error creating position data block when writing frame.");
            }
        }
        foundDataBlock->writeFrame(x, step, time);

        if (box != nullptr)
        {
            std::string wantedName = "/particles/" + systemOutputName_ + "/box/edges";
            foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
            if (foundDataBlock == dataBlocks_.end())
            {
                hid_t boxGroup = openOrCreateGroup(systemGroup, "box");
                setBoxGroupAttributes(boxGroup, PbcType::Xyz);
                GmxH5mdTimeDataBlock box(boxGroup,
                                         "edges",
                                         "nm",
                                         -1,
                                         numFramesPerChunk,
                                         DIM,
                                         DIM,
                                         datatype,
                                         CompressionAlgorithm::LosslessNoShuffle, // Never compress box output
                                         0);
                dataBlocks_.emplace_back(box);
                foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
                if (foundDataBlock == dataBlocks_.end())
                {
                    throw gmx::FileIOError("Error creating box data block when writing frame.");
                }
            }
            foundDataBlock->writeFrame(box, step, time);
        }
    }

    numFramesPerChunk    = 1;
    compressionAlgorithm = CompressionAlgorithm::LosslessWithShuffle;
    if (v != nullptr)
    {
        std::string wantedName = "/particles/" + systemOutputName_ + "/velocity";
        auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
        if (foundDataBlock == dataBlocks_.end())
        {
            GmxH5mdTimeDataBlock velocity(systemGroup,
                                          "velocity",
                                          "nm ps-1",
                                          -1,
                                          numFramesPerChunk,
                                          numParticles,
                                          DIM,
                                          datatype,
                                          compressionAlgorithm,
                                          0);
            dataBlocks_.emplace_back(velocity);
            foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
            if (foundDataBlock == dataBlocks_.end())
            {
                throw gmx::FileIOError("Error creating velocity data block when writing frame.");
            }
        }
        foundDataBlock->writeFrame(v, step, time);
    }
    if (f != nullptr)
    {
        std::string wantedName = "/particles/" + systemOutputName_ + "/force";
        auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
        if (foundDataBlock == dataBlocks_.end())
        {
            GmxH5mdTimeDataBlock force(systemGroup,
                                       "force",
                                       "kJ mol-1 nm-1",
                                       -1,
                                       numFramesPerChunk,
                                       numParticles,
                                       DIM,
                                       datatype,
                                       compressionAlgorithm,
                                       0);
            dataBlocks_.emplace_back(force);
            foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
            if (foundDataBlock == dataBlocks_.end())
            {
                throw gmx::FileIOError("Error creating force data block when writing frame.");
            }
        }
        foundDataBlock->writeFrame(f, step, time);
    }
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

bool GmxH5mdIo::readNextFrameOfStandardDataBlocks(int64_t* step,
                                                  real*    time,
                                                  rvec*    box,
                                                  rvec*    x,
                                                  rvec*    v,
                                                  rvec*    f,
                                                  real*    xCompressionError,
                                                  bool*    readBox,
                                                  bool*    readX,
                                                  bool*    readV,
                                                  bool*    readF)
{
    std::string                     nameStem = "/particles/" + systemOutputName_;
    std::list<std::string>          dataBlockNames{ nameStem + "/box/edges",
                                           nameStem + "/position",
                                           nameStem + "/force",
                                           nameStem + "/velocity" };
    std::list<GmxH5mdTimeDataBlock> dataBlocksNextFrame{};
    int64_t                         minStepNextFrame = std::numeric_limits<int64_t>::max();
    for (std::string name : dataBlockNames)
    {
        auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), name.c_str());
        if (foundDataBlock == dataBlocks_.end())
        {
            continue;
        }
        int frameIndex = foundDataBlock->readingFrameIndex();
        if (frameIndex >= foundDataBlock->writingFrameIndex())
        {
            continue;
        }
        int64_t frameStep = foundDataBlock->getStepOfFrame(frameIndex);
        /* Discard data sets that had a higher time stamp if an earlier data point has been found. */
        if (frameStep < minStepNextFrame)
        {
            dataBlocksNextFrame.clear();
            minStepNextFrame = frameStep;
            *time            = foundDataBlock->getTimeOfFrame(frameIndex);
        }
        if (frameStep <= minStepNextFrame)
        {
            dataBlocksNextFrame.emplace_back(*foundDataBlock);
        }
    }
    *step              = minStepNextFrame;
    bool didReadFrame  = false;
    *xCompressionError = -1;
    for (std::list<GmxH5mdTimeDataBlock>::iterator dataBlock = dataBlocks_.begin();
         dataBlock != dataBlocks_.end();
         ++dataBlock)
    {
        /* FIXME: Can this be done more elegantly? */
        if (box != nullptr && dataBlock->name() == "edges")
        {
            if (dataBlock->readNextFrame(static_cast<real*>(box[0])))
            {
                *readBox     = true;
                didReadFrame = true;
            }
        }
        else if (x != nullptr && dataBlock->name() == "position")
        {
            if (dataBlock->readNextFrame(static_cast<real*>(x[0])))
            {
                *readX             = true;
                didReadFrame       = true;
                *xCompressionError = dataBlock->getLossyCompressionError();
            }
        }
        else if (v != nullptr && dataBlock->name() == "velocity")
        {
            if (dataBlock->readNextFrame(static_cast<real*>(v[0])))
            {
                *readV       = true;
                didReadFrame = true;
            }
        }
        else if (f != nullptr && dataBlock->name() == "force")
        {
            if (dataBlock->readNextFrame(static_cast<real*>(f[0])))
            {
                *readF       = true;
                didReadFrame = true;
            }
        }
        else
        {
            throw gmx::FileIOError("Unexpected data type.");
        }
    }
    return didReadFrame;
}


int64_t GmxH5mdIo::getNumberOfFrames(const std::string dataBlockName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + systemOutputName_ + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->numberOfFrames();
}

int64_t GmxH5mdIo::getNumberOfParticles(const std::string dataBlockName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + systemOutputName_ + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getNumParticles();
}

real GmxH5mdIo::getFirstTime(const std::string dataBlockName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + systemOutputName_ + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getTimeOfFrame(0);
}

real GmxH5mdIo::getFirstTimeFromAllDataBlocks()
{
    real firstTime = std::numeric_limits<real>::max();
    bool foundAny  = false;
    for (const auto& dataBlock : dataBlocks_)
    {
        foundAny  = true;
        real time = dataBlock.getTimeOfFrame(0);
        firstTime = std::min(firstTime, time);
    }
    if (foundAny)
    {
        return firstTime;
    }
    return -1;
}

real GmxH5mdIo::getFinalTime(const std::string dataBlockName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + systemOutputName_ + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getTimeOfFrame(foundDataBlock->numberOfFrames() - 1);
}

real GmxH5mdIo::getFinalTimeFromAllDataBlocks()
{
    real finalTime = 0;
    bool foundAny  = false;
    for (auto& dataBlock : dataBlocks_)
    {
        int numFrames = dataBlock.numberOfFrames();
        if (numFrames < 1)
        {
            continue;
        }
        foundAny  = true;
        real time = dataBlock.getTimeOfFrame(numFrames - 1);
        finalTime = std::max(finalTime, time);
    }
    if (foundAny)
    {
        return finalTime;
    }
    return -1;
}


extern template hid_t
openOrCreateDataSet<1>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);
extern template hid_t
openOrCreateDataSet<2>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);

extern template void writeData<1, true>(hid_t, const void*, hsize_t);

extern template void setAttribute<int>(hid_t, const char*, int, hid_t);
extern template void setAttribute<float>(hid_t, const char*, float, hid_t);
extern template void setAttribute<double>(hid_t, const char*, double, hid_t);
extern template void setAttribute<char*>(hid_t, const char*, char*, hid_t);
