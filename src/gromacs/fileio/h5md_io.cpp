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

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include <sys/_types/_int64_t.h>

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

namespace
{

/*! \brief Iterates through groups with contents matching time dependent particles data blocks,
 * i.e., "step", "time" and "value". Then it creates corresponding H5MD data blocks.
 * Inspired by https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/C/H5G/h5ex_g_traverse.c
 */
herr_t iterativeSetupTimeDataBlocks(hid_t            locationId,
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
            if (gmx::h5mdio::objectExists(locationId, stepDataSetName.c_str())
                && gmx::h5mdio::objectExists(locationId, timeDataSetName.c_str())
                && gmx::h5mdio::objectExists(locationId, valueDataSetName.c_str()))
            {
                char containerFullName[gmx::h5mdio::c_maxFullNameLength];
                H5Iget_name(locationId, containerFullName, gmx::h5mdio::c_maxFullNameLength);
                gmx::h5mdio::GmxH5mdTimeDataBlock             dataBlock(locationId, name);
                std::list<gmx::h5mdio::GmxH5mdTimeDataBlock>* dataBlocks =
                        static_cast<std::list<gmx::h5mdio::GmxH5mdTimeDataBlock>*>(operatorData);

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

void setupSystemParticleProperties(gmx::h5mdio::GmxH5mdIo*  file,
                                   const t_atoms&           atoms,
                                   gmx::ArrayRef<const int> selectionIndices,
                                   std::string              selectionName)
{
    /* Vectors are used to keep the values in a continuous memory block. */
    std::vector<real> atomCharges;
    std::vector<real> atomMasses;
    std::vector<int>  atomElements;
    /* Since the system block contains all atoms it is not necessary to record the ID,
     * but we do that in order to allow changing the mapping or "remove" particles,
     * in order to enable grand canonical simulations. */
    std::vector<int> atomIds;

    const size_t numSelectedParticles = selectionIndices.size() > 0 ? selectionIndices.size() : atoms.nr;

    atomCharges.reserve(numSelectedParticles);
    atomMasses.reserve(numSelectedParticles);
    atomElements.reserve(numSelectedParticles);
    atomIds.reserve(numSelectedParticles);

    /* FIXME: Should use int64_t. Needs changes in atoms. */
    for (size_t i = 0; i < numSelectedParticles; i++)
    {
        size_t iParticle = selectionIndices.size() > 0 ? selectionIndices[i] : i;
        atomCharges.push_back(atoms.atom[iParticle].q);
        atomMasses.push_back(atoms.atom[iParticle].m);
        atomElements.push_back(atoms.atom[iParticle].atomnumber);
        atomIds.push_back(iParticle);
    }

    file->setNumericProperty("/particles/" + selectionName, "charge", atomCharges, false);
    file->setNumericProperty("/particles/" + selectionName, "mass", atomMasses, false);
    file->setNumericProperty("/particles/" + selectionName, "species", atomElements, false);
    file->setNumericProperty("/particles/" + selectionName, "id", atomIds, false);
}

/* Unused. May be useful later. */
/*
herr_t getGroupNamesInLocation(hid_t location, const char* name, const H5O_info_t* info, void* operatorData)
{
    if(info->type == H5O_TYPE_GROUP)
    {
        auto vec = static_cast<std::vector<std::string>*>(operatorData);
        vec->push_back(std::string(name));
    }
}
*/

} // namespace

namespace gmx
{
namespace h5mdio
{

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
        }
        else
        {
            file_ = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        /* Create H5MD group. They should already be there if appending to a valid H5MD file, but it's better to be on the safe side. */
        hid_t h5mdGroup = openOrCreateGroup(file_, "h5md");
        setVersionAttribute(h5mdGroup, c_h5mdMajorVersion, c_h5mdMinorVersion);
    }
    else
    {
        file_ = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    initGroupTimeDataBlocksFromFile("particles");
    initGroupTimeDataBlocksFromFile("observables");
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

int GmxH5mdIo::initGroupTimeDataBlocksFromFile(std::string groupName)
{
    int   numDataBlocksBefore = dataBlocks_.size();
    hid_t group               = H5Gopen(file_, groupName.c_str(), H5P_DEFAULT);
    if (group < 0)
    {
        if (debug)
        {
            fprintf(debug,
                    "Cannot find group %s when initializing particles data blocks. Invalid file?",
                    groupName.c_str());
        }
        return 0;
    }
    if (H5Literate(group,
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
    return dataBlocks_.size() - numDataBlocksBefore;
}

std::string GmxH5mdIo::getH5mdRootVersionNumber()
{
    int   majorVersion, minorVersion;
    hid_t h5mdGroup = H5Gopen(file_, "h5md", H5P_DEFAULT);

    if (getVersionAttribute(h5mdGroup, &majorVersion, &minorVersion))
    {
        return std::to_string(majorVersion) + "." + std::to_string(minorVersion);
    }
    return "";
}

void GmxH5mdIo::setAuthor(std::string authorName)
{
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    setAttribute(authorGroup, "name", authorName.c_str());
}

std::string GmxH5mdIo::getAuthor()
{
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    char* tmpName     = nullptr;
    getAttribute(authorGroup, "name", &tmpName);
    std::string name(tmpName);
    free(tmpName);
    return name;
}

void GmxH5mdIo::setCreatorProgramName(std::string creatorName)
{
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "name", creatorName.c_str());
}

std::string GmxH5mdIo::getCreatorProgramName()
{
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    char* tmpName      = nullptr;
    getAttribute(creatorGroup, "name", &tmpName);
    std::string name(tmpName);
    free(tmpName);
    return name;
}

void GmxH5mdIo::setCreatorProgramVersion(std::string version)
{
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "version", version.c_str());
}

std::string GmxH5mdIo::getCreatorProgramVersion()
{
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    char* tmpVersion   = nullptr;
    getAttribute(creatorGroup, "version", &tmpVersion);
    std::string version(tmpVersion);
    free(tmpVersion);
    return version;
}

hid_t GmxH5mdIo::getGroupId(const std::string& name)
{
    hid_t group = H5Gopen(file_, name.c_str(), H5P_DEFAULT);

    return group;
}

hid_t GmxH5mdIo::createGroup(const std::string& name)
{
    hid_t group = openOrCreateGroup(file_, name.c_str());

    return group;
}

void GmxH5mdIo::setStringProperty(const std::string&              containerName,
                                  const std::string&              propertyName,
                                  const std::vector<std::string>& propertyValues,
                                  bool                            replaceExisting,
                                  size_t                          maxStringLength)
{
    openOrCreateGroup(file_, containerName.c_str());
    std::string dataSetName(containerName + "/" + propertyName);

    if (!H5Lexists(file_, dataSetName.c_str(), H5P_DEFAULT) || replaceExisting == true)
    {
        hid_t stringDataType = H5Tcopy(H5T_C_S1);
        H5Tset_cset(stringDataType, H5T_CSET_UTF8);
        hsize_t chunkDims[1];
        chunkDims[0] = propertyValues.size();
        if (maxStringLength > 0)
        {
            /* FIXME: Is there a more convenient way to do this? std::string is nice above, but cannot be used for writing in HDF5. */
            char* propertyValuesChars;
            snew(propertyValuesChars, propertyValues.size() * maxStringLength);
            for (size_t i = 0; i < propertyValues.size(); i++)
            {
                strncpy(&propertyValuesChars[i * maxStringLength], propertyValues[i].c_str(), maxStringLength);
            }

            H5Tset_size(stringDataType, maxStringLength);
            hid_t dataSet = openOrCreateDataSet<1>(file_,
                                                   dataSetName.c_str(),
                                                   "",
                                                   stringDataType,
                                                   chunkDims,
                                                   CompressionAlgorithm::LosslessNoShuffle,
                                                   0);
            writeData<1, true>(dataSet, propertyValuesChars, 0);
            H5Dclose(dataSet);
            sfree(propertyValuesChars);
        }
        else
        {
            /* Is there a more convenient way to do this? std::string is nice above, but cannot be used for writing. */
            std::vector<const char*> propertyValuesChars(propertyValues.size());
            std::transform(propertyValues.begin(),
                           propertyValues.end(),
                           propertyValuesChars.begin(),
                           std::mem_fn(&std::string::c_str));

            H5Tset_size(stringDataType, H5T_VARIABLE);
            H5Tset_strpad(stringDataType, H5T_STR_NULLTERM);
            hid_t dataSet = openOrCreateDataSet<1>(file_,
                                                   dataSetName.c_str(),
                                                   "",
                                                   stringDataType,
                                                   chunkDims,
                                                   CompressionAlgorithm::LosslessNoShuffle,
                                                   0);
            writeData<1, true>(dataSet, propertyValuesChars.data(), 0);
            H5Dclose(dataSet);
        }
    }
}

template<typename T>
void GmxH5mdIo::setNumericProperty(const std::string&    containerName,
                                   const std::string&    propertyName,
                                   const std::vector<T>& propertyValues,
                                   bool                  replaceExisting)
{
    openOrCreateGroup(file_, containerName.c_str());
    std::string dataSetName(containerName + "/" + propertyName);

    if (!H5Lexists(file_, dataSetName.c_str(), H5P_DEFAULT) || replaceExisting == true)
    {
        hid_t dataType;
        if constexpr (std::is_same<T, float>::value)
        {
            dataType = H5Tcopy(H5T_NATIVE_FLOAT);
        }
        else if constexpr (std::is_same<T, double>::value)
        {
            dataType = H5Tcopy(H5T_NATIVE_DOUBLE);
        }
        else if constexpr (std::is_same<T, int>::value)
        {
            dataType = H5Tcopy(H5T_NATIVE_INT);
        }
        else if constexpr (std::is_same<T, std::int64_t>::value)
        {
            dataType = H5Tcopy(H5T_NATIVE_INT64);
        }

        hsize_t atomPropertiesChunkDims[1];
        atomPropertiesChunkDims[0] = propertyValues.size();

        hid_t dataSet = openOrCreateDataSet<1>(file_,
                                               dataSetName.c_str(),
                                               "",
                                               dataType,
                                               atomPropertiesChunkDims,
                                               CompressionAlgorithm::LosslessNoShuffle,
                                               0);
        writeData<1, true>(dataSet, propertyValues.data(), 0);
        H5Dclose(dataSet);
    }
}

std::vector<std::string> GmxH5mdIo::readStringProperty(const std::string& containerName,
                                                       const std::string& propertyName)
{
    std::string              dataSetName(containerName + "/" + propertyName);
    hid_t                    dataSet = H5Dopen(file_, dataSetName.c_str(), H5P_DEFAULT);
    std::vector<std::string> propertyValues;

    if (dataSet < 0)
    {
        return propertyValues;
    }


    hid_t   origDatatype        = H5Dget_type(dataSet);
    hid_t   nativeDatatype      = H5Tget_native_type(origDatatype, H5T_DIR_DEFAULT);
    hsize_t dataTypeSize        = H5Tget_size(nativeDatatype);
    char*   propertyValuesChars = nullptr;
    size_t  totalNumElements, varStringLengthMaxLength;
    readData<1, true>(
            dataSet, 0, reinterpret_cast<void**>(&propertyValuesChars), &totalNumElements, &varStringLengthMaxLength);
    propertyValues.reserve(totalNumElements);
    dataTypeSize = varStringLengthMaxLength != 0 ? varStringLengthMaxLength : dataTypeSize;
    for (size_t i = 0; i < totalNumElements; i++)
    {
        propertyValues.push_back(propertyValuesChars + i * dataTypeSize);
    }
    free(propertyValuesChars);

    return propertyValues;
}

template<typename T>
std::vector<T> GmxH5mdIo::readNumericProperty(const std::string& containerName, const std::string& propertyName)
{
    std::string    dataSetName(containerName + "/" + propertyName);
    hid_t          dataSet = H5Dopen(file_, dataSetName.c_str(), H5P_DEFAULT);
    std::vector<T> propertyValues;

    if (dataSet < 0)
    {
        return propertyValues;
    }

    size_t totalNumElements, dummy;
    void*  buffer = nullptr;
    readData<1, true>(dataSet, 0, &buffer, &totalNumElements, &dummy);
    propertyValues.reserve(totalNumElements);

    hid_t dataType       = H5Dget_type(dataSet);
    hid_t nativeDataType = H5Tget_native_type(dataType, H5T_DIR_DEFAULT);

    if (H5Tequal(nativeDataType, H5T_NATIVE_FLOAT))
    {
        for (size_t i = 0; i < totalNumElements; i++)
        {
            propertyValues.push_back(static_cast<float*>(buffer)[i]);
        }
    }
    else if (H5Tequal(nativeDataType, H5T_NATIVE_DOUBLE))
    {
        for (size_t i = 0; i < totalNumElements; i++)
        {
            propertyValues.push_back(static_cast<double*>(buffer)[i]);
        }
    }
    else if (H5Tequal(nativeDataType, H5T_NATIVE_INT))
    {
        for (size_t i = 0; i < totalNumElements; i++)
        {
            propertyValues.push_back(static_cast<int*>(buffer)[i]);
        }
    }
    else if (H5Tequal(nativeDataType, H5T_NATIVE_INT64))
    {
        for (size_t i = 0; i < totalNumElements; i++)
        {
            propertyValues.push_back(static_cast<std::int64_t*>(buffer)[i]);
        }
    }
    else
    {
        throw gmx::FileIOError("Unhandled numeric data type when reading.");
    }

    free(buffer);

    return propertyValues;
}

void GmxH5mdIo::writeDataFrame(int64_t              step,
                               real                 time,
                               std::string          dataBlockFullName,
                               int                  dataDimensionalityFirstDim,
                               int                  dataDimensionalitySecondDim,
                               const real*          data,
                               std::string          unit,
                               hsize_t              numberOfFramesPerChunk,
                               CompressionAlgorithm compressionAlgorithm,
                               double               lossyCompressionError)

{
    GMX_ASSERT(data != nullptr, "Needs valid data to write a data frame.");
    GMX_ASSERT(dataDimensionalityFirstDim > 0 && dataDimensionalitySecondDim > 0,
               "The data dimensionality must be at least 1 in both dimensions.");
    /* See if the data block (API container for time dependent data sets) exists, otherwise create it */
    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), dataBlockFullName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        std::size_t lastSeparatorPos = dataBlockFullName.find_last_of("/");
        std::string groupName        = dataBlockFullName.substr(0, lastSeparatorPos);
        std::string dataBlockName =
                dataBlockFullName.substr(lastSeparatorPos + 1, dataBlockFullName.length());
        hid_t group = openOrCreateGroup(file_, groupName.c_str());

#if GMX_DOUBLE
        const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
        const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif

        GmxH5mdTimeDataBlock dataBlock(group,
                                       dataBlockName,
                                       unit,
                                       numberOfFramesPerChunk,
                                       dataDimensionalityFirstDim,
                                       dataDimensionalitySecondDim,
                                       datatype,
                                       compressionAlgorithm,
                                       lossyCompressionError);
        dataBlocks_.emplace_back(dataBlock);
        foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), dataBlockFullName.c_str());
        if (foundDataBlock == dataBlocks_.end())
        {
            throw gmx::FileIOError("Error creating data block when writing frame.");
        }
        H5Gclose(group);
    }
    foundDataBlock->writeFrame(data, step, time);
}

bool GmxH5mdIo::readNextFrameOfDataBlock(std::string dataBlockFullName, real* data, int64_t stepToRead)
{
    for (auto& dataBlock : dataBlocks_)
    {
        if (dataBlock.fullName() == dataBlockFullName)
        {
            if (stepToRead < 0 || dataBlock.getStepOfNextReadingFrame() == stepToRead)
            {
                return dataBlock.readNextFrame(data);
            }
            return false;
        }
    }
    return false;
}

real GmxH5mdIo::getLossyCompressionErrorOfDataBlock(std::string dataBlockFullName)
{
    for (const auto& dataBlock : dataBlocks_)
    {
        if (dataBlock.fullName() == dataBlockFullName)
        {
            return dataBlock.getLossyCompressionError();
        }
    }
    return -1;
}

int64_t GmxH5mdIo::getNumberOfFrames(const std::string dataBlockName, std::string selectionName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->numberOfFrames();
}

int64_t GmxH5mdIo::getNumberOfParticles(const std::string dataBlockName, std::string selectionName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getNumParticles();
}

real GmxH5mdIo::getFirstTime(const std::string dataBlockName, std::string selectionName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

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

std::tuple<int64_t, real> GmxH5mdIo::getNextStepAndTimeToRead()
{
    int64_t minStepNextFrame = std::numeric_limits<int64_t>::max();
    real    minTime          = std::numeric_limits<real>::max();
    for (const auto& dataBlock : dataBlocks_)
    {
        int64_t frameStep = dataBlock.getStepOfNextReadingFrame();
        /* Discard data sets that had a higher time stamp if an earlier data point has been found. */
        if (frameStep >= 0 && frameStep < minStepNextFrame)
        {
            minStepNextFrame = frameStep;
            minTime          = dataBlock.getTimeOfFrame(dataBlock.readingFrameIndex());
        }
    }
    return std::tuple<int64_t, real>(minStepNextFrame, minTime);
}

real GmxH5mdIo::getFinalTime(const std::string dataBlockName, std::string selectionName)
{
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

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
        int64_t numFrames = dataBlock.numberOfFrames();
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

extern template void writeData<1, false>(hid_t, const void*, hsize_t);
extern template void writeData<1, true>(hid_t, const void*, hsize_t);

extern template void readData<1, true>(hid_t, hsize_t, void**, size_t*, size_t*);

extern template void setAttribute<int>(hid_t, const char*, int, hid_t);
extern template void setAttribute<int64_t>(hid_t, const char*, int64_t, hid_t);
extern template void setAttribute<float>(hid_t, const char*, float, hid_t);
extern template void setAttribute<double>(hid_t, const char*, double, hid_t);

template void GmxH5mdIo::setNumericProperty<float>(const std::string&,
                                                   const std::string&,
                                                   const std::vector<float>&,
                                                   bool);
template void GmxH5mdIo::setNumericProperty<double>(const std::string&,
                                                    const std::string&,
                                                    const std::vector<double>&,
                                                    bool);
template void GmxH5mdIo::setNumericProperty<int>(const std::string&,
                                                 const std::string&,
                                                 const std::vector<int>&,
                                                 bool);
template void GmxH5mdIo::setNumericProperty<std::int64_t>(const std::string&,
                                                          const std::string&,
                                                          const std::vector<std::int64_t>&,
                                                          bool);

template std::vector<float> GmxH5mdIo::readNumericProperty<float>(const std::string&, const std::string&);
template std::vector<double> GmxH5mdIo::readNumericProperty<double>(const std::string&, const std::string&);
template std::vector<int> GmxH5mdIo::readNumericProperty<int>(const std::string&, const std::string&);
template std::vector<std::int64_t> GmxH5mdIo::readNumericProperty<std::int64_t>(const std::string&,
                                                                                const std::string&);

} // namespace h5mdio

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
    bool separateSelection = false;
    if (index.ssize() > 0)
    {
        separateSelection = true;
    }
    else
    {
        /* FIXME: Should use int64_t. Needs changes in topology. */
        for (int i = 0; i < topology.natoms; i++)
        {
            if (getGroupType(topology.groups, SimulationAtomGroupType::CompressedPositionOutput, i) != 0)
            {
                separateSelection = true;
                break;
            }
        }
    }
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

h5mdio::hid_t addMoleculeType(h5mdio::GmxH5mdIo* file, const gmx_moltype_t& molType)
{
#if GMX_USE_HDF5
    std::string moleculesGroupName = h5mdio::s_gromacsTopologyGroupName + "/molecules";
    file->createGroup(moleculesGroupName);
    hid_t moleculeTypeGroup = file->createGroup(moleculesGroupName + "/" + (*molType.name));

    h5mdio::setAttribute(
            moleculeTypeGroup, "number_of_atoms", static_cast<int64_t>(molType.atoms.nr), H5T_NATIVE_INT64);

    hid_t stringDataType = H5Tcopy(H5T_C_S1);
    H5Tset_cset(stringDataType, H5T_CSET_UTF8);
    size_t maxNameStringLength = 17;
    H5Tset_size(stringDataType, maxNameStringLength);
    hsize_t chunkDims[1];
    chunkDims[0] = molType.atoms.nr;

    hid_t atomNameDataSet = h5mdio::openOrCreateDataSet<1>(moleculeTypeGroup,
                                                           "atom_name",
                                                           "",
                                                           stringDataType,
                                                           chunkDims,
                                                           h5mdio::CompressionAlgorithm::LosslessNoShuffle,
                                                           0);
    hid_t residueNameDataSet =
            h5mdio::openOrCreateDataSet<1>(moleculeTypeGroup,
                                           "residue_name",
                                           "",
                                           stringDataType,
                                           chunkDims,
                                           h5mdio::CompressionAlgorithm::LosslessNoShuffle,
                                           0);
    for (ssize_t i = 0; i < molType.atoms.nr; i++)
    {
        int residueIndex = molType.atoms.atom[i].resind;

        h5mdio::writeData<1, false>(atomNameDataSet, *(molType.atoms.atomname[i]), i);
        h5mdio::writeData<1, false>(residueNameDataSet, *(molType.atoms.resinfo[residueIndex].name), i);
    }

    return moleculeTypeGroup;
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void addBlockOfMoleculeType(hid_t  moleculeTypeGroup,
                            size_t molBlockIndex,
                            size_t moleculeIndexStart,
                            size_t numMol,
                            size_t globalAtomStart)
{

    hid_t   dataType                  = H5Tcopy(H5T_NATIVE_INT64);
    hsize_t chunkDims[1]              = { 1 };
    hid_t   moleculeIndexStartDataSet = h5mdio::openOrCreateDataSet<1>(
            moleculeTypeGroup, "molecule_index_start", "", dataType, chunkDims, h5mdio::CompressionAlgorithm::None, 0);
    h5mdio::writeData<1, false>(moleculeIndexStartDataSet, &moleculeIndexStart, molBlockIndex);

    hid_t numMolStartDataSet = h5mdio::openOrCreateDataSet<1>(
            moleculeTypeGroup, "number_of_molecules", "", dataType, chunkDims, h5mdio::CompressionAlgorithm::None, 0);
    h5mdio::writeData<1, false>(numMolStartDataSet, &numMol, molBlockIndex);

    hid_t globalAtomsStartDataSet = h5mdio::openOrCreateDataSet<1>(moleculeTypeGroup,
                                                                   "global_atoms_start_index",
                                                                   "",
                                                                   dataType,
                                                                   chunkDims,
                                                                   h5mdio::CompressionAlgorithm::None,
                                                                   0);
    h5mdio::writeData<1, false>(globalAtomsStartDataSet, &globalAtomStart, molBlockIndex);
}


void setupMolecularSystemTopology(h5mdio::GmxH5mdIo* file, const gmx_mtop_t& topology, bool abortIfPresent)
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

    for (size_t i = 0; i < numMolBlocks; i++)
    {
        const gmx_molblock_t&       molBlock      = topology.molblock[i];
        const MoleculeBlockIndices& molBlockIndex = topology.moleculeBlockIndices[i];
        const gmx_moltype_t&        molType       = topology.moltype[molBlock.type];
        const std::string           molName       = *molType.name;
        const size_t                numMol        = molBlock.nmol;
        hid_t                       molTypeId     = addMoleculeType(file, molType);
        addBlockOfMoleculeType(
                molTypeId, i, molBlockIndex.moleculeIndexStart, numMol, molBlockIndex.globalAtomStart);
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
}

} // namespace gmx
