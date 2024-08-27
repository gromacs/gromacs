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
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 */

#include "gmxpre.h"

#include "h5md.h"

#include "config.h"

#include <strings.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"

#include "h5md_low_level_util.h"
#include "h5md_time_datablock.h"

#if GMX_USE_HDF5
CLANG_DIAGNOSTIC_IGNORE("-Wold-style-cast")
#    include <hdf5.h>
#else
CLANG_DIAGNOSTIC_IGNORE("-Wmissing-noreturn")
#endif

#if GMX_USE_HDF5
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
#    if H5_VERS_MINOR < 12
    H5Oget_info_by_name(locationId, name, &infoBuffer, H5P_DEFAULT);
#    else
    H5Oget_info_by_name(locationId, name, &infoBuffer, H5O_INFO_BASIC, H5P_DEFAULT);
#    endif
    herr_t            returnVal        = 0;
    const std::string stepDataSetName  = std::string(name) + std::string("/step");
    const std::string timeDataSetName  = std::string(name) + std::string("/time");
    const std::string valueDataSetName = std::string(name) + std::string("/value");
    switch (infoBuffer.type)
    {
        case H5O_TYPE_GROUP:
            if (gmx::objectExists(locationId, stepDataSetName.c_str())
                && gmx::objectExists(locationId, timeDataSetName.c_str())
                && gmx::objectExists(locationId, valueDataSetName.c_str()))
            {
                char containerFullName[gmx::c_maxFullNameLength];
                H5Iget_name(locationId, containerFullName, gmx::c_maxFullNameLength - 1);
                gmx::GmxH5mdTimeDataBlock             dataBlock(locationId, name);
                std::list<gmx::GmxH5mdTimeDataBlock>* dataBlocks =
                        static_cast<std::list<gmx::GmxH5mdTimeDataBlock>*>(operatorData);

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

void setPluginPath()
{
    const gmx::IProgramContext&       programContext = gmx::getProgramContext();
    const gmx::InstallationPrefixInfo installPrefix  = programContext.installationPrefix();
    if (!installPrefix.sourceLayoutTreeLike_)
    {
        std::filesystem::path pluginDir = installPrefix.path_;
        pluginDir.append("lib");
        if (H5PLprepend(pluginDir.string().c_str()) < 0)
        {
            printf("Cannot set the HDF5 plugin path. Please set the HDF5_PLUGIN_PATH environment "
                   "variable manually.\n");
        }
    }
}

} // namespace
#endif // GMX_USE_HDF5

namespace gmx
{

H5md::H5md(const std::filesystem::path& fileName, const char mode)
{
#if GMX_USE_HDF5
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error output, e.g. when items are not found.

    if (debug)
    {
        fprintf(debug, "Opening H5MD file %s with mode %c\n", fileName.string().c_str(), mode);
    }
    if (mode == 'w' || mode == 'a')
    {
        bool fileExists = gmx_fexist(fileName);
        if (!fileExists || mode == 'w')
        {
            make_backup(fileName.string().c_str());
            hid_t createPropertyList = H5Pcreate(H5P_FILE_CREATE);
            file_ = H5Fcreate(fileName.string().c_str(), H5F_ACC_TRUNC, createPropertyList, H5P_DEFAULT);
            if (file_ == H5I_INVALID_HID)
            {
                throw gmx::FileIOError("Cannot create H5MD file.");
            }
        }
        else
        {
            file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        }
        /* Create H5MD group. It should already be there if appending to a valid H5MD file, but it's better to be on the safe side. */
        hid_t h5mdGroup = openOrCreateGroup(file_, "h5md");
        setVersionAttribute(h5mdGroup, c_h5mdMajorVersion, c_h5mdMinorVersion);
    }
    else
    {
        file_ = H5Fopen(fileName.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    filemode_ = mode;

    setPluginPath();

    initGroupTimeDataBlocksFromFile("particles");
    initGroupTimeDataBlocksFromFile("observables");
    if (file_ == H5I_INVALID_HID)
    {
        throw gmx::FileIOError("Cannot open H5MD file.");
    }

#else
    GMX_UNUSED_VALUE(fileName);
    GMX_UNUSED_VALUE(mode);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

H5md::~H5md()
{
#if GMX_USE_HDF5
    if (file_ != H5I_INVALID_HID)
    {
        flush();
        if (debug)
        {
            fprintf(debug, "Closing H5MD file.\n");
        }
        for (auto dataBlock : dataBlocks_)
        {
            dataBlock.closeAllDataSets();
        }
        H5Fclose(file_);
    }

    /* Do not throw, if GMX_USE_HDF5 is false, in the destructor. */

#endif
}

void H5md::flush()
{
#if GMX_USE_HDF5
    if (file_ != H5I_INVALID_HID)
    {
        if (debug)
        {
            fprintf(debug, "Flushing H5MD file.\n");
        }
        if (filemode_ == 'w' || filemode_ == 'a')
        {
            addToProvenanceRecord();
        }
        if (H5Fflush(file_, H5F_SCOPE_LOCAL) < 0)
        {
            throw gmx::FileIOError("Error flushing H5MD.");
        }
    }

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

int H5md::initGroupTimeDataBlocksFromFile(const std::string& groupName)
{
#if GMX_USE_HDF5
    int   numDataBlocksBefore = dataBlocks_.size();
    hid_t group               = H5Gopen(file_, groupName.c_str(), H5P_DEFAULT);
    if (group == H5I_INVALID_HID)
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
        throw gmx::FileIOError("Error iterating over particles data blocks.");
    }
    return dataBlocks_.size() - numDataBlocksBefore;

#else
    GMX_UNUSED_VALUE(groupName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::string H5md::getH5mdRootVersionNumber()
{
#if GMX_USE_HDF5
    int   majorVersion, minorVersion;
    hid_t h5mdGroup = H5Gopen(file_, "h5md", H5P_DEFAULT);

    if (getVersionAttribute(h5mdGroup, &majorVersion, &minorVersion))
    {
        return std::to_string(majorVersion) + "." + std::to_string(minorVersion);
    }
    return "";

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::setAuthor(const std::string& authorName)
{
#if GMX_USE_HDF5
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    setAttribute(authorGroup, "name", authorName.c_str());

#else
    GMX_UNUSED_VALUE(authorName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::string H5md::getAuthor()
{
#if GMX_USE_HDF5
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    char* tmpName     = nullptr;
    getAttribute(authorGroup, "name", &tmpName);
    std::string name(tmpName);
    free(tmpName);
    return name;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::setCreatorProgramName(const std::string& creatorName)
{
#if GMX_USE_HDF5
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "name", creatorName.c_str());
#else
    GMX_UNUSED_VALUE(creatorName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::string H5md::getCreatorProgramName()
{
#if GMX_USE_HDF5
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    char* tmpName      = nullptr;
    getAttribute(creatorGroup, "name", &tmpName);
    std::string name(tmpName);
    free(tmpName);
    return name;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::setCreatorProgramVersion(const std::string& version)
{
#if GMX_USE_HDF5
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "version", version.c_str());

#else
    GMX_UNUSED_VALUE(version);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::string H5md::getCreatorProgramVersion()
{
#if GMX_USE_HDF5
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    char* tmpVersion   = nullptr;
    getAttribute(creatorGroup, "version", &tmpVersion);
    std::string version(tmpVersion);
    free(tmpVersion);
    return version;

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

hid_t H5md::getGroupId(const std::string& fullName)
{
#if GMX_USE_HDF5
    hid_t group = H5Gopen(file_, fullName.c_str(), H5P_DEFAULT);

    return group;

#else
    GMX_UNUSED_VALUE(fullName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

hid_t H5md::createGroup(const std::string& fullName)
{
#if GMX_USE_HDF5
    hid_t group = openOrCreateGroup(file_, fullName.c_str());

    return group;
#else
    GMX_UNUSED_VALUE(fullName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

hid_t H5md::createGroup(hid_t container, const std::string& nameInContainer)
{
#if GMX_USE_HDF5
    hid_t group = openOrCreateGroup(container, nameInContainer.c_str());

    return group;

#else
    GMX_UNUSED_VALUE(container);
    GMX_UNUSED_VALUE(nameInContainer);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::setStringDataSet(const std::string&              containerName,
                            const std::string&              dataSetName,
                            const std::vector<std::string>& propertyValues,
                            bool                            replaceExisting,
                            size_t                          maxStringLength)
{
#if GMX_USE_HDF5
    openOrCreateGroup(file_, containerName.c_str());
    std::string fullDataSetName(containerName + "/" + dataSetName);

    if (!H5Lexists(file_, fullDataSetName.c_str(), H5P_DEFAULT) || replaceExisting == true)
    {
        hid_t stringDataType = H5Tcopy(H5T_C_S1);
        H5Tset_cset(stringDataType, H5T_CSET_UTF8);
        hsize_t chunkDims[1];
        chunkDims[0] = propertyValues.size();
        if (maxStringLength > 0)
        {
            char* propertyValuesChars = new char[propertyValues.size() * maxStringLength];
            for (size_t i = 0; i < propertyValues.size(); i++)
            {
                strncpy(&propertyValuesChars[i * maxStringLength], propertyValues[i].c_str(), maxStringLength);
            }

            H5Tset_size(stringDataType, maxStringLength);
            hid_t dataSet = openOrCreateDataSet<1>(file_,
                                                   fullDataSetName.c_str(),
                                                   nullptr,
                                                   stringDataType,
                                                   chunkDims,
                                                   CompressionAlgorithm::LosslessNoShuffle,
                                                   0);
            writeData<1, true>(dataSet, propertyValuesChars, 0);
            H5Dclose(dataSet);
            delete[] propertyValuesChars;
        }
        else
        {
            std::vector<const char*> propertyValuesChars(propertyValues.size());
            std::transform(propertyValues.begin(),
                           propertyValues.end(),
                           propertyValuesChars.begin(),
                           std::mem_fn(&std::string::c_str));

            H5Tset_size(stringDataType, H5T_VARIABLE);
            H5Tset_strpad(stringDataType, H5T_STR_NULLTERM);
            hid_t dataSet = openOrCreateDataSet<1>(file_,
                                                   fullDataSetName.c_str(),
                                                   nullptr,
                                                   stringDataType,
                                                   chunkDims,
                                                   CompressionAlgorithm::LosslessNoShuffle,
                                                   0);
            writeData<1, true>(dataSet, propertyValuesChars.data(), 0);
            H5Dclose(dataSet);
        }
    }

#else
    GMX_UNUSED_VALUE(containerName);
    GMX_UNUSED_VALUE(dataSetName);
    GMX_UNUSED_VALUE(propertyValues);
    GMX_UNUSED_VALUE(replaceExisting);
    GMX_UNUSED_VALUE(maxStringLength);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

template<typename T>
void H5md::setNumericDataSet(const std::string&    containerName,
                             const std::string&    dataSetName,
                             const std::vector<T>& propertyValues,
                             const std::string&    unit,
                             bool                  replaceExisting)
{
#if GMX_USE_HDF5
    openOrCreateGroup(file_, containerName.c_str());
    std::string fullDataSetName(containerName + "/" + dataSetName);

    if (!H5Lexists(file_, fullDataSetName.c_str(), H5P_DEFAULT) || replaceExisting == true)
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
        hid_t dataSet;
        if constexpr (std::is_same<T, std::pair<std::int64_t, std::int64_t>>::value)
        {
            dataType                           = H5Tcopy(H5T_NATIVE_INT64);
            hsize_t atomPropertiesChunkDims[2] = { propertyValues.size(), 2 };

            dataSet = openOrCreateDataSet<2>(file_,
                                             fullDataSetName.c_str(),
                                             unit.empty() ? nullptr : unit.c_str(),
                                             dataType,
                                             atomPropertiesChunkDims,
                                             CompressionAlgorithm::LosslessNoShuffle,
                                             0);
            writeData<2, true>(dataSet, propertyValues.data(), 0);
        }
        else
        {
            hsize_t atomPropertiesChunkDims[1];
            atomPropertiesChunkDims[0] = propertyValues.size();

            dataSet = openOrCreateDataSet<1>(file_,
                                             fullDataSetName.c_str(),
                                             unit.empty() ? nullptr : unit.c_str(),
                                             dataType,
                                             atomPropertiesChunkDims,
                                             CompressionAlgorithm::LosslessNoShuffle,
                                             0);
            writeData<1, true>(dataSet, propertyValues.data(), 0);
        }
        H5Dclose(dataSet);
    }

#else
    GMX_UNUSED_VALUE(containerName);
    GMX_UNUSED_VALUE(dataSetName);
    GMX_UNUSED_VALUE(propertyValues);
    GMX_UNUSED_VALUE(unit);
    GMX_UNUSED_VALUE(replaceExisting);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::vector<std::string> H5md::readStringDataSet(const std::string& containerName, const std::string& dataSetName)
{
#if GMX_USE_HDF5
    std::string              fullDataSetName(containerName + "/" + dataSetName);
    hid_t                    dataSet = H5Dopen(file_, fullDataSetName.c_str(), H5P_DEFAULT);
    std::vector<std::string> propertyValues;

    if (dataSet == H5I_INVALID_HID)
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

#else
    GMX_UNUSED_VALUE(containerName);
    GMX_UNUSED_VALUE(dataSetName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

template<typename T>
std::vector<T> H5md::readNumericDataSet(const std::string& containerName, const std::string& dataSetName)
{
#if GMX_USE_HDF5
    std::string    fullDataSetName(containerName + "/" + dataSetName);
    hid_t          dataSet = H5Dopen(file_, fullDataSetName.c_str(), H5P_DEFAULT);
    std::vector<T> propertyValues;

    if (dataSet == H5I_INVALID_HID)
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

#else
    GMX_UNUSED_VALUE(containerName);
    GMX_UNUSED_VALUE(dataSetName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::writeDataFrame(int64_t              step,
                          real                 time,
                          const std::string&   dataBlockFullName,
                          int                  dataDimensionalityFirstDim,
                          int                  dataDimensionalitySecondDim,
                          const real*          data,
                          const std::string&   unit,
                          hsize_t              numberOfFramesPerChunk,
                          CompressionAlgorithm compressionAlgorithm,
                          double               lossyCompressionError)

{
#if GMX_USE_HDF5
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

#    if GMX_DOUBLE
        const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#    else
        const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#    endif

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

#else
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(time);
    GMX_UNUSED_VALUE(dataBlockFullName);
    GMX_UNUSED_VALUE(dataDimensionalityFirstDim);
    GMX_UNUSED_VALUE(dataDimensionalitySecondDim);
    GMX_UNUSED_VALUE(data);
    GMX_UNUSED_VALUE(unit);
    GMX_UNUSED_VALUE(numberOfFramesPerChunk);
    GMX_UNUSED_VALUE(compressionAlgorithm);
    GMX_UNUSED_VALUE(lossyCompressionError);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

bool H5md::readNextFrameOfDataBlock(const std::string& dataBlockFullName, real* data, int64_t stepToRead)
{
#if GMX_USE_HDF5
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

#else
    GMX_UNUSED_VALUE(dataBlockFullName);
    GMX_UNUSED_VALUE(data);
    GMX_UNUSED_VALUE(stepToRead);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

double H5md::getLossyCompressionErrorOfDataBlock(const std::string& dataBlockFullName)
{
#if GMX_USE_HDF5
    for (const auto& dataBlock : dataBlocks_)
    {
        if (dataBlock.fullName() == dataBlockFullName)
        {
            return dataBlock.getLossyCompressionError();
        }
    }
    return -1;

#else
    GMX_UNUSED_VALUE(dataBlockFullName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

int64_t H5md::getNumberOfFrames(const std::string& dataBlockName, const std::string& selectionName)
{
#if GMX_USE_HDF5
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->numberOfFrames();

#else
    GMX_UNUSED_VALUE(dataBlockName);
    GMX_UNUSED_VALUE(selectionName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

int64_t H5md::getNumberOfParticles(const std::string& dataBlockName, const std::string& selectionName)
{
#if GMX_USE_HDF5
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getNumParticles();

#else
    GMX_UNUSED_VALUE(dataBlockName);
    GMX_UNUSED_VALUE(selectionName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

real H5md::getFirstTime(const std::string& dataBlockName, const std::string& selectionName)
{
#if GMX_USE_HDF5
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getTimeOfFrame(0);

#else
    GMX_UNUSED_VALUE(dataBlockName);
    GMX_UNUSED_VALUE(selectionName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

real H5md::getFirstTimeFromAllDataBlocks()
{
#if GMX_USE_HDF5
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

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

std::tuple<int64_t, real> H5md::getNextStepAndTimeToRead()
{
#if GMX_USE_HDF5
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

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

real H5md::getFinalTime(const std::string& dataBlockName, const std::string& selectionName)
{
#if GMX_USE_HDF5
    GMX_ASSERT(dataBlockName != "", "There must be a datablock name to look for.");

    std::string wantedName = "/particles/" + selectionName + "/" + dataBlockName;

    auto foundDataBlock = std::find(dataBlocks_.begin(), dataBlocks_.end(), wantedName.c_str());
    if (foundDataBlock == dataBlocks_.end())
    {
        return -1;
    }
    return foundDataBlock->getTimeOfFrame(foundDataBlock->numberOfFrames() - 1);

#else
    GMX_UNUSED_VALUE(dataBlockName);
    GMX_UNUSED_VALUE(selectionName);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

real H5md::getFinalTimeFromAllDataBlocks()
{
#if GMX_USE_HDF5
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

#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void H5md::addToProvenanceRecord(const std::string& commandLine,
                                 const std::string& programVersion,
                                 const std::string& comment)
{
#if GMX_USE_HDF5
    hid_t provenanceGroup = createGroup(s_provenanceGroupName);
    setVersionAttribute(
            provenanceGroup, c_gmxH5mdProvenanceGroupMajorVersion, c_gmxH5mdProvenanceGroupMinorVersion);

    hid_t stringDataType = H5Tcopy(H5T_C_S1);
    H5Tset_cset(stringDataType, H5T_CSET_UTF8);
    size_t recordStringLength = c_provenanceRecordStringLen;
    H5Tset_size(stringDataType, recordStringLength);
    hsize_t chunkDims[1] = { 1 };

    /* When creating a new data set the number of frames is 1 (there is a first empty record).
     * Therefore handle the number of frames (used to specify what records to write) differently
     * if the data set is created or if it already exists (and should have data entries). */
    hsize_t numFrames          = 0;
    hid_t   commandLineDataSet = H5Dopen(provenanceGroup, "command_line", H5P_DEFAULT);
    if (commandLineDataSet == H5I_INVALID_HID)
    {
        commandLineDataSet = openOrCreateDataSet<1>(provenanceGroup,
                                                    "command_line",
                                                    nullptr,
                                                    stringDataType,
                                                    chunkDims,
                                                    CompressionAlgorithm::LosslessNoShuffle,
                                                    0);
    }
    else
    {

        hid_t dataSpace = H5Dget_space(commandLineDataSet);
        if (dataSpace == H5I_INVALID_HID)
        {
            throw gmx::FileIOError("The main data block of the provenance record cannot be found.");
        }
        H5Sget_simple_extent_dims(dataSpace, &numFrames, nullptr);
    }

    char tmpString[c_provenanceRecordStringLen];
    snprintf(tmpString,
             c_provenanceRecordStringLen - 1,
             "%s",
             commandLine.empty() ? gmx::getProgramContext().commandLine() : commandLine.c_str());
    writeData<1, false>(commandLineDataSet, tmpString, numFrames);

    hid_t programVersionDataSet = openOrCreateDataSet<1>(provenanceGroup,
                                                         "program_version",
                                                         nullptr,
                                                         stringDataType,
                                                         chunkDims,
                                                         CompressionAlgorithm::LosslessNoShuffle,
                                                         0);
    snprintf(tmpString,
             c_provenanceRecordStringLen - 1,
             "%s",
             programVersion.empty() ? gmx_version() : programVersion.c_str());
    writeData<1, false>(programVersionDataSet, tmpString, numFrames);

    hid_t dataType    = H5Tcopy(H5T_NATIVE_INT64);
    hid_t timeDataSet = openOrCreateDataSet<1>(
            provenanceGroup, "time", "s", dataType, chunkDims, CompressionAlgorithm::LosslessNoShuffle, 0);
    const int64_t timeStamp = std::time(nullptr);
    writeData<1, false>(timeDataSet, &timeStamp, numFrames);

    hid_t commentDataSet = openOrCreateDataSet<1>(
            provenanceGroup, "comment", nullptr, stringDataType, chunkDims, CompressionAlgorithm::LosslessNoShuffle, 0);
    snprintf(tmpString, c_provenanceRecordStringLen - 1, "%s", comment.c_str());
    writeData<1, false>(commentDataSet, tmpString, numFrames);

#else
    GMX_UNUSED_VALUE(commandLine);
    GMX_UNUSED_VALUE(programVersion);
    GMX_UNUSED_VALUE(comment);
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}


extern template hid_t
openOrCreateDataSet<1>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);
extern template hid_t
openOrCreateDataSet<2>(hid_t, const char*, const char*, hid_t, const hsize_t*, CompressionAlgorithm, double);

extern template void writeData<1, false>(hid_t, const void*, hsize_t);
extern template void writeData<1, true>(hid_t, const void*, hsize_t);
extern template void writeData<2, false>(hid_t, const void*, hsize_t);
extern template void writeData<2, true>(hid_t, const void*, hsize_t);

extern template void readData<1, true>(hid_t, hsize_t, void**, size_t*, size_t*);
extern template void readData<1, false>(hid_t, hsize_t, void**, size_t*, size_t*);
extern template void readData<1>(hid_t, hsize_t, void**);

extern template void setAttribute<int>(hid_t, const char*, int, hid_t);
extern template void setAttribute<int64_t>(hid_t, const char*, int64_t, hid_t);
extern template void setAttribute<float>(hid_t, const char*, float, hid_t);
extern template void setAttribute<double>(hid_t, const char*, double, hid_t);

extern template bool getAttribute<int64_t>(hid_t, const char*, int64_t*);

template void H5md::setNumericDataSet<float>(const std::string&,
                                             const std::string&,
                                             const std::vector<float>&,
                                             const std::string&,
                                             bool);
template void H5md::setNumericDataSet<double>(const std::string&,
                                              const std::string&,
                                              const std::vector<double>&,
                                              const std::string&,
                                              bool);
template void H5md::setNumericDataSet<int>(const std::string&,
                                           const std::string&,
                                           const std::vector<int>&,
                                           const std::string&,
                                           bool);
template void H5md::setNumericDataSet<std::int64_t>(const std::string&,
                                                    const std::string&,
                                                    const std::vector<std::int64_t>&,
                                                    const std::string&,
                                                    bool);
template void H5md::setNumericDataSet<std::pair<std::int64_t, std::int64_t>>(
        const std::string&,
        const std::string&,
        const std::vector<std::pair<std::int64_t, std::int64_t>>&,
        const std::string&,
        bool);

template std::vector<float> H5md::readNumericDataSet<float>(const std::string&, const std::string&);
template std::vector<double> H5md::readNumericDataSet<double>(const std::string&, const std::string&);
template std::vector<int>    H5md::readNumericDataSet<int>(const std::string&, const std::string&);
template std::vector<std::int64_t> H5md::readNumericDataSet<std::int64_t>(const std::string&,
                                                                          const std::string&);

} // namespace gmx

CLANG_DIAGNOSTIC_RESET
