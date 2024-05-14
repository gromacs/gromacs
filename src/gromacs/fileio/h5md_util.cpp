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

#include "h5md_util.h"

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#    include <hdf5.h>

#    include "external/SZ3-bio/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif

namespace
{

/*! Set the fill value to -1 in a data set property list.
 * \param[in] dataSetCreatePropertyList The ID of the propery list to update.
 * \param[in] dataType The ID of the HDF5 data type of the data set.
 */
void setNumericFillValue(const hid_t dataSetCreatePropertyList, const hid_t dataType)
{
    if (H5Tequal(dataType, H5T_NATIVE_INT))
    {
        const int dataFill = -1;
        H5Pset_fill_value(dataSetCreatePropertyList, dataType, &dataFill);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_INT64))
    {
        const int64_t dataFill = -1;
        H5Pset_fill_value(dataSetCreatePropertyList, dataType, &dataFill);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_FLOAT))
    {
        const float dataFill = -1;
        H5Pset_fill_value(dataSetCreatePropertyList, dataType, &dataFill);
    }
    else if (H5Tequal(dataType, H5T_NATIVE_DOUBLE))
    {
        const double dataFill = -1;
        H5Pset_fill_value(dataSetCreatePropertyList, dataType, &dataFill);
    }
    /* No fill value for other data types */
}

/*! \brief Set the filter settings (of an SZ3 filter) in a property list
 * \param[in] propertyList The ID of the property list to update
 * \param[in] sz3Mode The compression error mode, 0 == absolute, 1 == relative
 * \param[in] compressionError The maximum error when compressing (compressionError = accuracy / 2)
 */
void setLossySz3CompressionProperties(const hid_t propertyList, const int sz3Mode, const double compressionError)
{
    size_t        numCompressionSettingsElements;
    unsigned int* compressionSettings = nullptr;
    SZ_errConfigToCdArray(
            &numCompressionSettingsElements, &compressionSettings, sz3Mode, compressionError, compressionError, 0, 0);
    herr_t status = H5Pset_filter(
            propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings);
    free(compressionSettings);
    if (status < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot set SZ3 compression.");
    }
    if (H5Zfilter_avail(H5Z_FILTER_SZ3) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("SZ3 filter not available.");
    }
}

} // namespace

namespace gmx
{
namespace h5mdio
{

hid_t openOrCreateGroup(const hid_t container, const char* name)
{
#if GMX_USE_HDF5
    hid_t group = H5Gopen(container, name, H5P_DEFAULT);
    if (group < 0)
    {
        hid_t linkPropertyList = H5Pcreate(H5P_LINK_CREATE); // create group creation property list
        if (linkPropertyList < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Cannot create linkPropertyList when creating group.");
        }
        H5Pset_create_intermediate_group(linkPropertyList, 1); // set intermediate link creation
        group = H5Gcreate(container, name, linkPropertyList, H5P_DEFAULT, H5P_DEFAULT);
        if (group < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Cannot create group.");
        }
    }
    return group;
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void registerSz3FilterImplicitly()
{
#if GMX_USE_HDF5
    hid_t         propertyList = H5Pcreate(H5P_DATASET_CREATE);
    int           sz3Mode      = 0; // 0: ABS, 1: REL
    size_t        numCompressionSettingsElements;
    unsigned int* compressionSettings = nullptr;
    SZ_errConfigToCdArray(&numCompressionSettingsElements, &compressionSettings, sz3Mode, 0, 0, 0, 0);
    herr_t status = H5Pset_filter(
            propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings);
    free(compressionSettings);

    /* For some reason status is 0 even if the filter could not be found. Check if any HDF5 errors have occured */
    ssize_t numHdf5Errors = H5Eget_num(H5E_DEFAULT);
    if (status < 0 || numHdf5Errors > 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError(
                "Cannot use SZ3 compression filter. Please check that the SZ3 filter is in "
                "HDF5_PLUGIN_PATH.");
    }
#endif
}

template<int numDims>
hid_t openOrCreateDataSet(const hid_t                container,
                          const char*                name,
                          const char*                unit,
                          const hid_t                dataType,
                          const hsize_t*             chunkDims,
                          const CompressionAlgorithm compression,
                          const double               compressionError)


{
    hid_t dataSet = H5Dopen(container, name, H5P_DEFAULT);

    if (dataSet < 0)
    {
        hsize_t maxDims[numDims];
        maxDims[0] = H5S_UNLIMITED;
        for (int i = 1; i < numDims; i++)
        {
            maxDims[i] = chunkDims[i];
        }
        hid_t dataSpace          = H5Screate_simple(numDims, chunkDims, maxDims);
        hid_t createPropertyList = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(createPropertyList, numDims, chunkDims);
        setNumericFillValue(createPropertyList, dataType);

        /* It would be nice to have an option not to write full incomplete edge chunks,
         * but the closest option is:
         * H5Pset_chunk_opts(createPropertyList, H5D_CHUNK_DONT_FILTER_PARTIAL_CHUNKS);
         * But that only avoids compressing/decompressing the edge chunks.
         * Keep an eye open for alternatives.
         * TODO: It is possible that it would be time-efficient
         * to avoid compressing edge chunks when writing checkpoints. Pros and cons for slightly
         * larger files vs slightly faster checkpoint writing must be evaluated.
         * Currently it seems like incomplete edge chunks are compressed even with this option. */

        switch (compression)
        {
            case CompressionAlgorithm::LossySz3:
            {
                const int sz3Mode = 0; // 0: ABS, 1: REL
                setLossySz3CompressionProperties(createPropertyList, sz3Mode, compressionError);
                break;
            }
            case CompressionAlgorithm::LosslessWithShuffle:
                if (H5Pset_shuffle(createPropertyList) < 0)
                {
                    throw gmx::FileIOError("Cannot set shuffle filter.");
                }
                /* Fall through */
            case CompressionAlgorithm::LosslessNoShuffle:
                if (H5Pset_deflate(createPropertyList, 1) < 0)
                {
                    H5Eprint2(H5E_DEFAULT, nullptr);
                    throw gmx::FileIOError("Cannot set GZIP compression.");
                }
                break;
            case CompressionAlgorithm::None: break;
            default: throw gmx::FileIOError("Unrecognized compression mode.");
        }
        /* Set a reasonable cache based on chunk sizes. The cache is not stored in file, so must be set when opening a dataset */
        size_t cacheSize = sizeof(real);
        for (int i = 0; i < numDims; i++)
        {
            cacheSize *= chunkDims[i];
        }
        hid_t accessPropertyList = H5Pcreate(H5P_DATASET_ACCESS);
        H5Pset_chunk_cache(
                accessPropertyList, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, cacheSize, H5D_CHUNK_CACHE_W0_DEFAULT);

        dataSet = H5Dcreate(
                container, name, dataType, dataSpace, H5P_DEFAULT, createPropertyList, accessPropertyList);
        if (dataSet < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Cannot create dataSet.");
        }

        if (unit != nullptr)
        {
            char unitElementString[] = "unit";
            setAttribute(dataSet, unitElementString, unit);
        }
    }
    return dataSet;
}

template<int numDims, bool writeFullDataSet>
void writeData(const hid_t dataSet, const void* data, const hsize_t frameToWrite)
{
#if GMX_USE_HDF5
    GMX_ASSERT(dataSet >= 0, "Needs a valid dataSet to write data.");
    GMX_ASSERT(data != nullptr_t, "Needs valid data to write.");
    GMX_ASSERT(!writeFullDataSet || frameToWrite == 0,
               "Must start writing from frame 0 if writing the whole data set.");

    hid_t   dataSpace          = H5Dget_space(dataSet);
    hid_t   createPropertyList = H5Dget_create_plist(dataSet);
    hsize_t currentDims[numDims], chunkDims[numDims];
    H5Sget_simple_extent_dims(dataSpace, currentDims, nullptr);
    H5Pget_chunk(createPropertyList, numDims, chunkDims);
    const hsize_t numFramesPerChunk = chunkDims[0];

    /* Resize the dataSet if needed. */
    if (frameToWrite >= currentDims[0])
    {
        hsize_t newDims[numDims];
        newDims[0] = (frameToWrite / numFramesPerChunk + 1) * numFramesPerChunk;
        for (int i = 1; i < numDims; i++)
        {
            newDims[i] = currentDims[i];
        }
        if (debug)
        {
            fprintf(debug, "Resizing dataSet from %llu to %llu\n", currentDims[0], newDims[0]);
        }
        H5Dset_extent(dataSet, newDims);
        dataSpace = H5Dget_space(dataSet);
    }
    hsize_t fileOffset[numDims];
    fileOffset[0] = frameToWrite;
    hsize_t outputBlockSize[numDims];
    if constexpr (writeFullDataSet)
    {
        outputBlockSize[0] = currentDims[0];
    }
    else
    {
        outputBlockSize[0] = 1;
    }

    for (int i = 1; i < numDims; i++)
    {
        fileOffset[i]      = 0;
        outputBlockSize[i] = currentDims[i];
    }
    if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, fileOffset, nullptr, outputBlockSize, nullptr) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot select the output region.");
    }

    hid_t memoryDataspace = H5Screate_simple(numDims, outputBlockSize, nullptr);
    hid_t dataType        = H5Dget_type(dataSet);
    if (H5Dwrite(dataSet, dataType, memoryDataspace, dataSpace, H5P_DEFAULT, data) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Error writing data.");
    }

    // It would be good to close the dataset here, but that means compressing and writing the whole chunk every time - very slow.
#else
    throw gmx::FileIOError(
            "GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

size_t getDataTypeSize(const hid_t dataSet)
{
    hid_t origDatatype   = H5Dget_type(dataSet);
    hid_t nativeDatatype = H5Tget_native_type(origDatatype, H5T_DIR_DEFAULT);

    return H5Tget_size(nativeDatatype);
}

template<int numDims, bool readFullDataSet>
void readData(const hid_t   dataSet,
              const hsize_t frameToRead,
              void**        buffer,
              size_t*       totalNumElements,
              size_t*       varLengthStringMaxLength)
{
    GMX_ASSERT(dataSet >= 0, "Needs a valid dataSet to read data.");
    GMX_ASSERT(!readFullDataSet || frameToRead == 0,
               "Must start reading from frame 0 if reading the whole data set.");

    hid_t dataSpace = H5Dget_space(dataSet);
    if (dataSpace < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError(
                "The main data block of the time dependent data set cannot be found.");
    }
    const int dataSpaceNumDims = H5Sget_simple_extent_ndims(dataSpace);
    if (numDims != dataSpaceNumDims)
    {
        throw gmx::FileIOError("The data set dimensions do not match what is expected.");
    }

    hsize_t dimExtents[numDims];
    H5Sget_simple_extent_dims(dataSpace, dimExtents, nullptr);
    hsize_t maxNumFrames = dimExtents[0];
    if (frameToRead >= maxNumFrames)
    {
        throw gmx::FileIOError("Trying to read outside the valid frame range.");
    }

    hsize_t fileOffset[numDims];
    fileOffset[0] = frameToRead;
    hsize_t inputBlockSize[numDims];
    if constexpr (readFullDataSet)
    {
        inputBlockSize[0] = dimExtents[0];
    }
    else
    {
        inputBlockSize[0] = 1;
    }
    size_t totalBlockSize = inputBlockSize[0];

    for (int i = 1; i < numDims; i++)
    {
        fileOffset[i]     = 0;
        inputBlockSize[i] = dimExtents[i];
        totalBlockSize *= inputBlockSize[i];
    }
    if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, fileOffset, nullptr, inputBlockSize, nullptr) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot select the input region.");
    }
    *totalNumElements = totalBlockSize;

    hid_t   memoryDataspace = H5Screate_simple(numDims, inputBlockSize, nullptr);
    hid_t   origDatatype    = H5Dget_type(dataSet);
    hid_t   nativeDatatype  = H5Tget_native_type(origDatatype, H5T_DIR_DEFAULT);
    hsize_t dataTypeSize    = H5Tget_size(nativeDatatype);
    if (H5Tis_variable_str(origDatatype))
    {
        char** tmpBuffer = static_cast<char**>(malloc(dataTypeSize * totalBlockSize));
        if (H5Dread(dataSet, nativeDatatype, memoryDataspace, dataSpace, H5P_DEFAULT, tmpBuffer) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error reading data set.");
        }
        size_t maxLength = 0;
        for (size_t i = 0; i < totalBlockSize; i++)
        {
            maxLength = std::max(strlen(tmpBuffer[i]), maxLength);
        }
        maxLength += 1;
        *varLengthStringMaxLength = maxLength;
        if (*buffer == nullptr)
        {
            *buffer = malloc(maxLength * totalBlockSize);
        }
        for (size_t i = 0; i < totalBlockSize; i++)
        {
            strncpy(static_cast<char*>(*buffer) + i * maxLength, tmpBuffer[i], maxLength);
        }
        H5Dvlen_reclaim(origDatatype, dataSpace, H5P_DEFAULT, tmpBuffer);
        free(tmpBuffer);
    }
    else
    {
        *varLengthStringMaxLength = 0;
        if (*buffer == nullptr)
        {
            *buffer = malloc(dataTypeSize * totalBlockSize);
        }
        if (H5Dread(dataSet, nativeDatatype, memoryDataspace, dataSpace, H5P_DEFAULT, *buffer) < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            throw gmx::FileIOError("Error reading data set.");
        }
    }
}

template<int numDims>
void readData(const hid_t dataSet, const hsize_t frameToRead, void** buffer)
{
    size_t totalNumElementsDummy, varLengthStringMaxLengthDummy;

    readData<numDims, false>(
            dataSet, frameToRead, buffer, &totalNumElementsDummy, &varLengthStringMaxLengthDummy);
}

void setBoxGroupAttributes(const hid_t boxGroup, const PbcType pbcType)
{
    setAttribute(boxGroup, "dimension", DIM, H5T_NATIVE_INT);
    static constexpr int c_pbcTypeStringLength                               = 9;
    char                 boundaryAttributeString[DIM][c_pbcTypeStringLength] = { "periodic",
                                                                 "periodic",
                                                                 "periodic" };
    switch (pbcType)
    {
        case PbcType::Xyz: break;
        case PbcType::XY: strcpy(boundaryAttributeString[2], "none"); break;
        default:
            for (int i = 0; i < DIM; i++)
            {
                strcpy(boundaryAttributeString[i], "none");
            }
            break;
    }
    setAttributeStringList<DIM, c_pbcTypeStringLength>(boxGroup, "boundary", boundaryAttributeString);
}

void setVersionAttribute(const hid_t group, const int majorVersion, const int minorVersion)
{
    char  name[]    = "version";
    hid_t attribute = H5Aopen(group, name, H5P_DEFAULT);
    hid_t dataType  = H5Tcopy(H5T_NATIVE_INT32);

    if (attribute < 0)
    {
        hsize_t dataSize[1] = { 2 };
        hid_t   dataSpace   = H5Screate_simple(1, dataSize, nullptr);
        attribute = H5Acreate2(group, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    }
    const int value[2] = { majorVersion, minorVersion };
    if (H5Awrite(attribute, dataType, &value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

bool getVersionAttribute(const hid_t group, int* majorVersion, int* minorVersion)
{
    char  name[]    = "version";
    hid_t attribute = H5Aopen(group, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        return false;
    }
    hid_t dataType = H5Aget_type(attribute);
    int   value[2];
    if (H5Aread(attribute, dataType, &value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot write attribute.");
    }
    *majorVersion = value[0];
    *minorVersion = value[1];

    H5Aclose(attribute);
    return true;
}

template<typename T>
void setAttribute(const hid_t dataSet, const char* name, const T value, const hid_t dataType)
{
    hid_t attribute = H5Aopen(dataSet, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hid_t dataSpace = H5Screate(H5S_SCALAR);
        attribute       = H5Acreate2(dataSet, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Awrite(attribute, dataType, &value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

void setAttribute(const hid_t dataSet, const char* name, const char* value)
{
    hid_t dataType = H5Tcopy(H5T_C_S1);
    H5Tset_size(dataType, strlen(value));
    H5Tset_strpad(dataType, H5T_STR_NULLTERM);
    H5Tset_cset(dataType, H5T_CSET_UTF8);

    hid_t attribute = H5Aopen(dataSet, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hid_t dataSpace = H5Screate(H5S_SCALAR);
        attribute       = H5Acreate2(dataSet, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Awrite(attribute, dataType, value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

template<typename T>
bool getAttribute(const hid_t dataSet, const char* name, T* value)
{
    hid_t attribute = H5Aopen(dataSet, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        return false;
    }
    hid_t dataType = H5Aget_type(attribute);
    if (H5Aread(attribute, dataType, value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot read attribute.");
    }

    H5Aclose(attribute);
    return true;
}

bool getAttribute(const hid_t dataSet, const char* name, char** value)
{
    if (!H5Aexists(dataSet, name))
    {
        return false;
    }
    hid_t attribute = H5Aopen(dataSet, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        return false;
    }
    hid_t dataType = H5Aget_type(attribute);
    /* Make room for string termination as well. */
    size_t allocationSize = H5Tget_size(dataType) + 1;
    *value                = reinterpret_cast<char*>(malloc(allocationSize));
    memset(*value, NULL, allocationSize);
    if (H5Aread(attribute, dataType, *value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot read attribute.");
    }

    H5Aclose(attribute);
    return true;
}

template<hid_t numEntries, size_t stringLength>
void setAttributeStringList(const hid_t dataSet, const char* name, const char value[numEntries][stringLength])
{
    hid_t dataType = H5Tcopy(H5T_C_S1);
    H5Tset_size(dataType, stringLength);
    H5Tset_strpad(dataType, H5T_STR_NULLTERM);
    H5Tset_cset(dataType, H5T_CSET_UTF8);
    hid_t attribute = H5Aopen(dataSet, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hsize_t dataSize[1] = { numEntries };
        hid_t   dataSpace   = H5Screate_simple(1, dataSize, nullptr);
        attribute = H5Acreate2(dataSet, name, dataType, dataSpace, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Awrite(attribute, dataType, &value[0]) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        throw gmx::FileIOError("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

real getDataSetSz3CompressionError(const hid_t dataSet)
{
    hid_t        propertyList                   = H5Dget_create_plist(dataSet);
    unsigned int flags                          = 0;
    size_t       numCompressionSettingsElements = 9;
    unsigned int compressionSettingsValues[16];
    real         compressionError = -1;
    if (H5Pget_filter_by_id(
                propertyList, H5Z_FILTER_SZ3, &flags, &numCompressionSettingsElements, compressionSettingsValues, 0, nullptr, nullptr)
        < 0)
    {
        return compressionError;
    }
    int    dimSize     = 0;
    int    dataType    = 0;
    int    errorMode   = 0;
    double absError    = 0;
    double relError    = 0;
    double l2normError = 0;
    double psNr        = 0;
    size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
    int withErrorInfo = checkCDValuesWithErrors(numCompressionSettingsElements, compressionSettingsValues);
    if (!withErrorInfo)
    {
        return -1;
    }
    SZ_cdArrayToMetaDataErr(numCompressionSettingsElements,
                            compressionSettingsValues,
                            &dimSize,
                            &dataType,
                            &r5,
                            &r4,
                            &r3,
                            &r2,
                            &r1,
                            &errorMode,
                            &absError,
                            &relError,
                            &l2normError,
                            &psNr);
    if (errorMode == 0)
    {
        return absError;
    }
    else if (errorMode == 1)
    {
        return relError;
    }
    return -1;
}

bool objectExists(const hid_t container, const char* name)
{
    return H5Lexists(container, name, H5P_DEFAULT) >= 0
           && H5Oexists_by_name(container, name, H5P_DEFAULT) >= 0;
}

} // namespace h5mdio
} // namespace gmx

template hid_t gmx::h5mdio::openOrCreateDataSet<1>(hid_t,
                                                   const char*,
                                                   const char*,
                                                   hid_t,
                                                   const hsize_t*,
                                                   gmx::h5mdio::CompressionAlgorithm,
                                                   double);
template hid_t gmx::h5mdio::openOrCreateDataSet<2>(hid_t,
                                                   const char*,
                                                   const char*,
                                                   hid_t,
                                                   const hsize_t*,
                                                   gmx::h5mdio::CompressionAlgorithm,
                                                   double);
template hid_t gmx::h5mdio::openOrCreateDataSet<3>(hid_t,
                                                   const char*,
                                                   const char*,
                                                   hid_t,
                                                   const hsize_t*,
                                                   gmx::h5mdio::CompressionAlgorithm,
                                                   double);

template void gmx::h5mdio::writeData<1, false>(hid_t, const void*, hsize_t);
template void gmx::h5mdio::writeData<1, true>(hid_t, const void*, hsize_t);
template void gmx::h5mdio::writeData<3, false>(hid_t, const void*, hsize_t);

template void gmx::h5mdio::readData<1, false>(hid_t, hsize_t, void**, size_t*, size_t*);
template void gmx::h5mdio::readData<1, true>(hid_t, hsize_t, void**, size_t*, size_t*);
template void gmx::h5mdio::readData<3, false>(hid_t, hsize_t, void**, size_t*, size_t*);
template void gmx::h5mdio::readData<1>(hid_t, hsize_t, void**);

template void gmx::h5mdio::setAttribute<int>(hid_t, const char*, int, hid_t);
template void gmx::h5mdio::setAttribute<int64_t>(hid_t, const char*, int64_t, hid_t);
template void gmx::h5mdio::setAttribute<float>(hid_t, const char*, float, hid_t);
template void gmx::h5mdio::setAttribute<double>(hid_t, const char*, double, hid_t);

template bool gmx::h5mdio::getAttribute<int64_t>(hid_t, const char*, int64_t*);
