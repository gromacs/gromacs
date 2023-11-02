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

#include "h5md_util.h"

#include "config.h"

#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/fatalerror.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#include <hdf5.h>
#include "external/SZ3/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif

static void setNumericFillValue(hid_t datasetCreatePropertyList, const hid_t datatype)
{
    if(H5Tequal(datatype, H5T_NATIVE_INT))
    {
        const int dataFill = -1;
        H5Pset_fill_value(datasetCreatePropertyList, datatype, &dataFill);
    }
    else if(H5Tequal(datatype, H5T_NATIVE_INT64))
    {
        const int64_t dataFill = -1;
        H5Pset_fill_value(datasetCreatePropertyList, datatype, &dataFill);
    }
    else if(H5Tequal(datatype, H5T_NATIVE_FLOAT))
    {
        const float dataFill = -1;
        H5Pset_fill_value(datasetCreatePropertyList, datatype, &dataFill);
    }
    else if(H5Tequal(datatype, H5T_NATIVE_DOUBLE))
    {
        const double dataFill = -1;
        H5Pset_fill_value(datasetCreatePropertyList, datatype, &dataFill);
    }
}

hid_t openOrCreateGroup(hid_t container, const char *name)
{
#if GMX_USE_HDF5
    hid_t group = H5Gopen(container, name, H5P_DEFAULT);
    if (group < 0)
    {
        hid_t linkPropertyList = H5Pcreate(H5P_LINK_CREATE);     // create group creation property list
        H5Pset_create_intermediate_group(linkPropertyList, 1);   // set intermediate link creation
        group = H5Gcreate(container, name, linkPropertyList, H5P_DEFAULT, H5P_DEFAULT);
        if( group < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            gmx_file("Cannot create group.");
        }
    }
    return group;
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void registerSz3FilterImplicitly()
{
#if GMX_USE_HDF5
    hid_t propertyList = H5Pcreate(H5P_DATASET_CREATE);
    int sz3_mode = 0; //0: ABS, 1: REL
    size_t numCompressionSettingsElements;
    unsigned int *compressionSettings = nullptr;
    SZ_errConfigToCdArray(&numCompressionSettingsElements, &compressionSettings, sz3_mode, 0, 0, 0, 0);
    if (H5Pset_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Cannot use SZ3 compression filter.");
    }
#endif
}

void writeData(hid_t container, const char* name, const char* unit, const void* data, hsize_t numFramesPerChunk, hsize_t numEntries, hsize_t numValuesPerEntry, hsize_t positionToWrite, hid_t datatype, CompressionAlgorithm compression, double compressionError)
{
#if GMX_USE_HDF5
    /* Set a reasonable cache based on chunk sizes. The cache is not stored in file, so must be set when opening a dataset */
    hsize_t chunkDims[3] = {numFramesPerChunk, numEntries, numValuesPerEntry};
    size_t cacheSize = sizeof(real) * chunkDims[0] * chunkDims[1] * chunkDims[2];
    hid_t accessPropertyList = H5Pcreate(H5P_DATASET_ACCESS);
    H5Pset_chunk_cache(accessPropertyList, H5D_CHUNK_CACHE_NSLOTS_DEFAULT, cacheSize, H5D_CHUNK_CACHE_W0_DEFAULT);
    hid_t dataset = H5Dopen(container, name, accessPropertyList);
    if (dataset < 0)
    {
        hsize_t dataSize[3] = {numFramesPerChunk, numEntries, numValuesPerEntry};
        hsize_t maxDims[3] = {H5S_UNLIMITED, numEntries, numValuesPerEntry};
        hid_t dataspace = H5Screate_simple(3, dataSize, maxDims);
        hid_t createPropertyList = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(createPropertyList, 3, chunkDims);
        setNumericFillValue(createPropertyList, datatype);

        switch(compression)
        {
        case CompressionAlgorithm::LossySz3:
        {
            // int sz3_mode = 0; //0: ABS, 1: REL
            // size_t numCompressionSettingsElements;
            // unsigned int *compressionSettings = nullptr;
            // SZ_errConfigToCdArray(&numCompressionSettingsElements, &compressionSettings, sz3_mode, compressionError, compressionError, 0, 0);
            // if (H5Pset_filter(createPropertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings) < 0)
            if (H5Pset_filter(createPropertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, 0, nullptr) < 0)
            {
                H5Eprint2(H5E_DEFAULT, nullptr);
                gmx_file("Cannot set SZ3 compression.");
            }
            if(H5Zfilter_avail(H5Z_FILTER_SZ3) < 0)
            {
                H5Eprint2(H5E_DEFAULT, nullptr);
                gmx_file("SZ3 filter not available.");
            }
            break;
        }
        case CompressionAlgorithm::LosslessWithShuffle:
            if (H5Pset_shuffle(createPropertyList) < 0)
            {
                gmx_file("Cannot set shuffle filter.");
            }
        case CompressionAlgorithm::LosslessNoShuffle:
            if (H5Pset_deflate(createPropertyList, 6) < 0)
            {
                H5Eprint2(H5E_DEFAULT, nullptr);
                gmx_file("Cannot set GZIP compression.");
            }
            break;
        case CompressionAlgorithm::None:
            break;
        }

        dataset = H5Dcreate(container, name, datatype, dataspace, H5P_DEFAULT, createPropertyList, accessPropertyList);
        if(dataset < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            gmx_file("Cannot create dataset.");
        }

        if (unit != nullptr)
        {
            char unitElementString[] = "unit";
            setAttribute(dataset, unitElementString, unit);
        }
    }
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t currentDims[3], maxDims[3];
    H5Sget_simple_extent_dims(dataspace, currentDims, maxDims);
    /* Resize the dataset if needed. */
    if(positionToWrite >= currentDims[0])
    {
        hsize_t newDims[3] = {(positionToWrite / numFramesPerChunk + 1) * numFramesPerChunk, currentDims[1], currentDims[2]};
        if (debug)
        {
            fprintf(debug, "Resizing dataset from %" PRId64" to %" PRId64 "\n", currentDims[0], newDims[0]);
        }
        H5Dset_extent(dataset, newDims);
        dataspace = H5Dget_space(dataset);
    }
    hsize_t fileOffset[3] = {positionToWrite, 0, 0};
    hsize_t outputBlockSize[3] = {1, numEntries, numValuesPerEntry};
    if (dataspace < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Cannot get dataspace of existing dataset.");
    }
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, fileOffset, nullptr, outputBlockSize, nullptr);

    hid_t memoryDataspace = H5Screate_simple(3, outputBlockSize, nullptr);
    if (H5Dwrite(dataset, datatype, memoryDataspace, dataspace, H5P_DEFAULT, data) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Error writing data.");
    }

    // It would be good to close the dataset here, but that means compressing and writing the whole chunk every time - very slow.
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

template <typename T>
void setAttribute(hid_t container, const char *name, const T value, hid_t dataType)
{
    hid_t attribute = H5Aopen(container, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hid_t dataspace = H5Screate(H5S_SCALAR);
        attribute = H5Acreate2(container, name, dataType, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Awrite(attribute, dataType, &value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

void setAttribute(hid_t container, const char *name, const char* value)
{
    hid_t dataType = H5Tcopy(H5T_C_S1);
    H5Tset_size(dataType, H5T_VARIABLE);
    H5Tset_strpad(dataType, H5T_STR_NULLTERM);
    H5Tset_cset(dataType, H5T_CSET_UTF8);

    setAttribute(container, name, value, dataType);
}

template void setAttribute<int>(hid_t, const char*, int, hid_t);
template void setAttribute<float>(hid_t, const char*, float, hid_t);
template void setAttribute<double>(hid_t, const char*, double, hid_t);
template void setAttribute<char *>(hid_t, const char*, char *, hid_t);
