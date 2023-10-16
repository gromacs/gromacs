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

#include "hdf5mdio.h"

#include "config.h"

#include <functional>
#include <string>
#include <sys/_types/_int64_t.h>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#if GMX_USE_HDF5
#    include <hdf5.h>
#endif

static unsigned convertModeStringToInt(const std::string& modeString)
{
#if GMX_USE_HDF5
    const gmx_unused int maxStringLength = 4; // 3 is standard max, but there may also be a "+", which is ignored
    GMX_ASSERT(modeString.length() <= maxStringLength, "The mode string is too long");
    std::string modeStringLower = gmx::toLowerCase(modeString);
    unsigned mode = H5F_ACC_RDONLY; // Reading is always enabled.

    std::size_t found = modeStringLower.find("w");
    if(found < std::string::npos)
    {
        mode |= H5F_ACC_TRUNC;
    }
    // Treat "a" (append) as write
    found = modeStringLower.find("a");
    if(found < std::string::npos)
    {
        mode |= H5F_ACC_TRUNC;
    }
    found = modeStringLower.find("t");
    if(found < std::string::npos)
    {
        mode |= H5F_ACC_TRUNC;
    }
    found = modeStringLower.find("e");
    if(found < std::string::npos)
    {
        mode |= H5F_ACC_EXCL;
    }

    return mode;
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdParticlesBox::GmxHdf5MdParticlesBox()
{
    strcpy(name_, "box");
    numFramesPerChunk_ = 1;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
#if GMX_DOUBLE
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
}

GmxHdf5MdParticlesBox::GmxHdf5MdParticlesBox(int numFramesPerChunk)
{
    strcpy(name_, "box");
    numFramesPerChunk_ = numFramesPerChunk;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
#if GMX_DOUBLE
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
}

void GmxHdf5MdParticlesBox::initNumFramesPerChunk(int numFramesPerChunk)
{
    if (numWrittenFrames_ != 0 || numWrittenFrames_ != 0)
    {
        gmx_file("Cannot change number of frames per chunk after writing.");
    }
    numFramesPerChunk_ = numFramesPerChunk;
}

void GmxHdf5MdParticlesBox::setupForWriting(int numFramesPerChunk)
{
    if (numWrittenFrames_ != 0 || numWrittenFrames_ != 0)
    {
        initNumFramesPerChunk(numFramesPerChunk);
    }
}

void GmxHdf5MdParticlesBox::writeFrame(int64_t          step,
                                       real             time,
                                       hid_t            container,
                                       const rvec*      box)
{
    hid_t boxDataset = H5Dopen(container, name_, H5P_DEFAULT);
    if (boxDataset < 0)
    {
        hsize_t dataSize[3] = {numFramesPerChunk_, DIM, DIM};
        hsize_t maxDims[3] = {H5S_UNLIMITED, DIM, DIM};
        hid_t boxDataspace = H5Screate_simple(3, dataSize, maxDims);
        if (boxDataspace < 0)
        {
            gmx_file("Cannot create box dataspace.");
        }
        hsize_t chunkDims[3] = {numFramesPerChunk_, DIM, DIM};
        printf("%s, ChunkDims: %d %d %d\n", name_, chunkDims[0], chunkDims[1], chunkDims[2]);
        hid_t boxPropertyList = H5Pcreate(H5P_DATASET_CREATE);
        if (boxPropertyList < 0)
        {
            gmx_file("Cannot create box property list.");
        }
        if (H5Pset_chunk(boxPropertyList, 3, chunkDims) < 0)
        {
            gmx_file("Cannot set box chunk dimensions.");
        }
        if (H5Pset_deflate(boxPropertyList, 6) < 0)
        {
            gmx_file("Cannot set box GZIP compression.");
        }

        boxDataset = H5Dcreate(container, name_, datatype_, boxDataspace, H5P_DEFAULT, boxPropertyList, H5P_DEFAULT);
        if (boxDataset < 0)
        {
            gmx_file("Cannot create box dataset");
        }
        printf("Created box dataset\n");
        numDatasetFrames_ = numFramesPerChunk_;
        numWrittenFrames_ = 0;
        // H5Dclose(boxPropertyList);
        // H5Dclose(boxDataspace);
    }
    if(numWrittenFrames_ >= numDatasetFrames_)
    {
        printf("Extending size from %d ", numDatasetFrames_);
        numDatasetFrames_ += numFramesPerChunk_;
        printf("to %d\n", numDatasetFrames_);
        hsize_t newDims[3] = {numDatasetFrames_, DIM, DIM};
        H5Dset_extent(boxDataset, newDims);
    }
    hsize_t fileOffset[3] = {numWrittenFrames_, 0, 0};
    hsize_t outputBlockSize[3] = {1, DIM, DIM};
    hid_t boxDataspace = H5Dget_space(boxDataset);
    if (boxDataspace < 0)
    {
        gmx_file("Cannot get dataspace of existing box dataset.");
    }
    H5Sselect_hyperslab(boxDataspace, H5S_SELECT_SET, fileOffset, NULL, outputBlockSize, NULL);

    hid_t boxMemoryDataspace = H5Screate_simple(3, outputBlockSize, NULL);
    printf("Writing box.\n");
    H5Dwrite(boxDataset, datatype_, boxMemoryDataspace, boxDataspace, H5P_DEFAULT, box);
    printf("Box written\n");
    numWrittenFrames_++;
    H5Dclose(boxDataset);
}

GmxHdf5MdParticlesProperties::GmxHdf5MdParticlesProperties()
{
    strcpy(name_, "");
    numFramesPerChunk_ = 1;
    numAtoms_ = 0;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
#if GMX_DOUBLE
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
}

GmxHdf5MdParticlesProperties::GmxHdf5MdParticlesProperties(const char* name, int numFramesPerChunk, int64_t numAtoms)
{
    strncpy(name_, name, 15);
    numFramesPerChunk_ = numFramesPerChunk;
    numAtoms_ = numAtoms;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
#if GMX_DOUBLE
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype_ = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
}

void GmxHdf5MdParticlesProperties::initNumFramesPerChunkAndNumAtoms(int numFramesPerChunk, int64_t numAtoms)
{
    if (numWrittenFrames_ != 0 || numWrittenFrames_ != 0)
    {
        gmx_file("Cannot change number of frames per chunk after writing.");
    }
    numFramesPerChunk_ = numFramesPerChunk;
}

void GmxHdf5MdParticlesProperties::setupForWriting(int numFramesPerChunk, int64_t numAtoms)
{
    if (numWrittenFrames_ != 0 || numWrittenFrames_ != 0)
    {
        initNumFramesPerChunkAndNumAtoms(numFramesPerChunk, numAtoms);
    }
}
void GmxHdf5MdParticlesProperties::writeFrame(int64_t          step,
                                              real             time,
                                              hid_t            container,
                                              const rvec*      data)
{
    hid_t dataset = H5Dopen(container, name_, H5P_DEFAULT);
    if (dataset < 0)
    {
        hsize_t dataSize[3] = {numFramesPerChunk_, numAtoms_, DIM};
        hsize_t maxDims[3] = {H5S_UNLIMITED, numAtoms_, DIM};
        hid_t dataspace = H5Screate_simple(3, dataSize, maxDims);
        if (dataspace < 0)
        {
            gmx_file("Cannot create dataspace.");
        }
        hsize_t chunkDims[3] = {numFramesPerChunk_, numAtoms_, DIM};
        printf("%s, ChunkDims: %d %d %d\n", name_, chunkDims[0], chunkDims[1], chunkDims[2]);
        hid_t propertyList = H5Pcreate(H5P_DATASET_CREATE);
        if (propertyList < 0)
        {
            gmx_file("Cannot create property list.");
        }
        if (H5Pset_chunk(propertyList, 3, chunkDims) < 0)
        {
            gmx_file("Cannot set chunk dimensions.");
        }
        if (H5Pset_deflate(propertyList, 6) < 0)
        {
            gmx_file("Cannot set GZIP compression.");
        }

        dataset = H5Dcreate(container, name_, datatype_, dataspace, H5P_DEFAULT, propertyList, H5P_DEFAULT);
        if (dataset < 0)
        {
            gmx_file("Cannot create dataset");
        }
        printf("Created %s dataset\n", name_);
        numDatasetFrames_ = numFramesPerChunk_;
        numWrittenFrames_ = 0;
    }
    /* Resize the dataset if needed. */
    if(numWrittenFrames_ >= numDatasetFrames_)
    {
        printf("Extending size from %d ", numDatasetFrames_);
        numDatasetFrames_ += numFramesPerChunk_;
        printf("to %d\n", numDatasetFrames_);
        hsize_t newDims[3] = {numDatasetFrames_, numAtoms_, DIM};
        H5Dset_extent(dataset, newDims);
    }
    hsize_t fileOffset[3] = {numWrittenFrames_, 0, 0};
    hsize_t outputBlockSize[3] = {1, numAtoms_, DIM};
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        gmx_file("Cannot get dataspace of existing dataset.");
    }
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, fileOffset, NULL, outputBlockSize, NULL);

    hid_t memoryDataspace = H5Screate_simple(3, outputBlockSize, NULL);
    printf("Writing box.\n");
    H5Dwrite(dataset, datatype_, memoryDataspace, dataspace, H5P_DEFAULT, data);
    printf("Box written\n");
    numWrittenFrames_++;
    H5Dclose(dataset);
}

GmxHdf5MdIo::GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5
    file_ = -1;
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo(const char* fileName, const char *modeString)
{
#ifdef GMX_USE_HDF5
    printf("In constructor. %s (%s)\n", fileName, modeString);
    file_ = -1;
    openFile(fileName, modeString);
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::~GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5
    if(file_ != -1)
    {
        closeFile();
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::openFile(const char* fileName, const char* modeString)
{
#ifdef GMX_USE_HDF5
    unsigned mode = convertModeStringToInt(modeString);
    printf("Opening %s with mode %s == %d\n", fileName, modeString, mode);
    file_ = H5Fcreate(fileName, mode, H5P_DEFAULT, H5P_DEFAULT);
    if (file_ < 0)
    {
        gmx_file("Cannot open file.");
    }
    printf("Opened file\n");
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::closeFile()
{
#ifdef GMX_USE_HDF5
    if (file_ >= 0)
    {
        printf("Closing file, flushing.\n");
        H5Fflush(file_, H5F_SCOPE_LOCAL);
        printf("Closing file.\n");
        H5Fclose(file_);
        file_ = -1;
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::flush()
{
#ifdef GMX_USE_HDF5
    if (file_ >= 0)
    {
        printf("Flushing.\n");
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}


void GmxHdf5MdIo::writeFrame(int64_t          step,
                             real             time,
                             real             lambda,
                             const rvec*      box,
                             int64_t          numAtoms,
                             const rvec*      x,
                             const rvec*      v,
                             const rvec*      f)
{
    if (file_ < 0)
    {
        gmx_file("No file open for writing");
    }
    const char particlesGroupName[] = "particles";
    hid_t particlesGroup = H5Gopen(file_, particlesGroupName, H5P_DEFAULT);
    if (particlesGroup < 0)
    {
        hid_t linkPropertyList = H5Pcreate(H5P_LINK_CREATE);     // create group creation property list
        H5Pset_create_intermediate_group(linkPropertyList, 1);   // set intermediate link creation
        particlesGroup = H5Gcreate(file_, particlesGroupName, linkPropertyList, H5P_DEFAULT, H5P_DEFAULT);
        if( particlesGroup < 0)
        {
            gmx_file("Cannot create particles group.");
        }
        printf("Created group. hid: %d\n", particlesGroup);
    }

    const int numFramesPerChunk = 5;

    // box_.setupForWriting(numFramesPerChunk);
    // x_.setupForWriting(numFramesPerChunk, numAtoms);
    //
    // box_.writeFrame(step, time, particlesGroup, box);
    // x_.writeFrame(step, time, particlesGroup, x);
    // writeMultiScalarToDataset(particlesGroup, boxDatasetName, box, dataSize);
    // H5Dclose(boxDataspace);
    // H5Dclose(boxMemoryDataspace);


    //
    // std::string positionDatasetName = "position";
    // if (!h5xx::exists_dataset(particlesGroup, positionDatasetName))
    // {
    //     std::vector<size_t> chunkDims{numFramesPerChunk, natoms * DIM};
    //     printf("ChunkDims: %d %d %d\n", chunkDims.data()[0], chunkDims.data()[1], chunkDims.data()[2]);
    //     h5xx::policy::storage::chunked storagePolicy(chunkDims);
    //     // storagePolicy.add(h5xx::policy::filter::deflate());
    //
    //     h5xx::create_dataset<real, h5xx::group, h5xx::policy::storage::chunked>
    //     (particlesGroup, positionDatasetName, storagePolicy);
    //     printf("Created position dataset\n");
    // }
    // writeMultiScalarToDataset(particlesGroup, positionDatasetName, x, natoms);
}
