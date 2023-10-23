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

#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#if GMX_USE_HDF5
#include <hdf5.h>
#include "external/SZ/hdf5-filter/H5Z-SZ/include/H5Z_SZ.h"
#include "external/SZ3/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif

GmxHdf5MdDataBlock::GmxHdf5MdDataBlock()
{
    strcpy(name_, "");
    numFramesPerChunk_ = 1;
    numEntries_ = 0;
    numValuesPerEntry_ = 1;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
}

GmxHdf5MdDataBlock::GmxHdf5MdDataBlock(const char* name, hsize_t numFramesPerChunk, hsize_t numEntries, hsize_t numValuesPerEntry)
{
    strncpy(name_, name, 15);
    numFramesPerChunk_ = numFramesPerChunk;
    numEntries_ = numEntries;
    numValuesPerEntry_ = numValuesPerEntry;
    numDatasetFrames_ = 0;
    numWrittenFrames_ = 0;
}

void GmxHdf5MdDataBlock::initDataProperties(hsize_t numFramesPerChunk, hsize_t numEntries, hsize_t numValuesPerEntry)
{
    if (numWrittenFrames_ != 0)
    {
        gmx_file("Cannot change number of frames per chunk after writing.");
    }
    printf("initNumFramesPerChunk\n");
    numFramesPerChunk_ = numFramesPerChunk;
    numEntries_ = numEntries;
    numValuesPerEntry_ = numValuesPerEntry;
}

void GmxHdf5MdDataBlock::setupForWriting(hsize_t numFramesPerChunk, hsize_t numEntries, hsize_t numValuesPerEntry)
{
    printf("Setup for writing\n");
    if (numWrittenFrames_ == 0)
    {
        initDataProperties(numFramesPerChunk, numEntries, numValuesPerEntry);
    }
}

template<typename T>
void GmxHdf5MdDataBlock::writeFrame(hid_t            container,
                                    const T*         data)
{
#if GMX_USE_HDF5
    hid_t dataset = H5Dopen(container, name_, H5P_DEFAULT);
#if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    if (dataset < 0)
    {
        hsize_t dataSize[3] = {numFramesPerChunk_, numEntries_, numValuesPerEntry_};
        hsize_t maxDims[3] = {H5S_UNLIMITED, numEntries_, numValuesPerEntry_};
        hid_t dataspace = H5Screate_simple(3, dataSize, maxDims);
        if (dataspace < 0)
        {
            gmx_file("Cannot create dataspace.");
        }
        hsize_t chunkDims[3] = {numFramesPerChunk_, numEntries_, numValuesPerEntry_};
        printf("%s, ChunkDims: %" PRId64 " %" PRId64 " %" PRId64 "\n", name_, chunkDims[0], chunkDims[1], chunkDims[2]);
        hid_t propertyList = H5Pcreate(H5P_DATASET_CREATE);
        if (propertyList < 0)
        {
            gmx_file("Cannot create property list.");
        }
        if (H5Pset_chunk(propertyList, 3, chunkDims) < 0)
        {
            gmx_file("Cannot set chunk dimensions.");
        }

        int sz_mode = 0; //0: ABS, 1: REL
        size_t numCompressionSettingsElements;
        unsigned int *compressionSettings = nullptr;
        printf("Setting SZ_errConfigToCdArray\n");
        SZ_errConfigToCdArray(&numCompressionSettingsElements, &compressionSettings, sz_mode, 0.005, 0.005, 0, 0);
        printf("Setting SZ3 filter\n");
        if (H5Pset_filter(propertyList, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings) < 0)
        {
            gmx_file("Cannot set SZ compression.");
        }

        // int sz3_mode = 0; //0: ABS, 1: REL
        // size_t numCompressionSettingsElements;
        // unsigned int *compressionSettings = nullptr;
        // printf("Setting SZ3_errConfigToCdArray\n");
        // SZ_errConfigToCdArray(&numCompressionSettingsElements, &compressionSettings, sz3_mode, 0.001, 0.001, 0, 0);
        // printf("Setting SZ3 filter\n");
        // if (H5Pset_filter(propertyList, H5Z_FILTER_SZ3, H5Z_FLAG_MANDATORY, numCompressionSettingsElements, compressionSettings) < 0)
        // {
        //     gmx_file("Cannot set SZ3 compression.");
        // }
        // if(H5Zfilter_avail(H5Z_FILTER_SZ3))
        // {
        //     unsigned filterConfig;
        //     H5Zget_filter_info(H5Z_FILTER_SZ3, &filterConfig);
        //
        //     if(filterConfig & H5Z_FILTER_CONFIG_ENCODE_ENABLED)
        //     {
        //         printf("SZ3 filter is available for encoding and decoding.\n");
        //     }
        //     else
        //     {
        //         printf("SZ3 filter is NOT available for encoding and decoding!\n");
        //     }
        // }

        // /* This gives a lossy compression with 0.001 precision. */
        // if (H5Pset_scaleoffset(propertyList, H5Z_SO_FLOAT_DSCALE, 3) < 0)
        // {
        //     gmx_file("Cannot set scale offset filter.");
        // }
        // if (H5Pset_shuffle(propertyList) < 0)
        // {
        //     gmx_file("Cannot set shuffle filter.");
        // }
        // if (H5Pset_deflate(propertyList, 6) < 0)
        // {
        //     gmx_file("Cannot set GZIP compression.");
        // }

        dataset = H5Dcreate(container, name_, datatype, dataspace, H5P_DEFAULT, propertyList, H5P_DEFAULT);
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
        printf("Extending size from %" PRId64 " ", numDatasetFrames_);
        numDatasetFrames_ += numFramesPerChunk_;
        printf("to % " PRId64 "\n", numDatasetFrames_);
        hsize_t newDims[3] = {numDatasetFrames_, numEntries_, numValuesPerEntry_};
        H5Dset_extent(dataset, newDims);
    }
    hsize_t fileOffset[3] = {numWrittenFrames_, 0, 0};
    hsize_t outputBlockSize[3] = {1, numEntries_, numValuesPerEntry_};
    hid_t dataspace = H5Dget_space(dataset);
    if (dataspace < 0)
    {
        gmx_file("Cannot get dataspace of existing dataset.");
    }
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, fileOffset, nullptr, outputBlockSize, nullptr);

    hid_t memoryDataspace = H5Screate_simple(3, outputBlockSize, nullptr);
    printf("Writing %s.\n", name_);
    H5Dwrite(dataset, datatype, memoryDataspace, dataspace, H5P_DEFAULT, data);
    printf("%s written\n", name_);
    numWrittenFrames_++;
    H5Dclose(dataset);
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo() :
box_("box", 1, 0),
x_("positions", 1, 0),
v_("velocities", 1, 0),
f_("forces", 1, 0),
charges_{"charges", 1, 0},
masses_{"masses", 1, 0}
{
    file_ = -1;
}

GmxHdf5MdIo::GmxHdf5MdIo(const char* fileName, const char *modeString) :
box_("box", 1, 0),
x_("positions", 1, 0),
v_("velocities", 1, 0),
f_("forces", 1, 0),
charges_{"charges", 1, 0},
masses_{"masses", 1, 0}
{
    printf("In constructor. %s (%s)\n", fileName, modeString);
    file_ = -1;
    openFile(fileName, modeString);
}

GmxHdf5MdIo::~GmxHdf5MdIo()
{
    if(file_ != -1)
    {
        closeFile();
    }
}

void GmxHdf5MdIo::openFile(const char* fileName, const char* modeString)
{
#if GMX_USE_HDF5
    std::string modeStringLower = gmx::toLowerCase(modeString);

    bool read = true;
    bool write = false;
    bool append = false;
    std::size_t found = modeStringLower.find("w");
    if(found < std::string::npos)
    {
        write = true;
    }
    found = modeStringLower.find("a");
    if(found < std::string::npos)
    {
        append = true;
        write = true;
    }
    closeFile();

    printf("Opening %s with mode %s\n", fileName, modeString);
    if (write)
    {
        if (append)
        {
            file_ = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        if (file_ < 0)
        {
            file_ = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
    }
    else
    {
        file_ = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    }
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
#if GMX_USE_HDF5
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
#if GMX_USE_HDF5
    if (file_ >= 0)
    {
        printf("Flushing.\n");
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::setupMolecularSystem(const gmx_mtop_t& topology)
{
#if GMX_USE_HDF5
    std::vector<real> atomCharges;
    std::vector<real> atomMasses;

    atomCharges.reserve(topology.natoms);
    atomMasses.reserve(topology.natoms);

    printf("Setup molsystem\n");
    for (const gmx_molblock_t& molBlock : topology.molblock)
    {
        const gmx_moltype_t* molType = &topology.moltype[molBlock.type];
        for (int atomCounter = 0; atomCounter < molType->atoms.nr; atomCounter++)
        {
            atomCharges.push_back(molType->atoms.atom[atomCounter].q);
            atomMasses.push_back(molType->atoms.atom[atomCounter].m);
        }
        for (int molCounter = 1; molCounter < molBlock.nmol; molCounter++)
        {
            std::copy_n(atomCharges.end() - molType->atoms.nr,
                        molType->atoms.nr,
                        std::back_inserter(atomCharges));
            std::copy_n(atomMasses.end() - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomMasses));
        }
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
    printf("Setting up for writing charges.");
    charges_.setupForWriting(1, topology.natoms, 1);
    charges_.writeFrame(particlesGroup, atomCharges.data());
    masses_.setupForWriting(1, topology.natoms, 1);
    masses_.writeFrame(particlesGroup, atomMasses.data());

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
#if GMX_USE_HDF5
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

    const int numFramesPerChunk = 10;

    if (box != nullptr)
    {
        box_.setupForWriting(numFramesPerChunk, DIM, DIM);
        box_.writeFrame(particlesGroup, box);
    }

    if (x != nullptr)
    {
        x_.setupForWriting(numFramesPerChunk, numAtoms, DIM);
        x_.writeFrame(particlesGroup, x);
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif

}
