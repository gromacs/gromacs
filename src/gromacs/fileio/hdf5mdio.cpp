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
#include "tngio.h"

#include "config.h"

#include <functional>
#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#if GMX_USE_HDF5
#    include <h5xx/h5xx.hpp>
#endif

static unsigned convertModeStringToInt(const std::string& modeString)
{
#if GMX_USE_HDF5
    const gmx_unused int maxStringLength = 4; // 3 is standard max, but there may also be a "+", which is ignored
    GMX_ASSERT(modeString.length() <= maxStringLength, "The mode string is too long");
    std::string modeStringLower = gmx::toLowerCase(modeString);
    unsigned mode = h5xx::file::in; // Reading is always enabled.

    std::size_t found = modeStringLower.find("w");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::out;
    }
    // Treat "a" (append) as write
    found = modeStringLower.find("a");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::out;
    }
    found = modeStringLower.find("t");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::trunc;
    }
    found = modeStringLower.find("e");
    if(found < std::string::npos)
    {
        mode |= h5xx::file::excl;
    }

    return mode;
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

static void inline writeMultiScalarToDataset(const h5xx::group& group, const std::string& datasetName, const rvec *data, const int numElements)
{
#if GMX_USE_HDF5
    for (int i = 0; i < numElements; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
                h5xx::write_dataset(group, datasetName, data[i][j]);
        }
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5
    file_ = new h5xx::file();
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::GmxHdf5MdIo(const std::string &fileName, const std::string &modeString)
{
#ifdef GMX_USE_HDF5
    printf("In constructor. %s (%s)\n", fileName.c_str(), modeString.c_str());
    file_ = new h5xx::file();
    openFile(fileName, modeString);
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

GmxHdf5MdIo::~GmxHdf5MdIo()
{
#ifdef GMX_USE_HDF5
    if(file_ != nullptr)
    {
        closeFile();
        delete file_;
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::openFile(const std::string &fileName, const std::string &modeString)
{
#ifdef GMX_USE_HDF5
    printf("Opening %s with mode %s\n", fileName.c_str(), modeString.c_str());
    unsigned mode = convertModeStringToInt(modeString);
    file_->open(fileName, mode);
    printf("Opened file\n");
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::closeFile()
{
#ifdef GMX_USE_HDF5
    file_->flush();
    file_->close();
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::flush()
{
#ifdef GMX_USE_HDF5
    file_->flush();
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxHdf5MdIo::writeFrame(int64_t          step,
                             real             time,
                             real             lambda,
                             const rvec*      box,
                             int              natoms,
                             const rvec*      x,
                             const rvec*      v,
                             const rvec*      f)
{
    std::string particlesGroupName = "particles";
    h5xx::group particlesGroup(*file_, particlesGroupName);
    printf("Created group. hid: %d\n", particlesGroup.hid());

    std::string boxDatasetName = "box";
    constexpr int c_boxDataSize = 3;
    const int numFramesPerChunk = 5;

    if (!h5xx::exists_dataset(particlesGroup, boxDatasetName))
    {
        std::vector<size_t> chunkDims{c_boxDataSize};
        const real* boxPointer = reinterpret_cast<const real*>(box);
        std::vector<real> boxVector (boxPointer, boxPointer + c_boxDataSize * DIM);
        // std::vector<real> boxVector {0,1,2,3,4,5,6,7,8};
        printf("%s, ChunkDims: %d %d\n", boxDatasetName.c_str(), chunkDims.data()[0], chunkDims.data()[1]);
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        // storagePolicy.add(h5xx::policy::filter::deflate());

        h5xx::dataset boxDataset = h5xx::create_dataset<real, h5xx::group, h5xx::policy::storage::chunked>
        (particlesGroup, boxDatasetName, storagePolicy);
        // h5xx::dataset boxDataset = h5xx::create_dataset<h5xx::group, std::vector<real>>(particlesGroup, boxDatasetName, boxVector);
        hid_t propertyListId(boxDataset.hid());
        printf("%d\n", propertyListId);
        storagePolicy.set_storage(propertyListId);
        printf("Created box dataset\n");
    }
    writeMultiScalarToDataset(particlesGroup, boxDatasetName, box, c_boxDataSize);

    std::string positionDatasetName = "position";
    if (!h5xx::exists_dataset(particlesGroup, positionDatasetName))
    {
        std::vector<size_t> chunkDims{numFramesPerChunk, natoms * DIM};
        printf("ChunkDims: %d %d %d\n", chunkDims.data()[0], chunkDims.data()[1], chunkDims.data()[2]);
        h5xx::policy::storage::chunked storagePolicy(chunkDims);
        // storagePolicy.add(h5xx::policy::filter::deflate());

        h5xx::create_dataset<real, h5xx::group, h5xx::policy::storage::chunked>
        (particlesGroup, positionDatasetName, storagePolicy);
        printf("Created position dataset\n");
    }
    writeMultiScalarToDataset(particlesGroup, positionDatasetName, x, natoms);
}
