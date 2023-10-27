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

#include <string>

#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

#include "h5md_datablock.h"
#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#include <hdf5.h>
#endif


GmxH5mdIo::GmxH5mdIo(const char* fileName, const char *modeString)
{
    file_ = -1;
    numFramesPerChunkCompressed_ = 0;
    if (strlen(fileName) > 0)
    {
        openFile(fileName, modeString);
    }
}

GmxH5mdIo::~GmxH5mdIo()
{
    if(file_ != -1)
    {
        closeFile();
    }
}

void GmxH5mdIo::openFile(const char* fileName, const char* modeString)
{
#if GMX_USE_HDF5
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error output, e.g. when items are not found.
    std::string modeStringLower = gmx::toLowerCase(modeString);

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
        /* Create H5MD groups */
        openOrCreateGroup(file_, "h5md");
        openOrCreateGroup(file_, "particles");
        openOrCreateGroup(file_, "particles/system");
        setAuthorAndCreator();
    }
    else
    {
        file_ = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    }
    if (file_ < 0)
    {
        gmx_file("Cannot open file.");
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::closeFile()
{
#if GMX_USE_HDF5
    if (file_ >= 0)
    {
        H5Fflush(file_, H5F_SCOPE_LOCAL);
        printf("Closing file.\n");
        H5Fclose(file_);
        file_ = -1;
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::flush()
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

void GmxH5mdIo::setUpParticlesDataBlocks(bool writeCoordinates, bool writeCoordinatesCompressed, bool writeForces, bool writeVelocities, int numParticles, int numParticlesCompressed, double compressionError)
{
    if (writeCoordinates || writeForces || writeVelocities)
    {
        openOrCreateGroup(file_, "particles/system");
        hid_t boxGroup = openOrCreateGroup(file_, "particles/system/box");
        setAttribute(boxGroup, "dimension", DIM, H5T_NATIVE_INT);
    }

#if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    if (writeCoordinates)
    {
        position_ = GmxH5mdDataBlock(file_, "/particles/system/position", "nm", 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeForces)
    {
        force_ = GmxH5mdDataBlock(file_, "/particles/system/force", "kJ/(mol*nm)", 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeVelocities)
    {
        velocity_ = GmxH5mdDataBlock(file_, "/particles/system/velocity", "nm/ps", 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeCoordinates || writeForces || writeVelocities)
    {
        box_ = GmxH5mdDataBlock(file_, "/particles/system/box/edges", "nm", 1, DIM, DIM, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    }
    if (writeCoordinatesCompressed)
    {
        openOrCreateGroup(file_, "particles/selection_compressed");
        hid_t boxGroup = openOrCreateGroup(file_, "particles/selection_compressed/box");
        setAttribute(boxGroup, "dimension", DIM, H5T_NATIVE_INT);
        positionLossy_ = GmxH5mdDataBlock(file_, "/particles/selection_compressed/position", "nm", numFramesPerChunkCompressed_, numParticlesCompressed, DIM, datatype, CompressionAlgorithm::LossySz3, compressionError);
        boxLossy_ = GmxH5mdDataBlock(file_, "/particles/selection_compressed/box/edges", "nm", numFramesPerChunkCompressed_, DIM, DIM, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    }
}

template <typename T>
void GmxH5mdIo::setAttribute(hid_t container, const char *name, const T value, hid_t dataType)
{
    hid_t attribute = H5Aopen(container, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hid_t propertyList = H5Pcreate(H5P_ATTRIBUTE_CREATE);
        H5Pset_char_encoding(propertyList, H5T_CSET_UTF8);
        hid_t dataspace = H5Screate(H5S_SCALAR);
        attribute = H5Acreate2(container, name, dataType, dataspace, propertyList, H5P_DEFAULT);
        if (attribute < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            gmx_file("Cannot create attribute.");
        }
    }
    if (H5Awrite(attribute, dataType, &value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

void GmxH5mdIo::setAttribute(hid_t container, const char *name, const char* value)
{
    size_t stringLength = strlen(value);
    hid_t dataType = H5Tcopy(H5T_C_S1);
    H5Tset_size(dataType, stringLength);
    H5Tset_strpad(dataType, H5T_STR_NULLTERM);

    hid_t attribute = H5Aopen(container, name, H5P_DEFAULT);
    if (attribute < 0)
    {
        hid_t propertyList = H5Pcreate(H5P_ATTRIBUTE_CREATE);
        H5Pset_char_encoding(propertyList, H5T_CSET_UTF8);
        hid_t dataspace = H5Screate(H5S_SCALAR);
        attribute = H5Acreate2(container, name, dataType, dataspace, propertyList, H5P_DEFAULT);
        if (attribute < 0)
        {
            H5Eprint2(H5E_DEFAULT, nullptr);
            gmx_file("Cannot create attribute.");
        }
    }
    if (H5Awrite(attribute, dataType, value) < 0)
    {
        H5Eprint2(H5E_DEFAULT, nullptr);
        gmx_file("Cannot write attribute.");
    }
    H5Aclose(attribute);
}

void GmxH5mdIo::setAuthorAndCreator()
{
    const char* precisionString = "";
#    if GMX_DOUBLE
        precisionString = " (double precision)";
#    endif

    char programInfo[128];
    sprintf(programInfo, "%.100s %.24s", gmx::getProgramContext().displayName(), precisionString);
    hid_t creatorGroup = openOrCreateGroup(file_, "h5md/creator");
    setAttribute(creatorGroup, "name", programInfo);
    const char* gmxVersion = gmx_version();
    setAttribute(creatorGroup, "version", gmxVersion);
    hid_t authorGroup = openOrCreateGroup(file_, "h5md/author");
    char username[256];
    if (!gmx_getusername(username, 256))
    {
        setAttribute(authorGroup, "name", username);
    }
}

void GmxH5mdIo::setupMolecularSystem(const gmx_mtop_t& topology)
{
#if GMX_USE_HDF5
    if (file_ < 0)
    {
        gmx_file("No file open for writing");
    }

    std::vector<real> atomCharges;
    std::vector<real> atomMasses;

    atomCharges.reserve(topology.natoms);
    atomMasses.reserve(topology.natoms);

    // TODO: Check this and set better values.
    numFramesPerChunkCompressed_ = std::min(11, int(std::ceil(1000000/topology.natoms)) + 1);

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

#if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    charge_ = GmxH5mdDataBlock(file_, "/particles/charge", "e", 1, topology.natoms, 1, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    charge_.writeFrame(atomCharges.data(), 0, 0);
    mass_ = GmxH5mdDataBlock(file_, "/particles/mass", "g/mol", 1, topology.natoms, 1, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    mass_.writeFrame(atomMasses.data(), 0, 0);

#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::writeFrame(int64_t        step,
                           real           time,
                           real           lambda,
                           const rvec*    box,
                           const rvec*    x,
                           const rvec*    v,
                           const rvec*    f,
                           const rvec*    xLossy)
{
#if GMX_USE_HDF5
    if (file_ < 0)
    {
        gmx_file("No file open for writing");
    }

    if (x != nullptr)
    {
        position_.writeFrame(x, step, time);
    }
    if (v != nullptr)
    {
        velocity_.writeFrame(v, step, time);
    }
    if (f != nullptr)
    {
        force_.writeFrame(f, step, time);
    }
    if ((x != nullptr || v != nullptr || f != nullptr) && box != nullptr)
    {
        box_.writeFrame(box, step, time);
    }
    if (xLossy != nullptr)
    {
        positionLossy_.writeFrame(xLossy, step, time);
        if(box != nullptr)
        {
            boxLossy_.writeFrame(box, step, time);
        }
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif

}
