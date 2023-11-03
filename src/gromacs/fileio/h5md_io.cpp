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

#include <functional>
#include <string>

#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/sysinfo.h"

#include "h5md_datablock.h"
#include "h5md_util.h"

#define GMX_USE_HDF5 1 // FIXME: Temporary just for the editor

#if GMX_USE_HDF5
#include <hdf5.h>
#include "external/SZ3/tools/H5Z-SZ3/include/H5Z_SZ3.hpp"
#endif


GmxH5mdIo::GmxH5mdIo(const char* fileName, const char mode)
{
    file_ = -1;
    if (strlen(fileName) > 0)
    {
        openFile(fileName, mode);
    }
}

GmxH5mdIo::~GmxH5mdIo()
{
    if(file_ != -1)
    {
        closeFile();
    }
}

void GmxH5mdIo::openFile(const char* fileName, const char mode)
{
#if GMX_USE_HDF5
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error output, e.g. when items are not found.

    closeFile();

    if (debug)
    {
        fprintf(debug, "Opening H5MD file %s with mode %c\n", fileName, mode);
    }
    if (mode == 'w' || mode == 'a')
    {
        if (mode == 'w')
        {
            make_backup(fileName);
            hid_t createPropertyList = H5Pcreate(H5P_FILE_CREATE);
            if(H5Pset_file_space_strategy(createPropertyList, H5F_FSPACE_STRATEGY_FSM_AGGR, 1, 1) < 0)
            {
                printf("Cannot set H5MD file space strategy.\n");
            }
            file_ = H5Fcreate(fileName, H5F_ACC_TRUNC, createPropertyList, H5P_DEFAULT);
            setAuthorAndCreator();
        }
        else
        {
            file_ = H5Fopen(fileName, H5F_ACC_RDWR, H5P_DEFAULT);
        }
        /* Create H5MD groups. They should already be there if appending to a valid H5MD file, but it's better to be on the safe side. */
        openOrCreateGroup(file_, "h5md");
        openOrCreateGroup(file_, "particles");
        openOrCreateGroup(file_, "particles/system");
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
        if (debug)
        {
            fprintf(debug, "Closing H5MD file.\n");
        }
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
        if (debug)
        {
            fprintf(debug, "Flushing H5MD file.\n");
        }
        H5Fflush(file_, H5F_SCOPE_LOCAL);
    }
#else
    gmx_file("GROMACS was compiled without HDF5 support, cannot handle this file type");
#endif
}

void GmxH5mdIo::setUpParticlesDataBlocks(int writeCoordinatesSteps, int writeCoordinatesCompressedSteps, int writeForcesSteps, int writeVelocitiesSteps, int numParticles, int numParticlesCompressed, double compressionError)
{
#if GMX_DOUBLE
    const hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t datatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    hid_t systemGroup = openOrCreateGroup(file_, "particles/system");
    hid_t boxGroup = openOrCreateGroup(systemGroup, "box");
    setAttribute(boxGroup, "dimension", DIM, H5T_NATIVE_INT);
    box_ = GmxH5mdDataBlock(boxGroup, "edges", "value", "nm", writeCoordinatesSteps, 1, DIM, DIM, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    // TODO: Write box 'boundary' attribute ('periodic' or 'none')

    if (writeCoordinatesSteps > 0)
    {
        position_ = GmxH5mdDataBlock(systemGroup, "position", "value", "nm", writeCoordinatesSteps, 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeForcesSteps > 0)
    {
        force_ = GmxH5mdDataBlock(systemGroup, "force", "value", "kJ mol-1 nm-1)", writeForcesSteps, 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeVelocitiesSteps > 0)
    {
        velocity_ = GmxH5mdDataBlock(systemGroup, "velocity", "value", "nm ps-1", writeVelocitiesSteps, 1, numParticles, DIM, datatype, CompressionAlgorithm::LosslessWithShuffle, 0);
    }
    if (writeCoordinatesCompressedSteps > 0)
    {
        /* Use no more than 20 frames per chunk (compression unit). Use fewer frames per chunk if there are many atoms. */
        hsize_t numFramesPerChunkCompressed = std::min(20, int(std::ceil(1e6/numParticlesCompressed)));

        hid_t compressedGroup = openOrCreateGroup(file_, "particles/selection_compressed");
        hid_t boxGroup = openOrCreateGroup(compressedGroup, "box");
        setAttribute(boxGroup, "dimension", DIM, H5T_NATIVE_INT);
        boxLossy_ = GmxH5mdDataBlock(boxGroup, "edges", "value", "nm", writeCoordinatesCompressedSteps, numFramesPerChunkCompressed, DIM, DIM, datatype, CompressionAlgorithm::LosslessNoShuffle, 0);
        // TODO: Write box 'boundary' attribute ('periodic' or 'none')

        /* Register the SZ3 filter. This is not necessary when creating a dataset with the filter, but must be done to append to an existing file (e.g. when restarting from checkpoint).*/
        registerSz3FilterImplicitly();
        positionLossy_ = GmxH5mdDataBlock(compressedGroup, "position", "value", "nm", writeCoordinatesCompressedSteps, numFramesPerChunkCompressed, numParticlesCompressed, DIM, datatype, CompressionAlgorithm::LossySz3, compressionError);
    }
}

void GmxH5mdIo::setAuthorAndCreator()
{
    const char* precisionString = "";
#if GMX_DOUBLE
    precisionString = " (double precision)";
#endif

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
    std::vector<std::string> atomNames;

    atomCharges.reserve(topology.natoms);
    atomMasses.reserve(topology.natoms);
    atomNames.reserve(topology.natoms);

    for (const gmx_molblock_t& molBlock : topology.molblock)
    {
        const gmx_moltype_t* molType = &topology.moltype[molBlock.type];
        for (int atomCounter = 0; atomCounter < molType->atoms.nr; atomCounter++)
        {
            atomCharges.push_back(molType->atoms.atom[atomCounter].q);
            atomMasses.push_back(molType->atoms.atom[atomCounter].m);
            atomNames.push_back(*(molType->atoms.atomname[atomCounter]));
        }
        for (int molCounter = 1; molCounter < molBlock.nmol; molCounter++)
        {
            std::copy_n(atomCharges.end() - molType->atoms.nr,
                        molType->atoms.nr,
                        std::back_inserter(atomCharges));
            std::copy_n(atomMasses.end() - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomMasses));
            std::copy_n(atomNames.end() - molType->atoms.nr, molType->atoms.nr, std::back_inserter(atomNames));
        }
    }
    /* Is there a more convenient way to do this? std::string is nice above, but cannot be used for writing. */
    std::vector<const char *> atomNamesChars(atomNames.size());
    std::transform(atomNames.begin(), atomNames.end(), atomNamesChars.begin(), std::mem_fn(&std::string::c_str));

    hid_t stringDataType = H5Tcopy(H5T_C_S1);
    H5Tset_size(stringDataType, H5T_VARIABLE);
    H5Tset_strpad(stringDataType, H5T_STR_NULLTERM);
    H5Tset_cset(stringDataType, H5T_CSET_UTF8);
    /* FIXME: Currently atom names cannot change during the simulation. */
    atomName_ = GmxH5mdDataBlock(file_, "/particles/system/atomname", "value", "", 0, 1, topology.natoms, 1, stringDataType, CompressionAlgorithm::LosslessNoShuffle, 0);
    atomName_.writeFrame(atomNamesChars.data(), 0, 0);

#if GMX_DOUBLE
    const hid_t floatDatatype = H5Tcopy(H5T_NATIVE_DOUBLE);
#else
    const hid_t floatDatatype = H5Tcopy(H5T_NATIVE_FLOAT);
#endif
    /* FIXME: Currently charges and masses cannot change during the simulation. */
    charge_ = GmxH5mdDataBlock(file_, "/particles/system/charge", "value", "e", 0, 1, topology.natoms, 1, floatDatatype, CompressionAlgorithm::LosslessNoShuffle, 0);
    charge_.writeFrame(atomCharges.data(), 0, 0);
    mass_ = GmxH5mdDataBlock(file_, "/particles/system/mass", "value", "g mol-1", 0, 1, topology.natoms, 1, floatDatatype, CompressionAlgorithm::LosslessNoShuffle, 0);
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
        box_.writeFrame(box, step, time);
    }
    if (v != nullptr)
    {
        velocity_.writeFrame(v, step, time);
    }
    if (f != nullptr)
    {
        force_.writeFrame(f, step, time);
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

extern template void setAttribute<int>(hid_t, const char*, int, hid_t);
extern template void setAttribute<float>(hid_t, const char*, float, hid_t);
extern template void setAttribute<double>(hid_t, const char*, double, hid_t);
extern template void setAttribute<char *>(hid_t, const char*, char *, hid_t);
