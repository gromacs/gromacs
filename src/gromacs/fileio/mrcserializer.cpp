/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements methods from mrcserializer.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "mrcserializer.h"

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/utility/inmemoryserializer.h"

namespace gmx
{

namespace
{

void serialize(ISerializer * serializer, std::array<int32_t, 3> * valueArray)
{
    serializer->doInt(&(*valueArray)[0]);
    serializer->doInt(&(*valueArray)[1]);
    serializer->doInt(&(*valueArray)[2]);
}

void serialize(ISerializer * serializer, std::array<float, 3> * valueArray)
{
    serializer->doFloat(&(*valueArray)[0]);
    serializer->doFloat(&(*valueArray)[1]);
    serializer->doFloat(&(*valueArray)[2]);
}

/*! \brief
 * Serializes an integer array with three integers, adding unity when writing and substracting unity when reading.
 */
void serializeOffsetPlusOne(ISerializer * serializer, std::array<int32_t, 3> * valueArray)
{
    std::array<int32_t, 3> offsetValues;
    if (!serializer->reading())
    {
        offsetValues = std::array<int32_t, 3>{{
                                                  (*valueArray)[0]+1, (*valueArray)[1]+1, (*valueArray)[2]+1
                                              }};
    }
    serialize(serializer, &offsetValues);
    if (serializer->reading())
    {
        *valueArray = std::array<int32_t, 3>{{
                                                 offsetValues[0]-1, offsetValues[1]-1, offsetValues[2]-1
                                             }};
    }
}

template <class EnumBasedOnInt32>
std::enable_if_t<std::is_same<std::underlying_type_t<EnumBasedOnInt32>, int32_t>::value, void>
serializeEnumBasedOnInt32(ISerializer *serializer, EnumBasedOnInt32 *enumValue)
{
    int32_t enumAsInt;
    if (!serializer->reading())
    {
        enumAsInt = static_cast<int32_t>(*enumValue);
    }
    serializer->doInt(&enumAsInt);
    if (serializer->reading())
    {
        *enumValue = static_cast<EnumBasedOnInt32>(enumAsInt);
    }
}

void serializeCrystallographicSkewData(ISerializer * serializer, MrcDensitySkewData * skewData)
{
    /* 25 | LSKFLG | signed int | 0,1
     * flag for skew matrix
     * emdb convention 0 */
    int32_t hasSkewMatrix;
    if (!serializer->reading())
    {
        if (skewData->valid_)
        {
            hasSkewMatrix = 1;
        }
        else
        {
            hasSkewMatrix = 0;
        }
    }

    serializer->doInt(&hasSkewMatrix);

    if (serializer->reading())
    {
        if (hasSkewMatrix == 1)
        {
            skewData->valid_ = true;
        }
        else
        {
            skewData->valid_ = false;
        }
    }
    /* 26-34 | SKWMAT | floating pt
     * skew matrix-S11, S12, S13, S21, S22, S23, S31, S32, S33
     * emdb convention: not set
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */

    /* 35-37 | SKWTRN | floating pt
     * skew translation-T1, T2, T3
     * emdb convention: not set
     *
     * SKWMAT, SKWTRN, and EXTRA fields are not currently used by EMDB. */
    for (auto && matrixEntry : skewData->matrix_)
    {
        float convertedValue;
        if (!serializer->reading())
        {
            convertedValue = matrixEntry * MrcDensityMapHeader::c_nmToMrcUnits;
        }
        serializer->doFloat(&convertedValue);
        if (serializer->reading())
        {
            matrixEntry = convertedValue / MrcDensityMapHeader::c_nmToMrcUnits;
        }
    }
    for (auto && translationEntry : skewData->translation_)
    {
        float convertedValue;
        if (!serializer->reading())
        {
            convertedValue = translationEntry * MrcDensityMapHeader::c_nmToMrcUnits;
        }
        serializer->doFloat(&convertedValue);
        if (serializer->reading())
        {
            translationEntry = convertedValue / MrcDensityMapHeader::c_nmToMrcUnits;
        }
    }
}

/* Supports reading according to
 *  ftp://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
 *  NOTE: http://www.ccpem.ac.uk/mrc_format/mrc2014.php differs slightly
 */
void doMrcDensityMapHeader(ISerializer * serializer, MrcDensityMapHeader * mrcFile)
{
    // 1-3 | NC, NR, NS | signed int >0 | emdb: NC=NR=NS
    // # of columns (fastest changing),rows, sections (slowest changing)
    serialize(serializer, &(mrcFile->numCrs_));

    // 4   | MODE | signed int | 0,1,2,3,4 | emdb: 2
    serializeEnumBasedOnInt32(serializer, &(mrcFile->dataMode_));

    // 5-7 | NCSTART, NRSTART, NSSTART | signed int
    // position of first column, first row, and first section
    serialize(serializer, &(mrcFile->crsStart_));

    // 8-10 | NX, NY, NZ | signed int >0 |  emdb: same as NC, NR, NS
    // intervals per unit cell repeat along X,Y Z
    serialize(serializer, &(mrcFile->extent_));

    // 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0
    // emdb Map lengths along X,Y,Z in Aangstrom
    // Length in Aangstroms for a single voxel is unit_cell_length/n_voxel
    serialize(serializer, &(mrcFile->cellLength_));

    // 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180 | emdb: 90, 90, 90
    // Unit Cell angles (degrees)following IUCr space group conventions for crystals
    serialize(serializer, &(mrcFile->cellAngles_));

    // 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)| emdb: 1, 2, 3
    // relationship of X,Y,Z axes to columns, rows, sections
    // NOTE beware offsetting: GROMACS uses X=0, whereas mrc uses X=1
    serializeOffsetPlusOne(serializer, &(mrcFile->crsToXyz_));

    // 20-22 | AMIN, AMAX, AMEAN | floating pt
    // Minimum, maximum, average density
    serializer->doFloat(&(mrcFile->dataStats_.min_));
    serializer->doFloat(&(mrcFile->dataStats_.max_));
    serializer->doFloat(&(mrcFile->dataStats_.mean_));
    // 23 | ISPG | signed int 1-230 |
    // space group number
    // For 3D volumes of single particle or tomogram entries, ISPG=1
    // For image stacks ISPG = 0
    serializeEnumBasedOnInt32(serializer, &(mrcFile->spaceGroup_));

    // 24 | NSYMBT | signed int | 80n | emdb: 0
    // # of bytes in symmetry table, expected to be multiple of 80
    int numBytesExtendedHeader;
    if (!serializer->reading())
    {
        numBytesExtendedHeader = mrcFile->extendedHeader_.size();
    }
    serializer->doInt(&numBytesExtendedHeader);
    if (!serializer->reading())
    {
        mrcFile->extendedHeader_.resize(numBytesExtendedHeader);
    }

    // 25-37 | SKEW DATA | 1 byte flag, 9 byte matrix, 3 byte translation
    serializeCrystallographicSkewData(serializer, &(mrcFile->skewData_));

    // 38-52 | EXTRA | 32 bit binary
    // 15 floats of user-defined metadata
    // EMDB might use fields 50,51 and 52 for setting the coordinate system origin
    for (auto && userFloat : mrcFile->userDefinedFloat_)
    {
        serializer->doFloat(&userFloat);
    }

    // 53 | MAP | ASCII char | emdb: "MAP "
    // MRC/CCP4 MAP format identifier
    std::string formatIdentifierString(mrcFile->formatIdentifier_.data(), mrcFile->formatIdentifier_.size());
    serializer->doString(&formatIdentifierString);
    for (size_t i = 0; i <  mrcFile->formatIdentifier_.size(); ++i)
    {
        mrcFile->formatIdentifier_[i] = formatIdentifierString[i];
    }
    // 54 | MACHST | 32 bit
    // binary machine stamp
    // MACHST written/read as 4 hex byte sequence
    // 0x44,0x41,0x00,0x00  for little endian machines
    // 0x11,0x11,0x00,0x00  for big endian machines
    serializeEnumBasedOnInt32(serializer, &(mrcFile->machineStamp_));

    // 55 | RMS | floating pt
    // Density root-mean-square deviation
    serializer->doFloat(&(mrcFile->dataStats_.rms_));

    // 56 | NLABL | signed int | 0-10
    // Number of labels
    serializer->doInt(&(mrcFile->labels_.numUsedLabels_));

    // 57-256 | LABEL_N | ASCII char | emdb : “::::EMDataBank.org::::EMD-1234::::”.
    // 10 user-defined labels each 80 characters
    for (auto && label : mrcFile->labels_.labels_)
    {
        if (serializer->reading())
        {
            label.resize(CrystallographicLables::c_labelSize);
        }
        serializer->doString(&label);
    }

    // 257-257+NSYMBT | user defined extended header information | emdb : none
    serializer->doString(&(mrcFile->extendedHeader_));
}
}   // namespace

void serializeMrcDensityMapHeader(ISerializer              * serializer,
                                  const MrcDensityMapHeader &mrcFile)
{
    MrcDensityMapHeader mrcHeaderCopy = mrcFile;
    doMrcDensityMapHeader(serializer, &mrcHeaderCopy);
}

MrcDensityMapHeader deserializeMrcDensityMapHeader(ISerializer * serializer)
{
    MrcDensityMapHeader mrcHeader;
    doMrcDensityMapHeader(serializer, &mrcHeader);
    return mrcHeader;
}

} // namespace gmx
