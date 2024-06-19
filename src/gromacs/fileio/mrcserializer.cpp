/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements methods from mrcserializer.h
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "mrcserializer.h"

#include "config.h"

#include <cstdint>

#include <array>
#include <type_traits>
#include <vector>

#include "gromacs/fileio/mrcdensitymapheader.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/iserializer.h"

namespace gmx
{

namespace
{

/*! \brief Machine stamp to indicate endianess of mrc/ccp4 file.
 * As named in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 * The first two bytes of a 4 byte unsigned int are used to indicate endianess:
 * 0x11 0x11 big endian
 * 0x44 0x44 small endian
 * Byte-swap data if appropriate, when transferring data files between machines.
 */
enum class MachineStamp : int32_t
{
    bigEndian   = 0x11110000, //!< big endian magic number 0x11 0x11 0x00 0x00 = 286,326,784
    smallEndian = 0x44440000, //!< small endian magic number 0x44 0x44 0x00 0x00 = 1,145,307,136
};

/*! \brief Serialize a container of int32_t values.
 * Serializes all containers with value_type int32_t that may looped over in a
 * range based for loop and have modifiable elements.
 *
 * \tparam ContainerType type of container to be serialized
 * \param[in,out] serializer the serializer
 * \param[in,out] valueContainer the array to be serialized
 */
template<typename ContainerType>
std::enable_if_t<std::is_same_v<typename ContainerType::value_type, int32_t>, void>
serialize(ISerializer* serializer, ContainerType* valueContainer)
{
    for (auto& value : *valueContainer)
    {
        serializer->doInt32(&value);
    }
}

/*! \brief Serialize a container of float values.
 * Serializes all containers with value_type float that may looped over in a
 * range based for loop and have modifiable elements.
 *
 * \tparam ContainerType type of container to be serialized
 * \param[in,out] serializer the serializer
 * \param[in,out] valueContainer the array to be serialized
 */
template<typename ContainerType>
std::enable_if_t<std::is_same_v<typename ContainerType::value_type, float>, void>
serialize(ISerializer* serializer, ContainerType* valueContainer)
{
    for (auto& value : *valueContainer)
    {
        serializer->doFloat(&value);
    }
}

//! Serialize and convert from FORTRAN 1-based to C 0-based indices when reading and vice versa when writing
void serializeIndex(ISerializer* serializer, int32_t* index)
{
    int32_t fortranIndex;
    if (!serializer->reading())
    {
        fortranIndex = *index + 1;
    }
    serializer->doInt32(&fortranIndex);
    if (serializer->reading())
    {
        *index = fortranIndex - 1;
    }
}

/*! \brief
 * Serializes an integer array and add unity when writing, substracting unity when reading.
 */
void serializeIndices(ISerializer* serializer, std::array<int32_t, 3>* valueArray)
{
    for (auto& value : *valueArray)
    {
        serializeIndex(serializer, &value);
    }
}

/*! \brief Serialize input as int32_t via static casting.
 * \tparam IntegralType type to be serialized as int32_t
 */
template<class IntegralType>
std::enable_if_t<(std::is_integral_v<IntegralType> || std::is_enum_v<IntegralType>), void>
serializeAsInt32(ISerializer* serializer, IntegralType* value)
{
    int32_t serializedValue;
    if (!serializer->reading())
    {
        serializedValue = static_cast<int32_t>(*value);
    }
    serializer->doInt32(&serializedValue);
    if (serializer->reading())
    {
        *value = static_cast<IntegralType>(serializedValue);
    }
}

//! conversion constant from nm to MRC distance units (Ångström)
constexpr float c_nmToMrcUnits = 10;

//! Convert MRC distances to nm
float mrcUnitsToNm(float mrcValue)
{
    return mrcValue / c_nmToMrcUnits;
}

//! Convert nm to MRC distances
float nmToMrcUnits(float nmValue)
{
    return nmValue * c_nmToMrcUnits;
}

//! Serialize and convert between MRC and GROMACS distance units
void serializeDistance(ISerializer* serializer, float* distance)
{
    float convertedDistance;
    if (!serializer->reading())
    {
        convertedDistance = nmToMrcUnits(*distance);
    }
    serializer->doFloat(&convertedDistance);
    if (serializer->reading())
    {
        *distance = mrcUnitsToNm(convertedDistance);
    }
}

//! Serialize the skew data, words 25-37 in an mrc file.
void serializeCrystallographicSkewData(ISerializer* serializer, MrcDensitySkewData* skewData)
{
    /* 25 | LSKFLG | signed int | emdb: 0 or 1
     * flag for skew matrix
     */
    serializeAsInt32(serializer, &(skewData->valid_));
    /* 26-34 | SKWMAT | floating pt | emdb: not set
     * skew matrix-S11, S12, S13, S21, S22, S23, S31, S32, S33
     * 35-37 | SKWTRN | floating pt | emdb: not set
     * skew translation-T1, T2, T3
     */
    for (auto& matrixEntry : skewData->matrix_)
    {
        serializeDistance(serializer, &matrixEntry);
    }
    for (auto& translationEntry : skewData->translation_)
    {
        serializeDistance(serializer, &translationEntry);
    }
}

/*! \brief Symmetrise and de-serializes the mrc density map header.
 *
 * Supports reading EMDB density maps in mrc format according to
 * ftp://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
 * \note format has small differences to http://www.ccpem.ac.uk/mrc_format/mrc2014.php
 */
void doMrcDensityMapHeader(ISerializer* serializer, MrcDensityMapHeader* mrcFile)
{
    // 1-3 | NC, NR, NS | signed int >0 | emdb: NC=NR=NS
    // # of columns (fastest changing),rows, sections (slowest changing)
    serialize(serializer, &(mrcFile->numColumnRowSection_));

    // 4   | MODE | signed int | 0,1,2,3,4 | emdb: 2
    serializeAsInt32(serializer, &(mrcFile->dataMode_));

    // 5-7 | NCSTART, NRSTART, NSSTART | signed int |  position of first column, first row, and first section
    serialize(serializer, &(mrcFile->columnRowSectionStart_));

    // 8-10 | NX, NY, NZ | signed int >0 |  emdb: same as NC, NR, NS | intervals per unit cell repeat along X,Y Z
    serialize(serializer, &(mrcFile->extent_));

    // 11-13 | X_LENGTH, Y_LENGTH, Z_LENGTH | floating pt >0 | emdb Map lengths along X,Y,Z in
    // Ångström Length in Ångström for a single voxel is unit_cell_length/n_voxel
    serialize(serializer, &(mrcFile->cellLength_));

    // 14-16 | ALPHA,BETA,GAMMA | floating pt >0, <180 | emdb: 90, 90, 90
    // Unit Cell angles (degrees)following IUCr space group conventions for crystals
    serialize(serializer, &(mrcFile->cellAngles_));

    // 17-19 | MAPC, MAPR, MAPS | signed int | 1 (=X) 2 (=Y) 3 (=Z)| emdb: 1, 2, 3
    // relationship of X,Y,Z axes to columns, rows, sections
    // GROMACS uses C-style X=0, whereas mrc uses FORTRAN-style X=1
    serializeIndices(serializer, &(mrcFile->columnRowSectionToXyz_));

    // 20-22 | AMIN, AMAX, AMEAN | floating pt |  Minimum, maximum, average density
    serializer->doFloat(&(mrcFile->dataStatistics_.min_));
    serializer->doFloat(&(mrcFile->dataStatistics_.max_));
    serializer->doFloat(&(mrcFile->dataStatistics_.mean_));

    // 23 | ISPG | signed int 1-230 | emdb: 1; 0 for image stacks | space group number
    serializeAsInt32(serializer, &(mrcFile->spaceGroup_));

    // 24 | NSYMBT | signed int | 80n | emdb: 0
    // # of bytes in symmetry table, expected to be multiple of 80
    int32_t numBytesExtendedHeader;
    if (!serializer->reading())
    {
        numBytesExtendedHeader = mrcFile->extendedHeader_.size();
    }
    serializer->doInt32(&numBytesExtendedHeader);
    if (serializer->reading())
    {
        mrcFile->extendedHeader_.resize(numBytesExtendedHeader);
    }

    // 25-37 | SKEW DATA | 1 byte flag, 9 byte matrix, 3 byte translation
    serializeCrystallographicSkewData(serializer, &(mrcFile->skewData_));

    // 38-52 | EXTRA | 32 bit binary
    // 15 floats of user-defined metadata
    // EMDB might use fields 50,51 and 52 for setting the coordinate system origin
    for (auto& userFloat : mrcFile->userDefinedFloat_)
    {
        serializer->doFloat(&userFloat);
    }

    // 53 | MAP | ASCII char | emdb: "MAP "
    // MRC/CCP4 MAP format identifier
    for (auto& formatIdentifierCharacter : mrcFile->formatIdentifier_)
    {
        serializer->doUChar(&formatIdentifierCharacter);
    }

    // 54 | MACHST | 32 bit | binary machine stamp written/read as 4 hex byte sequence
    MachineStamp machineStamp = GMX_INTEGER_BIG_ENDIAN ? MachineStamp::bigEndian : MachineStamp::smallEndian;
    serializeAsInt32(serializer, &machineStamp);

    // 55 | RMS | floating pt | Density root-mean-square deviation
    serializer->doFloat(&(mrcFile->dataStatistics_.rms_));

    // 56 | NLABL | signed int | 0-10 | Number of used user defined labels
    serializer->doInt32(&(mrcFile->labels_.numUsedLabels_));

    // 57-256 | LABEL_N | ASCII char | emdb : “::::EMDataBank.org::::EMD-1234::::”.
    // 10 user-defined labels each 80 characters
    for (auto&& label : mrcFile->labels_.labels_)
    {
        for (auto&& labelCharacter : label)
        {
            serializer->doUChar(&labelCharacter);
        }
    }

    // 257-257+NSYMBT | user defined extended header information | emdb : none
    for (auto& extendedHeaderCharacter : mrcFile->extendedHeader_)
    {
        serializer->doUChar(&extendedHeaderCharacter);
    }
}

} // namespace

void serializeMrcDensityMapHeader(ISerializer* serializer, const MrcDensityMapHeader& mrcFile)
{
    MrcDensityMapHeader mrcHeaderCopy = mrcFile;
    doMrcDensityMapHeader(serializer, &mrcHeaderCopy);
}

MrcDensityMapHeader deserializeMrcDensityMapHeader(ISerializer* serializer)
{
    MrcDensityMapHeader mrcHeader;
    doMrcDensityMapHeader(serializer, &mrcHeader);
    return mrcHeader;
}

} // namespace gmx
