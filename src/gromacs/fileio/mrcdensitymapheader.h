/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Implement mrc/ccp4-file metadata.
 *
 * \author Christian Blau <cblau@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_MRCDENSITYMAPHEADER_H
#define GMX_FILEIO_MRCDENSITYMAPHEADER_H
#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace gmx
{

/*! \libinternal \brief
 * Space group in three dimensions.
 *
 * Currently only "no symmetry" is supported, the complete enum class would hold
 * 230 symbols.
 *
 * Table 12.3.4.1 Standard space-group symbols, pages 824-831,
 * International Tables for Crystallography, Volume A, fifth edition
 */
enum class SpaceGroup : int32_t
{
    P1 = 1, //!< no symmetry
};

/*! \libinternal \brief
 * The type of density data stored in an mrc file.
 *
 * Modes 0-4 are defined by the standard.
 * NOTE only mode 2 is currently implemented and used.
 */
enum class MrcDataMode : int32_t
{
    uInt8          = 0, //!< compressed data mode, 8 bits, signed byte (range -128 to 127, ISO/IEC 10967)
    int16          = 1, //!< 16 bits, signed integer (range -32768 to 32767, ISO/IEC 10967)
    float32        = 2, //!< 32 bits, floating point number (IEEE 754)
    complexInt32   = 3, //!< 32 bits, complex signed integers (ISO/IEC 10967)
    complexFloat64 = 4, //!< 64 bits, complex floating point numbers (IEEE 754)
};

/*! \libinternal \brief
 * Machine stamp to indicate endianess of mrc/ccp4 file.
 *
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

/*! \libinternal 
 * \brief Statistics about data arrays.
 */
struct MrcDataStatistics
{
    MrcDataStatistics();
    float min_;  //!< Minimum voxel value scales values in (currently unsupported) compressed data mode.
    float max_;  //!< Maximum data value scales values in (currently unsupported) compressed data mode.
    float mean_; //!< mean of the data (evident from density)
    float rms_;  //!< rms of the data (evident from density)
};

/*! \libinternal 
 * \brief Skew matrix and translation.
 */
struct MrcDensitySkewData
{
    MrcDensitySkewData();
    bool                     valid_;       //!< True if skew matrix is stored.
    std::array<float, 3 * 3> matrix_;      //!< Skew matrix for crystallographic unit cell.
    std::array<float, 3>     translation_; //!< Translatation of crystallographic unit cell.
};

/*! \libinternal 
 * \brief Crystallographic labels for mrc data.
 */
struct MrcCrystallographicLables
{
    MrcCrystallographicLables();
    //! Wether used, or not, the mrc format always employs ten lines for crystallographic labels.
    static constexpr int c_numCrystallographicLabels = 10;
    static constexpr int c_labelSize                 = 80; //!< Length of crystallographic labels is eighty.

    //! Number of used crystallographic labels, 0 for imagestacks, 1 for emdb data
    int32_t                                              numUsedLabels_;
    //! Crystallographic labels or "::::EMDataBank.org::::EMD-1234::::" for EMDB entries
    std::array<std::string, c_numCrystallographicLabels> labels_;
};

/*! \libinternal \brief
 * A container for the data in mrc density map file formats.
 *
 * Mrc files are a widely used file format in crystallography and cryo electron
 * microscopy to represent volumetric data. Covers ccp4, map and imod.
 *
 * For a detailed decription see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcDensityMapHeader{

    //! Sets values to emdb defauts.
    explicit MrcDensityMapHeader();

    static constexpr float                     c_nmToMrcUnits = 10.0f;       //<! Conversion factor from nm to mrc units (Angstrom)

    SpaceGroup                                 spaceGroup_;                  //!< Space group of stored data.
    MrcDataMode                                dataMode_;                    //<! Data mode, currently only mode 2 is supported (32-bit float real values)
    MachineStamp                               machineStamp_;                //<! //! Machine stamp to determine endianess of map writing architecture

    static constexpr int                       c_formatIdentifierLength = 4; //!< Length of the format identifier tag
    std::array<char, c_formatIdentifierLength> formatIdentifier_;            //! Expected to be "MAP "

    static constexpr int                       c_numUserDefinedFloats = 15;  //!< Number of user defined float numbers
    std::array<float, c_numUserDefinedFloats>  userDefinedFloat_;            //!< 15 freely defined float numbers

    MrcCrystallographicLables                  labels_;                      //!< Labels for crystallographic data

    std::array<float, 3>                       cellLength_;                  //!< Length of the crystallographic unit cell
    std::array<float, 3>                       cellAngles_;                  //!< crystallographic unit cell angles

    //! Axis order of stored data. Data is stored with colums varying the fastest, rows second and sections the slowest.
    std::array<int32_t, 3> crsToXyz_;

    std::array<int32_t, 3> numCrs_;         //!< Column, row and section count determines the number of data points.
    std::array<int32_t, 3> crsStart_;       //!< Start of values in grid
    std::array<int32_t, 3> extent_;         //!< The number of grid points in the crystall cell

    MrcDataStatistics      dataStats_;      //!< Statistics about the data stored in the file.
    MrcDensitySkewData     skewData_;       //!< Data to perform crystallographic unit cell skewing
    std::string            extendedHeader_; //!< Extended header with symbol tables.
};

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_MRCDENSITYMAPHEADER_H */
