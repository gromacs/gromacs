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
/*! \file
 * \brief
 * Implement mrc/ccp4-file metadata.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_MRCFILE_H
#define GMX_FILEIO_MRCFILE_H
#include <array>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

namespace gmx
{

/*! \brief
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
    P1 = 1, //< no symmetry
};

/*! \brief
 * The type of density data stored in an mrc file.
 *
 * Modes 0-4 are defined by the standard.
 * NOTE only mode 2 is currently implemented and used.
 */
enum class MrcDataMode : int32_t
{
    uInt8          = 0, //< compressed data mode, 8 bits, signed byte (range -128 to 127, ISO/IEC 10967)
    int16          = 1, //< 16 bits, signed integer (range -32768 to 32767, ISO/IEC 10967)
    float32        = 2, //< 32 bits, floating point number (IEEE 754)
    complexInt32   = 3, //< 32 bits, complex signed integers (ISO/IEC 10967)
    complexFloat64 = 4, //< 64 bits, complex floating point numbers (IEEE 754)
};

/*! \brief
 * Machine stamp to indicate endianess of mrc/ccp4 file.
 *
 * The first two bytes of a 4 byte unsigned int are used to indicate endianess:
 * 0x11 0x11 big endian
 * 0x44 0x44 small endian
 * Byte-swap data if appropriate, when transferring data files between machines.
 */
enum class MachineStamp : int32_t
{
    bigEndian   = 0x11110000, //< big endian magic number 0x11 0x11 0x00 0x00 = 286,326,784
    smallEndian = 0x44440000, //< small endian magic number 0x44 0x44 0x00 0x00 = 1,145,307,136
};

/*! \brief
 * A container for the data in mrc density map file formats.
 *
 * Mrc files are a widely used file format in crystallography and cryo electron
 * microscopy to represent volumetric data. Covers ccp4, map and imod.
 *
 * For a detailed decription see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcDensityMapHeader{

    //! ensure 32 bit integers as demanded by the mrc format
    using mrc_integer = int32_t;
    //! Conversion factor from nm to mrc units (Angstrom)
    static constexpr float c_nmToMrcUnits = 10.0f;

    //! Sets values to emdb defauts.
    explicit MrcDensityMapHeader();
    //! Space group of stored data.
    SpaceGroup   spaceGroup;
    //! Data mode, currently only mode 2 is supported (32-bit float real values)
    MrcDataMode  dataMode;
    //! Machine stamp to determine endianess of map writing architecture. (big endian: 0x44410000 , little endian: 0x11110000)
    MachineStamp machineStamp;

    //! Length of the format identifier tag
    static constexpr int c_formatIdentifierLength = 4;

    //! expected for all density formats: four 1-byte chars reading "MAP "
    std::array<char, c_formatIdentifierLength> formatIdentifier;

    //! Number of user defined float numbers
    static constexpr int                      c_numUserDefinedFloats = 15;
    //! float numbers to be defined by users to their liking
    std::array<float, c_numUserDefinedFloats> userDefinedFloat;

    //! Number of used crystallographic labels, 0 for imagestacks, 1 for emdb data.
    mrc_integer          numUsedLabels;
    //! Wether used, or not, the mrc format always employs ten lines for crystallographic labels.
    static constexpr int c_numCrystallographicLabels = 10;

    //! Length of crystallographic labels is eighty.
    static constexpr int                                 c_labelSize = 80;
    //! crystallographic labels or ::::EMDataBank.org::::EMD-1234:::: for EMDB entries
    std::array<std::string, c_numCrystallographicLabels> labels;

    //! length of the crystallographic unit cell
    std::array<float, 3> cellLength;
    //! crystallographic unit cell angles
    std::array<float, 3> cellAngles;
    /*! \brief Axis order of stored data.
     * Data is stored with colums varying the fastest, rows second and sections the slowest.
     */
    std::array<mrc_integer, 3> crsToXyz;
    //! Column, row and section count determines the number of data points.
    std::array<mrc_integer, 3> numCrs;
    //! Start of values in grid, typically 0,0,0 for cryo-EM data.
    std::array<mrc_integer, 3> crsStart;
    //! the number of grid points in the crystall cell.
    std::array<mrc_integer, 3> extent;

    /*! \brief Statistics about data arrays.
     */
    struct DataStats{
        /*! \brief Minimum voxel value.
         * Scales values in (currently unsupported) compressed data mode.
         */
        float min;
        /*! \brief Maximum data value.
         * Scales values in (currently unsupported) compressed data mode.
         */
        float max;
        //! mean of the data (never required, evident from density)
        float mean;
        //! rms of the data (never required, evident from density)
        float rms;
    };
    //! Statistics about the data stored in the file.
    DataStats dataStats;

    /*! \brief Skew matrix and translation.
     */
    struct SkewData
    {
        //! True if skew matrix is stored.
        bool valid;

        //! Skew matrix for crystallographic unit cell.
        std::array<float, 3*3> matrix;

        //! Translatation of crystallographic unit cell.
        std::array<float, 3> translation;
    };
    /*! \brief Data to perform crystallographic unit cell skewing
     * NOTE the SkewData valid signalls if these values are to be ignored
     */
    SkewData skewData;

    //! Extended header with symbol tables.
    std::string extendedHeader;
};

/*! \brief Verbose result of checking the mrc file for consistency.
 */
struct MrcDiagnostic
{
    //! A fixed copy of the diagnosed header where
    MrcDensityMapHeader fixedHeader;
    //! A verbose message summarising all found inconsistencies in the MRC header.
    std::string         report;
};

/*! \brief
 * Check a mrc density map header for consistency and report warnings and a
 * fixed version with default values.
 *
 * \param[in] mrcFile the mrc header to be checked
 * \returns Diagnostic with fixed header and error report
 */
MrcDiagnostic mrcDiagnosticAndFixes(const MrcDensityMapHeader &mrcFile);

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_MRCFILE_H */
