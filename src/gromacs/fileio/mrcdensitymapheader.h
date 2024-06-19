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

#include <cstddef>
#include <cstdint>

#include <array>
#include <vector>

#include "gromacs/math/coordinatetransformation.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"

namespace gmx
{

/*! \brief Space group in three dimensions.
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

/*! \brief The type of density data stored in an mrc file.
 * As named in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 * Modes 0-4 are defined by the standard.
 * NOTE only mode 2 is currently implemented and used.
 */
enum class MrcDataMode : int32_t
{
    uInt8   = 0, //!< compressed data mode, 8 bits, signed byte (range -128 to 127, ISO/IEC 10967)
    int16   = 1, //!< 16 bits, signed integer (range -32768 to 32767, ISO/IEC 10967)
    float32 = 2, //!< 32 bits, floating point number (IEEE 754)
    complexInt32   = 3, //!< 32 bits, complex signed integers (ISO/IEC 10967)
    complexFloat64 = 4, //!< 64 bits, complex floating point numbers (IEEE 754)
};

/*! \libinternal
 * \brief Statistics about mrc data arrays.
 */
struct MrcDataStatistics
{
    float min_ = 0.; //!< Minimum data value scales values in (currently unsupported) compressed data mode.
    float max_ = 0.; //!< Maximum data value scales values in (currently unsupported) compressed data mode.
    float mean_ = 0.; //!< mean of the data
    float rms_  = 0.; //!< rms of the data
};

/*! \libinternal
 * \brief Skew matrix and translation.
 * As named in
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcDensitySkewData
{
    bool                        valid_ = false; //!< True if skew matrix is stored.
    std::array<float, DIM* DIM> matrix_ = {}; //!< Skew matrix for crystallographic unit cell in Ångström
    std::array<float, DIM> translation_ = {}; //!< Translation of crystallographic unit cell in Ångström
};

/*! \libinternal
 * \brief Crystallographic labels for mrc data.
 */
struct CrystallographicLabels
{
    static constexpr int c_labelSize = 80; //!< Length of crystallographic labels is eighty.
    //! Number of used crystallographic labels, 0 for imagestacks, 1 for emdb data
    int32_t numUsedLabels_ = 0;

    //! Crystallographic labels or "::::EMDataBank.org::::EMD-1234::::" for EMDB entries
    std::array<std::array<unsigned char, c_labelSize>, 10> labels_ = {};
};

/*! \libinternal
 * \brief A container for the data in mrc density map file formats.
 *
 * Mrc files are a widely used file format in crystallography and cryo electron
 * microscopy to represent volumetric data. Covers ccp4, map and imod.
 *
 * For a detailed description see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcDensityMapHeader
{

    //! Space group of stored data.
    SpaceGroup spaceGroup_ = SpaceGroup::P1;
    //! Data mode, currently only mode 2 is supported
    MrcDataMode dataMode_ = MrcDataMode::float32;
    //! Identifies file format, expected to be "MAP "
    std::array<unsigned char, 4> formatIdentifier_ = { { 'M', 'A', 'P', ' ' } };
    //! 15 unspecified float numbers
    std::array<float, 15> userDefinedFloat_ = {};

    //! Labels for crystallographic data
    CrystallographicLabels labels_ = {};

    //! Length of the crystallographic unit cell in Ångström
    std::array<float, DIM> cellLength_ = { { 1., 1., 1. } };
    //! crystallographic unit cell angles
    std::array<float, DIM> cellAngles_ = { { 90., 90., 90. } };

    //! Data axis order with columns varying the fastest, and sections the slowest.
    std::array<int32_t, DIM> columnRowSectionToXyz_ = { { 0, 1, 2 } };

    std::array<int32_t, DIM> numColumnRowSection_   = {}; //!< Column, row and section count
    std::array<int32_t, DIM> columnRowSectionStart_ = {}; //!< Start of values in grid
    std::array<int32_t, DIM> extent_ = {}; //!< The number of grid points in the crystall cell

    //! Statistics about the data stored in the file.
    MrcDataStatistics dataStatistics_ = {};
    //! Data to perform crystallographic unit cell skewing
    MrcDensitySkewData skewData_ = {};
    //! Extended header with symmetry tables
    std::vector<unsigned char> extendedHeader_ = {};
};

/*! \brief Return the number of density data items that are expected
 *         to follow this header.
 * \throws InternalError if the number of data items cannot be determined
 * \returns the number of voxels
 */
size_t numberOfExpectedDataItems(const MrcDensityMapHeader& header);

/*! \brief Extract the transformation into lattice coordinates.
 * \note Transformation into lattice coordinates is not treated uniformly
 *       in different implementations for the mrc format,e.g., vmd, pymol and
 *       chimera have different conventions. Following the vmd implementation here.
 *
 * In determining the density origin coordinates, explicit ORIGIN records
 * (also called origin2k) in the user defined floats 13 - 15, corresponding to
 * words 50,51 and 52 in the mrc header, precedence over ColumnRowSectionStart.
 * Only if above values are zero, using the column, row and section start to
 * determine the translation vector.
 *
 * \param[in] header from which the coordinate transformation is to be extracted
 * \returns a functor that transforms real space coordinates into the lattice
 */
TranslateAndScale getCoordinateTransformationToLattice(const MrcDensityMapHeader& header);

/*! \brief Extract the extents of the density data
 * \param[in] header from which the extents are to be extracted
 * \returns density data extents in three dimensions.
 */
dynamicExtents3D getDynamicExtents3D(const MrcDensityMapHeader& header);

/*! \brief Checks if the values in the header are sane.
 *
 * Checks extents and numbers of columns, rows and sections, as well as unit
 * cell angles for positivity and to be within bounds.
 *
 * Bounds are set generously not to hamper future creative uses of mrc files.
 *
 * \returns true if all header values are within resonable albeit generous bounds
 */
bool mrcHeaderIsSane(const MrcDensityMapHeader& header);

} // namespace gmx
#endif /* end of include guard: GMX_FILEIO_MRCDENSITYMAPHEADER_H */
