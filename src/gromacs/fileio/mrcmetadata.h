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
/*! \file
 * \brief
 * Data structure to hold the complete mrc-file metadata as specified in the mrc format.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_MRCMETADATA_H_
#define GMX_FILEIO_MRCMETADATA_H_
#include <array>
#include <string>
#include <vector>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/griddata.h"

namespace gmx
{

template<int N> class IGrid;

/*! \brief
 * A container for the metadata in mrc file formats (compatible with ccp4 and map and mostly imod).
 *
 * For a detailed decription see
 * "EMDB Map Distribution Format Description Version 1.01 (c) emdatabank.org 2014"
 */
struct MrcMetaData{
    bool                       swap_bytes;               //!< swap bytes upon reading/writing (applied, when endianess is different between file and machine architecture)
    int                        space_group;              //!< space group as defined by IUCr conventions (Table 12.3.4.1 Standard space-group symbols, pages 824-831, International Tables for Crystallography, Volume A, fifth edition)
    /*!\brief The mrc standard defines modes 0-4.
     *
     * MODE = 0: 8 bits, density stored as a signed byte (range -128 to 127, ISO/IEC 10967)
     * MODE = 1: 16 bits, density stored as a signed integer (range -32768 to 32767, ISO/IEC 10967)
     * MODE = 2: 32 bits, density stored as a floating point number (IEEE 754)
     * MODE = 3: 32 bits, Fourier transform stored as complex signed integers (ISO/IEC 10967)
     * MODE = 4: 64 bits, Fourier transform stored as complex floating point numbers (IEEE 754)
     */
    enum class                 MrcDataMode : int { uInt8 = 0, int16 = 1, float32 = 2, complexInt32 = 3, complexFloat64 = 4 };
    int                        mrc_data_mode;            //!< data mode, currently only mode 2 is supported (32-bit float real values)
    int                        machine_stamp;            //!< endianess of map writing architecture (big endian: 0x44410000 , little endian: 0x11110000)
    std::string                format_identifier;        //!< for all density formats: four 1-byte chars reading "MAP " (a little pointless, I know)

    int                        num_labels;               //!< number of used crystallographic labels, 0 for imagestacks, 1 for emdb data
    std::vector < std::string> labels;                   //!< crystallographic labels or ::::EMDataBank.org::::EMD-1234:::: for EMDB entries

    RVec                       cell_length;              //!< length of the crystallographic unit cell
    RVec                       cell_angles;              //!< crystallographic unit cell angles

    std::array<int, 3>         crs_to_xyz;               //!< Axis order
    std::array<int, 3>         num_crs;                  //!< redundand entry, we use the grid extend (NX,NY,NZ) from header words 8-10
    std::array<int, 3>         extend;                   //!< the grid extend, check against num_crs
    std::array<int, 3>         crs_start;                //!< Start of values in grid, typically 0,0,0

    float                      min_value;                //!< minimum voxel value. may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                      max_value;                //!< maximum voxel value. may be used to scale values in currently unsupported compressed data mode (mrc_data_mode=0)
    float                      mean_value;               //!< mean voxel value   (not always reported,as evident from density)
    float                      rms_value;                //!< rms of the density (not always reported,as evident from density)

    bool                       is_crystallographic;      //!< true if crystallographic data is to be read
    bool                       has_skew_matrix;          //!< only crystallographic data: true if skew matrix is stored
    std::array<float, 9>       skew_matrix;              //!< only crystallographic data: skew matrix or, if skew flag is zero, data in place of skew matrix
    RVec                       skew_translation;         //!< only crystallographic data: skew translatation or, if skew flag is zero, data in place of skew translation
    int                        num_bytes_extened_header; //!< only crystallographic data: the size of the symbol table in bytes
    std::vector<char>          extended_header;          //!< only crystallographic data: extended header, usually symbol tables

    std::array<float, 13>      extraskew;                //!< fields unused in EMDB standard, but used for skew matrix and translation in crystallogrphic data (skew flag, skew matrix and skew translation)
    std::array<float, 15>      extra;                    //!< extra data in header, currently unused

    const static float         nmToMrcUnits;             //!< Conversion factor from nm to mrc units (Angstrom)
    /*! \brief
     * Set values to emdb defauts.
     */
    void setEMDBDefaults();
    /*! \brief
     * Print contents of the metadata to string.
     *
     * \returns string describing metadata
     */
    std::string  to_string() const;
    MrcMetaData &fromGrid(const IGrid<DIM> &grid);
    std::unique_ptr < IGrid < DIM>> toGrid();

    std::array<int, 3> to_crs_order(const std::array<int, 3> &xyz_order);
    std::array<int, 3> to_xyz_order(const std::array<int, 3> &i_crs);

    MrcMetaData       &set_grid_stats(const GridDataReal3D::container_type &grid_data);

};

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_MRCMETADATA_H_ */
