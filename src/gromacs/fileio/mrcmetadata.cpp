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
 * Implements mrc header.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#include "mrcmetadata.h"

#include "gromacs/math/griddata/grid.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/container/containermeasure.h"

namespace gmx
{

const float MrcMetaData::nmToMrcUnits = 10.0f;

void MrcMetaData::setEMDBDefaults()
{
    {
        swap_bytes               = false;
        space_group              = 1;
        mrc_data_mode            = 2;
        num_bytes_extened_header = 0;
        has_skew_matrix          = false;
        crs_start                = {{0, 0, 0}};
        crs_to_xyz               = {{0, 1, 2}};
        skew_matrix              = {{0, 0, 0, 0, 0, 0, 0, 0, 0}};
        skew_translation         = {0, 0, 0};
        is_crystallographic      = false;
        extra                    = {};
        extraskew                = {};
        format_identifier        = "MAP ";

        // check endianess
        #ifdef GMX_INTEGER_BIG_ENDIAN
        machine_stamp            = 1145110528;
        #else
        machine_stamp            = 4369;
        #endif
        std::string empty80CharLabel = {"                                                                                "};
        std::string emdbCustomLabel  = {"::::EMDataBank.org::::EMD-xxxx::::Own Data Following EMDB convention::::::::::::"};
        labels[0] = emdbCustomLabel;
        for (int i = 1; i < 10; ++i)
        {
            labels[i] = empty80CharLabel;
        }
        num_labels               = 1;
        extended_header          = {};
    }
}

namespace
{

constexpr int lineWidth = 82;

std::string decorateCommentString(std::string commentString)
{
    const int         spacesWidth = lineWidth - commentString.length();
    const std::string spaces      = spacesWidth > 0 ? std::string(spacesWidth, ' ') : std::string();
    return "\n" +  spaces + commentString + "\n";
}

std::string decorateDataItemString(std::string dataItemText)
{
    return " " + dataItemText + " : ";
}

template <class T>
std::string dataItemString(T item)
{
    return std::to_string(item) + " ";
}

template <class T>
std::string dataVectorString(T dataVector, size_t linePosition)
{
    const std::string newlineString("\n     ");
    std::string       vectorString;
    for (const auto &vectorElement : dataVector)
    {
        const auto newItem = dataItemString(vectorElement);
        linePosition += newItem.length();
        if (linePosition > lineWidth)
        {
            vectorString += newlineString;
            linePosition  = newlineString.size()-1;
        }
        vectorString += newItem;
    }
    return vectorString;
}

std::string sectionString( std::string sectionTitle)
{
    const auto titleLength = sectionTitle.length();
    const int  leftWidth   = (lineWidth - titleLength)/2;
    const int  rightWidth  = lineWidth - leftWidth - titleLength;
    if (leftWidth > 0 && rightWidth > 0)
    {
        return std::string(leftWidth, '-') + sectionTitle + std::string(rightWidth, '-') + "\n";
    }
    return sectionTitle + "\n";
}

std::string printBoolItem(std::string dataItemText, bool item, std::string commentString)
{
    return decorateDataItemString(dataItemText) + (item ? std::string("yes") : std::string("no")) + decorateCommentString(commentString);
}

template<class T>
std::string printItem (std::string dataItemText, T item, std::string commentString)
{
    return decorateDataItemString(dataItemText) + dataItemString(item) + decorateCommentString(commentString);
}

template<class T>
std::string print3DVec (std::string dataItemText,  BasicVector<T> item, std::string commentString)
{
    return decorateDataItemString(dataItemText) + dataItemString(item[XX]) + dataItemString(item[YY])+dataItemString(item[ZZ])+ decorateCommentString(commentString);
}


template<class T>
std::string printVector (std::string dataItemText, T item, std::string commentString)
{
    auto dataItemString = decorateDataItemString(dataItemText);
    return dataItemString + dataVectorString(item, dataItemText.length()) + decorateCommentString(commentString);
}

}

std::string MrcMetaData::to_string() const
{
    std::string result;

    result += sectionString("Grid Properties");
    result += printVector("extend",
                          extend, "grid extend in x, y, z");
    result += print3DVec("cell length [nm]",
                         cell_length, "length of cell, not all of it need be filled with voxel data");
    result += print3DVec("cell angles [deg]",
                         cell_angles, "cell angles");
    result += printVector("start column row section",
                          crs_start, "start of values in grid, typically 0, 0, 0");
    result += printVector("used columns rows sections",
                          num_crs, "voxel data will have the size of columns * rows * sections");
    result += printVector("axis order",
                          crs_to_xyz, "column, row, section correspond to these dimensions x, y, z");

    result += sectionString("Grid Data Statistics");
    result += printItem("minimum",
                        min_value, "for scaling values in compressed data mode '0'");
    result += printItem("maximum",
                        max_value, "");
    result += printItem("mean",
                        mean_value, "mean value");
    result += printItem("rms",
                        rms_value, "root mean square deviation");

    result += sectionString("Meta Information");
    result += printBoolItem("swap bytes",
                            swap_bytes, "swap bytes upon reading/writing (file and machine endianess differs)");
    result += printItem("data mode",
                        mrc_data_mode, "may be 0-4, only 2 is supported (32-bit float)");
    result += printItem("machine stamp",
                        machine_stamp, "endianess of writing "
                        "architecture - big " + std::to_string(0x44411111) +
                        " little "+ std::to_string(0x11111444));
    result += printVector("format identifier",
                          format_identifier, "four 1-byte chars 77 65 80 32 ('MAP ')");
    result += printItem("number of labels",
                        num_labels, "number of crystallographic labels");
    if (!labels.empty())
    {
        result += decorateDataItemString("labels");
        for (const auto label : labels)
        {
            result += "\n' " + label + "'";
        }
        result += "\n";
        result +=  decorateCommentString("crystallographic labels or "
                                         "'::::EMDataBank.org::::EMD-1234::::' for EMDB entries");
    }
    result += printBoolItem("is crystallographic",
                            is_crystallographic, "true if additional crystallographic data is to be read");
    result += printBoolItem("has skew matrix",
                            has_skew_matrix, "true if skew matrix is stored (only crystallographic data)");
    result += printVector("skew matrix",
                          skew_matrix, "skew matrix or data in place");
    result += print3DVec("skew translation: ",
                         skew_translation, "skew translatation or, data in place");
    result += printItem("extended header size [bytes]",
                        num_bytes_extened_header, "the size of the symbol table in bytes" );
    result += printVector("textended header",
                          extended_header, "extended header, usually symbol tables");
    result += printVector("extra skew",
                          extraskew, "skew matrix and translation for crystallographic data, unused in EMDB");
    result += printVector("extra",
                          extra, "extra data in header, unused");

    return result;
}

offset<DIM> MrcMetaData::to_xyz_order(const std::array<int, 3> &i_crs)
{
    offset<DIM> i_xyz;

    i_xyz[crs_to_xyz[XX]] = i_crs[XX];
    i_xyz[crs_to_xyz[YY]] = i_crs[YY];
    i_xyz[crs_to_xyz[ZZ]] = i_crs[ZZ];

    return i_xyz;
}

std::array<int, 3> MrcMetaData::to_crs_order(const std::array<int, 3> &i_xyz)
{
    std::array<int, 3> i_crs;

    i_crs[XX] = i_xyz[crs_to_xyz[XX]];
    i_crs[YY] = i_xyz[crs_to_xyz[YY]];
    i_crs[ZZ] = i_xyz[crs_to_xyz[ZZ]];

    return i_crs;
}



MrcMetaData &MrcMetaData::fromGrid(const IGrid<DIM> &grid )
{
    setEMDBDefaults();
    auto index_of_origin = grid.coordinateToFloorMultiIndex({{1e-6, 1e-6, 1e-6}});
    crs_start   = {{static_cast<int>(-index_of_origin[XX]), static_cast<int>(-index_of_origin[YY]), static_cast<int>(-index_of_origin[ZZ])}};
    extend      = {{static_cast<int>(grid.lattice()[XX]), static_cast<int>(grid.lattice()[YY]), static_cast<int>(grid.lattice()[ZZ])}};
    num_crs     = to_crs_order(extend);
    cell_angles = { 90, 90, 90 };

    for (int dimension = 0; dimension <= ZZ; ++dimension)
    {
        cell_length[dimension] = nmToMrcUnits * grid.cell()[dimension];
    }
    num_labels = labels.size();
    return *this;
}



std::unique_ptr < IGrid < DIM>> MrcMetaData::toGrid()
{
    RVec                gridCellLength;
    const float         mrcUnitsToNm = 1/nmToMrcUnits;
    svmul(mrcUnitsToNm, cell_length, gridCellLength);
    Grid<DIM>::NdVector cellLength {{
                                        gridCellLength[XX], gridCellLength[YY], gridCellLength[ZZ]
                                    }};

    offset<DIM> gridExtend;
    if (std::all_of(std::begin(num_crs), std::end(num_crs), [](int i){return i > 0; }))
    {
        gridExtend = to_xyz_order(num_crs);
    }
    else
    {
        if (std::all_of(std::begin(extend), std::end(extend), [](int i){return i > 0; }))
        {
            gridExtend = to_xyz_order(extend);
        }
        else
        {
            return nullptr;
        }
    }

    Grid<DIM> gridWithoutTranslation(cellLength, {extend[XX], extend[YY], extend[ZZ]});
    gridWithoutTranslation.setLatticeAndRescaleCell({gridExtend[XX], gridExtend[YY], gridExtend[ZZ]});

    Grid<DIM>::MultiIndex indexOfStartVoxel {
        crs_start[crs_to_xyz[XX]], crs_start[crs_to_xyz[YY]], crs_start[crs_to_xyz[ZZ]]
    };

    auto translation =  gridWithoutTranslation.multiIndexToCoordinate(indexOfStartVoxel);

    GridWithTranslation<DIM>  result(gridWithoutTranslation.cell(), gridWithoutTranslation.lattice(), translation);
    /* If this the map origin is shifted, because the grid indexing starts at other values than zero,
     * values read here are ignored.
     * Convention is not clear at this point, whether the translation here should be treated as extra shift,
     * in observed cases this was not the case.
     *
     * Silently ignore if the map translation due to grid-index start shift
     * does not match the shift reported here.
     */
    if (!is_crystallographic && crs_start[XX] == 0 && crs_start[YY] == 0 && crs_start[ZZ] == 0)
    {
        result.setTranslation({{(float) mrcUnitsToNm * *(extra.end()-3), (float) mrcUnitsToNm * *(extra.end()-2), (float) mrcUnitsToNm * *(extra.end()-1)}});
    }

    return std::unique_ptr < IGrid < DIM>>(new GridWithTranslation<DIM>(result));
}

MrcMetaData &MrcMetaData::set_grid_stats(const GridDataReal3D::container_type &grid_data)
{
    min_value  = containermeasure::minIgnoreNaN(grid_data);
    max_value  = containermeasure::maxIgnoreNaN(grid_data);
    mean_value = containermeasure::meanIgnoreNaN(grid_data);
    rms_value  = containermeasure::rmsIgnoreNaN(grid_data);
    return *this;
}


}      // namespace gmx
