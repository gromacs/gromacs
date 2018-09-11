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
 * Implement file data to GridDataReal3D and back conversions.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_fileio
 */
#include "griddataview.h"

#include "gromacs/math/container/containermeasure.h"
#include "gromacs/math/griddata/griddata.h"
#include "xplor.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/mrcmetadata.h"

namespace gmx
{

GridDataReal3D gridDataFromXplor(const XplorData &xplor)
{
    RVec                                resolution;
    Grid<DIM>::NdVector                 translation;
    CanonicalVectorBasis<DIM>::NdVector cellLength;
    bounds<DIM> extend;
    for (size_t i = 0; i < DIM; i++)
    {
        extend[i]      = xplor.crs_end[i] - xplor.crs_start[i] + 1;
        resolution[i]  = xplor.cell_length[i] / xplor.extend[i];
        cellLength[i]  = (extend[i]) * resolution[i];
        // origin of the lattice
        translation[i] = xplor.crs_start[i] * resolution[i];
    }

    CanonicalVectorBasis<DIM> cell(cellLength);

    GridDataReal3D            result(GridWithTranslation<DIM>(cell, extend, translation));
    std::copy(std::begin(xplor.data), std::end(xplor.data), std::begin(result));

    return result;
}

Grid<DIM> cellFromXplor(const XplorData &xplor)
{
    CanonicalVectorBasis<DIM> cell(CanonicalVectorBasis<DIM>::NdVector({{xplor.cell_length[XX], xplor.cell_length[YY], xplor.cell_length[ZZ]}}));
    return Grid<DIM>(cell, {xplor.extend[XX], xplor.extend[YY], xplor.extend[ZZ]});
}

XplorData xplorFromGridData(const GridDataReal3D &map,
                            std::array<int, DIM> cellExtend)
{
    XplorData              xplor;
    bool                   cellExtendIsValid =
        std::all_of(std::begin(cellExtend), std::end(cellExtend),
                    [](int i) { return i != 0; });
    xplor.cell_angles = {90., 90., 90};

    /* With valid cell extend, set up the unit cell from grid resolution and unit
     * cell lengths
     * The cell and the griddata from gromacs may have different extends
     */
    if (cellExtendIsValid)
    {
        std::array<float, DIM> cellLength =
            map.getGrid().unitCell().basisVectorLengths();
        std::transform(
                std::begin(cellLength), std::end(cellLength), std::begin(cellExtend),
                std::begin(cellLength),
                [](float cellLength, int nExtend) { return cellLength * nExtend; });
        std::copy(std::begin(cellLength), std::end(cellLength), std::begin(xplor.cell_length.as_vec()));
        xplor.extend = cellExtend;
    }
    else
    {
        const auto &basisVectorLengths = map.getGrid().cell().basisVectorLengths();
        std::copy(std::begin(basisVectorLengths), std::end(basisVectorLengths),
                  std::begin(xplor.cell_length.as_vec()));
        xplor.extend = {{static_cast<int>(map.getGrid().lattice()[XX]), static_cast<int>(map.getGrid().lattice()[YY]), static_cast<int>(map.getGrid().lattice()[ZZ])}};
    }

    auto originCoordinate = map.getGrid().multiIndexToCoordinate({{0, 0, 0}});
    std::transform(std::begin(originCoordinate), std::end(originCoordinate),
                   std::begin(map.getGrid().unitCell().basisVectorLengths()),
                   std::begin(xplor.crs_start),
                   [](float origin, float unitVector) { return std::round(origin/unitVector); });

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        xplor.crs_end[dimension] = xplor.crs_start[dimension] +
            map.getGrid().lattice()[dimension] - 1;
    }

    xplor.mean_value = containermeasure::mean(map);
    xplor.rms_value  = containermeasure::rms(map);

    std::copy(std::begin(map), std::end(map), std::back_inserter(xplor.data));
    return xplor;
}


MapConverter::MapConverter(std::string filename)
{
    const auto &extension = filename.substr(filename.find_last_of(".") + 1);
    if (extension == "ccp4" || extension == "map" || extension == "mrc")
    {
        grid_ = MrcFile().read(filename);
    }
    if (extension == "xplor")
    {
        auto outputfile = gmx_fio_fopen(filename.c_str(), "r");
        grid_ = gridDataFromXplor(XplorData(outputfile));
        gmx_fio_fclose(outputfile);
    }
}

void MapConverter::MapConverter::to(std::string filename)
{
    const auto &extension = filename.substr(filename.find_last_of(".") + 1);
    if (extension == "ccp4" || extension == "map" || extension == "mrc")
    {
        MrcFile().write(filename, grid_);
    }
    if (extension == "xplor" || extension == "dat")
    {
        auto outputfile = gmx_fio_fopen(filename.c_str(), "w");
        if (extension == "xplor")
        {
            xplorFromGridData(grid_).write(outputfile);
        }
        if (extension == "dat")
        {
            const auto textDump = MrcMetaData().fromGrid(grid_.getGrid()).set_grid_stats(grid_).to_string();
            fprintf(outputfile, "%s", textDump.c_str());
        }
        gmx_fio_fclose(outputfile);
    }
}

MapConverter::MapConverter(const GridDataReal3D &grid)
{
    grid_ = grid;
}

MapConverter::MapConverter() = default;

const GridDataReal3D &MapConverter::map()
{
    return grid_;
}
}
