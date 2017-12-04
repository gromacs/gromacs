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
 * Data structure mirroring the content of an xplor-file.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 * \ingroup module_fileio
 */

#ifndef GMX_FILEIO_XPLOR_H_
#define GMX_FILEIO_XPLOR_H_
#include <array>
#include <string>
#include <vector>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/math/griddata/griddata.h"


namespace gmx
{

struct t_mapdata;
/*! \brief
 * A container for the metadata in xplor file formats.
 *
 * For a detailed decription see
 * http://www.scripps.edu/rc/softwaredocs/msi/xplor981/formats.html
 */
struct XplorData{
    XplorData() = default;
    XplorData(FILE* inputStream);

    void                  maskData(const std::vector<float> &mask, float threshold);
    void                  write(FILE * outputStream) const;

    RVec                  cell_length;      //!< length of the crystallographic unit cell
    RVec                  cell_angles;      //!< crystallographic unit cell angles

    std::array<int, DIM>  extend;           //!< the grid extend, check against num_crs
    std::array<int, DIM>  crs_start;        //!< Start of values in grid, typically 0,0,0
    std::array<int, DIM>  crs_end;          //!< End of values in grid, typically extend[XX],extend[YY],extend[ZZ]

    float                 mean_value;       //!< mean voxel value   (not always reported,as evident from density)
    float                 rms_value;        //!< rms of the density (not always reported,as evident from density)

    std::vector<float>    data;             //!< The density data

    static const float    nmToXplorUnits;   //!< Multiply by this factor to obtain nm from xplor units (Aangstrom)
    static const int      magicCheckNumber; //!< Typically -9999; establishes sucessful file reading
};

}      // namespace gmx
#endif /* end of include guard: GMX_FILEIO_XPLOR_H_ */
