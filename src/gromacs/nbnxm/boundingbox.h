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

/*! \file
 * \internal
 *
 * \brief Declares and defines the BoundingBox class.
 *
 * Can be used for computing rectangular bounding boxes of clusters of atoms
 * and for computing an underestimate of the distance between such clusters.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXM_BOUNDINGBOX_H
#define GMX_NBNXM_BOUNDINGBOX_H

#include <algorithm>

#include "gromacs/math/vec.h"

namespace gmx
{

/*! \internal
 * \brief Bounding box for a nbnxm atom cluster
 *
 * \note Should be aligned in memory to enable 4-wide SIMD operations.
 */
struct BoundingBox
{
    /*! \internal
     * \brief Corner for the bounding box, padded with one element to enable 4-wide SIMD operations
     */
    struct Corner
    {
        //! Returns a corner with the minimum coordinates along each dimension
        static Corner min(const Corner& c1, const Corner& c2)
        {
            Corner cMin;

            cMin.x = std::min(c1.x, c2.x);
            cMin.y = std::min(c1.y, c2.y);
            cMin.z = std::min(c1.z, c2.z);
            /* This value of the padding is irrelevant, as long as it
             * is initialized. We use min to allow auto-vectorization.
             */
            cMin.padding = std::min(c1.padding, c2.padding);

            return cMin;
        }

        //! Returns a corner with the maximum coordinates along each dimension
        static Corner max(const Corner& c1, const Corner& c2)
        {
            Corner cMax;

            cMax.x       = std::max(c1.x, c2.x);
            cMax.y       = std::max(c1.y, c2.y);
            cMax.z       = std::max(c1.z, c2.z);
            cMax.padding = std::max(c1.padding, c2.padding);

            return cMax;
        }

        //! Returns a pointer for SIMD loading of a Corner object
        const float* ptr() const { return &x; }

        //! Returns a pointer for SIMD storing of a Corner object
        float* ptr() { return &x; }

        //! x coordinate
        float x;
        //! y coordinate
        float y;
        //! z coordinate
        float z;
        //! padding, unused, but should be set to avoid operations on uninitialized data
        float padding;
    };

    //! lower, along x and y and z, corner
    Corner lower;
    //! upper, along x and y and z, corner
    Corner upper;
};

} // namespace gmx

#endif // GMX_NBNXM_BOUNDINGBOX_H
