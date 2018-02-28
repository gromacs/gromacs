/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017,2018, by the GROMACS development team, led by
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
 *
 *
 * \brief
 * This file contains datatypes for pull statistics history.
 *
 * \author Magnus Lundborg, Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDLIB_PULLHISTORY_H
#define GMX_MDLIB_PULLHISTORY_H

#include <vector>

//! \cond INTERNAL

//! \brief Pull statistics history, to allow output of average pull data.
class PullHistory
{
    public:
        int                 numCoordinates;         //!< The number of pull coordinates.
        int                 numValuesPerCoordinate; //!< The number of values per pull coordinate (depends on what is to be output).
        int                 numValuesInSum;         //!< Number of steps in the ener_ave and ener_sum.
        std::vector<double> sum;                    //!< Sum of pull force or coordinates (n=numCoordinates*numValuesPerCoordinate).

        //! Constructor
        PullHistory() : numCoordinates(0),
                        numValuesPerCoordinate(0),
                        numValuesInSum(0),
                        sum()
        {
        }
};

//! \endcond

#endif
