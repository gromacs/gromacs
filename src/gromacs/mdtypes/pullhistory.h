/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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

//! \brief Contains the sum of coordinate observables to enable calculation of the average of pull data.
class PullCoordinateHistory
{
public:
    double value;       //!< The sum of the current value of the coordinate, units of nm or rad.
    double valueRef;    //!< The sum of the reference value of the coordinate, units of nm or rad.
    double scalarForce; //!< The sum of the scalar force of the coordinate.
    dvec   dr01;        //!< The sum of the direction vector of group 1 relative to group 0.
    dvec   dr23;        //!< The sum of the direction vector of group 3 relative to group 2.
    dvec   dr45;        //!< The sum of the direction vector of group 5 relative to group 4.
    dvec   dynaX;       //!< The sum of the coordinate of the dynamic groups for geom=cylinder.

    //! Constructor
    PullCoordinateHistory() : value(0), valueRef(0), scalarForce(0), dr01(), dr23(), dr45(), dynaX()
    {
    }
};

//! \brief Contains the sum of group observables to enable calculation of the average of pull data.
class PullGroupHistory
{
public:
    dvec x; //!< The sum of the coordinates of the group.

    //! Constructor
    PullGroupHistory() : x() {}
};


//! \brief Pull statistics history, to allow output of average pull data.
class PullHistory
{
public:
    int numValuesInXSum; //!< Number of values of the coordinate values in the pull sums.
    int numValuesInFSum; //!< Number of values in the pull force sums.
    std::vector<PullCoordinateHistory> pullCoordinateSums; //!< The container of the sums of the values of the pull coordinate, also contains the scalar force.
    std::vector<PullGroupHistory> pullGroupSums; //!< The container of the sums of the values of the pull group.

    //! Constructor
    PullHistory() : numValuesInXSum(0), numValuesInFSum(0) {}
};

//! \endcond

#endif
