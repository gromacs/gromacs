/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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



#include <vector>

#ifndef GMX_CORRELATIONDATASET_H
#define GMX_CORRELATIONDATASET_H

#ifdef __cplusplus
extern "C" {
#endif


class CorrelationDataSet
{
    std::vector<real> values;

    int               nrLines;
    real              startTime;
    real              endTime;
    real              dt;
    int               nrColums;

    public:

/*! \brief
 * Constructor
 * \param[in] file containing fucntion to test. *.xvg
 * \param[in] the dimeson of the .xvg file, the nummber of colums to read
 */
        CorrelationDataSet(std::string fileName, int dim);

/*! \brief
 * Return a value at a index
 * \param[in] x the index of the value
 */
        real getValue(int x);


/*! \brief
 * Return the nummber of colums
 */
        int getNrColums();

/*! \brief
 * Return the nummber of Lines
 */
        int getNrLines();

/*! \brief
 * Return the time witch the function starts at
 */
        real getStartTime();

/*! \brief
 * Return the time the function ends at
 */
        real getEndTime();

/*! \brief
 * return delta time
 */
        real getDt();

/*! \brief
 * Destructor
 */
        ~CorrelationDataSet()
        {

        }


};
#ifdef __cplusplus
}
#endif

#endif
