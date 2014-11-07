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
/*! \internal \file
 * \brief
 * Declares helper class for autocorrelation tests
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \ingroup module_correlationfunctions
 */
#ifndef GMX_CORRELATIONDATASET_H
#define GMX_CORRELATIONDATASET_H

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

class CorrelationDataSet
{
    double         ** tempValues_;

    int               nrLines_;
    int               nrColumns_;
    double            startTime_;
    double            endTime_;
    double            dt_;

    public:

        /*! \brief
         * Constructor
         * \param[in] fileName containing function to test. *.xvg
         */
        explicit CorrelationDataSet(const std::string &fileName);

        /*! \brief
         * Return a value at an index
         * \param[in] set the set number
         * \param[in] t the time index of the value
         */
        real getValue(int set, int t) const;

        /*! \brief
         * Return the nummber of columns
         */
        int getNrColumns() const { return nrColumns_; }

        /*! \brief
         * Return the nummber of Lines
         */
        int getNrLines() const { return nrLines_; }

        /*! \brief
         * Return the time witch the function starts at
         */
        real getStartTime() const { return startTime_; }

        /*! \brief
         * Return the time the function ends at
         */
        real getEndTime() const { return endTime_; }

        /*! \brief
         * return delta time
         */
        real getDt() const { return dt_; }

        /*! \brief
         * Destructor
         */
        ~CorrelationDataSet();

    private:
        //! This class should not be copyable or assignable
        GMX_DISALLOW_COPY_AND_ASSIGN(CorrelationDataSet);
};

#endif
