/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_PLOTFILE_H
#define GMX_FILEIO_PLOTFILE_H

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;

namespace gmx
{

//! Forward declaration for PlotFileImpl
class PlotFileImpl;

/*! \brief Wrapper around xvgr routines
 *
 * This wrapper automates the closing of xvg files in a robust manner
 * when the variable goes out of scope.
 *
 * \todo Extend to a full class by implementing other functions
 * \todo Implement in analysis tools
 * \todo Implement error handling in case the file cannot be written
 */
class PlotFile
{
    public:
        /*! \brief Constructor
         *
         * \param[in] fn    The outpuf file name
         * \param[in] title The title of the graph
         * \param[in] xaxis The X-axis legend
         * \param[in] yaxis The Y-axis legend
         * \param[in] oenv  The Gromacs output environment
         */
        PlotFile(const char                    *fn,
                 const char                    *title,
                 const char                    *xaxis,
                 const char                    *yaxis,
                 const struct gmx_output_env_t *oenv);

        //! Disallow copy constructor
        PlotFile(const PlotFile &) = delete;

        ~PlotFile();

        /*! \brief Add a legend description for a set
         * The order for calling this function should correspond
         * to the data.
         * \param[in] legend The string describing the data
         */
        void addLegend(const std::string &legend)
        {
            legend_.push_back(legend);
        }

        //! Add a comment
        void addComment(const std::string &comment)
        {
            comment_.push_back(comment);
        }

        /*! \brief Add "Time" series
         *
         * This resets any previous stored time series
         * \param[in] t The time
         */
        void setTimeAxis(const std::vector<real> &t);

        /*! \brief Set the data array
         *
         * Store one or more y data series
         * Existing data arrays willbe removed
         * \param[in] data The 2D data array
         */
        void setDataSets(const std::vector<std::vector<real> > &data);

        /*! Plot the file at last
         *
         * \todo Implement flag to select format instead of xvg
         */
        void plot();

    private:
        //! Low level handling of files
        PlotFileImpl                   *impl_;

        //! Comments
        std::vector<std::string>        comment_;

        //! Legends
        std::vector<std::string>        legend_;

        //! "Time" series
        std::vector<real>               t_;

        //! Data series
        std::vector<std::vector<real> > data_;

};

}

#endif
