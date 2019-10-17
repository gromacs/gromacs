/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2018,2019, by the GROMACS development team, led by
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
 * This file contains datatypes for energy statistics history.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDLIB_ENERGYHISTORY_H
#define GMX_MDLIB_ENERGYHISTORY_H

#include <memory>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

//! \cond INTERNAL

//! \brief Energy history for delta_h histograms in between energy file frames
class delta_h_history_t
{
public:
    //! Vector (size number of intermediate data points) of vector of Hamiltonian differences for each foreign lambda
    std::vector<std::vector<real>> dh;
    //! The start time of these energy diff blocks
    double start_time;
    //! Lambda at start time
    double start_lambda;
    //! Whether the lambda value is set. Here for backward-compatibility.
    gmx_bool start_lambda_set;

    delta_h_history_t() : start_time(0), start_lambda(0), start_lambda_set(false) {}
};


//! \brief Energy statistics history, only used for output and reporting
class energyhistory_t
{
public:
    int64_t             nsteps;       //! The number of steps in the history
    int64_t             nsum;         //! Nr. of steps in the ener_ave and ener_sum
    std::vector<double> ener_ave;     //! Energy terms difference^2 sum to get fluctuations
    std::vector<double> ener_sum;     //! Energy terms sum
    int64_t             nsteps_sim;   //! The number of steps in ener_sum_sim
    int64_t             nsum_sim;     //! The number of frames in ener_sum_sim
    std::vector<double> ener_sum_sim; //! Energy term history sum of the whole sim

    //! History for energy difference for foreign lambdas (useful for BAR)
    std::unique_ptr<delta_h_history_t> deltaHForeignLambdas;

    energyhistory_t() :
        nsteps(0),
        nsum(0),

        nsteps_sim(0),
        nsum_sim(0),
        ener_sum_sim(0)
    {
    }
};

//! \endcond

#endif
