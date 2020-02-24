/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * \brief Declares the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_MULTISIM_H
#define GMX_MDRUNUTILITY_MULTISIM_H

#include <memory>
#include <string>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/mpiinplacebuffers.h"

/*! \libinternal
 * \brief Coordinate multi-simulation resources for mdrun
 *
 * \todo Change this to class
 */
struct gmx_multisim_t
{
    //! Default constructor
    gmx_multisim_t();
    /*! \brief Constructor useful for mdrun simulations
     *
     * Splits the communicator into multidirs.size() separate
     * simulations, if >1, and creates a communication structure
     * between the master these simulations.
     *
     * Valid to call regardless of build configuration, but \c
     * multidirs must be empty unless a real MPI build is used. */
    gmx_multisim_t(MPI_Comm comm, gmx::ArrayRef<const std::string> multidirs);
    //! Destructor
    ~gmx_multisim_t();

    //! The number of simulations in the set of multi-simulations
    int nsim = 1;
    //! The index of the simulation that owns this object within the set
    int sim = 0;
    //! The MPI Group between master ranks of simulations, valid only on master ranks.
    MPI_Group mpi_group_masters = MPI_GROUP_NULL;
    //! The MPI communicator between master ranks of simulations, valid only on master ranks.
    MPI_Comm mpi_comm_masters = MPI_COMM_NULL;
    //! Communication buffers needed if MPI_IN_PLACE isn't supported
    mpi_in_place_buf_t* mpb = nullptr;
};

//! Calculate the sum over the simulations of an array of ints
void gmx_sumi_sim(int nr, int r[], const gmx_multisim_t* ms);

//! Calculate the sum over the simulations of an array of large ints
void gmx_sumli_sim(int nr, int64_t r[], const gmx_multisim_t* ms);

//! Calculate the sum over the simulations of an array of floats
void gmx_sumf_sim(int nr, float r[], const gmx_multisim_t* ms);

//! Calculate the sum over the simulations of an array of doubles
void gmx_sumd_sim(int nr, double r[], const gmx_multisim_t* ms);

/*! \brief Return a vector containing the gathered values of \c
 * localValue found on the master rank of each simulation. */
std::vector<int> gatherIntFromMultiSimulation(const gmx_multisim_t* ms, int localValue);

/*! \brief Check if val is the same on all simulations for a mdrun
 * -multidir run
 *
 * The string name is used to print to the log file and in a fatal error
 * if the val's don't match. If bQuiet is true and the check passes,
 * no output is written. */
void check_multi_int(FILE* log, const gmx_multisim_t* ms, int val, const char* name, gmx_bool bQuiet);
/*! \copydoc check_multi_int() */
void check_multi_int64(FILE* log, const gmx_multisim_t* ms, int64_t val, const char* name, gmx_bool bQuiet);

#if GMX_DOUBLE
//! Convenience define for sum of reals
#    define gmx_sum_sim gmx_sumd_sim
#else
//! Convenience define for sum of reals
#    define gmx_sum_sim gmx_sumf_sim
#endif

//! Are we doing multiple independent simulations?
static bool inline isMultiSim(const gmx_multisim_t* ms)
{
    return ms != nullptr;
}

/*! \brief Return whether this rank is the master rank of a
 * simulation, using \c ms (if it is valid) and otherwise \c
 * communicator */
bool findIsSimulationMasterRank(const gmx_multisim_t* ms, MPI_Comm communicator);

//! Are we the master simulation of a possible multi-simulation?
bool isMasterSim(const gmx_multisim_t* ms);

/*! \brief Are we the master rank (of the master simulation, for a multi-sim).
 *
 * This rank prints the remaining run time etc. */
bool isMasterSimMasterRank(const gmx_multisim_t* ms, bool isMaster);

#endif
