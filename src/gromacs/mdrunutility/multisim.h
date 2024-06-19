/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares the multi-simulation support routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_MULTISIM_H
#define GMX_MDRUNUTILITY_MULTISIM_H

#include <cstdint>
#include <cstdio>

#include <memory>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"

namespace gmx
{
class MDLogger;
}

struct gmx_multisim_t;
struct t_commrec;

/*! \libinternal
 * \brief Builder function for gmx_multisim_t
 *
 * \param[in]  worldComm   MPI communicator to split when
 *                         multi-simulation is requested.
 * \param[in]  multidirs   Strings naming the subdirectories when
 *                         multi-simulation is requested, otherwise empty
 *
 * Splits \c worldComm into \c multidirs.size() separate
 * simulations, if >1, and creates a communication structure
 * between the main ranks of these simulations.
 *
 * Valid to call regardless of build configuration, but \c
 * multidirs must be empty unless a real MPI build is used.
 *
 * \throws NotImplementedError     when \c multidirs is non-empty unless using real MPI is true
 * \throws NotImplementedError     when \c multidirs has exactly one element
 * \throws InconsistentInputError  when the number of MPI ranks is not a multiple of the number of \c multidirs
 * \throws FileIOError             when the simulation cannot change to the working directory in \c multidirs
 */
std::unique_ptr<gmx_multisim_t> buildMultiSimulation(MPI_Comm                         worldComm,
                                                     gmx::ArrayRef<const std::string> multidirs);

/*! \libinternal
 * \brief Coordinate multi-simulation resources for mdrun
 *
 * \todo Change this to class
 */
struct gmx_multisim_t
{
    /*! \brief Constructor
     *
     * \param[in]  numSimulations   The number of simulations in the MPI world.
     * \param[in]  simulationIndex  The index of this simulation in the set of simulations.
     * \param[in]  mainsComm        On main ranks, the communicator among main ranks;
     *                              otherwise MPI_COMM_NULL.
     * \param[in]  simulationComm   The communicator among ranks of this simulation.
     *
     * Assumes ownership of the communicators if they are neither
     * MPI_COMM_WORLD nor MPI_COMM_NULL. If so, upon destruction will
     * call MPI_Comm_free on them.
     */
    gmx_multisim_t(int numSimulations, int simulationIndex, MPI_Comm mainsComm, MPI_Comm simulationComm);
    //! Destructor
    ~gmx_multisim_t();

    //! The number of simulations in the set of multi-simulations
    int numSimulations_ = 1;
    //! The index of the simulation that owns this object within the set
    int simulationIndex_ = 0;
    //! The MPI communicator between main ranks of simulations, valid only on main ranks.
    MPI_Comm mainRanksComm_ = MPI_COMM_NULL;
    //! The MPI communicator between ranks of this simulation.
    MPI_Comm simulationComm_ = MPI_COMM_NULL;
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
 * localValue found on the main rank of each simulation. */
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

/*! \brief Return whether this rank is the main rank of a
 * simulation, using \c ms (if it is valid) and otherwise \c
 * communicator */
bool findIsSimulationMainRank(const gmx_multisim_t* ms, MPI_Comm communicator);

//! Are we the main simulation of a possible multi-simulation?
bool isMainSim(const gmx_multisim_t* ms);

/*! \brief Are we the main rank (of the main simulation, for a multi-sim).
 *
 * This rank prints the remaining run time etc. */
bool isMainSimMainRank(const gmx_multisim_t* ms, bool isMain);

/*! \brief Log the initial state of the multi-sim
 *
 * The simulations may be at different steps, etc so we
 * report that.
 *
 * \param[in]  ms                     The multi-sum object
 * \param[in]  cr                     The commrec object
 * \param[in]  mdlog                  Logger
 * \param[in]  simulationsShareState  Whether the simulations share state
 * \param[in]  numSteps               The number of steps in this simulation
 * \param[in]  initialStep            The initial step for this simulation
 */
void logInitialMultisimStatus(const gmx_multisim_t* ms,
                              const t_commrec*      cr,
                              const gmx::MDLogger&  mdlog,
                              bool                  simulationsShareState,
                              int                   numSteps,
                              int                   initialStep);

#endif
