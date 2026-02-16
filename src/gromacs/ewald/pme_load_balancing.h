/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
 * \brief This file declares the PmeLoadBalancing class
 *
 * This class managing automatic load balance of PME calculations (Coulomb and
 * LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_ewald
 */

#ifndef GMX_EWALD_PME_LOAD_BALANCING_H
#define GMX_EWALD_PME_LOAD_BALANCING_H

#include <memory>

#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/vectypes.h"


enum class CoulombInteractionType;
struct gmx_domdec_t;
struct t_forcerec;
struct t_inputrec;
struct interaction_const_t;
struct gmx_pme_t;
class t_state;

namespace gmx
{
struct nonbonded_verlet_t;
class MDLogger;
template<typename T>
class ArrayRef;
class SimulationWorkload;

//! Returns whether PME tuning is supported
bool pmeTuningIsSupported(CoulombInteractionType    coulombInteractionType,
                          bool                      reproducibilityRequested,
                          const SimulationWorkload& simulationWork);

/*! \brief Object to manage PME load balancing */
class PmeLoadBalancing
{
public:
    /*! \brief Initialize the PP-PME load balancing data and infrastructure
     *
     * Initialize the PP-PME load balancing data and infrastructure.
     * The actual load balancing might start right away, later or never.
     * The PME grid in \p pmedata is reused for smaller grids to lower the memory
     * usage. The memory passed in \p pmedata needs to be freed after destructing
     * PmeLoadBalancing object.
     *
     * \note This constructor should only be called when \c pmeTuningIsSupported() returns true.
     */
    PmeLoadBalancing(gmx_domdec_t*              dd,
                     const MDLogger&            mdlog,
                     const t_inputrec&          ir,
                     const matrix               box,
                     const interaction_const_t& ic,
                     const nonbonded_verlet_t&  nbv,
                     gmx_pme_t*                 pmedata,
                     const SimulationWorkload&  simulationWork);

    ~PmeLoadBalancing();

    PmeLoadBalancing(const PmeLoadBalancing&)            = delete;
    PmeLoadBalancing& operator=(const PmeLoadBalancing&) = delete;
    PmeLoadBalancing(PmeLoadBalancing&&)                 = delete;
    PmeLoadBalancing& operator=(PmeLoadBalancing&&)      = delete;

    /*! \brief Return whether PME load balancing is active */
    bool isActive() const;

    /*! \brief Return whether PME load balancing is printing load numbers to the log file */
    bool isPrintingLoad() const;

    /*! \brief Process cycles and PME load balance when necessary
     *
     * Process the cycles measured over the last nstlist steps and then
     * either continue balancing or check if we need to trigger balancing.
     * Should be called after the WallCycleCounter::Step cycle counter has been stopped.
     */
    void addCycles(FILE*                fp_err,
                   t_forcerec*          fr,
                   const matrix         box,
                   ArrayRef<const RVec> x,
                   const gmx_wallcycle* wcycle,
                   int64_t              step,
                   int64_t              step_rel);

    /*! \brief Print the chosen PME settings to the mdlogger */
    void printSettings() const;

private:
    //! Implementation type
    class Impl;
    //! Implementation object
    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
