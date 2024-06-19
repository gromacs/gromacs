/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * \brief
 * Declares functions to calculate both essential dynamics constraints
 * as well as flooding potentials and forces.
 *
 * \authors Bert de Groot <bgroot@gwdg.de>, Oliver Lange <oliver.lange@tum.de>,
 * Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 */
#ifndef GMX_ESSENTIALDYNAMICS_EDSAM_H
#define GMX_ESSENTIALDYNAMICS_EDSAM_H

#include <cstdint>

#include <memory>

#include "gromacs/math/vectypes.h"

/*! \brief Abstract type for essential dynamics
 *
 * The main type is defined only in edsam.cpp
 */
struct gmx_edsam;
struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct ObservablesHistory;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
class t_state;

namespace gmx
{
enum class StartingBehavior;
class Constraints;
template<typename>
class ArrayRef;

class EssentialDynamics
{
public:
    EssentialDynamics();
    ~EssentialDynamics();

    /*! \brief Getter for working data
     *
     * This is needed while the module is still under
     * construction. */
    gmx_edsam* getLegacyED();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};
class MDLogger;
} // namespace gmx

/*! \brief Applies essential dynamics constrains as defined in the .edi input file.
 *
 * \param ir                MD input parameter record.
 * \param step              Number of the time step.
 * \param cr                Data needed for MPI communication.
 * \param coords            The local positions on this processor.
 * \param velocities        The local velocities.
 * \param box               The simulation box.
 * \param ed                The essential dynamics data.
 */
void do_edsam(const t_inputrec*        ir,
              int64_t                  step,
              const t_commrec*         cr,
              gmx::ArrayRef<gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec> velocities,
              const matrix             box,
              gmx_edsam*               ed);


/*! \brief Initializes the essential dynamics and flooding module.
 *
 * \param mdlog             Logger.
 * \param ediFileName       Essential dynamics input file.
 * \param edoFileName       Output file for essential dynamics data.
 * \param mtop              Molecular topology.
 * \param ir                MD input parameter record.
 * \param cr                Data needed for MPI communication.
 * \param constr            Data structure keeping the constraint information.
 * \param globalState       The global state, only used on the main rank.
 * \param oh                The observables history container.
 * \param oenv              The output environment information.
 * \param startingBehavior  Describes whether this is a restart appending to output files
 *
 * \returns                 A pointer to the ED data structure.
 */
std::unique_ptr<gmx::EssentialDynamics> init_edsam(const gmx::MDLogger&    mdlog,
                                                   const char*             ediFileName,
                                                   const char*             edoFileName,
                                                   const gmx_mtop_t&       mtop,
                                                   const t_inputrec&       ir,
                                                   const t_commrec*        cr,
                                                   gmx::Constraints*       constr,
                                                   const t_state*          globalState,
                                                   ObservablesHistory*     oh,
                                                   const gmx_output_env_t* oenv,
                                                   gmx::StartingBehavior   startingBehavior);

/*! \brief Make a selection of the home atoms for the ED groups.
 *
 * Should be called at every domain decomposition.
 *
 * \param dd                Domain decomposition data.
 * \param ed                Essential dynamics and flooding data.
 */
void dd_make_local_ed_indices(gmx_domdec_t* dd, gmx_edsam* ed);


/*! \brief Evaluate the flooding potential(s) and forces as requested in the .edi input file.
 *
 * \param cr                Data needed for MPI communication.
 * \param ir                MD input parameter record.
 * \param coords            Positions on the local processor.
 * \param force             Forcefield forces to which the flooding forces are added.
 * \param ed                The essential dynamics data.
 * \param box               The simulation box.
 * \param step              Number of the time step.
 * \param bNS               Are we in a neighbor searching step?
 */
void do_flood(const t_commrec*               cr,
              const t_inputrec&              ir,
              gmx::ArrayRef<const gmx::RVec> coords,
              gmx::ArrayRef<gmx::RVec>       force,
              gmx_edsam*                     ed,
              const matrix                   box,
              int64_t                        step,
              bool                           bNS);

#endif
