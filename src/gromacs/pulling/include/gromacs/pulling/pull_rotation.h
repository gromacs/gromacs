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
 * Declares functions to enforce rotational motion upon a group of particles.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 */

#ifndef GMX_PULLING_PULL_ROTATION_H
#define GMX_PULLING_PULL_ROTATION_H

#include <cstdio>

#include <memory>

#include "gromacs/math/vectypes.h"

struct gmx_domdec_t;
struct gmx_enfrot;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct t_commrec;
struct t_filenm;
struct t_inputrec;
struct t_rot;
class t_state;

namespace gmx
{
enum class StartingBehavior;
class LocalAtomSetManager;
struct MdrunOptions;
template<typename>
class ArrayRef;

class EnforcedRotation
{
public:
    EnforcedRotation();
    ~EnforcedRotation();

    /*! \brief Getter for working data
     *
     * This is needed while the module is still under
     * construction. */
    gmx_enfrot* getLegacyEnfrot();

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

/*! \brief Initializes the enforced rotation groups.
 *
 * This routine does the memory allocation for various helper arrays, opens
 * the output files etc.
 *
 * \param fplog    General output file, normally md.log.
 * \param ir       Struct containing MD input parameters, among those
 *                 also the enforced rotation parameters.
 * \param nfile    Number of entries in the fnm structure.
 * \param fnm      The filenames struct containing also the names
 *                 of the rotation output files.
 * \param atomSets Tracks indices of atoms subject to enforced rotation for each DD rank.
 * \param cr       Pointer to MPI communication data.
 * \param globalState  The global state, only used on the master rank.
 * \param mtop     Molecular topology.
 * \param oenv     Needed to open the rotation output xvgr file.
 * \param mdrunOptions  Options for mdrun.
 * \param startingBehavior  Describes whether this is a restart appending to output files
 * \return         An enforced rotation module.
 */
std::unique_ptr<gmx::EnforcedRotation> init_rot(FILE*                     fplog,
                                                t_inputrec*               ir,
                                                int                       nfile,
                                                const t_filenm            fnm[],
                                                const t_commrec*          cr,
                                                gmx::LocalAtomSetManager* atomSets,
                                                const t_state*            globalState,
                                                const gmx_mtop_t&         mtop,
                                                const gmx_output_env_t*   oenv,
                                                const gmx::MdrunOptions&  mdrunOptions,
                                                gmx::StartingBehavior     startingBehavior);

/*! \brief Calculates the enforced rotation potential(s).
 *
 * This is the main enforced rotation module which is called during every time
 * step. Here the rotation potential as well as the resulting forces are
 * calculated.
 *
 * \param cr      Pointer to MPI communication data.
 * \param er      Pointer to the enforced rotation working data.
 * \param box     Simulation box, needed to make group whole.
 * \param coords  The positions of all the local particles.
 * \param t       Time.
 * \param step    The time step.
 * \param bNS     After domain decomposition / neighbor searching several
 *                local arrays have to be updated (masses, shifts)
 */
void do_rotation(const t_commrec*               cr,
                 gmx_enfrot*                    er,
                 const matrix                   box,
                 gmx::ArrayRef<const gmx::RVec> coords,
                 real                           t,
                 int64_t                        step,
                 bool                           bNS);


/*! \brief Add the enforced rotation forces to the official force array.
 *
 * Adds the forces from enforced rotation potential to the local forces and
 * sums up the contributions to the rotation potential from all the nodes. Since
 * this needs communication, this routine should be called after the short range
 * forces have been evaluated (in order not to spoil cycle counts).
 * This routine also outputs data to the rotation output files (e.g.
 * the potential, the angle of the group(s), and torques).
 *
 * \param er      Pointer to the enforced rotation working data.
 * \param force   The local forces to which the rotational forces have
 *                to be added.
 * \param cr      Pointer to MPI communication data.
 * \param step    The time step, used for output.
 * \param t       Time, used for output.
 * \returns       The potential energy of the rotation potentials.
 */
real add_rot_forces(gmx_enfrot* er, gmx::ArrayRef<gmx::RVec> force, const t_commrec* cr, int64_t step, real t);


#endif
