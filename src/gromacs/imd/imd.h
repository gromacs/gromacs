/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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

/*! \libinternal
 * \defgroup module_imd Interactive molecular dynamics (IMD)
 * \ingroup group_mdrun
 *
 * \brief
 * Allows mdrun to interface with VMD via the interactive molecular dynamics
 * (IMD) protocol.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \todo Rename the directory, source and test files to
 * interactive_md, and prefer type names like
 * InteractiveMDSession. Avoid ambiguity with IMDModule.
 */

/*! \libinternal \file
 *
 * \brief This file declares the class that coordinates with VMD via
 * the Interactive Molecular Dynamics protocol, along with supporting
 * free functions.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_imd
 */

#ifndef GMX_IMD_IMD_H
#define GMX_IMD_IMD_H

#include <cstdint>

#include <memory>

#include "gromacs/math/vectypes.h"

struct gmx_domdec_t;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_filenm;
struct t_IMD;
struct t_inputrec;
class t_state;

namespace gmx
{
enum class StartingBehavior;
class IMDModule;
class ImdSession;
class InteractiveMolecularDynamics;
class MDLogger;
struct ImdOptions;
struct MdrunOptions;
template<typename>
class ArrayRef;

/*! \brief
 * Creates a module for interactive molecular dynamics.
 */
std::unique_ptr<IMDModule> createInteractiveMolecularDynamicsModule();

static const char IMDstr[] = "IMD:"; /**< Tag output from the IMD module with this string. */

/*! \brief Writes out the group of atoms selected for interactive manipulation.
 *
 * Called by grompp.
 * The resulting file has to be read in by VMD if one wants it to connect to mdrun.
 *
 * \param bIMD    Only springs into action if bIMD is TRUE. Otherwise returns directly.
 * \param ir      Structure containing MD input parameters, among those
 *                the IMD data structure.
 * \param state   The current state of the MD system.
 * \param sys     The global, complete system topology.
 * \param nfile   Number of files.
 * \param fnm     Filename struct.
 */
void write_IMDgroup_to_file(bool              bIMD,
                            t_inputrec*       ir,
                            const t_state*    state,
                            const gmx_mtop_t& sys,
                            int               nfile,
                            const t_filenm    fnm[]);


/*! \brief Makes and returns an initialized IMD session, which may be inactive.
 *
 * This function is called before the main MD loop over time steps.
 *
 * \param ir           The inputrec structure containing the MD input parameters
 * \param cr           Information structure for MPI communication.
 * \param wcycle       Count wallcycles of IMD routines for diagnostic output.
 * \param enerd        Contains the GROMACS energies for the different interaction types.
 * \param ms           Handler for multi-simulations.
 * \param top_global   The topology of the whole system.
 * \param mdlog        Logger
 * \param coords       The starting positions of the atoms.
 * \param nfile        Number of files.
 * \param fnm          Struct containing file names etc.
 * \param oenv         Output options.
 * \param options      Options for interactive MD.
 * \param startingBehavior  Describes whether this is a restart appending to output files
 */
std::unique_ptr<ImdSession> makeImdSession(const t_inputrec*              ir,
                                           const t_commrec*               cr,
                                           gmx_wallcycle*                 wcycle,
                                           gmx_enerdata_t*                enerd,
                                           const gmx_multisim_t*          ms,
                                           const gmx_mtop_t&              top_global,
                                           const MDLogger&                mdlog,
                                           gmx::ArrayRef<const gmx::RVec> coords,
                                           int                            nfile,
                                           const t_filenm                 fnm[],
                                           const gmx_output_env_t*        oenv,
                                           const ImdOptions&              options,
                                           StartingBehavior               startingBehavior);

class ImdSession
{
private:
    //! Private constructor, to force the use of makeImdSession()
    ImdSession(const MDLogger& mdlog);

public:
    ~ImdSession();

    /*! \brief Make a selection of the home atoms for the IMD group.
     *
     * Should be called at every domain decomposition.
     * Each node checks which of the atoms from "ind" are local and puts its local
     * atom numbers into the "ind_local" array. Furthermore, in "xa_ind" it is
     * stored at which position each local atom belongs in the assembled/collective
     * array, so that on the main node all positions can be merged into the
     * assembled array correctly.
     *
     * \param dd          Structure containing domain decomposition data.
     */
    void dd_make_local_IMD_atoms(const gmx_domdec_t* dd);

    /*! \brief Prepare to do IMD if required in this time step
     * Also checks for new IMD connection and syncs the nodes.
     *
     * \param step         The time step.
     * \param bNS          Is this a neighbor searching step?
     * \param box          The simulation box.
     * \param coords       The local atomic positions on this node.
     * \param t            The time.
     *
     * \returns            Whether or not we have to do IMD communication at this step.
     */
    bool run(int64_t step, bool bNS, const matrix box, gmx::ArrayRef<const gmx::RVec> coords, double t);

    /*! \brief Add external forces from a running interactive molecular dynamics session.
     *
     * \param force The forces.
     */
    void applyForces(gmx::ArrayRef<gmx::RVec> force);

    /*! \brief Copy energies and convert to float from enerdata to the IMD energy record.
     *
     * We do no conversion, so units in client are SI!
     *
     * \param step             The time step.
     * \param bHaveNewEnergies Only copy energies if we have done global summing of them before.
     */
    void fillEnergyRecord(int64_t step, bool bHaveNewEnergies);

    /*! \brief Send positions and energies to the client. */
    void sendPositionsAndEnergies();

    /*! \brief Updates the energy record sent to VMD if needed, and sends positions and energies.
     *
     * \param bIMDstep         If true, transfer the positions. Otherwise just update the time step and potentially the energy record.
     * \param step             The time step.
     * \param bHaveNewEnergies Update the energy record if we have done global summing of the energies.
     */
    void updateEnergyRecordAndSendPositionsAndEnergies(bool bIMDstep, int64_t step, bool bHaveNewEnergies);

private:
    //! Implementation type.
    class Impl;
    //! Implementation object.
    std::unique_ptr<Impl> impl_;

public:
    // Befriend the factory function.
    friend std::unique_ptr<ImdSession> makeImdSession(const t_inputrec*              ir,
                                                      const t_commrec*               cr,
                                                      gmx_wallcycle*                 wcycle,
                                                      gmx_enerdata_t*                enerd,
                                                      const gmx_multisim_t*          ms,
                                                      const gmx_mtop_t&              top_global,
                                                      const MDLogger&                mdlog,
                                                      gmx::ArrayRef<const gmx::RVec> coords,
                                                      int                            nfile,
                                                      const t_filenm                 fnm[],
                                                      const gmx_output_env_t*        oenv,
                                                      const ImdOptions&              options,
                                                      StartingBehavior startingBehavior);
};

} // namespace gmx

#endif
