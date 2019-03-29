/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
 * \brief
 * This file contains datatypes and function declarations necessary for mdrun
 * to interface with VMD via the Interactive Molecular Dynamics protocol.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_imd
 */

#ifndef GMX_IMD_IMD_H
#define GMX_IMD_IMD_H

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_domdec_t;
struct gmx_enerdata_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct t_commrec;
struct t_filenm;
struct t_gmx_IMD;
struct t_IMD;
struct t_inputrec;
class t_state;

namespace gmx
{
class MDLogger;
struct MdrunOptions;
}

static const char IMDstr[] = "IMD:";  /**< Tag output from the IMD module with this string. */

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
void write_IMDgroup_to_file(gmx_bool bIMD, t_inputrec *ir, const t_state *state,
                            const gmx_mtop_t *sys, int nfile, const t_filenm fnm[]);


/*! \brief Make a selection of the home atoms for the IMD group.
 *
 * Should be called at every domain decomposition.
 * Each node checks which of the atoms from "ind" are local and puts its local
 * atom numbers into the "ind_local" array. Furthermore, in "xa_ind" it is
 * stored at which position each local atom belongs in the assembled/collective
 * array, so that on the master node all positions can be merged into the
 * assembled array correctly.
 *
 * \param dd          Structure containing domain decomposition data.
 * \param imdSession  The IMD session
 */
void dd_make_local_IMD_atoms(const gmx_domdec_t *dd, t_gmx_IMD *imdSession);


/*! \brief Returns an initialized IMD session.
 *
 * This function is called before the main MD loop over time steps.
 *
 * \param ir           The inputrec structure containing the MD input parameters
 * \param cr           Information structure for MPI communication.
 * \param ms           Handler for multi-simulations.
 * \param top_global   The topology of the whole system.
 * \param mdlog        Logger
 * \param x            The starting positions of the atoms.
 * \param nfile        Number of files.
 * \param fnm          Struct containing file names etc.
 * \param oenv         Output options.
 * \param mdrunOptions Options for mdrun.
 *
 * \returns A pointer to an initialized IMD session.
 */
t_gmx_IMD *init_IMD(const t_inputrec *ir,
                    const t_commrec *cr,
                    const gmx_multisim_t *ms,
                    const gmx_mtop_t *top_global,
                    const gmx::MDLogger &mdlog,
                    const rvec x[],
                    int nfile, const t_filenm fnm[], const gmx_output_env_t *oenv,
                    const gmx::MdrunOptions &mdrunOptions);


/*! \brief IMD required in this time step?
 * Also checks for new IMD connection and syncs the nodes.
 *
 * \param IMDsetup     The IMD session.
 * \param step         The time step.
 * \param cr           Information structure for MPI communication.
 * \param bNS          Is this a neighbor searching step?
 * \param box          The simulation box.
 * \param x            The local atomic positions on this node.
 * \param t            The time.
 * \param wcycle       Count wallcycles of IMD routines for diagnostic output.
 *
 * \returns            Whether or not we have to do IMD communication at this step.
 */
gmx_bool do_IMD(t_gmx_IMD *IMDsetup, int64_t step, const t_commrec *cr,
                gmx_bool bNS,
                const matrix box, const rvec x[], double t,
                gmx_wallcycle *wcycle);


/*! \brief Add external forces from a running interactive molecular dynamics session.
 *
 * \param IMDsetup     The IMD session.
 * \param cr           Information structure for MPI communication.
 * \param f            The forces.
 * \param wcycle       Count wallcycles of IMD routines for diagnostic output.
 */
void IMD_apply_forces(t_gmx_IMD *IMDsetup,
                      const t_commrec *cr, rvec *f,
                      gmx_wallcycle *wcycle);


/*! \brief Copy energies and convert to float from enerdata to the IMD energy record.
 *
 * We do no conversion, so units in client are SI!
 *
 * \param IMDsetup         The IMD session.
 * \param enerd            Contains the GROMACS energies for the different interaction types.
 * \param step             The time step.
 * \param bHaveNewEnergies Only copy energies if we have done global summing of them before.
 *
 */
void IMD_fill_energy_record(t_gmx_IMD *IMDsetup, const gmx_enerdata_t *enerd,
                            int64_t step, gmx_bool bHaveNewEnergies);


/*! \brief Send positions and energies to the client.
 *
 * \param IMDsetup         The IMD session.
 */
void IMD_send_positions(t_gmx_IMD *IMDsetup);


/*! \brief Calls IMD_prepare_energies() and then IMD_send_positions().
 *
 * \param IMDsetup         The IMD session.
 * \param bIMDstep         If true, transfer the positions. Otherwise just update the time step and potentially the energy record.
 * \param enerd            Contains the GROMACS energies for the different interaction types.
 * \param step             The time step.
 * \param bHaveNewEnergies Only update the energy record if we have done global summing of the energies.
 * \param wcycle           Count wallcycles of IMD routines for diagnostic output.
 *
 */
void IMD_prep_energies_send_positions(t_gmx_IMD *IMDsetup, gmx_bool bIMDstep,
                                      const gmx_enerdata_t *enerd,
                                      int64_t step, gmx_bool bHaveNewEnergies,
                                      gmx_wallcycle *wcycle);

/*! \brief Finalize IMD and do some cleaning up.
 *
 * Currently, IMD finalize closes the force output file.
 *
 * \param IMDsetup     The IMD session.
 */
void IMD_finalize(t_gmx_IMD *IMDsetup);

#endif
