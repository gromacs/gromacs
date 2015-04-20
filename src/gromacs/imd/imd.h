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
 */

/*! \libinternal \file
 *
 * \brief
 * This file contains datatypes and function declarations necessary for mdrun
 * to interface with VMD via the interactive molecular dynamics protocol.
 *
 * \author Martin Hoefling, Carsten Kutzner <ckutzne@gwdg.de>
 *
 * \inlibraryapi
 * \ingroup module_imd
 */

#ifndef GMX_IMD_IMD_H
#define GMX_IMD_IMD_H

#include "config.h"

#include "gromacs/fileio/filenm.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/timing/wallcycle.h"

#ifdef GMX_NATIVE_WINDOWS
#include <Windows.h>
#define NOFLAGS 0
#endif


#ifdef __cplusplus
extern "C" {
#endif

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
extern void write_IMDgroup_to_file(gmx_bool bIMD, t_inputrec *ir, t_state *state,
                                   gmx_mtop_t *sys, int nfile, const t_filenm fnm[]);


/*! \brief Make a selection of the home atoms for the IMD group.
 *
 * Should be called at every domain decomposition.
 * Each node checks which of the atoms from "ind" are local and puts its local
 * atom numbers into the "ind_local" array. Furthermore, in "xa_ind" it is
 * stored at which position each local atom belongs in the assembled/collective
 * array, so that on the master node all positions can be merged into the
 * assembled array correctly.
 *
 * \param bIMD    Only springs into action if bIMD is TRUE. Otherwise returns directly.
 * \param dd      Structure containing domain decomposition data.
 * \param imd     The IMD group of atoms.
 */
extern void dd_make_local_IMD_atoms(gmx_bool bIMD, gmx_domdec_t *dd, t_IMD *imd);


/*! \brief Initializes (or disables) IMD.
 *
 * This function is called before the main MD loop over time steps,
 * and it must be called prior to any call to dd_partition_system if in parallel.
 *
 * \param ir           The inputrec structure containing the MD input parameters
 *                     including a pointer to the IMD data structure.
 * \param cr           Information structure for MPI communication.
 * \param top_global   The topology of the whole system.
 * \param fplog        General output file, normally md.log.
 * \param defnstimd    Default IMD update (=communication) frequency.
 * \param x            The starting positions of the atoms.
 * \param nfile        Number of files.
 * \param fnm          Struct containing file names etc.
 * \param oenv         Output options.
 * \param imdport      Port to use for IMD connections.
 * \param Flags        Flags passed over from main, used to determine
 *                     whether or not we are appending.
 */
extern void init_IMD(t_inputrec *ir, t_commrec *cr, gmx_mtop_t *top_global,
                     FILE *fplog, int defnstimd, rvec x[],
                     int nfile, const t_filenm fnm[], output_env_t oenv,
                     int imdport, unsigned long  Flags);


/*! \brief IMD required in this time step?
 * Also checks for new IMD connection and syncs the nodes.
 *
 * \param bIMD         Only springs into action if bIMD is TRUE. Otherwise returns directly.
 * \param step         The time step.
 * \param cr           Information structure for MPI communication.
 * \param bNS          Is this a neighbor searching step?
 * \param box          The simulation box.
 * \param x            The local atomic positions on this node.
 * \param ir           The inputrec structure containing the MD input parameters
 *                     including a pointer to the IMD data structure.
 * \param t            The time.
 * \param wcycle       Count wallcycles of IMD routines for diagnostic output.
 *
 * \returns            Whether or not we have to do IMD communication at this step.
 */
extern gmx_bool do_IMD(gmx_bool bIMD, gmx_int64_t step, t_commrec *cr, gmx_bool bNS,
                       matrix box, rvec x[], t_inputrec *ir, double t,
                       gmx_wallcycle_t wcycle);


/*! \brief Get the IMD update frequency.
 *
 * \param IMDsetup     Opaque pointer to IMD private data.
 *
 * \returns            The current IMD update/communication frequency
 */
extern int IMD_get_step(t_gmx_IMD IMDsetup);


/*! \brief Add external forces from a running interactive molecular dynamics session.
 *
 * \param bIMD         Returns directly if bIMD is FALSE.
 * \param imd          The IMD data structure.
 * \param cr           Information structure for MPI communication.
 * \param f            The forces.
 * \param wcycle       Count wallcycles of IMD routines for diagnostic output.
 */
extern void IMD_apply_forces(gmx_bool bIMD, t_IMD *imd, t_commrec *cr, rvec *f,
                             gmx_wallcycle_t wcycle);


/*! \brief Copy energies and convert to float from enerdata to the IMD energy record.
 *
 * We do no conversion, so units in client are SI!
 *
 * \param bIMD             Only springs into action if bIMD is TRUE. Otherwise returns directly.
 * \param imd              The IMD data structure.
 * \param enerd            Contains the GROMACS energies for the different interaction types.
 * \param step             The time step.
 * \param bHaveNewEnergies Only copy energies if we have done global summing of them before.
 *
 */
extern void IMD_fill_energy_record(gmx_bool bIMD, t_IMD *imd, gmx_enerdata_t *enerd,
                                   gmx_int64_t step, gmx_bool bHaveNewEnergies);


/*! \brief Send positions and energies to the client.
 *
 * \param imd              The IMD data structure.
 */
extern void IMD_send_positions(t_IMD *imd);


/*! \brief Calls IMD_prepare_energies() and then IMD_send_positions().
 *
 * \param bIMD             Returns directly if bIMD is FALSE.
 * \param bIMDstep         If true, transfer the positions. Otherwise just update the time step and potentially the energy record.
 * \param imd              The IMD data structure.
 * \param enerd            Contains the GROMACS energies for the different interaction types.
 * \param step             The time step.
 * \param bHaveNewEnergies Only update the energy record if we have done global summing of the energies.
 * \param wcycle           Count wallcycles of IMD routines for diagnostic output.
 *
 */
extern void IMD_prep_energies_send_positions(gmx_bool bIMD, gmx_bool bIMDstep,
                                             t_IMD *imd, gmx_enerdata_t *enerd,
                                             gmx_int64_t step, gmx_bool bHaveNewEnergies,
                                             gmx_wallcycle_t wcycle);

/*! \brief Finalize IMD and do some cleaning up.
 *
 * Currently, IMD finalize closes the force output file.
 *
 * \param bIMD         Returns directly if bIMD is FALSE.
 * \param imd          The IMD data structure.
 */
extern void IMD_finalize(gmx_bool bIMD, t_IMD *imd);


#ifdef __cplusplus
}
#endif

#endif
