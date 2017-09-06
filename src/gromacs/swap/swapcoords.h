/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * \defgroup module_swap "Computational Electrophysiology" position swapping (swap)
 * \ingroup group_mdrun
 * \brief
 * Implements the "Computational Electrophysiology" protocol.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 */
/*! \libinternal \file
 * \brief
 * The "Computational Electrophysiology" protocol for ion/water position swapping.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \inlibraryapi
 * \ingroup module_swap
 */
#ifndef GMX_SWAP_SWAPCOORDS_H
#define GMX_SWAP_SWAPCOORDS_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_domdec_t;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct gmx_wallcycle;
struct MdrunOptions;
struct swaphistory_t;
struct t_commrec;
struct t_inputrec;
class t_state;
struct t_swapcoords;
struct ObservablesHistory;


/*! \brief Initialize ion / water position swapping ("Computational Electrophysiology").
 *
 * This routine does the memory allocation for various helper arrays, opens
 * the output file, sets up swap data checkpoint writing, etc.
 *
 * \param[in] fplog         General output file, normally md.log.
 * \param[in] ir            Structure containing MD input parameters, among those
 *                          also the structure needed for position swapping.
 * \param[in] fn            Output file name for swap data.
 * \param[in] mtop          Molecular topology.
 * \param[in] globalState   The global state, only used on the master rank.
 * \param[in] oh            Contains struct with swap data that is read from or written to checkpoint.
 * \param[in] cr            Pointer to MPI communication data.
 * \param[in] oenv          Needed to open the swap output XVGR file.
 * \param[in] mdrunOptions  Options for mdrun.
 */
void init_swapcoords(
        FILE                   *fplog,
        t_inputrec             *ir,
        const char             *fn,
        gmx_mtop_t             *mtop,
        const t_state          *globalState,
        ObservablesHistory     *oh,
        t_commrec              *cr,
        const gmx_output_env_t *oenv,
        const MdrunOptions     &mdrunOptions);


/*! \brief Finalizes ion / water position swapping.
 *
 * \param[in] sc            Pointer to swap data.
 */
void finish_swapcoords(t_swapcoords *sc);


/*! \brief Make a selection of the home atoms for the swap groups. These are
 * the ions, the water, and the channels. This routine should be called at every
 * domain decomposition.
 *
 * \param[in] dd            Structure containing domain decomposition data.
 * \param[in] si_pub        Pointer to the swap data structure.
 */
void dd_make_local_swap_groups(gmx_domdec_t *dd, t_swapcoords *si_pub);


/*! \brief "Computational Electrophysiology" main routine within MD loop.
 *
 * \param[in] cr       Pointer to MPI communication data.
 * \param[in] step     The number of the MD time step.
 * \param[in] t        The time.
 * \param[in] ir       Structure containing MD input parameters, among those
 *                     also the structure needed for position swapping.
 * \param[in] wcycle   Count wallcycles of swap routines for diagnostic output.
 * \param[in] x        Positions of home particles this node owns.
 * \param[in] box      The simulation box.
 * \param[in] bVerbose Should we be quiet or verbose?
 * \param[in] bRerun   Are we doing a rerun?
 *
 * \returns Whether at least one pair of molecules was swapped.
 */
gmx_bool do_swapcoords(
        t_commrec        *cr,
        gmx_int64_t       step,
        double            t,
        t_inputrec       *ir,
        gmx_wallcycle    *wcycle,
        rvec              x[],
        matrix            box,
        gmx_bool          bVerbose,
        gmx_bool          bRerun);

#endif
