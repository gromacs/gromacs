/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "typedefs.h"
#include "types/commrec.h"
#include "enums.h" //should this be merged in here?
#include "gromacs/timing/wallcycle.h"

#ifdef __cplusplus
extern "C" {
#endif


/*! \brief Initialize ion / water position swapping ("Computational Electrophysiology").
 *
 * This routine does the memory allocation for various helper arrays, opens
 * the output file, sets up swap data checkpoint writing, etc.
 *
 * \param[in] fplog         General output file, normally md.log.
 * \param[in] bVerbose      Should we be quiet or verbose?
 * \param[in] ir            Structure containing MD input parameters, among those
 *                          also the structure needed for position swapping.
 * \param[in] fn            Output file name for swap data.
 * \param[in] mtop          Molecular topology.
 * \param[in] x             The initial positions of all particles.
 * \param[in] box           The simulation box.
 * \param[in] swapstate     Swap-related data that is read from or written to checkpoint.
 * \param[in] cr            Pointer to MPI communication data.
 * \param[in] oenv          Needed to open the swap output XVGR file.
 * \param[in] Flags         Flags passed over from main, used to determine
 *                          whether we are doing a rerun, appending, etc.
 */
extern void init_swapcoords(
        FILE              *fplog,
        gmx_bool           bVerbose,
        t_inputrec        *ir,
        const char        *fn,
        gmx_mtop_t        *mtop,
        rvec               x[],
        matrix             box,
        swapstate_t       *swapstate,
        t_commrec         *cr,
        const output_env_t oenv,
        unsigned long      Flags);


/*! \brief Make a selection of the home atoms for the swap groups. These are
 * the ions, the water, and the channels. This routine should be called at every
 * domain decomposition.
 *
 * \param[in] dd            Structure containing domain decomposition data.
 * \param[in] si_pub        Pointer to the swap data structure.
 */
extern void dd_make_local_swap_groups(gmx_domdec_t *dd, t_swapcoords *si_pub);


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
 * \param[in] mtop     Molecular topology.
 * \param[in] bVerbose Should we be quiet or verbose?
 * \param[in] bRerun   Are we doing a rerun?
 *
 * \returns Whether at least one pair of molecules was swapped.
 */
extern gmx_bool do_swapcoords(
        t_commrec        *cr,
        gmx_int64_t       step,
        double            t,
        t_inputrec       *ir,
        gmx_wallcycle_t   wcycle,
        rvec              x[],
        matrix            box,
        gmx_mtop_t       *mtop,
        gmx_bool          bVerbose,
        gmx_bool          bRerun);


struct t_swapcoords
{
    int              nstswap;           /* Every how many steps a swap is attempted?    */
    int              nat;               /* Number of atoms in the ion group             */
    int              nat_split[2];      /* Number of atoms in the split group           */
    int              nat_sol;           /* Number of atoms in the solvent group         */
    atom_id         *ind;               /* The global ion group atoms numbers           */
    atom_id         *ind_split[2];      /* Split groups for compartment partitioning    */
    atom_id         *ind_sol;           /* The global solvent group atom numbers        */
    gmx_bool         massw_split[2];    /* Use mass-weighted positions in split group?  */
    real             cyl0r, cyl1r;      /* Split cylinders defined by radius, upper and */
    real             cyl0u, cyl1u;      /* ... lower extension. The split cylinders de- */
    real             cyl0l, cyl1l;      /* ... fine the channels and are each anchored  */
                                        /* ... in the center of the split group         */
    int              nanions[eCompNR];  /* Requested number of anions and               */
    int              nAverage;          /* Coupling constant (nr of swap attempt steps) */
    real             threshold;         /* Ion counts may deviate from the requested
                                           values by +-threshold before a swap is done  */
    int              ncations[eCompNR]; /* ... cations for both compartments            */
    gmx_swapcoords_t si_priv;           /* swap private data accessible in
                                         * swapcoords.c                                 */
};

struct swapstate_t
{
    int        eSwapCoords;                         /* Swapping along x, y, or z-direction?      */
    int        nat_req[eCompNR][eIonNR];            /* Requested ion numbers per type an comp.   */
    int       *nat_req_p[eCompNR][eIonNR];          /* Pointer to this data (for .cpt writing)   */
    int        nAverage;                            /* Use average over this many swap attempt
                                                       steps when determining the ion counts     */
    int        inflow_netto[eCompNR][eIonNR];       /* Flux determined from the # of swaps       */
    int       *inflow_netto_p[eCompNR][eIonNR];     /* Pointer to this data                      */
    int       *nat_past[eCompNR][eIonNR];           /* Array with nAverage entries for history   */
    int       *nat_past_p[eCompNR][eIonNR];         /* Pointer points to the first entry only    */

    /* Channel flux detection, this is counting only and has no influence on whether swaps
     * are performed or not: */
    int            fluxfromAtoB[eCompNR][eIonNR];   /* Flux determined from the split cylinders  */
    int           *fluxfromAtoB_p[eCompNR][eIonNR]; /* Pointer to this data                      */
    int           *fluxleak;                        /* Flux not going through any channel        */
    int            nions;                           /* Size of the following arrays              */
    unsigned char *comp_from;                       /* Ion came from which compartment?          */
    unsigned char *channel_label;                   /* Through which channel did this ion pass?  */

    /* To also make multimeric channel proteins whole, we save the last whole configuration of
     * the channels in the checkpoint file. If we have no checkpoint file, we assume that the
     * starting configuration hast the correct PBC representation after making the individual
     * molecules whole */
    gmx_bool    bFromCpt;                           /* Did we started from a checkpoint file?    */
    int         nat[eChanNR];                       /* Size of xc_old_whole, i.e. the number of
                                                       atoms in each channel                     */
    rvec       *xc_old_whole[eChanNR];              /* Last known whole positions of the two
                                                       channels (important for multimeric ch.!)  */
    rvec      **xc_old_whole_p[eChanNR];            /* Pointer to these positions                */
};

#ifdef __cplusplus
}
#endif

#endif
