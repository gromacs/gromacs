/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */

#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-pp-communication.h"

/*! \brief Master PP-PME communication data structure */
struct gmx_pme_pp {
#if GMX_MPI
    MPI_Comm     mpi_comm_mysim; /**< MPI communicator for this simulation */
#endif
    int          nnode;          /**< The number of PP node to communicate with  */
    int         *node;           /**< The PP node ranks                          */
    int          node_peer;      /**< The peer PP node rank                      */
    int         *nat;            /**< The number of atom for each PP node        */
    //@{
    /**< Vectors of A- and B-state parameters used to transfer vectors to PME ranks  */
    real        *chargeA;
    real        *chargeB;
    real        *sqrt_c6A;
    real        *sqrt_c6B;
    real        *sigmaA;
    real        *sigmaB;
    //@}
    rvec        *x;             /**< Vector of atom coordinates to transfer to PME ranks */
    rvec        *f;             /**< Vector of atom forces received from PME ranks */
    int          nalloc;        /**< Allocation size of transfer vectors (>= \p nat) */
#if GMX_MPI
    //@{
    /**< Vectors of MPI objects used in non-blocking communication between multiple PP ranks per PME rank */
    MPI_Request *req;
    MPI_Status  *stat;
    //@}
#endif
};

/*! \brief Initialize the PME-only side of the PME <-> PP communication */
static gmx_pme_pp *gmx_pme_pp_init(t_commrec *cr)
{
    struct gmx_pme_pp *pme_pp;

    snew(pme_pp, 1);

#if GMX_MPI
    int rank;

    pme_pp->mpi_comm_mysim = cr->mpi_comm_mysim;
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);
    get_pme_ddnodes(cr, rank, &pme_pp->nnode, &pme_pp->node, &pme_pp->node_peer);
    snew(pme_pp->nat, pme_pp->nnode);
    snew(pme_pp->req, eCommType_NR*pme_pp->nnode);
    snew(pme_pp->stat, eCommType_NR*pme_pp->nnode);
    pme_pp->nalloc       = 0;
#else
    GMX_UNUSED_VALUE(cr);
#endif

    return pme_pp;
}

static void reset_pmeonly_counters(gmx_wallcycle_t wcycle,
                                   gmx_walltime_accounting_t walltime_accounting,
                                   t_nrnb *nrnb, t_inputrec *ir,
                                   gmx_int64_t step)
{
    /* Reset all the counters related to performance over the run */
    wallcycle_stop(wcycle, ewcRUN);
    wallcycle_reset_all(wcycle);
    init_nrnb(nrnb);
    if (ir->nsteps >= 0)
    {
        /* ir->nsteps is not used here, but we update it for consistency */
        ir->nsteps -= step - ir->init_step;
    }
    ir->init_step = step;
    wallcycle_start(wcycle, ewcRUN);
    walltime_accounting_start(walltime_accounting);
}


static void gmx_pmeonly_switch(int *npmedata, struct gmx_pme_t ***pmedata,
                               ivec grid_size,
                               real ewaldcoeff_q, real ewaldcoeff_lj,
                               t_commrec *cr, t_inputrec *ir,
                               struct gmx_pme_t **pme_ret)
{
    int               ind;
    struct gmx_pme_t *pme = nullptr;

    ind = 0;
    while (ind < *npmedata)
    {
        pme = (*pmedata)[ind];
        if (pme->nkx == grid_size[XX] &&
            pme->nky == grid_size[YY] &&
            pme->nkz == grid_size[ZZ])
        {
            /* Here we have found an existing PME data structure that suits us.
             * However, in the GPU case, we have to reinitialize it - there's only one GPU structure.
             * This should not cause actual GPU reallocations, at least (the allocated buffers are never shrunk).
             * So, just some grid size updates in the GPU kernel parameters.
             */
            gmx_pme_reinit(&((*pmedata)[ind]), cr, pme, ir, grid_size, ewaldcoeff_q, ewaldcoeff_lj);
            *pme_ret = pme;
            return;
        }

        ind++;
    }

    (*npmedata)++;
    srenew(*pmedata, *npmedata);

    /* Generate a new PME data structure, copying part of the old pointers */
    gmx_pme_reinit(&((*pmedata)[ind]), cr, pme, ir, grid_size, ewaldcoeff_q, ewaldcoeff_lj);

    *pme_ret = (*pmedata)[ind];
}

/*! \brief Called by PME-only ranks to receive coefficients and coordinates
 *
 * \param[in,out] pme_pp    PME-PP communication structure.
 * \param[out] natoms       Number of received atoms.
 * \param[out] box        System box, if received.
 * \param[out] maxshift_x        Maximum shift in X direction, if received.
 * \param[out] maxshift_y        Maximum shift in Y direction, if received.
 * \param[out] lambda_q         Free-energy lambda for electrostatics, if received.
 * \param[out] lambda_lj         Free-energy lambda for Lennard-Jones, if received.
 * \param[out] bEnerVir          Set to true if this is an energy/virial calculation step, otherwise set to false.
 * \param[out] step              MD integration step number.
 * \param[out] grid_size         PME grid size, if received.
 * \param[out] ewaldcoeff_q         Ewald cut-off parameter for electrostatics, if received.
 * \param[out] ewaldcoeff_lj         Ewald cut-off parameter for Lennard-Jones, if received.
 * \param[out] atomSetChanged    Set to true only if the local domain atom data (charges/coefficients)
 *                               has been received (after DD) and should be reinitialized. Otherwise not changed.
 *
 * \retval pmerecvqxX             All parameters were set, chargeA and chargeB can be NULL.
 * \retval pmerecvqxFINISH        No parameters were set.
 * \retval pmerecvqxSWITCHGRID    Only grid_size and *ewaldcoeff were set.
 * \retval pmerecvqxRESETCOUNTERS *step was set.
 */
static int gmx_pme_recv_coeffs_coords(gmx_pme_pp        *pme_pp,
                                      int               *natoms,
                                      matrix             box,
                                      int               *maxshift_x,
                                      int               *maxshift_y,
                                      real              *lambda_q,
                                      real              *lambda_lj,
                                      gmx_bool          *bEnerVir,
                                      gmx_int64_t       *step,
                                      ivec               grid_size,
                                      real              *ewaldcoeff_q,
                                      real              *ewaldcoeff_lj,
                                      bool              *atomSetChanged)
{
    int status = -1;
    int nat    = 0;

#if GMX_MPI
    unsigned int flags    = 0;
    int          messages = 0;

    do
    {
        gmx_pme_comm_n_box_t cnb;
        cnb.flags = 0;

        /* Receive the send count, box and time step from the peer PP node */
        MPI_Recv(&cnb, sizeof(cnb), MPI_BYTE,
                 pme_pp->node_peer, eCommType_CNB,
                 pme_pp->mpi_comm_mysim, MPI_STATUS_IGNORE);

        /* We accumulate all received flags */
        flags |= cnb.flags;

        *step  = cnb.step;

        if (debug)
        {
            fprintf(debug, "PME only rank receiving:%s%s%s%s%s\n",
                    (cnb.flags & PP_PME_CHARGE)        ? " charges" : "",
                    (cnb.flags & PP_PME_COORD )        ? " coordinates" : "",
                    (cnb.flags & PP_PME_FINISH)        ? " finish" : "",
                    (cnb.flags & PP_PME_SWITCHGRID)    ? " switch grid" : "",
                    (cnb.flags & PP_PME_RESETCOUNTERS) ? " reset counters" : "");
        }

        if (cnb.flags & PP_PME_FINISH)
        {
            status = pmerecvqxFINISH;
        }

        if (cnb.flags & PP_PME_SWITCHGRID)
        {
            /* Special case, receive the new parameters and return */
            copy_ivec(cnb.grid_size, grid_size);
            *ewaldcoeff_q  = cnb.ewaldcoeff_q;
            *ewaldcoeff_lj = cnb.ewaldcoeff_lj;

            status         = pmerecvqxSWITCHGRID;
        }

        if (cnb.flags & PP_PME_RESETCOUNTERS)
        {
            /* Special case, receive the step (set above) and return */
            status = pmerecvqxRESETCOUNTERS;
        }

        if (cnb.flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
        {
            *atomSetChanged = true;

            /* Receive the send counts from the other PP nodes */
            for (int sender = 0; sender < pme_pp->nnode; sender++)
            {
                if (pme_pp->node[sender] == pme_pp->node_peer)
                {
                    pme_pp->nat[sender] = cnb.natoms;
                }
                else
                {
                    MPI_Irecv(&(pme_pp->nat[sender]), sizeof(pme_pp->nat[0]),
                              MPI_BYTE,
                              pme_pp->node[sender], eCommType_CNB,
                              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                }
            }
            MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
            messages = 0;

            nat = 0;
            for (int sender = 0; sender < pme_pp->nnode; sender++)
            {
                nat += pme_pp->nat[sender];
            }

            if (nat > pme_pp->nalloc)
            {
                pme_pp->nalloc = over_alloc_dd(nat);
                if (cnb.flags & PP_PME_CHARGE)
                {
                    srenew(pme_pp->chargeA, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_CHARGEB)
                {
                    srenew(pme_pp->chargeB, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_SQRTC6)
                {
                    srenew(pme_pp->sqrt_c6A, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_SQRTC6B)
                {
                    srenew(pme_pp->sqrt_c6B, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_SIGMA)
                {
                    srenew(pme_pp->sigmaA, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_SIGMAB)
                {
                    srenew(pme_pp->sigmaB, pme_pp->nalloc);
                }
                srenew(pme_pp->x, pme_pp->nalloc);
                srenew(pme_pp->f, pme_pp->nalloc);
            }

            /* maxshift is sent when the charges are sent */
            *maxshift_x = cnb.maxshift_x;
            *maxshift_y = cnb.maxshift_y;

            /* Receive the charges in place */
            for (int q = 0; q < eCommType_NR; q++)
            {
                real *charge_pp;

                if (!(cnb.flags & (PP_PME_CHARGE<<q)))
                {
                    continue;
                }
                switch (q)
                {
                    case eCommType_ChargeA: charge_pp = pme_pp->chargeA;  break;
                    case eCommType_ChargeB: charge_pp = pme_pp->chargeB;  break;
                    case eCommType_SQRTC6A: charge_pp = pme_pp->sqrt_c6A; break;
                    case eCommType_SQRTC6B: charge_pp = pme_pp->sqrt_c6B; break;
                    case eCommType_SigmaA:  charge_pp = pme_pp->sigmaA;   break;
                    case eCommType_SigmaB:  charge_pp = pme_pp->sigmaB;   break;
                    default: gmx_incons("Wrong eCommType");
                }
                nat = 0;
                for (int sender = 0; sender < pme_pp->nnode; sender++)
                {
                    if (pme_pp->nat[sender] > 0)
                    {
                        MPI_Irecv(charge_pp+nat,
                                  pme_pp->nat[sender]*sizeof(real),
                                  MPI_BYTE,
                                  pme_pp->node[sender], q,
                                  pme_pp->mpi_comm_mysim,
                                  &pme_pp->req[messages++]);
                        nat += pme_pp->nat[sender];
                        if (debug)
                        {
                            fprintf(debug, "Received from PP rank %d: %d %s\n",
                                    pme_pp->node[sender], pme_pp->nat[sender],
                                    (q == eCommType_ChargeA ||
                                     q == eCommType_ChargeB) ? "charges" : "params");
                        }
                    }
                }
            }
        }

        if (cnb.flags & PP_PME_COORD)
        {
            /* The box, FE flag and lambda are sent along with the coordinates
             *  */
            copy_mat(cnb.box, box);
            *lambda_q       = cnb.lambda_q;
            *lambda_lj      = cnb.lambda_lj;
            *bEnerVir       = (cnb.flags & PP_PME_ENER_VIR);
            *step           = cnb.step;

            /* Receive the coordinates in place */
            nat = 0;
            for (int sender = 0; sender < pme_pp->nnode; sender++)
            {
                if (pme_pp->nat[sender] > 0)
                {
                    MPI_Irecv(pme_pp->x[nat], pme_pp->nat[sender]*sizeof(rvec),
                              MPI_BYTE,
                              pme_pp->node[sender], eCommType_COORD,
                              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                    nat += pme_pp->nat[sender];
                    if (debug)
                    {
                        fprintf(debug, "Received from PP rank %d: %d "
                                "coordinates\n",
                                pme_pp->node[sender], pme_pp->nat[sender]);
                    }
                }
            }

            status = pmerecvqxX;
        }

        /* Wait for the coordinates and/or charges to arrive */
        MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
        messages = 0;
    }
    while (status == -1);
#else
    GMX_UNUSED_VALUE(pme_pp);
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(maxshift_x);
    GMX_UNUSED_VALUE(maxshift_y);
    GMX_UNUSED_VALUE(lambda_q);
    GMX_UNUSED_VALUE(lambda_lj);
    GMX_UNUSED_VALUE(bEnerVir);
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(grid_size);
    GMX_UNUSED_VALUE(ewaldcoeff_q);
    GMX_UNUSED_VALUE(ewaldcoeff_lj);
    GMX_UNUSED_VALUE(atomSetChanged);

    status = pmerecvqxX;
#endif

    if (status == pmerecvqxX)
    {
        *natoms   = nat;
    }

    return status;
}

/*! \brief Send the PME mesh force, virial and energy to the PP-only ranks. */
static void gmx_pme_send_force_vir_ener(gmx_pme_pp *pme_pp,
                                        matrix vir_q, real energy_q,
                                        matrix vir_lj, real energy_lj,
                                        real dvdlambda_q, real dvdlambda_lj,
                                        float cycles)
{
#if GMX_MPI
    gmx_pme_comm_vir_ene_t cve;
    int                    messages, ind_start, ind_end;
    cve.cycles = cycles;

    /* Now the evaluated forces have to be transferred to the PP nodes */
    messages = 0;
    ind_end  = 0;
    for (int receiver = 0; receiver < pme_pp->nnode; receiver++)
    {
        ind_start = ind_end;
        ind_end   = ind_start + pme_pp->nat[receiver];
        if (MPI_Isend(pme_pp->f[ind_start], (ind_end-ind_start)*sizeof(rvec), MPI_BYTE,
                      pme_pp->node[receiver], 0,
                      pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]) != 0)
        {
            gmx_comm("MPI_Isend failed in do_pmeonly");
        }
    }

    /* send virial and energy to our last PP node */
    copy_mat(vir_q, cve.vir_q);
    copy_mat(vir_lj, cve.vir_lj);
    cve.energy_q     = energy_q;
    cve.energy_lj    = energy_lj;
    cve.dvdlambda_q  = dvdlambda_q;
    cve.dvdlambda_lj = dvdlambda_lj;
    /* check for the signals to send back to a PP node */
    cve.stop_cond = gmx_get_stop_condition();

    cve.cycles = cycles;

    if (debug)
    {
        fprintf(debug, "PME rank sending to PP rank %d: virial and energy\n",
                pme_pp->node_peer);
    }
    MPI_Isend(&cve, sizeof(cve), MPI_BYTE,
              pme_pp->node_peer, 1,
              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);

    /* Wait for the forces to arrive */
    MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
#else
    gmx_call("MPI not enabled");
    GMX_UNUSED_VALUE(pme_pp);
    GMX_UNUSED_VALUE(vir_q);
    GMX_UNUSED_VALUE(energy_q);
    GMX_UNUSED_VALUE(vir_lj);
    GMX_UNUSED_VALUE(energy_lj);
    GMX_UNUSED_VALUE(dvdlambda_q);
    GMX_UNUSED_VALUE(dvdlambda_lj);
    GMX_UNUSED_VALUE(cycles);
#endif
}

int gmx_pmeonly(struct gmx_pme_t *pme,
                t_commrec *cr,    t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                gmx_walltime_accounting_t walltime_accounting,
                t_inputrec *ir, PmeRunMode runMode)
{
    int                npmedata;
    struct gmx_pme_t **pmedata;
    gmx_pme_pp        *pme_pp;
    int                ret;
    int                natoms = 0;
    matrix             box;
    real               lambda_q   = 0;
    real               lambda_lj  = 0;
    int                maxshift_x = 0, maxshift_y = 0;
    real               energy_q, energy_lj, dvdlambda_q, dvdlambda_lj;
    matrix             vir_q, vir_lj;
    float              cycles;
    int                count;
    gmx_bool           bEnerVir = FALSE;
    gmx_int64_t        step;
    ivec               grid_switch;

    /* This data will only use with PME tuning, i.e. switching PME grids */
    npmedata = 1;
    snew(pmedata, npmedata);
    pmedata[0] = pme;

    pme_pp = gmx_pme_pp_init(cr);

    init_nrnb(mynrnb);

    count = 0;
    do /****** this is a quasi-loop over time steps! */
    {
        /* The reason for having a loop here is PME grid tuning/switching */
        do
        {
            /* Domain decomposition */
            bool atomSetChanged = false;
            real ewaldcoeff_q   = 0, ewaldcoeff_lj = 0;
            ret = gmx_pme_recv_coeffs_coords(pme_pp,
                                             &natoms,
                                             box,
                                             &maxshift_x, &maxshift_y,
                                             &lambda_q, &lambda_lj,
                                             &bEnerVir,
                                             &step,
                                             grid_switch,
                                             &ewaldcoeff_q,
                                             &ewaldcoeff_lj,
                                             &atomSetChanged);

            if (ret == pmerecvqxSWITCHGRID)
            {
                /* Switch the PME grid to grid_switch */
                gmx_pmeonly_switch(&npmedata, &pmedata, grid_switch, ewaldcoeff_q, ewaldcoeff_lj, cr, ir, &pme);
            }

            if (atomSetChanged)
            {
                gmx_pme_reinit_atoms(pme, natoms, pme_pp->chargeA);
            }

            if (ret == pmerecvqxRESETCOUNTERS)
            {
                /* Reset the cycle and flop counters */
                reset_pmeonly_counters(wcycle, walltime_accounting, mynrnb, ir, step);
            }
        }
        while (ret == pmerecvqxSWITCHGRID || ret == pmerecvqxRESETCOUNTERS);

        if (ret == pmerecvqxFINISH)
        {
            /* We should stop: break out of the loop */
            break;
        }

        if (count == 0)
        {
            wallcycle_start(wcycle, ewcRUN);
            walltime_accounting_start(walltime_accounting);
        }

        wallcycle_start(wcycle, ewcPMEMESH);

        dvdlambda_q  = 0;
        dvdlambda_lj = 0;
        clear_mat(vir_q);
        clear_mat(vir_lj);
        energy_q  = 0;
        energy_lj = 0;

        // TODO Make a struct of array refs onto these per-atom fields
        // of pme_pp (maybe box, energy and virial, too; and likewise
        // from mdatoms for the other call to gmx_pme_do), so we have
        // fewer lines of code and less parameter passing.
        const int pmeFlags = GMX_PME_DO_ALL_F | (bEnerVir ? GMX_PME_CALC_ENER_VIR : 0);
        if (runMode != PmeRunMode::CPU)
        {
            const bool boxChanged = true;
            //TODO this should be set properly by gmx_pme_recv_coeffs_coords,
            // or maybe use inputrecDynamicBox(ir), at the very least - change this when this codepath is tested!
            pme_gpu_prepare_step(pme, boxChanged, box, wcycle, pmeFlags);
            pme_gpu_launch_spread(pme, pme_pp->x, wcycle);
            pme_gpu_launch_complex_transforms(pme, wcycle);
            pme_gpu_launch_gather(pme, wcycle, pme_pp->f, PmeForceOutputHandling::Set);
            pme_gpu_wait_for_gpu(pme, wcycle, vir_q, &energy_q);
        }
        else
        {
            gmx_pme_do(pme, 0, natoms, pme_pp->x, pme_pp->f,
                       pme_pp->chargeA, pme_pp->chargeB,
                       pme_pp->sqrt_c6A, pme_pp->sqrt_c6B,
                       pme_pp->sigmaA, pme_pp->sigmaB, box,
                       cr, maxshift_x, maxshift_y, mynrnb, wcycle,
                       vir_q, vir_lj,
                       &energy_q, &energy_lj, lambda_q, lambda_lj, &dvdlambda_q, &dvdlambda_lj,
                       pmeFlags);
        }

        cycles = wallcycle_stop(wcycle, ewcPMEMESH);

        gmx_pme_send_force_vir_ener(pme_pp,
                                    vir_q, energy_q, vir_lj, energy_lj,
                                    dvdlambda_q, dvdlambda_lj, cycles);

        count++;
    } /***** end of quasi-loop, we stop with the break above */
    while (TRUE);

    walltime_accounting_end(walltime_accounting);

    return 0;
}
