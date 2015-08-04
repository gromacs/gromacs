/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file contains function definitions necessary for
 * managing the offload of long-ranged PME work to separate MPI rank,
 * for computing energies and forces (Coulomb and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/sighandler.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"

/*! \brief MPI Tags used to separate communication of different types of quantities */
enum {
    eCommType_ChargeA, eCommType_ChargeB, eCommType_SQRTC6A, eCommType_SQRTC6B,
    eCommType_SigmaA, eCommType_SigmaB, eCommType_NR, eCommType_COORD,
    eCommType_CNB
};

//@{
/*! \brief Flags used to coordinate PP-PME communication and computation phases
 *
 * Some parts of the code(gmx_pme_send_q, gmx_pme_recv_q_x) assume
 * that the six first flags are exactly in this order.
 * If more PP_PME_...-flags are to be introduced be aware of some of
 * the PME-specific flags in pme.h. Currently, they are also passed
 * through here.
 */

#define PP_PME_CHARGE         (1<<0)
#define PP_PME_CHARGEB        (1<<1)
#define PP_PME_SQRTC6         (1<<2)
#define PP_PME_SQRTC6B        (1<<3)
#define PP_PME_SIGMA          (1<<4)
#define PP_PME_SIGMAB         (1<<5)
#define PP_PME_COORD          (1<<6)
#define PP_PME_FEP_Q          (1<<7)
#define PP_PME_FEP_LJ         (1<<8)
#define PP_PME_ENER_VIR       (1<<9)
#define PP_PME_FINISH         (1<<10)
#define PP_PME_SWITCHGRID     (1<<11)
#define PP_PME_RESETCOUNTERS  (1<<12)

#define PME_PP_SIGSTOP        (1<<0)
#define PME_PP_SIGSTOPNSS     (1<<1)
//@}

/*! \brief Master PP-PME communication data structure */
struct gmx_pme_pp {
#ifdef GMX_MPI
    MPI_Comm     mpi_comm_mysim; /**< MPI communicator for this simulation */
#endif
    int          nnode;          /**< The number of PP node to communicate with  */
    int         *node;           /**< The PP node ranks                          */
    int          node_peer;      /**< The peer PP node rank                      */
    int         *nat;            /**< The number of atom for each PP node        */
    int          flags_charge;   /**< The flags sent along with the last charges */
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
#ifdef GMX_MPI
    //@{
    /**< Vectors of MPI objects used in non-blocking communication between multiple PP ranks per PME rank */
    MPI_Request *req;
    MPI_Status  *stat;
    //@}
#endif
};

/*! \brief Helper struct for PP-PME communication of parameters */
typedef struct gmx_pme_comm_n_box {
    int             natoms;     /**< Number of atoms */
    matrix          box;        /**< Box */
    int             maxshift_x; /**< Maximum shift in x direction */
    int             maxshift_y; /**< Maximum shift in y direction */
    real            lambda_q;   /**< Free-energy lambda for electrostatics */
    real            lambda_lj;  /**< Free-energy lambda for Lennard-Jones */
    int             flags;      /**< Control flags */
    gmx_int64_t     step;       /**< MD integration step number */
    //@{
    /*! \brief Used in PME grid tuning */
    ivec            grid_size;
    real            ewaldcoeff_q;
    real            ewaldcoeff_lj;
    //@}
} gmx_pme_comm_n_box_t;

/*! \brief Helper struct for PP-PME communication of virial and energy */
typedef struct {
    //@{
    /*! \brief Virial, energy, and derivative of potential w.r.t. lambda for charge and Lennard-Jones */
    matrix          vir_q;
    matrix          vir_lj;
    real            energy_q;
    real            energy_lj;
    real            dvdlambda_q;
    real            dvdlambda_lj;
    //@}
    float           cycles;     /**< Counter of CPU cycles used */
    gmx_stop_cond_t stop_cond;  /**< Flag used in responding to an external signal to terminate */
} gmx_pme_comm_vir_ene_t;

gmx_pme_pp_t gmx_pme_pp_init(t_commrec *cr)
{
    struct gmx_pme_pp *pme_pp;

    snew(pme_pp, 1);

#ifdef GMX_MPI
    int rank;

    pme_pp->mpi_comm_mysim = cr->mpi_comm_mysim;
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);
    get_pme_ddnodes(cr, rank, &pme_pp->nnode, &pme_pp->node, &pme_pp->node_peer);
    snew(pme_pp->nat, pme_pp->nnode);
    snew(pme_pp->req, eCommType_NR*pme_pp->nnode);
    snew(pme_pp->stat, eCommType_NR*pme_pp->nnode);
    pme_pp->nalloc       = 0;
    pme_pp->flags_charge = 0;
#else
    GMX_UNUSED_VALUE(cr);
#endif

    return pme_pp;
}

/*! \brief Block to wait for communication to PME ranks to complete
 *
 * This should be faster with a real non-blocking MPI implementation */
/* #define GMX_PME_DELAYED_WAIT */

static void gmx_pme_send_coeffs_coords_wait(gmx_domdec_t gmx_unused *dd)
{
#ifdef GMX_MPI
    if (dd->nreq_pme)
    {
        MPI_Waitall(dd->nreq_pme, dd->req_pme, MPI_STATUSES_IGNORE);
        dd->nreq_pme = 0;
    }
#endif
}

/*! \brief Send data to PME ranks */
static void gmx_pme_send_coeffs_coords(t_commrec *cr, int flags,
                                       real gmx_unused *chargeA, real gmx_unused *chargeB,
                                       real gmx_unused *c6A, real gmx_unused *c6B,
                                       real gmx_unused *sigmaA, real gmx_unused *sigmaB,
                                       matrix box, rvec gmx_unused *x,
                                       real lambda_q, real lambda_lj,
                                       int maxshift_x, int maxshift_y,
                                       gmx_int64_t step)
{
    gmx_domdec_t         *dd;
    gmx_pme_comm_n_box_t *cnb;
    int                   n;

    dd = cr->dd;
    n  = dd->nat_home;

    if (debug)
    {
        fprintf(debug, "PP rank %d sending to PME rank %d: %d%s%s%s%s\n",
                cr->sim_nodeid, dd->pme_nodeid, n,
                flags & PP_PME_CHARGE ? " charges" : "",
                flags & PP_PME_SQRTC6 ? " sqrtC6" : "",
                flags & PP_PME_SIGMA  ? " sigma" : "",
                flags & PP_PME_COORD  ? " coordinates" : "");
    }

#ifdef GMX_PME_DELAYED_WAIT
    /* When can not use cnb until pending communication has finished */
    gmx_pme_send_coeffs_coords_wait(dd);
#endif

    if (dd->pme_receive_vir_ener)
    {
        /* Peer PP node: communicate all data */
        if (dd->cnb == NULL)
        {
            snew(dd->cnb, 1);
        }
        cnb = dd->cnb;

        cnb->flags      = flags;
        cnb->natoms     = n;
        cnb->maxshift_x = maxshift_x;
        cnb->maxshift_y = maxshift_y;
        cnb->lambda_q   = lambda_q;
        cnb->lambda_lj  = lambda_lj;
        cnb->step       = step;
        if (flags & PP_PME_COORD)
        {
            copy_mat(box, cnb->box);
        }
#ifdef GMX_MPI
        MPI_Isend(cnb, sizeof(*cnb), MPI_BYTE,
                  dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }
    else if (flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
    {
#ifdef GMX_MPI
        /* Communicate only the number of atoms */
        MPI_Isend(&n, sizeof(n), MPI_BYTE,
                  dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }

#ifdef GMX_MPI
    if (n > 0)
    {
        if (flags & PP_PME_CHARGE)
        {
            MPI_Isend(chargeA, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_ChargeA, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_CHARGEB)
        {
            MPI_Isend(chargeB, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_ChargeB, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SQRTC6)
        {
            MPI_Isend(c6A, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_SQRTC6A, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SQRTC6B)
        {
            MPI_Isend(c6B, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_SQRTC6B, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMA)
        {
            MPI_Isend(sigmaA, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_SigmaA, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMAB)
        {
            MPI_Isend(sigmaB, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, eCommType_SigmaB, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_COORD)
        {
            MPI_Isend(x[0], n*sizeof(rvec), MPI_BYTE,
                      dd->pme_nodeid, eCommType_COORD, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
    }

#ifndef GMX_PME_DELAYED_WAIT
    /* Wait for the data to arrive */
    /* We can skip this wait as we are sure x and q will not be modified
     * before the next call to gmx_pme_send_x_q or gmx_pme_receive_f.
     */
    gmx_pme_send_coeffs_coords_wait(dd);
#endif
#endif
}

void gmx_pme_send_parameters(t_commrec *cr,
                             const interaction_const_t *ic,
                             gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                             real *chargeA, real *chargeB,
                             real *sqrt_c6A, real *sqrt_c6B,
                             real *sigmaA, real *sigmaB,
                             int maxshift_x, int maxshift_y)
{
    int flags;

    flags = 0;
    if (EEL_PME(ic->eeltype))
    {
        flags |= PP_PME_CHARGE;
    }
    if (EVDW_PME(ic->vdwtype))
    {
        flags |= (PP_PME_SQRTC6 | PP_PME_SIGMA);
    }
    if (bFreeEnergy_q || bFreeEnergy_lj)
    {
        /* Assumes that the B state flags are in the bits just above
         * the ones for the A state. */
        flags |= (flags << 1);
    }

    gmx_pme_send_coeffs_coords(cr, flags,
                               chargeA, chargeB,
                               sqrt_c6A, sqrt_c6B, sigmaA, sigmaB,
                               NULL, NULL, 0, 0, maxshift_x, maxshift_y, -1);
}

void gmx_pme_send_coordinates(t_commrec *cr, matrix box, rvec *x,
                              gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                              real lambda_q, real lambda_lj,
                              gmx_bool bEnerVir, int pme_flags,
                              gmx_int64_t step)
{
    int flags;

    flags = pme_flags | PP_PME_COORD;
    if (bFreeEnergy_q)
    {
        flags |= PP_PME_FEP_Q;
    }
    if (bFreeEnergy_lj)
    {
        flags |= PP_PME_FEP_LJ;
    }
    if (bEnerVir)
    {
        flags |= PP_PME_ENER_VIR;
    }
    gmx_pme_send_coeffs_coords(cr, flags, NULL, NULL, NULL, NULL, NULL, NULL,
                               box, x, lambda_q, lambda_lj, 0, 0, step);
}

void gmx_pme_send_finish(t_commrec *cr)
{
    int flags;

    flags = PP_PME_FINISH;

    gmx_pme_send_coeffs_coords(cr, flags, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, -1);
}

void gmx_pme_send_switchgrid(t_commrec gmx_unused *cr,
                             ivec gmx_unused       grid_size,
                             real gmx_unused       ewaldcoeff_q,
                             real gmx_unused       ewaldcoeff_lj)
{
#ifdef GMX_MPI
    gmx_pme_comm_n_box_t cnb;

    /* Only let one PP node signal each PME node */
    if (cr->dd->pme_receive_vir_ener)
    {
        cnb.flags = PP_PME_SWITCHGRID;
        copy_ivec(grid_size, cnb.grid_size);
        cnb.ewaldcoeff_q  = ewaldcoeff_q;
        cnb.ewaldcoeff_lj = ewaldcoeff_lj;

        /* We send this, uncommon, message blocking to simplify the code */
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE,
                 cr->dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim);
    }
#endif
}

void gmx_pme_send_resetcounters(t_commrec gmx_unused *cr, gmx_int64_t gmx_unused step)
{
#ifdef GMX_MPI
    gmx_pme_comm_n_box_t cnb;

    /* Only let one PP node signal each PME node */
    if (cr->dd->pme_receive_vir_ener)
    {
        cnb.flags = PP_PME_RESETCOUNTERS;
        cnb.step  = step;

        /* We send this, uncommon, message blocking to simplify the code */
        MPI_Send(&cnb, sizeof(cnb), MPI_BYTE,
                 cr->dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim);
    }
#endif
}

int gmx_pme_recv_coeffs_coords(struct gmx_pme_pp *pme_pp,
                               int               *natoms,
                               real             **chargeA,
                               real             **chargeB,
                               real             **sqrt_c6A,
                               real             **sqrt_c6B,
                               real             **sigmaA,
                               real             **sigmaB,
                               matrix             box,
                               rvec             **x,
                               rvec             **f,
                               int               *maxshift_x,
                               int               *maxshift_y,
                               gmx_bool          *bFreeEnergy_q,
                               gmx_bool          *bFreeEnergy_lj,
                               real              *lambda_q,
                               real              *lambda_lj,
                               gmx_bool          *bEnerVir,
                               int               *pme_flags,
                               gmx_int64_t       *step,
                               ivec               grid_size,
                               real              *ewaldcoeff_q,
                               real              *ewaldcoeff_lj)
{
    int                  nat = 0, status;

    *pme_flags = 0;
#ifdef GMX_MPI
    gmx_pme_comm_n_box_t cnb;
    int                  messages;

    cnb.flags  = 0;
    messages   = 0;
    do
    {

        /* Receive the send count, box and time step from the peer PP node */
        MPI_Recv(&cnb, sizeof(cnb), MPI_BYTE,
                 pme_pp->node_peer, eCommType_CNB,
                 pme_pp->mpi_comm_mysim, MPI_STATUS_IGNORE);

        if (debug)
        {
            fprintf(debug, "PME only rank receiving:%s%s%s%s%s\n",
                    (cnb.flags & PP_PME_CHARGE)        ? " charges" : "",
                    (cnb.flags & PP_PME_COORD )        ? " coordinates" : "",
                    (cnb.flags & PP_PME_FINISH)        ? " finish" : "",
                    (cnb.flags & PP_PME_SWITCHGRID)    ? " switch grid" : "",
                    (cnb.flags & PP_PME_RESETCOUNTERS) ? " reset counters" : "");
        }

        if (cnb.flags & PP_PME_SWITCHGRID)
        {
            /* Special case, receive the new parameters and return */
            copy_ivec(cnb.grid_size, grid_size);
            *ewaldcoeff_q  = cnb.ewaldcoeff_q;
            *ewaldcoeff_lj = cnb.ewaldcoeff_lj;
            return pmerecvqxSWITCHGRID;
        }

        if (cnb.flags & PP_PME_RESETCOUNTERS)
        {
            /* Special case, receive the step and return */
            *step = cnb.step;

            return pmerecvqxRESETCOUNTERS;
        }

        if (cnb.flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
        {
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

            pme_pp->flags_charge = cnb.flags;
        }

        if (cnb.flags & PP_PME_COORD)
        {
            if (!(pme_pp->flags_charge & (PP_PME_CHARGE | PP_PME_SQRTC6)))
            {
                gmx_incons("PME-only rank received coordinates before charges and/or C6-values"
                           );
            }

            /* The box, FE flag and lambda are sent along with the coordinates
             *  */
            copy_mat(cnb.box, box);
            *bFreeEnergy_q  = ((cnb.flags & GMX_PME_DO_COULOMB) &&
                               (cnb.flags & PP_PME_FEP_Q));
            *bFreeEnergy_lj = ((cnb.flags & GMX_PME_DO_LJ) &&
                               (cnb.flags & PP_PME_FEP_LJ));
            *lambda_q       = cnb.lambda_q;
            *lambda_lj      = cnb.lambda_lj;
            *bEnerVir       = (cnb.flags & PP_PME_ENER_VIR);
            *pme_flags      = cnb.flags;

            if (*bFreeEnergy_q && !(pme_pp->flags_charge & PP_PME_CHARGEB))
            {
                gmx_incons("PME-only rank received free energy request, but "
                           "did not receive B-state charges");
            }

            if (*bFreeEnergy_lj && !(pme_pp->flags_charge & PP_PME_SQRTC6B))
            {
                gmx_incons("PME-only rank received free energy request, but "
                           "did not receive B-state C6-values");
            }

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
        }

        /* Wait for the coordinates and/or charges to arrive */
        MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
        messages = 0;
    }
    while (!(cnb.flags & (PP_PME_COORD | PP_PME_FINISH)));
    status = ((cnb.flags & PP_PME_FINISH) ? pmerecvqxFINISH : pmerecvqxX);

    *step = cnb.step;
#else
    GMX_UNUSED_VALUE(box);
    GMX_UNUSED_VALUE(maxshift_x);
    GMX_UNUSED_VALUE(maxshift_y);
    GMX_UNUSED_VALUE(bFreeEnergy_q);
    GMX_UNUSED_VALUE(bFreeEnergy_lj);
    GMX_UNUSED_VALUE(lambda_q);
    GMX_UNUSED_VALUE(lambda_lj);
    GMX_UNUSED_VALUE(bEnerVir);
    GMX_UNUSED_VALUE(step);
    GMX_UNUSED_VALUE(grid_size);
    GMX_UNUSED_VALUE(ewaldcoeff_q);
    GMX_UNUSED_VALUE(ewaldcoeff_lj);

    status = pmerecvqxX;
#endif

    *natoms   = nat;
    *chargeA  = pme_pp->chargeA;
    *chargeB  = pme_pp->chargeB;
    *sqrt_c6A = pme_pp->sqrt_c6A;
    *sqrt_c6B = pme_pp->sqrt_c6B;
    *sigmaA   = pme_pp->sigmaA;
    *sigmaB   = pme_pp->sigmaB;
    *x        = pme_pp->x;
    *f        = pme_pp->f;

    return status;
}
/*! \brief Receive virial and energy from PME rank */
static void receive_virial_energy(t_commrec *cr,
                                  matrix vir_q, real *energy_q,
                                  matrix vir_lj, real *energy_lj,
                                  real *dvdlambda_q, real *dvdlambda_lj,
                                  float *pme_cycles)
{
    gmx_pme_comm_vir_ene_t cve;

    if (cr->dd->pme_receive_vir_ener)
    {
        if (debug)
        {
            fprintf(debug,
                    "PP rank %d receiving from PME rank %d: virial and energy\n",
                    cr->sim_nodeid, cr->dd->pme_nodeid);
        }
#ifdef GMX_MPI
        MPI_Recv(&cve, sizeof(cve), MPI_BYTE, cr->dd->pme_nodeid, 1, cr->mpi_comm_mysim,
                 MPI_STATUS_IGNORE);
#else
        memset(&cve, 0, sizeof(cve));
#endif

        m_add(vir_q, cve.vir_q, vir_q);
        m_add(vir_lj, cve.vir_lj, vir_lj);
        *energy_q      = cve.energy_q;
        *energy_lj     = cve.energy_lj;
        *dvdlambda_q  += cve.dvdlambda_q;
        *dvdlambda_lj += cve.dvdlambda_lj;
        *pme_cycles    = cve.cycles;

        if (cve.stop_cond != gmx_stop_cond_none)
        {
            gmx_set_stop_condition(cve.stop_cond);
        }
    }
    else
    {
        *energy_q   = 0;
        *energy_lj  = 0;
        *pme_cycles = 0;
    }
}

void gmx_pme_receive_f(t_commrec *cr,
                       rvec f[], matrix vir_q, real *energy_q,
                       matrix vir_lj, real *energy_lj,
                       real *dvdlambda_q, real *dvdlambda_lj,
                       float *pme_cycles)
{
    int natoms, i;

#ifdef GMX_PME_DELAYED_WAIT
    /* Wait for the x request to finish */
    gmx_pme_send_coeffs_coords_wait(cr->dd);
#endif

    natoms = cr->dd->nat_home;

    if (natoms > cr->dd->pme_recv_f_alloc)
    {
        cr->dd->pme_recv_f_alloc = over_alloc_dd(natoms);
        srenew(cr->dd->pme_recv_f_buf, cr->dd->pme_recv_f_alloc);
    }

#ifdef GMX_MPI
    MPI_Recv(cr->dd->pme_recv_f_buf[0],
             natoms*sizeof(rvec), MPI_BYTE,
             cr->dd->pme_nodeid, 0, cr->mpi_comm_mysim,
             MPI_STATUS_IGNORE);
#endif

    for (i = 0; i < natoms; i++)
    {
        rvec_inc(f[i], cr->dd->pme_recv_f_buf[i]);
    }


    receive_virial_energy(cr, vir_q, energy_q, vir_lj, energy_lj, dvdlambda_q, dvdlambda_lj, pme_cycles);
}

void gmx_pme_send_force_vir_ener(struct gmx_pme_pp *pme_pp,
                                 rvec gmx_unused *f,
                                 matrix vir_q, real energy_q,
                                 matrix vir_lj, real energy_lj,
                                 real dvdlambda_q, real dvdlambda_lj,
                                 float cycles)
{
#ifdef GMX_MPI
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
        if (MPI_Isend(f[ind_start], (ind_end-ind_start)*sizeof(rvec), MPI_BYTE,
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
    GMX_UNUSED_VALUE(f);
    GMX_UNUSED_VALUE(vir_q);
    GMX_UNUSED_VALUE(energy_q);
    GMX_UNUSED_VALUE(vir_lj);
    GMX_UNUSED_VALUE(energy_lj);
    GMX_UNUSED_VALUE(dvdlambda_q);
    GMX_UNUSED_VALUE(dvdlambda_lj);
    GMX_UNUSED_VALUE(cycles);
#endif
}
