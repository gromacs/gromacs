/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "pme.h"
#include "network.h"
#include "domdec.h"
#include "sighandler.h"

#include "gromacs/utility/gmxmpi.h"

/* Some parts of the code(gmx_pme_send_q, gmx_pme_recv_q_x) assume
 * that the six first flags are exactly in this order.
 * If more PP_PME_...-flags are to be introduced be aware of some of
 * the PME-specific flags in pme.h. Currently, they are also passed
 * through here.
 */

#define PP_PME_CHARGE         (1<<0)
#define PP_PME_CHARGEB        (1<<1)
#define PP_PME_C6             (1<<2)
#define PP_PME_C6B            (1<<3)
#define PP_PME_SIGMA          (1<<4)
#define PP_PME_SIGMAB         (1<<5)
#define PP_PME_COORD          (1<<6)
#define PP_PME_FEP            (1<<7)
#define PP_PME_ENER_VIR       (1<<8)
#define PP_PME_FINISH         (1<<9)
#define PP_PME_SWITCHGRID     (1<<10)
#define PP_PME_RESETCOUNTERS  (1<<11)

#define PME_PP_SIGSTOP        (1<<0)
#define PME_PP_SIGSTOPNSS     (1<<1)

typedef struct gmx_pme_pp {
#ifdef GMX_MPI
    MPI_Comm     mpi_comm_mysim;
#endif
    int          nnode;        /* The number of PP node to communicate with  */
    int         *node;         /* The PP node ranks                          */
    int          node_peer;    /* The peer PP node rank                      */
    int         *nat;          /* The number of atom for each PP node        */
    int          flags_charge; /* The flags sent along with the last charges */
    real        *chargeA;
    real        *chargeB;
    real        *c6A;
    real        *c6B;
    real        *sigmaA;
    real        *sigmaB;
    rvec        *x;
    rvec        *f;
    int          nalloc;
#ifdef GMX_MPI
    MPI_Request *req;
    MPI_Status  *stat;
#endif
} t_gmx_pme_pp;

typedef struct gmx_pme_comm_n_box {
    int             natoms;
    matrix          box;
    int             maxshift_x;
    int             maxshift_y;
    real            lambda_q;
    real            lambda_lj;
    int             flags;
    gmx_large_int_t step;
    ivec            grid_size;    /* For PME grid tuning */
    real            ewaldcoeff_q; /* For PME grid tuning */
    real            ewaldcoeff_lj;
} gmx_pme_comm_n_box_t;

typedef struct {
    matrix          vir_q;
    matrix          vir_lj;
    real            energy_q;
    real            energy_lj;
    real            dvdlambda_q;
    real            dvdlambda_lj;
    float           cycles;
    gmx_stop_cond_t stop_cond;
} gmx_pme_comm_vir_ene_t;




gmx_pme_pp_t gmx_pme_pp_init(t_commrec *cr)
{
    struct gmx_pme_pp *pme_pp;
    int                rank;

    snew(pme_pp, 1);

#ifdef GMX_MPI
    pme_pp->mpi_comm_mysim = cr->mpi_comm_mysim;
    MPI_Comm_rank(cr->mpi_comm_mygroup, &rank);
    get_pme_ddnodes(cr, rank, &pme_pp->nnode, &pme_pp->node, &pme_pp->node_peer);
    snew(pme_pp->nat, pme_pp->nnode);
    snew(pme_pp->req, 6*pme_pp->nnode);
    snew(pme_pp->stat, 6*pme_pp->nnode);
    pme_pp->nalloc       = 0;
    pme_pp->flags_charge = 0;
#endif

    return pme_pp;
}

/* This should be faster with a real non-blocking MPI implementation */
/* #define GMX_PME_DELAYED_WAIT */

static void gmx_pme_send_q_x_wait(gmx_domdec_t *dd)
{
#ifdef GMX_MPI
    if (dd->nreq_pme)
    {
        MPI_Waitall(dd->nreq_pme, dd->req_pme, MPI_STATUSES_IGNORE);
        dd->nreq_pme = 0;
    }
#endif
}

static void gmx_pme_send_q_x(t_commrec *cr, int flags,
                             real *chargeA, real *chargeB,
                             real *c6A, real *c6B,
                             real *sigmaA, real *sigmaB,
                             matrix box, rvec *x,
                             real lambda_q, real lambda_lj,
                             int maxshift_x, int maxshift_y,
                             gmx_large_int_t step)
{
    gmx_domdec_t         *dd;
    gmx_pme_comm_n_box_t *cnb;
    int                   n;

    dd = cr->dd;
    n  = dd->nat_home;

    if (debug)
    {
        fprintf(debug, "PP node %d sending to PME node %d: %d%s%s\n",
                cr->sim_nodeid, dd->pme_nodeid, n,
                flags & PP_PME_CHARGE ? " charges" : "",
                flags & PP_PME_COORD  ? " coordinates" : "");
    }

#ifdef GMX_PME_DELAYED_WAIT
    /* When can not use cnb until pending communication has finished */
    gmx_pme_send_x_q_wait(dd);
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
                  dd->pme_nodeid, 0, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }
    else if (flags & (PP_PME_CHARGE | PP_PME_C6 | PP_PME_SIGMA))
    {
#ifdef GMX_MPI
        /* Communicate only the number of atoms */
        MPI_Isend(&n, sizeof(n), MPI_BYTE,
                  dd->pme_nodeid, 0, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }

#ifdef GMX_MPI
    if (n > 0)
    {
        if (flags & PP_PME_CHARGE)
        {
            MPI_Isend(chargeA, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 1, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_CHARGEB)
        {
            MPI_Isend(chargeB, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 2, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_C6)
        {
            MPI_Isend(c6A, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 3, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_C6B)
        {
            MPI_Isend(c6B, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 4, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMA)
        {
            MPI_Isend(sigmaA, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 5, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_SIGMAB)
        {
            MPI_Isend(sigmaB, n*sizeof(real), MPI_BYTE,
                      dd->pme_nodeid, 6, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
        if (flags & PP_PME_COORD)
        {
            MPI_Isend(x[0], n*sizeof(rvec), MPI_BYTE,
                      dd->pme_nodeid, 7, cr->mpi_comm_mysim,
                      &dd->req_pme[dd->nreq_pme++]);
        }
    }

#ifndef GMX_PME_DELAYED_WAIT
    /* Wait for the data to arrive */
    /* We can skip this wait as we are sure x and q will not be modified
     * before the next call to gmx_pme_send_x_q or gmx_pme_receive_f.
     */
    gmx_pme_send_q_x_wait(dd);
#endif
#endif
}

void gmx_pme_send_q(t_commrec *cr,
                    gmx_bool bFreeEnergy, real *chargeA, real *chargeB,
                    real *c6A, real *c6B, real *sigmaA, real *sigmaB,
                    int maxshift_x, int maxshift_y)
{
    int flags;

    flags = PP_PME_CHARGE | PP_PME_C6 | PP_PME_SIGMA;
    if (bFreeEnergy)
    {
        /* Assumes that the B state flags are in the bits just above
         * the ones for the A state. */
        flags |= (flags << 1);
    }

    gmx_pme_send_q_x(cr, flags, chargeA, chargeB, c6A, c6B, sigmaA, sigmaB,
                     NULL, NULL, 0, 0, maxshift_x, maxshift_y, -1);
}

void gmx_pme_send_x(t_commrec *cr, matrix box, rvec *x,
                    gmx_bool bFreeEnergy,
                    real lambda_q, real lambda_lj,
                    gmx_bool bEnerVir, int pme_flags,
                    gmx_large_int_t step)
{
    int flags;

    flags = pme_flags | PP_PME_COORD;
    if (bFreeEnergy)
    {
        flags |= PP_PME_FEP;
    }
    if (bEnerVir)
    {
        flags |= PP_PME_ENER_VIR;
    }
    gmx_pme_send_q_x(cr, flags, NULL, NULL, NULL, NULL, NULL, NULL,
                     box, x, lambda_q, lambda_lj, 0, 0, step);
}

void gmx_pme_send_finish(t_commrec *cr)
{
    int flags;

    flags = PP_PME_FINISH;

    gmx_pme_send_q_x(cr, flags, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0, -1);
}

void gmx_pme_send_switchgrid(t_commrec *cr, ivec grid_size, real ewaldcoeff_q, real ewaldcoeff_lj)
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
                 cr->dd->pme_nodeid, 0, cr->mpi_comm_mysim);
    }
#endif
}

void gmx_pme_send_resetcounters(t_commrec *cr, gmx_large_int_t step)
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
                 cr->dd->pme_nodeid, 0, cr->mpi_comm_mysim);
    }
#endif
}

int gmx_pme_recv_q_x(struct gmx_pme_pp *pme_pp,
                     int *natoms,
                     real **chargeA, real **chargeB,
                     real **c6A, real **c6B,
                     real **sigmaA, real **sigmaB,
                     matrix box, rvec **x, rvec **f,
                     int *maxshift_x, int *maxshift_y,
                     gmx_bool *bFreeEnergy,
                     real *lambda_q, real *lambda_lj,
                     gmx_bool *bEnerVir, int *pme_flags,
                     gmx_large_int_t *step,
                     ivec grid_size, real *ewaldcoeff_q, real *ewaldcoeff_lj)
{
    gmx_pme_comm_n_box_t cnb;
    int                  nat = 0, q, messages, sender;
    real                *charge_pp;

    messages = 0;

    /* avoid compiler warning about unused variable without MPI support */
    cnb.flags  = 0;
    *pme_flags = 0;
#ifdef GMX_MPI
    do
    {
        /* Receive the send count, box and time step from the peer PP node */
        MPI_Recv(&cnb, sizeof(cnb), MPI_BYTE,
                 pme_pp->node_peer, 0,
                 pme_pp->mpi_comm_mysim, MPI_STATUS_IGNORE);

        if (debug)
        {
            fprintf(debug, "PME only node receiving:%s%s%s%s%s\n",
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

        if (cnb.flags & (PP_PME_CHARGE | PP_PME_C6 | PP_PME_SIGMA))
        {
            /* Receive the send counts from the other PP nodes */
            for (sender = 0; sender < pme_pp->nnode; sender++)
            {
                if (pme_pp->node[sender] == pme_pp->node_peer)
                {
                    pme_pp->nat[sender] = cnb.natoms;
                }
                else
                {
                    MPI_Irecv(&(pme_pp->nat[sender]), sizeof(pme_pp->nat[0]),
                              MPI_BYTE,
                              pme_pp->node[sender], 0,
                              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                }
            }
            MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
            messages = 0;

            nat = 0;
            for (sender = 0; sender < pme_pp->nnode; sender++)
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
                if (cnb.flags & PP_PME_C6)
                {
                    srenew(pme_pp->c6A, pme_pp->nalloc);
                }
                if (cnb.flags & PP_PME_C6B)
                {
                    srenew(pme_pp->c6B, pme_pp->nalloc);
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
            for (q = 0; q < 6; q++)
            {
                if (!(cnb.flags & (PP_PME_CHARGE<<q)))
                {
                    continue;
                }
                switch (q)
                {
                    case 0: charge_pp = pme_pp->chargeA; break;
                    case 1: charge_pp = pme_pp->chargeB; break;
                    case 2: charge_pp = pme_pp->c6A;     break;
                    case 3: charge_pp = pme_pp->c6B;     break;
                    case 4: charge_pp = pme_pp->sigmaA;  break;
                    case 5: charge_pp = pme_pp->sigmaB;  break;
                }
                nat = 0;
                for (sender = 0; sender < pme_pp->nnode; sender++)
                {
                    if (pme_pp->nat[sender] > 0)
                    {
                        MPI_Irecv(charge_pp+nat,
                                  pme_pp->nat[sender]*sizeof(real),
                                  MPI_BYTE,
                                  pme_pp->node[sender], 1+q,
                                  pme_pp->mpi_comm_mysim,
                                  &pme_pp->req[messages++]);
                        nat += pme_pp->nat[sender];
                        if (debug)
                        {
                            fprintf(debug, "Received from PP node %d: %d "
                                    "charges\n",
                                    pme_pp->node[sender], pme_pp->nat[sender]);
                        }
                    }
                }
            }

            pme_pp->flags_charge = cnb.flags;
        }

        if (cnb.flags & PP_PME_COORD)
        {
            if (!(pme_pp->flags_charge & (PP_PME_CHARGE | PP_PME_C6)))
            {
                gmx_incons("PME-only node received coordinates before charges and/or C6-values"
                           );
            }

            /* The box, FE flag and lambda are sent along with the coordinates
             *  */
            copy_mat(cnb.box, box);
            *bFreeEnergy = (cnb.flags & PP_PME_FEP);
            *lambda_q    = cnb.lambda_q;
            *lambda_lj   = cnb.lambda_lj;
            *bEnerVir    = (cnb.flags & PP_PME_ENER_VIR);
            *pme_flags   = cnb.flags;

            if (*bFreeEnergy && !(pme_pp->flags_charge & PP_PME_CHARGEB))
            {
                gmx_incons("PME-only node received free energy request, but "
                           "did not receive B-state charges");
            }

            /* Receive the coordinates in place */
            nat = 0;
            for (sender = 0; sender < pme_pp->nnode; sender++)
            {
                if (pme_pp->nat[sender] > 0)
                {
                    MPI_Irecv(pme_pp->x[nat], pme_pp->nat[sender]*sizeof(rvec),
                              MPI_BYTE,
                              pme_pp->node[sender], 7,
                              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);
                    nat += pme_pp->nat[sender];
                    if (debug)
                    {
                        fprintf(debug, "Received from PP node %d: %d "
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

    *step = cnb.step;
#endif

    *natoms  = nat;
    *chargeA = pme_pp->chargeA;
    *chargeB = pme_pp->chargeB;
    *c6A     = pme_pp->c6A;
    *c6B     = pme_pp->c6B;
    *sigmaA  = pme_pp->sigmaA;
    *sigmaB  = pme_pp->sigmaB;
    *x       = pme_pp->x;
    *f       = pme_pp->f;

    return ((cnb.flags & PP_PME_FINISH) ? pmerecvqxFINISH : pmerecvqxX);
}

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
                    "PP node %d receiving from PME node %d: virial and energy\n",
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
    gmx_pme_send_q_x_wait(cr->dd);
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
                                 rvec *f, matrix vir_q, real energy_q,
                                 matrix vir_lj, real energy_lj,
                                 real dvdlambda_q, real dvdlambda_lj,
                                 float cycles)
{
    gmx_pme_comm_vir_ene_t cve;
    int                    messages, ind_start, ind_end, receiver;

    cve.cycles = cycles;

    /* Now the evaluated forces have to be transferred to the PP nodes */
    messages = 0;
    ind_end  = 0;
    for (receiver = 0; receiver < pme_pp->nnode; receiver++)
    {
        ind_start = ind_end;
        ind_end   = ind_start + pme_pp->nat[receiver];
#ifdef GMX_MPI
        if (MPI_Isend(f[ind_start], (ind_end-ind_start)*sizeof(rvec), MPI_BYTE,
                      pme_pp->node[receiver], 0,
                      pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]) != 0)
        {
            gmx_comm("MPI_Isend failed in do_pmeonly");
        }
#endif
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
        fprintf(debug, "PME node sending to PP node %d: virial and energy\n",
                pme_pp->node_peer);
    }
#ifdef GMX_MPI
    MPI_Isend(&cve, sizeof(cve), MPI_BYTE,
              pme_pp->node_peer, 1,
              pme_pp->mpi_comm_mysim, &pme_pp->req[messages++]);

    /* Wait for the forces to arrive */
    MPI_Waitall(messages, pme_pp->req, pme_pp->stat);
#endif
}
