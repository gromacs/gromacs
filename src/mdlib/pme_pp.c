/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#include "mpelogging.h"

#define PP_PME_CHARGE         (1<<0)
#define PP_PME_CHARGEB        (1<<1)
#define PP_PME_COORD          (1<<2)
#define PP_PME_FEP            (1<<3)
#define PP_PME_ENER_VIR       (1<<4)
#define PP_PME_FINISH         (1<<5)
#define PP_PME_SWITCHGRID     (1<<6)
#define PP_PME_RESETCOUNTERS  (1<<7)


#define PME_PP_SIGSTOP     (1<<0)
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
    real            lambda;
    int             flags;
    gmx_large_int_t step;
    ivec            grid_size;  /* For PME grid tuning */
    real            ewaldcoeff; /* For PME grid tuning */
} gmx_pme_comm_n_box_t;

typedef struct {
    matrix          vir;
    real            energy;
    real            dvdlambda;
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
    snew(pme_pp->req, 2*pme_pp->nnode);
    snew(pme_pp->stat, 2*pme_pp->nnode);
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
                             matrix box, rvec *x,
                             real lambda,
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
        cnb->lambda     = lambda;
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
    else if (flags & PP_PME_CHARGE)
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
        if (flags & PP_PME_COORD)
        {
            MPI_Isend(x[0], n*sizeof(rvec), MPI_BYTE,
                      dd->pme_nodeid, 3, cr->mpi_comm_mysim,
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
                    int maxshift_x, int maxshift_y)
{
    int flags;

    flags = PP_PME_CHARGE;
    if (bFreeEnergy)
    {
        flags |= PP_PME_CHARGEB;
    }

    gmx_pme_send_q_x(cr, flags,
                     chargeA, chargeB, NULL, NULL, 0, maxshift_x, maxshift_y, -1);
}

void gmx_pme_send_x(t_commrec *cr, matrix box, rvec *x,
                    gmx_bool bFreeEnergy, real lambda,
                    gmx_bool bEnerVir,
                    gmx_large_int_t step)
{
    int flags;

    flags = PP_PME_COORD;
    if (bFreeEnergy)
    {
        flags |= PP_PME_FEP;
    }
    if (bEnerVir)
    {
        flags |= PP_PME_ENER_VIR;
    }

    gmx_pme_send_q_x(cr, flags, NULL, NULL, box, x, lambda, 0, 0, step);
}

void gmx_pme_send_finish(t_commrec *cr)
{
    int flags;

    flags = PP_PME_FINISH;

    gmx_pme_send_q_x(cr, flags, NULL, NULL, NULL, NULL, 0, 0, 0, -1);
}

void gmx_pme_send_switchgrid(t_commrec *cr, ivec grid_size, real ewaldcoeff)
{
#ifdef GMX_MPI
    gmx_pme_comm_n_box_t cnb;

    /* Only let one PP node signal each PME node */
    if (cr->dd->pme_receive_vir_ener)
    {
        cnb.flags = PP_PME_SWITCHGRID;
        copy_ivec(grid_size, cnb.grid_size);
        cnb.ewaldcoeff = ewaldcoeff;

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
                     matrix box, rvec **x, rvec **f,
                     int *maxshift_x, int *maxshift_y,
                     gmx_bool *bFreeEnergy, real *lambda,
                     gmx_bool *bEnerVir,
                     gmx_large_int_t *step,
                     ivec grid_size, real *ewaldcoeff)
{
    gmx_pme_comm_n_box_t cnb;
    int                  nat = 0, q, messages, sender;
    real                *charge_pp;

    messages = 0;

    /* avoid compiler warning about unused variable without MPI support */
    cnb.flags = 0;
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
            *ewaldcoeff = cnb.ewaldcoeff;

            return pmerecvqxSWITCHGRID;
        }

        if (cnb.flags & PP_PME_RESETCOUNTERS)
        {
            /* Special case, receive the step and return */
            *step = cnb.step;

            return pmerecvqxRESETCOUNTERS;
        }

        if (cnb.flags & PP_PME_CHARGE)
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
                srenew(pme_pp->chargeA, pme_pp->nalloc);
                if (cnb.flags & PP_PME_CHARGEB)
                {
                    srenew(pme_pp->chargeB, pme_pp->nalloc);
                }
                srenew(pme_pp->x, pme_pp->nalloc);
                srenew(pme_pp->f, pme_pp->nalloc);
            }

            /* maxshift is sent when the charges are sent */
            *maxshift_x = cnb.maxshift_x;
            *maxshift_y = cnb.maxshift_y;

            /* Receive the charges in place */
            for (q = 0; q < ((cnb.flags & PP_PME_CHARGEB) ? 2 : 1); q++)
            {
                if (q == 0)
                {
                    charge_pp = pme_pp->chargeA;
                }
                else
                {
                    charge_pp = pme_pp->chargeB;
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
            if (!(pme_pp->flags_charge & PP_PME_CHARGE))
            {
                gmx_incons("PME-only node received coordinates before charges"
                           );
            }

            /* The box, FE flag and lambda are sent along with the coordinates
             *  */
            copy_mat(cnb.box, box);
            *bFreeEnergy = (cnb.flags & PP_PME_FEP);
            *lambda      = cnb.lambda;
            *bEnerVir    = (cnb.flags & PP_PME_ENER_VIR);

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
                              pme_pp->node[sender], 3,
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
    *x       = pme_pp->x;
    *f       = pme_pp->f;

    return ((cnb.flags & PP_PME_FINISH) ? pmerecvqxFINISH : pmerecvqxX);
}

static void receive_virial_energy(t_commrec *cr,
                                  matrix vir, real *energy, real *dvdlambda,
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

        m_add(vir, cve.vir, vir);
        *energy     = cve.energy;
        *dvdlambda += cve.dvdlambda;
        *pme_cycles = cve.cycles;

        if (cve.stop_cond != gmx_stop_cond_none)
        {
            gmx_set_stop_condition(cve.stop_cond);
        }
    }
    else
    {
        *energy     = 0;
        *pme_cycles = 0;
    }
}

void gmx_pme_receive_f(t_commrec *cr,
                       rvec f[], matrix vir,
                       real *energy, real *dvdlambda,
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


    receive_virial_energy(cr, vir, energy, dvdlambda, pme_cycles);
}

void gmx_pme_send_force_vir_ener(struct gmx_pme_pp *pme_pp,
                                 rvec *f, matrix vir,
                                 real energy, real dvdlambda,
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
    copy_mat(vir, cve.vir);
    cve.energy    = energy;
    cve.dvdlambda = dvdlambda;
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
