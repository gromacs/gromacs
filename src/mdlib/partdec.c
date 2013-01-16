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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "invblock.h"
#include "macros.h"
#include "main.h"
#include "ns.h"
#include "partdec.h"
#include "splitter.h"
#include "gmx_random.h"
#include "mtop_util.h"
#include "mvdata.h"
#include "vec.h"

typedef struct gmx_partdec_constraint
{
    int                  left_range_receive;
    int                  right_range_receive;
    int                  left_range_send;
    int                  right_range_send;
    int                  nconstraints;
    int *                nlocalatoms;
    rvec *               sendbuf;
    rvec *               recvbuf;
}
gmx_partdec_constraint_t;


typedef struct gmx_partdec {
    int   neighbor[2];                         /* The nodeids of left and right neighb */
    int  *cgindex;                             /* The charge group boundaries,         */
                                               /* size nnodes+1,                       */
                                               /* only allocated with particle decomp. */
    int  *index;                               /* The home particle boundaries,        */
                                               /* size nnodes+1,                       */
                                               /* only allocated with particle decomp. */
    int  shift, bshift;                        /* Coordinates are shifted left for     */
                                               /* 'shift' systolic pulses, and right   */
                                               /* for 'bshift' pulses. Forces are      */
                                               /* shifted right for 'shift' pulses     */
                                               /* and left for 'bshift' pulses         */
                                               /* This way is not necessary to shift   */
                                               /* the coordinates over the entire ring */
    rvec                          *vbuf;       /* Buffer for summing the forces        */
#ifdef GMX_MPI
    MPI_Request                    mpi_req_rx; /* MPI reqs for async transfers */
    MPI_Request                    mpi_req_tx;
#endif
    gmx_partdec_constraint_t *     constraints;
} gmx_partdec_t;


void gmx_tx(const t_commrec *cr, int dir, void *buf, int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_tx");
#else
    int        nodeid;
    int        tag, flag;
    MPI_Status status;

    nodeid = cr->pd->neighbor[dir];

#ifdef DEBUG
    fprintf(stderr, "gmx_tx: nodeid=%d, buf=%x, bufsize=%d\n",
            nodeid, buf, bufsize);
#endif
#ifdef MPI_TEST
    /* workaround for crashes encountered with MPI on IRIX 6.5 */
    if (cr->pd->mpi_req_tx != MPI_REQUEST_NULL)
    {
        MPI_Test(&cr->pd->mpi_req_tx, &flag, &status);
        if (flag == FALSE)
        {
            fprintf(stdlog, "gmx_tx called before previous send was complete: nodeid=%d, buf=%x, bufsize=%d\n",
                    nodeid, buf, bufsize);
            gmx_tx_wait(nodeid);
        }
    }
#endif
    tag = 0;
    if (MPI_Isend(buf, bufsize, MPI_BYTE, RANK(cr, nodeid), tag, cr->mpi_comm_mygroup, &cr->pd->mpi_req_tx) != 0)
    {
        gmx_comm("MPI_Isend Failed");
    }
#endif
}

void gmx_tx_wait(const t_commrec *cr, int dir)
{
#ifndef GMX_MPI
    gmx_call("gmx_tx_wait");
#else
    MPI_Status  status;
    int         mpi_result;

    if ((mpi_result = MPI_Wait(&cr->pd->mpi_req_tx, &status)) != 0)
    {
        gmx_fatal(FARGS, "MPI_Wait: result=%d", mpi_result);
    }
#endif
}

void gmx_rx(const t_commrec *cr, int dir, void *buf, int bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_rx");
#else
    int nodeid;
    int tag;

    nodeid = cr->pd->neighbor[dir];
#ifdef DEBUG
    fprintf(stderr, "gmx_rx: nodeid=%d, buf=%x, bufsize=%d\n",
            nodeid, buf, bufsize);
#endif
    tag = 0;
    if (MPI_Irecv( buf, bufsize, MPI_BYTE, RANK(cr, nodeid), tag, cr->mpi_comm_mygroup, &cr->pd->mpi_req_rx) != 0)
    {
        gmx_comm("MPI_Recv Failed");
    }
#endif
}

void gmx_rx_wait(const t_commrec *cr, int nodeid)
{
#ifndef GMX_MPI
    gmx_call("gmx_rx_wait");
#else
    MPI_Status  status;
    int         mpi_result;

    if ((mpi_result = MPI_Wait(&cr->pd->mpi_req_rx, &status)) != 0)
    {
        gmx_fatal(FARGS, "MPI_Wait: result=%d", mpi_result);
    }
#endif
}

void gmx_tx_rx_real(const t_commrec *cr,
                    int send_dir, real *send_buf, int send_bufsize,
                    int recv_dir, real *recv_buf, int recv_bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_tx_rx_real");
#else
    int        send_nodeid, recv_nodeid;
    int        tx_tag = 0, rx_tag = 0;
    MPI_Status stat;

    send_nodeid = cr->pd->neighbor[send_dir];
    recv_nodeid = cr->pd->neighbor[recv_dir];

#ifdef GMX_DOUBLE
#define mpi_type MPI_DOUBLE
#else
#define mpi_type MPI_FLOAT
#endif

    if (send_bufsize > 0 && recv_bufsize > 0)
    {
        MPI_Sendrecv(send_buf, send_bufsize, mpi_type, RANK(cr, send_nodeid), tx_tag,
                     recv_buf, recv_bufsize, mpi_type, RANK(cr, recv_nodeid), rx_tag,
                     cr->mpi_comm_mygroup, &stat);
    }
    else if (send_bufsize > 0)
    {
        MPI_Send(send_buf, send_bufsize, mpi_type, RANK(cr, send_nodeid), tx_tag,
                 cr->mpi_comm_mygroup);
    }
    else if (recv_bufsize > 0)
    {
        MPI_Recv(recv_buf, recv_bufsize, mpi_type, RANK(cr, recv_nodeid), rx_tag,
                 cr->mpi_comm_mygroup, &stat);
    }
#undef mpi_type
#endif
}


void gmx_tx_rx_void(const t_commrec *cr,
                    int send_dir, void *send_buf, int send_bufsize,
                    int recv_dir, void *recv_buf, int recv_bufsize)
{
#ifndef GMX_MPI
    gmx_call("gmx_tx_rx_void");
#else
    int        send_nodeid, recv_nodeid;
    int        tx_tag = 0, rx_tag = 0;
    MPI_Status stat;

    send_nodeid = cr->pd->neighbor[send_dir];
    recv_nodeid = cr->pd->neighbor[recv_dir];


    MPI_Sendrecv(send_buf, send_bufsize, MPI_BYTE, RANK(cr, send_nodeid), tx_tag,
                 recv_buf, recv_bufsize, MPI_BYTE, RANK(cr, recv_nodeid), rx_tag,
                 cr->mpi_comm_mygroup, &stat);

#endif
}


/*void gmx_wait(int dir_send,int dir_recv)*/

void gmx_wait(const t_commrec *cr, int dir_send, int dir_recv)
{
#ifndef GMX_MPI
    gmx_call("gmx_wait");
#else
    gmx_tx_wait(cr, dir_send);
    gmx_rx_wait(cr, dir_recv);
#endif
}

static void set_left_right(t_commrec *cr)
{
    cr->pd->neighbor[GMX_LEFT]  = (cr->nnodes + cr->nodeid - 1) % cr->nnodes;
    cr->pd->neighbor[GMX_RIGHT] = (cr->nodeid + 1) % cr->nnodes;
}

void pd_move_f(const t_commrec *cr, rvec f[], t_nrnb *nrnb)
{
    move_f(NULL, cr, GMX_LEFT, GMX_RIGHT, f, cr->pd->vbuf, nrnb);
}

int *pd_cgindex(const t_commrec *cr)
{
    return cr->pd->cgindex;
}

int *pd_index(const t_commrec *cr)
{
    return cr->pd->index;
}

int pd_shift(const t_commrec *cr)
{
    return cr->pd->shift;
}

int pd_bshift(const t_commrec *cr)
{
    return cr->pd->bshift;
}

void pd_cg_range(const t_commrec *cr, int *cg0, int *cg1)
{
    *cg0 = cr->pd->cgindex[cr->nodeid];
    *cg1 = cr->pd->cgindex[cr->nodeid+1];
}

void pd_at_range(const t_commrec *cr, int *at0, int *at1)
{
    *at0 = cr->pd->index[cr->nodeid];
    *at1 = cr->pd->index[cr->nodeid+1];
}

void
pd_get_constraint_range(gmx_partdec_p_t pd, int *start, int *natoms)
{
    *start  = pd->constraints->left_range_receive;
    *natoms = pd->constraints->right_range_receive-pd->constraints->left_range_receive;
}

int *
pd_constraints_nlocalatoms(gmx_partdec_p_t pd)
{
    int *rc;

    if (NULL != pd && NULL != pd->constraints)
    {
        rc = pd->constraints->nlocalatoms;
    }
    else
    {
        rc = NULL;
    }
    return rc;
}




/* This routine is used to communicate the non-home-atoms needed for constrains.
 * We have already calculated this range of atoms during setup, and stored in the
 * partdec constraints structure.
 *
 * When called, we send/receive left_range_send/receive atoms to our left (lower)
 * node neighbor, and similar to the right (higher) node.
 *
 * This has not been tested for periodic molecules...
 */
void
pd_move_x_constraints(t_commrec *  cr,
                      rvec *       x0,
                      rvec *       x1)
{
#ifdef GMX_MPI
    gmx_partdec_t            *pd;
    gmx_partdec_constraint_t *pdc;

    rvec                *     sendptr;
    rvec                *     recvptr;
    int                       thisnode;
    int                       i;
    int                       cnt;
    int                       sendcnt, recvcnt;

    pd  = cr->pd;
    pdc = pd->constraints;

    if (pdc == NULL)
    {
        return;
    }

    thisnode  = cr->nodeid;

    /* First pulse to right */

    recvcnt = 3*(pd->index[thisnode]-pdc->left_range_receive);
    sendcnt = 3*(cr->pd->index[thisnode+1]-cr->pd->constraints->right_range_send);

    if (x1 != NULL)
    {
        /* Assemble temporary array with both x0 & x1 */
        recvptr = pdc->recvbuf;
        sendptr = pdc->sendbuf;

        cnt = 0;
        for (i = pdc->right_range_send; i < pd->index[thisnode+1]; i++)
        {
            copy_rvec(x0[i], sendptr[cnt++]);
        }
        for (i = pdc->right_range_send; i < pd->index[thisnode+1]; i++)
        {
            copy_rvec(x1[i], sendptr[cnt++]);
        }
        recvcnt *= 2;
        sendcnt *= 2;
    }
    else
    {
        recvptr = x0 + pdc->left_range_receive;
        sendptr = x0 + pdc->right_range_send;
    }

    gmx_tx_rx_real(cr,
                   GMX_RIGHT, (real *)sendptr, sendcnt,
                   GMX_LEFT, (real *)recvptr, recvcnt);

    if (x1 != NULL)
    {
        /* copy back to x0/x1 */
        cnt = 0;
        for (i = pdc->left_range_receive; i < pd->index[thisnode]; i++)
        {
            copy_rvec(recvptr[cnt++], x0[i]);
        }
        for (i = pdc->left_range_receive; i < pd->index[thisnode]; i++)
        {
            copy_rvec(recvptr[cnt++], x1[i]);
        }
    }

    /* And pulse to left */
    sendcnt = 3*(pdc->left_range_send-pd->index[thisnode]);
    recvcnt = 3*(pdc->right_range_receive-pd->index[thisnode+1]);

    if (x1 != NULL)
    {
        cnt = 0;
        for (i = cr->pd->index[thisnode]; i < pdc->left_range_send; i++)
        {
            copy_rvec(x0[i], sendptr[cnt++]);
        }
        for (i = cr->pd->index[thisnode]; i < pdc->left_range_send; i++)
        {
            copy_rvec(x1[i], sendptr[cnt++]);
        }
        recvcnt *= 2;
        sendcnt *= 2;
    }
    else
    {
        sendptr = x0 + pd->index[thisnode];
        recvptr = x0 + pd->index[thisnode+1];
    }

    gmx_tx_rx_real(cr,
                   GMX_LEFT, (real *)sendptr, sendcnt,
                   GMX_RIGHT, (real *)recvptr, recvcnt);

    /* Final copy back from buffers */
    if (x1 != NULL)
    {
        /* First copy received data back into x0 & x1 */
        cnt = 0;
        for (i = pd->index[thisnode+1]; i < pdc->right_range_receive; i++)
        {
            copy_rvec(recvptr[cnt++], x0[i]);
        }
        for (i = pd->index[thisnode+1]; i < pdc->right_range_receive; i++)
        {
            copy_rvec(recvptr[cnt++], x1[i]);
        }
    }
#endif
}

static int home_cpu(int nnodes, gmx_partdec_t *pd, int atomid)
{
    int k;

    for (k = 0; (k < nnodes); k++)
    {
        if (atomid < pd->index[k+1])
        {
            return k;
        }
    }
    gmx_fatal(FARGS, "Atomid %d is larger than number of atoms (%d)",
              atomid+1, pd->index[nnodes]+1);

    return -1;
}

static void calc_nsbshift(FILE *fp, int nnodes, gmx_partdec_t *pd, t_idef *idef)
{
    int i, j, k;
    int lastcg, targetcg, nshift, naaj;
    int homeid[32];

    pd->bshift = 0;
    for (i = 1; (i < nnodes); i++)
    {
        targetcg = pd->cgindex[i];
        for (nshift = i; (nshift > 0) && (pd->cgindex[nshift] > targetcg); nshift--)
        {
            ;
        }
        pd->bshift = max(pd->bshift, i-nshift);
    }

    pd->shift = (nnodes + 1)/2;
    for (i = 0; (i < nnodes); i++)
    {
        lastcg   = pd->cgindex[i+1]-1;
        naaj     = calc_naaj(lastcg, pd->cgindex[nnodes]);
        targetcg = (lastcg+naaj) % pd->cgindex[nnodes];

        /* Search until we find the target charge group */
        for (nshift = 0;
             (nshift < nnodes) && (targetcg > pd->cgindex[nshift+1]);
             nshift++)
        {
            ;
        }
        /* Now compute the shift, that is the difference in node index */
        nshift = ((nshift - i + nnodes) % nnodes);

        if (fp)
        {
            fprintf(fp, "CPU=%3d, lastcg=%5d, targetcg=%5d, myshift=%5d\n",
                    i, lastcg, targetcg, nshift);
        }

        /* It's the largest shift that matters */
        pd->shift = max(nshift, pd->shift);
    }
    /* Now for the bonded forces */
    for (i = 0; (i < F_NRE); i++)
    {
        if (interaction_function[i].flags & IF_BOND)
        {
            int nratoms = interaction_function[i].nratoms;
            for (j = 0; (j < idef->il[i].nr); j += nratoms+1)
            {
                for (k = 1; (k <= nratoms); k++)
                {
                    homeid[k-1] = home_cpu(nnodes, pd, idef->il[i].iatoms[j+k]);
                }
                for (k = 1; (k < nratoms); k++)
                {
                    pd->shift = max(pd->shift, abs(homeid[k]-homeid[0]));
                }
            }
        }
    }
    if (fp)
    {
        fprintf(fp, "pd->shift = %3d, pd->bshift=%3d\n",
                pd->shift, pd->bshift);
    }
}


static void
init_partdec_constraint(t_commrec *cr,
                        t_idef *   idef,
                        int       *left_range,
                        int       *right_range)
{
    gmx_partdec_t *            pd = cr->pd;
    gmx_partdec_constraint_t  *pdc;
    int i, cnt, k;
    int ai, aj, nodei, nodej;
    int nratoms;
    int nodeid;

    snew(pdc, 1);
    cr->pd->constraints = pdc;


    /* Who am I? */
    nodeid = cr->nodeid;

    /* Setup LINCS communication ranges */
    pdc->left_range_receive   = left_range[nodeid];
    pdc->right_range_receive  = right_range[nodeid]+1;
    pdc->left_range_send      = (nodeid > 0) ? right_range[nodeid-1]+1 : 0;
    pdc->right_range_send     = (nodeid < cr->nnodes-1) ? left_range[nodeid+1] : pd->index[cr->nnodes];

    snew(pdc->nlocalatoms, idef->il[F_CONSTR].nr);
    nratoms = interaction_function[F_CONSTR].nratoms;

    for (i = 0, cnt = 0; i < idef->il[F_CONSTR].nr; i += nratoms+1, cnt++)
    {
        ai    = idef->il[F_CONSTR].iatoms[i+1];
        aj    = idef->il[F_CONSTR].iatoms[i+2];
        nodei = 0;
        while (ai >= pd->index[nodei+1])
        {
            nodei++;
        }
        nodej = 0;
        while (aj >= pd->index[nodej+1])
        {
            nodej++;
        }
        pdc->nlocalatoms[cnt] = 0;
        if (nodei == nodeid)
        {
            pdc->nlocalatoms[cnt]++;
        }
        if (nodej == nodeid)
        {
            pdc->nlocalatoms[cnt]++;
        }
    }
    pdc->nconstraints = cnt;

    snew(pdc->sendbuf, max(6*(pd->index[cr->nodeid+1]-pd->constraints->right_range_send), 6*(pdc->left_range_send-pd->index[cr->nodeid])));
    snew(pdc->recvbuf, max(6*(pd->index[cr->nodeid]-pdc->left_range_receive), 6*(pdc->right_range_receive-pd->index[cr->nodeid+1])));

}

static void init_partdec(FILE *fp, t_commrec *cr, t_block *cgs, int *multinr,
                         t_idef *idef)
{
    int            i, nodeid;
    gmx_partdec_t *pd;

    snew(pd, 1);
    cr->pd = pd;

    set_left_right(cr);

    if (cr->nnodes > 1)
    {
        if (multinr == NULL)
        {
            gmx_fatal(FARGS, "Internal error in init_partdec: multinr = NULL");
        }
        snew(pd->index, cr->nnodes+1);
        snew(pd->cgindex, cr->nnodes+1);
        pd->cgindex[0] = 0;
        pd->index[0]   = 0;
        for (i = 0; (i < cr->nnodes); i++)
        {
            pd->cgindex[i+1] = multinr[i];
            pd->index[i+1]   = cgs->index[multinr[i]];
        }
        calc_nsbshift(fp, cr->nnodes, pd, idef);
        /* This is a hack to do with bugzilla 148 */
        /*pd->shift = cr->nnodes-1;
           pd->bshift = 0;*/

        /* Allocate a buffer of size natoms of the whole system
         * for summing the forces over the nodes.
         */
        snew(pd->vbuf, cgs->index[cgs->nr]);
        pd->constraints = NULL;
    }
#ifdef GMX_MPI
    pd->mpi_req_tx = MPI_REQUEST_NULL;
    pd->mpi_req_rx = MPI_REQUEST_NULL;
#endif
}

static void print_partdec(FILE *fp, const char *title,
                          int nnodes, gmx_partdec_t *pd)
{
    int i;

    fprintf(fp, "%s\n", title);
    fprintf(fp, "nnodes:       %5d\n", nnodes);
    fprintf(fp, "pd->shift:    %5d\n", pd->shift);
    fprintf(fp, "pd->bshift:   %5d\n", pd->bshift);

    fprintf(fp, "Nodeid   atom0   #atom     cg0       #cg\n");
    for (i = 0; (i < nnodes); i++)
    {
        fprintf(fp, "%6d%8d%8d%8d%10d\n",
                i,
                pd->index[i], pd->index[i+1]-pd->index[i],
                pd->cgindex[i], pd->cgindex[i+1]-pd->cgindex[i]);
    }
    fprintf(fp, "\n");
}

static void pr_idef_division(FILE *fp, t_idef *idef, int nnodes, int **multinr)
{
    int i, ftype, nr, nra, m0, m1;

    fprintf(fp, "Division of bonded forces over processors\n");
    fprintf(fp, "%-12s", "CPU");
    for (i = 0; (i < nnodes); i++)
    {
        fprintf(fp, " %5d", i);
    }
    fprintf(fp, "\n");

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (idef->il[ftype].nr > 0)
        {
            nr  = idef->il[ftype].nr;
            nra = 1+interaction_function[ftype].nratoms;
            fprintf(fp, "%-12s", interaction_function[ftype].name);
            /* Loop over processors */
            for (i = 0; (i < nnodes); i++)
            {
                m0 = (i == 0) ? 0 : multinr[ftype][i-1]/nra;
                m1 = multinr[ftype][i]/nra;
                fprintf(fp, " %5d", m1-m0);
            }
            fprintf(fp, "\n");
        }
    }
}

static void select_my_ilist(FILE *log, t_ilist *il, int *multinr, t_commrec *cr)
{
    t_iatom *ia;
    int      i, start, end, nr;

    if (cr->nodeid == 0)
    {
        start = 0;
    }
    else
    {
        start = multinr[cr->nodeid-1];
    }
    end = multinr[cr->nodeid];

    nr = end-start;
    if (nr < 0)
    {
        gmx_fatal(FARGS, "Negative number of atoms (%d) on node %d\n"
                  "You have probably not used the same value for -np with grompp"
                  " and mdrun",
                  nr, cr->nodeid);
    }
    snew(ia, nr);

    for (i = 0; (i < nr); i++)
    {
        ia[i] = il->iatoms[start+i];
    }

    sfree(il->iatoms);
    il->iatoms = ia;

    il->nr = nr;
}

static void select_my_idef(FILE *log, t_idef *idef, int **multinr, t_commrec *cr)
{
    int i;

    for (i = 0; (i < F_NRE); i++)
    {
        select_my_ilist(log, &idef->il[i], multinr[i], cr);
    }
}

gmx_localtop_t *split_system(FILE *log,
                             gmx_mtop_t *mtop, t_inputrec *inputrec,
                             t_commrec *cr)
{
    int             i, npp, n, nn;
    real           *capacity;
    double          tcap = 0, cap;
    int            *multinr_cgs, **multinr_nre;
    char           *cap_env;
    gmx_localtop_t *top;
    int            *left_range;
    int            *right_range;

    /* Time to setup the division of charge groups over processors */
    npp = cr->nnodes-cr->npmenodes;
    snew(capacity, npp);
    cap_env = getenv("GMX_CAPACITY");
    if (cap_env == NULL)
    {
        for (i = 0; (i < npp-1); i++)
        {
            capacity[i] = 1.0/(double)npp;
            tcap       += capacity[i];
        }
        /* Take care that the sum of capacities is 1.0 */
        capacity[npp-1] = 1.0 - tcap;
    }
    else
    {
        tcap = 0;
        nn   = 0;
        for (i = 0; i < npp; i++)
        {
            cap = 0;
            sscanf(cap_env+nn, "%lf%n", &cap, &n);
            if (cap == 0)
            {
                gmx_fatal(FARGS, "Incorrect entry or number of entries in GMX_CAPACITY='%s'", cap_env);
            }
            capacity[i] = cap;
            tcap       += cap;
            nn         += n;
        }
        for (i = 0; i < npp; i++)
        {
            capacity[i] /= tcap;
        }
    }

    /* Convert the molecular topology to a fully listed topology */
    top = gmx_mtop_generate_local_top(mtop, inputrec);

    snew(multinr_cgs, npp);
    snew(multinr_nre, F_NRE);
    for (i = 0; i < F_NRE; i++)
    {
        snew(multinr_nre[i], npp);
    }


    snew(left_range, cr->nnodes);
    snew(right_range, cr->nnodes);

    /* This computes which entities can be placed on processors */
    split_top(log, npp, top, inputrec, &mtop->mols, capacity, multinr_cgs, multinr_nre, left_range, right_range);

    sfree(capacity);
    init_partdec(log, cr, &(top->cgs), multinr_cgs, &(top->idef));

    /* This should be fine */
    /*split_idef(&(top->idef),cr->nnodes-cr->npmenodes);*/

    select_my_idef(log, &(top->idef), multinr_nre, cr);

    if (top->idef.il[F_CONSTR].nr > 0)
    {
        init_partdec_constraint(cr, &(top->idef), left_range, right_range);
    }

    if (log)
    {
        pr_idef_division(log, &(top->idef), npp, multinr_nre);
    }

    for (i = 0; i < F_NRE; i++)
    {
        sfree(multinr_nre[i]);
    }
    sfree(multinr_nre);
    sfree(multinr_cgs);

    sfree(left_range);
    sfree(right_range);

    if (log)
    {
        print_partdec(log, "Workload division", cr->nnodes, cr->pd);
    }

    return top;
}

static void
add_to_vsitelist(int **list, int *nitem, int *nalloc, int newitem)
{
    int      i, idx;
    gmx_bool found;

    found = FALSE;
    idx   = *nitem;
    for (i = 0; i < idx && !found; i++)
    {
        found = (newitem == (*list)[i]);
    }
    if (!found)
    {
        *nalloc += 100;
        srenew(*list, *nalloc);
        (*list)[idx++] = newitem;
        *nitem         = idx;
    }
}

gmx_bool setup_parallel_vsites(t_idef *idef, t_commrec *cr,
                               t_comm_vsites *vsitecomm)
{
    int            i, j, ftype;
    int            nra;
    gmx_bool       do_comm;
    t_iatom       *ia;
    gmx_partdec_t *pd;
    int            iconstruct;
    int            i0, i1;
    int            nalloc_left_construct, nalloc_right_construct;
    int            sendbuf[2], recvbuf[2];
    int            bufsize, leftbuf, rightbuf;

    pd = cr->pd;

    i0 = pd->index[cr->nodeid];
    i1 = pd->index[cr->nodeid+1];

    vsitecomm->left_import_construct  = NULL;
    vsitecomm->left_import_nconstruct = 0;
    nalloc_left_construct             = 0;

    vsitecomm->right_import_construct  = NULL;
    vsitecomm->right_import_nconstruct = 0;
    nalloc_right_construct             = 0;

    for (ftype = 0; (ftype < F_NRE); ftype++)
    {
        if (!(interaction_function[ftype].flags & IF_VSITE) )
        {
            continue;
        }

        nra    = interaction_function[ftype].nratoms;
        ia     = idef->il[ftype].iatoms;

        for (i = 0; i < idef->il[ftype].nr; i += nra+1)
        {
            for (j = 2; j < 1+nra; j++)
            {
                iconstruct = ia[i+j];
                if (iconstruct < i0)
                {
                    add_to_vsitelist(&vsitecomm->left_import_construct,
                                     &vsitecomm->left_import_nconstruct,
                                     &nalloc_left_construct, iconstruct);
                }
                else if (iconstruct >= i1)
                {
                    add_to_vsitelist(&vsitecomm->right_import_construct,
                                     &vsitecomm->right_import_nconstruct,
                                     &nalloc_right_construct, iconstruct);
                }
            }
        }
    }

    /* Pre-communicate the array lengths */
    gmx_tx_rx_void(cr,
                   GMX_RIGHT, (void *)&vsitecomm->right_import_nconstruct, sizeof(int),
                   GMX_LEFT, (void *)&vsitecomm->left_export_nconstruct, sizeof(int));
    gmx_tx_rx_void(cr,
                   GMX_LEFT, (void *)&vsitecomm->left_import_nconstruct, sizeof(int),
                   GMX_RIGHT, (void *)&vsitecomm->right_export_nconstruct, sizeof(int));

    snew(vsitecomm->left_export_construct, vsitecomm->left_export_nconstruct);
    snew(vsitecomm->right_export_construct, vsitecomm->right_export_nconstruct);

    /* Communicate the construcing atom index arrays */
    gmx_tx_rx_void(cr,
                   GMX_RIGHT, (void *)vsitecomm->right_import_construct, vsitecomm->right_import_nconstruct*sizeof(int),
                   GMX_LEFT, (void *)vsitecomm->left_export_construct, vsitecomm->left_export_nconstruct*sizeof(int));

    /* Communicate the construcing atom index arrays */
    gmx_tx_rx_void(cr,
                   GMX_LEFT, (void *)vsitecomm->left_import_construct, vsitecomm->left_import_nconstruct*sizeof(int),
                   GMX_RIGHT, (void *)vsitecomm->right_export_construct, vsitecomm->right_export_nconstruct*sizeof(int));

    leftbuf  = max(vsitecomm->left_export_nconstruct, vsitecomm->left_import_nconstruct);
    rightbuf = max(vsitecomm->right_export_nconstruct, vsitecomm->right_import_nconstruct);

    bufsize  = max(leftbuf, rightbuf);

    do_comm = (bufsize > 0);

    snew(vsitecomm->send_buf, 2*bufsize);
    snew(vsitecomm->recv_buf, 2*bufsize);

    return do_comm;
}

t_state *partdec_init_local_state(t_commrec *cr, t_state *state_global)
{
    int      i;
    t_state *state_local;

    snew(state_local, 1);

    /* Copy all the contents */
    *state_local = *state_global;
    snew(state_local->lambda, efptNR);
    /* local storage for lambda */
    for (i = 0; i < efptNR; i++)
    {
        state_local->lambda[i] = state_global->lambda[i];
    }
    if (state_global->nrngi > 1)
    {
        /* With stochastic dynamics we need local storage for the random state */
        if (state_local->flags & (1<<estLD_RNG))
        {
            state_local->nrng = gmx_rng_n();
            snew(state_local->ld_rng, state_local->nrng);
#ifdef GMX_MPI
            if (PAR(cr))
            {
                MPI_Scatter(state_global->ld_rng,
                            state_local->nrng*sizeof(state_local->ld_rng[0]), MPI_BYTE,
                            state_local->ld_rng,
                            state_local->nrng*sizeof(state_local->ld_rng[0]), MPI_BYTE,
                            MASTERRANK(cr), cr->mpi_comm_mygroup);
            }
#endif
        }
        if (state_local->flags & (1<<estLD_RNGI))
        {
            snew(state_local->ld_rngi, 1);
#ifdef GMX_MPI
            if (PAR(cr))
            {
                MPI_Scatter(state_global->ld_rngi,
                            sizeof(state_local->ld_rngi[0]), MPI_BYTE,
                            state_local->ld_rngi,
                            sizeof(state_local->ld_rngi[0]), MPI_BYTE,
                            MASTERRANK(cr), cr->mpi_comm_mygroup);
            }
#endif
        }
    }

    return state_local;
}
