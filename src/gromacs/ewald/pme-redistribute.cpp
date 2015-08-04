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
#include "gmxpre.h"

#include "pme-redistribute.h"

#include "config.h"

#include <algorithm>

#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-simd.h"

static void pme_calc_pidx(int start, int end,
                          matrix recipbox, rvec x[],
                          pme_atomcomm_t *atc, int *count)
{
    int   nslab, i;
    int   si;
    real *xptr, s;
    real  rxx, ryx, rzx, ryy, rzy;
    int  *pd;

    /* Calculate PME task index (pidx) for each grid index.
     * Here we always assign equally sized slabs to each node
     * for load balancing reasons (the PME grid spacing is not used).
     */

    nslab = atc->nslab;
    pd    = atc->pd;

    /* Reset the count */
    for (i = 0; i < nslab; i++)
    {
        count[i] = 0;
    }

    if (atc->dimind == 0)
    {
        rxx = recipbox[XX][XX];
        ryx = recipbox[YY][XX];
        rzx = recipbox[ZZ][XX];
        /* Calculate the node index in x-dimension */
        for (i = start; i < end; i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab*(xptr[XX]*rxx + xptr[YY]*ryx + xptr[ZZ]*rzx);
            si    = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
    else
    {
        ryy = recipbox[YY][YY];
        rzy = recipbox[ZZ][YY];
        /* Calculate the node index in y-dimension */
        for (i = start; i < end; i++)
        {
            xptr   = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab*(xptr[YY]*ryy + xptr[ZZ]*rzy);
            si    = (int)(s + 2*nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
}

static void pme_calc_pidx_wrapper(int natoms, matrix recipbox, rvec x[],
                                  pme_atomcomm_t *atc)
{
    int nthread, thread, slab;

    nthread = atc->nthread;

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (thread = 0; thread < nthread; thread++)
    {
        pme_calc_pidx(natoms* thread   /nthread,
                      natoms*(thread+1)/nthread,
                      recipbox, x, atc, atc->count_thread[thread]);
    }
    /* Non-parallel reduction, since nslab is small */

    for (thread = 1; thread < nthread; thread++)
    {
        for (slab = 0; slab < atc->nslab; slab++)
        {
            atc->count_thread[0][slab] += atc->count_thread[thread][slab];
        }
    }
}

static void realloc_splinevec(splinevec th, real **ptr_z, int nalloc)
{
    const int padding = 4;
    int       i;

    srenew(th[XX], nalloc);
    srenew(th[YY], nalloc);
    /* In z we add padding, this is only required for the aligned SIMD code */
    sfree_aligned(*ptr_z);
    snew_aligned(*ptr_z, nalloc+2*padding, SIMD4_ALIGNMENT);
    th[ZZ] = *ptr_z + padding;

    for (i = 0; i < padding; i++)
    {
        (*ptr_z)[               i] = 0;
        (*ptr_z)[padding+nalloc+i] = 0;
    }
}

static void pme_realloc_splinedata(splinedata_t *spline, pme_atomcomm_t *atc)
{
    int i;

    srenew(spline->ind, atc->nalloc);
    /* Initialize the index to identity so it works without threads */
    for (i = 0; i < atc->nalloc; i++)
    {
        spline->ind[i] = i;
    }

    realloc_splinevec(spline->theta, &spline->ptr_theta_z,
                      atc->pme_order*atc->nalloc);
    realloc_splinevec(spline->dtheta, &spline->ptr_dtheta_z,
                      atc->pme_order*atc->nalloc);
}

void pme_realloc_atomcomm_things(pme_atomcomm_t *atc)
{
    int nalloc_old, i;

    /* We have to avoid a NULL pointer for atc->x to avoid
     * possible fatal errors in MPI routines.
     */
    if (atc->n > atc->nalloc || atc->nalloc == 0)
    {
        nalloc_old  = atc->nalloc;
        atc->nalloc = over_alloc_dd(std::max(atc->n, 1));

        if (atc->nslab > 1)
        {
            srenew(atc->x, atc->nalloc);
            srenew(atc->coefficient, atc->nalloc);
            srenew(atc->f, atc->nalloc);
            for (i = nalloc_old; i < atc->nalloc; i++)
            {
                clear_rvec(atc->f[i]);
            }
        }
        if (atc->bSpread)
        {
            srenew(atc->fractx, atc->nalloc);
            srenew(atc->idx, atc->nalloc);

            if (atc->nthread > 1)
            {
                srenew(atc->thread_idx, atc->nalloc);
            }

            for (i = 0; i < atc->nthread; i++)
            {
                pme_realloc_splinedata(&atc->spline[i], atc);
            }
        }
    }
}

static void pme_dd_sendrecv(pme_atomcomm_t gmx_unused *atc,
                            gmx_bool gmx_unused bBackward, int gmx_unused shift,
                            void gmx_unused *buf_s, int gmx_unused nbyte_s,
                            void gmx_unused *buf_r, int gmx_unused nbyte_r)
{
#ifdef GMX_MPI
    int        dest, src;
    MPI_Status stat;

    if (bBackward == FALSE)
    {
        dest = atc->node_dest[shift];
        src  = atc->node_src[shift];
    }
    else
    {
        dest = atc->node_src[shift];
        src  = atc->node_dest[shift];
    }

    if (nbyte_s > 0 && nbyte_r > 0)
    {
        MPI_Sendrecv(buf_s, nbyte_s, MPI_BYTE,
                     dest, shift,
                     buf_r, nbyte_r, MPI_BYTE,
                     src, shift,
                     atc->mpi_comm, &stat);
    }
    else if (nbyte_s > 0)
    {
        MPI_Send(buf_s, nbyte_s, MPI_BYTE,
                 dest, shift,
                 atc->mpi_comm);
    }
    else if (nbyte_r > 0)
    {
        MPI_Recv(buf_r, nbyte_r, MPI_BYTE,
                 src, shift,
                 atc->mpi_comm, &stat);
    }
#endif
}

static void dd_pmeredist_pos_coeffs(struct gmx_pme_t *pme,
                                    int n, gmx_bool bX, rvec *x, real *data,
                                    pme_atomcomm_t *atc)
{
    int *commnode, *buf_index;
    int  nnodes_comm, i, nsend, local_pos, buf_pos, node, scount, rcount;

    commnode  = atc->node_dest;
    buf_index = atc->buf_index;

    nnodes_comm = std::min(2*atc->maxshift, atc->nslab-1);

    nsend = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        buf_index[commnode[i]] = nsend;
        nsend                 += atc->count[commnode[i]];
    }
    if (bX)
    {
        if (atc->count[atc->nodeid] + nsend != n)
        {
            gmx_fatal(FARGS, "%d particles communicated to PME rank %d are more than 2/3 times the cut-off out of the domain decomposition cell of their charge group in dimension %c.\n"
                      "This usually means that your system is not well equilibrated.",
                      n - (atc->count[atc->nodeid] + nsend),
                      pme->nodeid, 'x'+atc->dimind);
        }

        if (nsend > pme->buf_nalloc)
        {
            pme->buf_nalloc = over_alloc_dd(nsend);
            srenew(pme->bufv, pme->buf_nalloc);
            srenew(pme->bufr, pme->buf_nalloc);
        }

        atc->n = atc->count[atc->nodeid];
        for (i = 0; i < nnodes_comm; i++)
        {
            scount = atc->count[commnode[i]];
            /* Communicate the count */
            if (debug)
            {
                fprintf(debug, "dimind %d PME rank %d send to rank %d: %d\n",
                        atc->dimind, atc->nodeid, commnode[i], scount);
            }
            pme_dd_sendrecv(atc, FALSE, i,
                            &scount, sizeof(int),
                            &atc->rcount[i], sizeof(int));
            atc->n += atc->rcount[i];
        }

        pme_realloc_atomcomm_things(atc);
    }

    local_pos = 0;
    for (i = 0; i < n; i++)
    {
        node = atc->pd[i];
        if (node == atc->nodeid)
        {
            /* Copy direct to the receive buffer */
            if (bX)
            {
                copy_rvec(x[i], atc->x[local_pos]);
            }
            atc->coefficient[local_pos] = data[i];
            local_pos++;
        }
        else
        {
            /* Copy to the send buffer */
            if (bX)
            {
                copy_rvec(x[i], pme->bufv[buf_index[node]]);
            }
            pme->bufr[buf_index[node]] = data[i];
            buf_index[node]++;
        }
    }

    buf_pos = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        scount = atc->count[commnode[i]];
        rcount = atc->rcount[i];
        if (scount > 0 || rcount > 0)
        {
            if (bX)
            {
                /* Communicate the coordinates */
                pme_dd_sendrecv(atc, FALSE, i,
                                pme->bufv[buf_pos], scount*sizeof(rvec),
                                atc->x[local_pos], rcount*sizeof(rvec));
            }
            /* Communicate the coefficients */
            pme_dd_sendrecv(atc, FALSE, i,
                            pme->bufr+buf_pos, scount*sizeof(real),
                            atc->coefficient+local_pos, rcount*sizeof(real));
            buf_pos   += scount;
            local_pos += atc->rcount[i];
        }
    }
}

void dd_pmeredist_f(struct gmx_pme_t *pme, pme_atomcomm_t *atc,
                    int n, rvec *f,
                    gmx_bool bAddF)
{
    int *commnode, *buf_index;
    int  nnodes_comm, local_pos, buf_pos, i, scount, rcount, node;

    commnode  = atc->node_dest;
    buf_index = atc->buf_index;

    nnodes_comm = std::min(2*atc->maxshift, atc->nslab-1);

    local_pos = atc->count[atc->nodeid];
    buf_pos   = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        scount = atc->rcount[i];
        rcount = atc->count[commnode[i]];
        if (scount > 0 || rcount > 0)
        {
            /* Communicate the forces */
            pme_dd_sendrecv(atc, TRUE, i,
                            atc->f[local_pos], scount*sizeof(rvec),
                            pme->bufv[buf_pos], rcount*sizeof(rvec));
            local_pos += scount;
        }
        buf_index[commnode[i]] = buf_pos;
        buf_pos               += rcount;
    }

    local_pos = 0;
    if (bAddF)
    {
        for (i = 0; i < n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Add from the local force array */
                rvec_inc(f[i], atc->f[local_pos]);
                local_pos++;
            }
            else
            {
                /* Add from the receive buffer */
                rvec_inc(f[i], pme->bufv[buf_index[node]]);
                buf_index[node]++;
            }
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            node = atc->pd[i];
            if (node == atc->nodeid)
            {
                /* Copy from the local force array */
                copy_rvec(atc->f[local_pos], f[i]);
                local_pos++;
            }
            else
            {
                /* Copy from the receive buffer */
                copy_rvec(pme->bufv[buf_index[node]], f[i]);
                buf_index[node]++;
            }
        }
    }
}

void
do_redist_pos_coeffs(struct gmx_pme_t *pme, t_commrec *cr, int start, int homenr,
                     gmx_bool bFirst, rvec x[], real *data)
{
    int             d;
    pme_atomcomm_t *atc;
    atc = &pme->atc[0];

    for (d = pme->ndecompdim - 1; d >= 0; d--)
    {
        int             n_d;
        rvec           *x_d;
        real           *param_d;

        if (d == pme->ndecompdim - 1)
        {
            n_d     = homenr;
            x_d     = x + start;
            param_d = data;
        }
        else
        {
            n_d     = pme->atc[d + 1].n;
            x_d     = atc->x;
            param_d = atc->coefficient;
        }
        atc      = &pme->atc[d];
        atc->npd = n_d;
        if (atc->npd > atc->pd_nalloc)
        {
            atc->pd_nalloc = over_alloc_dd(atc->npd);
            srenew(atc->pd, atc->pd_nalloc);
        }
        pme_calc_pidx_wrapper(n_d, pme->recipbox, x_d, atc);
        where();
        /* Redistribute x (only once) and qA/c6A or qB/c6B */
        if (DOMAINDECOMP(cr))
        {
            dd_pmeredist_pos_coeffs(pme, n_d, bFirst, x_d, param_d, atc);
        }
    }
}
