/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2008
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
#include <assert.h>

#include "smalloc.h"
#include "vec.h"
#include "constr.h"
#include "domdec.h"
#include "domdec_network.h"
#include "mtop_util.h"
#include "gmx_ga2la.h"
#include "gmx_hash.h"
#include "gmx_omp_nthreads.h"

typedef struct {
    int  nsend;
    int *a;
    int  a_nalloc;
    int  nrecv;
} gmx_specatsend_t;

typedef struct {
    int *ind;
    int  nalloc;
    int  n;
} ind_req_t;

typedef struct gmx_domdec_specat_comm {
    /* The number of indices to receive during the setup */
    int              nreq[DIM][2][2];
    /* The atoms to send */
    gmx_specatsend_t spas[DIM][2];
    gmx_bool        *bSendAtom;
    int              bSendAtom_nalloc;
    /* Send buffers */
    int             *ibuf;
    int              ibuf_nalloc;
    rvec            *vbuf;
    int              vbuf_nalloc;
    rvec            *vbuf2;
    int              vbuf2_nalloc;
    /* The range in the local buffer(s) for received atoms */
    int              at_start;
    int              at_end;

    /* The atom indices we need from the surrounding cells.
     * We can gather the indices over nthread threads.
     */
    int        nthread;
    ind_req_t *ireq;
} gmx_domdec_specat_comm_t;

typedef struct gmx_domdec_constraints {
    int       *molb_con_offset;
    int       *molb_ncon_mol;
    /* The fully local and connected constraints */
    int        ncon;
    /* The global constraint number, only required for clearing gc_req */
    int       *con_gl;
    int       *con_nlocat;
    int        con_nalloc;
    /* Boolean that tells if a global constraint index has been requested */
    char      *gc_req;
    /* Global to local communicated constraint atom only index */
    gmx_hash_t ga2la;

    /* Multi-threading stuff */
    int      nthread;
    t_ilist *ils;
} gmx_domdec_constraints_t;


static void dd_move_f_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                             rvec *f, rvec *fshift)
{
    gmx_specatsend_t *spas;
    rvec             *vbuf;
    int               n, n0, n1, d, dim, dir, i;
    ivec              vis;
    int               is;
    gmx_bool          bPBC, bScrew;

    n = spac->at_end;
    for (d = dd->ndim-1; d >= 0; d--)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            spas = spac->spas[d];
            n0   = spas[0].nrecv;
            n1   = spas[1].nrecv;
            n   -= n1 + n0;
            vbuf = spac->vbuf;
            /* Send and receive the coordinates */
            dd_sendrecv2_rvec(dd, d,
                              f+n+n1, n0, vbuf, spas[0].nsend,
                              f+n, n1, vbuf+spas[0].nsend, spas[1].nsend);
            for (dir = 0; dir < 2; dir++)
            {
                bPBC   = ((dir == 0 && dd->ci[dim] == 0) ||
                          (dir == 1 && dd->ci[dim] == dd->nc[dim]-1));
                bScrew = (bPBC && dd->bScrewPBC && dim == XX);

                spas = &spac->spas[d][dir];
                /* Sum the buffer into the required forces */
                if (!bPBC || (!bScrew && fshift == NULL))
                {
                    for (i = 0; i < spas->nsend; i++)
                    {
                        rvec_inc(f[spas->a[i]], *vbuf);
                        vbuf++;
                    }
                }
                else
                {
                    clear_ivec(vis);
                    vis[dim] = (dir == 0 ? 1 : -1);
                    is       = IVEC2IS(vis);
                    if (!bScrew)
                    {
                        /* Sum and add to shift forces */
                        for (i = 0; i < spas->nsend; i++)
                        {
                            rvec_inc(f[spas->a[i]], *vbuf);
                            rvec_inc(fshift[is], *vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Rotate the forces */
                        for (i = 0; i < spas->nsend; i++)
                        {
                            f[spas->a[i]][XX] += (*vbuf)[XX];
                            f[spas->a[i]][YY] -= (*vbuf)[YY];
                            f[spas->a[i]][ZZ] -= (*vbuf)[ZZ];
                            if (fshift)
                            {
                                rvec_inc(fshift[is], *vbuf);
                            }
                            vbuf++;
                        }
                    }
                }
            }
        }
        else
        {
            /* Two cells, so we only need to communicate one way */
            spas = &spac->spas[d][0];
            n   -= spas->nrecv;
            /* Send and receive the coordinates */
            dd_sendrecv_rvec(dd, d, dddirForward,
                             f+n, spas->nrecv, spac->vbuf, spas->nsend);
            /* Sum the buffer into the required forces */
            if (dd->bScrewPBC && dim == XX &&
                (dd->ci[dim] == 0 ||
                 dd->ci[dim] == dd->nc[dim]-1))
            {
                for (i = 0; i < spas->nsend; i++)
                {
                    /* Rotate the force */
                    f[spas->a[i]][XX] += spac->vbuf[i][XX];
                    f[spas->a[i]][YY] -= spac->vbuf[i][YY];
                    f[spas->a[i]][ZZ] -= spac->vbuf[i][ZZ];
                }
            }
            else
            {
                for (i = 0; i < spas->nsend; i++)
                {
                    rvec_inc(f[spas->a[i]], spac->vbuf[i]);
                }
            }
        }
    }
}

void dd_move_f_vsites(gmx_domdec_t *dd, rvec *f, rvec *fshift)
{
    if (dd->vsite_comm)
    {
        dd_move_f_specat(dd, dd->vsite_comm, f, fshift);
    }
}

void dd_clear_f_vsites(gmx_domdec_t *dd, rvec *f)
{
    int i;

    if (dd->vsite_comm)
    {
        for (i = dd->vsite_comm->at_start; i < dd->vsite_comm->at_end; i++)
        {
            clear_rvec(f[i]);
        }
    }
}

static void dd_move_x_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                             matrix box,
                             rvec *x0,
                             rvec *x1, gmx_bool bX1IsCoord)
{
    gmx_specatsend_t *spas;
    rvec             *x, *vbuf, *rbuf;
    int               nvec, v, n, nn, ns0, ns1, nr0, nr1, nr, d, dim, dir, i;
    gmx_bool          bPBC, bScrew = FALSE;
    rvec              shift = {0, 0, 0};

    nvec = 1;
    if (x1 != NULL)
    {
        nvec++;
    }

    n = spac->at_start;
    for (d = 0; d < dd->ndim; d++)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            vbuf = spac->vbuf;
            for (dir = 0; dir < 2; dir++)
            {
                if (dir == 0 && dd->ci[dim] == 0)
                {
                    bPBC   = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    copy_rvec(box[dim], shift);
                }
                else if (dir == 1 && dd->ci[dim] == dd->nc[dim]-1)
                {
                    bPBC   = TRUE;
                    bScrew = (dd->bScrewPBC && dim == XX);
                    for (i = 0; i < DIM; i++)
                    {
                        shift[i] = -box[dim][i];
                    }
                }
                else
                {
                    bPBC   = FALSE;
                    bScrew = FALSE;
                }
                spas = &spac->spas[d][dir];
                for (v = 0; v < nvec; v++)
                {
                    x = (v == 0 ? x0 : x1);
                    /* Copy the required coordinates to the send buffer */
                    if (!bPBC || (v == 1 && !bX1IsCoord))
                    {
                        /* Only copy */
                        for (i = 0; i < spas->nsend; i++)
                        {
                            copy_rvec(x[spas->a[i]], *vbuf);
                            vbuf++;
                        }
                    }
                    else if (!bScrew)
                    {
                        /* Shift coordinates */
                        for (i = 0; i < spas->nsend; i++)
                        {
                            rvec_add(x[spas->a[i]], shift, *vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Shift and rotate coordinates */
                        for (i = 0; i < spas->nsend; i++)
                        {
                            (*vbuf)[XX] =               x[spas->a[i]][XX] + shift[XX];
                            (*vbuf)[YY] = box[YY][YY] - x[spas->a[i]][YY] + shift[YY];
                            (*vbuf)[ZZ] = box[ZZ][ZZ] - x[spas->a[i]][ZZ] + shift[ZZ];
                            vbuf++;
                        }
                    }
                }
            }
            /* Send and receive the coordinates */
            spas = spac->spas[d];
            ns0  = spas[0].nsend;
            nr0  = spas[0].nrecv;
            ns1  = spas[1].nsend;
            nr1  = spas[1].nrecv;
            if (nvec == 1)
            {
                dd_sendrecv2_rvec(dd, d,
                                  spac->vbuf+ns0, ns1, x0+n, nr1,
                                  spac->vbuf, ns0, x0+n+nr1, nr0);
            }
            else
            {
                /* Communicate both vectors in one buffer */
                rbuf = spac->vbuf2;
                dd_sendrecv2_rvec(dd, d,
                                  spac->vbuf+2*ns0, 2*ns1, rbuf, 2*nr1,
                                  spac->vbuf, 2*ns0, rbuf+2*nr1, 2*nr0);
                /* Split the buffer into the two vectors */
                nn = n;
                for (dir = 1; dir >= 0; dir--)
                {
                    nr = spas[dir].nrecv;
                    for (v = 0; v < 2; v++)
                    {
                        x = (v == 0 ? x0 : x1);
                        for (i = 0; i < nr; i++)
                        {
                            copy_rvec(*rbuf, x[nn+i]);
                            rbuf++;
                        }
                    }
                    nn += nr;
                }
            }
            n += nr0 + nr1;
        }
        else
        {
            spas = &spac->spas[d][0];
            /* Copy the required coordinates to the send buffer */
            vbuf = spac->vbuf;
            for (v = 0; v < nvec; v++)
            {
                x = (v == 0 ? x0 : x1);
                if (dd->bScrewPBC && dim == XX &&
                    (dd->ci[XX] == 0 || dd->ci[XX] == dd->nc[XX]-1))
                {
                    /* Here we only perform the rotation, the rest of the pbc
                     * is handled in the constraint or viste routines.
                     */
                    for (i = 0; i < spas->nsend; i++)
                    {
                        (*vbuf)[XX] =               x[spas->a[i]][XX];
                        (*vbuf)[YY] = box[YY][YY] - x[spas->a[i]][YY];
                        (*vbuf)[ZZ] = box[ZZ][ZZ] - x[spas->a[i]][ZZ];
                        vbuf++;
                    }
                }
                else
                {
                    for (i = 0; i < spas->nsend; i++)
                    {
                        copy_rvec(x[spas->a[i]], *vbuf);
                        vbuf++;
                    }
                }
            }
            /* Send and receive the coordinates */
            if (nvec == 1)
            {
                dd_sendrecv_rvec(dd, d, dddirBackward,
                                 spac->vbuf, spas->nsend, x0+n, spas->nrecv);
            }
            else
            {
                /* Communicate both vectors in one buffer */
                rbuf = spac->vbuf2;
                dd_sendrecv_rvec(dd, d, dddirBackward,
                                 spac->vbuf, 2*spas->nsend, rbuf, 2*spas->nrecv);
                /* Split the buffer into the two vectors */
                nr = spas[0].nrecv;
                for (v = 0; v < 2; v++)
                {
                    x = (v == 0 ? x0 : x1);
                    for (i = 0; i < nr; i++)
                    {
                        copy_rvec(*rbuf, x[n+i]);
                        rbuf++;
                    }
                }
            }
            n += spas->nrecv;
        }
    }
}

void dd_move_x_constraints(gmx_domdec_t *dd, matrix box,
                           rvec *x0, rvec *x1, gmx_bool bX1IsCoord)
{
    if (dd->constraint_comm)
    {
        dd_move_x_specat(dd, dd->constraint_comm, box, x0, x1, bX1IsCoord);
    }
}

void dd_move_x_vsites(gmx_domdec_t *dd, matrix box, rvec *x)
{
    if (dd->vsite_comm)
    {
        dd_move_x_specat(dd, dd->vsite_comm, box, x, NULL, FALSE);
    }
}

int *dd_constraints_nlocalatoms(gmx_domdec_t *dd)
{
    if (dd->constraints)
    {
        return dd->constraints->con_nlocat;
    }
    else
    {
        return NULL;
    }
}

void dd_clear_local_constraint_indices(gmx_domdec_t *dd)
{
    gmx_domdec_constraints_t *dc;
    int i;

    dc = dd->constraints;

    for (i = 0; i < dc->ncon; i++)
    {
        dc->gc_req[dc->con_gl[i]] = 0;
    }

    if (dd->constraint_comm)
    {
        gmx_hash_clear_and_optimize(dc->ga2la);
    }
}

void dd_clear_local_vsite_indices(gmx_domdec_t *dd)
{
    int i;

    if (dd->vsite_comm)
    {
        gmx_hash_clear_and_optimize(dd->ga2la_vsite);
    }
}

static int setup_specat_communication(gmx_domdec_t             *dd,
                                      ind_req_t                *ireq,
                                      gmx_domdec_specat_comm_t *spac,
                                      gmx_hash_t                ga2la_specat,
                                      int                       at_start,
                                      int                       vbuf_fac,
                                      const char               *specat_type,
                                      const char               *add_err)
{
    int               nsend[2], nlast, nsend_zero[2] = {0, 0}, *nsend_ptr;
    int               d, dim, ndir, dir, nr, ns, i, nrecv_local, n0, start, indr, ind, buf[2];
    int               nat_tot_specat, nat_tot_prev, nalloc_old;
    gmx_bool          bPBC, bFirst;
    gmx_specatsend_t *spas;

    if (debug)
    {
        fprintf(debug, "Begin setup_specat_communication for %s\n", specat_type);
    }

    /* nsend[0]: the number of atoms requested by this node only,
     *           we communicate this for more efficients checks
     * nsend[1]: the total number of requested atoms
     */
    nsend[0] = ireq->n;
    nsend[1] = nsend[0];
    nlast    = nsend[1];
    for (d = dd->ndim-1; d >= 0; d--)
    {
        /* Pulse the grid forward and backward */
        dim  = dd->dim[d];
        bPBC = (dim < dd->npbcdim);
        if (dd->nc[dim] == 2)
        {
            /* Only 2 cells, so we only need to communicate once */
            ndir = 1;
        }
        else
        {
            ndir = 2;
        }
        for (dir = 0; dir < ndir; dir++)
        {
            if (!bPBC &&
                dd->nc[dim] > 2 &&
                ((dir == 0 && dd->ci[dim] == dd->nc[dim] - 1) ||
                 (dir == 1 && dd->ci[dim] == 0)))
            {
                /* No pbc: the fist/last cell should not request atoms */
                nsend_ptr = nsend_zero;
            }
            else
            {
                nsend_ptr = nsend;
            }
            /* Communicate the number of indices */
            dd_sendrecv_int(dd, d, dir == 0 ? dddirForward : dddirBackward,
                            nsend_ptr, 2, spac->nreq[d][dir], 2);
            nr = spac->nreq[d][dir][1];
            if (nlast+nr > ireq->nalloc)
            {
                ireq->nalloc = over_alloc_dd(nlast+nr);
                srenew(ireq->ind, ireq->nalloc);
            }
            /* Communicate the indices */
            dd_sendrecv_int(dd, d, dir == 0 ? dddirForward : dddirBackward,
                            ireq->ind, nsend_ptr[1], ireq->ind+nlast, nr);
            nlast += nr;
        }
        nsend[1] = nlast;
    }
    if (debug)
    {
        fprintf(debug, "Communicated the counts\n");
    }

    /* Search for the requested atoms and communicate the indices we have */
    nat_tot_specat = at_start;
    nrecv_local    = 0;
    for (d = 0; d < dd->ndim; d++)
    {
        bFirst = (d == 0);
        /* Pulse the grid forward and backward */
        if (dd->dim[d] >= dd->npbcdim || dd->nc[dd->dim[d]] > 2)
        {
            ndir = 2;
        }
        else
        {
            ndir = 1;
        }
        nat_tot_prev = nat_tot_specat;
        for (dir = ndir-1; dir >= 0; dir--)
        {
            if (nat_tot_specat > spac->bSendAtom_nalloc)
            {
                nalloc_old             = spac->bSendAtom_nalloc;
                spac->bSendAtom_nalloc = over_alloc_dd(nat_tot_specat);
                srenew(spac->bSendAtom, spac->bSendAtom_nalloc);
                for (i = nalloc_old; i < spac->bSendAtom_nalloc; i++)
                {
                    spac->bSendAtom[i] = FALSE;
                }
            }
            spas = &spac->spas[d][dir];
            n0   = spac->nreq[d][dir][0];
            nr   = spac->nreq[d][dir][1];
            if (debug)
            {
                fprintf(debug, "dim=%d, dir=%d, searching for %d atoms\n",
                        d, dir, nr);
            }
            start       = nlast - nr;
            spas->nsend = 0;
            nsend[0]    = 0;
            for (i = 0; i < nr; i++)
            {
                indr = ireq->ind[start+i];
                ind  = -1;
                /* Check if this is a home atom and if so ind will be set */
                if (!ga2la_get_home(dd->ga2la, indr, &ind))
                {
                    /* Search in the communicated atoms */
                    ind = gmx_hash_get_minone(ga2la_specat, indr);
                }
                if (ind >= 0)
                {
                    if (i < n0 || !spac->bSendAtom[ind])
                    {
                        if (spas->nsend+1 > spas->a_nalloc)
                        {
                            spas->a_nalloc = over_alloc_large(spas->nsend+1);
                            srenew(spas->a, spas->a_nalloc);
                        }
                        /* Store the local index so we know which coordinates
                         * to send out later.
                         */
                        spas->a[spas->nsend] = ind;
                        spac->bSendAtom[ind] = TRUE;
                        if (spas->nsend+1 > spac->ibuf_nalloc)
                        {
                            spac->ibuf_nalloc = over_alloc_large(spas->nsend+1);
                            srenew(spac->ibuf, spac->ibuf_nalloc);
                        }
                        /* Store the global index so we can send it now */
                        spac->ibuf[spas->nsend] = indr;
                        if (i < n0)
                        {
                            nsend[0]++;
                        }
                        spas->nsend++;
                    }
                }
            }
            nlast = start;
            /* Clear the local flags */
            for (i = 0; i < spas->nsend; i++)
            {
                spac->bSendAtom[spas->a[i]] = FALSE;
            }
            /* Send and receive the number of indices to communicate */
            nsend[1] = spas->nsend;
            dd_sendrecv_int(dd, d, dir == 0 ? dddirBackward : dddirForward,
                            nsend, 2, buf, 2);
            if (debug)
            {
                fprintf(debug, "Send to node %d, %d (%d) indices, "
                        "receive from node %d, %d (%d) indices\n",
                        dd->neighbor[d][1-dir], nsend[1], nsend[0],
                        dd->neighbor[d][dir], buf[1], buf[0]);
                if (gmx_debug_at)
                {
                    for (i = 0; i < spas->nsend; i++)
                    {
                        fprintf(debug, " %d", spac->ibuf[i]+1);
                    }
                    fprintf(debug, "\n");
                }
            }
            nrecv_local += buf[0];
            spas->nrecv  = buf[1];
            if (nat_tot_specat + spas->nrecv > dd->gatindex_nalloc)
            {
                dd->gatindex_nalloc =
                    over_alloc_dd(nat_tot_specat + spas->nrecv);
                srenew(dd->gatindex, dd->gatindex_nalloc);
            }
            /* Send and receive the indices */
            dd_sendrecv_int(dd, d, dir == 0 ? dddirBackward : dddirForward,
                            spac->ibuf, spas->nsend,
                            dd->gatindex+nat_tot_specat, spas->nrecv);
            nat_tot_specat += spas->nrecv;
        }

        /* Allocate the x/f communication buffers */
        ns = spac->spas[d][0].nsend;
        nr = spac->spas[d][0].nrecv;
        if (ndir == 2)
        {
            ns += spac->spas[d][1].nsend;
            nr += spac->spas[d][1].nrecv;
        }
        if (vbuf_fac*ns > spac->vbuf_nalloc)
        {
            spac->vbuf_nalloc = over_alloc_dd(vbuf_fac*ns);
            srenew(spac->vbuf, spac->vbuf_nalloc);
        }
        if (vbuf_fac == 2 && vbuf_fac*nr > spac->vbuf2_nalloc)
        {
            spac->vbuf2_nalloc = over_alloc_dd(vbuf_fac*nr);
            srenew(spac->vbuf2, spac->vbuf2_nalloc);
        }

        /* Make a global to local index for the communication atoms */
        for (i = nat_tot_prev; i < nat_tot_specat; i++)
        {
            gmx_hash_change_or_set(ga2la_specat, dd->gatindex[i], i);
        }
    }

    /* Check that in the end we got the number of atoms we asked for */
    if (nrecv_local != ireq->n)
    {
        if (debug)
        {
            fprintf(debug, "Requested %d, received %d (tot recv %d)\n",
                    ireq->n, nrecv_local, nat_tot_specat-at_start);
            if (gmx_debug_at)
            {
                for (i = 0; i < ireq->n; i++)
                {
                    ind = gmx_hash_get_minone(ga2la_specat, ireq->ind[i]);
                    fprintf(debug, " %s%d",
                            (ind >= 0) ? "" : "!",
                            ireq->ind[i]+1);
                }
                fprintf(debug, "\n");
            }
        }
        fprintf(stderr, "\nDD cell %d %d %d: Neighboring cells do not have atoms:",
                dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
        for (i = 0; i < ireq->n; i++)
        {
            if (gmx_hash_get_minone(ga2la_specat, ireq->ind[i]) < 0)
            {
                fprintf(stderr, " %d", ireq->ind[i]+1);
            }
        }
        fprintf(stderr, "\n");
        gmx_fatal(FARGS, "DD cell %d %d %d could only obtain %d of the %d atoms that are connected via %ss from the neighboring cells. This probably means your %s lengths are too long compared to the domain decomposition cell size. Decrease the number of domain decomposition grid cells%s%s.",
                  dd->ci[XX], dd->ci[YY], dd->ci[ZZ],
                  nrecv_local, ireq->n, specat_type,
                  specat_type, add_err,
                  dd->bGridJump ? " or use the -rcon option of mdrun" : "");
    }

    spac->at_start = at_start;
    spac->at_end   = nat_tot_specat;

    if (debug)
    {
        fprintf(debug, "Done setup_specat_communication\n");
    }

    return nat_tot_specat;
}

static void walk_out(int con, int con_offset, int a, int offset, int nrec,
                     int ncon1, const t_iatom *ia1, const t_iatom *ia2,
                     const t_blocka *at2con,
                     const gmx_ga2la_t ga2la, gmx_bool bHomeConnect,
                     gmx_domdec_constraints_t *dc,
                     gmx_domdec_specat_comm_t *dcc,
                     t_ilist *il_local,
                     ind_req_t *ireq)
{
    int            a1_gl, a2_gl, a_loc, i, coni, b;
    const t_iatom *iap;

    if (dc->gc_req[con_offset+con] == 0)
    {
        /* Add this non-home constraint to the list */
        if (dc->ncon+1 > dc->con_nalloc)
        {
            dc->con_nalloc = over_alloc_large(dc->ncon+1);
            srenew(dc->con_gl, dc->con_nalloc);
            srenew(dc->con_nlocat, dc->con_nalloc);
        }
        dc->con_gl[dc->ncon]       = con_offset + con;
        dc->con_nlocat[dc->ncon]   = (bHomeConnect ? 1 : 0);
        dc->gc_req[con_offset+con] = 1;
        if (il_local->nr + 3 > il_local->nalloc)
        {
            il_local->nalloc = over_alloc_dd(il_local->nr+3);
            srenew(il_local->iatoms, il_local->nalloc);
        }
        iap = constr_iatomptr(ncon1, ia1, ia2, con);
        il_local->iatoms[il_local->nr++] = iap[0];
        a1_gl = offset + iap[1];
        a2_gl = offset + iap[2];
        /* The following indexing code can probably be optizimed */
        if (ga2la_get_home(ga2la, a1_gl, &a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a1_gl - 1;
        }
        if (ga2la_get_home(ga2la, a2_gl, &a_loc))
        {
            il_local->iatoms[il_local->nr++] = a_loc;
        }
        else
        {
            /* We set this index later */
            il_local->iatoms[il_local->nr++] = -a2_gl - 1;
        }
        dc->ncon++;
    }
    /* Check to not ask for the same atom more than once */
    if (gmx_hash_get_minone(dc->ga2la, offset+a) == -1)
    {
        assert(dcc);
        /* Add this non-home atom to the list */
        if (ireq->n+1 > ireq->nalloc)
        {
            ireq->nalloc = over_alloc_large(ireq->n+1);
            srenew(ireq->ind, ireq->nalloc);
        }
        ireq->ind[ireq->n++] = offset + a;
        /* Temporarily mark with -2, we get the index later */
        gmx_hash_set(dc->ga2la, offset+a, -2);
    }

    if (nrec > 0)
    {
        for (i = at2con->index[a]; i < at2con->index[a+1]; i++)
        {
            coni = at2con->a[i];
            if (coni != con)
            {
                /* Walk further */
                iap = constr_iatomptr(ncon1, ia1, ia2, coni);
                if (a == iap[1])
                {
                    b = iap[2];
                }
                else
                {
                    b = iap[1];
                }
                if (!ga2la_get_home(ga2la, offset+b, &a_loc))
                {
                    walk_out(coni, con_offset, b, offset, nrec-1,
                             ncon1, ia1, ia2, at2con,
                             ga2la, FALSE, dc, dcc, il_local, ireq);
                }
            }
        }
    }
}

static void atoms_to_settles(gmx_domdec_t *dd,
                             const gmx_mtop_t *mtop,
                             const int *cginfo,
                             const int **at2settle_mt,
                             int cg_start, int cg_end,
                             t_ilist *ils_local,
                             ind_req_t *ireq)
{
    gmx_ga2la_t           ga2la;
    gmx_mtop_atomlookup_t alook;
    int                   settle;
    int                   nral, sa;
    int                   cg, a, a_gl, a_glsa, a_gls[3], a_locs[3];
    int                   mb, molnr, a_mol, offset;
    const gmx_molblock_t *molb;
    const t_iatom        *ia1;
    gmx_bool              a_home[3];
    int                   nlocal;
    gmx_bool              bAssign;

    ga2la  = dd->ga2la;

    alook = gmx_mtop_atomlookup_settle_init(mtop);

    nral = NRAL(F_SETTLE);

    for (cg = cg_start; cg < cg_end; cg++)
    {
        if (GET_CGINFO_SETTLE(cginfo[cg]))
        {
            for (a = dd->cgindex[cg]; a < dd->cgindex[cg+1]; a++)
            {
                a_gl = dd->gatindex[a];

                gmx_mtop_atomnr_to_molblock_ind(alook, a_gl, &mb, &molnr, &a_mol);
                molb = &mtop->molblock[mb];

                settle = at2settle_mt[molb->type][a_mol];

                if (settle >= 0)
                {
                    offset = a_gl - a_mol;

                    ia1 = mtop->moltype[molb->type].ilist[F_SETTLE].iatoms;

                    bAssign = FALSE;
                    nlocal  = 0;
                    for (sa = 0; sa < nral; sa++)
                    {
                        a_glsa     = offset + ia1[settle*(1+nral)+1+sa];
                        a_gls[sa]  = a_glsa;
                        a_home[sa] = ga2la_get_home(ga2la, a_glsa, &a_locs[sa]);
                        if (a_home[sa])
                        {
                            if (nlocal == 0 && a_gl == a_glsa)
                            {
                                bAssign = TRUE;
                            }
                            nlocal++;
                        }
                    }

                    if (bAssign)
                    {
                        if (ils_local->nr+1+nral > ils_local->nalloc)
                        {
                            ils_local->nalloc = over_alloc_dd(ils_local->nr+1+nral);
                            srenew(ils_local->iatoms, ils_local->nalloc);
                        }

                        ils_local->iatoms[ils_local->nr++] = ia1[settle*4];

                        for (sa = 0; sa < nral; sa++)
                        {
                            if (ga2la_get_home(ga2la, a_gls[sa], &a_locs[sa]))
                            {
                                ils_local->iatoms[ils_local->nr++] = a_locs[sa];
                            }
                            else
                            {
                                ils_local->iatoms[ils_local->nr++] = -a_gls[sa] - 1;
                                /* Add this non-home atom to the list */
                                if (ireq->n+1 > ireq->nalloc)
                                {
                                    ireq->nalloc = over_alloc_large(ireq->n+1);
                                    srenew(ireq->ind, ireq->nalloc);
                                }
                                ireq->ind[ireq->n++] = a_gls[sa];
                                /* A check on double atom requests is
                                 * not required for settle.
                                 */
                            }
                        }
                    }
                }
            }
        }
    }

    gmx_mtop_atomlookup_destroy(alook);
}

static void atoms_to_constraints(gmx_domdec_t *dd,
                                 const gmx_mtop_t *mtop,
                                 const int *cginfo,
                                 const t_blocka *at2con_mt, int nrec,
                                 t_ilist *ilc_local,
                                 ind_req_t *ireq)
{
    const t_blocka           *at2con;
    gmx_ga2la_t               ga2la;
    gmx_mtop_atomlookup_t     alook;
    int                       ncon1;
    gmx_molblock_t           *molb;
    t_iatom                  *ia1, *ia2, *iap;
    int                       nhome, cg, a, a_gl, a_mol, a_loc, b_lo, offset, mb, molnr, b_mol, i, con, con_offset;
    gmx_domdec_constraints_t *dc;
    gmx_domdec_specat_comm_t *dcc;

    dc  = dd->constraints;
    dcc = dd->constraint_comm;

    ga2la  = dd->ga2la;

    alook = gmx_mtop_atomlookup_init(mtop);

    nhome = 0;
    for (cg = 0; cg < dd->ncg_home; cg++)
    {
        if (GET_CGINFO_CONSTR(cginfo[cg]))
        {
            for (a = dd->cgindex[cg]; a < dd->cgindex[cg+1]; a++)
            {
                a_gl = dd->gatindex[a];

                gmx_mtop_atomnr_to_molblock_ind(alook, a_gl, &mb, &molnr, &a_mol);
                molb = &mtop->molblock[mb];

                ncon1 = mtop->moltype[molb->type].ilist[F_CONSTR].nr/NRAL(F_SETTLE);

                ia1 = mtop->moltype[molb->type].ilist[F_CONSTR].iatoms;
                ia2 = mtop->moltype[molb->type].ilist[F_CONSTRNC].iatoms;

                /* Calculate the global constraint number offset for the molecule.
                 * This is only required for the global index to make sure
                 * that we use each constraint only once.
                 */
                con_offset =
                    dc->molb_con_offset[mb] + molnr*dc->molb_ncon_mol[mb];

                /* The global atom number offset for this molecule */
                offset = a_gl - a_mol;
                at2con = &at2con_mt[molb->type];
                for (i = at2con->index[a_mol]; i < at2con->index[a_mol+1]; i++)
                {
                    con = at2con->a[i];
                    iap = constr_iatomptr(ncon1, ia1, ia2, con);
                    if (a_mol == iap[1])
                    {
                        b_mol = iap[2];
                    }
                    else
                    {
                        b_mol = iap[1];
                    }
                    if (ga2la_get_home(ga2la, offset+b_mol, &a_loc))
                    {
                        /* Add this fully home constraint at the first atom */
                        if (a_mol < b_mol)
                        {
                            if (dc->ncon+1 > dc->con_nalloc)
                            {
                                dc->con_nalloc = over_alloc_large(dc->ncon+1);
                                srenew(dc->con_gl, dc->con_nalloc);
                                srenew(dc->con_nlocat, dc->con_nalloc);
                            }
                            dc->con_gl[dc->ncon]     = con_offset + con;
                            dc->con_nlocat[dc->ncon] = 2;
                            if (ilc_local->nr + 3 > ilc_local->nalloc)
                            {
                                ilc_local->nalloc = over_alloc_dd(ilc_local->nr + 3);
                                srenew(ilc_local->iatoms, ilc_local->nalloc);
                            }
                            b_lo = a_loc;
                            ilc_local->iatoms[ilc_local->nr++] = iap[0];
                            ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? a    : b_lo);
                            ilc_local->iatoms[ilc_local->nr++] = (a_gl == iap[1] ? b_lo : a   );
                            dc->ncon++;
                            nhome++;
                        }
                    }
                    else
                    {
                        /* We need the nrec constraints coupled to this constraint,
                         * so we need to walk out of the home cell by nrec+1 atoms,
                         * since already atom bg is not locally present.
                         * Therefore we call walk_out with nrec recursions to go
                         * after this first call.
                         */
                        walk_out(con, con_offset, b_mol, offset, nrec,
                                 ncon1, ia1, ia2, at2con,
                                 dd->ga2la, TRUE, dc, dcc, ilc_local, ireq);
                    }
                }
            }
        }
    }

    gmx_mtop_atomlookup_destroy(alook);

    if (debug)
    {
        fprintf(debug,
                "Constraints: home %3d border %3d atoms: %3d\n",
                nhome, dc->ncon-nhome,
                dd->constraint_comm ? ireq->n : 0);
    }
}

int dd_make_local_constraints(gmx_domdec_t *dd, int at_start,
                              const gmx_mtop_t *mtop,
                              const int *cginfo,
                              gmx_constr_t constr, int nrec,
                              t_ilist *il_local)
{
    gmx_domdec_constraints_t *dc;
    t_ilist                  *ilc_local, *ils_local;
    ind_req_t                *ireq;
    const t_blocka           *at2con_mt;
    const int               **at2settle_mt;
    gmx_hash_t                ga2la_specat;
    int at_end, i, j;
    t_iatom                  *iap;

    dc = dd->constraints;

    ilc_local = &il_local[F_CONSTR];
    ils_local = &il_local[F_SETTLE];

    dc->ncon      = 0;
    ilc_local->nr = 0;
    if (dd->constraint_comm)
    {
        at2con_mt = atom2constraints_moltype(constr);
        ireq      = &dd->constraint_comm->ireq[0];
        ireq->n   = 0;
    }
    else
    {
        at2con_mt = NULL;
        ireq      = NULL;
    }

    if (dd->bInterCGsettles)
    {
        at2settle_mt  = atom2settle_moltype(constr);
        ils_local->nr = 0;
    }
    else
    {
        /* Settle works inside charge groups, we assigned them already */
        at2settle_mt = NULL;
    }

    if (at2settle_mt == NULL)
    {
        atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec,
                             ilc_local, ireq);
    }
    else
    {
        int t0_set;
        int thread;

        /* Do the constraints, if present, on the first thread.
         * Do the settles on all other threads.
         */
        t0_set = ((at2con_mt != NULL && dc->nthread > 1) ? 1 : 0);

#pragma omp parallel for num_threads(dc->nthread) schedule(static)
        for (thread = 0; thread < dc->nthread; thread++)
        {
            if (at2con_mt && thread == 0)
            {
                atoms_to_constraints(dd, mtop, cginfo, at2con_mt, nrec,
                                     ilc_local, ireq);
            }

            if (thread >= t0_set)
            {
                int        cg0, cg1;
                t_ilist   *ilst;
                ind_req_t *ireqt;

                /* Distribute the settle check+assignments over
                 * dc->nthread or dc->nthread-1 threads.
                 */
                cg0 = (dd->ncg_home*(thread-t0_set  ))/(dc->nthread-t0_set);
                cg1 = (dd->ncg_home*(thread-t0_set+1))/(dc->nthread-t0_set);

                if (thread == t0_set)
                {
                    ilst = ils_local;
                }
                else
                {
                    ilst = &dc->ils[thread];
                }
                ilst->nr = 0;

                ireqt = &dd->constraint_comm->ireq[thread];
                if (thread > 0)
                {
                    ireqt->n = 0;
                }

                atoms_to_settles(dd, mtop, cginfo, at2settle_mt,
                                 cg0, cg1,
                                 ilst, ireqt);
            }
        }

        /* Combine the generate settles and requested indices */
        for (thread = 1; thread < dc->nthread; thread++)
        {
            t_ilist   *ilst;
            ind_req_t *ireqt;
            int        ia;

            if (thread > t0_set)
            {
                ilst = &dc->ils[thread];
                if (ils_local->nr + ilst->nr > ils_local->nalloc)
                {
                    ils_local->nalloc = over_alloc_large(ils_local->nr + ilst->nr);
                    srenew(ils_local->iatoms, ils_local->nalloc);
                }
                for (ia = 0; ia < ilst->nr; ia++)
                {
                    ils_local->iatoms[ils_local->nr+ia] = ilst->iatoms[ia];
                }
                ils_local->nr += ilst->nr;
            }

            ireqt = &dd->constraint_comm->ireq[thread];
            if (ireq->n+ireqt->n > ireq->nalloc)
            {
                ireq->nalloc = over_alloc_large(ireq->n+ireqt->n);
                srenew(ireq->ind, ireq->nalloc);
            }
            for (ia = 0; ia < ireqt->n; ia++)
            {
                ireq->ind[ireq->n+ia] = ireqt->ind[ia];
            }
            ireq->n += ireqt->n;
        }

        if (debug)
        {
            fprintf(debug, "Settles: total %3d\n", ils_local->nr/4);
        }
    }

    if (dd->constraint_comm)
    {
        int nral1;

        at_end =
            setup_specat_communication(dd, ireq, dd->constraint_comm,
                                       dd->constraints->ga2la,
                                       at_start, 2,
                                       "constraint", " or lincs-order");

        /* Fill in the missing indices */
        ga2la_specat = dd->constraints->ga2la;

        nral1 = 1 + NRAL(F_CONSTR);
        for (i = 0; i < ilc_local->nr; i += nral1)
        {
            iap = ilc_local->iatoms + i;
            for (j = 1; j < nral1; j++)
            {
                if (iap[j] < 0)
                {
                    iap[j] = gmx_hash_get_minone(ga2la_specat, -iap[j]-1);
                }
            }
        }

        nral1 = 1 + NRAL(F_SETTLE);
        for (i = 0; i < ils_local->nr; i += nral1)
        {
            iap = ils_local->iatoms + i;
            for (j = 1; j < nral1; j++)
            {
                if (iap[j] < 0)
                {
                    iap[j] = gmx_hash_get_minone(ga2la_specat, -iap[j]-1);
                }
            }
        }
    }
    else
    {
        at_end = at_start;
    }

    return at_end;
}

int dd_make_local_vsites(gmx_domdec_t *dd, int at_start, t_ilist *lil)
{
    gmx_domdec_specat_comm_t *spac;
    ind_req_t                *ireq;
    gmx_hash_t                ga2la_specat;
    int  ftype, nral, i, j, gat, a;
    t_ilist                  *lilf;
    t_iatom                  *iatoms;
    int  at_end;

    spac         = dd->vsite_comm;
    ireq         = &spac->ireq[0];
    ga2la_specat = dd->ga2la_vsite;

    ireq->n = 0;
    /* Loop over all the home vsites */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for (i = 0; i < lilf->nr; i += 1+nral)
            {
                iatoms = lilf->iatoms + i;
                /* Check if we have the other atoms */
                for (j = 1; j < 1+nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        /* This is not a home atom,
                         * we need to ask our neighbors.
                         */
                        a = -iatoms[j] - 1;
                        /* Check to not ask for the same atom more than once */
                        if (gmx_hash_get_minone(dd->ga2la_vsite, a) == -1)
                        {
                            /* Add this non-home atom to the list */
                            if (ireq->n+1 > ireq->nalloc)
                            {
                                ireq->nalloc = over_alloc_large(ireq->n+1);
                                srenew(ireq->ind, ireq->nalloc);
                            }
                            ireq->ind[ireq->n++] = a;
                            /* Temporarily mark with -2,
                             * we get the index later.
                             */
                            gmx_hash_set(ga2la_specat, a, -2);
                        }
                    }
                }
            }
        }
    }

    at_end = setup_specat_communication(dd, ireq, dd->vsite_comm, ga2la_specat,
                                        at_start, 1, "vsite", "");

    /* Fill in the missing indices */
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            nral = NRAL(ftype);
            lilf = &lil[ftype];
            for (i = 0; i < lilf->nr; i += 1+nral)
            {
                iatoms = lilf->iatoms + i;
                for (j = 1; j < 1+nral; j++)
                {
                    if (iatoms[j] < 0)
                    {
                        iatoms[j] = gmx_hash_get_minone(ga2la_specat, -iatoms[j]-1);
                    }
                }
            }
        }
    }

    return at_end;
}

static gmx_domdec_specat_comm_t *specat_comm_init(int nthread)
{
    gmx_domdec_specat_comm_t *spac;

    snew(spac, 1);
    spac->nthread = nthread;
    snew(spac->ireq, spac->nthread);

    return spac;
}

void init_domdec_constraints(gmx_domdec_t *dd,
                             gmx_mtop_t   *mtop,
                             gmx_constr_t  constr)
{
    gmx_domdec_constraints_t *dc;
    gmx_molblock_t           *molb;
    int mb, ncon, c, a;

    if (debug)
    {
        fprintf(debug, "Begin init_domdec_constraints\n");
    }

    snew(dd->constraints, 1);
    dc = dd->constraints;

    snew(dc->molb_con_offset, mtop->nmolblock);
    snew(dc->molb_ncon_mol, mtop->nmolblock);

    ncon = 0;
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        molb                    = &mtop->molblock[mb];
        dc->molb_con_offset[mb] = ncon;
        dc->molb_ncon_mol[mb]   =
            mtop->moltype[molb->type].ilist[F_CONSTR].nr/3 +
            mtop->moltype[molb->type].ilist[F_CONSTRNC].nr/3;
        ncon += molb->nmol*dc->molb_ncon_mol[mb];
    }

    if (ncon > 0)
    {
        snew(dc->gc_req, ncon);
        for (c = 0; c < ncon; c++)
        {
            dc->gc_req[c] = 0;
        }
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    dc->ga2la = gmx_hash_init(min(mtop->natoms/20,
                                  mtop->natoms/(2*dd->nnodes)));

    dc->nthread = gmx_omp_nthreads_get(emntDomdec);
    snew(dc->ils, dc->nthread);

    dd->constraint_comm = specat_comm_init(dc->nthread);
}

void init_domdec_vsites(gmx_domdec_t *dd, int n_intercg_vsite)
{
    int i;
    gmx_domdec_constraints_t *dc;

    if (debug)
    {
        fprintf(debug, "Begin init_domdec_vsites\n");
    }

    /* Use a hash table for the global to local index.
     * The number of keys is a rough estimate, it will be optimized later.
     */
    dd->ga2la_vsite = gmx_hash_init(min(n_intercg_vsite/20,
                                        n_intercg_vsite/(2*dd->nnodes)));

    dd->vsite_comm = specat_comm_init(1);
}
