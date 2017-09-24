/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
 * \brief This file implements functions for domdec to use
 * while managing inter-atomic constraints.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_specatomcomm.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "hash.h"

void dd_move_f_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
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
                if (!bPBC || (!bScrew && fshift == nullptr))
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

void dd_move_x_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                      const matrix box,
                      rvec *x0,
                      rvec *x1, gmx_bool bX1IsCoord)
{
    gmx_specatsend_t *spas;
    rvec             *x, *vbuf, *rbuf;
    int               nvec, v, n, nn, ns0, ns1, nr0, nr1, nr, d, dim, dir, i;
    gmx_bool          bPBC, bScrew = FALSE;
    rvec              shift = {0, 0, 0};

    nvec = 1;
    if (x1 != nullptr)
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

int setup_specat_communication(gmx_domdec_t               *dd,
                               ind_req_t                  *ireq,
                               gmx_domdec_specat_comm_t   *spac,
                               gmx_hash_t                 *ga2la_specat,
                               int                         at_start,
                               int                         vbuf_fac,
                               const char                 *specat_type,
                               const char                 *add_err)
{
    int               nsend[2], nlast, nsend_zero[2] = {0, 0}, *nsend_ptr;
    int               d, dim, ndir, dir, nr, ns, i, nrecv_local, n0, start, indr, ind, buf[2];
    int               nat_tot_specat, nat_tot_prev, nalloc_old;
    gmx_bool          bPBC;
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
                fprintf(debug, "Send to rank %d, %d (%d) indices, "
                        "receive from rank %d, %d (%d) indices\n",
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
                  dd_dlb_is_on(dd) ? " or use the -rcon option of mdrun" : "");
    }

    spac->at_start = at_start;
    spac->at_end   = nat_tot_specat;

    if (debug)
    {
        fprintf(debug, "Done setup_specat_communication\n");
    }

    return nat_tot_specat;
}

gmx_domdec_specat_comm_t *specat_comm_init(int nthread)
{
    gmx_domdec_specat_comm_t *spac;

    snew(spac, 1);
    spac->nthread = nthread;
    snew(spac->ireq, spac->nthread);

    return spac;
}
