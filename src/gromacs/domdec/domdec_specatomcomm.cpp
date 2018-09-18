/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006,2007,2008,2009,2010,2012,2013,2014,2015,2017,2018, by the GROMACS development team, led by
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

#include <cassert>

#include <algorithm>

#include "gromacs/domdec/dlb.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/hashedmap.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

void dd_move_f_specat(gmx_domdec_t *dd, gmx_domdec_specat_comm_t *spac,
                      rvec *f, rvec *fshift)
{
    gmx_specatsend_t *spas;
    rvec             *vbuf;
    int               n, n0, n1, dim, dir;
    ivec              vis;
    int               is;
    gmx_bool          bPBC, bScrew;

    n = spac->at_end;
    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        dim = dd->dim[d];
        if (dd->nc[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            spas = spac->spas[d];
            n0   = spas[0].nrecv;
            n1   = spas[1].nrecv;
            n   -= n1 + n0;
            vbuf = as_rvec_array(spac->vbuf.data());
            /* Send and receive the coordinates */
            dd_sendrecv2_rvec(dd, d,
                              f + n + n1, n0, vbuf, spas[0].a.size(),
                              f + n, n1, vbuf + spas[0].a.size(),
                              spas[1].a.size());
            for (dir = 0; dir < 2; dir++)
            {
                bPBC   = ((dir == 0 && dd->ci[dim] == 0) ||
                          (dir == 1 && dd->ci[dim] == dd->nc[dim]-1));
                bScrew = (bPBC && dd->bScrewPBC && dim == XX);

                spas = &spac->spas[d][dir];
                /* Sum the buffer into the required forces */
                if (!bPBC || (!bScrew && fshift == nullptr))
                {
                    for (int a : spas->a)
                    {
                        rvec_inc(f[a], *vbuf);
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
                        for (int a : spas->a)
                        {
                            rvec_inc(f[a], *vbuf);
                            rvec_inc(fshift[is], *vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Rotate the forces */
                        for (int a : spas->a)
                        {
                            f[a][XX] += (*vbuf)[XX];
                            f[a][YY] -= (*vbuf)[YY];
                            f[a][ZZ] -= (*vbuf)[ZZ];
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
            ddSendrecv(dd, d, dddirForward,
                       f + n, spas->nrecv,
                       as_rvec_array(spac->vbuf.data()), spas->a.size());
            /* Sum the buffer into the required forces */
            if (dd->bScrewPBC && dim == XX &&
                (dd->ci[dim] == 0 ||
                 dd->ci[dim] == dd->nc[dim]-1))
            {
                int i = 0;
                for (int a : spas->a)
                {
                    /* Rotate the force */
                    f[a][XX] += spac->vbuf[i][XX];
                    f[a][YY] -= spac->vbuf[i][YY];
                    f[a][ZZ] -= spac->vbuf[i][ZZ];
                    i++;
                }
            }
            else
            {
                int i = 0;
                for (int a : spas->a)
                {
                    rvec_inc(f[a], spac->vbuf[i]);
                    i++;
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
            rvec *vbuf = as_rvec_array(spac->vbuf.data());
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
                    rvec *x = (v == 0 ? x0 : x1);
                    /* Copy the required coordinates to the send buffer */
                    if (!bPBC || (v == 1 && !bX1IsCoord))
                    {
                        /* Only copy */
                        for (int a : spas->a)
                        {
                            copy_rvec(x[a], *vbuf);
                            vbuf++;
                        }
                    }
                    else if (!bScrew)
                    {
                        /* Shift coordinates */
                        for (int a : spas->a)
                        {
                            rvec_add(x[a], shift, *vbuf);
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Shift and rotate coordinates */
                        for (int a : spas->a)
                        {
                            (*vbuf)[XX] =               x[a][XX] + shift[XX];
                            (*vbuf)[YY] = box[YY][YY] - x[a][YY] + shift[YY];
                            (*vbuf)[ZZ] = box[ZZ][ZZ] - x[a][ZZ] + shift[ZZ];
                            vbuf++;
                        }
                    }
                }
            }
            /* Send and receive the coordinates */
            spas = spac->spas[d];
            ns0  = spas[0].a.size();
            nr0  = spas[0].nrecv;
            ns1  = spas[1].a.size();
            nr1  = spas[1].nrecv;
            if (nvec == 1)
            {
                rvec *vbuf = as_rvec_array(spac->vbuf.data());
                dd_sendrecv2_rvec(dd, d,
                                  vbuf + ns0, ns1, x0 + n, nr1,
                                  vbuf, ns0, x0 + n + nr1, nr0);
            }
            else
            {
                rvec *vbuf = as_rvec_array(spac->vbuf.data());
                /* Communicate both vectors in one buffer */
                rvec *rbuf = as_rvec_array(spac->vbuf2.data());
                dd_sendrecv2_rvec(dd, d,
                                  vbuf + 2*ns0, 2*ns1, rbuf, 2*nr1,
                                  vbuf, 2*ns0, rbuf + 2*nr1, 2*nr0);
                /* Split the buffer into the two vectors */
                nn = n;
                for (dir = 1; dir >= 0; dir--)
                {
                    nr = spas[dir].nrecv;
                    for (v = 0; v < 2; v++)
                    {
                        rvec *x = (v == 0 ? x0 : x1);
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
            rvec *vbuf = as_rvec_array(spac->vbuf.data());
            for (v = 0; v < nvec; v++)
            {
                rvec *x = (v == 0 ? x0 : x1);
                if (dd->bScrewPBC && dim == XX &&
                    (dd->ci[XX] == 0 || dd->ci[XX] == dd->nc[XX]-1))
                {
                    /* Here we only perform the rotation, the rest of the pbc
                     * is handled in the constraint or viste routines.
                     */
                    for (int a : spas->a)
                    {
                        (*vbuf)[XX] =               x[a][XX];
                        (*vbuf)[YY] = box[YY][YY] - x[a][YY];
                        (*vbuf)[ZZ] = box[ZZ][ZZ] - x[a][ZZ];
                        vbuf++;
                    }
                }
                else
                {
                    for (int a : spas->a)
                    {
                        copy_rvec(x[a], *vbuf);
                        vbuf++;
                    }
                }
            }
            /* Send and receive the coordinates */
            if (nvec == 1)
            {
                rvec *vbuf = as_rvec_array(spac->vbuf.data());
                ddSendrecv(dd, d, dddirBackward,
                           vbuf, spas->a.size(), x0 + n, spas->nrecv);
            }
            else
            {
                rvec *vbuf = as_rvec_array(spac->vbuf.data());
                /* Communicate both vectors in one buffer */
                rvec *rbuf = as_rvec_array(spac->vbuf2.data());
                ddSendrecv(dd, d, dddirBackward,
                           vbuf, 2*spas->a.size(), rbuf, 2*spas->nrecv);
                /* Split the buffer into the two vectors */
                nr = spas[0].nrecv;
                for (v = 0; v < 2; v++)
                {
                    rvec *x = (v == 0 ? x0 : x1);
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
                               std::vector<int>           *ireq,
                               gmx_domdec_specat_comm_t   *spac,
                               gmx::HashedMap<int>        *ga2la_specat,
                               int                         at_start,
                               int                         vbuf_fac,
                               const char                 *specat_type,
                               const char                 *add_err)
{
    int               nsend[2], nlast, nsend_zero[2] = {0, 0}, *nsend_ptr;
    int               dim, ndir, nr, ns, nrecv_local, n0, start, buf[2];
    int               nat_tot_specat, nat_tot_prev;
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
    const int numRequested = ireq->size();
    nsend[0]               = ireq->size();
    nsend[1]               = nsend[0];
    nlast                  = nsend[1];
    for (int d = dd->ndim-1; d >= 0; d--)
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
        for (int dir = 0; dir < ndir; dir++)
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
            ddSendrecv(dd, d, dir == 0 ? dddirForward : dddirBackward,
                       nsend_ptr, 2, spac->nreq[d][dir], 2);
            nr = spac->nreq[d][dir][1];
            ireq->resize(nlast + nr);
            /* Communicate the indices */
            ddSendrecv(dd, d, dir == 0 ? dddirForward : dddirBackward,
                       ireq->data(), nsend_ptr[1], ireq->data() + nlast, nr);
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
    for (int d = 0; d < dd->ndim; d++)
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
        for (int dir = ndir - 1; dir >= 0; dir--)
        {
            /* To avoid cost of clearing by resize(), we only increase size */
            if (static_cast<size_t>(nat_tot_specat) > spac->sendAtom.size())
            {
                /* Note: resize initializes new elements to false, which is actually needed here */
                spac->sendAtom.resize(nat_tot_specat);
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
            spas->a.clear();
            spac->ibuf.clear();
            nsend[0]    = 0;
            for (int i = 0; i < nr; i++)
            {
                const int indr = (*ireq)[start + i];
                int       ind;
                /* Check if this is a home atom and if so ind will be set */
                if (const int *homeIndex = dd->ga2la->findHome(indr))
                {
                    ind = *homeIndex;
                }
                else
                {
                    /* Search in the communicated atoms */
                    if (const int *a = ga2la_specat->find(indr))
                    {
                        ind = *a;
                    }
                    else
                    {
                        ind = -1;
                    }
                }
                if (ind >= 0)
                {
                    if (i < n0 || !spac->sendAtom[ind])
                    {
                        /* Store the local index so we know which coordinates
                         * to send out later.
                         */
                        spas->a.push_back(ind);
                        spac->sendAtom[ind] = true;
                        /* Store the global index so we can send it now */
                        spac->ibuf.push_back(indr);
                        if (i < n0)
                        {
                            nsend[0]++;
                        }
                    }
                }
            }
            nlast = start;
            /* Clear the local flags */
            for (int a : spas->a)
            {
                spac->sendAtom[a] = false;
            }
            /* Send and receive the number of indices to communicate */
            nsend[1] = spas->a.size();
            ddSendrecv(dd, d, dir == 0 ? dddirBackward : dddirForward,
                       nsend, 2, buf, 2);
            if (debug)
            {
                fprintf(debug, "Send to rank %d, %d (%d) indices, "
                        "receive from rank %d, %d (%d) indices\n",
                        dd->neighbor[d][1-dir], nsend[1], nsend[0],
                        dd->neighbor[d][dir], buf[1], buf[0]);
                if (gmx_debug_at)
                {
                    for (int i : spac->ibuf)
                    {
                        fprintf(debug, " %d", i + 1);
                    }
                    fprintf(debug, "\n");
                }
            }
            nrecv_local += buf[0];
            spas->nrecv  = buf[1];
            dd->globalAtomIndices.resize(nat_tot_specat + spas->nrecv);
            /* Send and receive the indices */
            ddSendrecv(dd, d, dir == 0 ? dddirBackward : dddirForward,
                       spac->ibuf.data(), spac->ibuf.size(),
                       dd->globalAtomIndices.data() + nat_tot_specat, spas->nrecv);
            nat_tot_specat += spas->nrecv;
        }

        /* Increase the x/f communication buffer sizes, when necessary */
        ns = spac->spas[d][0].a.size();
        nr = spac->spas[d][0].nrecv;
        if (ndir == 2)
        {
            ns += spac->spas[d][1].a.size();
            nr += spac->spas[d][1].nrecv;
        }
        if (vbuf_fac*ns > gmx::index(spac->vbuf.size()))
        {
            spac->vbuf.resize(vbuf_fac*ns);
        }
        if (vbuf_fac == 2 && vbuf_fac*nr > gmx::index(spac->vbuf2.size()))
        {
            spac->vbuf2.resize(vbuf_fac*nr);
        }

        /* Make a global to local index for the communication atoms */
        for (int i = nat_tot_prev; i < nat_tot_specat; i++)
        {
            ga2la_specat->insert_or_assign(dd->globalAtomIndices[i], i);
        }
    }

    /* Check that in the end we got the number of atoms we asked for */
    if (nrecv_local != numRequested)
    {
        if (debug)
        {
            fprintf(debug, "Requested %d, received %d (tot recv %d)\n",
                    numRequested, nrecv_local, nat_tot_specat - at_start);
            if (gmx_debug_at)
            {
                for (int i = 0; i < numRequested; i++)
                {
                    const int *ind = ga2la_specat->find((*ireq)[i]);
                    fprintf(debug, " %s%d",
                            ind ? "" : "!",
                            (*ireq)[i] + 1);
                }
                fprintf(debug, "\n");
            }
        }
        fprintf(stderr, "\nDD cell %d %d %d: Neighboring cells do not have atoms:",
                dd->ci[XX], dd->ci[YY], dd->ci[ZZ]);
        for (int i = 0; i < numRequested; i++)
        {
            if (!ga2la_specat->find((*ireq)[i]))
            {
                fprintf(stderr, " %d", (*ireq)[i] + 1);
            }
        }
        fprintf(stderr, "\n");
        gmx_fatal(FARGS, "DD cell %d %d %d could only obtain %d of the %d atoms that are connected via %ss from the neighboring cells. This probably means your %s lengths are too long compared to the domain decomposition cell size. Decrease the number of domain decomposition grid cells%s%s.",
                  dd->ci[XX], dd->ci[YY], dd->ci[ZZ],
                  nrecv_local, numRequested, specat_type,
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
