/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2006- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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
#include <cstdio>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>

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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

void dd_move_f_specat(const gmx_domdec_t* dd, gmx_domdec_specat_comm_t* spac, gmx::RVec* f, gmx::RVec* fshift)
{
    ivec vis;

    int n = spac->at_end;
    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        const int dim = dd->dim[d];
        if (dd->numCells[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            const gmx_specatsend_t* spas = spac->spas[d];
            int                     n0   = spas[0].nrecv;
            int                     n1   = spas[1].nrecv;
            n -= n1 + n0;
            rvec* vbuf = as_rvec_array(spac->vbuf.data());
            /* Send and receive the coordinates */
            dd_sendrecv2_rvec(dd,
                              d,
                              as_rvec_array(f + n + n1),
                              n0,
                              vbuf,
                              spas[0].a.size(),
                              as_rvec_array(f + n),
                              n1,
                              vbuf + spas[0].a.size(),
                              spas[1].a.size());
            for (int dir = 0; dir < 2; dir++)
            {
                bool bPBC   = ((dir == 0 && dd->ci[dim] == 0)
                             || (dir == 1 && dd->ci[dim] == dd->numCells[dim] - 1));
                bool bScrew = (bPBC && dd->unitCellInfo.haveScrewPBC && dim == XX);

                const gmx_specatsend_t* spas = &spac->spas[d][dir];
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
                    int is   = gmx::ivecToShiftIndex(vis);
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
            const gmx_specatsend_t* spas = &spac->spas[d][0];
            n -= spas->nrecv;
            /* Send and receive the coordinates */
            ddSendrecv(dd,
                       d,
                       dddirForward,
                       gmx::arrayRefFromArray(f + n, spas->nrecv),
                       gmx::arrayRefFromArray(spac->vbuf.data(), spas->a.size()));
            /* Sum the buffer into the required forces */
            if (dd->unitCellInfo.haveScrewPBC && dim == XX
                && (dd->ci[dim] == 0 || dd->ci[dim] == dd->numCells[dim] - 1))
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

void dd_move_x_specat(const gmx_domdec_t*       dd,
                      gmx_domdec_specat_comm_t* spac,
                      const matrix              box,
                      gmx::RVec*                x0,
                      gmx::RVec*                x1,
                      bool                      bX1IsCoord)
{
    gmx::RVec shift = { 0, 0, 0 };

    int nvec = 1;
    if (x1 != nullptr)
    {
        nvec++;
    }

    int n = spac->at_start;
    for (int d = 0; d < dd->ndim; d++)
    {
        int dim = dd->dim[d];
        if (dd->numCells[dim] > 2)
        {
            /* Pulse the grid forward and backward */
            gmx::RVec* vbuf = spac->vbuf.data();
            for (int dir = 0; dir < 2; dir++)
            {
                bool bPBC   = false;
                bool bScrew = false;
                if (dir == 0 && dd->ci[dim] == 0)
                {
                    bPBC   = true;
                    bScrew = (dd->unitCellInfo.haveScrewPBC && dim == XX);
                    copy_rvec(box[dim], shift);
                }
                else if (dir == 1 && dd->ci[dim] == dd->numCells[dim] - 1)
                {
                    bPBC   = true;
                    bScrew = (dd->unitCellInfo.haveScrewPBC && dim == XX);
                    for (int i = 0; i < DIM; i++)
                    {
                        shift[i] = -box[dim][i];
                    }
                }
                const gmx_specatsend_t* spas = &spac->spas[d][dir];
                for (int v = 0; v < nvec; v++)
                {
                    gmx::RVec* x = (v == 0 ? x0 : x1);
                    /* Copy the required coordinates to the send buffer */
                    if (!bPBC || (v == 1 && !bX1IsCoord))
                    {
                        /* Only copy */
                        for (int a : spas->a)
                        {
                            *vbuf = x[a];
                            vbuf++;
                        }
                    }
                    else if (!bScrew)
                    {
                        /* Shift coordinates */
                        for (int a : spas->a)
                        {
                            *vbuf = x[a] + shift;
                            vbuf++;
                        }
                    }
                    else
                    {
                        /* Shift and rotate coordinates */
                        for (int a : spas->a)
                        {
                            (*vbuf)[XX] = x[a][XX] + shift[XX];
                            (*vbuf)[YY] = box[YY][YY] - x[a][YY] + shift[YY];
                            (*vbuf)[ZZ] = box[ZZ][ZZ] - x[a][ZZ] + shift[ZZ];
                            vbuf++;
                        }
                    }
                }
            }
            /* Send and receive the coordinates */
            const gmx_specatsend_t* spas = spac->spas[d];

            int ns0 = spas[0].a.size();
            int nr0 = spas[0].nrecv;
            int ns1 = spas[1].a.size();
            int nr1 = spas[1].nrecv;
            if (nvec == 1)
            {
                rvec* vbuf = as_rvec_array(spac->vbuf.data());
                dd_sendrecv2_rvec(
                        dd, d, vbuf + ns0, ns1, as_rvec_array(x0 + n), nr1, vbuf, ns0, as_rvec_array(x0 + n + nr1), nr0);
            }
            else
            {
                rvec* vbuf = as_rvec_array(spac->vbuf.data());
                /* Communicate both vectors in one buffer */
                rvec* rbuf = as_rvec_array(spac->vbuf2.data());
                dd_sendrecv2_rvec(
                        dd, d, vbuf + 2 * ns0, 2 * ns1, rbuf, 2 * nr1, vbuf, 2 * ns0, rbuf + 2 * nr1, 2 * nr0);
                /* Split the buffer into the two vectors */
                int nn = n;
                for (int dir = 1; dir >= 0; dir--)
                {
                    int nr = spas[dir].nrecv;
                    for (int v = 0; v < 2; v++)
                    {
                        gmx::RVec* x = (v == 0 ? x0 : x1);
                        for (int i = 0; i < nr; i++)
                        {
                            x[nn + i] = *rbuf;
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
            const gmx_specatsend_t* spas = &spac->spas[d][0];
            /* Copy the required coordinates to the send buffer */
            gmx::RVec* vbuf = spac->vbuf.data();
            for (int v = 0; v < nvec; v++)
            {
                const gmx::RVec* x = (v == 0 ? x0 : x1);
                if (dd->unitCellInfo.haveScrewPBC && dim == XX
                    && (dd->ci[XX] == 0 || dd->ci[XX] == dd->numCells[XX] - 1))
                {
                    /* Here we only perform the rotation, the rest of the pbc
                     * is handled in the constraint or viste routines.
                     */
                    for (int a : spas->a)
                    {
                        (*vbuf)[XX] = x[a][XX];
                        (*vbuf)[YY] = box[YY][YY] - x[a][YY];
                        (*vbuf)[ZZ] = box[ZZ][ZZ] - x[a][ZZ];
                        vbuf++;
                    }
                }
                else
                {
                    for (int a : spas->a)
                    {
                        *vbuf = x[a];
                        vbuf++;
                    }
                }
            }
            /* Send and receive the coordinates */
            if (nvec == 1)
            {
                ddSendrecv(dd,
                           d,
                           dddirBackward,
                           gmx::arrayRefFromArray(spac->vbuf.data(), spas->a.size()),
                           gmx::arrayRefFromArray(x0 + n, spas->nrecv));
            }
            else
            {
                gmx::RVec* vbuf = spac->vbuf.data();
                /* Communicate both vectors in one buffer */
                gmx::RVec* rbuf = spac->vbuf2.data();
                ddSendrecv(dd,
                           d,
                           dddirBackward,
                           gmx::arrayRefFromArray(vbuf, 2 * spas->a.size()),
                           gmx::arrayRefFromArray(rbuf, 2 * spas->nrecv));
                /* Split the buffer into the two vectors */
                int nr = spas[0].nrecv;
                for (int v = 0; v < 2; v++)
                {
                    gmx::RVec* x = (v == 0 ? x0 : x1);
                    for (int i = 0; i < nr; i++)
                    {
                        x[n + i] = *rbuf;
                        rbuf++;
                    }
                }
            }
            n += spas->nrecv;
        }
    }
}

int setup_specat_communication(gmx_domdec_t*             dd,
                               std::vector<int>*         ireq,
                               gmx_domdec_specat_comm_t* spac,
                               gmx::HashedMap<int>*      ga2la_specat,
                               int                       at_start,
                               int                       vbuf_fac,
                               const char*               specat_type,
                               const char*               add_err)
{
    std::array<int, 2> nsend;
    std::array<int, 2> nsendZero = { 0, 0 };
    std::array<int, 2> buf;

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
    int nlast              = nsend[1];
    for (int d = dd->ndim - 1; d >= 0; d--)
    {
        /* Pulse the grid forward and backward */
        int                 dim  = dd->dim[d];
        bool                bPBC = (dim < dd->unitCellInfo.npbcdim);
        const int           ndir = (dd->numCells[dim] == 2)
                                           ?
                                           /* Only 2 cells, so we only need to communicate once */
                                 1
                                           : 2;
        std::array<int, 2>* nsendPtr;
        for (int dir = 0; dir < ndir; dir++)
        {
            if (!bPBC && dd->numCells[dim] > 2
                && ((dir == 0 && dd->ci[dim] == dd->numCells[dim] - 1) || (dir == 1 && dd->ci[dim] == 0)))
            {
                /* No pbc: the fist/last cell should not request atoms */
                nsendPtr = &nsendZero;
            }
            else
            {
                nsendPtr = &nsend;
            }
            /* Communicate the number of indices */
            ddSendrecv(dd,
                       d,
                       dir == 0 ? dddirForward : dddirBackward,
                       gmx::makeArrayRef(*nsendPtr),
                       gmx::arrayRefFromArray(spac->nreq[d][dir], 2));
            int nr = spac->nreq[d][dir][1];
            ireq->resize(nlast + nr);
            /* Communicate the indices */
            ddSendrecv(dd,
                       d,
                       dir == 0 ? dddirForward : dddirBackward,
                       gmx::arrayRefFromArray(ireq->data(), (*nsendPtr)[1]),
                       gmx::arrayRefFromArray(ireq->data() + nlast, nr));
            nlast += nr;
        }
        nsend[1] = nlast;
    }
    if (debug)
    {
        fprintf(debug, "Communicated the counts\n");
    }

    /* Search for the requested atoms and communicate the indices we have */
    int nat_tot_specat = at_start;
    int nrecv_local    = 0;
    for (int d = 0; d < dd->ndim; d++)
    {
        /* Pulse the grid forward and backward */
        const int ndir = (dd->dim[d] >= dd->unitCellInfo.npbcdim || dd->numCells[dd->dim[d]] > 2) ? 2 : 1;
        int       nat_tot_prev = nat_tot_specat;
        for (int dir = ndir - 1; dir >= 0; dir--)
        {
            /* To avoid cost of clearing by resize(), we only increase size */
            if (static_cast<size_t>(nat_tot_specat) > spac->sendAtom.size())
            {
                /* Note: resize initializes new elements to false, which is actually needed here */
                spac->sendAtom.resize(nat_tot_specat);
            }
            gmx_specatsend_t* spas = &spac->spas[d][dir];
            const int         n0   = spac->nreq[d][dir][0];
            const int         nr   = spac->nreq[d][dir][1];
            if (debug)
            {
                fprintf(debug, "dim=%d, dir=%d, searching for %d atoms\n", d, dir, nr);
            }
            const int start = nlast - nr;
            spas->a.clear();
            spac->ibuf.clear();
            nsend[0] = 0;
            for (int i = 0; i < nr; i++)
            {
                const int indr = (*ireq)[start + i];
                int       ind  = 0;
                /* Check if this is a home atom and if so ind will be set */
                if (const int* homeIndex = dd->ga2la->findHome(indr))
                {
                    ind = *homeIndex;
                }
                else
                {
                    /* Search in the communicated atoms */
                    if (const int* a = ga2la_specat->find(indr))
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
            ddSendrecv(dd,
                       d,
                       dir == 0 ? dddirBackward : dddirForward,
                       gmx::makeArrayRef(nsend),
                       gmx::makeArrayRef(buf));
            if (debug)
            {
                fprintf(debug,
                        "Send to rank %d, %d (%d) indices, "
                        "receive from rank %d, %d (%d) indices\n",
                        dd->neighbor[d][1 - dir],
                        nsend[1],
                        nsend[0],
                        dd->neighbor[d][dir],
                        buf[1],
                        buf[0]);
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
            spas->nrecv = buf[1];
            dd->globalAtomIndices.resize(nat_tot_specat + spas->nrecv);
            /* Send and receive the indices */
            ddSendrecv(dd,
                       d,
                       dir == 0 ? dddirBackward : dddirForward,
                       gmx::makeArrayRef(spac->ibuf),
                       gmx::arrayRefFromArray(dd->globalAtomIndices.data() + nat_tot_specat, spas->nrecv));
            nat_tot_specat += spas->nrecv;
        }

        /* Increase the x/f communication buffer sizes, when necessary */
        int ns = spac->spas[d][0].a.size();
        int nr = spac->spas[d][0].nrecv;
        if (ndir == 2)
        {
            ns += spac->spas[d][1].a.size();
            nr += spac->spas[d][1].nrecv;
        }
        if (vbuf_fac * ns > gmx::Index(spac->vbuf.size()))
        {
            spac->vbuf.resize(vbuf_fac * ns);
        }
        if (vbuf_fac == 2 && vbuf_fac * nr > gmx::Index(spac->vbuf2.size()))
        {
            spac->vbuf2.resize(vbuf_fac * nr);
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
            fprintf(debug,
                    "Requested %d, received %d (tot recv %d)\n",
                    numRequested,
                    nrecv_local,
                    nat_tot_specat - at_start);
            if (gmx_debug_at)
            {
                for (int i = 0; i < numRequested; i++)
                {
                    const int* ind = ga2la_specat->find((*ireq)[i]);
                    fprintf(debug, " %s%d", ind ? "" : "!", (*ireq)[i] + 1);
                }
                fprintf(debug, "\n");
            }
        }
        fprintf(stderr,
                "\nDD cell %d %d %d: Neighboring cells do not have atoms:",
                dd->ci[XX],
                dd->ci[YY],
                dd->ci[ZZ]);
        for (int i = 0; i < numRequested; i++)
        {
            if (!ga2la_specat->find((*ireq)[i]))
            {
                fprintf(stderr, " %d", (*ireq)[i] + 1);
            }
        }
        fprintf(stderr, "\n");
        gmx_fatal(FARGS,
                  "DD cell %d %d %d could only obtain %d of the %d atoms that are connected via "
                  "%ss from the neighboring cells. This probably means your %s lengths are too "
                  "long compared to the domain decomposition cell size. Decrease the number of "
                  "domain decomposition grid cells%s%s.",
                  dd->ci[XX],
                  dd->ci[YY],
                  dd->ci[ZZ],
                  nrecv_local,
                  numRequested,
                  specat_type,
                  specat_type,
                  add_err,
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
