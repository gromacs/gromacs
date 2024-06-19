/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * \brief This file contains function definitions for redistributing
 * atoms over the PME domains
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */

#include "gmxpre.h"

#include "pme_redistribute.h"

#include "config.h"

#include <cstdio>

#include <algorithm>
#include <array>
#include <filesystem>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxmpi.h"

#include "pme_internal.h"

//! Calculate the slab indices and store in \p atc, store counts in \p count
static void pme_calc_pidx(int                            start,
                          int                            end,
                          const matrix                   recipbox,
                          gmx::ArrayRef<const gmx::RVec> x,
                          PmeAtomComm*                   atc,
                          int*                           count)
{
    int         nslab, i;
    int         si;
    const real* xptr;
    real        s;
    real        rxx, ryx, rzx, ryy, rzy;
    int*        pd;

    /* Calculate PME task index (pidx) for each grid index.
     * Here we always assign equally sized slabs to each rank
     * for load balancing reasons (the PME grid spacing is not used).
     */

    nslab = atc->nslab;
    pd    = atc->pd.data();

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
        /* Calculate the slab index in x-dimension */
        for (i = start; i < end; i++)
        {
            xptr = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab * (xptr[XX] * rxx + xptr[YY] * ryx + xptr[ZZ] * rzx);
            si    = static_cast<int>(s + 2 * nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
    else
    {
        ryy = recipbox[YY][YY];
        rzy = recipbox[ZZ][YY];
        /* Calculate the slab index in y-dimension */
        for (i = start; i < end; i++)
        {
            xptr = x[i];
            /* Fractional coordinates along box vectors */
            s     = nslab * (xptr[YY] * ryy + xptr[ZZ] * rzy);
            si    = static_cast<int>(s + 2 * nslab) % nslab;
            pd[i] = si;
            count[si]++;
        }
    }
}

//! Wrapper function for calculating slab indices, stored in \p atc
static void pme_calc_pidx_wrapper(gmx::ArrayRef<const gmx::RVec> x, const matrix recipbox, PmeAtomComm* atc)
{
    int nthread = atc->nthread;

#pragma omp parallel for num_threads(nthread) schedule(static)
    for (int thread = 0; thread < nthread; thread++)
    {
        try
        {
            const int natoms = x.ssize();
            pme_calc_pidx(natoms * thread / nthread,
                          natoms * (thread + 1) / nthread,
                          recipbox,
                          x,
                          atc,
                          atc->count_thread[thread].data());
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    /* Non-parallel reduction, since nslab is small */

    for (int thread = 1; thread < nthread; thread++)
    {
        for (int slab = 0; slab < atc->nslab; slab++)
        {
            atc->count_thread[0][slab] += atc->count_thread[thread][slab];
        }
    }
}

#ifndef DOXYGEN

void SplineCoefficients::realloc(const int nalloc)
{
    const int padding = 4;

    bufferX_.resize(nalloc);
    coefficients[XX] = bufferX_.data();
    bufferY_.resize(nalloc);
    coefficients[YY] = bufferY_.data();
    /* In z we add padding, this is only required for the aligned 4-wide SIMD code */
    bufferZ_.resize(nalloc + 2 * padding);
    coefficients[ZZ] = bufferZ_.data() + padding;
}

#endif // !DOXYGEN

//! Reallocates all buffers in \p spline to fit atoms in \p atc
static void pme_realloc_splinedata(splinedata_t* spline, const PmeAtomComm* atc)
{
    if (spline->nalloc >= atc->x.ssize() && spline->nalloc >= atc->numAtoms())
    {
        return;
    }

    spline->nalloc = std::max(atc->x.capacity(), static_cast<size_t>(atc->numAtoms()));
    spline->ind.resize(spline->nalloc);
    /* Initialize the index to identity so it works without threads */
    for (int i = 0; i < spline->nalloc; i++)
    {
        spline->ind[i] = i;
    }

    spline->theta.realloc(atc->pme_order * spline->nalloc);
    spline->dtheta.realloc(atc->pme_order * spline->nalloc);
}

#ifndef DOXYGEN

void PmeAtomComm::setNumAtoms(const int numAtoms)
{
    numAtoms_ = numAtoms;

    if (nslab > 1)
    {
        /* We have to avoid a NULL pointer for atc->x to avoid
         * possible fatal errors in MPI routines.
         */
        xBuffer.resize(numAtoms_);
        if (xBuffer.capacity() == 0)
        {
            xBuffer.reserve(1);
        }
        x = xBuffer;
        coefficientBuffer.resize(numAtoms_);
        if (coefficientBuffer.capacity() == 0)
        {
            coefficientBuffer.reserve(1);
        }
        coefficient          = coefficientBuffer;
        const int nalloc_old = fBuffer.size();
        fBuffer.resize(numAtoms_);
        for (int i = nalloc_old; i < numAtoms_; i++)
        {
            clear_rvec(fBuffer[i]);
        }
        f = fBuffer;
    }
    if (bSpread)
    {
        fractx.resize(numAtoms_);
        idx.resize(numAtoms_);

        if (nthread > 1)
        {
            thread_idx.resize(numAtoms_);
        }

        for (int i = 0; i < nthread; i++)
        {
            pme_realloc_splinedata(&spline[i], this);
        }
    }
}

#endif // !DOXYGEN

//! Communicates buffers between rank separated by \p shift slabs
static void pme_dd_sendrecv(PmeAtomComm gmx_unused* atc,
                            gmx_bool gmx_unused     bBackward,
                            int gmx_unused          shift,
                            void gmx_unused* buf_s,
                            int gmx_unused   nbyte_s,
                            void gmx_unused* buf_r,
                            int gmx_unused   nbyte_r)
{
#if GMX_MPI
    int        dest, src;
    MPI_Status stat;

    if (!bBackward)
    {
        dest = atc->slabCommSetup[shift].node_dest;
        src  = atc->slabCommSetup[shift].node_src;
    }
    else
    {
        dest = atc->slabCommSetup[shift].node_src;
        src  = atc->slabCommSetup[shift].node_dest;
    }

    if (nbyte_s > 0 && nbyte_r > 0)
    {
        MPI_Sendrecv(
                buf_s, nbyte_s, MPI_BYTE, dest, shift, buf_r, nbyte_r, MPI_BYTE, src, shift, atc->mpi_comm, &stat);
    }
    else if (nbyte_s > 0)
    {
        MPI_Send(buf_s, nbyte_s, MPI_BYTE, dest, shift, atc->mpi_comm);
    }
    else if (nbyte_r > 0)
    {
        MPI_Recv(buf_r, nbyte_r, MPI_BYTE, src, shift, atc->mpi_comm, &stat);
    }
#endif
}

//! Redistristributes \p data and optionally coordinates between MPI ranks
static void dd_pmeredist_pos_coeffs(gmx_pme_t*                     pme,
                                    const gmx_bool                 bX,
                                    gmx::ArrayRef<const gmx::RVec> x,
                                    gmx::ArrayRef<const real>      data,
                                    PmeAtomComm*                   atc)
{
    int nnodes_comm, i, local_pos, buf_pos;

    nnodes_comm = std::min(2 * atc->maxshift, atc->nslab - 1);

    auto sendCount = atc->sendCount();
    int  nsend     = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        const int commnode           = atc->slabCommSetup[i].node_dest;
        atc->bufferIndices[commnode] = nsend;
        nsend += sendCount[commnode];
    }
    if (bX)
    {
        if (sendCount[atc->slabIndex] + nsend != x.ssize())
        {
            gmx_fatal(
                    FARGS,
                    "%zd particles communicated to PME rank %d are more than 2/3 times the cut-off "
                    "out of the domain decomposition cell of their charge group in dimension %c.\n"
                    "This usually means that your system is not well equilibrated.",
                    x.ssize() - (sendCount[atc->slabIndex] + nsend),
                    pme->nodeid,
                    'x' + atc->dimind);
        }

        if (nsend > gmx::ssize(pme->bufr) || nsend > gmx::ssize(pme->bufv))
        {
            pme->bufv.resize(nsend);
            pme->bufr.resize(nsend);
        }

        int numAtoms = sendCount[atc->slabIndex];
        for (i = 0; i < nnodes_comm; i++)
        {
            const int commnode = atc->slabCommSetup[i].node_dest;
            int       scount   = sendCount[commnode];
            /* Communicate the count */
            if (debug)
            {
                fprintf(debug,
                        "dimind %d PME rank %d send to rank %d: %d\n",
                        atc->dimind,
                        atc->slabIndex,
                        commnode,
                        scount);
            }
            pme_dd_sendrecv(
                    atc, FALSE, i, &scount, sizeof(int), &atc->slabCommSetup[i].rcount, sizeof(int));
            numAtoms += atc->slabCommSetup[i].rcount;
        }

        atc->setNumAtoms(numAtoms);
    }

    local_pos = 0;
    for (gmx::Index i = 0; i < x.ssize(); i++)
    {
        const int slabIndex = atc->pd[i];
        if (slabIndex == atc->slabIndex)
        {
            /* Copy direct to the receive buffer */
            if (bX)
            {
                atc->xBuffer[local_pos] = x[i];
            }
            atc->coefficientBuffer[local_pos] = data[i];
            local_pos++;
        }
        else
        {
            /* Copy to the send buffer */
            int& buf_index = atc->bufferIndices[slabIndex];
            if (bX)
            {
                pme->bufv[buf_index] = x[i];
            }
            pme->bufr[buf_index] = data[i];
            buf_index++;
        }
    }

    buf_pos = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        const int scount = atc->sendCount()[atc->slabCommSetup[i].node_dest];
        const int rcount = atc->slabCommSetup[i].rcount;
        if (scount > 0 || rcount > 0)
        {
            if (bX)
            {
                /* Communicate the coordinates */
                pme_dd_sendrecv(atc,
                                FALSE,
                                i,
                                pme->bufv.data() + buf_pos,
                                scount * sizeof(rvec),
                                atc->xBuffer.data() + local_pos,
                                rcount * sizeof(rvec));
            }
            /* Communicate the coefficients */
            pme_dd_sendrecv(atc,
                            FALSE,
                            i,
                            pme->bufr.data() + buf_pos,
                            scount * sizeof(real),
                            atc->coefficientBuffer.data() + local_pos,
                            rcount * sizeof(real));
            buf_pos += scount;
            local_pos += atc->slabCommSetup[i].rcount;
        }
    }
    GMX_ASSERT(local_pos == atc->numAtoms(), "After receiving we should have numAtoms coordinates");
}

void dd_pmeredist_f(struct gmx_pme_t* pme, PmeAtomComm* atc, gmx::ArrayRef<gmx::RVec> f, gmx_bool bAddF)
{
    int nnodes_comm, local_pos, buf_pos, i;

    nnodes_comm = std::min(2 * atc->maxshift, atc->nslab - 1);

    local_pos = atc->sendCount()[atc->slabIndex];
    buf_pos   = 0;
    for (i = 0; i < nnodes_comm; i++)
    {
        const int commnode = atc->slabCommSetup[i].node_dest;
        const int scount   = atc->slabCommSetup[i].rcount;
        const int rcount   = atc->sendCount()[commnode];
        if (scount > 0 || rcount > 0)
        {
            /* Communicate the forces */
            pme_dd_sendrecv(atc,
                            TRUE,
                            i,
                            atc->f.data() + local_pos,
                            scount * sizeof(rvec),
                            pme->bufv.data() + buf_pos,
                            rcount * sizeof(rvec));
            local_pos += scount;
        }
        atc->bufferIndices[commnode] = buf_pos;
        buf_pos += rcount;
    }

    local_pos = 0;
    if (bAddF)
    {
        for (gmx::Index i = 0; i < f.ssize(); i++)
        {
            const int slabIndex = atc->pd[i];
            if (slabIndex == atc->slabIndex)
            {
                /* Add from the local force array */
                f[i] += atc->f[local_pos];
                local_pos++;
            }
            else
            {
                /* Add from the receive buffer */
                f[i] += pme->bufv[atc->bufferIndices[slabIndex]];
                atc->bufferIndices[slabIndex]++;
            }
        }
    }
    else
    {
        for (gmx::Index i = 0; i < f.ssize(); i++)
        {
            const int slabIndex = atc->pd[i];
            if (slabIndex == atc->slabIndex)
            {
                /* Copy from the local force array */
                f[i] = atc->f[local_pos];
                local_pos++;
            }
            else
            {
                /* Copy from the receive buffer */
                f[i] = pme->bufv[atc->bufferIndices[slabIndex]];
                atc->bufferIndices[slabIndex]++;
            }
        }
    }
}

void do_redist_pos_coeffs(struct gmx_pme_t*              pme,
                          const t_commrec*               cr,
                          gmx_bool                       bFirst,
                          gmx::ArrayRef<const gmx::RVec> x,
                          gmx::ArrayRef<const real>      data)
{
    for (int d = pme->ndecompdim - 1; d >= 0; d--)
    {
        gmx::ArrayRef<const gmx::RVec> xRef;
        gmx::ArrayRef<const real>      param_d;
        if (d == pme->ndecompdim - 1)
        {
            /* Start out with the local coordinates and charges */
            xRef    = x;
            param_d = data;
        }
        else
        {
            /* Redistribute the data collected along the previous dimension */
            const PmeAtomComm& atc = pme->atc[d + 1];
            xRef                   = atc.x;
            param_d                = atc.coefficient;
        }
        PmeAtomComm& atc = pme->atc[d];
        atc.pd.resize(xRef.size());
        pme_calc_pidx_wrapper(xRef, pme->recipbox, &atc);
        /* Redistribute x (only once) and qA/c6A or qB/c6B */
        if (haveDDAtomOrdering(*cr))
        {
            dd_pmeredist_pos_coeffs(pme, bFirst, xRef, param_d, &atc);
        }
    }
}
