/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <stdio.h>

#include <cstring>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

#include "pme-internal.h"
#include "pme-pp-communication.h"

/*! \brief Block to wait for communication to PME ranks to complete
 *
 * This should be faster with a real non-blocking MPI implementation
 */
static constexpr bool c_useDelayedWait = false;

/*! \brief Wait for the pending data send requests to PME ranks to complete */
static void gmx_pme_send_coeffs_coords_wait(gmx_domdec_t *dd)
{
    if (dd->nreq_pme)
    {
#if GMX_MPI
        MPI_Waitall(dd->nreq_pme, dd->req_pme, MPI_STATUSES_IGNORE);
#endif
        dd->nreq_pme = 0;
    }
}

/*! \brief Send data to PME ranks */
static void gmx_pme_send_coeffs_coords(t_commrec *cr, unsigned int flags,
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
                (flags & PP_PME_CHARGE) ? " charges" : "",
                (flags & PP_PME_SQRTC6) ? " sqrtC6" : "",
                (flags & PP_PME_SIGMA)  ? " sigma" : "",
                (flags & PP_PME_COORD)  ? " coordinates" : "");
    }

    if (c_useDelayedWait)
    {
        /* We can not use cnb until pending communication has finished */
        gmx_pme_send_coeffs_coords_wait(dd);
    }

    if (dd->pme_receive_vir_ener)
    {
        /* Peer PP node: communicate all data */
        if (dd->cnb == nullptr)
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
#if GMX_MPI
        MPI_Isend(cnb, sizeof(*cnb), MPI_BYTE,
                  dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }
    else if (flags & (PP_PME_CHARGE | PP_PME_SQRTC6 | PP_PME_SIGMA))
    {
#if GMX_MPI
        /* Communicate only the number of atoms */
        MPI_Isend(&n, sizeof(n), MPI_BYTE,
                  dd->pme_nodeid, eCommType_CNB, cr->mpi_comm_mysim,
                  &dd->req_pme[dd->nreq_pme++]);
#endif
    }

#if GMX_MPI
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
#endif
    if (!c_useDelayedWait)
    {
        /* Wait for the data to arrive */
        /* We can skip this wait as we are sure x and q will not be modified
         * before the next call to gmx_pme_send_x_q or gmx_pme_receive_f.
         */
        gmx_pme_send_coeffs_coords_wait(dd);
    }
}

void gmx_pme_send_parameters(t_commrec *cr,
                             const interaction_const_t *ic,
                             gmx_bool bFreeEnergy_q, gmx_bool bFreeEnergy_lj,
                             real *chargeA, real *chargeB,
                             real *sqrt_c6A, real *sqrt_c6B,
                             real *sigmaA, real *sigmaB,
                             int maxshift_x, int maxshift_y)
{
    unsigned int flags = 0;

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
                               nullptr, nullptr, 0, 0, maxshift_x, maxshift_y, -1);
}

void gmx_pme_send_coordinates(t_commrec *cr, matrix box, rvec *x,
                              real lambda_q, real lambda_lj,
                              gmx_bool bEnerVir,
                              gmx_int64_t step)
{
    unsigned int flags = PP_PME_COORD;
    if (bEnerVir)
    {
        flags |= PP_PME_ENER_VIR;
    }
    gmx_pme_send_coeffs_coords(cr, flags, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                               box, x, lambda_q, lambda_lj, 0, 0, step);
}

void gmx_pme_send_finish(t_commrec *cr)
{
    unsigned int flags = PP_PME_FINISH;

    gmx_pme_send_coeffs_coords(cr, flags, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 0, 0, 0, 0, -1);
}

void gmx_pme_send_switchgrid(t_commrec gmx_unused *cr,
                             ivec gmx_unused       grid_size,
                             real gmx_unused       ewaldcoeff_q,
                             real gmx_unused       ewaldcoeff_lj)
{
#if GMX_MPI
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
#if GMX_MPI
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

/*! \brief Receive virial and energy from PME rank */
static void receive_virial_energy(t_commrec *cr,
                                  gmx::ForceWithVirial *forceWithVirial,
                                  real *energy_q, real *energy_lj,
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
#if GMX_MPI
        MPI_Recv(&cve, sizeof(cve), MPI_BYTE, cr->dd->pme_nodeid, 1, cr->mpi_comm_mysim,
                 MPI_STATUS_IGNORE);
#else
        memset(&cve, 0, sizeof(cve));
#endif

        forceWithVirial->addVirialContribution(cve.vir_q);
        forceWithVirial->addVirialContribution(cve.vir_lj);
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
                       gmx::ForceWithVirial *forceWithVirial,
                       real *energy_q, real *energy_lj,
                       real *dvdlambda_q, real *dvdlambda_lj,
                       float *pme_cycles)
{
    if (c_useDelayedWait)
    {
        /* Wait for the x request to finish */
        gmx_pme_send_coeffs_coords_wait(cr->dd);
    }

    int natoms = cr->dd->nat_home;

    if (natoms > cr->dd->pme_recv_f_alloc)
    {
        cr->dd->pme_recv_f_alloc = over_alloc_dd(natoms);
        srenew(cr->dd->pme_recv_f_buf, cr->dd->pme_recv_f_alloc);
    }

#if GMX_MPI
    MPI_Recv(cr->dd->pme_recv_f_buf[0],
             natoms*sizeof(rvec), MPI_BYTE,
             cr->dd->pme_nodeid, 0, cr->mpi_comm_mysim,
             MPI_STATUS_IGNORE);
#endif

    int nt = gmx_omp_nthreads_get_simple_rvec_task(emntDefault, natoms);

    gmx::ArrayRef<gmx::RVec> f = forceWithVirial->force_;

    /* Note that we would like to avoid this conditional by putting it
     * into the omp pragma instead, but then we still take the full
     * omp parallel for overhead (at least with gcc5).
     */
    if (nt == 1)
    {
        for (int i = 0; i < natoms; i++)
        {
            rvec_inc(f[i], cr->dd->pme_recv_f_buf[i]);
        }
    }
    else
    {
#pragma omp parallel for num_threads(nt) schedule(static)
        for (int i = 0; i < natoms; i++)
        {
            rvec_inc(f[i], cr->dd->pme_recv_f_buf[i]);
        }
    }

    receive_virial_energy(cr, forceWithVirial, energy_q, energy_lj, dvdlambda_q, dvdlambda_lj, pme_cycles);
}
