/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 The GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
 * computing energies and forces for the PME long-ranged part (Coulomb
 * and LJ).
 *
 * \author Erik Lindahl <erik@kth.se>
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_ewald
 */
/* IMPORTANT FOR DEVELOPERS:
 *
 * Triclinic pme stuff isn't entirely trivial, and we've experienced
 * some bugs during development (many of them due to me). To avoid
 * this in the future, please check the following things if you make
 * changes in this file:
 *
 * 1. You should obtain identical (at least to the PME precision)
 *    energies, forces, and virial for
 *    a rectangular box and a triclinic one where the z (or y) axis is
 *    tilted a whole box side. For instance you could use these boxes:
 *
 *    rectangular       triclinic
 *     2  0  0           2  0  0
 *     0  2  0           0  2  0
 *     0  0  6           2  2  6
 *
 * 2. You should check the energy conservation in a triclinic box.
 *
 * It might seem an overkill, but better safe than sorry.
 * /Erik 001109
 */

#include "gmxpre.h"

#include "pme.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <list>

#include "gromacs/domdec/domdec.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "calculate_spline_moduli.h"
#include "pme_gather.h"
#include "pme_gpu_internal.h"
#include "pme_grid.h"
#include "pme_internal.h"
#include "pme_redistribute.h"
#include "pme_solve.h"
#include "pme_spline_work.h"
#include "pme_spread.h"

/*! \brief Help build a descriptive message in \c error if there are
 * \c errorReasons why PME on GPU is not supported.
 *
 * \returns Whether the lack of errorReasons indicate there is support. */
static bool addMessageIfNotSupported(const std::list<std::string>& errorReasons, std::string* error)
{
    bool isSupported = errorReasons.empty();
    if (!isSupported && error)
    {
        std::string regressionTestMarker = "PME GPU does not support";
        // this prefix is tested for in the regression tests script gmxtest.pl
        *error = regressionTestMarker;
        if (errorReasons.size() == 1)
        {
            *error += " " + errorReasons.back();
        }
        else
        {
            *error += ": " + gmx::joinStrings(errorReasons, "; ");
        }
        *error += ".";
    }
    return isSupported;
}

bool pme_gpu_supports_build(std::string* error)
{
    std::list<std::string> errorReasons;
    if (GMX_DOUBLE)
    {
        errorReasons.emplace_back("a double-precision build");
    }
    if (GMX_GPU == GMX_GPU_NONE)
    {
        errorReasons.emplace_back("a non-GPU build");
    }
    return addMessageIfNotSupported(errorReasons, error);
}

bool pme_gpu_supports_hardware(const gmx_hw_info_t gmx_unused& hwinfo, std::string* error)
{
    std::list<std::string> errorReasons;

    if (GMX_GPU == GMX_GPU_OPENCL)
    {
#ifdef __APPLE__
        errorReasons.emplace_back("Apple OS X operating system");
#endif
    }
    return addMessageIfNotSupported(errorReasons, error);
}

bool pme_gpu_supports_input(const t_inputrec& ir, const gmx_mtop_t& mtop, std::string* error)
{
    std::list<std::string> errorReasons;
    if (!EEL_PME(ir.coulombtype))
    {
        errorReasons.emplace_back("systems that do not use PME for electrostatics");
    }
    if (ir.pme_order != 4)
    {
        errorReasons.emplace_back("interpolation orders other than 4");
    }
    if (ir.efep != efepNO)
    {
        if (gmx_mtop_has_perturbed_charges(mtop))
        {
            errorReasons.emplace_back(
                    "free energy calculations with perturbed charges (multiple grids)");
        }
    }
    if (EVDW_PME(ir.vdwtype))
    {
        errorReasons.emplace_back("Lennard-Jones PME");
    }
    if (!EI_DYNAMICS(ir.eI))
    {
        errorReasons.emplace_back("not a dynamical integrator");
    }
    return addMessageIfNotSupported(errorReasons, error);
}

/*! \brief \libinternal
 * Finds out if PME with given inputs is possible to run on GPU.
 * This function is an internal final check, validating the whole PME structure on creation,
 * but it still duplicates the preliminary checks from the above (externally exposed) pme_gpu_supports_input() - just in case.
 *
 * \param[in]  pme          The PME structure.
 * \param[out] error        The error message if the input is not supported on GPU.
 * \returns                 True if this PME input is possible to run on GPU, false otherwise.
 */
static bool pme_gpu_check_restrictions(const gmx_pme_t* pme, std::string* error)
{
    std::list<std::string> errorReasons;
    if (pme->nnodes != 1)
    {
        errorReasons.emplace_back("PME decomposition");
    }
    if (pme->pme_order != 4)
    {
        errorReasons.emplace_back("interpolation orders other than 4");
    }
    if (pme->bFEP)
    {
        errorReasons.emplace_back("free energy calculations (multiple grids)");
    }
    if (pme->doLJ)
    {
        errorReasons.emplace_back("Lennard-Jones PME");
    }
    if (GMX_DOUBLE)
    {
        errorReasons.emplace_back("double precision");
    }
    if (GMX_GPU == GMX_GPU_NONE)
    {
        errorReasons.emplace_back("non-GPU build of GROMACS");
    }

    return addMessageIfNotSupported(errorReasons, error);
}

PmeRunMode pme_run_mode(const gmx_pme_t* pme)
{
    GMX_ASSERT(pme != nullptr, "Expecting valid PME data pointer");
    return pme->runMode;
}

gmx::PinningPolicy pme_get_pinning_policy()
{
    return gmx::PinningPolicy::PinnedIfSupported;
}

/*! \brief Number of bytes in a cache line.
 *
 * Must also be a multiple of the SIMD and SIMD4 register size, to
 * preserve alignment.
 */
const int gmxCacheLineSize = 64;

//! Set up coordinate communication
static void setup_coordinate_communication(PmeAtomComm* atc)
{
    int nslab, n, i;
    int fw, bw;

    nslab = atc->nslab;

    n = 0;
    for (i = 1; i <= nslab / 2; i++)
    {
        fw = (atc->nodeid + i) % nslab;
        bw = (atc->nodeid - i + nslab) % nslab;
        if (n < nslab - 1)
        {
            atc->slabCommSetup[n].node_dest = fw;
            atc->slabCommSetup[n].node_src  = bw;
            n++;
        }
        if (n < nslab - 1)
        {
            atc->slabCommSetup[n].node_dest = bw;
            atc->slabCommSetup[n].node_src  = fw;
            n++;
        }
    }
}

/*! \brief Round \p n up to the next multiple of \p f */
static int mult_up(int n, int f)
{
    return ((n + f - 1) / f) * f;
}

/*! \brief Return estimate of the load imbalance from the PME grid not being a good match for the number of PME ranks */
static double estimate_pme_load_imbalance(struct gmx_pme_t* pme)
{
    int    nma, nmi;
    double n1, n2, n3;

    nma = pme->nnodes_major;
    nmi = pme->nnodes_minor;

    n1 = mult_up(pme->nkx, nma) * mult_up(pme->nky, nmi) * pme->nkz;
    n2 = mult_up(pme->nkx, nma) * mult_up(pme->nkz, nmi) * pme->nky;
    n3 = mult_up(pme->nky, nma) * mult_up(pme->nkz, nmi) * pme->nkx;

    /* pme_solve is roughly double the cost of an fft */

    return (n1 + n2 + 3 * n3) / static_cast<double>(6 * pme->nkx * pme->nky * pme->nkz);
}

#ifndef DOXYGEN

PmeAtomComm::PmeAtomComm(MPI_Comm   PmeMpiCommunicator,
                         const int  numThreads,
                         const int  pmeOrder,
                         const int  dimIndex,
                         const bool doSpread) :
    dimind(dimIndex),
    bSpread(doSpread),
    pme_order(pmeOrder),
    nthread(numThreads),
    spline(nthread)
{
    if (PmeMpiCommunicator != MPI_COMM_NULL)
    {
        mpi_comm = PmeMpiCommunicator;
#    if GMX_MPI
        MPI_Comm_size(mpi_comm, &nslab);
        MPI_Comm_rank(mpi_comm, &nodeid);
#    endif
    }
    if (debug)
    {
        fprintf(debug, "For PME atom communication in dimind %d: nslab %d rank %d\n", dimind, nslab, nodeid);
    }

    if (nslab > 1)
    {
        slabCommSetup.resize(nslab);
        setup_coordinate_communication(this);

        count_thread.resize(nthread);
        for (auto& countThread : count_thread)
        {
            countThread.resize(nslab);
        }
    }

    if (nthread > 1)
    {
        threadMap.resize(nthread);

#    pragma omp parallel for num_threads(nthread) schedule(static)
        for (int thread = 0; thread < nthread; thread++)
        {
            try
            {
                /* Allocate buffer with padding to avoid cache polution */
                threadMap[thread].nBuffer.resize(nthread + 2 * gmxCacheLineSize);
                threadMap[thread].n = threadMap[thread].nBuffer.data() + gmxCacheLineSize;
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
    }
}

#endif // !DOXYGEN

/*! \brief Initialize data structure for communication */
static void init_overlap_comm(pme_overlap_t* ol, int norder, MPI_Comm comm, int nnodes, int nodeid, int ndata, int commplainsize)
{
    gmx_bool bCont;

    ol->mpi_comm = comm;
    ol->nnodes   = nnodes;
    ol->nodeid   = nodeid;

    /* Linear translation of the PME grid won't affect reciprocal space
     * calculations, so to optimize we only interpolate "upwards",
     * which also means we only have to consider overlap in one direction.
     * I.e., particles on this node might also be spread to grid indices
     * that belong to higher nodes (modulo nnodes)
     */

    ol->s2g0.resize(ol->nnodes + 1);
    ol->s2g1.resize(ol->nnodes);
    if (debug)
    {
        fprintf(debug, "PME slab boundaries:");
    }
    for (int i = 0; i < nnodes; i++)
    {
        /* s2g0 the local interpolation grid start.
         * s2g1 the local interpolation grid end.
         * Since in calc_pidx we divide particles, and not grid lines,
         * spatially uniform along dimension x or y, we need to round
         * s2g0 down and s2g1 up.
         */
        ol->s2g0[i] = (i * ndata + 0) / nnodes;
        ol->s2g1[i] = ((i + 1) * ndata + nnodes - 1) / nnodes + norder - 1;

        if (debug)
        {
            fprintf(debug, "  %3d %3d", ol->s2g0[i], ol->s2g1[i]);
        }
    }
    ol->s2g0[nnodes] = ndata;
    if (debug)
    {
        fprintf(debug, "\n");
    }

    /* Determine with how many nodes we need to communicate the grid overlap */
    int testRankCount = 0;
    do
    {
        testRankCount++;
        bCont = FALSE;
        for (int i = 0; i < nnodes; i++)
        {
            if ((i + testRankCount < nnodes && ol->s2g1[i] > ol->s2g0[i + testRankCount])
                || (i + testRankCount >= nnodes && ol->s2g1[i] > ol->s2g0[i + testRankCount - nnodes] + ndata))
            {
                bCont = TRUE;
            }
        }
    } while (bCont && testRankCount < nnodes);

    ol->comm_data.resize(testRankCount - 1);
    ol->send_size = 0;

    for (size_t b = 0; b < ol->comm_data.size(); b++)
    {
        pme_grid_comm_t* pgc = &ol->comm_data[b];

        /* Send */
        pgc->send_id  = (ol->nodeid + (b + 1)) % ol->nnodes;
        int fft_start = ol->s2g0[pgc->send_id];
        int fft_end   = ol->s2g0[pgc->send_id + 1];
        if (pgc->send_id < nodeid)
        {
            fft_start += ndata;
            fft_end += ndata;
        }
        int send_index1  = ol->s2g1[nodeid];
        send_index1      = std::min(send_index1, fft_end);
        pgc->send_index0 = fft_start;
        pgc->send_nindex = std::max(0, send_index1 - pgc->send_index0);
        ol->send_size += pgc->send_nindex;

        /* We always start receiving to the first index of our slab */
        pgc->recv_id    = (ol->nodeid - (b + 1) + ol->nnodes) % ol->nnodes;
        fft_start       = ol->s2g0[ol->nodeid];
        fft_end         = ol->s2g0[ol->nodeid + 1];
        int recv_index1 = ol->s2g1[pgc->recv_id];
        if (pgc->recv_id > nodeid)
        {
            recv_index1 -= ndata;
        }
        recv_index1      = std::min(recv_index1, fft_end);
        pgc->recv_index0 = fft_start;
        pgc->recv_nindex = std::max(0, recv_index1 - pgc->recv_index0);
    }

#if GMX_MPI
    /* Communicate the buffer sizes to receive */
    MPI_Status stat;
    for (size_t b = 0; b < ol->comm_data.size(); b++)
    {
        MPI_Sendrecv(&ol->send_size, 1, MPI_INT, ol->comm_data[b].send_id, b, &ol->comm_data[b].recv_size,
                     1, MPI_INT, ol->comm_data[b].recv_id, b, ol->mpi_comm, &stat);
    }
#endif

    /* For non-divisible grid we need pme_order iso pme_order-1 */
    ol->sendbuf.resize(norder * commplainsize);
    ol->recvbuf.resize(norder * commplainsize);
}

int minimalPmeGridSize(int pmeOrder)
{
    /* The actual grid size limitations are:
     *   serial:        >= pme_order
     *   DD, no OpenMP: >= 2*(pme_order - 1)
     *   DD, OpenMP:    >= pme_order + 1
     * But we use the maximum for simplicity since in practice there is not
     * much performance difference between pme_order and 2*(pme_order -1).
     */
    int minimalSize = 2 * (pmeOrder - 1);

    GMX_RELEASE_ASSERT(pmeOrder >= 3, "pmeOrder has to be >= 3");
    GMX_RELEASE_ASSERT(minimalSize >= pmeOrder + 1, "The grid size should be >= pmeOrder + 1");

    return minimalSize;
}

bool gmx_pme_check_restrictions(int pme_order, int nkx, int nky, int nkz, int numPmeDomainsAlongX, bool useThreads, bool errorsAreFatal)
{
    if (pme_order > PME_ORDER_MAX)
    {
        if (!errorsAreFatal)
        {
            return false;
        }

        std::string message = gmx::formatString(
                "pme_order (%d) is larger than the maximum allowed value (%d). Modify and "
                "recompile the code if you really need such a high order.",
                pme_order, PME_ORDER_MAX);
        GMX_THROW(gmx::InconsistentInputError(message));
    }

    const int minGridSize = minimalPmeGridSize(pme_order);
    if (nkx < minGridSize || nky < minGridSize || nkz < minGridSize)
    {
        if (!errorsAreFatal)
        {
            return false;
        }
        std::string message =
                gmx::formatString("The PME grid sizes need to be >= 2*(pme_order-1) (%d)", minGridSize);
        GMX_THROW(gmx::InconsistentInputError(message));
    }

    /* Check for a limitation of the (current) sum_fftgrid_dd code.
     * We only allow multiple communication pulses in dim 1, not in dim 0.
     */
    if (useThreads
        && (nkx < numPmeDomainsAlongX * pme_order && nkx != numPmeDomainsAlongX * (pme_order - 1)))
    {
        if (!errorsAreFatal)
        {
            return false;
        }
        gmx_fatal(FARGS,
                  "The number of PME grid lines per rank along x is %g. But when using OpenMP "
                  "threads, the number of grid lines per rank along x should be >= pme_order (%d) "
                  "or = pmeorder-1. To resolve this issue, use fewer ranks along x (and possibly "
                  "more along y and/or z) by specifying -dd manually.",
                  nkx / static_cast<double>(numPmeDomainsAlongX), pme_order);
    }

    return true;
}

/*! \brief Round \p enumerator */
static int div_round_up(int enumerator, int denominator)
{
    return (enumerator + denominator - 1) / denominator;
}

gmx_pme_t* gmx_pme_init(const t_commrec*     cr,
                        const NumPmeDomains& numPmeDomains,
                        const t_inputrec*    ir,
                        gmx_bool             bFreeEnergy_q,
                        gmx_bool             bFreeEnergy_lj,
                        gmx_bool             bReproducible,
                        real                 ewaldcoeff_q,
                        real                 ewaldcoeff_lj,
                        int                  nthread,
                        PmeRunMode           runMode,
                        PmeGpu*              pmeGpu,
                        const DeviceContext* deviceContext,
                        const DeviceStream*  deviceStream,
                        const PmeGpuProgram* pmeGpuProgram,
                        const gmx::MDLogger& /*mdlog*/)
{
    int  use_threads, sum_use_threads, i;
    ivec ndata;

    if (debug)
    {
        fprintf(debug, "Creating PME data structures.\n");
    }

    gmx::unique_cptr<gmx_pme_t, gmx_pme_destroy> pme(new gmx_pme_t());

    pme->sum_qgrid_tmp    = nullptr;
    pme->sum_qgrid_dd_tmp = nullptr;

    pme->buf_nalloc = 0;

    pme->nnodes  = 1;
    pme->bPPnode = TRUE;

    pme->nnodes_major = numPmeDomains.x;
    pme->nnodes_minor = numPmeDomains.y;

    if (numPmeDomains.x * numPmeDomains.y > 1)
    {
        pme->mpi_comm = cr->mpi_comm_mygroup;

#if GMX_MPI
        MPI_Comm_rank(pme->mpi_comm, &pme->nodeid);
        MPI_Comm_size(pme->mpi_comm, &pme->nnodes);
#endif
        if (pme->nnodes != numPmeDomains.x * numPmeDomains.y)
        {
            gmx_incons("PME rank count mismatch");
        }
    }
    else
    {
        pme->mpi_comm = MPI_COMM_NULL;
    }

    if (pme->nnodes == 1)
    {
        pme->mpi_comm_d[0] = MPI_COMM_NULL;
        pme->mpi_comm_d[1] = MPI_COMM_NULL;
        pme->ndecompdim    = 0;
        pme->nodeid_major  = 0;
        pme->nodeid_minor  = 0;
    }
    else
    {
        if (numPmeDomains.y == 1)
        {
            pme->mpi_comm_d[0] = pme->mpi_comm;
            pme->mpi_comm_d[1] = MPI_COMM_NULL;
            pme->ndecompdim    = 1;
            pme->nodeid_major  = pme->nodeid;
            pme->nodeid_minor  = 0;
        }
        else if (numPmeDomains.x == 1)
        {
            pme->mpi_comm_d[0] = MPI_COMM_NULL;
            pme->mpi_comm_d[1] = pme->mpi_comm;
            pme->ndecompdim    = 1;
            pme->nodeid_major  = 0;
            pme->nodeid_minor  = pme->nodeid;
        }
        else
        {
            if (pme->nnodes % numPmeDomains.x != 0)
            {
                gmx_incons(
                        "For 2D PME decomposition, #PME ranks must be divisible by the number of "
                        "domains along x");
            }
            pme->ndecompdim = 2;

#if GMX_MPI
            MPI_Comm_split(pme->mpi_comm, pme->nodeid % numPmeDomains.y, pme->nodeid,
                           &pme->mpi_comm_d[0]); /* My communicator along major dimension */
            MPI_Comm_split(pme->mpi_comm, pme->nodeid / numPmeDomains.y, pme->nodeid,
                           &pme->mpi_comm_d[1]); /* My communicator along minor dimension */

            MPI_Comm_rank(pme->mpi_comm_d[0], &pme->nodeid_major);
            MPI_Comm_size(pme->mpi_comm_d[0], &pme->nnodes_major);
            MPI_Comm_rank(pme->mpi_comm_d[1], &pme->nodeid_minor);
            MPI_Comm_size(pme->mpi_comm_d[1], &pme->nnodes_minor);
#endif
        }
    }
    // cr is always initialized if there is a a PP rank, so we can safely assume
    // that when it is not, like in ewald tests, we not on a PP rank.
    pme->bPPnode = ((cr != nullptr && cr->duty != 0) && thisRankHasDuty(cr, DUTY_PP));

    pme->nthread = nthread;

    /* Check if any of the PME MPI ranks uses threads */
    use_threads = (pme->nthread > 1 ? 1 : 0);
#if GMX_MPI
    if (pme->nnodes > 1)
    {
        MPI_Allreduce(&use_threads, &sum_use_threads, 1, MPI_INT, MPI_SUM, pme->mpi_comm);
    }
    else
#endif
    {
        sum_use_threads = use_threads;
    }
    pme->bUseThreads = (sum_use_threads > 0);

    if (ir->pbcType == PbcType::Screw)
    {
        gmx_fatal(FARGS, "pme does not (yet) work with pbc = screw");
    }

    /* NOTE:
     * It is likely that the current gmx_pme_do() routine supports calculating
     * only Coulomb or LJ while gmx_pme_init() configures for both,
     * but that has never been tested.
     * It is likely that the current gmx_pme_do() routine supports calculating,
     * not calculating free-energy for Coulomb and/or LJ while gmx_pme_init()
     * configures with free-energy, but that has never been tested.
     */
    pme->doCoulomb     = EEL_PME(ir->coulombtype);
    pme->doLJ          = EVDW_PME(ir->vdwtype);
    pme->bFEP_q        = ((ir->efep != efepNO) && bFreeEnergy_q);
    pme->bFEP_lj       = ((ir->efep != efepNO) && bFreeEnergy_lj);
    pme->bFEP          = (pme->bFEP_q || pme->bFEP_lj);
    pme->nkx           = ir->nkx;
    pme->nky           = ir->nky;
    pme->nkz           = ir->nkz;
    pme->bP3M          = (ir->coulombtype == eelP3M_AD || getenv("GMX_PME_P3M") != nullptr);
    pme->pme_order     = ir->pme_order;
    pme->ewaldcoeff_q  = ewaldcoeff_q;
    pme->ewaldcoeff_lj = ewaldcoeff_lj;

    /* Always constant electrostatics coefficients */
    pme->epsilon_r = ir->epsilon_r;

    /* Always constant LJ coefficients */
    pme->ljpme_combination_rule = ir->ljpme_combination_rule;

    // The box requires scaling with nwalls = 2, we store that condition as well
    // as the scaling factor
    delete pme->boxScaler;
    pme->boxScaler = new EwaldBoxZScaler(*ir);

    /* If we violate restrictions, generate a fatal error here */
    gmx_pme_check_restrictions(pme->pme_order, pme->nkx, pme->nky, pme->nkz, pme->nnodes_major,
                               pme->bUseThreads, true);

    if (pme->nnodes > 1)
    {
        double imbal;

#if GMX_MPI
        MPI_Type_contiguous(DIM, GMX_MPI_REAL, &(pme->rvec_mpi));
        MPI_Type_commit(&(pme->rvec_mpi));
#endif

        /* Note that the coefficient spreading and force gathering, which usually
         * takes about the same amount of time as FFT+solve_pme,
         * is always fully load balanced
         * (unless the coefficient distribution is inhomogeneous).
         */

        imbal = estimate_pme_load_imbalance(pme.get());
        if (imbal >= 1.2 && pme->nodeid_major == 0 && pme->nodeid_minor == 0)
        {
            fprintf(stderr,
                    "\n"
                    "NOTE: The load imbalance in PME FFT and solve is %d%%.\n"
                    "      For optimal PME load balancing\n"
                    "      PME grid_x (%d) and grid_y (%d) should be divisible by #PME_ranks_x "
                    "(%d)\n"
                    "      and PME grid_y (%d) and grid_z (%d) should be divisible by #PME_ranks_y "
                    "(%d)\n"
                    "\n",
                    gmx::roundToInt((imbal - 1) * 100), pme->nkx, pme->nky, pme->nnodes_major,
                    pme->nky, pme->nkz, pme->nnodes_minor);
        }
    }

    /* For non-divisible grid we need pme_order iso pme_order-1 */
    /* In sum_qgrid_dd x overlap is copied in place: take padding into account.
     * y is always copied through a buffer: we don't need padding in z,
     * but we do need the overlap in x because of the communication order.
     */
    init_overlap_comm(&pme->overlap[0], pme->pme_order, pme->mpi_comm_d[0], pme->nnodes_major,
                      pme->nodeid_major, pme->nkx,
                      (div_round_up(pme->nky, pme->nnodes_minor) + pme->pme_order)
                              * (pme->nkz + pme->pme_order - 1));

    /* Along overlap dim 1 we can send in multiple pulses in sum_fftgrid_dd.
     * We do this with an offset buffer of equal size, so we need to allocate
     * extra for the offset. That's what the (+1)*pme->nkz is for.
     */
    init_overlap_comm(&pme->overlap[1], pme->pme_order, pme->mpi_comm_d[1], pme->nnodes_minor,
                      pme->nodeid_minor, pme->nky,
                      (div_round_up(pme->nkx, pme->nnodes_major) + pme->pme_order + 1) * pme->nkz);

    /* Double-check for a limitation of the (current) sum_fftgrid_dd code.
     * Note that gmx_pme_check_restrictions checked for this already.
     */
    if (pme->bUseThreads && (pme->overlap[0].comm_data.size() > 1))
    {
        gmx_incons(
                "More than one communication pulse required for grid overlap communication along "
                "the major dimension while using threads");
    }

    snew(pme->bsp_mod[XX], pme->nkx);
    snew(pme->bsp_mod[YY], pme->nky);
    snew(pme->bsp_mod[ZZ], pme->nkz);

    pme->gpu     = pmeGpu; /* Carrying over the single GPU structure */
    pme->runMode = runMode;

    /* The required size of the interpolation grid, including overlap.
     * The allocated size (pmegrid_n?) might be slightly larger.
     */
    pme->pmegrid_nx = pme->overlap[0].s2g1[pme->nodeid_major] - pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_ny = pme->overlap[1].s2g1[pme->nodeid_minor] - pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_nz_base = pme->nkz;
    pme->pmegrid_nz      = pme->pmegrid_nz_base + pme->pme_order - 1;
    set_grid_alignment(&pme->pmegrid_nz, pme->pme_order);
    pme->pmegrid_start_ix = pme->overlap[0].s2g0[pme->nodeid_major];
    pme->pmegrid_start_iy = pme->overlap[1].s2g0[pme->nodeid_minor];
    pme->pmegrid_start_iz = 0;

    make_gridindex_to_localindex(pme->nkx, pme->pmegrid_start_ix,
                                 pme->pmegrid_nx - (pme->pme_order - 1), &pme->nnx, &pme->fshx);
    make_gridindex_to_localindex(pme->nky, pme->pmegrid_start_iy,
                                 pme->pmegrid_ny - (pme->pme_order - 1), &pme->nny, &pme->fshy);
    make_gridindex_to_localindex(pme->nkz, pme->pmegrid_start_iz, pme->pmegrid_nz_base, &pme->nnz,
                                 &pme->fshz);

    pme->spline_work = make_pme_spline_work(pme->pme_order);

    ndata[0] = pme->nkx;
    ndata[1] = pme->nky;
    ndata[2] = pme->nkz;
    /* It doesn't matter if we allocate too many grids here,
     * we only allocate and use the ones we need.
     */
    if (pme->doLJ)
    {
        pme->ngrids = ((ir->ljpme_combination_rule == eljpmeLB) ? DO_Q_AND_LJ_LB : DO_Q_AND_LJ);
    }
    else
    {
        pme->ngrids = DO_Q;
    }
    snew(pme->fftgrid, pme->ngrids);
    snew(pme->cfftgrid, pme->ngrids);
    snew(pme->pfft_setup, pme->ngrids);

    for (i = 0; i < pme->ngrids; ++i)
    {
        if ((i < DO_Q && pme->doCoulomb && (i == 0 || bFreeEnergy_q))
            || (i >= DO_Q && pme->doLJ
                && (i == 2 || bFreeEnergy_lj || ir->ljpme_combination_rule == eljpmeLB)))
        {
            pmegrids_init(&pme->pmegrid[i], pme->pmegrid_nx, pme->pmegrid_ny, pme->pmegrid_nz,
                          pme->pmegrid_nz_base, pme->pme_order, pme->bUseThreads, pme->nthread,
                          pme->overlap[0].s2g1[pme->nodeid_major]
                                  - pme->overlap[0].s2g0[pme->nodeid_major + 1],
                          pme->overlap[1].s2g1[pme->nodeid_minor]
                                  - pme->overlap[1].s2g0[pme->nodeid_minor + 1]);
            /* This routine will allocate the grid data to fit the FFTs */
            const auto allocateRealGridForGpu = (pme->runMode == PmeRunMode::Mixed)
                                                        ? gmx::PinningPolicy::PinnedIfSupported
                                                        : gmx::PinningPolicy::CannotBePinned;
            gmx_parallel_3dfft_init(&pme->pfft_setup[i], ndata, &pme->fftgrid[i], &pme->cfftgrid[i],
                                    pme->mpi_comm_d, bReproducible, pme->nthread, allocateRealGridForGpu);
        }
    }

    if (!pme->bP3M)
    {
        /* Use plain SPME B-spline interpolation */
        make_bspline_moduli(pme->bsp_mod, pme->nkx, pme->nky, pme->nkz, pme->pme_order);
    }
    else
    {
        /* Use the P3M grid-optimized influence function */
        make_p3m_bspline_moduli(pme->bsp_mod, pme->nkx, pme->nky, pme->nkz, pme->pme_order);
    }

    /* Use atc[0] for spreading */
    const int firstDimIndex   = (numPmeDomains.x > 1 ? 0 : 1);
    MPI_Comm  mpiCommFirstDim = (pme->nnodes > 1 ? pme->mpi_comm_d[firstDimIndex] : MPI_COMM_NULL);
    bool      doSpread        = true;
    pme->atc.emplace_back(mpiCommFirstDim, pme->nthread, pme->pme_order, firstDimIndex, doSpread);
    if (pme->ndecompdim >= 2)
    {
        const int secondDimIndex = 1;
        doSpread                 = false;
        pme->atc.emplace_back(pme->mpi_comm_d[1], pme->nthread, pme->pme_order, secondDimIndex, doSpread);
    }

    // Initial check of validity of the input for running on the GPU
    if (pme->runMode != PmeRunMode::CPU)
    {
        std::string errorString;
        bool        canRunOnGpu = pme_gpu_check_restrictions(pme.get(), &errorString);
        if (!canRunOnGpu)
        {
            GMX_THROW(gmx::NotImplementedError(errorString));
        }
        pme_gpu_reinit(pme.get(), deviceContext, deviceStream, pmeGpuProgram);
    }
    else
    {
        GMX_ASSERT(pme->gpu == nullptr, "Should not have PME GPU object when PME is on a CPU.");
    }


    pme_init_all_work(&pme->solve_work, pme->nthread, pme->nkx);

    // no exception was thrown during the init, so we hand over the PME structure handle
    return pme.release();
}

void gmx_pme_reinit(struct gmx_pme_t** pmedata,
                    const t_commrec*   cr,
                    struct gmx_pme_t*  pme_src,
                    const t_inputrec*  ir,
                    const ivec         grid_size,
                    real               ewaldcoeff_q,
                    real               ewaldcoeff_lj)
{
    // Create a copy of t_inputrec fields that are used in gmx_pme_init().
    // TODO: This would be better as just copying a sub-structure that contains
    // all the PME parameters and nothing else.
    t_inputrec irc;
    irc.pbcType                = ir->pbcType;
    irc.coulombtype            = ir->coulombtype;
    irc.vdwtype                = ir->vdwtype;
    irc.efep                   = ir->efep;
    irc.pme_order              = ir->pme_order;
    irc.epsilon_r              = ir->epsilon_r;
    irc.ljpme_combination_rule = ir->ljpme_combination_rule;
    irc.nkx                    = grid_size[XX];
    irc.nky                    = grid_size[YY];
    irc.nkz                    = grid_size[ZZ];

    try
    {
        const gmx::MDLogger dummyLogger;
        // This is reinit which is currently only changing grid size/coefficients,
        // so we don't expect the actual logging.
        // TODO: when PME is an object, it should take reference to mdlog on construction and save it.
        GMX_ASSERT(pmedata, "Invalid PME pointer");
        NumPmeDomains numPmeDomains = { pme_src->nnodes_major, pme_src->nnodes_minor };
        *pmedata = gmx_pme_init(cr, numPmeDomains, &irc, pme_src->bFEP_q, pme_src->bFEP_lj, FALSE,
                                ewaldcoeff_q, ewaldcoeff_lj, pme_src->nthread, pme_src->runMode,
                                pme_src->gpu, nullptr, nullptr, nullptr, dummyLogger);
        /* When running PME on the CPU not using domain decomposition,
         * the atom data is allocated once only in gmx_pme_(re)init().
         */
        if (!pme_src->gpu && pme_src->nnodes == 1)
        {
            gmx_pme_reinit_atoms(*pmedata, pme_src->atc[0].numAtoms(), nullptr);
        }
        // TODO this is mostly passing around current values
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    /* We can easily reuse the allocated pme grids in pme_src */
    reuse_pmegrids(&pme_src->pmegrid[PME_GRID_QA], &(*pmedata)->pmegrid[PME_GRID_QA]);
    /* We would like to reuse the fft grids, but that's harder */
}

void gmx_pme_calc_energy(gmx_pme_t* pme, gmx::ArrayRef<const gmx::RVec> x, gmx::ArrayRef<const real> q, real* V)
{
    pmegrids_t* grid;

    if (pme->nnodes > 1)
    {
        gmx_incons("gmx_pme_calc_energy called in parallel");
    }
    if (pme->bFEP_q)
    {
        gmx_incons("gmx_pme_calc_energy with free energy");
    }

    if (!pme->atc_energy)
    {
        pme->atc_energy = std::make_unique<PmeAtomComm>(MPI_COMM_NULL, 1, pme->pme_order, 0, true);
    }
    PmeAtomComm* atc = pme->atc_energy.get();
    atc->setNumAtoms(x.ssize());
    atc->x           = x;
    atc->coefficient = q;

    /* We only use the A-charges grid */
    grid = &pme->pmegrid[PME_GRID_QA];

    /* Only calculate the spline coefficients, don't actually spread */
    spread_on_grid(pme, atc, nullptr, TRUE, FALSE, pme->fftgrid[PME_GRID_QA], FALSE, PME_GRID_QA);

    *V = gather_energy_bsplines(pme, grid->grid.grid, atc);
}

/*! \brief Calculate initial Lorentz-Berthelot coefficients for LJ-PME */
static void calc_initial_lb_coeffs(gmx::ArrayRef<real> coefficient, const real* local_c6, const real* local_sigma)
{
    for (gmx::index i = 0; i < coefficient.ssize(); ++i)
    {
        real sigma4    = local_sigma[i];
        sigma4         = sigma4 * sigma4;
        sigma4         = sigma4 * sigma4;
        coefficient[i] = local_c6[i] / sigma4;
    }
}

/*! \brief Calculate next Lorentz-Berthelot coefficients for LJ-PME */
static void calc_next_lb_coeffs(gmx::ArrayRef<real> coefficient, const real* local_sigma)
{
    for (gmx::index i = 0; i < coefficient.ssize(); ++i)
    {
        coefficient[i] *= local_sigma[i];
    }
}

int gmx_pme_do(struct gmx_pme_t*              pme,
               gmx::ArrayRef<const gmx::RVec> coordinates,
               gmx::ArrayRef<gmx::RVec>       forces,
               real                           chargeA[],
               real                           chargeB[],
               real                           c6A[],
               real                           c6B[],
               real                           sigmaA[],
               real                           sigmaB[],
               const matrix                   box,
               const t_commrec*               cr,
               int                            maxshift_x,
               int                            maxshift_y,
               t_nrnb*                        nrnb,
               gmx_wallcycle*                 wcycle,
               matrix                         vir_q,
               matrix                         vir_lj,
               real*                          energy_q,
               real*                          energy_lj,
               real                           lambda_q,
               real                           lambda_lj,
               real*                          dvdlambda_q,
               real*                          dvdlambda_lj,
               const gmx::StepWorkload&       stepWork)
{
    GMX_ASSERT(pme->runMode == PmeRunMode::CPU,
               "gmx_pme_do should not be called on the GPU PME run.");

    int                  d, npme, grid_index, max_grid_index;
    PmeAtomComm&         atc         = pme->atc[0];
    pmegrids_t*          pmegrid     = nullptr;
    real*                grid        = nullptr;
    real*                coefficient = nullptr;
    PmeOutput            output[2]; // The second is used for the B state with FEP
    real                 scale, lambda;
    gmx_bool             bClearF;
    gmx_parallel_3dfft_t pfft_setup;
    real*                fftgrid;
    t_complex*           cfftgrid;
    int                  thread;
    gmx_bool             bFirst, bDoSplines;
    int                  fep_state;
    int                  fep_states_lj = pme->bFEP_lj ? 2 : 1;
    // There's no support for computing energy without virial, or vice versa
    const bool computeEnergyAndVirial = (stepWork.computeEnergy || stepWork.computeVirial);

    /* We could be passing lambda!=0 while no q or LJ is actually perturbed */
    if (!pme->bFEP_q)
    {
        lambda_q = 0;
    }
    if (!pme->bFEP_lj)
    {
        lambda_lj = 0;
    }

    assert(pme->nnodes > 0);
    assert(pme->nnodes == 1 || pme->ndecompdim > 0);

    if (pme->nnodes > 1)
    {
        atc.pd.resize(coordinates.ssize());
        for (int d = pme->ndecompdim - 1; d >= 0; d--)
        {
            PmeAtomComm& atc = pme->atc[d];
            atc.maxshift     = (atc.dimind == 0 ? maxshift_x : maxshift_y);
        }
    }
    else
    {
        GMX_ASSERT(coordinates.ssize() == atc.numAtoms(), "We expect atc.numAtoms() coordinates");
        GMX_ASSERT(forces.ssize() >= atc.numAtoms(),
                   "We need a force buffer with at least atc.numAtoms() elements");

        atc.x = coordinates;
        atc.f = forces;
    }

    matrix scaledBox;
    pme->boxScaler->scaleBox(box, scaledBox);

    gmx::invertBoxMatrix(scaledBox, pme->recipbox);
    bFirst = TRUE;

    /* For simplicity, we construct the splines for all particles if
     * more than one PME calculations is needed. Some optimization
     * could be done by keeping track of which atoms have splines
     * constructed, and construct new splines on each pass for atoms
     * that don't yet have them.
     */

    bDoSplines = pme->bFEP || (pme->doCoulomb && pme->doLJ);

    /* We need a maximum of four separate PME calculations:
     * grid_index=0: Coulomb PME with charges from state A
     * grid_index=1: Coulomb PME with charges from state B
     * grid_index=2: LJ PME with C6 from state A
     * grid_index=3: LJ PME with C6 from state B
     * For Lorentz-Berthelot combination rules, a separate loop is used to
     * calculate all the terms
     */

    /* If we are doing LJ-PME with LB, we only do Q here */
    max_grid_index = (pme->ljpme_combination_rule == eljpmeLB) ? DO_Q : DO_Q_AND_LJ;

    for (grid_index = 0; grid_index < max_grid_index; ++grid_index)
    {
        /* Check if we should do calculations at this grid_index
         * If grid_index is odd we should be doing FEP
         * If grid_index < 2 we should be doing electrostatic PME
         * If grid_index >= 2 we should be doing LJ-PME
         */
        if ((grid_index < DO_Q && (!pme->doCoulomb || (grid_index == 1 && !pme->bFEP_q)))
            || (grid_index >= DO_Q && (!pme->doLJ || (grid_index == 3 && !pme->bFEP_lj))))
        {
            continue;
        }
        /* Unpack structure */
        pmegrid    = &pme->pmegrid[grid_index];
        fftgrid    = pme->fftgrid[grid_index];
        cfftgrid   = pme->cfftgrid[grid_index];
        pfft_setup = pme->pfft_setup[grid_index];
        switch (grid_index)
        {
            case 0: coefficient = chargeA; break;
            case 1: coefficient = chargeB; break;
            case 2: coefficient = c6A; break;
            case 3: coefficient = c6B; break;
        }

        grid = pmegrid->grid.grid;

        if (debug)
        {
            fprintf(debug, "PME: number of ranks = %d, rank = %d\n", cr->nnodes, cr->nodeid);
            fprintf(debug, "Grid = %p\n", static_cast<void*>(grid));
            if (grid == nullptr)
            {
                gmx_fatal(FARGS, "No grid!");
            }
        }

        if (pme->nnodes == 1)
        {
            atc.coefficient = gmx::arrayRefFromArray(coefficient, coordinates.size());
        }
        else
        {
            wallcycle_start(wcycle, ewcPME_REDISTXF);
            do_redist_pos_coeffs(pme, cr, bFirst, coordinates, coefficient);

            wallcycle_stop(wcycle, ewcPME_REDISTXF);
        }

        if (debug)
        {
            fprintf(debug, "Rank= %6d, pme local particles=%6d\n", cr->nodeid, atc.numAtoms());
        }

        wallcycle_start(wcycle, ewcPME_SPREAD);

        /* Spread the coefficients on a grid */
        spread_on_grid(pme, &atc, pmegrid, bFirst, TRUE, fftgrid, bDoSplines, grid_index);

        if (bFirst)
        {
            inc_nrnb(nrnb, eNR_WEIGHTS, DIM * atc.numAtoms());
        }
        inc_nrnb(nrnb, eNR_SPREADBSP, pme->pme_order * pme->pme_order * pme->pme_order * atc.numAtoms());

        if (!pme->bUseThreads)
        {
            wrap_periodic_pmegrid(pme, grid);

            /* sum contributions to local grid from other nodes */
            if (pme->nnodes > 1)
            {
                gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_FORWARD);
            }

            copy_pmegrid_to_fftgrid(pme, grid, fftgrid, grid_index);
        }

        wallcycle_stop(wcycle, ewcPME_SPREAD);

        /* TODO If the OpenMP and single-threaded implementations
           converge, then spread_on_grid() and
           copy_pmegrid_to_fftgrid() will perhaps live in the same
           source file.
        */

        /* Here we start a large thread parallel region */
#pragma omp parallel num_threads(pme->nthread) private(thread)
        {
            try
            {
                thread = gmx_omp_get_thread_num();
                int loop_count;

                /* do 3d-fft */
                if (thread == 0)
                {
                    wallcycle_start(wcycle, ewcPME_FFT);
                }
                gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_REAL_TO_COMPLEX, thread, wcycle);
                if (thread == 0)
                {
                    wallcycle_stop(wcycle, ewcPME_FFT);
                }

                /* solve in k-space for our local cells */
                if (thread == 0)
                {
                    wallcycle_start(wcycle, (grid_index < DO_Q ? ewcPME_SOLVE : ewcLJPME));
                }
                if (grid_index < DO_Q)
                {
                    loop_count = solve_pme_yzx(
                            pme, cfftgrid, scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ],
                            computeEnergyAndVirial, pme->nthread, thread);
                }
                else
                {
                    loop_count =
                            solve_pme_lj_yzx(pme, &cfftgrid, FALSE,
                                             scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ],
                                             computeEnergyAndVirial, pme->nthread, thread);
                }

                if (thread == 0)
                {
                    wallcycle_stop(wcycle, (grid_index < DO_Q ? ewcPME_SOLVE : ewcLJPME));
                    inc_nrnb(nrnb, eNR_SOLVEPME, loop_count);
                }

                /* do 3d-invfft */
                if (thread == 0)
                {
                    wallcycle_start(wcycle, ewcPME_FFT);
                }
                gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_COMPLEX_TO_REAL, thread, wcycle);
                if (thread == 0)
                {
                    wallcycle_stop(wcycle, ewcPME_FFT);


                    if (pme->nodeid == 0)
                    {
                        real ntot = pme->nkx * pme->nky * pme->nkz;
                        npme      = static_cast<int>(ntot * std::log(ntot) / std::log(2.0));
                        inc_nrnb(nrnb, eNR_FFT, 2 * npme);
                    }

                    /* Note: this wallcycle region is closed below
                       outside an OpenMP region, so take care if
                       refactoring code here. */
                    wallcycle_start(wcycle, ewcPME_GATHER);
                }

                copy_fftgrid_to_pmegrid(pme, fftgrid, grid, grid_index, pme->nthread, thread);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }
        /* End of thread parallel section.
         * With MPI we have to synchronize here before gmx_sum_qgrid_dd.
         */

        /* distribute local grid to all nodes */
        if (pme->nnodes > 1)
        {
            gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_BACKWARD);
        }

        unwrap_periodic_pmegrid(pme, grid);

        if (stepWork.computeForces)
        {
            /* interpolate forces for our local atoms */


            /* If we are running without parallelization,
             * atc->f is the actual force array, not a buffer,
             * therefore we should not clear it.
             */
            lambda  = grid_index < DO_Q ? lambda_q : lambda_lj;
            bClearF = (bFirst && PAR(cr));
#pragma omp parallel for num_threads(pme->nthread) schedule(static)
            for (thread = 0; thread < pme->nthread; thread++)
            {
                try
                {
                    gather_f_bsplines(pme, grid, bClearF, &atc, &atc.spline[thread],
                                      pme->bFEP ? (grid_index % 2 == 0 ? 1.0 - lambda : lambda) : 1.0);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }


            inc_nrnb(nrnb, eNR_GATHERFBSP,
                     pme->pme_order * pme->pme_order * pme->pme_order * atc.numAtoms());
            /* Note: this wallcycle region is opened above inside an OpenMP
               region, so take care if refactoring code here. */
            wallcycle_stop(wcycle, ewcPME_GATHER);
        }

        if (computeEnergyAndVirial)
        {
            /* This should only be called on the master thread
             * and after the threads have synchronized.
             */
            if (grid_index < 2)
            {
                get_pme_ener_vir_q(pme->solve_work, pme->nthread, &output[grid_index % 2]);
            }
            else
            {
                get_pme_ener_vir_lj(pme->solve_work, pme->nthread, &output[grid_index % 2]);
            }
        }
        bFirst = FALSE;
    } /* of grid_index-loop */

    /* For Lorentz-Berthelot combination rules in LJ-PME, we need to calculate
     * seven terms. */

    if (pme->doLJ && pme->ljpme_combination_rule == eljpmeLB)
    {
        /* Loop over A- and B-state if we are doing FEP */
        for (fep_state = 0; fep_state < fep_states_lj; ++fep_state)
        {
            real *local_c6 = nullptr, *local_sigma = nullptr, *RedistC6 = nullptr, *RedistSigma = nullptr;
            gmx::ArrayRef<real> coefficientBuffer;
            if (pme->nnodes == 1)
            {
                pme->lb_buf1.resize(atc.numAtoms());
                coefficientBuffer = pme->lb_buf1;
                switch (fep_state)
                {
                    case 0:
                        local_c6    = c6A;
                        local_sigma = sigmaA;
                        break;
                    case 1:
                        local_c6    = c6B;
                        local_sigma = sigmaB;
                        break;
                    default: gmx_incons("Trying to access wrong FEP-state in LJ-PME routine");
                }
            }
            else
            {
                coefficientBuffer = atc.coefficientBuffer;
                switch (fep_state)
                {
                    case 0:
                        RedistC6    = c6A;
                        RedistSigma = sigmaA;
                        break;
                    case 1:
                        RedistC6    = c6B;
                        RedistSigma = sigmaB;
                        break;
                    default: gmx_incons("Trying to access wrong FEP-state in LJ-PME routine");
                }
                wallcycle_start(wcycle, ewcPME_REDISTXF);

                do_redist_pos_coeffs(pme, cr, bFirst, coordinates, RedistC6);
                pme->lb_buf1.resize(atc.numAtoms());
                pme->lb_buf2.resize(atc.numAtoms());
                local_c6 = pme->lb_buf1.data();
                for (int i = 0; i < atc.numAtoms(); ++i)
                {
                    local_c6[i] = atc.coefficient[i];
                }

                do_redist_pos_coeffs(pme, cr, FALSE, coordinates, RedistSigma);
                local_sigma = pme->lb_buf2.data();
                for (int i = 0; i < atc.numAtoms(); ++i)
                {
                    local_sigma[i] = atc.coefficient[i];
                }

                wallcycle_stop(wcycle, ewcPME_REDISTXF);
            }
            atc.coefficient = coefficientBuffer;
            calc_initial_lb_coeffs(coefficientBuffer, local_c6, local_sigma);

            /*Seven terms in LJ-PME with LB, grid_index < 2 reserved for electrostatics*/
            for (grid_index = 2; grid_index < 9; ++grid_index)
            {
                /* Unpack structure */
                pmegrid    = &pme->pmegrid[grid_index];
                fftgrid    = pme->fftgrid[grid_index];
                pfft_setup = pme->pfft_setup[grid_index];
                calc_next_lb_coeffs(coefficientBuffer, local_sigma);
                grid = pmegrid->grid.grid;

                wallcycle_start(wcycle, ewcPME_SPREAD);
                /* Spread the c6 on a grid */
                spread_on_grid(pme, &atc, pmegrid, bFirst, TRUE, fftgrid, bDoSplines, grid_index);

                if (bFirst)
                {
                    inc_nrnb(nrnb, eNR_WEIGHTS, DIM * atc.numAtoms());
                }

                inc_nrnb(nrnb, eNR_SPREADBSP,
                         pme->pme_order * pme->pme_order * pme->pme_order * atc.numAtoms());
                if (pme->nthread == 1)
                {
                    wrap_periodic_pmegrid(pme, grid);
                    /* sum contributions to local grid from other nodes */
                    if (pme->nnodes > 1)
                    {
                        gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_FORWARD);
                    }
                    copy_pmegrid_to_fftgrid(pme, grid, fftgrid, grid_index);
                }
                wallcycle_stop(wcycle, ewcPME_SPREAD);

                /*Here we start a large thread parallel region*/
#pragma omp parallel num_threads(pme->nthread) private(thread)
                {
                    try
                    {
                        thread = gmx_omp_get_thread_num();
                        /* do 3d-fft */
                        if (thread == 0)
                        {
                            wallcycle_start(wcycle, ewcPME_FFT);
                        }

                        gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_REAL_TO_COMPLEX, thread, wcycle);
                        if (thread == 0)
                        {
                            wallcycle_stop(wcycle, ewcPME_FFT);
                        }
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                }
                bFirst = FALSE;
            }
            /* solve in k-space for our local cells */
#pragma omp parallel num_threads(pme->nthread) private(thread)
            {
                try
                {
                    int loop_count;
                    thread = gmx_omp_get_thread_num();
                    if (thread == 0)
                    {
                        wallcycle_start(wcycle, ewcLJPME);
                    }

                    loop_count =
                            solve_pme_lj_yzx(pme, &pme->cfftgrid[2], TRUE,
                                             scaledBox[XX][XX] * scaledBox[YY][YY] * scaledBox[ZZ][ZZ],
                                             computeEnergyAndVirial, pme->nthread, thread);
                    if (thread == 0)
                    {
                        wallcycle_stop(wcycle, ewcLJPME);
                        inc_nrnb(nrnb, eNR_SOLVEPME, loop_count);
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
            }

            if (computeEnergyAndVirial)
            {
                /* This should only be called on the master thread and
                 * after the threads have synchronized.
                 */
                get_pme_ener_vir_lj(pme->solve_work, pme->nthread, &output[fep_state]);
            }

            bFirst = !pme->doCoulomb;
            calc_initial_lb_coeffs(coefficientBuffer, local_c6, local_sigma);
            for (grid_index = 8; grid_index >= 2; --grid_index)
            {
                /* Unpack structure */
                pmegrid    = &pme->pmegrid[grid_index];
                fftgrid    = pme->fftgrid[grid_index];
                pfft_setup = pme->pfft_setup[grid_index];
                grid       = pmegrid->grid.grid;
                calc_next_lb_coeffs(coefficientBuffer, local_sigma);
#pragma omp parallel num_threads(pme->nthread) private(thread)
                {
                    try
                    {
                        thread = gmx_omp_get_thread_num();
                        /* do 3d-invfft */
                        if (thread == 0)
                        {
                            wallcycle_start(wcycle, ewcPME_FFT);
                        }

                        gmx_parallel_3dfft_execute(pfft_setup, GMX_FFT_COMPLEX_TO_REAL, thread, wcycle);
                        if (thread == 0)
                        {
                            wallcycle_stop(wcycle, ewcPME_FFT);


                            if (pme->nodeid == 0)
                            {
                                real ntot = pme->nkx * pme->nky * pme->nkz;
                                npme      = static_cast<int>(ntot * std::log(ntot) / std::log(2.0));
                                inc_nrnb(nrnb, eNR_FFT, 2 * npme);
                            }
                            wallcycle_start(wcycle, ewcPME_GATHER);
                        }

                        copy_fftgrid_to_pmegrid(pme, fftgrid, grid, grid_index, pme->nthread, thread);
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                } /*#pragma omp parallel*/

                /* distribute local grid to all nodes */
                if (pme->nnodes > 1)
                {
                    gmx_sum_qgrid_dd(pme, grid, GMX_SUM_GRID_BACKWARD);
                }

                unwrap_periodic_pmegrid(pme, grid);

                if (stepWork.computeForces)
                {
                    /* interpolate forces for our local atoms */
                    bClearF = (bFirst && PAR(cr));
                    scale   = pme->bFEP ? (fep_state < 1 ? 1.0 - lambda_lj : lambda_lj) : 1.0;
                    scale *= lb_scale_factor[grid_index - 2];

#pragma omp parallel for num_threads(pme->nthread) schedule(static)
                    for (thread = 0; thread < pme->nthread; thread++)
                    {
                        try
                        {
                            gather_f_bsplines(pme, grid, bClearF, &pme->atc[0],
                                              &pme->atc[0].spline[thread], scale);
                        }
                        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
                    }


                    inc_nrnb(nrnb, eNR_GATHERFBSP,
                             pme->pme_order * pme->pme_order * pme->pme_order * pme->atc[0].numAtoms());
                }
                wallcycle_stop(wcycle, ewcPME_GATHER);

                bFirst = FALSE;
            } /* for (grid_index = 8; grid_index >= 2; --grid_index) */
        }     /* for (fep_state = 0; fep_state < fep_states_lj; ++fep_state) */
    }         /* if (pme->doLJ && pme->ljpme_combination_rule == eljpmeLB) */

    if (stepWork.computeForces && pme->nnodes > 1)
    {
        wallcycle_start(wcycle, ewcPME_REDISTXF);
        for (d = 0; d < pme->ndecompdim; d++)
        {
            gmx::ArrayRef<gmx::RVec> forcesRef;
            if (d == pme->ndecompdim - 1)
            {
                const size_t numAtoms = coordinates.size();
                GMX_ASSERT(forces.size() >= numAtoms, "Need at least numAtoms forces");
                forcesRef = forces.subArray(0, numAtoms);
            }
            else
            {
                forcesRef = pme->atc[d + 1].f;
            }
            if (DOMAINDECOMP(cr))
            {
                dd_pmeredist_f(pme, &pme->atc[d], forcesRef, d == pme->ndecompdim - 1 && pme->bPPnode);
            }
        }

        wallcycle_stop(wcycle, ewcPME_REDISTXF);
    }

    if (computeEnergyAndVirial)
    {
        if (pme->doCoulomb)
        {
            if (!pme->bFEP_q)
            {
                *energy_q = output[0].coulombEnergy_;
                m_add(vir_q, output[0].coulombVirial_, vir_q);
            }
            else
            {
                *energy_q = (1.0 - lambda_q) * output[0].coulombEnergy_ + lambda_q * output[1].coulombEnergy_;
                *dvdlambda_q += output[1].coulombEnergy_ - output[0].coulombEnergy_;
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        vir_q[i][j] += (1.0 - lambda_q) * output[0].coulombVirial_[i][j]
                                       + lambda_q * output[1].coulombVirial_[i][j];
                    }
                }
            }
            if (debug)
            {
                fprintf(debug, "Electrostatic PME mesh energy: %g\n", *energy_q);
            }
        }
        else
        {
            *energy_q = 0;
        }

        if (pme->doLJ)
        {
            if (!pme->bFEP_lj)
            {
                *energy_lj = output[0].lennardJonesEnergy_;
                m_add(vir_lj, output[0].lennardJonesVirial_, vir_lj);
            }
            else
            {
                *energy_lj = (1.0 - lambda_lj) * output[0].lennardJonesEnergy_
                             + lambda_lj * output[1].lennardJonesEnergy_;
                *dvdlambda_lj += output[1].lennardJonesEnergy_ - output[0].lennardJonesEnergy_;
                for (int i = 0; i < DIM; i++)
                {
                    for (int j = 0; j < DIM; j++)
                    {
                        vir_lj[i][j] += (1.0 - lambda_lj) * output[0].lennardJonesVirial_[i][j]
                                        + lambda_lj * output[1].lennardJonesVirial_[i][j];
                    }
                }
            }
            if (debug)
            {
                fprintf(debug, "Lennard-Jones PME mesh energy: %g\n", *energy_lj);
            }
        }
        else
        {
            *energy_lj = 0;
        }
    }
    return 0;
}

void gmx_pme_destroy(gmx_pme_t* pme)
{
    if (!pme)
    {
        return;
    }

    delete pme->boxScaler;

    sfree(pme->nnx);
    sfree(pme->nny);
    sfree(pme->nnz);
    sfree(pme->fshx);
    sfree(pme->fshy);
    sfree(pme->fshz);

    for (int i = 0; i < pme->ngrids; ++i)
    {
        pmegrids_destroy(&pme->pmegrid[i]);
    }
    if (pme->pfft_setup)
    {
        for (int i = 0; i < pme->ngrids; ++i)
        {
            gmx_parallel_3dfft_destroy(pme->pfft_setup[i]);
        }
    }
    sfree(pme->fftgrid);
    sfree(pme->cfftgrid);
    sfree(pme->pfft_setup);

    for (int i = 0; i < DIM; i++)
    {
        sfree(pme->bsp_mod[i]);
    }

    sfree(pme->bufv);
    sfree(pme->bufr);

    if (pme->solve_work)
    {
        pme_free_all_work(&pme->solve_work, pme->nthread);
    }

    sfree(pme->sum_qgrid_tmp);
    sfree(pme->sum_qgrid_dd_tmp);

    destroy_pme_spline_work(pme->spline_work);

    if (pme->gpu != nullptr)
    {
        pme_gpu_destroy(pme->gpu);
    }

    delete pme;
}

void gmx_pme_reinit_atoms(gmx_pme_t* pme, const int numAtoms, const real* charges)
{
    if (pme->gpu != nullptr)
    {
        pme_gpu_reinit_atoms(pme->gpu, numAtoms, charges);
    }
    else
    {
        pme->atc[0].setNumAtoms(numAtoms);
        // TODO: set the charges here as well
    }
}

bool gmx_pme_grid_matches(const gmx_pme_t& pme, const ivec grid_size)
{
    return (pme.nkx == grid_size[XX] && pme.nky == grid_size[YY] && pme.nkz == grid_size[ZZ]);
}
