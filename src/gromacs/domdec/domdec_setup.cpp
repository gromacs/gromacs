/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2008- The GROMACS Authors
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
 * \brief This file defines functions used by the domdec module
 * in its initial setup phase.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "domdec_setup.h"

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdio>

#include <numeric>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/options.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/perf_est.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "atomdistribution.h"
#include "box.h"
#include "domdec_internal.h"
#include "utility.h"

// TODO remove this when moving domdec into gmx namespace
using gmx::DomdecOptions;

/*! \brief Margin for setting up the DD grid */
#define DD_GRID_MARGIN_PRES_SCALE 1.05

/*! \brief Factorize \p n.
 *
 * \param[in]    n     Value to factorize
 * \param[out]   fac   Vector of factors (to be allocated in this function)
 * \param[out]   mfac  Vector with the number of times each factor repeats in the factorization (to be allocated in this function)
 */
static void factorize(int n, std::vector<int>* fac, std::vector<int>* mfac)
{
    if (n <= 0)
    {
        gmx_fatal(FARGS, "Can only factorize positive integers.");
    }

    /* Decompose n in factors */
    fac->clear();
    mfac->clear();
    int d = 2;
    while (n > 1)
    {
        while (n % d == 0)
        {
            if (fac->empty() || fac->back() != d)
            {
                fac->push_back(d);
                mfac->push_back(1);
            }
            else
            {
                mfac->back()++;
            }
            n /= d;
        }
        d++;
    }
}

/*! \brief Find largest divisor of \p n smaller than \p n*/
static int largest_divisor(int n)
{
    std::vector<int> div;
    std::vector<int> mdiv;
    factorize(n, &div, &mdiv);

    return div.back();
}

/*! \brief Returns TRUE when there are enough PME ranks for the ratio */
static gmx_bool fits_pme_ratio(int nrank_tot, int nrank_pme, float ratio)
{
    return (static_cast<double>(nrank_pme) / static_cast<double>(nrank_tot) > 0.95 * ratio);
}

/*! \brief Returns TRUE when npme out of ntot ranks doing PME is expected to give reasonable performance */
static gmx_bool fits_pp_pme_perf(int ntot, int npme, float ratio)
{
    std::vector<int> div;
    std::vector<int> mdiv;
    factorize(ntot - npme, &div, &mdiv);

    int npp_root3  = gmx::roundToInt(std::cbrt(ntot - npme));
    int npme_root2 = gmx::roundToInt(std::sqrt(static_cast<double>(npme)));

    /* The check below gives a reasonable division:
     * factor 5 allowed at 5 or more PP ranks,
     * factor 7 allowed at 49 or more PP ranks.
     */
    if (div.back() > 3 + npp_root3)
    {
        return FALSE;
    }

    /* Check if the number of PP and PME ranks have a reasonable sized
     * denominator in common, such that we can use 2D PME decomposition
     * when required (which requires nx_pp == nx_pme).
     * The factor of 2 allows for a maximum ratio of 2^2=4
     * between nx_pme and ny_pme.
     */
    if (std::gcd(ntot - npme, npme) * 2 < npme_root2)
    {
        return FALSE;
    }

    /* Does this division gives a reasonable PME load? */
    return fits_pme_ratio(ntot, npme, ratio);
}

/*! \brief Make a guess for the number of PME ranks to use. */
static int guess_npme(const gmx::MDLogger& mdlog,
                      const gmx_mtop_t&    mtop,
                      const t_inputrec&    ir,
                      const matrix         box,
                      int                  nrank_tot)
{
    float ratio = pme_load_estimate(mtop, ir, box);

    GMX_LOG(mdlog.info).appendTextFormatted("Guess for relative PME load: %.2f", ratio);

    /* We assume the optimal rank ratio is close to the load ratio.
     * The communication load is neglected,
     * but (hopefully) this will balance out between PP and PME.
     */

    if (!fits_pme_ratio(nrank_tot, nrank_tot / 2, ratio))
    {
        /* We would need more than nrank_tot/2 PME only nodes,
         * which is not possible. Since the PME load is very high,
         * we will not loose much performance when all ranks do PME.
         */

        return 0;
    }

    /* First try to find npme as a factor of nrank_tot up to nrank_tot/3.
     * We start with a minimum PME node fraction of 1/16
     * and avoid ratios which lead to large prime factors in nnodes-npme.
     */
    int npme = (nrank_tot + 15) / 16;
    while (npme <= nrank_tot / 3)
    {
        if (nrank_tot % npme == 0)
        {
            /* Note that fits_perf might change the PME grid,
             * in the current implementation it does not.
             */
            if (fits_pp_pme_perf(nrank_tot, npme, ratio))
            {
                break;
            }
        }
        npme++;
    }
    if (npme > nrank_tot / 3)
    {
        /* Try any possible number for npme */
        npme = 1;
        while (npme <= nrank_tot / 2)
        {
            /* Note that fits_perf may change the PME grid */
            if (fits_pp_pme_perf(nrank_tot, npme, ratio))
            {
                break;
            }
            npme++;
        }
    }
    if (npme > nrank_tot / 2)
    {
        gmx_fatal(FARGS,
                  "Could not find an appropriate number of separate PME ranks. i.e. >= %5f*#ranks "
                  "(%d) and <= #ranks/2 (%d) and reasonable performance wise (grid_x=%d, "
                  "grid_y=%d).\n"
                  "Use the -npme option of mdrun or change the number of ranks or the PME grid "
                  "dimensions, see the manual for details.",
                  ratio,
                  gmx::roundToInt(0.95 * ratio * nrank_tot),
                  nrank_tot / 2,
                  ir.nkx,
                  ir.nky);
    }
    else
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Will use %d particle-particle and %d PME only ranks\n"
                        "This is a guess, check the performance at the end of the log file",
                        nrank_tot - npme,
                        npme);
    }

    return npme;
}

/*! \brief Return \p n divided by \p f rounded up to the next integer. */
static int div_up(int n, int f)
{
    return (n + f - 1) / f;
}

real comm_box_frac(const gmx::IVec& dd_nc, real cutoff, const gmx_ddbox_t& ddbox)
{
    rvec nw;

    for (int i = 0; i < DIM; i++)
    {
        real bt = ddbox.box_size[i] * ddbox.skew_fac[i];
        nw[i]   = dd_nc[i] * cutoff / bt;
    }

    real comm_vol = 0;
    for (int i = 0; i < DIM; i++)
    {
        if (dd_nc[i] > 1)
        {
            comm_vol += nw[i];
            for (int j = i + 1; j < DIM; j++)
            {
                if (dd_nc[j] > 1)
                {
                    comm_vol += nw[i] * nw[j] * M_PI / 4;
                    for (int k = j + 1; k < DIM; k++)
                    {
                        if (dd_nc[k] > 1)
                        {
                            comm_vol += nw[i] * nw[j] * nw[k] * M_PI / 6;
                        }
                    }
                }
            }
        }
    }

    return comm_vol;
}

/*! \brief Return whether the DD inhomogeneous in the z direction */
static gmx_bool inhomogeneous_z(const t_inputrec& ir)
{
    return ((usingPme(ir.coulombtype) || ir.coulombtype == CoulombInteractionType::Ewald)
            && ir.pbcType == PbcType::Xyz && ir.ewald_geometry == EwaldGeometry::ThreeDC);
}

/*! \brief Estimate cost of PME FFT communication
 *
 * This only takes the communication into account and not imbalance
 * in the calculation. But the imbalance in communication and calculation
 * are similar and therefore these formulas also prefer load balance
 * in the FFT and pme_solve calculation.
 */
static float comm_pme_cost_vol(int npme, int a, int b, int c)
{
    /* We use a float here, since an integer might overflow */
    float comm_vol = npme - 1;
    comm_vol *= npme;
    comm_vol *= div_up(a, npme);
    comm_vol *= div_up(b, npme);
    comm_vol *= c;

    return comm_vol;
}

/*! \brief Estimate cost of communication for a possible domain decomposition. */
static float comm_cost_est(real               limit,
                           real               cutoff,
                           const matrix       box,
                           const gmx_ddbox_t& ddbox,
                           const int64_t      natoms,
                           const t_inputrec&  ir,
                           float              pbcdxr,
                           int                npme_tot,
                           const gmx::IVec&   nc)
{
    gmx::IVec npme = { 1, 1, 1 };
    rvec      bt;
    /* This is the cost of a pbc_dx call relative to the cost
     * of communicating the coordinate and force of an atom.
     * This will be machine dependent.
     * These factors are for x86 with SMP or Infiniband.
     */
    float pbcdx_rect_fac = 0.1;
    float pbcdx_tric_fac = 0.2;

    /* Check the DD algorithm restrictions */
    if ((ir.pbcType == PbcType::XY && ir.nwall < 2 && nc[ZZ] > 1)
        || (ir.pbcType == PbcType::Screw && (nc[XX] == 1 || nc[YY] > 1 || nc[ZZ] > 1)))
    {
        return -1;
    }

    if (inhomogeneous_z(ir) && nc[ZZ] > 1)
    {
        return -1;
    }

    assert(ddbox.npbcdim <= DIM);

    /* Check if the triclinic requirements are met */
    for (int i = 0; i < DIM; i++)
    {
        for (int j = i + 1; j < ddbox.npbcdim; j++)
        {
            if (box[j][i] != 0 || ir.deform[j][i] != 0
                || (ir.pressureCouplingOptions.epc != PressureCoupling::No
                    && ir.pressureCouplingOptions.compress[j][i] != 0))
            {
                if (nc[j] > 1 && nc[i] == 1)
                {
                    return -1;
                }
            }
        }
    }

    for (int i = 0; i < DIM; i++)
    {
        bt[i] = ddbox.box_size[i] * ddbox.skew_fac[i];

        /* Without PBC and with 2 cells, there are no lower limits on the cell size */
        if (!(i >= ddbox.npbcdim && nc[i] <= 2) && bt[i] < nc[i] * limit)
        {
            return -1;
        }
        /* With PBC, check if the cut-off fits in nc[i]-1 cells */
        if (i < ddbox.npbcdim && nc[i] > 1 && (nc[i] - 1) * bt[i] < nc[i] * cutoff)
        {
            return -1;
        }
    }

    if (npme_tot > 1)
    {
        /* The following choices should match those
         * in init_domain_decomposition in domdec.c.
         */
        if (nc[XX] == 1 && nc[YY] > 1)
        {
            npme[XX] = 1;
            npme[YY] = npme_tot;
        }
        else if (nc[YY] == 1)
        {
            npme[XX] = npme_tot;
            npme[YY] = 1;
        }
        else
        {
            /* Will we use 1D or 2D PME decomposition? */
            npme[XX] = (npme_tot % nc[XX] == 0) ? nc[XX] : npme_tot;
            npme[YY] = npme_tot / npme[XX];
        }
    }

    if (usingPme(ir.coulombtype) || usingLJPme(ir.vdwtype))
    {
        /* Check the PME grid restrictions.
         * Currently these can only be invalid here with too few grid lines
         * along the x dimension per rank doing PME.
         */
        int npme_x = (npme_tot > 1 ? npme[XX] : nc[XX]);
        int npme_y = (npme_tot > 1 ? npme[YY] : nc[YY]);

        /* Currently we don't have the OpenMP thread count available here.
         * But with threads we have only tighter restrictions and it's
         * probably better anyhow to avoid settings where we need to reduce
         * grid lines over multiple ranks, as the thread check will do.
         *
         * extendedHaloRegion (used for PME GPU decomposition runs) is also not known
         * at this point. Just ignore it at this point.
         */
        bool useThreads         = true;
        int  extendedHaloRegion = 0;
        bool useGpuPme          = false;
        bool errorsAreFatal     = false;
        if (!gmx_pme_check_restrictions(
                    ir.pme_order, ir.nkx, ir.nky, ir.nkz, npme_x, npme_y, extendedHaloRegion, useGpuPme, useThreads, errorsAreFatal))
        {
            return -1;
        }
    }

    /* When two dimensions are (nearly) equal, use more cells
     * for the smallest index, so the decomposition does not
     * depend sensitively on the rounding of the box elements.
     */
    for (int i = 0; i < DIM; i++)
    {
        for (int j = i + 1; j < DIM; j++)
        {
            /* Check if the box size is nearly identical,
             * in that case we prefer nx > ny  and ny > nz.
             */
            if (std::fabs(bt[j] - bt[i]) < 0.01 * bt[i] && nc[j] > nc[i])
            {
                /* The XX/YY check is a bit compact. If nc[YY]==npme[YY]
                 * this means the swapped nc has nc[XX]==npme[XX],
                 * and we can also swap X and Y for PME.
                 */
                /* Check if dimension i and j are equivalent for PME.
                 * For x/y: if nc[YY]!=npme[YY], we can not swap x/y
                 * For y/z: we can not have PME decomposition in z
                 */
                if (npme_tot <= 1
                    || !((i == XX && j == YY && nc[YY] != npme[YY]) || (i == YY && j == ZZ && npme[YY] > 1)))
                {
                    return -1;
                }
            }
        }
    }

    /* This function determines only half of the communication cost.
     * All PP, PME and PP-PME communication is symmetric
     * and the "back"-communication cost is identical to the forward cost.
     */

    float comm_vol = comm_box_frac(nc, cutoff, ddbox);

    float comm_pme = 0;
    for (int i = 0; i < 2; i++)
    {
        /* Determine the largest volume for PME x/f redistribution */
        if (nc[i] % npme[i] != 0)
        {
            float comm_vol_xf =
                    (nc[i] > npme[i]) ? (npme[i] == 2 ? 1.0 / 3.0 : 0.5)
                                      : (1.0 - std::gcd(nc[i], npme[i]) / static_cast<double>(npme[i]));
            comm_pme += 3 * natoms * comm_vol_xf;
        }

        /* Grid overlap communication */
        if (npme[i] > 1)
        {
            const int nk      = (i == 0 ? ir.nkx : ir.nky);
            const int overlap = (nk % npme[i] == 0 ? ir.pme_order - 1 : ir.pme_order);
            float     temp    = npme[i];
            temp *= overlap;
            temp *= ir.nkx;
            temp *= ir.nky;
            temp *= ir.nkz;
            temp /= nk;
            comm_pme += temp;
            /* Old line comm_pme += npme[i]*overlap*ir.nkx*ir.nky*ir.nkz/nk; */
        }
    }

    comm_pme += comm_pme_cost_vol(npme[YY], ir.nky, ir.nkz, ir.nkx);
    comm_pme += comm_pme_cost_vol(npme[XX], ir.nkx, ir.nky, ir.nkz);

    /* Add cost of pbc_dx for bondeds */
    float cost_pbcdx = 0;
    if ((nc[XX] == 1 || nc[YY] == 1) || (nc[ZZ] == 1 && ir.pbcType != PbcType::XY))
    {
        if ((ddbox.tric_dir[XX] && nc[XX] == 1) || (ddbox.tric_dir[YY] && nc[YY] == 1))
        {
            cost_pbcdx = pbcdxr * pbcdx_tric_fac;
        }
        else
        {
            cost_pbcdx = pbcdxr * pbcdx_rect_fac;
        }
    }

    if (debug)
    {
        fprintf(debug,
                "nc %2d %2d %2d %2d %2d vol pp %6.4f pbcdx %6.4f pme %9.3e tot %9.3e\n",
                nc[XX],
                nc[YY],
                nc[ZZ],
                npme[XX],
                npme[YY],
                comm_vol,
                cost_pbcdx,
                comm_pme / (3 * natoms),
                comm_vol + cost_pbcdx + comm_pme / (3 * natoms));
    }

    return 3 * natoms * (comm_vol + cost_pbcdx) + comm_pme;
}

/*! \brief Assign penalty factors to possible domain decompositions,
 * based on the estimated communication costs. */
static void assign_factors(const real         limit,
                           const real         cutoff,
                           const matrix       box,
                           const gmx_ddbox_t& ddbox,
                           int                natoms,
                           const t_inputrec&  ir,
                           float              pbcdxr,
                           int                npme,
                           int                ndiv,
                           const int*         div,
                           const int*         mdiv,
                           gmx::IVec*         irTryPtr,
                           gmx::IVec*         opt)
{
    gmx::IVec& ir_try = *irTryPtr;

    if (ndiv == 0)
    {

        const float ce = comm_cost_est(limit, cutoff, box, ddbox, natoms, ir, pbcdxr, npme, ir_try);
        if (ce >= 0
            && ((*opt)[XX] == 0
                || ce < comm_cost_est(limit, cutoff, box, ddbox, natoms, ir, pbcdxr, npme, *opt)))
        {
            *opt = ir_try;
        }

        return;
    }

    for (int x = mdiv[0]; x >= 0; x--)
    {
        for (int i = 0; i < x; i++)
        {
            ir_try[XX] *= div[0];
        }
        for (int y = mdiv[0] - x; y >= 0; y--)
        {
            for (int i = 0; i < y; i++)
            {
                ir_try[YY] *= div[0];
            }
            for (int i = 0; i < mdiv[0] - x - y; i++)
            {
                ir_try[ZZ] *= div[0];
            }

            /* recurse */
            assign_factors(
                    limit, cutoff, box, ddbox, natoms, ir, pbcdxr, npme, ndiv - 1, div + 1, mdiv + 1, irTryPtr, opt);

            for (int i = 0; i < mdiv[0] - x - y; i++)
            {
                ir_try[ZZ] /= div[0];
            }
            for (int i = 0; i < y; i++)
            {
                ir_try[YY] /= div[0];
            }
        }
        for (int i = 0; i < x; i++)
        {
            ir_try[XX] /= div[0];
        }
    }
}

/*! \brief Determine the optimal distribution of DD cells for the
 * simulation system and number of MPI ranks
 *
 * \returns The optimal grid cell choice. The latter will contain all
 *          zeros if no valid cell choice exists. */
static gmx::IVec optimizeDDCells(const gmx::MDLogger& mdlog,
                                 const int            numRanksRequested,
                                 const int            numPmeOnlyRanks,
                                 const real           cellSizeLimit,
                                 const gmx_mtop_t&    mtop,
                                 const matrix         box,
                                 const gmx_ddbox_t&   ddbox,
                                 const t_inputrec&    ir,
                                 const DDSystemInfo&  systemInfo)
{
    double pbcdxr = 0;

    const int numPPRanks = numRanksRequested - numPmeOnlyRanks;

    GMX_LOG(mdlog.info)
            .appendTextFormatted(
                    "Optimizing the DD grid for %d cells with a minimum initial size of %.3f nm",
                    numPPRanks,
                    cellSizeLimit);
    if (inhomogeneous_z(ir))
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Ewald_geometry=%s: assuming inhomogeneous particle distribution in z, "
                        "will not decompose in z.",
                        enumValueToString(ir.ewald_geometry));
    }


    // For cost estimates, we need the number of ranks doing PME work,
    // which is the number of PP ranks when not using separate
    // PME-only ranks.
    const int numRanksDoingPmeWork =
            (usingPme(ir.coulombtype) ? ((numPmeOnlyRanks > 0) ? numPmeOnlyRanks : numPPRanks) : 0);

    if (systemInfo.haveInterDomainBondeds)
    {
        /* If we can skip PBC for distance calculations in plain-C bondeds,
         * we can save some time (e.g. 3D DD with pbc=xyz).
         * Here we ignore SIMD bondeds as they always do (fast) PBC.
         */
        count_bonded_distances(mtop, ir, &pbcdxr, nullptr);
        pbcdxr /= static_cast<double>(mtop.natoms);
    }
    else
    {
        /* Every molecule is a single charge group: no pbc required */
        pbcdxr = 0;
    }

    if (cellSizeLimit > 0)
    {
        std::string maximumCells = "The maximum allowed number of cells is:";
        for (int d = 0; d < DIM; d++)
        {
            int nmax = static_cast<int>(ddbox.box_size[d] * ddbox.skew_fac[d] / cellSizeLimit);
            if (d >= ddbox.npbcdim && nmax < 2)
            {
                nmax = 2;
            }
            if (d == ZZ && inhomogeneous_z(ir))
            {
                nmax = 1;
            }
            maximumCells += gmx::formatString(" %c %d", 'X' + d, nmax);
        }
        GMX_LOG(mdlog.info).appendText(maximumCells);
    }

    if (debug)
    {
        fprintf(debug, "Average nr of pbc_dx calls per atom %.2f\n", pbcdxr);
    }

    /* Decompose numPPRanks in factors */
    std::vector<int> div;
    std::vector<int> mdiv;
    factorize(numPPRanks, &div, &mdiv);

    gmx::IVec itry       = { 1, 1, 1 };
    gmx::IVec numDomains = { 0, 0, 0 };
    assign_factors(cellSizeLimit,
                   systemInfo.cutoff,
                   box,
                   ddbox,
                   mtop.natoms,
                   ir,
                   pbcdxr,
                   numRanksDoingPmeWork,
                   div.size(),
                   div.data(),
                   mdiv.data(),
                   &itry,
                   &numDomains);

    return numDomains;
}

real getDDGridSetupCellSizeLimit(const gmx::MDLogger& mdlog,
                                 const bool           bDynLoadBal,
                                 const real           dlb_scale,
                                 const t_inputrec&    ir,
                                 real                 systemInfoCellSizeLimit,
                                 int                  numRanksRequested)
{
    // do not impose DLB-related limits with single-domain
    if (numRanksRequested == 1)
    {
        return systemInfoCellSizeLimit;
    }

    real cellSizeLimit = systemInfoCellSizeLimit;

    /* Add a margin for DLB and/or pressure scaling */
    if (bDynLoadBal)
    {
        if (dlb_scale >= 1.0)
        {
            gmx_fatal(FARGS, "The value for option -dds should be smaller than 1");
        }
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Scaling the initial minimum size with 1/%g (option -dds) = %g", dlb_scale, 1 / dlb_scale);
        cellSizeLimit /= dlb_scale;
    }
    else if (ir.pressureCouplingOptions.epc != PressureCoupling::No)
    {
        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "To account for pressure scaling, scaling the initial minimum size with %g",
                        DD_GRID_MARGIN_PRES_SCALE);
        cellSizeLimit *= DD_GRID_MARGIN_PRES_SCALE;
    }

    return cellSizeLimit;
}

gmx::SeparatePmeRanksPermitted checkForSeparatePmeRanks(const gmx::MDModulesNotifiers& notifiers,
                                                        const DomdecOptions&           options,
                                                        int  numRanksRequested,
                                                        bool useGpuForNonbonded,
                                                        bool useGpuForPme,
                                                        bool canUseGpuPmeDecomposition)
{
    gmx::SeparatePmeRanksPermitted separatePmeRanksPermitted;

    /* Permit MDModules to notify whether they want to use PME-only ranks */
    notifiers.simulationSetupNotifier_.notify(&separatePmeRanksPermitted);

    /* With NB GPUs we don't automatically use PME-only ranks. PME-only CPU ranks can
     * improve performance with many ranks per GPU, since our OpenMP
     * scaling is bad, but it's difficult to automate the setup.
     */
    if (useGpuForNonbonded && options.numPmeRanks < 0)
    {
        separatePmeRanksPermitted.disablePmeRanks(
                "PME-only ranks are not automatically used when "
                "non-bonded interactions are computed on GPUs");
    }

    /* If more than one PME ranks requested, check if PME decomposition is supported */
    if (useGpuForPme && !canUseGpuPmeDecomposition && (options.numPmeRanks < 0 || options.numPmeRanks > 1))
    {
        separatePmeRanksPermitted.disablePmeRanks(
                "PME GPU decomposition is not supported for current build configuration, only one "
                "separate PME-only GPU rank "
                "can be used");
    }

    // With explicit DD grid setup we do not automatically use PME-only ranks
    if (options.numCells[XX] > 0 && options.numPmeRanks < 0)
    {
        separatePmeRanksPermitted.disablePmeRanks("explicit DD grid requested");
    }

    // Controls the automated choice of when to use separate PME-only ranks.
    const int minRankCountToDefaultToSeparatePmeRanks = 19;

    // Check if we have enough total ranks to use PME-only ranks
    if (numRanksRequested < minRankCountToDefaultToSeparatePmeRanks && options.numPmeRanks < 0)
    {
        separatePmeRanksPermitted.disablePmeRanks(
                "there are too few total ranks for efficient splitting");
    }

    return separatePmeRanksPermitted;
}

void checkForValidRankCountRequests(const int                             numRanksRequested,
                                    const bool                            usingPme,
                                    const int                             numPmeRanksRequested,
                                    const gmx::SeparatePmeRanksPermitted& separatePmeRanksPermitted,
                                    const bool                            checkForLargePrimeFactors)
{
    int numPPRanksRequested = numRanksRequested;
    if (usingPme && numPmeRanksRequested > 0)
    {
        numPPRanksRequested -= numPmeRanksRequested;
        if (numPmeRanksRequested > numPPRanksRequested)
        {
            gmx_fatal(FARGS,
                      "Cannot have %d separate PME ranks with only %d PP ranks, choose fewer or no "
                      "separate PME ranks",
                      numPmeRanksRequested,
                      numPPRanksRequested);
        }
    }

    /* If simulation is not using PME but -npme > 0 then inform user and throw an error */
    if (!usingPme && numPmeRanksRequested > 0)
    {
        gmx_fatal(FARGS,
                  "The system does not use PME for electrostatics or LJ. Requested "
                  "-npme %d option is not viable.",
                  numPmeRanksRequested);
    }

    /* If some parts of the code could not use PME-only ranks and
     * user explicitly used mdrun -npme > 0 option then throw an error */
    if (!separatePmeRanksPermitted.permitSeparatePmeRanks() && numPmeRanksRequested > 0)
    {
        gmx_fatal(FARGS,
                  "Cannot have %d separate PME ranks because: %s",
                  numPmeRanksRequested,
                  separatePmeRanksPermitted.reasonsWhyDisabled().c_str());
    }

    // Once the rank count is large enough, it becomes worth
    // suggesting improvements to the user.
    const int minPPRankCountToCheckForLargePrimeFactors = 13;
    if (checkForLargePrimeFactors && numPPRanksRequested >= minPPRankCountToCheckForLargePrimeFactors)
    {
        const int largestDivisor = largest_divisor(numPPRanksRequested);
        // Don't abort if the largest divisor is <= this value
        const int c_maxAllowedPrimeFactor = 7;
        /* Check if the largest divisor is more than numPPRanks ^ (2/3) */
        if (largestDivisor * largestDivisor * largestDivisor > numPPRanksRequested * numPPRanksRequested
            && largestDivisor > c_maxAllowedPrimeFactor)
        {
            gmx_fatal(FARGS,
                      "The number of ranks selected for particle-particle work (%d) "
                      "contains a large prime factor %d. In most cases this will lead to "
                      "bad performance. Choose a number with smaller prime factors or "
                      "set the decomposition (option -dd) manually.",
                      numPPRanksRequested,
                      largestDivisor);
        }
    }
}

/*! \brief Return the number of PME-only ranks used by the simulation
 *
 * If the user did not choose a number, then decide for them. */
static int getNumPmeOnlyRanksToUse(const gmx::MDLogger&                  mdlog,
                                   const gmx::DomdecOptions&             options,
                                   const gmx_mtop_t&                     mtop,
                                   const t_inputrec&                     ir,
                                   const gmx::SeparatePmeRanksPermitted& separatePmeRanksPermitted,
                                   const matrix                          box,
                                   const int                             numRanksRequested)
{
    int numPmeOnlyRanks = 0;

    if (!(usingPme(ir.coulombtype) || usingLJPme(ir.vdwtype)))
    {
        // System does not use PME for Electrostatics or LJ
        numPmeOnlyRanks = 0;
        GMX_LOG(mdlog.info)
                .appendTextFormatted("The system does not use PME for electrostatics or LJ");
    }
    else
    {
        std::string extraMessage;

        if (options.numPmeRanks >= 0)
        {
            numPmeOnlyRanks = options.numPmeRanks;
            extraMessage += ", as requested with -npme option";
        }
        else
        {
            /* Disable PME-only ranks if some parts of the code requested so and it's up to GROMACS to decide */
            if (!separatePmeRanksPermitted.permitSeparatePmeRanks())
            {
                numPmeOnlyRanks = 0;
                extraMessage += " because: " + separatePmeRanksPermitted.reasonsWhyDisabled();
            }
            else
            {
                numPmeOnlyRanks = guess_npme(mdlog, mtop, ir, box, numRanksRequested);
                extraMessage += ", as guessed by mdrun";
            }
        }

        GMX_RELEASE_ASSERT(numPmeOnlyRanks <= numRanksRequested,
                           "Cannot have more PME ranks than total ranks");

        GMX_LOG(mdlog.info)
                .appendTextFormatted(
                        "Using %d separate PME ranks%s", numPmeOnlyRanks, extraMessage.c_str());
    }

    return numPmeOnlyRanks;
}

/*! \brief Sets the order of the DD dimensions, returns the number of DD dimensions */
static int set_dd_dim(const gmx::IVec& numDDCells, const DDSettings& ddSettings, ivec* dims)
{
    int ndim = 0;
    if (ddSettings.useDDOrderZYX)
    {
        /* Decomposition order z,y,x */
        for (int dim = DIM - 1; dim >= 0; dim--)
        {
            if (numDDCells[dim] > 1)
            {
                (*dims)[ndim++] = dim;
            }
        }
    }
    else
    {
        /* Decomposition order x,y,z */
        for (int dim = 0; dim < DIM; dim++)
        {
            if (numDDCells[dim] > 1)
            {
                (*dims)[ndim++] = dim;
            }
        }
    }

    if (ndim == 0)
    {
        /* Set dim[0] to avoid extra checks on ndim in several places */
        (*dims)[0] = XX;
    }

    return ndim;
}

DDGridSetup getDDGridSetup(const gmx::MDLogger&                  mdlog,
                           DDRole                                ddRole,
                           MPI_Comm                              communicator,
                           const int                             numRanksRequested,
                           const DomdecOptions&                  options,
                           const DDSettings&                     ddSettings,
                           const DDSystemInfo&                   systemInfo,
                           const real                            cellSizeLimit,
                           const gmx_mtop_t&                     mtop,
                           const t_inputrec&                     ir,
                           const gmx::SeparatePmeRanksPermitted& separatePmeRanksPermitted,
                           const matrix                          box,
                           gmx::ArrayRef<const gmx::RVec>        xGlobal,
                           gmx_ddbox_t*                          ddbox)
{
    int numPmeOnlyRanks = getNumPmeOnlyRanksToUse(
            mdlog, options, mtop, ir, separatePmeRanksPermitted, box, numRanksRequested);

    gmx::IVec numDomains;
    if (options.numCells[XX] > 0)
    {
        numDomains                      = gmx::IVec(options.numCells);
        const ivec numDomainsLegacyIvec = { numDomains[XX], numDomains[YY], numDomains[ZZ] };
        set_ddbox_cr(ddRole, communicator, &numDomainsLegacyIvec, ir, box, xGlobal, ddbox);
    }
    else
    {
        set_ddbox_cr(ddRole, communicator, nullptr, ir, box, xGlobal, ddbox);

        if (ddRole == DDRole::Main)
        {
            numDomains = optimizeDDCells(
                    mdlog, numRanksRequested, numPmeOnlyRanks, cellSizeLimit, mtop, box, *ddbox, ir, systemInfo);
        }
    }

    /* Communicate the information set by the coordinator to all ranks */
    gmx_bcast(sizeof(numDomains), numDomains, communicator);
    if (usingPme(ir.coulombtype))
    {
        gmx_bcast(sizeof(numPmeOnlyRanks), &numPmeOnlyRanks, communicator);
    }

    DDGridSetup ddGridSetup;
    ddGridSetup.numPmeOnlyRanks = numPmeOnlyRanks;
    ddGridSetup.numDomains[XX]  = numDomains[XX];
    ddGridSetup.numDomains[YY]  = numDomains[YY];
    ddGridSetup.numDomains[ZZ]  = numDomains[ZZ];
    ddGridSetup.numDDDimensions = set_dd_dim(numDomains, ddSettings, &ddGridSetup.ddDimensions);

    return ddGridSetup;
}
