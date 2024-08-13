/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#include "gmxpre.h"

#include "kernel_gpu_ref.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/atomdata.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

static constexpr int c_clSize = c_nbnxnGpuClusterSize;

void nbnxn_kernel_gpu_ref(const NbnxnPairlistGpu*    nbl,
                          const nbnxn_atomdata_t*    nbat,
                          const interaction_const_t* iconst,
                          ArrayRef<const RVec>       shiftvec,
                          const StepWorkload&        stepWork,
                          int                        clearF,
                          ArrayRef<real>             f,
                          real*                      fshift,
                          real*                      Vc,
                          real*                      Vvdw)
{
    real                fscal = NAN;
    real                vcoul = 0;
    const nbnxn_excl_t* excl[c_nbnxnGpuClusterpairSplit];

    if (nbl->na_ci != c_clSize)
    {
        gmx_fatal(FARGS,
                  "The neighborlist cluster size in the GPU reference kernel is %d, expected it to "
                  "be %d",
                  nbl->na_ci,
                  c_clSize);
    }

    if (clearF == enbvClearFYes)
    {
        for (real& elem : f)
        {
            elem = 0;
        }
    }

    const bool bEwald = usingFullElectrostatics(iconst->eeltype);

    const real rcut2 = iconst->rcoulomb * iconst->rcoulomb;
    const real rvdw2 = iconst->rvdw * iconst->rvdw;

    const real rlist2 = nbl->rlist * nbl->rlist;

    const int*  type     = nbat->params().type.data();
    const real  facel    = iconst->epsfac;
    const real* vdwparam = nbat->params().nbfp.data();
    const int   ntype    = nbat->params().numTypes;

    const real* x = nbat->x().data();

    int npair_tot   = 0;
    int nhwu        = 0;
    int nhwu_pruned = 0;

    for (const nbnxn_sci_t& nbln : nbl->sci)
    {
        const int  ish           = nbln.shift;
        const int  ish3          = DIM * ish;
        const real shX           = shiftvec[ish][XX];
        const real shY           = shiftvec[ish][YY];
        const real shZ           = shiftvec[ish][ZZ];
        const int  cjPackedBegin = nbln.cjPackedBegin;
        const int  cjPackedEnd   = nbln.cjPackedEnd;
        const int  sci           = nbln.sci;
        real       vctot         = 0;
        real       Vvdwtot       = 0;

        if (nbln.shift == c_centralShiftIndex
            && nbl->cjPacked.list_[cjPackedBegin].cj[0] == sci * c_nbnxnGpuNumClusterPerSupercluster)
        {
            /* we have the diagonal:
             * add the charge self interaction energy term
             */
            for (int im = 0; im < c_nbnxnGpuNumClusterPerSupercluster; im++)
            {
                const int ci = sci * c_nbnxnGpuNumClusterPerSupercluster + im;
                for (int ic = 0; ic < c_clSize; ic++)
                {
                    const int ia = ci * c_clSize + ic;
                    real      iq = x[ia * nbat->xstride + 3];
                    vctot += iq * iq;
                }
            }
            if (!bEwald)
            {
                vctot *= -facel * 0.5 * iconst->reactionFieldShift;
            }
            else
            {
                /* last factor 1/sqrt(pi) */
                vctot *= -facel * iconst->ewaldcoeff_q * M_1_SQRTPI;
            }
        }

        for (int cjPacked = cjPackedBegin; (cjPacked < cjPackedEnd); cjPacked++)
        {
            for (int splitIdx = 0; splitIdx < c_nbnxnGpuClusterpairSplit; splitIdx++)
            {
                excl[splitIdx] = &nbl->excl[nbl->cjPacked.list_[cjPacked].imei[splitIdx].excl_ind];
            }

            for (int jm = 0; jm < c_nbnxnGpuJgroupSize; jm++)
            {
                const int cj = nbl->cjPacked.list_[cjPacked].cj[jm];

                for (int im = 0; im < c_nbnxnGpuNumClusterPerSupercluster; im++)
                {
                    /* We're only using the first imask,
                     * but here imei[1].imask is identical.
                     */
                    if ((nbl->cjPacked.list_[cjPacked].imei[0].imask
                         >> (jm * c_nbnxnGpuNumClusterPerSupercluster + im))
                        & 1)
                    {
                        const int ci = sci * c_nbnxnGpuNumClusterPerSupercluster + im;

                        bool within_rlist = false;
                        int  npair        = 0;
                        for (int ic = 0; ic < c_clSize; ic++)
                        {
                            const int ia = ci * c_clSize + ic;

                            const int  is  = ia * nbat->xstride;
                            const int  ifs = ia * nbat->fstride;
                            const real ix  = shX + x[is + 0];
                            const real iy  = shY + x[is + 1];
                            const real iz  = shZ + x[is + 2];
                            const real iq  = facel * x[is + 3];
                            const int  nti = ntype * 2 * type[ia];

                            real fix = 0;
                            real fiy = 0;
                            real fiz = 0;

                            for (int jc = 0; jc < c_clSize; jc++)
                            {
                                const int ja = cj * c_clSize + jc;

                                if (nbln.shift == c_centralShiftIndex && ci == cj && ja <= ia)
                                {
                                    continue;
                                }

                                constexpr int clusterPerSplit =
                                        c_nbnxnGpuClusterSize / c_nbnxnGpuClusterpairSplit;
                                const real int_bit = static_cast<real>(
                                        (excl[jc / clusterPerSplit]->pair[(jc & (clusterPerSplit - 1)) * c_clSize + ic]
                                         >> (jm * c_nbnxnGpuNumClusterPerSupercluster + im))
                                        & 1);

                                int        js  = ja * nbat->xstride;
                                int        jfs = ja * nbat->fstride;
                                const real jx  = x[js + 0];
                                const real jy  = x[js + 1];
                                const real jz  = x[js + 2];
                                const real dx  = ix - jx;
                                const real dy  = iy - jy;
                                const real dz  = iz - jz;
                                real       rsq = dx * dx + dy * dy + dz * dz;
                                if (rsq < rlist2)
                                {
                                    within_rlist = true;
                                }
                                if (rsq >= rcut2)
                                {
                                    continue;
                                }

                                if (type[ia] != ntype - 1 && type[ja] != ntype - 1)
                                {
                                    npair++;
                                }

                                // Ensure distance do not become so small that r^-12 overflows
                                rsq = std::max(rsq, c_nbnxnMinDistanceSquared);

                                const real rinv   = invsqrt(rsq);
                                const real rinvsq = rinv * rinv;

                                const real qq = iq * x[js + 3];
                                if (!bEwald)
                                {
                                    /* Reaction-field */
                                    const real krsq = iconst->reactionFieldCoefficient * rsq;
                                    fscal           = qq * (int_bit * rinv - 2 * krsq) * rinvsq;
                                    if (stepWork.computeEnergy)
                                    {
                                        vcoul = qq * (int_bit * rinv + krsq - iconst->reactionFieldShift);
                                    }
                                }
                                else
                                {
                                    const real  r    = rsq * rinv;
                                    const real  rt   = r * iconst->coulombEwaldTables->scale;
                                    const int   n0   = static_cast<int>(rt);
                                    const real  eps  = rt - static_cast<real>(n0);
                                    const real* Ftab = iconst->coulombEwaldTables->tableF.data();

                                    const real fexcl = (1 - eps) * Ftab[n0] + eps * Ftab[n0 + 1];

                                    fscal = qq * (int_bit * rinvsq - fexcl) * rinv;

                                    if (stepWork.computeEnergy)
                                    {
                                        vcoul = qq
                                                * ((int_bit - std::erf(iconst->ewaldcoeff_q * r)) * rinv
                                                   - int_bit * iconst->sh_ewald);
                                    }
                                }

                                if (rsq < rvdw2)
                                {
                                    const int tj = nti + 2 * type[ja];

                                    /* Vanilla Lennard-Jones cutoff */
                                    const real c6  = vdwparam[tj];
                                    const real c12 = vdwparam[tj + 1];

                                    const real rinvsix   = int_bit * rinvsq * rinvsq * rinvsq;
                                    const real Vvdw_disp = c6 * rinvsix;
                                    const real Vvdw_rep  = c12 * rinvsix * rinvsix;
                                    fscal += (Vvdw_rep - Vvdw_disp) * rinvsq;

                                    if (stepWork.computeEnergy)
                                    {
                                        vctot += vcoul;

                                        Vvdwtot +=
                                                (Vvdw_rep + int_bit * c12 * iconst->repulsion_shift.cpot) / 12
                                                - (Vvdw_disp
                                                   + int_bit * c6 * iconst->dispersion_shift.cpot)
                                                          / 6;
                                    }
                                }

                                real tx = fscal * dx;
                                real ty = fscal * dy;
                                real tz = fscal * dz;
                                fix     = fix + tx;
                                fiy     = fiy + ty;
                                fiz     = fiz + tz;
                                f[jfs + 0] -= tx;
                                f[jfs + 1] -= ty;
                                f[jfs + 2] -= tz;
                            }

                            f[ifs + 0] += fix;
                            f[ifs + 1] += fiy;
                            f[ifs + 2] += fiz;
                            fshift[ish3]     = fshift[ish3] + fix;
                            fshift[ish3 + 1] = fshift[ish3 + 1] + fiy;
                            fshift[ish3 + 2] = fshift[ish3 + 2] + fiz;

                            /* Count in half work-units.
                             * In CUDA one work-unit is 2 warps.
                             */
                            if ((ic + 1) % (c_clSize / c_nbnxnGpuClusterpairSplit) == 0)
                            {
                                npair_tot += npair;

                                nhwu++;
                                if (within_rlist)
                                {
                                    nhwu_pruned++;
                                }

                                within_rlist = false;
                                npair        = 0;
                            }
                        }
                    }
                }
            }
        }

        if (stepWork.computeEnergy)
        {
            const int ggid = 0;
            Vc[ggid]       = Vc[ggid] + vctot;
            Vvdw[ggid]     = Vvdw[ggid] + Vvdwtot;
        }
    }

    if (debug)
    {
        fprintf(debug,
                "number of half %dx%d atom pairs: %d after pruning: %d fraction %4.2f\n",
                nbl->na_ci,
                nbl->na_ci,
                nhwu,
                nhwu_pruned,
                nhwu_pruned / static_cast<double>(nhwu));
        fprintf(debug, "generic kernel pair interactions:            %d\n", nhwu * nbl->na_ci / 2 * nbl->na_ci);
        fprintf(debug,
                "generic kernel post-prune pair interactions: %d\n",
                nhwu_pruned * nbl->na_ci / 2 * nbl->na_ci);
        fprintf(debug, "generic kernel non-zero pair interactions:   %d\n", npair_tot);
        fprintf(debug,
                "ratio non-zero/post-prune pair interactions: %4.2f\n",
                npair_tot / static_cast<double>(nhwu_pruned * exactDiv(nbl->na_ci, 2) * nbl->na_ci));
    }
}

} // namespace gmx
