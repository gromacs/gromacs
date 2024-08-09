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
#include "gmxpre.h"

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/pulling/pull_internal.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

using gmx::ArrayRef;
using gmx::RVec;

#if GMX_MPI

// Helper function to deduce MPI datatype from the type of data
gmx_unused static MPI_Datatype mpiDatatype(const float gmx_unused* data)
{
    return MPI_FLOAT;
}

// Helper function to deduce MPI datatype from the type of data
gmx_unused static MPI_Datatype mpiDatatype(const double gmx_unused* data)
{
    return MPI_DOUBLE;
}

#endif // GMX_MPI

#if !GMX_DOUBLE
// Helper function; note that gmx_sum(d) should actually be templated
gmx_unused static void gmxAllReduce(int n, real* data, const t_commrec* cr)
{
    gmx_sum(n, data, cr);
}
#endif

// Helper function; note that gmx_sum(d) should actually be templated
gmx_unused static void gmxAllReduce(int n, double* data, const t_commrec* cr)
{
    gmx_sumd(n, data, cr);
}

// Reduce data of n elements over all ranks currently participating in pull
template<typename T>
static void pullAllReduce(const t_commrec* cr, pull_comm_t* comm, int n, T* data)
{
    if (cr != nullptr && PAR(cr))
    {
        if (comm->bParticipateAll)
        {
            /* Sum the contributions over all DD ranks */
            gmxAllReduce(n, data, cr);
        }
        else
        {
            /* Separate branch because gmx_sum uses cr->mpi_comm_mygroup */
#if GMX_MPI
            MPI_Allreduce(MPI_IN_PLACE, data, n, mpiDatatype(data), MPI_SUM, comm->mpi_comm_com);
#else
            gmx_incons("comm->bParticipateAll=FALSE without GMX_MPI");
#endif
        }
    }
}

/* Copies the coordinates of the PBC atom of pgrp to x_pbc.
 * When those coordinates are not available on this rank, clears x_pbc.
 */
static void setPbcAtomCoords(const pull_group_work_t& pgrp, ArrayRef<const RVec> x, rvec x_pbc)
{
    if (pgrp.pbcAtomSet != nullptr)
    {
        if (pgrp.pbcAtomSet->numAtomsLocal() > 0)
        {
            /* We have the atom locally, copy its coordinates */
            copy_rvec(x[pgrp.pbcAtomSet->localIndex()[0]], x_pbc);
        }
        else
        {
            /* Another rank has it, clear the coordinates for MPI_Allreduce */
            clear_rvec(x_pbc);
        }
    }
    else
    {
        copy_rvec(x[pgrp.params_.pbcatom], x_pbc);
    }
}

static void pull_set_pbcatoms(const t_commrec* cr, struct pull_t* pull, ArrayRef<const RVec> x, ArrayRef<RVec> x_pbc)
{
    int numPbcAtoms = 0;
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        const pull_group_work_t& group = pull->group[g];
        if (group.needToCalcCom && (group.epgrppbc == epgrppbcREFAT || group.epgrppbc == epgrppbcPREVSTEPCOM))
        {
            setPbcAtomCoords(pull->group[g], x, x_pbc[g]);
            numPbcAtoms++;
        }
        else
        {
            clear_rvec(x_pbc[g]);
        }
    }

    if (cr && PAR(cr) && numPbcAtoms > 0)
    {
        /* Sum over participating ranks to get x_pbc from the home ranks.
         * This can be very expensive at high parallelization, so we only
         * do this after each DD repartitioning.
         */
        pullAllReduce(cr, &pull->comm, pull->group.size() * DIM, static_cast<real*>(x_pbc[0]));
    }
}

static void make_cyl_refgrps(const t_commrec*     cr,
                             pull_t*              pull,
                             ArrayRef<const real> masses,
                             const t_pbc&         pbc,
                             double               t,
                             ArrayRef<const RVec> x)
{
    pull_comm_t* comm = &pull->comm;

    GMX_ASSERT(comm->cylinderBuffer.size() == pull->coord.size() * c_cylinderBufferStride,
               "cylinderBuffer should have the correct size");

    double inv_cyl_r2 = 1.0 / gmx::square(pull->params.cylinder_r);

    /* loop over all groups to make a reference group for each*/
    int bufferOffset = 0;
    for (pull_coord_work_t& pcrd : pull->coord)
    {
        double sum_a, wmass, wwmass;
        dvec   radf_fac0, radf_fac1;

        sum_a  = 0;
        wmass  = 0;
        wwmass = 0;
        clear_dvec(radf_fac0);
        clear_dvec(radf_fac1);

        if (pcrd.params_.eGeom == PullGroupGeometry::Cylinder)
        {
            /* pref will be the same group for all pull coordinates */
            const pull_group_work_t& pref  = pull->group[pcrd.params_.group[0]];
            const pull_group_work_t& pgrp  = pull->group[pcrd.params_.group[1]];
            pull_group_work_t&       pdyna = *pcrd.dynamicGroup0;
            rvec                     direction;
            copy_dvec_to_rvec(pcrd.spatialData.vec, direction);

            /* Since we have not calculated the COM of the cylinder group yet,
             * we calculate distances with respect to location of the pull
             * group minus the reference position along the vector.
             * here we already have the COM of the pull group. This resolves
             * any PBC issues and we don't need to use a PBC-atom here.
             */
            if (pcrd.params_.rate != 0)
            {
                /* With rate=0, value_ref is set initially */
                pcrd.value_ref = pcrd.params_.init + pcrd.params_.rate * t;
            }
            rvec reference;
            for (int m = 0; m < DIM; m++)
            {
                reference[m] = pgrp.x[m] - pcrd.spatialData.vec[m] * pcrd.value_ref;
            }

            auto localAtomIndices = pref.atomSet_.localIndex();

            /* This actually only needs to be done at init or DD time,
             * but resizing with the same size does not cause much overhead.
             */
            pdyna.localWeights.resize(localAtomIndices.size());
            pdyna.mdw.resize(localAtomIndices.size());
            pdyna.dv.resize(localAtomIndices.size());

            /* loop over all atoms in the main ref group */
            for (gmx::Index indexInSet = 0; indexInSet < localAtomIndices.ssize(); indexInSet++)
            {
                int  atomIndex = localAtomIndices[indexInSet];
                rvec dx;
                pbc_dx_aiuc(&pbc, x[atomIndex], reference, dx);
                double axialLocation = iprod(direction, dx);
                dvec   radialLocation;
                double dr2 = 0;
                for (int m = 0; m < DIM; m++)
                {
                    /* Determine the radial components */
                    radialLocation[m] = dx[m] - axialLocation * direction[m];
                    dr2 += gmx::square(radialLocation[m]);
                }
                double dr2_rel = dr2 * inv_cyl_r2;

                if (dr2_rel < 1)
                {
                    /* add atom to sum of COM and to weight array */

                    double mass = masses[atomIndex];
                    /* The radial weight function is 1-2x^2+x^4,
                     * where x=r/cylinder_r. Since this function depends
                     * on the radial component, we also get radial forces
                     * on both groups.
                     */
                    double weight                  = 1 + (-2 + dr2_rel) * dr2_rel;
                    double dweight_r               = (-4 + 4 * dr2_rel) * inv_cyl_r2;
                    pdyna.localWeights[indexInSet] = weight;
                    sum_a += mass * weight * axialLocation;
                    wmass += mass * weight;
                    wwmass += mass * weight * weight;
                    dvec mdw;
                    dsvmul(mass * dweight_r, radialLocation, mdw);
                    copy_dvec(mdw, pdyna.mdw[indexInSet]);
                    /* Currently we only have the axial component of the
                     * offset from the cylinder COM up to an unkown offset.
                     * We add this offset after the reduction needed
                     * for determining the COM of the cylinder group.
                     */
                    pdyna.dv[indexInSet] = axialLocation;
                    for (int m = 0; m < DIM; m++)
                    {
                        radf_fac0[m] += mdw[m];
                        radf_fac1[m] += mdw[m] * axialLocation;
                    }
                }
                else
                {
                    pdyna.localWeights[indexInSet] = 0;
                }
            }
        }

        auto buffer = gmx::arrayRefFromArray(comm->cylinderBuffer.data() + bufferOffset,
                                             c_cylinderBufferStride);
        bufferOffset += c_cylinderBufferStride;

        buffer[0] = wmass;
        buffer[1] = wwmass;
        buffer[2] = sum_a;

        buffer[3] = radf_fac0[XX];
        buffer[4] = radf_fac0[YY];
        buffer[5] = radf_fac0[ZZ];

        buffer[6] = radf_fac1[XX];
        buffer[7] = radf_fac1[YY];
        buffer[8] = radf_fac1[ZZ];
    }

    if (cr != nullptr && PAR(cr))
    {
        /* Sum the contributions over the ranks */
        pullAllReduce(cr, comm, pull->coord.size() * c_cylinderBufferStride, comm->cylinderBuffer.data());
    }

    int pullCoordIndex = 0;
    bufferOffset       = 0;
    for (pull_coord_work_t& pcrd : pull->coord)
    {
        if (pcrd.params_.eGeom == PullGroupGeometry::Cylinder)
        {
            pull_group_work_t&    dynamicGroup0 = *pcrd.dynamicGroup0;
            pull_group_work_t&    group1        = pull->group[pcrd.params_.group[1]];
            PullCoordSpatialData& spatialData   = pcrd.spatialData;

            auto buffer = gmx::constArrayRefFromArray(comm->cylinderBuffer.data() + bufferOffset,
                                                      c_cylinderBufferStride);
            bufferOffset += c_cylinderBufferStride;
            const double wmass  = buffer[0];
            const double wwmass = buffer[1];
            if (wmass == 0)
            {
                gmx_fatal(FARGS,
                          "The mass of the cylinder selection of pull group %d is zero. "
                          "This means that group %d does not cover the whole area. "
                          "Maybe you did not put the correct group as group 1 in the pull "
                          "coordinate?",
                          pullCoordIndex + 1,
                          pullCoordIndex + 1);
            }
            dynamicGroup0.mwscale = 1.0 / wmass;
            /* Cylinder pulling can't be used with constraints, but we set
             * wscale and invtm anyhow, in case someone would like to use them.
             */
            dynamicGroup0.wscale = wmass / wwmass;
            dynamicGroup0.invtm  = wwmass / (wmass * wmass);

            /* We store the deviation of the COM from the reference location
             * used above, since we need it when we apply the radial forces
             * to the atoms in the cylinder group.
             */
            spatialData.cyl_dev = 0;
            for (int m = 0; m < DIM; m++)
            {
                double reference   = group1.x[m] - spatialData.vec[m] * pcrd.value_ref;
                double dist        = -spatialData.vec[m] * buffer[2] * dynamicGroup0.mwscale;
                dynamicGroup0.x[m] = reference - dist;
                spatialData.cyl_dev += dist;
            }
            /* Now we know the exact COM of the cylinder reference group,
             * we can determine the radial force factor (ffrad) that when
             * multiplied with the axial pull force will give the radial
             * force on the pulled (non-cylinder) group.
             */
            for (int m = 0; m < DIM; m++)
            {
                spatialData.ffrad[m] = (buffer[6 + m] + buffer[3 + m] * spatialData.cyl_dev) / wmass;
            }
        }

        pullCoordIndex++;
    }
}

static double atan2_0_2pi(double y, double x)
{
    double a;

    a = std::atan2(y, x);
    if (a < 0)
    {
        a += 2.0 * M_PI;
    }
    return a;
}

static void sum_com_part(const pull_group_work_t* pgrp,
                         int                      ind_start,
                         int                      ind_end,
                         ArrayRef<const RVec>     x,
                         ArrayRef<const RVec>     xp,
                         ArrayRef<const real>     mass,
                         const t_pbc&             pbc,
                         const rvec               x_pbc,
                         ComSums*                 sum_com)
{
    double sum_wm   = 0;
    double sum_wwm  = 0;
    dvec   sum_wmx  = { 0, 0, 0 };
    dvec   sum_wmxp = { 0, 0, 0 };

    auto localAtomIndices = pgrp->atomSet_.localIndex();
    for (int i = ind_start; i < ind_end; i++)
    {
        int  ii = localAtomIndices[i];
        real wm;
        if (pgrp->localWeights.empty())
        {
            wm = mass[ii];
            sum_wm += wm;
        }
        else
        {
            real w;

            w  = pgrp->localWeights[i];
            wm = w * mass[ii];
            sum_wm += wm;
            sum_wwm += wm * w;
        }
        if (pgrp->epgrppbc == epgrppbcNONE)
        {
            /* Plain COM: sum the coordinates */
            for (int d = 0; d < DIM; d++)
            {
                sum_wmx[d] += wm * x[ii][d];
            }
            if (!xp.empty())
            {
                for (int d = 0; d < DIM; d++)
                {
                    sum_wmxp[d] += wm * xp[ii][d];
                }
            }
        }
        else
        {
            rvec dx;

            /* Sum the difference with the reference atom */
            pbc_dx(&pbc, x[ii], x_pbc, dx);
            for (int d = 0; d < DIM; d++)
            {
                sum_wmx[d] += wm * dx[d];
            }
            if (!xp.empty())
            {
                /* For xp add the difference between xp and x to dx,
                 * such that we use the same periodic image,
                 * also when xp has a large displacement.
                 */
                for (int d = 0; d < DIM; d++)
                {
                    sum_wmxp[d] += wm * (dx[d] + xp[ii][d] - x[ii][d]);
                }
            }
        }
    }

    sum_com->sum_wm  = sum_wm;
    sum_com->sum_wwm = sum_wwm;
    copy_dvec(sum_wmx, sum_com->sum_wmx);
    if (!xp.empty())
    {
        copy_dvec(sum_wmxp, sum_com->sum_wmxp);
    }
}

static void sum_com_part_cosweight(const pull_group_work_t* pgrp,
                                   int                      ind_start,
                                   int                      ind_end,
                                   int                      cosdim,
                                   real                     twopi_box,
                                   ArrayRef<const RVec>     x,
                                   ArrayRef<const RVec>     xp,
                                   ArrayRef<const real>     mass,
                                   ComSums*                 sum_com)
{
    /* Cosine weighting geometry */
    double sum_cm  = 0;
    double sum_sm  = 0;
    double sum_ccm = 0;
    double sum_csm = 0;
    double sum_ssm = 0;
    double sum_cmp = 0;
    double sum_smp = 0;

    auto localAtomIndices = pgrp->atomSet_.localIndex();

    for (int i = ind_start; i < ind_end; i++)
    {
        int  ii = localAtomIndices[i];
        real m  = mass[ii];
        /* Determine cos and sin sums */
        real cw = std::cos(x[ii][cosdim] * twopi_box);
        real sw = std::sin(x[ii][cosdim] * twopi_box);
        sum_cm += static_cast<double>(cw * m);
        sum_sm += static_cast<double>(sw * m);
        sum_ccm += static_cast<double>(cw * cw * m);
        sum_csm += static_cast<double>(cw * sw * m);
        sum_ssm += static_cast<double>(sw * sw * m);

        if (!xp.empty())
        {
            real cw = std::cos(xp[ii][cosdim] * twopi_box);
            real sw = std::sin(xp[ii][cosdim] * twopi_box);
            sum_cmp += static_cast<double>(cw * m);
            sum_smp += static_cast<double>(sw * m);
        }
    }

    sum_com->sum_cm  = sum_cm;
    sum_com->sum_sm  = sum_sm;
    sum_com->sum_ccm = sum_ccm;
    sum_com->sum_csm = sum_csm;
    sum_com->sum_ssm = sum_ssm;
    sum_com->sum_cmp = sum_cmp;
    sum_com->sum_smp = sum_smp;
}

/* calculates center of mass of selection index from all coordinates x */
void pull_calc_coms(const t_commrec*     cr,
                    pull_t*              pull,
                    ArrayRef<const real> masses,
                    const t_pbc&         pbc,
                    double               t,
                    ArrayRef<const RVec> x,
                    ArrayRef<RVec>       xp)
{
    real         twopi_box = 0;
    pull_comm_t* comm;

    comm = &pull->comm;

    GMX_ASSERT(comm->pbcAtomBuffer.size() == pull->group.size(),
               "pbcAtomBuffer should have size number of groups");
    GMX_ASSERT(comm->comBuffer.size() == pull->group.size() * c_comBufferStride,
               "comBuffer should have size #group*c_comBufferStride");

    if (pull->bRefAt && pull->bSetPBCatoms)
    {
        pull_set_pbcatoms(cr, pull, x, comm->pbcAtomBuffer);

        if (cr != nullptr && haveDDAtomOrdering(*cr))
        {
            /* We can keep these PBC reference coordinates fixed for nstlist
             * steps, since atoms won't jump over PBC.
             * This avoids a global reduction at the next nstlist-1 steps.
             * Note that the exact values of the pbc reference coordinates
             * are irrelevant, as long all atoms in the group are within
             * half a box distance of the reference coordinate.
             */
            pull->bSetPBCatoms = FALSE;
        }
    }

    if (pull->cosdim >= 0)
    {
        int m;

        assert(pull->npbcdim <= DIM);

        for (m = pull->cosdim + 1; m < pull->npbcdim; m++)
        {
            if (pbc.box[m][pull->cosdim] != 0)
            {
                gmx_fatal(FARGS, "Can not do cosine weighting for trilinic dimensions");
            }
        }
        twopi_box = 2.0 * M_PI / pbc.box[pull->cosdim][pull->cosdim];
    }

    for (size_t g = 0; g < pull->group.size(); g++)
    {
        pull_group_work_t* pgrp = &pull->group[g];

        /* Cosine-weighted COMs behave different from all other weighted COMs
         * in the sense that the weights depend on instantaneous coordinates,
         * not on pre-set weights. Thus we resize the local weight buffer here.
         */
        if (pgrp->epgrppbc == epgrppbcCOS)
        {
            pgrp->localWeights.resize(pgrp->atomSet_.localIndex().size());
        }

        auto comBuffer = gmx::arrayRefFromArray(comm->comBuffer.data() + g * c_comBufferStride,
                                                c_comBufferStride);

        if (pgrp->needToCalcCom)
        {
            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                rvec x_pbc = { 0, 0, 0 };

                switch (pgrp->epgrppbc)
                {
                    case epgrppbcREFAT:
                        /* Set the pbc atom */
                        copy_rvec(comm->pbcAtomBuffer[g], x_pbc);
                        break;
                    case epgrppbcPREVSTEPCOM:
                        /* Set the pbc reference to the COM of the group of the last step */
                        copy_dvec_to_rvec(pgrp->x_prev_step, comm->pbcAtomBuffer[g]);
                        copy_dvec_to_rvec(pgrp->x_prev_step, x_pbc);
                }

                /* The final sums should end up in comSums[0] */
                ComSums& comSumsTotal = pull->comSums[0];

                /* If we have a single-atom group the mass is irrelevant, so
                 * we can remove the mass factor to avoid division by zero.
                 * Note that with constraint pulling the mass does matter, but
                 * in that case a check group mass != 0 has been done before.
                 */
                if (pgrp->params_.ind.size() == 1 && pgrp->atomSet_.numAtomsLocal() == 1
                    && masses[pgrp->atomSet_.localIndex()[0]] == 0)
                {
                    GMX_ASSERT(xp.empty(),
                               "We should not have groups with zero mass with constraints, i.e. "
                               "xp not empty");

                    /* Copy the single atom coordinate */
                    for (int d = 0; d < DIM; d++)
                    {
                        comSumsTotal.sum_wmx[d] = x[pgrp->atomSet_.localIndex()[0]][d];
                    }
                    /* Set all mass factors to 1 to get the correct COM */
                    comSumsTotal.sum_wm  = 1;
                    comSumsTotal.sum_wwm = 1;
                }
                else if (pgrp->atomSet_.numAtomsLocal() <= c_pullMaxNumLocalAtomsSingleThreaded)
                {
                    sum_com_part(pgrp, 0, pgrp->atomSet_.numAtomsLocal(), x, xp, masses, pbc, x_pbc, &comSumsTotal);
                }
                else
                {
                    const int numThreads = pgrp->numThreads();
#pragma omp parallel for num_threads(numThreads) schedule(static)
                    for (int thread = 0; thread < numThreads; thread++)
                    {
                        int ind_start = (pgrp->atomSet_.numAtomsLocal() * (thread + 0)) / numThreads;
                        int ind_end = (pgrp->atomSet_.numAtomsLocal() * (thread + 1)) / numThreads;
                        sum_com_part(
                                pgrp, ind_start, ind_end, x, xp, masses, pbc, x_pbc, &pull->comSums[thread]);
                    }

                    /* Reduce the thread contributions to sum_com[0] */
                    for (int thread = 1; thread < numThreads; thread++)
                    {
                        comSumsTotal.sum_wm += pull->comSums[thread].sum_wm;
                        comSumsTotal.sum_wwm += pull->comSums[thread].sum_wwm;
                        dvec_inc(comSumsTotal.sum_wmx, pull->comSums[thread].sum_wmx);
                        dvec_inc(comSumsTotal.sum_wmxp, pull->comSums[thread].sum_wmxp);
                    }
                }

                if (pgrp->localWeights.empty())
                {
                    comSumsTotal.sum_wwm = comSumsTotal.sum_wm;
                }

                /* Copy local sums to a buffer for global summing */
                copy_dvec(comSumsTotal.sum_wmx, comBuffer[0]);

                copy_dvec(comSumsTotal.sum_wmxp, comBuffer[1]);

                comBuffer[2][0] = comSumsTotal.sum_wm;
                comBuffer[2][1] = comSumsTotal.sum_wwm;
                comBuffer[2][2] = 0;
            }
            else
            {
                /* Cosine weighting geometry.
                 * This uses a slab of the system, thus we always have many
                 * atoms in the pull groups. Therefore, always use threads.
                 */
                const int numThreads = pgrp->numThreads();
#pragma omp parallel for num_threads(numThreads) schedule(static)
                for (int thread = 0; thread < numThreads; thread++)
                {
                    int ind_start = (pgrp->atomSet_.numAtomsLocal() * (thread + 0)) / numThreads;
                    int ind_end   = (pgrp->atomSet_.numAtomsLocal() * (thread + 1)) / numThreads;
                    sum_com_part_cosweight(
                            pgrp, ind_start, ind_end, pull->cosdim, twopi_box, x, xp, masses, &pull->comSums[thread]);
                }

                /* Reduce the thread contributions to comSums[0] */
                ComSums& comSumsTotal = pull->comSums[0];
                for (int thread = 1; thread < numThreads; thread++)
                {
                    comSumsTotal.sum_cm += pull->comSums[thread].sum_cm;
                    comSumsTotal.sum_sm += pull->comSums[thread].sum_sm;
                    comSumsTotal.sum_ccm += pull->comSums[thread].sum_ccm;
                    comSumsTotal.sum_csm += pull->comSums[thread].sum_csm;
                    comSumsTotal.sum_ssm += pull->comSums[thread].sum_ssm;
                    comSumsTotal.sum_cmp += pull->comSums[thread].sum_cmp;
                    comSumsTotal.sum_smp += pull->comSums[thread].sum_smp;
                }

                /* Copy local sums to a buffer for global summing */
                comBuffer[0][0] = comSumsTotal.sum_cm;
                comBuffer[0][1] = comSumsTotal.sum_sm;
                comBuffer[0][2] = 0;
                comBuffer[1][0] = comSumsTotal.sum_ccm;
                comBuffer[1][1] = comSumsTotal.sum_csm;
                comBuffer[1][2] = comSumsTotal.sum_ssm;
                comBuffer[2][0] = comSumsTotal.sum_cmp;
                comBuffer[2][1] = comSumsTotal.sum_smp;
                comBuffer[2][2] = 0;
            }
        }
        else
        {
            clear_dvec(comBuffer[0]);
            clear_dvec(comBuffer[1]);
            clear_dvec(comBuffer[2]);
        }
    }

    pullAllReduce(cr,
                  comm,
                  pull->group.size() * c_comBufferStride * DIM,
                  static_cast<double*>(comm->comBuffer[0]));

    for (size_t g = 0; g < pull->group.size(); g++)
    {
        pull_group_work_t* pgrp;

        pgrp = &pull->group[g];
        if (pgrp->needToCalcCom)
        {
            GMX_ASSERT(!pgrp->params_.ind.empty(),
                       "Normal pull groups should have atoms, only group 0, which should have "
                       "bCalcCom=FALSE has nat=0");

            const auto comBuffer = gmx::constArrayRefFromArray(
                    comm->comBuffer.data() + g * c_comBufferStride, c_comBufferStride);

            if (pgrp->epgrppbc != epgrppbcCOS)
            {
                double wmass, wwmass;
                int    m;

                /* Determine the inverse mass */
                wmass         = comBuffer[2][0];
                wwmass        = comBuffer[2][1];
                pgrp->mwscale = 1.0 / wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm != 0)
                {
                    pgrp->wscale = wmass / wwmass;
                    pgrp->invtm  = wwmass / (wmass * wmass);
                }
                /* Divide by the total mass */
                for (m = 0; m < DIM; m++)
                {
                    pgrp->x[m] = comBuffer[0][m] * pgrp->mwscale;
                    if (!xp.empty())
                    {
                        pgrp->xp[m] = comBuffer[1][m] * pgrp->mwscale;
                    }
                    if (pgrp->epgrppbc == epgrppbcREFAT || pgrp->epgrppbc == epgrppbcPREVSTEPCOM)
                    {
                        pgrp->x[m] += comm->pbcAtomBuffer[g][m];
                        if (!xp.empty())
                        {
                            pgrp->xp[m] += comm->pbcAtomBuffer[g][m];
                        }
                    }
                }
            }
            else
            {
                /* Cosine weighting geometry */
                double csw, snw, wmass, wwmass;

                /* Determine the optimal location of the cosine weight */
                csw                   = comBuffer[0][0];
                snw                   = comBuffer[0][1];
                pgrp->x[pull->cosdim] = atan2_0_2pi(snw, csw) / twopi_box;
                /* Set the weights for the local atoms */
                wmass  = std::sqrt(csw * csw + snw * snw);
                wwmass = (comBuffer[1][0] * csw * csw + comBuffer[1][1] * csw * snw
                          + comBuffer[1][2] * snw * snw)
                         / (wmass * wmass);

                pgrp->mwscale = 1.0 / wmass;
                pgrp->wscale  = wmass / wwmass;
                pgrp->invtm   = wwmass / (wmass * wmass);
                /* Set the weights for the local atoms */
                csw *= pgrp->invtm;
                snw *= pgrp->invtm;
                for (size_t i = 0; i < pgrp->atomSet_.numAtomsLocal(); i++)
                {
                    int ii                = pgrp->atomSet_.localIndex()[i];
                    pgrp->localWeights[i] = csw * std::cos(twopi_box * x[ii][pull->cosdim])
                                            + snw * std::sin(twopi_box * x[ii][pull->cosdim]);
                }
                if (!xp.empty())
                {
                    csw                    = comBuffer[2][0];
                    snw                    = comBuffer[2][1];
                    pgrp->xp[pull->cosdim] = atan2_0_2pi(snw, csw) / twopi_box;
                }
            }
        }
    }

    if (pull->bCylinder)
    {
        /* Calculate the COMs for the cyclinder reference groups */
        make_cyl_refgrps(cr, pull, masses, pbc, t, x);
    }
}

using BoolVec = gmx::BasicVector<bool>;

/* Returns whether the pull group obeys the PBC restrictions */
static bool pullGroupObeysPbcRestrictions(const pull_group_work_t& group,
                                          const BoolVec&           dimUsed,
                                          ArrayRef<const RVec>     x,
                                          const t_pbc&             pbc,
                                          const RVec&              x_pbc,
                                          const real               pbcMargin)
{
    /* Determine which dimensions are relevant for PBC */
    BoolVec dimUsesPbc       = { false, false, false };
    bool    pbcIsRectangular = true;
    for (int d = 0; d < pbc.ndim_ePBC; d++)
    {
        if (dimUsed[d])
        {
            dimUsesPbc[d] = true;
            /* All non-zero dimensions of vector v are involved in PBC */
            for (int d2 = d + 1; d2 < pbc.ndim_ePBC; d2++)
            {
                assert(d2 < DIM);
                if (pbc.box[d2][d] != 0)
                {
                    dimUsesPbc[d2]   = true;
                    pbcIsRectangular = false;
                }
            }
        }
    }

    rvec marginPerDim    = {};
    real marginDistance2 = 0;
    if (pbcIsRectangular)
    {
        /* Use margins for dimensions independently */
        for (int d = 0; d < pbc.ndim_ePBC; d++)
        {
            marginPerDim[d] = pbcMargin * pbc.hbox_diag[d];
        }
    }
    else
    {
        /* Check the total distance along the relevant dimensions */
        for (int d = 0; d < pbc.ndim_ePBC; d++)
        {
            if (dimUsesPbc[d])
            {
                marginDistance2 += pbcMargin * gmx::square(0.5) * norm2(pbc.box[d]);
            }
        }
    }

    auto localAtomIndices = group.atomSet_.localIndex();
    for (gmx::Index indexInSet = 0; indexInSet < localAtomIndices.ssize(); indexInSet++)
    {
        rvec dx;
        pbc_dx(&pbc, x[localAtomIndices[indexInSet]], x_pbc, dx);

        bool atomIsTooFar = false;
        if (pbcIsRectangular)
        {
            for (int d = 0; d < pbc.ndim_ePBC; d++)
            {
                if (dimUsesPbc[d] && (dx[d] < -marginPerDim[d] || dx[d] > marginPerDim[d]))
                {
                    atomIsTooFar = true;
                }
            }
        }
        else
        {
            real pbcDistance2 = 0;
            for (int d = 0; d < pbc.ndim_ePBC; d++)
            {
                if (dimUsesPbc[d])
                {
                    pbcDistance2 += gmx::square(dx[d]);
                }
            }
            atomIsTooFar = (pbcDistance2 > marginDistance2);
        }
        if (atomIsTooFar)
        {
            return false;
        }
    }

    return true;
}

int pullCheckPbcWithinGroups(const pull_t& pull, ArrayRef<const RVec> x, const t_pbc& pbc, real pbcMargin)
{
    if (pbc.pbcType == PbcType::No)
    {
        return -1;
    }

    /* Determine what dimensions are used for each group by pull coordinates */
    std::vector<BoolVec> dimUsed(pull.group.size(), { false, false, false });
    for (const pull_coord_work_t& pcrd : pull.coord)
    {
        const t_pull_coord& coordParams = pcrd.params_;
        for (int groupIndex = 0; groupIndex < coordParams.ngroup; groupIndex++)
        {
            for (int d = 0; d < DIM; d++)
            {
                if (coordParams.dim[d]
                    && !(coordParams.eGeom == PullGroupGeometry::Cylinder && groupIndex == 0))
                {
                    dimUsed[coordParams.group[groupIndex]][d] = true;
                }
            }
        }
    }

    /* Check PBC for every group that uses a PBC reference atom treatment */
    for (size_t g = 0; g < pull.group.size(); g++)
    {
        const pull_group_work_t& group = pull.group[g];
        if ((group.epgrppbc == epgrppbcREFAT || group.epgrppbc == epgrppbcPREVSTEPCOM)
            && !pullGroupObeysPbcRestrictions(group, dimUsed[g], x, pbc, pull.comm.pbcAtomBuffer[g], pbcMargin))
        {
            return g;
        }
    }

    return -1;
}

bool pullCheckPbcWithinGroup(const pull_t& pull, ArrayRef<const RVec> x, const t_pbc& pbc, int groupNr, real pbcMargin)
{
    if (pbc.pbcType == PbcType::No)
    {
        return true;
    }
    GMX_ASSERT(groupNr < gmx::ssize(pull.group), "groupNr is out of range");

    /* Check PBC if the group uses a PBC reference atom treatment. */
    const pull_group_work_t& group = pull.group[groupNr];
    if (group.epgrppbc != epgrppbcREFAT && group.epgrppbc != epgrppbcPREVSTEPCOM)
    {
        return true;
    }

    /* Determine what dimensions are used for each group by pull coordinates */
    BoolVec dimUsed = { false, false, false };
    for (const pull_coord_work_t& pcrd : pull.coord)
    {
        const t_pull_coord& coordParams = pcrd.params_;
        for (int groupIndex = 0; groupIndex < coordParams.ngroup; groupIndex++)
        {
            if (coordParams.group[groupIndex] == groupNr)
            {
                for (int d = 0; d < DIM; d++)
                {
                    if (coordParams.dim[d]
                        && !(coordParams.eGeom == PullGroupGeometry::Cylinder && groupIndex == 0))
                    {
                        dimUsed[d] = true;
                    }
                }
            }
        }
    }

    return (pullGroupObeysPbcRestrictions(
            group, dimUsed, x, pbc, pull.comm.pbcAtomBuffer[groupNr], pbcMargin));
}

void setPrevStepPullComFromState(struct pull_t* pull, const t_state* state)
{
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        for (int j = 0; j < DIM; j++)
        {
            pull->group[g].x_prev_step[j] = state->pull_com_prev_step[g * DIM + j];
        }
    }
}

/*! \brief Whether pull functions save a backup to the t_state object
 *
 * Saving to the state object is only used for checkpointing in the legacy simulator.
 * Modular simulator doesn't use the t_state object for checkpointing.
 */
enum class PullBackupCOM
{
    Yes, //<! Save a copy of the previous step COM to state
    No,  //<! Don't save a copy of the previous step COM to state
};

/*! \brief Sets the previous step COM in pull to the current COM and optionally
 *         stores it in the provided ArrayRef
 *
 * \tparam     pullBackupToState  Whether we're storing the previous COM to state
 * \param[in]  pull  The COM pull force calculation data structure
 * \param[in]   comPreviousStep  The COM of the previous step of each pull group
 */
template<PullBackupCOM pullBackupToState>
static void updatePrevStepPullComImpl(pull_t* pull, gmx::ArrayRef<double> comPreviousStep)
{
    for (gmx::Index g = 0; g < gmx::ssize(pull->group); g++)
    {
        if (pull->group[g].needToCalcCom)
        {
            for (int j = 0; j < DIM; j++)
            {
                pull->group[g].x_prev_step[j] = pull->group[g].x[j];
                if (pullBackupToState == PullBackupCOM::Yes)
                {
                    comPreviousStep[g * DIM + j] = pull->group[g].x[j];
                }
            }
        }
    }
}

void updatePrevStepPullCom(pull_t* pull, std::optional<gmx::ArrayRef<double>> comPreviousStep)
{
    if (comPreviousStep.has_value())
    {
        updatePrevStepPullComImpl<PullBackupCOM::Yes>(pull, comPreviousStep.value());
    }
    else
    {
        updatePrevStepPullComImpl<PullBackupCOM::No>(pull, gmx::ArrayRef<double>());
    }
}

std::vector<double> prevStepPullCom(const pull_t* pull)
{
    std::vector<double> pullCom(pull->group.size() * DIM, 0.0);
    for (gmx::Index g = 0; g < gmx::ssize(pull->group); g++)
    {
        for (int j = 0; j < DIM; j++)
        {
            pullCom[g * DIM + j] = pull->group[g].x_prev_step[j];
        }
    }
    return pullCom;
}

void setPrevStepPullCom(pull_t* pull, gmx::ArrayRef<const double> prevStepPullCom)
{
    GMX_RELEASE_ASSERT(prevStepPullCom.size() >= pull->group.size() * DIM,
                       "Pull COM vector size mismatch.");
    for (gmx::Index g = 0; g < gmx::ssize(pull->group); g++)
    {
        for (int j = 0; j < DIM; j++)
        {
            pull->group[g].x_prev_step[j] = prevStepPullCom[g * DIM + j];
        }
    }
}

void allocStatePrevStepPullCom(t_state* state, const pull_t* pull)
{
    if (!pull)
    {
        state->pull_com_prev_step.clear();
        return;
    }
    size_t ngroup = pull->group.size();
    if (state->pull_com_prev_step.size() / DIM != ngroup)
    {
        state->pull_com_prev_step.resize(ngroup * DIM, NAN);
    }
}

void initPullComFromPrevStep(const t_commrec*     cr,
                             pull_t*              pull,
                             ArrayRef<const real> masses,
                             const t_pbc&         pbc,
                             ArrayRef<const RVec> x)
{
    pull_comm_t* comm   = &pull->comm;
    size_t       ngroup = pull->group.size();

    if (!comm->bParticipate)
    {
        return;
    }

    GMX_ASSERT(comm->pbcAtomBuffer.size() == pull->group.size(),
               "pbcAtomBuffer should have size number of groups");
    GMX_ASSERT(comm->comBuffer.size() == pull->group.size() * c_comBufferStride,
               "comBuffer should have size #group*c_comBufferStride");

    pull_set_pbcatoms(cr, pull, x, comm->pbcAtomBuffer);

    for (size_t g = 0; g < ngroup; g++)
    {
        pull_group_work_t* pgrp;

        pgrp = &pull->group[g];

        if (pgrp->needToCalcCom && pgrp->epgrppbc == epgrppbcPREVSTEPCOM)
        {
            GMX_ASSERT(pgrp->params_.ind.size() > 1,
                       "Groups with no atoms, or only one atom, should not "
                       "use the COM from the previous step as reference.");

            RVec x_pbc = { 0, 0, 0 };
            copy_rvec(comm->pbcAtomBuffer[g], x_pbc);

            /* The following is to a large extent similar to pull_calc_coms() */

            /* The final sums should end up in sum_com[0] */
            ComSums& comSumsTotal = pull->comSums[0];

            if (pgrp->atomSet_.numAtomsLocal() <= c_pullMaxNumLocalAtomsSingleThreaded)
            {
                sum_com_part(pgrp, 0, pgrp->atomSet_.numAtomsLocal(), x, {}, masses, pbc, x_pbc, &comSumsTotal);
            }
            else
            {
                const int numThreads = pgrp->numThreads();
#pragma omp parallel for num_threads(numThreads) schedule(static)
                for (int t = 0; t < numThreads; t++)
                {
                    int ind_start = (pgrp->atomSet_.numAtomsLocal() * (t + 0)) / numThreads;
                    int ind_end   = (pgrp->atomSet_.numAtomsLocal() * (t + 1)) / numThreads;
                    sum_com_part(pgrp, ind_start, ind_end, x, {}, masses, pbc, x_pbc, &pull->comSums[t]);
                }

                /* Reduce the thread contributions to sum_com[0] */
                for (int t = 1; t < numThreads; t++)
                {
                    comSumsTotal.sum_wm += pull->comSums[t].sum_wm;
                    comSumsTotal.sum_wwm += pull->comSums[t].sum_wwm;
                    dvec_inc(comSumsTotal.sum_wmx, pull->comSums[t].sum_wmx);
                    dvec_inc(comSumsTotal.sum_wmxp, pull->comSums[t].sum_wmxp);
                }
            }

            if (pgrp->localWeights.empty())
            {
                comSumsTotal.sum_wwm = comSumsTotal.sum_wm;
            }

            /* Copy local sums to a buffer for global summing */
            auto localSums = gmx::arrayRefFromArray(comm->comBuffer.data() + g * c_comBufferStride,
                                                    c_comBufferStride);

            localSums[0]    = comSumsTotal.sum_wmx;
            localSums[1]    = comSumsTotal.sum_wmxp;
            localSums[2][0] = comSumsTotal.sum_wm;
            localSums[2][1] = comSumsTotal.sum_wwm;
            localSums[2][2] = 0;
        }
    }

    pullAllReduce(cr, comm, ngroup * c_comBufferStride * DIM, static_cast<double*>(comm->comBuffer[0]));

    for (size_t g = 0; g < ngroup; g++)
    {
        pull_group_work_t* pgrp;

        pgrp = &pull->group[g];
        if (pgrp->needToCalcCom)
        {
            if (pgrp->epgrppbc == epgrppbcPREVSTEPCOM)
            {
                auto localSums = gmx::arrayRefFromArray(
                        comm->comBuffer.data() + g * c_comBufferStride, c_comBufferStride);
                double wmass, wwmass;

                /* Determine the inverse mass */
                wmass         = localSums[2][0];
                wwmass        = localSums[2][1];
                pgrp->mwscale = 1.0 / wmass;
                /* invtm==0 signals a frozen group, so then we should keep it zero */
                if (pgrp->invtm != 0)
                {
                    pgrp->wscale = wmass / wwmass;
                    pgrp->invtm  = wwmass / (wmass * wmass);
                }
                /* Divide by the total mass */
                for (int m = 0; m < DIM; m++)
                {
                    pgrp->x[m] = localSums[0][m] * pgrp->mwscale;
                    pgrp->x[m] += comm->pbcAtomBuffer[g][m];
                }
                copy_dvec(pgrp->x, pgrp->x_prev_step);
            }
        }
    }
}
