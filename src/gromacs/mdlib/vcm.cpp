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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "vcm.h"

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/smalloc.h"

t_vcm *init_vcm(FILE *fp, gmx_groups_t *groups, const t_inputrec *ir)
{
    t_vcm *vcm;
    int    g;

    snew(vcm, 1);

    vcm->mode     = (ir->nstcomm > 0) ? ir->comm_mode : ecmNO;
    vcm->ndim     = ndof_com(ir);
    vcm->timeStep = ir->nstcomm*ir->delta_t;

    if (vcm->mode == ecmANGULAR && vcm->ndim < 3)
    {
        gmx_fatal(FARGS, "Can not have angular comm removal with pbc=%s",
                  epbc_names[ir->ePBC]);
    }

    if (vcm->mode != ecmNO)
    {
        vcm->nr = groups->grps[egcVCM].nr;
        /* Allocate one extra for a possible rest group */
        vcm->size = vcm->nr + 1;
        /* We need vcm->nr+1 elements per thread, but to avoid cache
         * invalidation we add 2 elements to get a 152 byte separation.
         */
        vcm->stride = vcm->nr + 3;
        if (vcm->mode == ecmANGULAR)
        {
            snew(vcm->group_j, vcm->size);
            snew(vcm->group_x, vcm->size);
            snew(vcm->group_i, vcm->size);
            snew(vcm->group_w, vcm->size);
        }
        snew(vcm->group_p, vcm->size);
        snew(vcm->group_v, vcm->size);
        snew(vcm->group_mass, vcm->size);
        snew(vcm->group_name, vcm->size);
        snew(vcm->group_ndf, vcm->size);
        for (g = 0; (g < vcm->nr); g++)
        {
            vcm->group_ndf[g] = ir->opts.nrdf[g];
        }

        /* Copy pointer to group names and print it. */
        if (fp)
        {
            fprintf(fp, "Center of mass motion removal mode is %s\n",
                    ECOM(vcm->mode));
            fprintf(fp, "We have the following groups for center of"
                    " mass motion removal:\n");
        }
        for (g = 0; (g < vcm->nr); g++)
        {
            vcm->group_name[g] = *groups->grpname[groups->grps[egcVCM].nm_ind[g]];
            if (fp)
            {
                fprintf(fp, "%3d:  %s\n", g, vcm->group_name[g]);
            }
        }

        snew(vcm->thread_vcm, gmx_omp_nthreads_get(emntDefault) * vcm->stride);
    }

    return vcm;
}

static void update_tensor(rvec x, real m0, tensor I)
{
    real xy, xz, yz;

    /* Compute inertia tensor contribution due to this atom */
    xy         = x[XX]*x[YY]*m0;
    xz         = x[XX]*x[ZZ]*m0;
    yz         = x[YY]*x[ZZ]*m0;
    I[XX][XX] += x[XX]*x[XX]*m0;
    I[YY][YY] += x[YY]*x[YY]*m0;
    I[ZZ][ZZ] += x[ZZ]*x[ZZ]*m0;
    I[XX][YY] += xy;
    I[YY][XX] += xy;
    I[XX][ZZ] += xz;
    I[ZZ][XX] += xz;
    I[YY][ZZ] += yz;
    I[ZZ][YY] += yz;
}

/* Center of mass code for groups */
void calc_vcm_grp(int start, int homenr, t_mdatoms *md,
                  rvec x[], rvec v[], t_vcm *vcm)
{
    int nthreads = gmx_omp_nthreads_get(emntDefault);
    if (vcm->mode != ecmNO)
    {
#pragma omp parallel num_threads(nthreads)
        {
            int t = gmx_omp_get_thread_num();
            for (int g = 0; g < vcm->size; g++)
            {
                /* Reset linear momentum */
                t_vcm_thread *vcm_t = &vcm->thread_vcm[t*vcm->stride + g];
                vcm_t->mass = 0;
                clear_rvec(vcm_t->p);
                if (vcm->mode == ecmANGULAR)
                {
                    /* Reset angular momentum */
                    clear_rvec(vcm_t->j);
                    clear_rvec(vcm_t->x);
                    clear_mat(vcm_t->i);
                }
            }

#pragma omp for schedule(static)
            for (int i = start; i < start+homenr; i++)
            {
                int  g  = 0;
                real m0 = md->massT[i];
                if (md->cVCM)
                {
                    g = md->cVCM[i];
                }
                t_vcm_thread *vcm_t = &vcm->thread_vcm[t*vcm->stride + g];
                /* Calculate linear momentum */
                vcm_t->mass  += m0;
                int m;
                for (m = 0; (m < DIM); m++)
                {
                    vcm_t->p[m] += m0*v[i][m];
                }

                if (vcm->mode == ecmANGULAR)
                {
                    /* Calculate angular momentum */
                    rvec   j0;
                    cprod(x[i], v[i], j0);

                    for (m = 0; (m < DIM); m++)
                    {
                        vcm_t->j[m] += m0*j0[m];
                        vcm_t->x[m] += m0*x[i][m];
                    }
                    /* Update inertia tensor */
                    update_tensor(x[i], m0, vcm_t->i);
                }
            }
        }
        for (int g = 0; g < vcm->size; g++)
        {
            /* Reset linear momentum */
            vcm->group_mass[g] = 0;
            clear_rvec(vcm->group_p[g]);
            if (vcm->mode == ecmANGULAR)
            {
                /* Reset angular momentum */
                clear_rvec(vcm->group_j[g]);
                clear_rvec(vcm->group_x[g]);
                clear_rvec(vcm->group_w[g]);
                clear_mat(vcm->group_i[g]);
            }

            for (int t = 0; t < nthreads; t++)
            {
                t_vcm_thread *vcm_t = &vcm->thread_vcm[t*vcm->stride + g];
                vcm->group_mass[g] += vcm_t->mass;
                rvec_inc(vcm->group_p[g], vcm_t->p);
                if (vcm->mode == ecmANGULAR)
                {
                    rvec_inc(vcm->group_j[g], vcm_t->j);
                    rvec_inc(vcm->group_x[g], vcm_t->x);
                    m_add(vcm_t->i, vcm->group_i[g], vcm->group_i[g]);
                }
            }
        }

    }
}

/*! \brief Remove the COM motion velocity from the velocities
 *
 * \note This routine should be called from within an OpenMP parallel region.
 *
 * \tparam numDimensions    Correct dimensions 0 to \p numDimensions-1
 * \param[in]     homenr    The number of atoms to correct
 * \param[in]     group_id  List of VCM group ids, when nullptr is passed all atoms are assumed to be in group 0
 * \param[in,out] v         The velocities to correct
 * \param[in]     vcm       VCM data
 */
template<int numDimensions>
static void
doStopComMotionLinear(int                   homenr,
                      const unsigned short *group_id,
                      rvec                 *v,
                      const t_vcm          &vcm)
{
    if (group_id == nullptr)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            for (int d = 0; d < numDimensions; d++)
            {
                v[i][d] -= vcm.group_v[0][d];
            }
        }
    }
    else
    {
#pragma omp for schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            const int g = group_id[i];
            for (int d = 0; d < numDimensions; d++)
            {
                v[i][d] -= vcm.group_v[g][d];
            }
        }
    }
}

/*! \brief Remove the COM motion velocity from the velocities, correct the coordinates assuming constant acceleration
 *
 * \note This routine should be called from within an OpenMP parallel region.
 *
 * \tparam numDimensions    Correct dimensions 0 to \p numDimensions-1
 * \param[in]     homenr    The number of atoms to correct
 * \param[in]     group_id  List of VCM group ids, when nullptr is passed all atoms are assumed to be in group 0
 * \param[in,out] x         The coordinates to correct
 * \param[in,out] v         The velocities to correct
 * \param[in]     vcm       VCM data
 */
template<int numDimensions>
static void
doStopComMotionAccelerationCorrection(int                   homenr,
                                      const unsigned short *group_id,
                                      rvec * gmx_restrict   x,
                                      rvec * gmx_restrict   v,
                                      const t_vcm          &vcm)
{
    const real xCorrectionFactor = 0.5*vcm.timeStep;

    if (group_id == nullptr)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            for (int d = 0; d < numDimensions; d++)
            {
                x[i][d] -= vcm.group_v[0][d]*xCorrectionFactor;
                v[i][d] -= vcm.group_v[0][d];
            }
        }
    }
    else
    {
#pragma omp for schedule(static)
        for (int i = 0; i < homenr; i++)
        {
            const int g = group_id[i];
            for (int d = 0; d < numDimensions; d++)
            {
                x[i][d] -= vcm.group_v[g][d]*xCorrectionFactor;
                v[i][d] -= vcm.group_v[g][d];
            }
        }
    }
}

void do_stopcm_grp(int homenr, const unsigned short *group_id,
                   rvec x[], rvec v[], const t_vcm &vcm)
{
    if (vcm.mode != ecmNO)
    {
        // cppcheck-suppress unreadVariable
        int gmx_unused nth = gmx_omp_nthreads_get(emntDefault);
#pragma omp parallel num_threads(nth)
        {
            if (vcm.mode == ecmLINEAR ||
                vcm.mode == ecmANGULAR ||
                (vcm.mode == ecmLINEAR_ACCELERATION_CORRECTION && x == nullptr))
            {
                /* Subtract linear momentum for v */
                switch (vcm.ndim)
                {
                    case 1:
                        doStopComMotionLinear<1>(homenr, group_id, v, vcm);
                        break;
                    case 2:
                        doStopComMotionLinear<2>(homenr, group_id, v, vcm);
                        break;
                    case 3:
                        doStopComMotionLinear<3>(homenr, group_id, v, vcm);
                        break;
                }
            }
            else
            {
                GMX_ASSERT(vcm.mode == ecmLINEAR_ACCELERATION_CORRECTION, "When the mode is not linear or angular, it should be acceleration correction");
                /* Subtract linear momentum for v and x*/
                switch (vcm.ndim)
                {
                    case 1:
                        doStopComMotionAccelerationCorrection<1>(homenr, group_id, x, v, vcm);
                        break;
                    case 2:
                        doStopComMotionAccelerationCorrection<2>(homenr, group_id, x, v, vcm);
                        break;
                    case 3:
                        doStopComMotionAccelerationCorrection<3>(homenr, group_id, x, v, vcm);
                        break;
                }

            }
            if (vcm.mode == ecmANGULAR)
            {
                /* Subtract angular momentum */
                GMX_ASSERT(x, "Need x to compute angular momentum correction");

                int g = 0;
#pragma omp for schedule(static)
                for (int i = 0; i < homenr; i++)
                {
                    if (group_id)
                    {
                        g = group_id[i];
                    }
                    /* Compute the correction to the velocity for each atom */
                    rvec dv, dx;
                    rvec_sub(x[i], vcm.group_x[g], dx);
                    cprod(vcm.group_w[g], dx, dv);
                    rvec_dec(v[i], dv);
                }
            }
        }
    }
}

static void get_minv(tensor A, tensor B)
{
    int    m, n;
    double fac, rfac;
    tensor tmp;

    tmp[XX][XX] =  A[YY][YY] + A[ZZ][ZZ];
    tmp[YY][XX] = -A[XX][YY];
    tmp[ZZ][XX] = -A[XX][ZZ];
    tmp[XX][YY] = -A[XX][YY];
    tmp[YY][YY] =  A[XX][XX] + A[ZZ][ZZ];
    tmp[ZZ][YY] = -A[YY][ZZ];
    tmp[XX][ZZ] = -A[XX][ZZ];
    tmp[YY][ZZ] = -A[YY][ZZ];
    tmp[ZZ][ZZ] =  A[XX][XX] + A[YY][YY];

    /* This is a hack to prevent very large determinants */
    rfac  = (tmp[XX][XX]+tmp[YY][YY]+tmp[ZZ][ZZ])/3;
    if (rfac == 0.0)
    {
        gmx_fatal(FARGS, "Can not stop center of mass: maybe 2dimensional system");
    }
    fac = 1.0/rfac;
    for (m = 0; (m < DIM); m++)
    {
        for (n = 0; (n < DIM); n++)
        {
            tmp[m][n] *= fac;
        }
    }
    gmx::invertMatrix(tmp, B);
    for (m = 0; (m < DIM); m++)
    {
        for (n = 0; (n < DIM); n++)
        {
            B[m][n] *= fac;
        }
    }
}

void check_cm_grp(FILE *fp, t_vcm *vcm, t_inputrec *ir, real Temp_Max)
{
    int    m, g;
    real   ekcm, ekrot, tm, tm_1, Temp_cm;
    rvec   jcm;
    tensor Icm;

    /* First analyse the total results */
    if (vcm->mode != ecmNO)
    {
        for (g = 0; (g < vcm->nr); g++)
        {
            tm = vcm->group_mass[g];
            if (tm != 0)
            {
                tm_1 = 1.0/tm;
                svmul(tm_1, vcm->group_p[g], vcm->group_v[g]);
            }
            /* Else it's zero anyway! */
        }
        if (vcm->mode == ecmANGULAR)
        {
            for (g = 0; (g < vcm->nr); g++)
            {
                tm = vcm->group_mass[g];
                if (tm != 0)
                {
                    tm_1 = 1.0/tm;

                    /* Compute center of mass for this group */
                    for (m = 0; (m < DIM); m++)
                    {
                        vcm->group_x[g][m] *= tm_1;
                    }

                    /* Subtract the center of mass contribution to the
                     * angular momentum
                     */
                    cprod(vcm->group_x[g], vcm->group_v[g], jcm);
                    for (m = 0; (m < DIM); m++)
                    {
                        vcm->group_j[g][m] -= tm*jcm[m];
                    }

                    /* Subtract the center of mass contribution from the inertia
                     * tensor (this is not as trivial as it seems, but due to
                     * some cancellation we can still do it, even in parallel).
                     */
                    clear_mat(Icm);
                    update_tensor(vcm->group_x[g], tm, Icm);
                    m_sub(vcm->group_i[g], Icm, vcm->group_i[g]);

                    /* Compute angular velocity, using matrix operation
                     * Since J = I w
                     * we have
                     * w = I^-1 J
                     */
                    get_minv(vcm->group_i[g], Icm);
                    mvmul(Icm, vcm->group_j[g], vcm->group_w[g]);
                }
                /* Else it's zero anyway! */
            }
        }
    }
    for (g = 0; (g < vcm->nr); g++)
    {
        ekcm    = 0;
        if (vcm->group_mass[g] != 0 && vcm->group_ndf[g] > 0)
        {
            for (m = 0; m < vcm->ndim; m++)
            {
                ekcm += gmx::square(vcm->group_v[g][m]);
            }
            ekcm   *= 0.5*vcm->group_mass[g];
            Temp_cm = 2*ekcm/vcm->group_ndf[g];

            if ((Temp_cm > Temp_Max) && fp)
            {
                fprintf(fp, "Large VCM(group %s): %12.5f, %12.5f, %12.5f, Temp-cm: %12.5e\n",
                        vcm->group_name[g], vcm->group_v[g][XX],
                        vcm->group_v[g][YY], vcm->group_v[g][ZZ], Temp_cm);
            }

            if (vcm->mode == ecmANGULAR)
            {
                ekrot = 0.5*iprod(vcm->group_j[g], vcm->group_w[g]);
                if ((ekrot > 1) && fp && !EI_RANDOM(ir->eI))
                {
                    /* if we have an integrator that may not conserve momenta, skip */
                    tm    = vcm->group_mass[g];
                    fprintf(fp, "Group %s with mass %12.5e, Ekrot %12.5e Det(I) = %12.5e\n",
                            vcm->group_name[g], tm, ekrot, det(vcm->group_i[g]));
                    fprintf(fp, "  COM: %12.5f  %12.5f  %12.5f\n",
                            vcm->group_x[g][XX], vcm->group_x[g][YY], vcm->group_x[g][ZZ]);
                    fprintf(fp, "  P:   %12.5f  %12.5f  %12.5f\n",
                            vcm->group_p[g][XX], vcm->group_p[g][YY], vcm->group_p[g][ZZ]);
                    fprintf(fp, "  V:   %12.5f  %12.5f  %12.5f\n",
                            vcm->group_v[g][XX], vcm->group_v[g][YY], vcm->group_v[g][ZZ]);
                    fprintf(fp, "  J:   %12.5f  %12.5f  %12.5f\n",
                            vcm->group_j[g][XX], vcm->group_j[g][YY], vcm->group_j[g][ZZ]);
                    fprintf(fp, "  w:   %12.5f  %12.5f  %12.5f\n",
                            vcm->group_w[g][XX], vcm->group_w[g][YY], vcm->group_w[g][ZZ]);
                    pr_rvecs(fp, 0, "Inertia tensor", vcm->group_i[g], DIM);
                }
            }
        }
    }
}
