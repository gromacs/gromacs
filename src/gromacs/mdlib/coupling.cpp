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
#include "gmxpre.h"

#include <assert.h>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/sim_util.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/gammadistribution.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/trajectory/energy.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define NTROTTERPARTS 3

/* Suzuki-Yoshida Constants, for n=3 and n=5, for symplectic integration  */
/* for n=1, w0 = 1 */
/* for n=3, w0 = w2 = 1/(2-2^-(1/3)), w1 = 1-2*w0 */
/* for n=5, w0 = w1 = w3 = w4 = 1/(4-4^-(1/3)), w1 = 1-4*w0 */

#define MAX_SUZUKI_YOSHIDA_NUM 5
#define SUZUKI_YOSHIDA_NUM  5

static const double  sy_const_1[] = { 1. };
static const double  sy_const_3[] = { 0.828981543588751, -0.657963087177502, 0.828981543588751 };
static const double  sy_const_5[] = { 0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065, 0.2967324292201065 };

static const double* sy_const[] = {
    nullptr,
    sy_const_1,
    nullptr,
    sy_const_3,
    nullptr,
    sy_const_5
};

/*
   static const double sy_const[MAX_SUZUKI_YOSHIDA_NUM+1][MAX_SUZUKI_YOSHIDA_NUM+1] = {
    {},
    {1},
    {},
    {0.828981543588751,-0.657963087177502,0.828981543588751},
    {},
    {0.2967324292201065,0.2967324292201065,-0.186929716880426,0.2967324292201065,0.2967324292201065}
   };*/

/* these integration routines are only referenced inside this file */
static void NHC_trotter(t_grpopts *opts, int nvar, gmx_ekindata_t *ekind, real dtfull,
                        double xi[], double vxi[], double scalefac[], real *veta, t_extmass *MassQ, gmx_bool bEkinAveVel)

{
    /* general routine for both barostat and thermostat nose hoover chains */

    int           i, j, mi, mj;
    double        Ekin, Efac, reft, kT, nd;
    double        dt;
    t_grp_tcstat *tcstat;
    double       *ivxi, *ixi;
    double       *iQinv;
    double       *GQ;
    gmx_bool      bBarostat;
    int           mstepsi, mstepsj;
    int           ns = SUZUKI_YOSHIDA_NUM; /* set the degree of integration in the types/state.h file */
    int           nh = opts->nhchainlength;

    snew(GQ, nh);
    mstepsi = mstepsj = ns;

/* if scalefac is NULL, we are doing the NHC of the barostat */

    bBarostat = FALSE;
    if (scalefac == nullptr)
    {
        bBarostat = TRUE;
    }

    for (i = 0; i < nvar; i++)
    {

        /* make it easier to iterate by selecting
           out the sub-array that corresponds to this T group */

        ivxi = &vxi[i*nh];
        ixi  = &xi[i*nh];
        if (bBarostat)
        {
            iQinv = &(MassQ->QPinv[i*nh]);
            nd    = 1.0; /* THIS WILL CHANGE IF NOT ISOTROPIC */
            reft  = std::max<real>(0, opts->ref_t[0]);
            Ekin  = gmx::square(*veta)/MassQ->Winv;
        }
        else
        {
            iQinv  = &(MassQ->Qinv[i*nh]);
            tcstat = &ekind->tcstat[i];
            nd     = opts->nrdf[i];
            reft   = std::max<real>(0, opts->ref_t[i]);
            if (bEkinAveVel)
            {
                Ekin = 2*trace(tcstat->ekinf)*tcstat->ekinscalef_nhc;
            }
            else
            {
                Ekin = 2*trace(tcstat->ekinh)*tcstat->ekinscaleh_nhc;
            }
        }
        kT = BOLTZ*reft;

        for (mi = 0; mi < mstepsi; mi++)
        {
            for (mj = 0; mj < mstepsj; mj++)
            {
                /* weighting for this step using Suzuki-Yoshida integration - fixed at 5 */
                dt = sy_const[ns][mj] * dtfull / mstepsi;

                /* compute the thermal forces */
                GQ[0] = iQinv[0]*(Ekin - nd*kT);

                for (j = 0; j < nh-1; j++)
                {
                    if (iQinv[j+1] > 0)
                    {
                        /* we actually don't need to update here if we save the
                           state of the GQ, but it's easier to just recompute*/
                        GQ[j+1] = iQinv[j+1]*((gmx::square(ivxi[j])/iQinv[j])-kT);
                    }
                    else
                    {
                        GQ[j+1] = 0;
                    }
                }

                ivxi[nh-1] += 0.25*dt*GQ[nh-1];
                for (j = nh-1; j > 0; j--)
                {
                    Efac      = exp(-0.125*dt*ivxi[j]);
                    ivxi[j-1] = Efac*(ivxi[j-1]*Efac + 0.25*dt*GQ[j-1]);
                }

                Efac = exp(-0.5*dt*ivxi[0]);
                if (bBarostat)
                {
                    *veta *= Efac;
                }
                else
                {
                    scalefac[i] *= Efac;
                }
                Ekin *= (Efac*Efac);

                /* Issue - if the KE is an average of the last and the current temperatures, then we might not be
                   able to scale the kinetic energy directly with this factor.  Might take more bookkeeping -- have to
                   think about this a bit more . . . */

                GQ[0] = iQinv[0]*(Ekin - nd*kT);

                /* update thermostat positions */
                for (j = 0; j < nh; j++)
                {
                    ixi[j] += 0.5*dt*ivxi[j];
                }

                for (j = 0; j < nh-1; j++)
                {
                    Efac    = exp(-0.125*dt*ivxi[j+1]);
                    ivxi[j] = Efac*(ivxi[j]*Efac + 0.25*dt*GQ[j]);
                    if (iQinv[j+1] > 0)
                    {
                        GQ[j+1] = iQinv[j+1]*((gmx::square(ivxi[j])/iQinv[j])-kT);
                    }
                    else
                    {
                        GQ[j+1] = 0;
                    }
                }
                ivxi[nh-1] += 0.25*dt*GQ[nh-1];
            }
        }
    }
    sfree(GQ);
}

static void boxv_trotter(t_inputrec *ir, real *veta, real dt, tensor box,
                         gmx_ekindata_t *ekind, tensor vir, real pcorr, t_extmass *MassQ)
{

    real   pscal;
    double alpha;
    int    nwall;
    real   GW, vol;
    tensor ekinmod, localpres;

    /* The heat bath is coupled to a separate barostat, the last temperature group.  In the
       2006 Tuckerman et al paper., the order is iL_{T_baro} iL {T_part}
     */

    if (ir->epct == epctSEMIISOTROPIC)
    {
        nwall = 2;
    }
    else
    {
        nwall = 3;
    }

    /* eta is in pure units.  veta is in units of ps^-1. GW is in
       units of ps^-2.  However, eta has a reference of 1 nm^3, so care must be
       taken to use only RATIOS of eta in updating the volume. */

    /* we take the partial pressure tensors, modify the
       kinetic energy tensor, and recovert to pressure */

    if (ir->opts.nrdf[0] == 0)
    {
        gmx_fatal(FARGS, "Barostat is coupled to a T-group with no degrees of freedom\n");
    }
    /* alpha factor for phase space volume, then multiply by the ekin scaling factor.  */
    alpha  = 1.0 + DIM/((double)ir->opts.nrdf[0]);
    alpha *= ekind->tcstat[0].ekinscalef_nhc;
    msmul(ekind->ekin, alpha, ekinmod);
    /* for now, we use Elr = 0, because if you want to get it right, you
       really should be using PME. Maybe print a warning? */

    pscal   = calc_pres(ir->ePBC, nwall, box, ekinmod, vir, localpres)+pcorr;

    vol = det(box);
    GW  = (vol*(MassQ->Winv/PRESFAC))*(DIM*pscal - trace(ir->ref_p));  /* W is in ps^2 * bar * nm^3 */

    *veta += 0.5*dt*GW;
}

/*
 * This file implements temperature and pressure coupling algorithms:
 * For now only the Weak coupling and the modified weak coupling.
 *
 * Furthermore computation of pressure and temperature is done here
 *
 */

real calc_pres(int ePBC, int nwall, matrix box, tensor ekin, tensor vir,
               tensor pres)
{
    int  n, m;
    real fac;

    if (ePBC == epbcNONE || (ePBC == epbcXY && nwall != 2))
    {
        clear_mat(pres);
    }
    else
    {
        /* Uitzoeken welke ekin hier van toepassing is, zie Evans & Morris - E.
         * Wrs. moet de druktensor gecorrigeerd worden voor de netto stroom in
         * het systeem...
         */

        fac = PRESFAC*2.0/det(box);
        for (n = 0; (n < DIM); n++)
        {
            for (m = 0; (m < DIM); m++)
            {
                pres[n][m] = (ekin[n][m] - vir[n][m])*fac;
            }
        }

        if (debug)
        {
            pr_rvecs(debug, 0, "PC: pres", pres, DIM);
            pr_rvecs(debug, 0, "PC: ekin", ekin, DIM);
            pr_rvecs(debug, 0, "PC: vir ", vir, DIM);
            pr_rvecs(debug, 0, "PC: box ", box, DIM);
        }
    }
    return trace(pres)/DIM;
}

real calc_temp(real ekin, real nrdf)
{
    if (nrdf > 0)
    {
        return (2.0*ekin)/(nrdf*BOLTZ);
    }
    else
    {
        return 0;
    }
}

/*! \brief Sets 1/mass for Parrinello-Rahman in wInv; NOTE: PRESFAC is not included, so not in GROMACS units! */
static void calcParrinelloRahmanInvMass(const t_inputrec *ir, const matrix box,
                                        tensor wInv)
{
    real maxBoxLength;

    /* TODO: See if we can make the mass independent of the box size */
    maxBoxLength = std::max(box[XX][XX], box[YY][YY]);
    maxBoxLength = std::max(maxBoxLength, box[ZZ][ZZ]);

    for (int d = 0; d < DIM; d++)
    {
        for (int n = 0; n < DIM; n++)
        {
            wInv[d][n] = (4*M_PI*M_PI*ir->compress[d][n])/(3*ir->tau_p*ir->tau_p*maxBoxLength);
        }
    }
}

void parrinellorahman_pcoupl(FILE *fplog, gmx_int64_t step,
                             const t_inputrec *ir, real dt, const tensor pres,
                             tensor box, tensor box_rel, tensor boxv,
                             tensor M, matrix mu, gmx_bool bFirstStep)
{
    /* This doesn't do any coordinate updating. It just
     * integrates the box vector equations from the calculated
     * acceleration due to pressure difference. We also compute
     * the tensor M which is used in update to couple the particle
     * coordinates to the box vectors.
     *
     * In Nose and Klein (Mol.Phys 50 (1983) no 5., p 1055) this is
     * given as
     *            -1    .           .     -1
     * M_nk = (h')   * (h' * h + h' h) * h
     *
     * with the dots denoting time derivatives and h is the transformation from
     * the scaled frame to the real frame, i.e. the TRANSPOSE of the box.
     * This also goes for the pressure and M tensors - they are transposed relative
     * to ours. Our equation thus becomes:
     *
     *                  -1       .    .           -1
     * M_gmx = M_nk' = b  * (b * b' + b * b') * b'
     *
     * where b is the gromacs box matrix.
     * Our box accelerations are given by
     *   ..                                    ..
     *   b = vol/W inv(box') * (P-ref_P)     (=h')
     */

    real   vol = box[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
    real   atot, arel, change;
    tensor invbox, pdiff, t1, t2;

    gmx::invertBoxMatrix(box, invbox);

    if (!bFirstStep)
    {
        /* Note that PRESFAC does not occur here.
         * The pressure and compressibility always occur as a product,
         * therefore the pressure unit drops out.
         */
        tensor winv;
        calcParrinelloRahmanInvMass(ir, box, winv);

        m_sub(pres, ir->ref_p, pdiff);

        if (ir->epct == epctSURFACETENSION)
        {
            /* Unlike Berendsen coupling it might not be trivial to include a z
             * pressure correction here? On the other hand we don't scale the
             * box momentarily, but change accelerations, so it might not be crucial.
             */
            real xy_pressure = 0.5*(pres[XX][XX]+pres[YY][YY]);
            for (int d = 0; d < ZZ; d++)
            {
                pdiff[d][d] = (xy_pressure-(pres[ZZ][ZZ]-ir->ref_p[d][d]/box[d][d]));
            }
        }

        tmmul(invbox, pdiff, t1);
        /* Move the off-diagonal elements of the 'force' to one side to ensure
         * that we obey the box constraints.
         */
        for (int d = 0; d < DIM; d++)
        {
            for (int n = 0; n < d; n++)
            {
                t1[d][n] += t1[n][d];
                t1[n][d]  = 0;
            }
        }

        switch (ir->epct)
        {
            case epctANISOTROPIC:
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] *= winv[d][n]*vol;
                    }
                }
                break;
            case epctISOTROPIC:
                /* calculate total volume acceleration */
                atot = box[XX][XX]*box[YY][YY]*t1[ZZ][ZZ]+
                    box[XX][XX]*t1[YY][YY]*box[ZZ][ZZ]+
                    t1[XX][XX]*box[YY][YY]*box[ZZ][ZZ];
                arel = atot/(3*vol);
                /* set all RELATIVE box accelerations equal, and maintain total V
                 * change speed */
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] = winv[0][0]*vol*arel*box[d][n];
                    }
                }
                break;
            case epctSEMIISOTROPIC:
            case epctSURFACETENSION:
                /* Note the correction to pdiff above for surftens. coupling  */

                /* calculate total XY volume acceleration */
                atot = box[XX][XX]*t1[YY][YY]+t1[XX][XX]*box[YY][YY];
                arel = atot/(2*box[XX][XX]*box[YY][YY]);
                /* set RELATIVE XY box accelerations equal, and maintain total V
                 * change speed. Dont change the third box vector accelerations */
                for (int d = 0; d < ZZ; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] = winv[d][n]*vol*arel*box[d][n];
                    }
                }
                for (int n = 0; n < DIM; n++)
                {
                    t1[ZZ][n] *= winv[ZZ][n]*vol;
                }
                break;
            default:
                gmx_fatal(FARGS, "Parrinello-Rahman pressure coupling type %s "
                          "not supported yet\n", EPCOUPLTYPETYPE(ir->epct));
                break;
        }

        real maxchange = 0;
        for (int d = 0; d < DIM; d++)
        {
            for (int n = 0; n <= d; n++)
            {
                boxv[d][n] += dt*t1[d][n];

                /* We do NOT update the box vectors themselves here, since
                 * we need them for shifting later. It is instead done last
                 * in the update() routine.
                 */

                /* Calculate the change relative to diagonal elements-
                   since it's perfectly ok for the off-diagonal ones to
                   be zero it doesn't make sense to check the change relative
                   to its current size.
                 */

                change = fabs(dt*boxv[d][n]/box[d][d]);

                if (change > maxchange)
                {
                    maxchange = change;
                }
            }
        }

        if (maxchange > 0.01 && fplog)
        {
            char buf[22];
            fprintf(fplog,
                    "\nStep %s  Warning: Pressure scaling more than 1%%. "
                    "This may mean your system\n is not yet equilibrated. "
                    "Use of Parrinello-Rahman pressure coupling during\n"
                    "equilibration can lead to simulation instability, "
                    "and is discouraged.\n",
                    gmx_step_str(step, buf));
        }
    }

    preserve_box_shape(ir, box_rel, boxv);

    mtmul(boxv, box, t1);   /* t1=boxv * b' */
    mmul(invbox, t1, t2);
    mtmul(t2, invbox, M);

    /* Determine the scaling matrix mu for the coordinates */
    for (int d = 0; d < DIM; d++)
    {
        for (int n = 0; n <= d; n++)
        {
            t1[d][n] = box[d][n] + dt*boxv[d][n];
        }
    }
    preserve_box_shape(ir, box_rel, t1);
    /* t1 is the box at t+dt, determine mu as the relative change */
    mmul_ur0(invbox, t1, mu);
}

void berendsen_pcoupl(FILE *fplog, gmx_int64_t step,
                      const t_inputrec *ir, real dt,
                      const tensor pres, const matrix box,
                      const matrix force_vir, const matrix constraint_vir,
                      matrix mu, double *baros_integral)
{
    int     d, n;
    real    scalar_pressure, xy_pressure, p_corr_z;
    char    buf[STRLEN];

    /*
     *  Calculate the scaling matrix mu
     */
    scalar_pressure = 0;
    xy_pressure     = 0;
    for (d = 0; d < DIM; d++)
    {
        scalar_pressure += pres[d][d]/DIM;
        if (d != ZZ)
        {
            xy_pressure += pres[d][d]/(DIM-1);
        }
    }
    /* Pressure is now in bar, everywhere. */
#define factor(d, m) (ir->compress[d][m]*dt/ir->tau_p)

    /* mu has been changed from pow(1+...,1/3) to 1+.../3, since this is
     * necessary for triclinic scaling
     */
    clear_mat(mu);
    switch (ir->epct)
    {
        case epctISOTROPIC:
            for (d = 0; d < DIM; d++)
            {
                mu[d][d] = 1.0 - factor(d, d)*(ir->ref_p[d][d] - scalar_pressure) /DIM;
            }
            break;
        case epctSEMIISOTROPIC:
            for (d = 0; d < ZZ; d++)
            {
                mu[d][d] = 1.0 - factor(d, d)*(ir->ref_p[d][d]-xy_pressure)/DIM;
            }
            mu[ZZ][ZZ] =
                1.0 - factor(ZZ, ZZ)*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ])/DIM;
            break;
        case epctANISOTROPIC:
            for (d = 0; d < DIM; d++)
            {
                for (n = 0; n < DIM; n++)
                {
                    mu[d][n] = (d == n ? 1.0 : 0.0)
                        -factor(d, n)*(ir->ref_p[d][n] - pres[d][n])/DIM;
                }
            }
            break;
        case epctSURFACETENSION:
            /* ir->ref_p[0/1] is the reference surface-tension times *
             * the number of surfaces                                */
            if (ir->compress[ZZ][ZZ])
            {
                p_corr_z = dt/ir->tau_p*(ir->ref_p[ZZ][ZZ] - pres[ZZ][ZZ]);
            }
            else
            {
                /* when the compressibity is zero, set the pressure correction   *
                 * in the z-direction to zero to get the correct surface tension */
                p_corr_z = 0;
            }
            mu[ZZ][ZZ] = 1.0 - ir->compress[ZZ][ZZ]*p_corr_z;
            for (d = 0; d < DIM-1; d++)
            {
                mu[d][d] = 1.0 + factor(d, d)*(ir->ref_p[d][d]/(mu[ZZ][ZZ]*box[ZZ][ZZ])
                                               - (pres[ZZ][ZZ]+p_corr_z - xy_pressure))/(DIM-1);
            }
            break;
        default:
            gmx_fatal(FARGS, "Berendsen pressure coupling type %s not supported yet\n",
                      EPCOUPLTYPETYPE(ir->epct));
            break;
    }
    /* To fullfill the orientation restrictions on triclinic boxes
     * we will set mu_yx, mu_zx and mu_zy to 0 and correct
     * the other elements of mu to first order.
     */
    mu[YY][XX] += mu[XX][YY];
    mu[ZZ][XX] += mu[XX][ZZ];
    mu[ZZ][YY] += mu[YY][ZZ];
    mu[XX][YY]  = 0;
    mu[XX][ZZ]  = 0;
    mu[YY][ZZ]  = 0;

    /* Keep track of the work the barostat applies on the system.
     * Without constraints force_vir tells us how Epot changes when scaling.
     * With constraints constraint_vir gives us the constraint contribution
     * to both Epot and Ekin. Although we are not scaling velocities, scaling
     * the coordinates leads to scaling of distances involved in constraints.
     * This in turn changes the angular momentum (even if the constrained
     * distances are corrected at the next step). The kinetic component
     * of the constraint virial captures the angular momentum change.
     */
    for (int d = 0; d < DIM; d++)
    {
        for (int n = 0; n <= d; n++)
        {
            *baros_integral -= 2*(mu[d][n] - (n == d ? 1 : 0))*(force_vir[d][n] + constraint_vir[d][n]);
        }
    }

    if (debug)
    {
        pr_rvecs(debug, 0, "PC: pres ", pres, 3);
        pr_rvecs(debug, 0, "PC: mu   ", mu, 3);
    }

    if (mu[XX][XX] < 0.99 || mu[XX][XX] > 1.01 ||
        mu[YY][YY] < 0.99 || mu[YY][YY] > 1.01 ||
        mu[ZZ][ZZ] < 0.99 || mu[ZZ][ZZ] > 1.01)
    {
        char buf2[22];
        sprintf(buf, "\nStep %s  Warning: pressure scaling more than 1%%, "
                "mu: %g %g %g\n",
                gmx_step_str(step, buf2), mu[XX][XX], mu[YY][YY], mu[ZZ][ZZ]);
        if (fplog)
        {
            fprintf(fplog, "%s", buf);
        }
        fprintf(stderr, "%s", buf);
    }
}

void berendsen_pscale(const t_inputrec *ir, const matrix mu,
                      matrix box, matrix box_rel,
                      int start, int nr_atoms,
                      rvec x[], const unsigned short cFREEZE[],
                      t_nrnb *nrnb)
{
    ivec   *nFreeze = ir->opts.nFreeze;
    int     n, d;
    int     nthreads gmx_unused;

#ifndef __clang_analyzer__
    // cppcheck-suppress unreadVariable
    nthreads = gmx_omp_nthreads_get(emntUpdate);
#endif

    /* Scale the positions */
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (n = start; n < start+nr_atoms; n++)
    {
        // Trivial OpenMP region that does not throw
        int g;

        if (cFREEZE == nullptr)
        {
            g = 0;
        }
        else
        {
            g = cFREEZE[n];
        }

        if (!nFreeze[g][XX])
        {
            x[n][XX] = mu[XX][XX]*x[n][XX]+mu[YY][XX]*x[n][YY]+mu[ZZ][XX]*x[n][ZZ];
        }
        if (!nFreeze[g][YY])
        {
            x[n][YY] = mu[YY][YY]*x[n][YY]+mu[ZZ][YY]*x[n][ZZ];
        }
        if (!nFreeze[g][ZZ])
        {
            x[n][ZZ] = mu[ZZ][ZZ]*x[n][ZZ];
        }
    }
    /* compute final boxlengths */
    for (d = 0; d < DIM; d++)
    {
        box[d][XX] = mu[XX][XX]*box[d][XX]+mu[YY][XX]*box[d][YY]+mu[ZZ][XX]*box[d][ZZ];
        box[d][YY] = mu[YY][YY]*box[d][YY]+mu[ZZ][YY]*box[d][ZZ];
        box[d][ZZ] = mu[ZZ][ZZ]*box[d][ZZ];
    }

    preserve_box_shape(ir, box_rel, box);

    /* (un)shifting should NOT be done after this,
     * since the box vectors might have changed
     */
    inc_nrnb(nrnb, eNR_PCOUPL, nr_atoms);
}

void berendsen_tcoupl(const t_inputrec *ir, const gmx_ekindata_t *ekind, real dt,
                      std::vector<double> &therm_integral)
{
    const t_grpopts *opts = &ir->opts;

    for (int i = 0; (i < opts->ngtc); i++)
    {
        real Ek, T;

        if (ir->eI == eiVV)
        {
            Ek = trace(ekind->tcstat[i].ekinf);
            T  = ekind->tcstat[i].T;
        }
        else
        {
            Ek = trace(ekind->tcstat[i].ekinh);
            T  = ekind->tcstat[i].Th;
        }

        if ((opts->tau_t[i] > 0) && (T > 0.0))
        {
            real reft               = std::max<real>(0, opts->ref_t[i]);
            real lll                = std::sqrt(1.0 + (dt/opts->tau_t[i])*(reft/T-1.0));
            ekind->tcstat[i].lambda = std::max<real>(std::min<real>(lll, 1.25), 0.8);
        }
        else
        {
            ekind->tcstat[i].lambda = 1.0;
        }

        /* Keep track of the amount of energy we are adding to the system */
        therm_integral[i] -= (gmx::square(ekind->tcstat[i].lambda) - 1)*Ek;

        if (debug)
        {
            fprintf(debug, "TC: group %d: T: %g, Lambda: %g\n",
                    i, T, ekind->tcstat[i].lambda);
        }
    }
}

void andersen_tcoupl(t_inputrec *ir, gmx_int64_t step,
                     const t_commrec *cr, const t_mdatoms *md, t_state *state, real rate, const gmx_bool *randomize, const real *boltzfac)
{
    const int                                 *gatindex = (DOMAINDECOMP(cr) ? cr->dd->gatindex : nullptr);
    int                                        i;
    int                                        gc = 0;
    gmx::ThreeFry2x64<0>                       rng(ir->andersen_seed, gmx::RandomDomain::Thermostat);
    gmx::UniformRealDistribution<real>         uniformDist;
    gmx::TabulatedNormalDistribution<real, 14> normalDist;

    /* randomize the velocities of the selected particles */

    for (i = 0; i < md->homenr; i++)  /* now loop over the list of atoms */
    {
        int      ng = gatindex ? gatindex[i] : i;
        gmx_bool bRandomize;

        rng.restart(step, ng);

        if (md->cTC)
        {
            gc = md->cTC[i];  /* assign the atom to a temperature group if there are more than one */
        }
        if (randomize[gc])
        {
            if (ir->etc == etcANDERSENMASSIVE)
            {
                /* Randomize particle always */
                bRandomize = TRUE;
            }
            else
            {
                /* Randomize particle probabilistically */
                uniformDist.reset();
                bRandomize = uniformDist(rng) < rate;
            }
            if (bRandomize)
            {
                real scal;
                int  d;

                scal = std::sqrt(boltzfac[gc]*md->invmass[i]);

                normalDist.reset();

                for (d = 0; d < DIM; d++)
                {
                    state->v[i][d] = scal*normalDist(rng);
                }
            }
        }
    }
}


void nosehoover_tcoupl(t_grpopts *opts, gmx_ekindata_t *ekind, real dt,
                       double xi[], double vxi[], t_extmass *MassQ)
{
    int   i;
    real  reft, oldvxi;

    /* note that this routine does not include Nose-hoover chains yet. Should be easy to add. */

    for (i = 0; (i < opts->ngtc); i++)
    {
        reft     = std::max<real>(0, opts->ref_t[i]);
        oldvxi   = vxi[i];
        vxi[i]  += dt*MassQ->Qinv[i]*(ekind->tcstat[i].Th - reft);
        xi[i]   += dt*(oldvxi + vxi[i])*0.5;
    }
}

void trotter_update(t_inputrec *ir, gmx_int64_t step, gmx_ekindata_t *ekind,
                    gmx_enerdata_t *enerd, t_state *state,
                    tensor vir, t_mdatoms *md,
                    t_extmass *MassQ, int **trotter_seqlist, int trotter_seqno)
{

    int             n, i, d, ngtc, gc = 0, t;
    t_grp_tcstat   *tcstat;
    t_grpopts      *opts;
    gmx_int64_t     step_eff;
    real            dt;
    double         *scalefac, dtc;
    int            *trotter_seq;
    rvec            sumv = {0, 0, 0};
    gmx_bool        bCouple;

    if (trotter_seqno <= ettTSEQ2)
    {
        step_eff = step-1;  /* the velocity verlet calls are actually out of order -- the first half step
                               is actually the last half step from the previous step.  Thus the first half step
                               actually corresponds to the n-1 step*/

    }
    else
    {
        step_eff = step;
    }

    bCouple = (ir->nsttcouple == 1 ||
               do_per_step(step_eff+ir->nsttcouple, ir->nsttcouple));

    trotter_seq = trotter_seqlist[trotter_seqno];

    if ((trotter_seq[0] == etrtSKIPALL) || (!bCouple))
    {
        return;
    }
    dtc  = ir->nsttcouple*ir->delta_t; /* This is OK for NPT, because nsttcouple == nstpcouple is enforcesd */
    opts = &(ir->opts);                /* just for ease of referencing */
    ngtc = opts->ngtc;
    assert(ngtc > 0);
    snew(scalefac, opts->ngtc);
    for (i = 0; i < ngtc; i++)
    {
        scalefac[i] = 1;
    }
    /* execute the series of trotter updates specified in the trotterpart array */

    for (i = 0; i < NTROTTERPARTS; i++)
    {
        /* allow for doubled intgrators by doubling dt instead of making 2 calls */
        if ((trotter_seq[i] == etrtBAROV2) || (trotter_seq[i] == etrtBARONHC2) || (trotter_seq[i] == etrtNHC2))
        {
            dt = 2 * dtc;
        }
        else
        {
            dt = dtc;
        }

        switch (trotter_seq[i])
        {
            case etrtBAROV:
            case etrtBAROV2:
                boxv_trotter(ir, &(state->veta), dt, state->box, ekind, vir,
                             enerd->term[F_PDISPCORR], MassQ);
                break;
            case etrtBARONHC:
            case etrtBARONHC2:
                NHC_trotter(opts, state->nnhpres, ekind, dt, state->nhpres_xi.data(),
                            state->nhpres_vxi.data(), nullptr, &(state->veta), MassQ, FALSE);
                break;
            case etrtNHC:
            case etrtNHC2:
                NHC_trotter(opts, opts->ngtc, ekind, dt, state->nosehoover_xi.data(),
                            state->nosehoover_vxi.data(), scalefac, nullptr, MassQ, (ir->eI == eiVV));
                /* need to rescale the kinetic energies and velocities here.  Could
                   scale the velocities later, but we need them scaled in order to
                   produce the correct outputs, so we'll scale them here. */

                for (t = 0; t < ngtc; t++)
                {
                    tcstat                  = &ekind->tcstat[t];
                    tcstat->vscale_nhc      = scalefac[t];
                    tcstat->ekinscaleh_nhc *= (scalefac[t]*scalefac[t]);
                    tcstat->ekinscalef_nhc *= (scalefac[t]*scalefac[t]);
                }
                /* now that we've scaled the groupwise velocities, we can add them up to get the total */
                /* but do we actually need the total? */

                /* modify the velocities as well */
                for (n = 0; n < md->homenr; n++)
                {
                    if (md->cTC) /* does this conditional need to be here? is this always true?*/
                    {
                        gc = md->cTC[n];
                    }
                    for (d = 0; d < DIM; d++)
                    {
                        state->v[n][d] *= scalefac[gc];
                    }

                    if (debug)
                    {
                        for (d = 0; d < DIM; d++)
                        {
                            sumv[d] += (state->v[n][d])/md->invmass[n];
                        }
                    }
                }
                break;
            default:
                break;
        }
    }
    /* check for conserved momentum -- worth looking at this again eventually, but not working right now.*/
#if 0
    if (debug)
    {
        if (bFirstHalf)
        {
            for (d = 0; d < DIM; d++)
            {
                consk[d] = sumv[d]*exp((1 + 1.0/opts->nrdf[0])*((1.0/DIM)*log(det(state->box)/state->vol0)) + state->nosehoover_xi[0]);
            }
            fprintf(debug, "Conserved kappa: %15.8f %15.8f %15.8f\n", consk[0], consk[1], consk[2]);
        }
    }
#endif
    sfree(scalefac);
}


extern void init_npt_masses(t_inputrec *ir, t_state *state, t_extmass *MassQ, gmx_bool bInit)
{
    int           n, i, j, d, ngtc, nh;
    t_grpopts    *opts;
    real          reft, kT, ndj, nd;

    opts    = &(ir->opts); /* just for ease of referencing */
    ngtc    = ir->opts.ngtc;
    nh      = state->nhchainlength;

    if (ir->eI == eiMD)
    {
        if (bInit)
        {
            snew(MassQ->Qinv, ngtc);
        }
        for (i = 0; (i < ngtc); i++)
        {
            if ((opts->tau_t[i] > 0) && (opts->ref_t[i] > 0))
            {
                MassQ->Qinv[i] = 1.0/(gmx::square(opts->tau_t[i]/M_2PI)*opts->ref_t[i]);
            }
            else
            {
                MassQ->Qinv[i] = 0.0;
            }
        }
    }
    else if (EI_VV(ir->eI))
    {
        /* Set pressure variables */

        if (bInit)
        {
            if (state->vol0 == 0)
            {
                state->vol0 = det(state->box);
                /* because we start by defining a fixed
                   compressibility, we need the volume at this
                   compressibility to solve the problem. */
            }
        }

        /* units are nm^3 * ns^2 / (nm^3 * bar / kJ/mol) = kJ/mol  */
        /* Consider evaluating eventually if this the right mass to use.  All are correct, some might be more stable  */
        MassQ->Winv = (PRESFAC*trace(ir->compress)*BOLTZ*opts->ref_t[0])/(DIM*state->vol0*gmx::square(ir->tau_p/M_2PI));
        /* An alternate mass definition, from Tuckerman et al. */
        /* MassQ->Winv = 1.0/(gmx::square(ir->tau_p/M_2PI)*(opts->nrdf[0]+DIM)*BOLTZ*opts->ref_t[0]); */
        for (d = 0; d < DIM; d++)
        {
            for (n = 0; n < DIM; n++)
            {
                MassQ->Winvm[d][n] = PRESFAC*ir->compress[d][n]/(state->vol0*gmx::square(ir->tau_p/M_2PI));
                /* not clear this is correct yet for the anisotropic case. Will need to reevaluate
                   before using MTTK for anisotropic states.*/
            }
        }
        /* Allocate space for thermostat variables */
        if (bInit)
        {
            snew(MassQ->Qinv, ngtc*nh);
        }

        /* now, set temperature variables */
        for (i = 0; i < ngtc; i++)
        {
            if (opts->tau_t[i] > 0 && opts->ref_t[i] > 0 && opts->nrdf[i] > 0)
            {
                reft = std::max<real>(0, opts->ref_t[i]);
                nd   = opts->nrdf[i];
                kT   = BOLTZ*reft;
                for (j = 0; j < nh; j++)
                {
                    if (j == 0)
                    {
                        ndj = nd;
                    }
                    else
                    {
                        ndj = 1;
                    }
                    MassQ->Qinv[i*nh+j]   = 1.0/(gmx::square(opts->tau_t[i]/M_2PI)*ndj*kT);
                }
            }
            else
            {
                for (j = 0; j < nh; j++)
                {
                    MassQ->Qinv[i*nh+j] = 0.0;
                }
            }
        }
    }
}

int **init_npt_vars(t_inputrec *ir, t_state *state, t_extmass *MassQ, gmx_bool bTrotter)
{
    int           i, j, nnhpres, nh;
    t_grpopts    *opts;
    real          bmass, qmass, reft, kT;
    int         **trotter_seq;

    opts    = &(ir->opts); /* just for ease of referencing */
    nnhpres = state->nnhpres;
    nh      = state->nhchainlength;

    if (EI_VV(ir->eI) && (ir->epc == epcMTTK) && (ir->etc != etcNOSEHOOVER))
    {
        gmx_fatal(FARGS, "Cannot do MTTK pressure coupling without Nose-Hoover temperature control");
    }

    init_npt_masses(ir, state, MassQ, TRUE);

    /* first, initialize clear all the trotter calls */
    snew(trotter_seq, ettTSEQMAX);
    for (i = 0; i < ettTSEQMAX; i++)
    {
        snew(trotter_seq[i], NTROTTERPARTS);
        for (j = 0; j < NTROTTERPARTS; j++)
        {
            trotter_seq[i][j] = etrtNONE;
        }
        trotter_seq[i][0] = etrtSKIPALL;
    }

    if (!bTrotter)
    {
        /* no trotter calls, so we never use the values in the array.
         * We access them (so we need to define them, but ignore
         * then.*/

        return trotter_seq;
    }

    /* compute the kinetic energy by using the half step velocities or
     * the kinetic energies, depending on the order of the trotter calls */

    if (ir->eI == eiVV)
    {
        if (inputrecNptTrotter(ir))
        {
            /* This is the complicated version - there are 4 possible calls, depending on ordering.
               We start with the initial one. */
            /* first, a round that estimates veta. */
            trotter_seq[0][0] = etrtBAROV;

            /* trotter_seq[1] is etrtNHC for 1/2 step velocities - leave zero */

            /* The first half trotter update */
            trotter_seq[2][0] = etrtBAROV;
            trotter_seq[2][1] = etrtNHC;
            trotter_seq[2][2] = etrtBARONHC;

            /* The second half trotter update */
            trotter_seq[3][0] = etrtBARONHC;
            trotter_seq[3][1] = etrtNHC;
            trotter_seq[3][2] = etrtBAROV;

            /* trotter_seq[4] is etrtNHC for second 1/2 step velocities - leave zero */

        }
        else if (inputrecNvtTrotter(ir))
        {
            /* This is the easy version - there are only two calls, both the same.
               Otherwise, even easier -- no calls  */
            trotter_seq[2][0] = etrtNHC;
            trotter_seq[3][0] = etrtNHC;
        }
        else if (inputrecNphTrotter(ir))
        {
            /* This is the complicated version - there are 4 possible calls, depending on ordering.
               We start with the initial one. */
            /* first, a round that estimates veta. */
            trotter_seq[0][0] = etrtBAROV;

            /* trotter_seq[1] is etrtNHC for 1/2 step velocities - leave zero */

            /* The first half trotter update */
            trotter_seq[2][0] = etrtBAROV;
            trotter_seq[2][1] = etrtBARONHC;

            /* The second half trotter update */
            trotter_seq[3][0] = etrtBARONHC;
            trotter_seq[3][1] = etrtBAROV;

            /* trotter_seq[4] is etrtNHC for second 1/2 step velocities - leave zero */
        }
    }
    else if (ir->eI == eiVVAK)
    {
        if (inputrecNptTrotter(ir))
        {
            /* This is the complicated version - there are 4 possible calls, depending on ordering.
               We start with the initial one. */
            /* first, a round that estimates veta. */
            trotter_seq[0][0] = etrtBAROV;

            /* The first half trotter update, part 1 -- double update, because it commutes */
            trotter_seq[1][0] = etrtNHC;

            /* The first half trotter update, part 2 */
            trotter_seq[2][0] = etrtBAROV;
            trotter_seq[2][1] = etrtBARONHC;

            /* The second half trotter update, part 1 */
            trotter_seq[3][0] = etrtBARONHC;
            trotter_seq[3][1] = etrtBAROV;

            /* The second half trotter update */
            trotter_seq[4][0] = etrtNHC;
        }
        else if (inputrecNvtTrotter(ir))
        {
            /* This is the easy version - there is only one call, both the same.
               Otherwise, even easier -- no calls  */
            trotter_seq[1][0] = etrtNHC;
            trotter_seq[4][0] = etrtNHC;
        }
        else if (inputrecNphTrotter(ir))
        {
            /* This is the complicated version - there are 4 possible calls, depending on ordering.
               We start with the initial one. */
            /* first, a round that estimates veta. */
            trotter_seq[0][0] = etrtBAROV;

            /* The first half trotter update, part 1 -- leave zero */
            trotter_seq[1][0] = etrtNHC;

            /* The first half trotter update, part 2 */
            trotter_seq[2][0] = etrtBAROV;
            trotter_seq[2][1] = etrtBARONHC;

            /* The second half trotter update, part 1 */
            trotter_seq[3][0] = etrtBARONHC;
            trotter_seq[3][1] = etrtBAROV;

            /* The second half trotter update -- blank for now */
        }
    }

    switch (ir->epct)
    {
        case epctISOTROPIC:
        default:
            bmass = DIM*DIM; /* recommended mass parameters for isotropic barostat */
    }

    snew(MassQ->QPinv, nnhpres*opts->nhchainlength);

    /* barostat temperature */
    if ((ir->tau_p > 0) && (opts->ref_t[0] > 0))
    {
        reft = std::max<real>(0, opts->ref_t[0]);
        kT   = BOLTZ*reft;
        for (i = 0; i < nnhpres; i++)
        {
            for (j = 0; j < nh; j++)
            {
                if (j == 0)
                {
                    qmass = bmass;
                }
                else
                {
                    qmass = 1;
                }
                MassQ->QPinv[i*opts->nhchainlength+j]   = 1.0/(gmx::square(opts->tau_t[0]/M_2PI)*qmass*kT);
            }
        }
    }
    else
    {
        for (i = 0; i < nnhpres; i++)
        {
            for (j = 0; j < nh; j++)
            {
                MassQ->QPinv[i*nh+j] = 0.0;
            }
        }
    }
    return trotter_seq;
}

static real energyNoseHoover(const t_inputrec *ir, const t_state *state, const t_extmass *MassQ)
{
    real energy = 0;

    int  nh     = state->nhchainlength;

    for (int i = 0; i < ir->opts.ngtc; i++)
    {
        const double *ixi   = &state->nosehoover_xi[i*nh];
        const double *ivxi  = &state->nosehoover_vxi[i*nh];
        const double *iQinv = &(MassQ->Qinv[i*nh]);

        int           nd    = ir->opts.nrdf[i];
        real          reft  = std::max<real>(ir->opts.ref_t[i], 0);
        real          kT    = BOLTZ * reft;

        if (nd > 0.0)
        {
            if (inputrecNvtTrotter(ir))
            {
                /* contribution from the thermal momenta of the NH chain */
                for (int j = 0; j < nh; j++)
                {
                    if (iQinv[j] > 0)
                    {
                        energy += 0.5*gmx::square(ivxi[j])/iQinv[j];
                        /* contribution from the thermal variable of the NH chain */
                        int ndj;
                        if (j == 0)
                        {
                            ndj = nd;
                        }
                        else
                        {
                            ndj = 1.0;
                        }
                        energy += ndj*ixi[j]*kT;
                    }
                }
            }
            else  /* Other non Trotter temperature NH control  -- no chains yet. */
            {
                energy += 0.5*BOLTZ*nd*gmx::square(ivxi[0])/iQinv[0];
                energy += nd*ixi[0]*kT;
            }
        }
    }

    return energy;
}

/* Returns the energy from the barostat thermostat chain */
static real energyPressureMTTK(const t_inputrec *ir, const t_state *state, const t_extmass *MassQ)
{
    real energy = 0;

    int  nh     = state->nhchainlength;

    for (int i = 0; i < state->nnhpres; i++)
    {
        /* note -- assumes only one degree of freedom that is thermostatted in barostat */
        real    reft  = std::max<real>(ir->opts.ref_t[0], 0.0); /* using 'System' temperature */
        real    kT    = BOLTZ * reft;

        for (int j = 0; j < nh; j++)
        {
            double iQinv = MassQ->QPinv[i*nh + j];
            if (iQinv > 0)
            {
                energy += 0.5*gmx::square(state->nhpres_vxi[i*nh + j]/iQinv);
                /* contribution from the thermal variable of the NH chain */
                energy += state->nhpres_xi[i*nh + j]*kT;
            }
            if (debug)
            {
                fprintf(debug, "P-T-group: %10d Chain %4d ThermV: %15.8f ThermX: %15.8f", i, j, state->nhpres_vxi[i*nh + j], state->nhpres_xi[i*nh + j]);
            }
        }
    }

    return energy;
}

/* Returns the energy accumulated by the V-rescale or Berendsen thermostat */
static real energyVrescale(const t_inputrec *ir, const t_state *state)
{
    real energy = 0;
    for (int i = 0; i < ir->opts.ngtc; i++)
    {
        energy += state->therm_integral[i];
    }

    return energy;
}

real NPT_energy(const t_inputrec *ir, const t_state *state, const t_extmass *MassQ)
{
    real energyNPT = 0;

    if (ir->epc != epcNO)
    {
        /* Compute the contribution of the pressure to the conserved quantity*/

        real vol  = det(state->box);

        switch (ir->epc)
        {
            case epcPARRINELLORAHMAN:
            {
                /* contribution from the pressure momenta */
                tensor invMass;
                calcParrinelloRahmanInvMass(ir, state->box, invMass);
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        if (invMass[d][n] > 0)
                        {
                            energyNPT += 0.5*gmx::square(state->boxv[d][n])/(invMass[d][n]*PRESFAC);
                        }
                    }
                }

                /* Contribution from the PV term.
                 * Not that with non-zero off-diagonal reference pressures,
                 * i.e. applied shear stresses, there are additional terms.
                 * We don't support this here, since that requires keeping
                 * track of unwrapped box diagonal elements. This case is
                 * excluded in integratorHasConservedEnergyQuantity().
                 */
                energyNPT += vol*trace(ir->ref_p)/(DIM*PRESFAC);
                break;
            }
            case epcMTTK:
                /* contribution from the pressure momenta */
                energyNPT += 0.5*gmx::square(state->veta)/MassQ->Winv;

                /* contribution from the PV term */
                energyNPT += vol*trace(ir->ref_p)/(DIM*PRESFAC);

                if (ir->epc == epcMTTK)
                {
                    /* contribution from the MTTK chain */
                    energyNPT += energyPressureMTTK(ir, state, MassQ);
                }
                break;
            case epcBERENDSEN:
                energyNPT += state->baros_integral;
                break;
            default:
                GMX_RELEASE_ASSERT(false, "Conserved energy quantity for pressure coupling is not handled. A case should be added with either the conserved quantity added or nothing added and an exclusion added to integratorHasConservedEnergyQuantity().");
        }
    }

    switch (ir->etc)
    {
        case etcNO:
            break;
        case etcVRESCALE:
        case etcBERENDSEN:
            energyNPT += energyVrescale(ir, state);
            break;
        case etcNOSEHOOVER:
            energyNPT += energyNoseHoover(ir, state, MassQ);
            break;
        case etcANDERSEN:
        case etcANDERSENMASSIVE:
            // Not supported, excluded in integratorHasConservedEnergyQuantity()
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Conserved energy quantity for temperature coupling is not handled. A case should be added with either the conserved quantity added or nothing added and an exclusion added to integratorHasConservedEnergyQuantity().");
    }

    return energyNPT;
}


static real vrescale_sumnoises(real                            nn,
                               gmx::ThreeFry2x64<>            *rng,
                               gmx::NormalDistribution<real>  *normalDist)
{
/*
 * Returns the sum of nn independent gaussian noises squared
 * (i.e. equivalent to summing the square of the return values
 * of nn calls to a normal distribution).
 */
    const real                     ndeg_tol = 0.0001;
    real                           r;
    gmx::GammaDistribution<real>   gammaDist(0.5*nn, 1.0);

    if (nn < 2 + ndeg_tol)
    {
        int  nn_int, i;
        real gauss;

        nn_int = (int)(nn + 0.5);

        if (nn - nn_int < -ndeg_tol || nn - nn_int > ndeg_tol)
        {
            gmx_fatal(FARGS, "The v-rescale thermostat was called with a group with #DOF=%f, but for #DOF<3 only integer #DOF are supported", nn + 1);
        }

        r = 0;
        for (i = 0; i < nn_int; i++)
        {
            gauss = (*normalDist)(*rng);
            r    += gauss*gauss;
        }
    }
    else
    {
        /* Use a gamma distribution for any real nn > 2 */
        r = 2.0*gammaDist(*rng);
    }

    return r;
}

static real vrescale_resamplekin(real kk, real sigma, real ndeg, real taut,
                                 gmx_int64_t step, gmx_int64_t seed)
{
/*
 * Generates a new value for the kinetic energy,
 * according to Bussi et al JCP (2007), Eq. (A7)
 * kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
 * sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
 * ndeg:  number of degrees of freedom of the atoms to be thermalized
 * taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
 */
    real                           factor, rr, ekin_new;
    gmx::ThreeFry2x64<64>          rng(seed, gmx::RandomDomain::Thermostat);
    gmx::NormalDistribution<real>  normalDist;

    if (taut > 0.1)
    {
        factor = exp(-1.0/taut);
    }
    else
    {
        factor = 0.0;
    }

    rng.restart(step, 0);

    rr = normalDist(rng);

    ekin_new =
        kk +
        (1.0 - factor)*(sigma*(vrescale_sumnoises(ndeg-1, &rng, &normalDist) + rr*rr)/ndeg - kk) +
        2.0*rr*std::sqrt(kk*sigma/ndeg*(1.0 - factor)*factor);

    return ekin_new;
}

void vrescale_tcoupl(t_inputrec *ir, gmx_int64_t step,
                     gmx_ekindata_t *ekind, real dt,
                     double therm_integral[])
{
    t_grpopts *opts;
    int        i;
    real       Ek, Ek_ref1, Ek_ref, Ek_new;

    opts = &ir->opts;

    for (i = 0; (i < opts->ngtc); i++)
    {
        if (ir->eI == eiVV)
        {
            Ek = trace(ekind->tcstat[i].ekinf);
        }
        else
        {
            Ek = trace(ekind->tcstat[i].ekinh);
        }

        if (opts->tau_t[i] >= 0 && opts->nrdf[i] > 0 && Ek > 0)
        {
            Ek_ref1 = 0.5*opts->ref_t[i]*BOLTZ;
            Ek_ref  = Ek_ref1*opts->nrdf[i];

            Ek_new  = vrescale_resamplekin(Ek, Ek_ref, opts->nrdf[i],
                                           opts->tau_t[i]/dt,
                                           step, ir->ld_seed);

            /* Analytically Ek_new>=0, but we check for rounding errors */
            if (Ek_new <= 0)
            {
                ekind->tcstat[i].lambda = 0.0;
            }
            else
            {
                ekind->tcstat[i].lambda = std::sqrt(Ek_new/Ek);
            }

            therm_integral[i] -= Ek_new - Ek;

            if (debug)
            {
                fprintf(debug, "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n",
                        i, Ek_ref, Ek, Ek_new, ekind->tcstat[i].lambda);
            }
        }
        else
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

void rescale_velocities(gmx_ekindata_t *ekind, t_mdatoms *mdatoms,
                        int start, int end, rvec v[])
{
    t_grp_acc      *gstat;
    t_grp_tcstat   *tcstat;
    unsigned short *cACC, *cTC;
    int             ga, gt, n, d;
    real            lg;
    rvec            vrel;

    tcstat = ekind->tcstat;
    cTC    = mdatoms->cTC;

    if (ekind->bNEMD)
    {
        gstat  = ekind->grpstat;
        cACC   = mdatoms->cACC;

        ga = 0;
        gt = 0;
        for (n = start; n < end; n++)
        {
            if (cACC)
            {
                ga   = cACC[n];
            }
            if (cTC)
            {
                gt   = cTC[n];
            }
            /* Only scale the velocity component relative to the COM velocity */
            rvec_sub(v[n], gstat[ga].u, vrel);
            lg = tcstat[gt].lambda;
            for (d = 0; d < DIM; d++)
            {
                v[n][d] = gstat[ga].u[d] + lg*vrel[d];
            }
        }
    }
    else
    {
        gt = 0;
        for (n = start; n < end; n++)
        {
            if (cTC)
            {
                gt   = cTC[n];
            }
            lg = tcstat[gt].lambda;
            for (d = 0; d < DIM; d++)
            {
                v[n][d] *= lg;
            }
        }
    }
}


/* set target temperatures if we are annealing */
void update_annealing_target_temp(t_inputrec *ir, real t, gmx_update_t *upd)
{
    int  i, j, n, npoints;
    real pert, thist = 0, x;

    for (i = 0; i < ir->opts.ngtc; i++)
    {
        npoints = ir->opts.anneal_npoints[i];
        switch (ir->opts.annealing[i])
        {
            case eannNO:
                continue;
            case  eannPERIODIC:
                /* calculate time modulo the period */
                pert  = ir->opts.anneal_time[i][npoints-1];
                n     = static_cast<int>(t / pert);
                thist = t - n*pert; /* modulo time */
                /* Make sure rounding didn't get us outside the interval */
                if (fabs(thist-pert) < GMX_REAL_EPS*100)
                {
                    thist = 0;
                }
                break;
            case eannSINGLE:
                thist = t;
                break;
            default:
                gmx_fatal(FARGS, "Death horror in update_annealing_target_temp (i=%d/%d npoints=%d)", i, ir->opts.ngtc, npoints);
        }
        /* We are doing annealing for this group if we got here,
         * and we have the (relative) time as thist.
         * calculate target temp */
        j = 0;
        while ((j < npoints-1) && (thist > (ir->opts.anneal_time[i][j+1])))
        {
            j++;
        }
        if (j < npoints-1)
        {
            /* Found our position between points j and j+1.
             * Interpolate: x is the amount from j+1, (1-x) from point j
             * First treat possible jumps in temperature as a special case.
             */
            if ((ir->opts.anneal_time[i][j+1]-ir->opts.anneal_time[i][j]) < GMX_REAL_EPS*100)
            {
                ir->opts.ref_t[i] = ir->opts.anneal_temp[i][j+1];
            }
            else
            {
                x = ((thist-ir->opts.anneal_time[i][j])/
                     (ir->opts.anneal_time[i][j+1]-ir->opts.anneal_time[i][j]));
                ir->opts.ref_t[i] = x*ir->opts.anneal_temp[i][j+1]+(1-x)*ir->opts.anneal_temp[i][j];
            }
        }
        else
        {
            ir->opts.ref_t[i] = ir->opts.anneal_temp[i][npoints-1];
        }
    }

    update_temperature_constants(upd, ir);
}
