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

#include "coupling.h"

#include <cassert>
#include <cinttypes>
#include <cmath>

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <string>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/boxmatrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/boxdeformation.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/boxutilities.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/gammadistribution.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#define NTROTTERPARTS 3

/* Suzuki-Yoshida Constants, for n=3 and n=5, for symplectic integration  */
/* for n=1, w0 = 1 */
/* for n=3, w0 = w2 = 1/(2-2^-(1/3)), w1 = 1-2*w0 */
/* for n=5, w0 = w1 = w3 = w4 = 1/(4-4^-(1/3)), w1 = 1-4*w0 */

using gmx::Matrix3x3;

#define SUZUKI_YOSHIDA_NUM 5

static const double sy_const_1[] = { 1. };
static const double sy_const_3[] = { 0.828981543588751, -0.657963087177502, 0.828981543588751 };
static const double sy_const_5[] = { 0.2967324292201065,
                                     0.2967324292201065,
                                     -0.186929716880426,
                                     0.2967324292201065,
                                     0.2967324292201065 };

static constexpr std::array<const double*, 6> sy_const = { nullptr,    sy_const_1, nullptr,
                                                           sy_const_3, nullptr,    sy_const_5 };

static void nosehoover_tcoupl(const gmx_ekindata_t& ekind,
                              real                  dt,
                              gmx::ArrayRef<double> xi,
                              gmx::ArrayRef<double> vxi,
                              const t_extmass&      MassQ);

void update_tcouple(int64_t                             step,
                    const t_inputrec*                   inputrec,
                    t_state*                            state,
                    gmx_ekindata_t*                     ekind,
                    const t_extmass*                    MassQ,
                    int                                 homenr,
                    gmx::ArrayRef<const unsigned short> cTC)

{
    // This condition was explicitly checked in previous version, but should have never been satisfied
    GMX_ASSERT(!(EI_VV(inputrec->eI)
                 && (inputrecNvtTrotter(inputrec) || inputrecNptTrotter(inputrec)
                     || inputrecNphTrotter(inputrec))),
               "Temperature coupling was requested with velocity verlet and trotter");

    bool doTemperatureCoupling = false;

    // For VV temperature coupling parameters are updated on the current
    // step, for the others - one step before.
    if (inputrec->etc == TemperatureCoupling::No)
    {
        doTemperatureCoupling = false;
    }
    else if (EI_VV(inputrec->eI))
    {
        doTemperatureCoupling = do_per_step(step, inputrec->nsttcouple);
    }
    else
    {
        doTemperatureCoupling = do_per_step(step + inputrec->nsttcouple - 1, inputrec->nsttcouple);
    }

    if (doTemperatureCoupling)
    {
        real dttc = inputrec->nsttcouple * inputrec->delta_t;

        // TODO: berendsen_tcoupl(...), nosehoover_tcoupl(...) and vrescale_tcoupl(...) update
        //      temperature coupling parameters, which should be reflected in the name of these
        //      subroutines
        switch (inputrec->etc)
        {
            case TemperatureCoupling::No:
            case TemperatureCoupling::Andersen:
            case TemperatureCoupling::AndersenMassive: break;
            case TemperatureCoupling::Berendsen:
                berendsen_tcoupl(inputrec, ekind, dttc, state->therm_integral);
                break;
            case TemperatureCoupling::NoseHoover:
                nosehoover_tcoupl(*ekind, dttc, state->nosehoover_xi, state->nosehoover_vxi, *MassQ);
                break;
            case TemperatureCoupling::VRescale:
                vrescale_tcoupl(inputrec, step, ekind, dttc, state->therm_integral);
                break;
            default: gmx_fatal(FARGS, "Unknown temperature coupling algorithm");
        }
        /* rescale in place here */
        if (EI_VV(inputrec->eI))
        {
            rescale_velocities(ekind, cTC, 0, homenr, state->v);
        }
    }
    else
    {
        // Set the T scaling lambda to 1 to have no scaling
        // TODO: Do we have to do it on every non-t-couple step?
        for (int i = 0; (i < inputrec->opts.ngtc); i++)
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

void update_pcouple_before_coordinates(const gmx::MDLogger&           mdlog,
                                       int64_t                        step,
                                       const PressureCouplingOptions& pressureCouplingOptions,
                                       const tensor                   deform,
                                       const real                     delta_t,
                                       t_state*                       state,
                                       Matrix3x3*                     parrinellorahmanMu,
                                       Matrix3x3*                     M)
{
    /* Berendsen P-coupling is completely handled after the coordinate update.
     * Trotter P-coupling is handled by separate calls to trotter_update().
     */
    if (pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
        && do_per_step(step + pressureCouplingOptions.nstpcouple - 1, pressureCouplingOptions.nstpcouple))
    {
        const real couplingTimePeriod = pressureCouplingOptions.nstpcouple * delta_t;

        parrinellorahman_pcoupl(mdlog,
                                step,
                                pressureCouplingOptions,
                                deform,
                                couplingTimePeriod,
                                state->pres_prev,
                                state->box,
                                state->box_rel,
                                state->boxv,
                                M,
                                parrinellorahmanMu);
    }
}

void update_pcouple_after_coordinates(FILE*                               fplog,
                                      int64_t                             step,
                                      const PressureCouplingOptions&      pressureCouplingOptions,
                                      const int64_t                       ld_seed,
                                      const real                          ensembleTemperature,
                                      const ivec*                         nFreeze,
                                      const tensor                        deform,
                                      const real                          delta_t,
                                      const int                           homenr,
                                      gmx::ArrayRef<const unsigned short> cFREEZE,
                                      const matrix                        pressure,
                                      const matrix                        forceVirial,
                                      const matrix                        constraintVirial,
                                      Matrix3x3*                          pressureCouplingMu,
                                      t_state*                            state,
                                      t_nrnb*                             nrnb,
                                      gmx::BoxDeformation*                boxDeformation,
                                      const bool                          scaleCoordinates)
{
    int start = 0;

    // Note that delta_t is expressed as a real, which causes no loss
    // in precision in a double-precision build. The only reason for
    // wanting delta_t as a double is to get accurate values for
    // t=delta_t*step when step is larger than float precision. For
    // integration, using delta_t with the accuracy of real suffices,
    // since with integral += delta_t*integrand the increment is
    // nearly always (much) smaller than the integral (and the
    // integrand has real precision).

    // Now update boxes, perhaps after calculating the box-scaling matrix mu
    switch (pressureCouplingOptions.epc)
    {
        case (PressureCoupling::No): break;
        case (PressureCoupling::Berendsen):
            if (do_per_step(step, pressureCouplingOptions.nstpcouple))
            {
                real couplingTimePeriod = pressureCouplingOptions.nstpcouple * delta_t;
                pressureCouplingCalculateScalingMatrix<PressureCoupling::Berendsen>(
                        fplog,
                        step,
                        pressureCouplingOptions,
                        ld_seed,
                        ensembleTemperature,
                        couplingTimePeriod,
                        pressure,
                        state->box,
                        forceVirial,
                        constraintVirial,
                        pressureCouplingMu,
                        &state->baros_integral);
                pressureCouplingScaleBoxAndCoordinates<PressureCoupling::Berendsen>(
                        pressureCouplingOptions,
                        deform,
                        nFreeze,
                        *pressureCouplingMu,
                        state->box,
                        state->box_rel,
                        start,
                        homenr,
                        state->x,
                        gmx::ArrayRef<gmx::RVec>(),
                        cFREEZE,
                        nrnb,
                        scaleCoordinates);
            }
            break;
        case (PressureCoupling::CRescale):
            if (do_per_step(step, pressureCouplingOptions.nstpcouple))
            {
                real couplingTimePeriod = pressureCouplingOptions.nstpcouple * delta_t;
                pressureCouplingCalculateScalingMatrix<PressureCoupling::CRescale>(
                        fplog,
                        step,
                        pressureCouplingOptions,
                        ld_seed,
                        ensembleTemperature,
                        couplingTimePeriod,
                        pressure,
                        state->box,
                        forceVirial,
                        constraintVirial,
                        pressureCouplingMu,
                        &state->baros_integral);
                pressureCouplingScaleBoxAndCoordinates<PressureCoupling::CRescale>(pressureCouplingOptions,
                                                                                   deform,
                                                                                   nFreeze,
                                                                                   *pressureCouplingMu,
                                                                                   state->box,
                                                                                   state->box_rel,
                                                                                   start,
                                                                                   homenr,
                                                                                   state->x,
                                                                                   state->v,
                                                                                   cFREEZE,
                                                                                   nrnb,
                                                                                   scaleCoordinates);
            }
            break;
        case (PressureCoupling::ParrinelloRahman):
            // Note that pressureCouplingMu and the box-velocity
            // update was calculated before the coordinates update.
            if (do_per_step(step + pressureCouplingOptions.nstpcouple - 1, pressureCouplingOptions.nstpcouple))
            {
                /* The box velocities were updated in parrinellorahman_pcoupl,
                 * but we dont change the box vectors until we get here
                 * since we need to be able to shift/unshift above.
                 */
                real couplingTimePeriod = pressureCouplingOptions.nstpcouple * delta_t;
                for (int i = 0; i < DIM; i++)
                {
                    for (int m = 0; m <= i; m++)
                    {
                        state->box[i][m] += couplingTimePeriod * state->boxv[i][m];
                    }
                }
                preserveBoxShape(pressureCouplingOptions, deform, state->box_rel, state->box);

                /* Scale the coordinates */
                if (scaleCoordinates)
                {
                    auto* x = state->x.rvec_array();
                    for (int n = start; n < start + homenr; n++)
                    {
                        if (cFREEZE.empty())
                        {
                            multiplyVectorByTransposeOfBoxMatrix(*pressureCouplingMu, x[n], x[n]);
                        }
                        else
                        {
                            int g = cFREEZE[n];
                            if (!nFreeze[g][XX])
                            {
                                x[n][XX] = (*pressureCouplingMu)(XX, XX) * x[n][XX]
                                           + (*pressureCouplingMu)(YY, XX) * x[n][YY]
                                           + (*pressureCouplingMu)(ZZ, XX) * x[n][ZZ];
                            }
                            if (!nFreeze[g][YY])
                            {
                                x[n][YY] = (*pressureCouplingMu)(YY, YY) * x[n][YY]
                                           + (*pressureCouplingMu)(ZZ, YY) * x[n][ZZ];
                            }
                            if (!nFreeze[g][ZZ])
                            {
                                x[n][ZZ] = (*pressureCouplingMu)(ZZ, ZZ) * x[n][ZZ];
                            }
                        }
                    }
                }
            }
            break;
        case (PressureCoupling::Mttk):
            switch (pressureCouplingOptions.epct)
            {
                case (PressureCouplingType::Isotropic):
                    /* DIM * eta = ln V.  so DIM*eta_new = DIM*eta_old + DIM*dt*veta =>
                       ln V_new = ln V_old + 3*dt*veta => V_new = V_old*exp(3*dt*veta) =>
                       Side length scales as exp(veta*dt) */

                    msmul(state->box, std::exp(state->veta * delta_t), state->box);

                    /* Relate veta to boxv.  veta = d(eta)/dT = (1/DIM)*1/V dV/dT.
                       o               If we assume isotropic scaling, and box length scaling
                       factor L, then V = L^DIM (det(M)).  So dV/dt = DIM
                       L^(DIM-1) dL/dt det(M), and veta = (1/L) dL/dt.  The
                       determinant of B is L^DIM det(M), and the determinant
                       of dB/dt is (dL/dT)^DIM det (M).  veta will be
                       (det(dB/dT)/det(B))^(1/3).  Then since M =
                       B_new*(vol_new)^(1/3), dB/dT_new = (veta_new)*B(new). */

                    msmul(state->box, state->veta, state->boxv);
                    break;
                default: break;
            }
            break;
        default: break;
    }

    if (boxDeformation)
    {
        auto localBox = gmx::createMatrix3x3FromLegacyMatrix(state->box);
        boxDeformation->apply(&localBox, step);
        fillLegacyMatrix(localBox, state->box);
    }
}

extern bool update_randomize_velocities(const t_inputrec*                   ir,
                                        int64_t                             step,
                                        const t_commrec*                    cr,
                                        int                                 homenr,
                                        gmx::ArrayRef<const unsigned short> cTC,
                                        gmx::ArrayRef<const real>           invMass,
                                        gmx::ArrayRef<gmx::RVec>            v,
                                        const gmx::Update*                  upd,
                                        const gmx::Constraints*             constr)
{

    real rate = (ir->delta_t) / ir->opts.tau_t[0];

    if (ir->etc == TemperatureCoupling::Andersen && constr != nullptr)
    {
        /* Currently, Andersen thermostat does not support constrained
           systems. Functionality exists in the andersen_tcoupl
           function in GROMACS 4.5.7 to allow this combination. That
           code could be ported to the current random-number
           generation approach, but has not yet been done because of
           lack of time and resources. */
        gmx_fatal(FARGS,
                  "Normal Andersen is currently not supported with constraints, use massive "
                  "Andersen instead");
    }

    /* proceed with andersen if 1) it's fixed probability per
       particle andersen or 2) it's massive andersen and it's tau_t/dt */
    if ((ir->etc == TemperatureCoupling::Andersen) || do_per_step(step, gmx::roundToInt(1.0 / rate)))
    {
        andersen_tcoupl(
                ir, step, cr, homenr, cTC, invMass, v, rate, upd->getAndersenRandomizeGroup(), upd->getBoltzmanFactor());
        return TRUE;
    }
    return FALSE;
}

/* these integration routines are only referenced inside this file */
static void NHC_trotter(const t_grpopts*      opts,
                        int                   nvar,
                        const gmx_ekindata_t* ekind,
                        real                  dtfull,
                        double                xi[],
                        double                vxi[],
                        double                scalefac[],
                        real*                 veta,
                        const t_extmass*      MassQ,
                        bool                  bEkinAveVel)

{
    /* general routine for both barostat and thermostat nose hoover chains */

    int     i, j, mi, mj;
    double  Ekin, Efac, reft, kT, nd;
    double  dt;
    double *ivxi, *ixi;
    double* GQ;
    bool    bBarostat;
    int     mstepsi, mstepsj;
    int     ns = SUZUKI_YOSHIDA_NUM; /* set the degree of integration in the types/state.h file */
    int     nh = opts->nhchainlength;

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

        ivxi = &vxi[i * nh];
        ixi  = &xi[i * nh];
        gmx::ArrayRef<const double> iQinv;
        if (bBarostat)
        {
            iQinv = gmx::arrayRefFromArray(&MassQ->QPinv[i * nh], nh);
            nd    = 1.0; /* THIS WILL CHANGE IF NOT ISOTROPIC */
            reft  = std::max<real>(0, ekind->currentEnsembleTemperature());
            Ekin  = gmx::square(*veta) / MassQ->Winv;
        }
        else
        {
            iQinv                      = gmx::arrayRefFromArray(&MassQ->Qinv[i * nh], nh);
            const t_grp_tcstat* tcstat = &ekind->tcstat[i];
            nd                         = opts->nrdf[i];
            reft                       = std::max<real>(0, ekind->currentReferenceTemperature(i));
            if (bEkinAveVel)
            {
                Ekin = 2 * trace(tcstat->ekinf) * tcstat->ekinscalef_nhc;
            }
            else
            {
                Ekin = 2 * trace(tcstat->ekinh) * tcstat->ekinscaleh_nhc;
            }
        }
        kT = gmx::c_boltz * reft;

        for (mi = 0; mi < mstepsi; mi++)
        {
            for (mj = 0; mj < mstepsj; mj++)
            {
                /* weighting for this step using Suzuki-Yoshida integration - fixed at 5 */
                dt = sy_const[ns][mj] * dtfull / mstepsi;

                /* compute the thermal forces */
                GQ[0] = iQinv[0] * (Ekin - nd * kT);

                for (j = 0; j < nh - 1; j++)
                {
                    if (iQinv[j + 1] > 0)
                    {
                        /* we actually don't need to update here if we save the
                           state of the GQ, but it's easier to just recompute*/
                        GQ[j + 1] = iQinv[j + 1] * ((gmx::square(ivxi[j]) / iQinv[j]) - kT);
                    }
                    else
                    {
                        GQ[j + 1] = 0;
                    }
                }

                ivxi[nh - 1] += 0.25 * dt * GQ[nh - 1];
                for (j = nh - 1; j > 0; j--)
                {
                    Efac        = std::exp(-0.125 * dt * ivxi[j]);
                    ivxi[j - 1] = Efac * (ivxi[j - 1] * Efac + 0.25 * dt * GQ[j - 1]);
                }

                Efac = std::exp(-0.5 * dt * ivxi[0]);
                if (bBarostat)
                {
                    *veta *= Efac;
                }
                else
                {
                    scalefac[i] *= Efac;
                }
                Ekin *= (Efac * Efac);

                /* Issue - if the KE is an average of the last and the current temperatures, then we
                   might not be able to scale the kinetic energy directly with this factor.  Might
                   take more bookkeeping -- have to think about this a bit more . . . */

                GQ[0] = iQinv[0] * (Ekin - nd * kT);

                /* update thermostat positions */
                for (j = 0; j < nh; j++)
                {
                    ixi[j] += 0.5 * dt * ivxi[j];
                }

                for (j = 0; j < nh - 1; j++)
                {
                    Efac    = std::exp(-0.125 * dt * ivxi[j + 1]);
                    ivxi[j] = Efac * (ivxi[j] * Efac + 0.25 * dt * GQ[j]);
                    if (iQinv[j + 1] > 0)
                    {
                        GQ[j + 1] = iQinv[j + 1] * ((gmx::square(ivxi[j]) / iQinv[j]) - kT);
                    }
                    else
                    {
                        GQ[j + 1] = 0;
                    }
                }
                ivxi[nh - 1] += 0.25 * dt * GQ[nh - 1];
            }
        }
    }
    sfree(GQ);
}

static void boxv_trotter(const t_inputrec*     ir,
                         real*                 veta,
                         real                  delta_t,
                         const tensor          box,
                         const gmx_ekindata_t* ekind,
                         const tensor          vir,
                         const t_extmass*      MassQ)
{

    real   pscal;
    double alpha;
    int    nwall;
    real   GW, vol;
    tensor ekinmod, localpres;

    /* The heat bath is coupled to a separate barostat, the last temperature group.  In the
       2006 Tuckerman et al paper., the order is iL_{T_baro} iL {T_part}
     */

    if (ir->pressureCouplingOptions.epct == PressureCouplingType::SemiIsotropic)
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
    alpha = 1.0 + DIM / (static_cast<double>(ir->opts.nrdf[0]));
    alpha *= ekind->tcstat[0].ekinscalef_nhc;
    msmul(ekind->ekin, alpha, ekinmod);
    /* for now, we use Elr = 0, because if you want to get it right, you
       really should be using PME. Maybe print a warning? */

    pscal = calc_pres(ir->pbcType, nwall, box, ekinmod, vir, localpres);

    vol = det(box);
    GW  = (vol * (MassQ->Winv / gmx::c_presfac))
         * (DIM * pscal - trace(ir->pressureCouplingOptions.ref_p)); /* W is in ps^2 * bar * nm^3 */

    *veta += 0.5 * delta_t * GW;
}

/*
 * This file implements temperature and pressure coupling algorithms:
 * For now only the Weak coupling and the modified weak coupling.
 *
 * Furthermore computation of pressure and temperature is done here
 *
 */

real calc_pres(PbcType pbcType, int nwall, const matrix box, const tensor ekin, const tensor vir, tensor pres)
{
    int  n, m;
    real fac;

    if (pbcType == PbcType::No || (pbcType == PbcType::XY && nwall != 2))
    {
        clear_mat(pres);
    }
    else
    {
        /* Uitzoeken welke ekin hier van toepassing is, zie Evans & Morris - E.
         * Wrs. moet de druktensor gecorrigeerd worden voor de netto stroom in
         * het systeem...
         */

        fac = gmx::c_presfac * 2.0 / det(box);
        for (n = 0; (n < DIM); n++)
        {
            for (m = 0; (m < DIM); m++)
            {
                pres[n][m] = (ekin[n][m] - vir[n][m]) * fac;
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
    return trace(pres) / DIM;
}

real calc_temp(real ekin, real nrdf)
{
    if (nrdf > 0)
    {
        return (2.0 * ekin) / (nrdf * gmx::c_boltz);
    }
    else
    {
        return 0;
    }
}

/*! \brief Sets 1/mass for Parrinello-Rahman in wInv; NOTE: c_presfac is not included, so not in GROMACS units! */
static void calcParrinelloRahmanInvMass(const PressureCouplingOptions& pressureCouplingOptions,
                                        const matrix                   box,
                                        tensor                         wInv)
{
    real maxBoxLength;

    /* TODO: See if we can make the mass independent of the box size */
    maxBoxLength = std::max(box[XX][XX], box[YY][YY]);
    maxBoxLength = std::max(maxBoxLength, box[ZZ][ZZ]);

    for (int d = 0; d < DIM; d++)
    {
        for (int n = 0; n < DIM; n++)
        {
            wInv[d][n] = (4 * M_PI * M_PI * pressureCouplingOptions.compress[d][n])
                         / (3 * pressureCouplingOptions.tau_p * pressureCouplingOptions.tau_p * maxBoxLength);
        }
    }
}

/*! \brief Returns the product of an inverse box matrix with a box or box velocity matrix
 *
 * The result has off-diagonal elements exactly zero when the compressibility
 * matrix has off-diagonal elements that are all exactly zero.
 */
static Matrix3x3 productOfInvBoxAndBoxMatrix(const PressureCouplingOptions& pressureCouplingOptions,
                                             const Matrix3x3&               invbox,
                                             const Matrix3x3&               box)

{
    Matrix3x3 M;

    if (TRICLINIC(pressureCouplingOptions.compress))
    {
        M = multiplyBoxMatrices(invbox, box);
    }
    else
    {
        // Compute only the diagonal elements to avoid non-zero off-diagonal due to rounding
        M = { { 0._real } };
        for (int d = 0; d < DIM; d++)
        {
            M(d, d) = invbox(d, d) * box(d, d);
        }
    }

    return M;
}

/*! \brief Calculate the mu tensor for Parrinello-Rahman pressure coupling
 *
 * mu describes the relative change in box vectors introduced by the
 * coupling and applied at each step. This ensures there is no large
 * jump in box-vector values.
 */
static Matrix3x3 calculateMu(const PressureCouplingOptions& pressureCouplingOptions,
                             const tensor                   deform,
                             tensor                         box_rel,
                             const matrix                   box,
                             const Matrix3x3&               invbox,
                             const matrix                   boxv,
                             const real                     couplingTimePeriod)
{
    tensor temp{ { 0._real } };
    for (int d = 0; d < DIM; d++)
    {
        for (int n = 0; n <= d; n++)
        {
            temp[d][n] = box[d][n] + couplingTimePeriod * boxv[d][n];
        }
    }
    preserveBoxShape(pressureCouplingOptions, deform, box_rel, temp);
    // Now temp is the box at t+dt, determine mu as the relative
    // change in the box.
    return productOfInvBoxAndBoxMatrix(
            pressureCouplingOptions, invbox, gmx::createMatrix3x3FromLegacyMatrix(temp));
}

void init_parrinellorahman(const PressureCouplingOptions& pressureCouplingOptions,
                           const tensor                   deform,
                           const real                     couplingTimePeriod,
                           const tensor                   box,
                           tensor                         box_rel,
                           tensor                         boxv,
                           Matrix3x3*                     M,
                           Matrix3x3*                     mu)
{
    const Matrix3x3 inverseBox = gmx::invertBoxMatrix(gmx::createMatrix3x3FromLegacyMatrix(box));
    preserveBoxShape(pressureCouplingOptions, deform, box_rel, boxv);
    *M = productOfInvBoxAndBoxMatrix(
            pressureCouplingOptions, inverseBox, gmx::createMatrix3x3FromLegacyMatrix(boxv));
    *mu = calculateMu(pressureCouplingOptions, deform, box_rel, box, inverseBox, boxv, couplingTimePeriod);
}

void parrinellorahman_pcoupl(const gmx::MDLogger&           mdlog,
                             int64_t                        step,
                             const PressureCouplingOptions& pressureCouplingOptions,
                             const tensor                   deform,
                             const real                     couplingTimePeriod,
                             const tensor                   pres,
                             const tensor                   box,
                             tensor                         box_rel,
                             tensor                         boxv,
                             Matrix3x3*                     M,
                             Matrix3x3*                     mu)
{
    real   vol = box[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
    real   atot, arel, change;
    tensor invbox, pdiff, t1;

    gmx::invertBoxMatrix(box, invbox);

    // Indentation preserved for review convenience, can be
    // removed later
    {
        /* First, calculate the acceleration of the box vectors.
         *
         * Note that c_presfac does not occur here.
         * The pressure and compressibility always occur as a product,
         * therefore the pressure unit drops out.
         */
        tensor winv;
        calcParrinelloRahmanInvMass(pressureCouplingOptions, box, winv);

        m_sub(pres, pressureCouplingOptions.ref_p, pdiff);

        if (pressureCouplingOptions.epct == PressureCouplingType::SurfaceTension)
        {
            /* Unlike Berendsen coupling it might not be trivial to include a z
             * pressure correction here? On the other hand we don't scale the
             * box momentarily, but change accelerations, so it might not be crucial.
             */
            real xy_pressure = 0.5 * (pres[XX][XX] + pres[YY][YY]);
            for (int d = 0; d < ZZ; d++)
            {
                pdiff[d][d] =
                        (xy_pressure - (pres[ZZ][ZZ] - pressureCouplingOptions.ref_p[d][d] / box[d][d]));
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
                t1[n][d] = 0;
            }
        }

        switch (pressureCouplingOptions.epct)
        {
            case PressureCouplingType::Anisotropic:
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] *= winv[d][n] * vol;
                    }
                }
                break;
            case PressureCouplingType::Isotropic:
                /* calculate total volume acceleration */
                atot = box[XX][XX] * box[YY][YY] * t1[ZZ][ZZ] + box[XX][XX] * t1[YY][YY] * box[ZZ][ZZ]
                       + t1[XX][XX] * box[YY][YY] * box[ZZ][ZZ];
                arel = atot / (3 * vol);
                /* set all RELATIVE box accelerations equal, and maintain total V
                 * change speed */
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] = winv[0][0] * vol * arel * box[d][n];
                    }
                }
                break;
            case PressureCouplingType::SemiIsotropic:
            case PressureCouplingType::SurfaceTension:
                /* Note the correction to pdiff above for surftens. coupling  */

                /* calculate total XY volume acceleration */
                atot = box[XX][XX] * t1[YY][YY] + t1[XX][XX] * box[YY][YY];
                arel = atot / (2 * box[XX][XX] * box[YY][YY]);
                /* set RELATIVE XY box accelerations equal, and maintain total V
                 * change speed. Dont change the third box vector accelerations */
                for (int d = 0; d < ZZ; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        t1[d][n] = winv[d][n] * vol * arel * box[d][n];
                    }
                }
                for (int n = 0; n < DIM; n++)
                {
                    t1[ZZ][n] *= winv[ZZ][n] * vol;
                }
                break;
            default:
                gmx_fatal(FARGS,
                          "Parrinello-Rahman pressure coupling type %s "
                          "not supported yet\n",
                          enumValueToString(pressureCouplingOptions.epct));
        }

        // Update the box velocities from the box accelerations, and
        // prepare to log a warning about large changes, if needed.
        real maxchange = 0;
        for (int d = 0; d < DIM; d++)
        {
            for (int n = 0; n <= d; n++)
            {
                boxv[d][n] += couplingTimePeriod * t1[d][n];

                /* Calculate the change relative to diagonal elements-
                   since it's perfectly ok for the off-diagonal ones to
                   be zero it doesn't make sense to check the change relative
                   to its current size.
                 */

                change = std::fabs(couplingTimePeriod * boxv[d][n] / box[d][d]);

                if (change > maxchange)
                {
                    maxchange = change;
                }
            }
        }

        if (maxchange > 0.01)
        {
            GMX_LOG(mdlog.warning)
                    .asParagraph()
                    .appendTextFormatted("Step %" PRId64
                                         " Pressure scaling more than 1%%. "
                                         "This may mean your system is not yet equilibrated. "
                                         "Use of Parrinello-Rahman pressure coupling during "
                                         "equilibration can lead to simulation instability, "
                                         "and is discouraged.",
                                         step);
        }
    }

    // The new box velocity has been calculated, but might need
    // correcting for box shape, e.g. when rounding error has
    // accumulated.
    preserveBoxShape(pressureCouplingOptions, deform, box_rel, boxv);

    const Matrix3x3 inverseBox = gmx::createMatrix3x3FromLegacyMatrix(invbox);
    *M                         = productOfInvBoxAndBoxMatrix(
            pressureCouplingOptions, inverseBox, gmx::createMatrix3x3FromLegacyMatrix(boxv));

    *mu = calculateMu(pressureCouplingOptions, deform, box_rel, box, inverseBox, boxv, couplingTimePeriod);
}

//! Return compressibility factor for entry (i,j) of Berendsen / C-rescale scaling matrix
static inline real compressibilityFactor(int                            i,
                                         int                            j,
                                         const PressureCouplingOptions& pressureCouplingOptions,
                                         const real                     couplingTimePeriod)
{
    return pressureCouplingOptions.compress[i][j] * couplingTimePeriod / pressureCouplingOptions.tau_p;
}

//! Details of Berendsen / C-rescale scaling matrix calculation
template<PressureCoupling pressureCouplingType>
static void calculateScalingMatrixImplDetail(const PressureCouplingOptions& pressureCouplingOptions,
                                             int64_t                        ld_seed,
                                             real                           ensembleTemperature,
                                             Matrix3x3*                     mu,
                                             real                           couplingTimePeriod,
                                             const matrix                   pres,
                                             const matrix                   box,
                                             real                           scalar_pressure,
                                             real                           xy_pressure,
                                             int64_t                        step);

//! Calculate Berendsen / C-rescale scaling matrix
template<PressureCoupling pressureCouplingType>
static void calculateScalingMatrixImpl(const PressureCouplingOptions& pressureCouplingOptions,
                                       const int64_t                  ld_seed,
                                       const real                     ensembleTemperature,
                                       Matrix3x3*                     mu,
                                       const real                     couplingTimePeriod,
                                       const matrix                   pres,
                                       const matrix                   box,
                                       int64_t                        step)
{
    real scalar_pressure = 0;
    real xy_pressure     = 0;
    for (int d = 0; d < DIM; d++)
    {
        scalar_pressure += pres[d][d] / DIM;
        if (d != ZZ)
        {
            xy_pressure += pres[d][d] / (DIM - 1);
        }
    }
    *mu = gmx::Matrix3x3{ { 0._real } };
    calculateScalingMatrixImplDetail<pressureCouplingType>(pressureCouplingOptions,
                                                           ld_seed,
                                                           ensembleTemperature,
                                                           mu,
                                                           couplingTimePeriod,
                                                           pres,
                                                           box,
                                                           scalar_pressure,
                                                           xy_pressure,
                                                           step);
}

template<>
void calculateScalingMatrixImplDetail<PressureCoupling::Berendsen>(const PressureCouplingOptions& pressureCouplingOptions,
                                                                   int64_t /*ld_seed*/,
                                                                   real /*ensembleTemperature*/,
                                                                   Matrix3x3*   mu,
                                                                   const real   couplingTimePeriod,
                                                                   const matrix pres,
                                                                   const matrix box,
                                                                   real         scalar_pressure,
                                                                   real         xy_pressure,
                                                                   int64_t gmx_unused step)
{
    real p_corr_z = 0;
    switch (pressureCouplingOptions.epct)
    {
        case PressureCouplingType::Isotropic:
            for (int d = 0; d < DIM; d++)
            {
                (*mu)(d, d) = 1.0
                              - compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod)
                                        * (pressureCouplingOptions.ref_p[d][d] - scalar_pressure) / DIM;
            }
            break;
        case PressureCouplingType::SemiIsotropic:
            for (int d = 0; d < ZZ; d++)
            {
                (*mu)(d, d) = 1.0
                              - compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod)
                                        * (pressureCouplingOptions.ref_p[d][d] - xy_pressure) / DIM;
            }
            (*mu)(ZZ, ZZ) = 1.0
                            - compressibilityFactor(ZZ, ZZ, pressureCouplingOptions, couplingTimePeriod)
                                      * (pressureCouplingOptions.ref_p[ZZ][ZZ] - pres[ZZ][ZZ]) / DIM;
            break;
        case PressureCouplingType::Anisotropic:
            for (int d = 0; d < DIM; d++)
            {
                for (int n = 0; n < DIM; n++)
                {
                    (*mu)(d, n) = (d == n ? 1.0 : 0.0)
                                  - compressibilityFactor(d, n, pressureCouplingOptions, couplingTimePeriod)
                                            * (pressureCouplingOptions.ref_p[d][n] - pres[d][n]) / DIM;
                }
            }
            break;
        case PressureCouplingType::SurfaceTension:
            /* pressureCouplingOptions.ref_p[0/1] is the reference surface-tension times *
             * the number of surfaces                                */
            if (pressureCouplingOptions.compress[ZZ][ZZ] != 0.0F)
            {
                p_corr_z = couplingTimePeriod / pressureCouplingOptions.tau_p
                           * (pressureCouplingOptions.ref_p[ZZ][ZZ] - pres[ZZ][ZZ]);
            }
            else
            {
                /* when the compressibity is zero, set the pressure correction   *
                 * in the z-direction to zero to get the correct surface tension */
                p_corr_z = 0;
            }
            (*mu)(ZZ, ZZ) = 1.0 - pressureCouplingOptions.compress[ZZ][ZZ] * p_corr_z;
            for (int d = 0; d < DIM - 1; d++)
            {
                (*mu)(d, d) =
                        1.0
                        + compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod)
                                  * (pressureCouplingOptions.ref_p[d][d] / ((*mu)(ZZ, ZZ) * box[ZZ][ZZ])
                                     - (pres[ZZ][ZZ] + p_corr_z - xy_pressure))
                                  / (DIM - 1);
            }
            break;
        default:
            gmx_fatal(FARGS,
                      "Berendsen pressure coupling type %s not supported yet\n",
                      enumValueToString(pressureCouplingOptions.epct));
    }
}

template<>
void calculateScalingMatrixImplDetail<PressureCoupling::CRescale>(const PressureCouplingOptions& pressureCouplingOptions,
                                                                  const int64_t ld_seed,
                                                                  const real    ensembleTemperature,
                                                                  Matrix3x3*    mu,
                                                                  const real    couplingTimePeriod,
                                                                  const matrix  pres,
                                                                  const matrix  box,
                                                                  real          scalar_pressure,
                                                                  real          xy_pressure,
                                                                  int64_t       step)
{
    gmx::ThreeFry2x64<64>         rng(ld_seed, gmx::RandomDomain::Barostat);
    gmx::NormalDistribution<real> normalDist;
    rng.restart(step, 0);
    real vol = 1.0;
    for (int d = 0; d < DIM; d++)
    {
        vol *= box[d][d];
    }
    real gauss  = 0;
    real gauss2 = 0;
    real kt     = ensembleTemperature * gmx::c_boltz;
    if (kt < 0.0)
    {
        kt = 0.0;
    }

    switch (pressureCouplingOptions.epct)
    {
        case PressureCouplingType::Isotropic:
            gauss = normalDist(rng);
            for (int d = 0; d < DIM; d++)
            {
                const real factor =
                        compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod);
                (*mu)(d, d) =
                        std::exp(-factor * (pressureCouplingOptions.ref_p[d][d] - scalar_pressure) / DIM
                                 + std::sqrt(2.0 * kt * factor * gmx::c_presfac / vol) * gauss / DIM);
            }
            break;
        case PressureCouplingType::SemiIsotropic:
            gauss  = normalDist(rng);
            gauss2 = normalDist(rng);
            for (int d = 0; d < ZZ; d++)
            {
                const real factor =
                        compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod);
                (*mu)(d, d) =
                        std::exp(-factor * (pressureCouplingOptions.ref_p[d][d] - xy_pressure) / DIM
                                 + std::sqrt((DIM - 1) * 2.0 * kt * factor * gmx::c_presfac / vol / DIM)
                                           / (DIM - 1) * gauss);
            }
            {
                const real factor =
                        compressibilityFactor(ZZ, ZZ, pressureCouplingOptions, couplingTimePeriod);
                (*mu)(ZZ, ZZ) =
                        std::exp(-factor * (pressureCouplingOptions.ref_p[ZZ][ZZ] - pres[ZZ][ZZ]) / DIM
                                 + std::sqrt(2.0 * kt * factor * gmx::c_presfac / vol / DIM) * gauss2);
            }
            break;
        case PressureCouplingType::SurfaceTension:
            gauss  = normalDist(rng);
            gauss2 = normalDist(rng);
            for (int d = 0; d < ZZ; d++)
            {
                const real factor =
                        compressibilityFactor(d, d, pressureCouplingOptions, couplingTimePeriod);
                /* Notice: we here use ref_p[ZZ][ZZ] as isotropic pressure and pressureCouplingOptions.ref_p[d][d] as surface tension */
                (*mu)(d, d) = std::exp(
                        -factor
                                * (pressureCouplingOptions.ref_p[ZZ][ZZ]
                                   - pressureCouplingOptions.ref_p[d][d] / box[ZZ][ZZ] - xy_pressure)
                                / DIM
                        + std::sqrt(4.0 / 3.0 * kt * factor * gmx::c_presfac / vol) / (DIM - 1) * gauss);
            }
            {
                const real factor =
                        compressibilityFactor(ZZ, ZZ, pressureCouplingOptions, couplingTimePeriod);
                (*mu)(ZZ, ZZ) =
                        std::exp(-factor * (pressureCouplingOptions.ref_p[ZZ][ZZ] - pres[ZZ][ZZ]) / DIM
                                 + std::sqrt(2.0 / 3.0 * kt * factor * gmx::c_presfac / vol) * gauss2);
            }
            break;
        default:
            gmx_fatal(FARGS,
                      "C-rescale pressure coupling type %s not supported yet\n",
                      enumValueToString(pressureCouplingOptions.epct));
    }
}

template<PressureCoupling pressureCouplingType>
void pressureCouplingCalculateScalingMatrix(FILE*                          fplog,
                                            int64_t                        step,
                                            const PressureCouplingOptions& pressureCouplingOptions,
                                            int64_t                        ld_seed,
                                            real                           ensembleTemperature,
                                            const real                     couplingTimePeriod,
                                            const tensor                   pres,
                                            const matrix                   box,
                                            const matrix                   force_vir,
                                            const matrix                   constraint_vir,
                                            Matrix3x3*                     mu,
                                            double*                        baros_integral)
{
    // If the support here increases, we need to also add the template declarations
    // for new cases below.
    static_assert(pressureCouplingType == PressureCoupling::Berendsen
                          || pressureCouplingType == PressureCoupling::CRescale,
                  "pressureCouplingCalculateScalingMatrix is only implemented for Berendsen and "
                  "C-rescale pressure coupling");

    calculateScalingMatrixImpl<pressureCouplingType>(
            pressureCouplingOptions, ld_seed, ensembleTemperature, mu, couplingTimePeriod, pres, box, step);

    /* To fulfill the orientation restrictions on triclinic boxes
     * we will set mu_yx, mu_zx and mu_zy to 0 and correct
     * the other elements of mu to first order.
     */
    (*mu)(YY, XX) += (*mu)(XX, YY);
    (*mu)(ZZ, XX) += (*mu)(XX, ZZ);
    (*mu)(ZZ, YY) += (*mu)(YY, ZZ);
    (*mu)(XX, YY) = 0;
    (*mu)(XX, ZZ) = 0;
    (*mu)(YY, ZZ) = 0;

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
            *baros_integral -=
                    2 * ((*mu)(d, n) - (n == d ? 1 : 0)) * (force_vir[d][n] + constraint_vir[d][n]);
        }
    }

    if ((*mu)(XX, XX) < 0.99 || (*mu)(XX, XX) > 1.01 || (*mu)(YY, YY) < 0.99 || (*mu)(YY, YY) > 1.01
        || (*mu)(ZZ, ZZ) < 0.99 || (*mu)(ZZ, ZZ) > 1.01)
    {
        char buf[STRLEN];
        char buf2[22];
        sprintf(buf,
                "\nStep %s  Warning: pressure scaling more than 1%%, "
                "mu: %g %g %g\n",
                gmx_step_str(step, buf2),
                (*mu)(XX, XX),
                (*mu)(YY, YY),
                (*mu)(ZZ, ZZ));
        if (fplog)
        {
            fprintf(fplog, "%s", buf);
        }
        fprintf(stderr, "%s", buf);
    }
}

template void pressureCouplingCalculateScalingMatrix<PressureCoupling::CRescale>(FILE*,
                                                                                 int64_t,
                                                                                 const PressureCouplingOptions&,
                                                                                 int64_t,
                                                                                 real,
                                                                                 real,
                                                                                 const tensor,
                                                                                 const matrix,
                                                                                 const matrix,
                                                                                 const matrix,
                                                                                 Matrix3x3*,
                                                                                 double*);

template void pressureCouplingCalculateScalingMatrix<PressureCoupling::Berendsen>(FILE*,
                                                                                  int64_t,
                                                                                  const PressureCouplingOptions&,
                                                                                  int64_t,
                                                                                  real,
                                                                                  real,
                                                                                  const tensor,
                                                                                  const matrix,
                                                                                  const matrix,
                                                                                  const matrix,
                                                                                  Matrix3x3*,
                                                                                  double*);

template<PressureCoupling pressureCouplingType>
void pressureCouplingScaleBoxAndCoordinates(const PressureCouplingOptions& pressureCouplingOptions,
                                            const tensor                   deform,
                                            const ivec*                    nFreeze,
                                            const Matrix3x3&               mu,
                                            matrix                         box,
                                            matrix                         box_rel,
                                            int                            start,
                                            int                            nr_atoms,
                                            gmx::ArrayRef<gmx::RVec>       x,
                                            gmx::ArrayRef<gmx::RVec>       v,
                                            gmx::ArrayRef<const unsigned short> cFREEZE,
                                            t_nrnb*                             nrnb,
                                            const bool                          scaleCoordinates)
{
    // If the support here increases, we need to also add the template declarations
    // for new cases below.
    static_assert(pressureCouplingType == PressureCoupling::Berendsen
                          || pressureCouplingType == PressureCoupling::CRescale,
                  "pressureCouplingScaleBoxAndCoordinates is only implemented for Berendsen and "
                  "C-rescale pressure coupling");

    Matrix3x3 inverseMu;
    if (pressureCouplingType == PressureCoupling::CRescale)
    {
        inverseMu = gmx::invertBoxMatrix(mu);
    }

    /* Scale the positions (for Berendsen and c-rescale) and perhaps the velocities (for c-rescale only) */
    if (scaleCoordinates)
    {
        const int gmx_unused numThreads = gmx_omp_nthreads_get(ModuleMultiThread::Update);
#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int n = start; n < start + nr_atoms; n++)
        {
            // Trivial OpenMP region that does not throw
            int g = 0;
            if (!cFREEZE.empty())
            {
                g = cFREEZE[n];
            }

            if (!nFreeze[g][XX])
            {
                x[n][XX] = mu(XX, XX) * x[n][XX] + mu(YY, XX) * x[n][YY] + mu(ZZ, XX) * x[n][ZZ];
                if (pressureCouplingType == PressureCoupling::CRescale)
                {
                    v[n][XX] = inverseMu(XX, XX) * v[n][XX] + inverseMu(YY, XX) * v[n][YY]
                               + inverseMu(ZZ, XX) * v[n][ZZ];
                }
            }
            if (!nFreeze[g][YY])
            {
                x[n][YY] = mu(YY, YY) * x[n][YY] + mu(ZZ, YY) * x[n][ZZ];
                if (pressureCouplingType == PressureCoupling::CRescale)
                {
                    v[n][YY] = inverseMu(YY, YY) * v[n][YY] + inverseMu(ZZ, YY) * v[n][ZZ];
                }
            }
            if (!nFreeze[g][ZZ])
            {
                x[n][ZZ] = mu(ZZ, ZZ) * x[n][ZZ];
                if (pressureCouplingType == PressureCoupling::CRescale)
                {
                    v[n][ZZ] = inverseMu(ZZ, ZZ) * v[n][ZZ];
                }
            }
        }
    }
    /* compute final boxlengths */
    for (int d = 0; d < DIM; d++)
    {
        box[d][XX] = mu(XX, XX) * box[d][XX] + mu(YY, XX) * box[d][YY] + mu(ZZ, XX) * box[d][ZZ];
        box[d][YY] = mu(YY, YY) * box[d][YY] + mu(ZZ, YY) * box[d][ZZ];
        box[d][ZZ] = mu(ZZ, ZZ) * box[d][ZZ];
    }

    preserveBoxShape(pressureCouplingOptions, deform, box_rel, box);

    /* (un)shifting should NOT be done after this,
     * since the box vectors might have changed
     */
    inc_nrnb(nrnb, eNR_PCOUPL, nr_atoms);
}

template void
pressureCouplingScaleBoxAndCoordinates<PressureCoupling::Berendsen>(const PressureCouplingOptions&,
                                                                    const tensor,
                                                                    const ivec*,
                                                                    const Matrix3x3&,
                                                                    matrix,
                                                                    matrix,
                                                                    int,
                                                                    int,
                                                                    gmx::ArrayRef<gmx::RVec>,
                                                                    gmx::ArrayRef<gmx::RVec>,
                                                                    gmx::ArrayRef<const unsigned short>,
                                                                    t_nrnb*,
                                                                    const bool);


template void
pressureCouplingScaleBoxAndCoordinates<PressureCoupling::CRescale>(const PressureCouplingOptions&,
                                                                   const tensor,
                                                                   const ivec*,
                                                                   const Matrix3x3&,
                                                                   matrix,
                                                                   matrix,
                                                                   int,
                                                                   int,
                                                                   gmx::ArrayRef<gmx::RVec>,
                                                                   gmx::ArrayRef<gmx::RVec>,
                                                                   gmx::ArrayRef<const unsigned short>,
                                                                   t_nrnb*,
                                                                   const bool);

void berendsen_tcoupl(const t_inputrec* ir, gmx_ekindata_t* ekind, real dt, std::vector<double>& therm_integral)
{
    const t_grpopts* opts = &ir->opts;

    for (int i = 0; (i < opts->ngtc); i++)
    {
        real Ek, T;

        if (ir->eI == IntegrationAlgorithm::VV)
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
            real reft               = std::max<real>(0, ekind->currentReferenceTemperature(i));
            real lll                = std::sqrt(1.0 + (dt / opts->tau_t[i]) * (reft / T - 1.0));
            ekind->tcstat[i].lambda = std::max<real>(std::min<real>(lll, 1.25), 0.8);
        }
        else
        {
            ekind->tcstat[i].lambda = 1.0;
        }

        /* Keep track of the amount of energy we are adding to the system */
        therm_integral[i] -= (gmx::square(ekind->tcstat[i].lambda) - 1) * Ek;

        if (debug)
        {
            fprintf(debug, "TC: group %d: T: %g, Lambda: %g\n", i, T, ekind->tcstat[i].lambda);
        }
    }
}

void andersen_tcoupl(const t_inputrec*                   ir,
                     int64_t                             step,
                     const t_commrec*                    cr,
                     const int                           homenr,
                     gmx::ArrayRef<const unsigned short> cTC,
                     gmx::ArrayRef<const real>           invMass,
                     gmx::ArrayRef<gmx::RVec>            v,
                     real                                rate,
                     const std::vector<bool>&            randomize,
                     gmx::ArrayRef<const real>           boltzfac)
{
    const int* gatindex = (haveDDAtomOrdering(*cr) ? cr->dd->globalAtomIndices.data() : nullptr);
    int        i;
    int        gc = 0;
    gmx::ThreeFry2x64<0>               rng(ir->andersen_seed, gmx::RandomDomain::Thermostat);
    gmx::UniformRealDistribution<real> uniformDist;
    gmx::TabulatedNormalDistribution<real, 14> normalDist;

    /* randomize the velocities of the selected particles */

    for (i = 0; i < homenr; i++) /* now loop over the list of atoms */
    {
        int  ng = gatindex ? gatindex[i] : i;
        bool bRandomize;

        rng.restart(step, ng);

        if (!cTC.empty())
        {
            gc = cTC[i]; /* assign the atom to a temperature group if there are more than one */
        }
        if (randomize[gc])
        {
            if (ir->etc == TemperatureCoupling::AndersenMassive)
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

                scal = std::sqrt(boltzfac[gc] * invMass[i]);

                normalDist.reset();

                for (d = 0; d < DIM; d++)
                {
                    v[i][d] = scal * normalDist(rng);
                }
            }
        }
    }
}


static void nosehoover_tcoupl(const gmx_ekindata_t& ekind,
                              real                  dt,
                              gmx::ArrayRef<double> xi,
                              gmx::ArrayRef<double> vxi,
                              const t_extmass&      MassQ)
{
    /* note that this routine does not include Nose-hoover chains yet. Should be easy to add. */

    for (int i = 0; i < ekind.numTemperatureCouplingGroups(); i++)
    {
        const real reft   = std::max<real>(0, ekind.currentReferenceTemperature(i));
        const real oldvxi = vxi[i];
        vxi[i] += dt * MassQ.Qinv[i] * (ekind.tcstat[i].Th - reft);
        xi[i] += dt * (oldvxi + vxi[i]) * 0.5;
    }
}

void trotter_update(const t_inputrec*                   ir,
                    int64_t                             step,
                    gmx_ekindata_t*                     ekind,
                    t_state*                            state,
                    const tensor                        vir,
                    int                                 homenr,
                    gmx::ArrayRef<const unsigned short> cTC,
                    gmx::ArrayRef<const real>           invMass,
                    const t_extmass*                    MassQ,
                    gmx::ArrayRef<std::vector<int>>     trotter_seqlist,
                    TrotterSequence                     trotter_seqno)
{

    int              n, i, d, ngtc, gc = 0, t;
    t_grp_tcstat*    tcstat;
    const t_grpopts* opts;
    int64_t          step_eff;
    real             dt;
    double *         scalefac, dtc;
    rvec             sumv = { 0, 0, 0 };
    bool             bCouple;

    if (trotter_seqno <= TrotterSequence::Two)
    {
        step_eff = step - 1; /* the velocity verlet calls are actually out of order -- the first
                                half step is actually the last half step from the previous step.
                                Thus the first half step actually corresponds to the n-1 step*/
    }
    else
    {
        step_eff = step;
    }

    bCouple = (ir->nsttcouple == 1 || do_per_step(step_eff + ir->nsttcouple, ir->nsttcouple));

    const gmx::ArrayRef<const int> trotter_seq = trotter_seqlist[static_cast<int>(trotter_seqno)];

    if ((trotter_seq[0] == etrtSKIPALL) || (!bCouple))
    {
        return;
    }
    dtc  = ir->nsttcouple * ir->delta_t; /* This is OK for NPT, because nsttcouple == nstpcouple is enforcesd */
    opts = &(ir->opts);                  /* just for ease of referencing */
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
        if ((trotter_seq[i] == etrtBAROV2) || (trotter_seq[i] == etrtBARONHC2)
            || (trotter_seq[i] == etrtNHC2))
        {
            dt = 2 * dtc;
        }
        else
        {
            dt = dtc;
        }

        auto v = makeArrayRef(state->v);
        switch (trotter_seq[i])
        {
            case etrtBAROV:
            case etrtBAROV2:
                boxv_trotter(ir, &(state->veta), dt, state->box, ekind, vir, MassQ);
                break;
            case etrtBARONHC:
            case etrtBARONHC2:
                NHC_trotter(opts,
                            state->nnhpres,
                            ekind,
                            dt,
                            state->nhpres_xi.data(),
                            state->nhpres_vxi.data(),
                            nullptr,
                            &(state->veta),
                            MassQ,
                            FALSE);
                break;
            case etrtNHC:
            case etrtNHC2:
                NHC_trotter(opts,
                            opts->ngtc,
                            ekind,
                            dt,
                            state->nosehoover_xi.data(),
                            state->nosehoover_vxi.data(),
                            scalefac,
                            nullptr,
                            MassQ,
                            (ir->eI == IntegrationAlgorithm::VV));
                /* need to rescale the kinetic energies and velocities here.  Could
                   scale the velocities later, but we need them scaled in order to
                   produce the correct outputs, so we'll scale them here. */

                for (t = 0; t < ngtc; t++)
                {
                    tcstat             = &ekind->tcstat[t];
                    tcstat->vscale_nhc = scalefac[t];
                    tcstat->ekinscaleh_nhc *= (scalefac[t] * scalefac[t]);
                    tcstat->ekinscalef_nhc *= (scalefac[t] * scalefac[t]);
                }
                /* now that we've scaled the groupwise velocities, we can add them up to get the total */
                /* but do we actually need the total? */

                /* modify the velocities as well */
                for (n = 0; n < homenr; n++)
                {
                    if (!cTC.empty()) /* does this conditional need to be here? is this always true?*/
                    {
                        gc = cTC[n];
                    }
                    for (d = 0; d < DIM; d++)
                    {
                        v[n][d] *= scalefac[gc];
                    }

                    if (debug)
                    {
                        for (d = 0; d < DIM; d++)
                        {
                            sumv[d] += (v[n][d]) / invMass[n];
                        }
                    }
                }
                break;
            default: break;
        }
    }
    /* check for conserved momentum -- worth looking at this again eventually, but not working right now.*/
    sfree(scalefac);
}

void init_npt_masses(const t_inputrec& ir, const gmx_ekindata_t& ekind, t_state* state, t_extmass* MassQ, bool bInit)
{
    int              i, j, ngtc, nh;
    const t_grpopts* opts;
    real             reft, kT, ndj, nd;

    opts = &ir.opts; /* just for ease of referencing */
    ngtc = opts->ngtc;
    nh   = state->nhchainlength;

    if (ir.eI == IntegrationAlgorithm::MD)
    {
        if (bInit)
        {
            MassQ->Qinv.resize(ngtc);
        }
        for (i = 0; (i < ngtc); i++)
        {
            if (opts->tau_t[i] > 0 && ekind.currentReferenceTemperature(i) > 0)
            {
                MassQ->Qinv[i] =
                        1.0 / (gmx::square(opts->tau_t[i] / M_2PI) * ekind.currentReferenceTemperature(i));
            }
            else
            {
                MassQ->Qinv[i] = 0.0;
            }
        }
    }
    else if (EI_VV(ir.eI))
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
        MassQ->Winv = (gmx::c_presfac * trace(ir.pressureCouplingOptions.compress) * gmx::c_boltz
                       * ekind.currentEnsembleTemperature())
                      / (DIM * state->vol0 * gmx::square(ir.pressureCouplingOptions.tau_p / M_2PI));
        /* Allocate space for thermostat variables */
        if (bInit)
        {
            MassQ->Qinv.resize(ngtc * nh);
        }

        /* now, set temperature variables */
        for (i = 0; i < ngtc; i++)
        {
            if (opts->tau_t[i] > 0 && ekind.currentReferenceTemperature(i) > 0 && opts->nrdf[i] > 0)
            {
                reft = ekind.currentReferenceTemperature(i);
                nd   = opts->nrdf[i];
                kT   = gmx::c_boltz * reft;
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
                    MassQ->Qinv[i * nh + j] = 1.0 / (gmx::square(opts->tau_t[i] / M_2PI) * ndj * kT);
                }
            }
            else
            {
                for (j = 0; j < nh; j++)
                {
                    MassQ->Qinv[i * nh + j] = 0.0;
                }
            }
        }
    }
}

gmx::EnumerationArray<TrotterSequence, std::vector<int>> init_npt_vars(const t_inputrec*     ir,
                                                                       const gmx_ekindata_t& ekind,
                                                                       t_state*              state,
                                                                       t_extmass*            MassQ,
                                                                       bool bTrotter)
{
    int              i, j, nnhpres, nh;
    const t_grpopts* opts;
    real             bmass, qmass, reft, kT;

    opts    = &(ir->opts); /* just for ease of referencing */
    nnhpres = state->nnhpres;
    nh      = state->nhchainlength;

    if (EI_VV(ir->eI) && (ir->pressureCouplingOptions.epc == PressureCoupling::Mttk)
        && (ir->etc != TemperatureCoupling::NoseHoover))
    {
        gmx_fatal(FARGS, "Cannot do MTTK pressure coupling without Nose-Hoover temperature control");
    }

    /* If we have second order coupling algorithms, initialize their masses here */
    if (ir->pressureCouplingOptions.epc == PressureCoupling::ParrinelloRahman
        || ir->pressureCouplingOptions.epc == PressureCoupling::Mttk
        || ir->etc == TemperatureCoupling::NoseHoover)
    {
        init_npt_masses(*ir, ekind, state, MassQ, TRUE);
    }

    /* first, initialize clear all the trotter calls */
    gmx::EnumerationArray<TrotterSequence, std::vector<int>> trotter_seq;
    for (i = 0; i < static_cast<int>(TrotterSequence::Count); i++)
    {
        trotter_seq[i].resize(NTROTTERPARTS, etrtNONE);
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

    if (ir->eI == IntegrationAlgorithm::VV)
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
    else if (ir->eI == IntegrationAlgorithm::VVAK)
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

    switch (ir->pressureCouplingOptions.epct)
    {
        case PressureCouplingType::Isotropic:
        default: bmass = DIM * DIM; /* recommended mass parameters for isotropic barostat */
    }

    MassQ->QPinv.resize(nnhpres * opts->nhchainlength);

    /* barostat temperature */
    if ((ir->pressureCouplingOptions.tau_p > 0) && (constantEnsembleTemperature(*ir) > 0))
    {
        reft = constantEnsembleTemperature(*ir);
        kT   = gmx::c_boltz * reft;
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
                MassQ->QPinv[i * opts->nhchainlength + j] =
                        1.0 / (gmx::square(opts->tau_t[0] / M_2PI) * qmass * kT);
            }
        }
    }
    else
    {
        for (i = 0; i < nnhpres; i++)
        {
            for (j = 0; j < nh; j++)
            {
                MassQ->QPinv[i * nh + j] = 0.0;
            }
        }
    }
    return trotter_seq;
}

static real energyNoseHoover(const gmx::ArrayRef<const real> degreesOfFreedom,
                             const gmx_ekindata_t&           ekind,
                             const bool                      isTrotterWithConstantTemperature,
                             const t_state*                  state,
                             const t_extmass*                MassQ)
{
    real energy = 0;

    int nh = state->nhchainlength;

    for (int i = 0; i < state->ngtc; i++)
    {
        const double* ixi   = &state->nosehoover_xi[i * nh];
        const double* ivxi  = &state->nosehoover_vxi[i * nh];
        const double* iQinv = &(MassQ->Qinv[i * nh]);

        real nd   = degreesOfFreedom[i];
        real reft = std::max<real>(ekind.currentReferenceTemperature(i), 0);
        real kT   = gmx::c_boltz * reft;

        if (nd > 0.0)
        {
            if (isTrotterWithConstantTemperature)
            {
                /* contribution from the thermal momenta of the NH chain */
                for (int j = 0; j < nh; j++)
                {
                    if (iQinv[j] > 0)
                    {
                        energy += 0.5 * gmx::square(ivxi[j]) / iQinv[j];
                        /* contribution from the thermal variable of the NH chain */
                        real ndj = 0;
                        if (j == 0)
                        {
                            ndj = nd;
                        }
                        else
                        {
                            ndj = 1;
                        }
                        energy += ndj * ixi[j] * kT;
                    }
                }
            }
            else /* Other non Trotter temperature NH control  -- no chains yet. */
            {
                energy += 0.5 * gmx::c_boltz * nd * gmx::square(ivxi[0]) / iQinv[0];
                energy += nd * ixi[0] * kT;
            }
        }
    }

    return energy;
}

/* Returns the energy from the barostat thermostat chain */
static real energyPressureMTTK(const real ensembleTemperature, const t_state* state, const t_extmass* MassQ)
{
    real energy = 0;

    int nh = state->nhchainlength;

    for (int i = 0; i < state->nnhpres; i++)
    {
        /* note -- assumes only one degree of freedom that is thermostatted in barostat */
        real reft = std::max<real>(ensembleTemperature, 0.0); /* using 'System' temperature */
        real kT   = gmx::c_boltz * reft;

        for (int j = 0; j < nh; j++)
        {
            double iQinv = MassQ->QPinv[i * nh + j];
            if (iQinv > 0)
            {
                energy += 0.5 * gmx::square(state->nhpres_vxi[i * nh + j]) / iQinv;
                /* contribution from the thermal variable of the NH chain */
                energy += state->nhpres_xi[i * nh + j] * kT;
            }
            if (debug)
            {
                fprintf(debug,
                        "P-T-group: %10d Chain %4d ThermV: %15.8f ThermX: %15.8f",
                        i,
                        j,
                        state->nhpres_vxi[i * nh + j],
                        state->nhpres_xi[i * nh + j]);
            }
        }
    }

    return energy;
}

/* Returns the energy accumulated by the V-rescale or Berendsen thermostat */
static real energyVrescale(const t_state* state)
{
    return std::accumulate(state->therm_integral.begin(), state->therm_integral.end(), 0.);
}

real NPT_energy(const PressureCouplingOptions&  pressureCouplingOptions,
                const TemperatureCoupling       etc,
                const gmx::ArrayRef<const real> degreesOfFreedom,
                const gmx_ekindata_t&           ekind,
                const bool                      isTrotterWithConstantTemperature,
                const t_state*                  state,
                const t_extmass*                MassQ)
{
    real energyNPT = 0;

    if (pressureCouplingOptions.epc != PressureCoupling::No)
    {
        /* Compute the contribution of the pressure to the conserved quantity*/

        real vol = det(state->box);

        switch (pressureCouplingOptions.epc)
        {
            case PressureCoupling::ParrinelloRahman:
            {
                /* contribution from the pressure momenta */
                tensor invMass;
                calcParrinelloRahmanInvMass(pressureCouplingOptions, state->box, invMass);
                for (int d = 0; d < DIM; d++)
                {
                    for (int n = 0; n <= d; n++)
                    {
                        if (invMass[d][n] > 0)
                        {
                            energyNPT += 0.5 * gmx::square(state->boxv[d][n])
                                         / (invMass[d][n] * gmx::c_presfac);
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
                energyNPT += vol * trace(pressureCouplingOptions.ref_p) / (DIM * gmx::c_presfac);
                break;
            }
            case PressureCoupling::Mttk:
                /* contribution from the pressure momenta */
                energyNPT += 0.5 * gmx::square(state->veta) / MassQ->Winv;

                /* contribution from the PV term */
                energyNPT += vol * trace(pressureCouplingOptions.ref_p) / (DIM * gmx::c_presfac);

                if (pressureCouplingOptions.epc == PressureCoupling::Mttk)
                {
                    /* contribution from the MTTK chain */
                    energyNPT += energyPressureMTTK(ekind.currentEnsembleTemperature(), state, MassQ);
                }
                break;
            case PressureCoupling::Berendsen:
            case PressureCoupling::CRescale: energyNPT += state->baros_integral; break;
            default:
                GMX_RELEASE_ASSERT(
                        false,
                        "Conserved energy quantity for pressure coupling is not handled. A case "
                        "should be added with either the conserved quantity added or nothing added "
                        "and an exclusion added to integratorHasConservedEnergyQuantity().");
        }
    }

    switch (etc)
    {
        case TemperatureCoupling::No: break;
        case TemperatureCoupling::VRescale:
        case TemperatureCoupling::Berendsen: energyNPT += energyVrescale(state); break;
        case TemperatureCoupling::NoseHoover:
            energyNPT += energyNoseHoover(
                    degreesOfFreedom, ekind, isTrotterWithConstantTemperature, state, MassQ);
            break;
        case TemperatureCoupling::Andersen:
        case TemperatureCoupling::AndersenMassive:
            // Not supported, excluded in integratorHasConservedEnergyQuantity()
            break;
        default:
            GMX_RELEASE_ASSERT(
                    false,
                    "Conserved energy quantity for temperature coupling is not handled. A case "
                    "should be added with either the conserved quantity added or nothing added and "
                    "an exclusion added to integratorHasConservedEnergyQuantity().");
    }

    return energyNPT;
}


static real vrescale_sumnoises(real nn, gmx::ThreeFry2x64<>* rng, gmx::NormalDistribution<real>* normalDist)
{
    /*
     * Returns the sum of nn independent gaussian noises squared
     * (i.e. equivalent to summing the square of the return values
     * of nn calls to a normal distribution).
     */
    const real                   ndeg_tol = 0.0001;
    real                         r;
    gmx::GammaDistribution<real> gammaDist(0.5 * nn, 1.0);

    if (nn < 2 + ndeg_tol)
    {
        int  nn_int, i;
        real gauss;

        nn_int = gmx::roundToInt(nn);

        if (nn - nn_int < -ndeg_tol || nn - nn_int > ndeg_tol)
        {
            gmx_fatal(FARGS,
                      "The v-rescale thermostat was called with a group with #DOF=%f, but for "
                      "#DOF<3 only integer #DOF are supported",
                      nn + 1);
        }

        r = 0;
        for (i = 0; i < nn_int; i++)
        {
            gauss = (*normalDist)(*rng);
            r += gauss * gauss;
        }
    }
    else
    {
        /* Use a gamma distribution for any real nn > 2 */
        r = 2.0 * gammaDist(*rng);
    }

    return r;
}

real vrescale_resamplekin(real kk, real sigma, real ndeg, real taut, int64_t step, int64_t seed)
{
    /*
     * Generates a new value for the kinetic energy,
     * according to Bussi et al JCP (2007), Eq. (A7)
     * kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
     * sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
     * ndeg:  number of degrees of freedom of the atoms to be thermalized
     * taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
     */
    real                          factor, rr, ekin_new;
    gmx::ThreeFry2x64<64>         rng(seed, gmx::RandomDomain::Thermostat);
    gmx::NormalDistribution<real> normalDist;

    if (taut > 0.1)
    {
        factor = std::exp(-1.0 / taut);
    }
    else
    {
        factor = 0.0;
    }

    rng.restart(step, 0);

    rr = normalDist(rng);

    ekin_new = kk
               + (1.0 - factor)
                         * (sigma * (vrescale_sumnoises(ndeg - 1, &rng, &normalDist) + rr * rr) / ndeg - kk)
               + 2.0 * rr * std::sqrt(kk * sigma / ndeg * (1.0 - factor) * factor);

    return ekin_new;
}

void vrescale_tcoupl(const t_inputrec* ir, int64_t step, gmx_ekindata_t* ekind, real dt, gmx::ArrayRef<double> therm_integral)
{
    const t_grpopts* opts;
    int              i;
    real             Ek, Ek_ref1, Ek_ref, Ek_new;

    opts = &ir->opts;

    for (i = 0; (i < opts->ngtc); i++)
    {
        if (ir->eI == IntegrationAlgorithm::VV)
        {
            Ek = trace(ekind->tcstat[i].ekinf);
        }
        else
        {
            Ek = trace(ekind->tcstat[i].ekinh);
        }

        if (opts->tau_t[i] >= 0 && opts->nrdf[i] > 0 && Ek > 0)
        {
            Ek_ref1 = 0.5 * ekind->currentReferenceTemperature(i) * gmx::c_boltz;
            Ek_ref  = Ek_ref1 * opts->nrdf[i];

            Ek_new = vrescale_resamplekin(Ek, Ek_ref, opts->nrdf[i], opts->tau_t[i] / dt, step, ir->ld_seed);

            /* Analytically Ek_new>=0, but we check for rounding errors */
            if (Ek_new <= 0)
            {
                ekind->tcstat[i].lambda = 0.0;
            }
            else
            {
                ekind->tcstat[i].lambda = std::sqrt(Ek_new / Ek);
            }

            therm_integral[i] -= Ek_new - Ek;

            if (debug)
            {
                fprintf(debug,
                        "TC: group %d: Ekr %g, Ek %g, Ek_new %g, Lambda: %g\n",
                        i,
                        Ek_ref,
                        Ek,
                        Ek_new,
                        ekind->tcstat[i].lambda);
            }
        }
        else
        {
            ekind->tcstat[i].lambda = 1.0;
        }
    }
}

void rescale_velocities(const gmx_ekindata_t*               ekind,
                        gmx::ArrayRef<const unsigned short> cTC,
                        int                                 start,
                        int                                 end,
                        gmx::ArrayRef<gmx::RVec>            v)
{
    gmx::ArrayRef<const t_grp_tcstat> tcstat = ekind->tcstat;

    for (int n = start; n < end; n++)
    {
        int gt = 0;
        if (!cTC.empty())
        {
            gt = cTC[n];
        }
        const real lg = tcstat[gt].lambda;
        for (int d = 0; d < DIM; d++)
        {
            v[n][d] *= lg;
        }
    }
}

//! Initialize simulated annealing
bool initSimulatedAnnealing(const t_inputrec& ir, gmx_ekindata_t* ekind, gmx::Update* upd)
{
    bool doSimAnnealing = doSimulatedAnnealing(ir);
    if (doSimAnnealing)
    {
        update_annealing_target_temp(ir, ir.init_t, ekind, upd);
    }
    return doSimAnnealing;
}

/*!
 * \brief Compute the new annealing temperature for a temperature group
 *
 * \param inputrec          The input record
 * \param temperatureGroup  The temperature group
 * \param time              The current time
 * \return  The new reference temperature for the group
 */
static real computeAnnealingTargetTemperature(const t_inputrec& inputrec, int temperatureGroup, real time)
{
    GMX_RELEASE_ASSERT(temperatureGroup >= 0 && temperatureGroup < inputrec.opts.ngtc,
                       "Invalid temperature group.");
    GMX_RELEASE_ASSERT(inputrec.opts.annealing[temperatureGroup] != SimulatedAnnealing::No,
                       "Should only compute temperature of annealed groups");
    GMX_RELEASE_ASSERT(
            inputrec.opts.annealing[temperatureGroup] == SimulatedAnnealing::Single
                    || inputrec.opts.annealing[temperatureGroup] == SimulatedAnnealing::Periodic,
            gmx::formatString("Unknown simulated annealing algorithm for temperature group %d", temperatureGroup)
                    .c_str());
    real       thist   = 0;
    const auto npoints = inputrec.opts.anneal_npoints[temperatureGroup];
    if (inputrec.opts.annealing[temperatureGroup] == SimulatedAnnealing::Periodic)
    {
        /* calculate time modulo the period */
        const auto pert = inputrec.opts.anneal_time[temperatureGroup][npoints - 1];
        const auto n    = static_cast<int>(time / pert);
        thist           = time - n * pert; /* modulo time */
        /* Make sure rounding didn't get us outside the interval */
        if (std::fabs(thist - pert) < GMX_REAL_EPS * 100)
        {
            thist = 0;
        }
    }
    else if (inputrec.opts.annealing[temperatureGroup] == SimulatedAnnealing::Single)
    {
        thist = time;
    }
    /* We are doing annealing for this group if we got here,
     * and we have the (relative) time as thist.
     * calculate target temp */
    int j = 0;
    while ((j < npoints - 1) && (thist > (inputrec.opts.anneal_time[temperatureGroup][j + 1])))
    {
        j++;
    }
    if (j < npoints - 1)
    {
        /* Found our position between points j and j+1.
         * Interpolate: x is the amount from j+1, (1-x) from point j
         * First treat possible jumps in temperature as a special case.
         */
        if ((inputrec.opts.anneal_time[temperatureGroup][j + 1]
             - inputrec.opts.anneal_time[temperatureGroup][j])
            < GMX_REAL_EPS * 100)
        {
            return inputrec.opts.anneal_temp[temperatureGroup][j + 1];
        }
        else
        {
            const real x = ((thist - inputrec.opts.anneal_time[temperatureGroup][j])
                            / (inputrec.opts.anneal_time[temperatureGroup][j + 1]
                               - inputrec.opts.anneal_time[temperatureGroup][j]));
            return x * inputrec.opts.anneal_temp[temperatureGroup][j + 1]
                   + (1 - x) * inputrec.opts.anneal_temp[temperatureGroup][j];
        }
    }
    else
    {
        return inputrec.opts.anneal_temp[temperatureGroup][npoints - 1];
    }
}

/* set target temperatures if we are annealing */
void update_annealing_target_temp(const t_inputrec& ir, const real t, gmx_ekindata_t* ekind, gmx::Update* upd)
{
    for (int temperatureGroup = 0; temperatureGroup < ir.opts.ngtc; temperatureGroup++)
    {
        if (ir.opts.annealing[temperatureGroup] != SimulatedAnnealing::No)
        {
            ekind->setCurrentReferenceTemperature(
                    temperatureGroup, computeAnnealingTargetTemperature(ir, temperatureGroup, t));
        }
    }

    upd->update_temperature_constants(ir, *ekind);
}

void pleaseCiteCouplingAlgorithms(FILE* fplog, const t_inputrec& ir)
{
    if (EI_DYNAMICS(ir.eI))
    {
        if (ir.etc == TemperatureCoupling::Berendsen)
        {
            please_cite(fplog, "Berendsen84a");
        }
        if (ir.etc == TemperatureCoupling::VRescale)
        {
            please_cite(fplog, "Bussi2007a");
        }
        if (ir.pressureCouplingOptions.epc == PressureCoupling::CRescale)
        {
            please_cite(fplog, "Bernetti2020");
        }
        // TODO this is actually an integrator, not a coupling algorithm
        if (ir.eI == IntegrationAlgorithm::SD1)
        {
            please_cite(fplog, "Goga2012");
        }
    }
}
