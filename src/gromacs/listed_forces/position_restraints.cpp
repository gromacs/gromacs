/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
 * \brief This file defines low-level functions necessary for
 * computing energies and forces for position restraints.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_listed_forces
 */

#include "gmxpre.h"

#include "position_restraints.h"

#include <cassert>
#include <cmath>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"

struct gmx_wallcycle;

namespace
{

/*! \brief returns dx, rdist, and dpdl for functions posres() and fbposres()
 */
void posres_dx(const rvec      x,
               const rvec      pos0A,
               const rvec      pos0B,
               const rvec      comA_sc,
               const rvec      comB_sc,
               real            lambda,
               const t_pbc*    pbc,
               RefCoordScaling refcoord_scaling,
               int             npbcdim,
               rvec            dx,
               rvec            rdist,
               rvec            dpdl)
{
    int  m, d;
    real posA, posB, L1, ref = 0.;
    rvec pos;

    L1 = 1.0 - lambda;

    for (m = 0; m < DIM; m++)
    {
        posA = pos0A[m];
        posB = pos0B[m];
        if (m < npbcdim)
        {
            switch (refcoord_scaling)
            {
                case RefCoordScaling::No:
                    ref      = 0;
                    rdist[m] = L1 * posA + lambda * posB;
                    dpdl[m]  = posB - posA;
                    break;
                case RefCoordScaling::All:
                    /* Box relative coordinates are stored for dimensions with pbc */
                    posA *= pbc->box[m][m];
                    posB *= pbc->box[m][m];
                    assert(npbcdim <= DIM);
                    for (d = m + 1; d < npbcdim && d < DIM; d++)
                    {
                        posA += pos0A[d] * pbc->box[d][m];
                        posB += pos0B[d] * pbc->box[d][m];
                    }
                    ref      = L1 * posA + lambda * posB;
                    rdist[m] = 0;
                    dpdl[m]  = posB - posA;
                    break;
                case RefCoordScaling::Com:
                    ref      = L1 * comA_sc[m] + lambda * comB_sc[m];
                    rdist[m] = L1 * posA + lambda * posB;
                    dpdl[m]  = comB_sc[m] - comA_sc[m] + posB - posA;
                    break;
                default: gmx_fatal(FARGS, "No such scaling method implemented");
            }
        }
        else
        {
            ref      = L1 * posA + lambda * posB;
            rdist[m] = 0;
            dpdl[m]  = posB - posA;
        }

        /* We do pbc_dx with ref+rdist,
         * since with only ref we can be up to half a box vector wrong.
         */
        pos[m] = ref + rdist[m];
    }

    if (pbc)
    {
        pbc_dx(pbc, x, pos, dx);
    }
    else
    {
        rvec_sub(x, pos, dx);
    }
}

/*! \brief Computes forces and potential for flat-bottom cylindrical restraints.
 *         Returns the flat-bottom potential. */
real do_fbposres_cylinder(int fbdim, rvec fm, rvec dx, real rfb, real kk, gmx_bool bInvert)
{
    int  d;
    real dr, dr2, invdr, v, rfb2;

    dr2  = 0.0;
    rfb2 = gmx::square(rfb);
    v    = 0.0;

    for (d = 0; d < DIM; d++)
    {
        if (d != fbdim)
        {
            dr2 += gmx::square(dx[d]);
        }
    }

    if (dr2 > 0.0 && ((dr2 > rfb2 && !bInvert) || (dr2 < rfb2 && bInvert)))
    {
        dr    = std::sqrt(dr2);
        invdr = 1. / dr;
        v     = 0.5 * kk * gmx::square(dr - rfb);
        for (d = 0; d < DIM; d++)
        {
            if (d != fbdim)
            {
                fm[d] = -kk * (dr - rfb) * dx[d] * invdr; /* Force pointing to the center */
            }
        }
    }

    return v;
}

/*! \brief Compute energies and forces for flat-bottomed position restraints
 *
 * Returns the flat-bottomed potential. Same PBC treatment as in
 * normal position restraints */
real fbposres(int                   nbonds,
              const t_iatom         forceatoms[],
              const t_iparams       forceparams[],
              const rvec            x[],
              gmx::ForceWithVirial* forceWithVirial,
              const t_pbc*          pbc,
              RefCoordScaling       refcoord_scaling,
              PbcType               pbcType,
              const rvec            com)
/* compute flat-bottomed positions restraints */
{
    int              i, ai, m, d, type, npbcdim = 0, fbdim;
    const t_iparams* pr;
    real             kk, v;
    real             dr, dr2, rfb, rfb2, fact;
    rvec             com_sc, rdist, dx, dpdl, fm;
    gmx_bool         bInvert;

    npbcdim = numPbcDimensions(pbcType);
    GMX_ASSERT((pbcType == PbcType::No) == (npbcdim == 0), "");
    if (refcoord_scaling == RefCoordScaling::Com)
    {
        clear_rvec(com_sc);
        for (m = 0; m < npbcdim; m++)
        {
            assert(npbcdim <= DIM);
            for (d = m; d < npbcdim; d++)
            {
                com_sc[m] += com[d] * pbc->box[d][m];
            }
        }
    }

    rvec* f      = as_rvec_array(forceWithVirial->force_.data());
    real  vtot   = 0.0;
    rvec  virial = { 0 };
    for (i = 0; (i < nbonds);)
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        pr   = &forceparams[type];

        /* same calculation as for normal posres, but with identical A and B states, and lambda==0 */
        posres_dx(x[ai],
                  forceparams[type].fbposres.pos0,
                  forceparams[type].fbposres.pos0,
                  com_sc,
                  com_sc,
                  0.0,
                  pbc,
                  refcoord_scaling,
                  npbcdim,
                  dx,
                  rdist,
                  dpdl);

        clear_rvec(fm);
        v = 0.0;

        kk   = pr->fbposres.k;
        rfb  = pr->fbposres.r;
        rfb2 = gmx::square(rfb);

        /* with rfb<0, push particle out of the sphere/cylinder/layer */
        bInvert = FALSE;
        if (rfb < 0.)
        {
            bInvert = TRUE;
            rfb     = -rfb;
        }

        switch (pr->fbposres.geom)
        {
            case efbposresSPHERE:
                /* spherical flat-bottom posres */
                dr2 = norm2(dx);
                if (dr2 > 0.0 && ((dr2 > rfb2 && !bInvert) || (dr2 < rfb2 && bInvert)))
                {
                    dr   = std::sqrt(dr2);
                    v    = 0.5 * kk * gmx::square(dr - rfb);
                    fact = -kk * (dr - rfb) / dr; /* Force pointing to the center pos0 */
                    svmul(fact, dx, fm);
                }
                break;
            case efbposresCYLINDERX:
                /* cylindrical flat-bottom posres in y-z plane. fm[XX] = 0. */
                fbdim = XX;
                v     = do_fbposres_cylinder(fbdim, fm, dx, rfb, kk, bInvert);
                break;
            case efbposresCYLINDERY:
                /* cylindrical flat-bottom posres in x-z plane. fm[YY] = 0. */
                fbdim = YY;
                v     = do_fbposres_cylinder(fbdim, fm, dx, rfb, kk, bInvert);
                break;
            case efbposresCYLINDER:
            /* equivalent to efbposresCYLINDERZ for backwards compatibility */
            case efbposresCYLINDERZ:
                /* cylindrical flat-bottom posres in x-y plane. fm[ZZ] = 0. */
                fbdim = ZZ;
                v     = do_fbposres_cylinder(fbdim, fm, dx, rfb, kk, bInvert);
                break;
            case efbposresX: /* fbdim=XX */
            case efbposresY: /* fbdim=YY */
            case efbposresZ: /* fbdim=ZZ */
                /* 1D flat-bottom potential */
                fbdim = pr->fbposres.geom - efbposresX;
                dr    = dx[fbdim];
                if ((dr > rfb && !bInvert) || (0 < dr && dr < rfb && bInvert))
                {
                    v         = 0.5 * kk * gmx::square(dr - rfb);
                    fm[fbdim] = -kk * (dr - rfb);
                }
                else if ((dr < (-rfb) && !bInvert) || ((-rfb) < dr && dr < 0 && bInvert))
                {
                    v         = 0.5 * kk * gmx::square(dr + rfb);
                    fm[fbdim] = -kk * (dr + rfb);
                }
                break;
        }

        vtot += v;

        for (m = 0; (m < DIM); m++)
        {
            f[ai][m] += fm[m];
            /* Here we correct for the pbc_dx which included rdist */
            virial[m] -= 0.5 * (dx[m] + rdist[m]) * fm[m];
        }
    }

    forceWithVirial->addVirialContribution(virial);

    return vtot;
}


/*! \brief Compute energies and forces, when requested, for position restraints
 *
 * Note that position restraints require a different pbc treatment
 * from other bondeds */
template<bool computeForce>
real posres(int                   nbonds,
            const t_iatom         forceatoms[],
            const t_iparams       forceparams[],
            const rvec            x[],
            gmx::ForceWithVirial* forceWithVirial,
            const struct t_pbc*   pbc,
            real                  lambda,
            real*                 dvdlambda,
            RefCoordScaling       refcoord_scaling,
            PbcType               pbcType,
            const rvec            comA,
            const rvec            comB)
{
    int              i, ai, m, d, type, npbcdim = 0;
    const t_iparams* pr;
    real             kk, fm;
    rvec             comA_sc, comB_sc, rdist, dpdl, dx;

    npbcdim = numPbcDimensions(pbcType);
    GMX_ASSERT((pbcType == PbcType::No) == (npbcdim == 0), "");
    if (refcoord_scaling == RefCoordScaling::Com)
    {
        clear_rvec(comA_sc);
        clear_rvec(comB_sc);
        for (m = 0; m < npbcdim; m++)
        {
            assert(npbcdim <= DIM);
            for (d = m; d < npbcdim; d++)
            {
                comA_sc[m] += comA[d] * pbc->box[d][m];
                comB_sc[m] += comB[d] * pbc->box[d][m];
            }
        }
    }

    const real L1 = 1.0 - lambda;

    rvec* f;
    if (computeForce)
    {
        GMX_ASSERT(forceWithVirial != nullptr, "When forces are requested we need a force object");
        f = as_rvec_array(forceWithVirial->force_.data());
    }
    real vtot = 0.0;
    /* Use intermediate virial buffer to reduce reduction rounding errors */
    rvec virial = { 0 };
    for (i = 0; (i < nbonds);)
    {
        type = forceatoms[i++];
        ai   = forceatoms[i++];
        pr   = &forceparams[type];

        /* return dx, rdist, and dpdl */
        posres_dx(x[ai],
                  forceparams[type].posres.pos0A,
                  forceparams[type].posres.pos0B,
                  comA_sc,
                  comB_sc,
                  lambda,
                  pbc,
                  refcoord_scaling,
                  npbcdim,
                  dx,
                  rdist,
                  dpdl);

        for (m = 0; (m < DIM); m++)
        {
            kk = L1 * pr->posres.fcA[m] + lambda * pr->posres.fcB[m];
            fm = -kk * dx[m];
            vtot += 0.5 * kk * dx[m] * dx[m];
            *dvdlambda += 0.5 * (pr->posres.fcB[m] - pr->posres.fcA[m]) * dx[m] * dx[m] + fm * dpdl[m];

            /* Here we correct for the pbc_dx which included rdist */
            if (computeForce)
            {
                f[ai][m] += fm;
                virial[m] -= 0.5 * (dx[m] + rdist[m]) * fm;
            }
        }
    }

    if (computeForce)
    {
        forceWithVirial->addVirialContribution(virial);
    }

    return vtot;
}

} // namespace

void posres_wrapper(t_nrnb*                       nrnb,
                    const InteractionDefinitions& idef,
                    const struct t_pbc*           pbc,
                    const rvec*                   x,
                    gmx_enerdata_t*               enerd,
                    gmx::ArrayRef<const real>     lambda,
                    const t_forcerec*             fr,
                    gmx::ForceWithVirial*         forceWithVirial)
{
    real v, dvdl;

    dvdl = 0;
    v    = posres<true>(idef.il[F_POSRES].size(),
                     idef.il[F_POSRES].iatoms.data(),
                     idef.iparams_posres.data(),
                     x,
                     forceWithVirial,
                     fr->pbcType == PbcType::No ? nullptr : pbc,
                     lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Restraint)],
                     &dvdl,
                     fr->rc_scaling,
                     fr->pbcType,
                     fr->posres_com,
                     fr->posres_comB);
    enerd->term[F_POSRES] += v;
    /* If just the force constant changes, the FEP term is linear,
     * but if k changes, it is not.
     */
    enerd->dvdl_nonlin[FreeEnergyPerturbationCouplingType::Restraint] += dvdl;
    inc_nrnb(nrnb, eNR_POSRES, gmx::exactDiv(idef.il[F_POSRES].size(), 2));
}

void posres_wrapper_lambda(struct gmx_wallcycle*         wcycle,
                           const InteractionDefinitions& idef,
                           const struct t_pbc*           pbc,
                           const rvec                    x[],
                           gmx_enerdata_t*               enerd,
                           gmx::ArrayRef<const real>     lambda,
                           const t_forcerec*             fr)
{
    wallcycle_sub_start_nocount(wcycle, WallCycleSubCounter::Restraints);

    auto& foreignTerms = enerd->foreignLambdaTerms;
    for (int i = 0; i < 1 + foreignTerms.numLambdas(); i++)
    {
        real dvdl = 0;

        const real lambda_dum =
                (i == 0 ? lambda[static_cast<int>(FreeEnergyPerturbationCouplingType::Restraint)]
                        : enerd->foreignLambdaTerms.foreignLambdas(
                                FreeEnergyPerturbationCouplingType::Restraint)[i - 1]);
        const real v = posres<false>(idef.il[F_POSRES].size(),
                                     idef.il[F_POSRES].iatoms.data(),
                                     idef.iparams_posres.data(),
                                     x,
                                     nullptr,
                                     fr->pbcType == PbcType::No ? nullptr : pbc,
                                     lambda_dum,
                                     &dvdl,
                                     fr->rc_scaling,
                                     fr->pbcType,
                                     fr->posres_com,
                                     fr->posres_comB);
        foreignTerms.accumulate(i, FreeEnergyPerturbationCouplingType::Restraint, v, dvdl);
    }
    wallcycle_sub_stop(wcycle, WallCycleSubCounter::Restraints);
}

/*! \brief Helper function that wraps calls to fbposres for
    free-energy perturbation */
void fbposres_wrapper(t_nrnb*                       nrnb,
                      const InteractionDefinitions& idef,
                      const struct t_pbc*           pbc,
                      const rvec*                   x,
                      gmx_enerdata_t*               enerd,
                      const t_forcerec*             fr,
                      gmx::ForceWithVirial*         forceWithVirial)
{
    real v;

    v = fbposres(idef.il[F_FBPOSRES].size(),
                 idef.il[F_FBPOSRES].iatoms.data(),
                 idef.iparams_fbposres.data(),
                 x,
                 forceWithVirial,
                 fr->pbcType == PbcType::No ? nullptr : pbc,
                 fr->rc_scaling,
                 fr->pbcType,
                 fr->posres_com);
    enerd->term[F_FBPOSRES] += v;
    inc_nrnb(nrnb, eNR_FBPOSRES, gmx::exactDiv(idef.il[F_FBPOSRES].size(), 2));
}
