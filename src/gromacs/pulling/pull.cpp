/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
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
#include "gmxpre.h"

#include "pull.h"

#include "config.h"

#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <memory>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/mutex.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "pull_internal.h"

namespace gmx
{
extern template LocalAtomSet LocalAtomSetManager::add<void, void>(ArrayRef<const int> globalAtomIndex);
} // namespace gmx

static int groupPbcFromParams(const t_pull_group& params, bool setPbcRefToPrevStepCOM)
{
    if (params.nat <= 1)
    {
        /* no PBC required */
        return epgrppbcNONE;
    }
    else if (params.pbcatom >= 0)
    {
        if (setPbcRefToPrevStepCOM)
        {
            return epgrppbcPREVSTEPCOM;
        }
        else
        {
            return epgrppbcREFAT;
        }
    }
    else
    {
        return epgrppbcCOS;
    }
}

/* NOTE: The params initialization currently copies pointers.
 *       So the lifetime of the source, currently always inputrec,
 *       should not end before that of this object.
 *       This will be fixed when the pointers are replacted by std::vector.
 */
pull_group_work_t::pull_group_work_t(const t_pull_group& params,
                                     gmx::LocalAtomSet   atomSet,
                                     bool                bSetPbcRefToPrevStepCOM) :
    params(params),
    epgrppbc(groupPbcFromParams(params, bSetPbcRefToPrevStepCOM)),
    needToCalcCom(false),
    atomSet(atomSet),
    mwscale(0),
    wscale(0),
    invtm(0)
{
    clear_dvec(x);
    clear_dvec(xp);
};

bool pull_coordinate_is_angletype(const t_pull_coord* pcrd)
{
    return (pcrd->eGeom == epullgANGLE || pcrd->eGeom == epullgDIHEDRAL || pcrd->eGeom == epullgANGLEAXIS);
}

static bool pull_coordinate_is_directional(const t_pull_coord* pcrd)
{
    return (pcrd->eGeom == epullgDIR || pcrd->eGeom == epullgDIRPBC
            || pcrd->eGeom == epullgDIRRELATIVE || pcrd->eGeom == epullgCYL);
}

const char* pull_coordinate_units(const t_pull_coord* pcrd)
{
    return pull_coordinate_is_angletype(pcrd) ? "deg" : "nm";
}

double pull_conversion_factor_userinput2internal(const t_pull_coord* pcrd)
{
    if (pull_coordinate_is_angletype(pcrd))
    {
        return DEG2RAD;
    }
    else
    {
        return 1.0;
    }
}

double pull_conversion_factor_internal2userinput(const t_pull_coord* pcrd)
{
    if (pull_coordinate_is_angletype(pcrd))
    {
        return RAD2DEG;
    }
    else
    {
        return 1.0;
    }
}

/* Apply forces in a mass weighted fashion for part of the pull group */
static void apply_forces_grp_part(const pull_group_work_t* pgrp,
                                  int                      ind_start,
                                  int                      ind_end,
                                  const t_mdatoms*         md,
                                  const dvec               f_pull,
                                  int                      sign,
                                  rvec*                    f)
{
    double inv_wm = pgrp->mwscale;

    auto localAtomIndices = pgrp->atomSet.localIndex();
    for (int i = ind_start; i < ind_end; i++)
    {
        int    ii    = localAtomIndices[i];
        double wmass = md->massT[ii];
        if (!pgrp->localWeights.empty())
        {
            wmass *= pgrp->localWeights[i];
        }

        for (int d = 0; d < DIM; d++)
        {
            f[ii][d] += sign * wmass * f_pull[d] * inv_wm;
        }
    }
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_grp(const pull_group_work_t* pgrp,
                             const t_mdatoms*         md,
                             const dvec               f_pull,
                             int                      sign,
                             rvec*                    f,
                             int                      nthreads)
{
    auto localAtomIndices = pgrp->atomSet.localIndex();

    if (pgrp->params.nat == 1 && pgrp->atomSet.numAtomsLocal() == 1)
    {
        /* Only one atom and our rank has this atom: we can skip
         * the mass weighting, which means that this code also works
         * for mass=0, e.g. with a virtual site.
         */
        for (int d = 0; d < DIM; d++)
        {
            f[localAtomIndices[0]][d] += sign * f_pull[d];
        }
    }
    else
    {
        if (localAtomIndices.size() <= c_pullMaxNumLocalAtomsSingleThreaded)
        {
            apply_forces_grp_part(pgrp, 0, localAtomIndices.size(), md, f_pull, sign, f);
        }
        else
        {
#pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int th = 0; th < nthreads; th++)
            {
                int ind_start = (localAtomIndices.size() * (th + 0)) / nthreads;
                int ind_end   = (localAtomIndices.size() * (th + 1)) / nthreads;
                apply_forces_grp_part(pgrp, ind_start, ind_end, md, f_pull, sign, f);
            }
        }
    }
}

/* Apply forces in a mass weighted fashion to a cylinder group */
static void apply_forces_cyl_grp(const pull_group_work_t* pgrp,
                                 const double             dv_corr,
                                 const t_mdatoms*         md,
                                 const dvec               f_pull,
                                 double                   f_scal,
                                 int                      sign,
                                 rvec*                    f,
                                 int gmx_unused nthreads)
{
    double inv_wm = pgrp->mwscale;

    auto localAtomIndices = pgrp->atomSet.localIndex();

    /* The cylinder group is always a slab in the system, thus large.
     * Therefore we always thread-parallelize this group.
     */
    int numAtomsLocal = localAtomIndices.size();
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i = 0; i < numAtomsLocal; i++)
    {
        double weight = pgrp->localWeights[i];
        if (weight == 0)
        {
            continue;
        }
        int    ii   = localAtomIndices[i];
        double mass = md->massT[ii];
        /* The stored axial distance from the cylinder center (dv) needs
         * to be corrected for an offset (dv_corr), which was unknown when
         * we calculated dv.
         */
        double dv_com = pgrp->dv[i] + dv_corr;

        /* Here we not only add the pull force working along vec (f_pull),
         * but also a radial component, due to the dependence of the weights
         * on the radial distance.
         */
        for (int m = 0; m < DIM; m++)
        {
            f[ii][m] += sign * inv_wm * (mass * weight * f_pull[m] + pgrp->mdw[i][m] * dv_com * f_scal);
        }
    }
}

/* Apply torque forces in a mass weighted fashion to the groups that define
 * the pull vector direction for pull coordinate pcrd.
 */
static void apply_forces_vec_torque(const struct pull_t*     pull,
                                    const pull_coord_work_t* pcrd,
                                    const t_mdatoms*         md,
                                    rvec*                    f)
{
    const PullCoordSpatialData& spatialData = pcrd->spatialData;

    /* The component inpr along the pull vector is accounted for in the usual
     * way. Here we account for the component perpendicular to vec.
     */
    double inpr = 0;
    for (int m = 0; m < DIM; m++)
    {
        inpr += spatialData.dr01[m] * spatialData.vec[m];
    }
    /* The torque force works along the component of the distance vector
     * of between the two "usual" pull groups that is perpendicular to
     * the pull vector. The magnitude of this force is the "usual" scale force
     * multiplied with the ratio of the distance between the two "usual" pull
     * groups and the distance between the two groups that define the vector.
     */
    dvec f_perp;
    for (int m = 0; m < DIM; m++)
    {
        f_perp[m] = (spatialData.dr01[m] - inpr * spatialData.vec[m]) / spatialData.vec_len * pcrd->scalarForce;
    }

    /* Apply the force to the groups defining the vector using opposite signs */
    apply_forces_grp(&pull->group[pcrd->params.group[2]], md, f_perp, -1, f, pull->nthreads);
    apply_forces_grp(&pull->group[pcrd->params.group[3]], md, f_perp, 1, f, pull->nthreads);
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_coord(struct pull_t*               pull,
                               int                          coord,
                               const PullCoordVectorForces& forces,
                               const t_mdatoms*             md,
                               rvec*                        f)
{
    /* Here it would be more efficient to use one large thread-parallel
     * region instead of potential parallel regions within apply_forces_grp.
     * But there could be overlap between pull groups and this would lead
     * to data races.
     */

    const pull_coord_work_t& pcrd = pull->coord[coord];

    if (pcrd.params.eGeom == epullgCYL)
    {
        apply_forces_cyl_grp(&pull->dyna[coord], pcrd.spatialData.cyl_dev, md, forces.force01,
                             pcrd.scalarForce, -1, f, pull->nthreads);

        /* Sum the force along the vector and the radial force */
        dvec f_tot;
        for (int m = 0; m < DIM; m++)
        {
            f_tot[m] = forces.force01[m] + pcrd.scalarForce * pcrd.spatialData.ffrad[m];
        }
        apply_forces_grp(&pull->group[pcrd.params.group[1]], md, f_tot, 1, f, pull->nthreads);
    }
    else
    {
        if (pcrd.params.eGeom == epullgDIRRELATIVE)
        {
            /* We need to apply the torque forces to the pull groups
             * that define the pull vector.
             */
            apply_forces_vec_torque(pull, &pcrd, md, f);
        }

        if (pull->group[pcrd.params.group[0]].params.nat > 0)
        {
            apply_forces_grp(&pull->group[pcrd.params.group[0]], md, forces.force01, -1, f, pull->nthreads);
        }
        apply_forces_grp(&pull->group[pcrd.params.group[1]], md, forces.force01, 1, f, pull->nthreads);

        if (pcrd.params.ngroup >= 4)
        {
            apply_forces_grp(&pull->group[pcrd.params.group[2]], md, forces.force23, -1, f, pull->nthreads);
            apply_forces_grp(&pull->group[pcrd.params.group[3]], md, forces.force23, 1, f, pull->nthreads);
        }
        if (pcrd.params.ngroup >= 6)
        {
            apply_forces_grp(&pull->group[pcrd.params.group[4]], md, forces.force45, -1, f, pull->nthreads);
            apply_forces_grp(&pull->group[pcrd.params.group[5]], md, forces.force45, 1, f, pull->nthreads);
        }
    }
}

real max_pull_distance2(const pull_coord_work_t* pcrd, const t_pbc* pbc)
{
    /* Note that this maximum distance calculation is more complex than
     * most other cases in GROMACS, since here we have to take care of
     * distance calculations that don't involve all three dimensions.
     * For example, we can use distances that are larger than the
     * box X and Y dimensions for a box that is elongated along Z.
     */

    real max_d2 = GMX_REAL_MAX;

    if (pull_coordinate_is_directional(&pcrd->params))
    {
        /* Directional pulling along along direction pcrd->vec.
         * Calculating the exact maximum distance is complex and bug-prone.
         * So we take a safe approach by not allowing distances that
         * are larger than half the distance between unit cell faces
         * along dimensions involved in pcrd->vec.
         */
        for (int m = 0; m < DIM; m++)
        {
            if (m < pbc->ndim_ePBC && pcrd->spatialData.vec[m] != 0)
            {
                real imageDistance2 = gmx::square(pbc->box[m][m]);
                for (int d = m + 1; d < DIM; d++)
                {
                    imageDistance2 -= gmx::square(pbc->box[d][m]);
                }
                max_d2 = std::min(max_d2, imageDistance2);
            }
        }
    }
    else
    {
        /* Distance pulling along dimensions with pcrd->params.dim[d]==1.
         * We use half the minimum box vector length of the dimensions involved.
         * This is correct for all cases, except for corner cases with
         * triclinic boxes where e.g. dim[XX]=1 and dim[YY]=0 with
         * box[YY][XX]!=0 and box[YY][YY] < box[XX][XX]. But since even
         * in such corner cases the user could get correct results,
         * depending on the details of the setup, we avoid further
         * code complications.
         */
        for (int m = 0; m < DIM; m++)
        {
            if (m < pbc->ndim_ePBC && pcrd->params.dim[m] != 0)
            {
                real imageDistance2 = gmx::square(pbc->box[m][m]);
                for (int d = 0; d < m; d++)
                {
                    if (pcrd->params.dim[d] != 0)
                    {
                        imageDistance2 += gmx::square(pbc->box[m][d]);
                    }
                }
                max_d2 = std::min(max_d2, imageDistance2);
            }
        }
    }

    return 0.25 * max_d2;
}

/* This function returns the distance based on coordinates xg and xref.
 * Note that the pull coordinate struct pcrd is not modified.
 */
static void low_get_pull_coord_dr(const struct pull_t*     pull,
                                  const pull_coord_work_t* pcrd,
                                  const t_pbc*             pbc,
                                  dvec                     xg,
                                  dvec                     xref,
                                  double                   max_dist2,
                                  dvec                     dr)
{
    const pull_group_work_t* pgrp0 = &pull->group[pcrd->params.group[0]];

    /* Only the first group can be an absolute reference, in that case nat=0 */
    if (pgrp0->params.nat == 0)
    {
        for (int m = 0; m < DIM; m++)
        {
            xref[m] = pcrd->params.origin[m];
        }
    }

    dvec xrefr;
    copy_dvec(xref, xrefr);

    dvec dref = { 0, 0, 0 };
    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        for (int m = 0; m < DIM; m++)
        {
            dref[m] = pcrd->value_ref * pcrd->spatialData.vec[m];
        }
        /* Add the reference position, so we use the correct periodic image */
        dvec_inc(xrefr, dref);
    }

    pbc_dx_d(pbc, xg, xrefr, dr);

    bool   directional = pull_coordinate_is_directional(&pcrd->params);
    double dr2         = 0;
    for (int m = 0; m < DIM; m++)
    {
        dr[m] *= pcrd->params.dim[m];
        if (pcrd->params.dim[m] && !(directional && pcrd->spatialData.vec[m] == 0))
        {
            dr2 += dr[m] * dr[m];
        }
    }
    /* Check if we are close to switching to another periodic image */
    if (max_dist2 > 0 && dr2 > 0.98 * 0.98 * max_dist2)
    {
        /* Note that technically there is no issue with switching periodic
         * image, as pbc_dx_d returns the distance to the closest periodic
         * image. However in all cases where periodic image switches occur,
         * the pull results will be useless in practice.
         */
        gmx_fatal(FARGS,
                  "Distance between pull groups %d and %d (%f nm) is larger than 0.49 times the "
                  "box size (%f).\n%s",
                  pcrd->params.group[0], pcrd->params.group[1], sqrt(dr2),
                  sqrt(0.98 * 0.98 * max_dist2), pcrd->params.eGeom == epullgDIR ? "You might want to consider using \"pull-geometry = direction-periodic\" instead.\n" : "");
    }

    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        dvec_inc(dr, dref);
    }
}

/* This function returns the distance based on the contents of the pull struct.
 * pull->coord[coord_ind].dr, and potentially vector, are updated.
 */
static void get_pull_coord_dr(struct pull_t* pull, int coord_ind, const t_pbc* pbc)
{
    pull_coord_work_t*    pcrd        = &pull->coord[coord_ind];
    PullCoordSpatialData& spatialData = pcrd->spatialData;

    double md2;
    /* With AWH pulling we allow for periodic pulling with geometry=direction.
     * TODO: Store a periodicity flag instead of checking for external pull provider.
     */
    if (pcrd->params.eGeom == epullgDIRPBC
        || (pcrd->params.eGeom == epullgDIR && pcrd->params.eType == epullEXTERNAL))
    {
        md2 = -1;
    }
    else
    {
        md2 = static_cast<double>(max_pull_distance2(pcrd, pbc));
    }

    if (pcrd->params.eGeom == epullgDIRRELATIVE)
    {
        /* We need to determine the pull vector */
        dvec vec;
        int  m;

        const pull_group_work_t& pgrp2 = pull->group[pcrd->params.group[2]];
        const pull_group_work_t& pgrp3 = pull->group[pcrd->params.group[3]];

        pbc_dx_d(pbc, pgrp3.x, pgrp2.x, vec);

        for (m = 0; m < DIM; m++)
        {
            vec[m] *= pcrd->params.dim[m];
        }
        spatialData.vec_len = dnorm(vec);
        for (m = 0; m < DIM; m++)
        {
            spatialData.vec[m] = vec[m] / spatialData.vec_len;
        }
        if (debug)
        {
            fprintf(debug, "pull coord %d vector: %6.3f %6.3f %6.3f normalized: %6.3f %6.3f %6.3f\n",
                    coord_ind, vec[XX], vec[YY], vec[ZZ], spatialData.vec[XX], spatialData.vec[YY],
                    spatialData.vec[ZZ]);
        }
    }

    pull_group_work_t* pgrp0 = &pull->group[pcrd->params.group[0]];
    pull_group_work_t* pgrp1 = &pull->group[pcrd->params.group[1]];

    low_get_pull_coord_dr(pull, pcrd, pbc, pgrp1->x,
                          pcrd->params.eGeom == epullgCYL ? pull->dyna[coord_ind].x : pgrp0->x, md2,
                          spatialData.dr01);

    if (pcrd->params.ngroup >= 4)
    {
        pull_group_work_t *pgrp2, *pgrp3;
        pgrp2 = &pull->group[pcrd->params.group[2]];
        pgrp3 = &pull->group[pcrd->params.group[3]];

        low_get_pull_coord_dr(pull, pcrd, pbc, pgrp3->x, pgrp2->x, md2, spatialData.dr23);
    }
    if (pcrd->params.ngroup >= 6)
    {
        pull_group_work_t *pgrp4, *pgrp5;
        pgrp4 = &pull->group[pcrd->params.group[4]];
        pgrp5 = &pull->group[pcrd->params.group[5]];

        low_get_pull_coord_dr(pull, pcrd, pbc, pgrp5->x, pgrp4->x, md2, spatialData.dr45);
    }
}

/* Modify x so that it is periodic in [-pi, pi)
 * It is assumed that x is in [-3pi, 3pi) so that x
 * needs to be shifted by at most one period.
 */
static void make_periodic_2pi(double* x)
{
    if (*x >= M_PI)
    {
        *x -= M_2PI;
    }
    else if (*x < -M_PI)
    {
        *x += M_2PI;
    }
}

/* This function should always be used to modify pcrd->value_ref */
static void low_set_pull_coord_reference_value(pull_coord_work_t* pcrd, int coord_ind, double value_ref)
{
    GMX_ASSERT(pcrd->params.eType != epullEXTERNAL,
               "The pull coord reference value should not be used with type external-potential");

    if (pcrd->params.eGeom == epullgDIST)
    {
        if (value_ref < 0)
        {
            gmx_fatal(FARGS,
                      "Pull reference distance for coordinate %d (%f) needs to be non-negative",
                      coord_ind + 1, value_ref);
        }
    }
    else if (pcrd->params.eGeom == epullgANGLE || pcrd->params.eGeom == epullgANGLEAXIS)
    {
        if (value_ref < 0 || value_ref > M_PI)
        {
            gmx_fatal(FARGS,
                      "Pull reference angle for coordinate %d (%f) needs to be in the allowed "
                      "interval [0,180] deg",
                      coord_ind + 1,
                      value_ref * pull_conversion_factor_internal2userinput(&pcrd->params));
        }
    }
    else if (pcrd->params.eGeom == epullgDIHEDRAL)
    {
        /* Allow pulling to be periodic for dihedral angles by remapping the reference value to the interval [-pi, pi). */
        make_periodic_2pi(&value_ref);
    }

    pcrd->value_ref = value_ref;
}

static void update_pull_coord_reference_value(pull_coord_work_t* pcrd, int coord_ind, double t)
{
    /* With zero rate the reference value is set initially and doesn't change */
    if (pcrd->params.rate != 0)
    {
        double value_ref = (pcrd->params.init + pcrd->params.rate * t)
                           * pull_conversion_factor_userinput2internal(&pcrd->params);
        low_set_pull_coord_reference_value(pcrd, coord_ind, value_ref);
    }
}

/* Returns the dihedral angle. Updates the plane normal vectors m, n. */
static double get_dihedral_angle_coord(PullCoordSpatialData* spatialData)
{
    double phi, sign;
    dvec   dr32; /* store instead of dr23? */

    dsvmul(-1, spatialData->dr23, dr32);
    dcprod(spatialData->dr01, dr32, spatialData->planevec_m); /* Normal of first plane */
    dcprod(dr32, spatialData->dr45, spatialData->planevec_n); /* Normal of second plane */
    phi = gmx_angle_between_dvecs(spatialData->planevec_m, spatialData->planevec_n);

    /* Note 1: the sign below is opposite of that in the bondeds or Bekker 1994
     * because there r_ij = ri - rj, while here dr01 = r_1 - r_0
     * Note 2: the angle between the plane normal vectors equals pi only when
     * both planes coincide. Thus, when phi = pi dr01 will lie in both planes and
     * we get a positive sign below. Thus, the range of the dihedral angle is (-180, 180].
     */
    sign = (diprod(spatialData->dr01, spatialData->planevec_n) < 0.0) ? 1.0 : -1.0;
    return sign * phi;
}

/* Calculates pull->coord[coord_ind].value.
 * This function also updates pull->coord[coord_ind].dr.
 */
static void get_pull_coord_distance(struct pull_t* pull, int coord_ind, const t_pbc* pbc)
{
    pull_coord_work_t* pcrd;
    int                m;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_dr(pull, coord_ind, pbc);

    PullCoordSpatialData& spatialData = pcrd->spatialData;

    switch (pcrd->params.eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            spatialData.value = dnorm(spatialData.dr01);
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgDIRRELATIVE:
        case epullgCYL:
            /* Pull along vec */
            spatialData.value = 0;
            for (m = 0; m < DIM; m++)
            {
                spatialData.value += spatialData.vec[m] * spatialData.dr01[m];
            }
            break;
        case epullgANGLE:
            spatialData.value = gmx_angle_between_dvecs(spatialData.dr01, spatialData.dr23);
            break;
        case epullgDIHEDRAL: spatialData.value = get_dihedral_angle_coord(&spatialData); break;
        case epullgANGLEAXIS:
            spatialData.value = gmx_angle_between_dvecs(spatialData.dr01, spatialData.vec);
            break;
        default: gmx_incons("Unsupported pull type in get_pull_coord_distance");
    }
}

/* Returns the deviation from the reference value.
 * Updates pull->coord[coord_ind].dr, .value and .value_ref.
 */
static double get_pull_coord_deviation(struct pull_t* pull, int coord_ind, const t_pbc* pbc, double t)
{
    pull_coord_work_t* pcrd;
    double             dev = 0;

    pcrd = &pull->coord[coord_ind];

    /* Update the reference value before computing the distance,
     * since it is used in the distance computation with periodic pulling.
     */
    update_pull_coord_reference_value(pcrd, coord_ind, t);

    get_pull_coord_distance(pull, coord_ind, pbc);

    get_pull_coord_distance(pull, coord_ind, pbc);

    /* Determine the deviation */
    dev = pcrd->spatialData.value - pcrd->value_ref;

    /* Check that values are allowed */
    if (pcrd->params.eGeom == epullgDIST && pcrd->spatialData.value == 0)
    {
        /* With no vector we can not determine the direction for the force,
         * so we set the force to zero.
         */
        dev = 0;
    }
    else if (pcrd->params.eGeom == epullgDIHEDRAL)
    {
        /* The reference value is in [-pi, pi). The coordinate value is in (-pi, pi].
           Thus, the unwrapped deviation is here in (-2pi, 2pi].
           After making it periodic, the deviation will be in [-pi, pi). */
        make_periodic_2pi(&dev);
    }

    return dev;
}

double get_pull_coord_value(struct pull_t* pull, int coord_ind, const t_pbc* pbc)
{
    get_pull_coord_distance(pull, coord_ind, pbc);

    return pull->coord[coord_ind].spatialData.value;
}

void clear_pull_forces(pull_t* pull)
{
    /* Zeroing the forces is only required for constraint pulling.
     * It can happen that multiple constraint steps need to be applied
     * and therefore the constraint forces need to be accumulated.
     */
    for (pull_coord_work_t& coord : pull->coord)
    {
        coord.scalarForce = 0;
    }
}

/* Apply constraint using SHAKE */
static void
do_constraint(struct pull_t* pull, t_pbc* pbc, rvec* x, rvec* v, gmx_bool bMaster, tensor vir, double dt, double t)
{

    dvec*    r_ij;   /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
    dvec     unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */
    dvec*    rnew;   /* current 'new' positions of the groups */
    double*  dr_tot; /* the total update of the coords */
    dvec     vec;
    double   inpr;
    double   lambda, rm, invdt = 0;
    gmx_bool bConverged_all, bConverged = FALSE;
    int      niter = 0, ii, j, m, max_iter = 100;
    double   a;
    dvec     tmp, tmp3;

    snew(r_ij, pull->coord.size());
    snew(dr_tot, pull->coord.size());

    snew(rnew, pull->group.size());

    /* copy the current unconstrained positions for use in iterations. We
       iterate until rinew[i] and rjnew[j] obey the constraints. Then
       rinew - pull.x_unc[i] is the correction dr to group i */
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        copy_dvec(pull->group[g].xp, rnew[g]);
    }

    /* Determine the constraint directions from the old positions */
    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        pull_coord_work_t* pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType != epullCONSTRAINT)
        {
            continue;
        }

        /* Note that get_pull_coord_distance also sets pcrd->dr and pcrd->value.
         * We don't modify dr and value anymore, so these values are also used
         * for printing.
         */
        get_pull_coord_distance(pull, c, pbc);

        const PullCoordSpatialData& spatialData = pcrd->spatialData;
        if (debug)
        {
            fprintf(debug, "Pull coord %zu dr %f %f %f\n", c, spatialData.dr01[XX],
                    spatialData.dr01[YY], spatialData.dr01[ZZ]);
        }

        if (pcrd->params.eGeom == epullgDIR || pcrd->params.eGeom == epullgDIRPBC)
        {
            /* Select the component along vec */
            a = 0;
            for (m = 0; m < DIM; m++)
            {
                a += spatialData.vec[m] * spatialData.dr01[m];
            }
            for (m = 0; m < DIM; m++)
            {
                r_ij[c][m] = a * spatialData.vec[m];
            }
        }
        else
        {
            copy_dvec(spatialData.dr01, r_ij[c]);
        }

        if (dnorm2(r_ij[c]) == 0)
        {
            gmx_fatal(FARGS,
                      "Distance for pull coordinate %zu is zero with constraint pulling, which is "
                      "not allowed.",
                      c + 1);
        }
    }

    bConverged_all = FALSE;
    while (!bConverged_all && niter < max_iter)
    {
        bConverged_all = TRUE;

        /* loop over all constraints */
        for (size_t c = 0; c < pull->coord.size(); c++)
        {
            pull_coord_work_t* pcrd;
            pull_group_work_t *pgrp0, *pgrp1;
            dvec               dr0, dr1;

            pcrd = &pull->coord[c];

            if (pcrd->params.eType != epullCONSTRAINT)
            {
                continue;
            }

            update_pull_coord_reference_value(pcrd, c, t);

            pgrp0 = &pull->group[pcrd->params.group[0]];
            pgrp1 = &pull->group[pcrd->params.group[1]];

            /* Get the current difference vector */
            low_get_pull_coord_dr(pull, pcrd, pbc, rnew[pcrd->params.group[1]],
                                  rnew[pcrd->params.group[0]], -1, unc_ij);

            if (debug)
            {
                fprintf(debug, "Pull coord %zu, iteration %d\n", c, niter);
            }

            rm = 1.0 / (pgrp0->invtm + pgrp1->invtm);

            switch (pcrd->params.eGeom)
            {
                case epullgDIST:
                    if (pcrd->value_ref <= 0)
                    {
                        gmx_fatal(
                                FARGS,
                                "The pull constraint reference distance for group %zu is <= 0 (%f)",
                                c, pcrd->value_ref);
                    }

                    {
                        double q, c_a, c_b, c_c;

                        c_a = diprod(r_ij[c], r_ij[c]);
                        c_b = diprod(unc_ij, r_ij[c]) * 2;
                        c_c = diprod(unc_ij, unc_ij) - gmx::square(pcrd->value_ref);

                        if (c_b < 0)
                        {
                            q      = -0.5 * (c_b - std::sqrt(c_b * c_b - 4 * c_a * c_c));
                            lambda = -q / c_a;
                        }
                        else
                        {
                            q      = -0.5 * (c_b + std::sqrt(c_b * c_b - 4 * c_a * c_c));
                            lambda = -c_c / q;
                        }

                        if (debug)
                        {
                            fprintf(debug, "Pull ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n", c_a, c_b,
                                    c_c, lambda);
                        }
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda * rm * pgrp1->invtm, r_ij[c], dr1);
                    dsvmul(lambda * rm * pgrp0->invtm, r_ij[c], dr0);
                    dr_tot[c] += -lambda * dnorm(r_ij[c]);
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    /* A 1-dimensional constraint along a vector */
                    a = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pcrd->spatialData.vec[m];
                        a += unc_ij[m] * vec[m];
                    }
                    /* Select only the component along the vector */
                    dsvmul(a, vec, unc_ij);
                    lambda = a - pcrd->value_ref;
                    if (debug)
                    {
                        fprintf(debug, "Pull inpr %e lambda: %e\n", a, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda * rm * pgrp1->invtm, vec, dr1);
                    dsvmul(lambda * rm * pgrp0->invtm, vec, dr0);
                    dr_tot[c] += -lambda;
                    break;
                default: gmx_incons("Invalid enumeration value for eGeom");
            }

            /* DEBUG */
            if (debug)
            {
                int g0, g1;

                g0 = pcrd->params.group[0];
                g1 = pcrd->params.group[1];
                low_get_pull_coord_dr(pull, pcrd, pbc, rnew[g1], rnew[g0], -1, tmp);
                low_get_pull_coord_dr(pull, pcrd, pbc, dr1, dr0, -1, tmp3);
                fprintf(debug, "Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n", rnew[g0][0],
                        rnew[g0][1], rnew[g0][2], rnew[g1][0], rnew[g1][1], rnew[g1][2], dnorm(tmp));
                fprintf(debug, "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f\n", "", "", "", "", "",
                        "", pcrd->value_ref);
                fprintf(debug, "Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n", dr0[0],
                        dr0[1], dr0[2], dr1[0], dr1[1], dr1[2], dnorm(tmp3));
            } /* END DEBUG */

            /* Update the COMs with dr */
            dvec_inc(rnew[pcrd->params.group[1]], dr1);
            dvec_inc(rnew[pcrd->params.group[0]], dr0);
        }

        /* Check if all constraints are fullfilled now */
        for (pull_coord_work_t& coord : pull->coord)
        {
            if (coord.params.eType != epullCONSTRAINT)
            {
                continue;
            }

            low_get_pull_coord_dr(pull, &coord, pbc, rnew[coord.params.group[1]],
                                  rnew[coord.params.group[0]], -1, unc_ij);

            switch (coord.params.eGeom)
            {
                case epullgDIST:
                    bConverged = fabs(dnorm(unc_ij) - coord.value_ref) < pull->params.constr_tol;
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = coord.spatialData.vec[m];
                    }
                    inpr = diprod(unc_ij, vec);
                    dsvmul(inpr, vec, unc_ij);
                    bConverged = fabs(diprod(unc_ij, vec) - coord.value_ref) < pull->params.constr_tol;
                    break;
            }

            if (!bConverged)
            {
                if (debug)
                {
                    fprintf(debug,
                            "Pull constraint not converged: "
                            "groups %d %d,"
                            "d_ref = %f, current d = %f\n",
                            coord.params.group[0], coord.params.group[1], coord.value_ref, dnorm(unc_ij));
                }

                bConverged_all = FALSE;
            }
        }

        niter++;
        /* if after all constraints are dealt with and bConverged is still TRUE
           we're finished, if not we do another iteration */
    }
    if (niter > max_iter)
    {
        gmx_fatal(FARGS, "Too many iterations for constraint run: %d", niter);
    }

    /* DONE ITERATING, NOW UPDATE COORDINATES AND CALC. CONSTRAINT FORCES */

    if (v)
    {
        invdt = 1 / dt;
    }

    /* update atoms in the groups */
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        const pull_group_work_t* pgrp;
        dvec                     dr;

        pgrp = &pull->group[g];

        /* get the final constraint displacement dr for group g */
        dvec_sub(rnew[g], pgrp->xp, dr);

        if (dnorm2(dr) == 0)
        {
            /* No displacement, no update necessary */
            continue;
        }

        /* update the atom positions */
        auto localAtomIndices = pgrp->atomSet.localIndex();
        copy_dvec(dr, tmp);
        for (gmx::index j = 0; j < localAtomIndices.ssize(); j++)
        {
            ii = localAtomIndices[j];
            if (!pgrp->localWeights.empty())
            {
                dsvmul(pgrp->wscale * pgrp->localWeights[j], dr, tmp);
            }
            for (m = 0; m < DIM; m++)
            {
                x[ii][m] += tmp[m];
            }
            if (v)
            {
                for (m = 0; m < DIM; m++)
                {
                    v[ii][m] += invdt * tmp[m];
                }
            }
        }
    }

    /* calculate the constraint forces, used for output and virial only */
    for (size_t c = 0; c < pull->coord.size(); c++)
    {
        pull_coord_work_t* pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType != epullCONSTRAINT)
        {
            continue;
        }

        /* Accumulate the forces, in case we have multiple constraint steps */
        double force =
                dr_tot[c]
                / ((pull->group[pcrd->params.group[0]].invtm + pull->group[pcrd->params.group[1]].invtm)
                   * dt * dt);
        pcrd->scalarForce += force;

        if (vir != nullptr && pcrd->params.eGeom != epullgDIRPBC && bMaster)
        {
            double f_invr;

            /* Add the pull contribution to the virial */
            /* We have already checked above that r_ij[c] != 0 */
            f_invr = pcrd->scalarForce / dnorm(r_ij[c]);

            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5 * f_invr * r_ij[c][j] * r_ij[c][m];
                }
            }
        }
    }

    /* finished! I hope. Give back some memory */
    sfree(r_ij);
    sfree(dr_tot);
    sfree(rnew);
}

static void add_virial_coord_dr(tensor vir, const dvec dr, const dvec f)
{
    for (int j = 0; j < DIM; j++)
    {
        for (int m = 0; m < DIM; m++)
        {
            vir[j][m] -= 0.5 * f[j] * dr[m];
        }
    }
}

/* Adds the pull contribution to the virial */
static void add_virial_coord(tensor vir, const pull_coord_work_t& pcrd, const PullCoordVectorForces& forces)
{
    if (vir != nullptr && pcrd.params.eGeom != epullgDIRPBC)
    {
        /* Add the pull contribution for each distance vector to the virial. */
        add_virial_coord_dr(vir, pcrd.spatialData.dr01, forces.force01);
        if (pcrd.params.ngroup >= 4)
        {
            add_virial_coord_dr(vir, pcrd.spatialData.dr23, forces.force23);
        }
        if (pcrd.params.ngroup >= 6)
        {
            add_virial_coord_dr(vir, pcrd.spatialData.dr45, forces.force45);
        }
    }
}

static void calc_pull_coord_scalar_force_and_potential(pull_coord_work_t* pcrd,
                                                       double             dev,
                                                       real               lambda,
                                                       real*              V,
                                                       real*              dVdl)
{
    real k, dkdl;

    k    = (1.0 - lambda) * pcrd->params.k + lambda * pcrd->params.kB;
    dkdl = pcrd->params.kB - pcrd->params.k;

    switch (pcrd->params.eType)
    {
        case epullUMBRELLA:
        case epullFLATBOTTOM:
        case epullFLATBOTTOMHIGH:
            /* The only difference between an umbrella and a flat-bottom
             * potential is that a flat-bottom is zero above or below
               the reference value.
             */
            if ((pcrd->params.eType == epullFLATBOTTOM && dev < 0)
                || (pcrd->params.eType == epullFLATBOTTOMHIGH && dev > 0))
            {
                dev = 0;
            }

            pcrd->scalarForce = -k * dev;
            *V += 0.5 * k * gmx::square(dev);
            *dVdl += 0.5 * dkdl * gmx::square(dev);
            break;
        case epullCONST_F:
            pcrd->scalarForce = -k;
            *V += k * pcrd->spatialData.value;
            *dVdl += dkdl * pcrd->spatialData.value;
            break;
        case epullEXTERNAL:
            gmx_incons(
                    "the scalar pull force should not be calculated internally for pull type "
                    "external");
        default: gmx_incons("Unsupported pull type in do_pull_pot");
    }
}

static PullCoordVectorForces calculateVectorForces(const pull_coord_work_t& pcrd)
{
    const t_pull_coord&         params      = pcrd.params;
    const PullCoordSpatialData& spatialData = pcrd.spatialData;

    /* The geometry of the coordinate determines how the scalar force relates to the force on each group */
    PullCoordVectorForces forces;

    if (params.eGeom == epullgDIST)
    {
        double invdr01 = spatialData.value > 0 ? 1. / spatialData.value : 0.;
        for (int m = 0; m < DIM; m++)
        {
            forces.force01[m] = pcrd.scalarForce * spatialData.dr01[m] * invdr01;
        }
    }
    else if (params.eGeom == epullgANGLE)
    {

        double cos_theta, cos_theta2;

        cos_theta  = cos(spatialData.value);
        cos_theta2 = gmx::square(cos_theta);

        /* The force at theta = 0, pi is undefined so we should not apply any force.
         * cos(theta) = 1 for theta = 0, pi so this is what we check for below.
         * Note that dr01 = 0 or dr23 = 0 evaluates to theta = 0 so we do not
         * have to check for this before dividing by their norm below.
         */
        if (cos_theta2 < 1)
        {
            /* The forces f01 and f23 can be decomposed into one vector parallel to dr01
             * and another vector parallel to dr23:
             * f01 = -dV/dtheta*(a*dr23/|dr23| - b*dr01*|dr01|)/|dr01|,
             * f23 = -dV/dtheta*(a*dr01/|dr01| - b*dr23*|dr23|)/|dr23|,
             */
            double a       = -gmx::invsqrt(1 - cos_theta2); /* comes from d/dx acos(x) */
            double b       = a * cos_theta;
            double invdr01 = 1. / dnorm(spatialData.dr01);
            double invdr23 = 1. / dnorm(spatialData.dr23);
            dvec   normalized_dr01, normalized_dr23;
            dsvmul(invdr01, spatialData.dr01, normalized_dr01);
            dsvmul(invdr23, spatialData.dr23, normalized_dr23);

            for (int m = 0; m < DIM; m++)
            {
                /* Here, f_scal is -dV/dtheta */
                forces.force01[m] =
                        pcrd.scalarForce * invdr01 * (a * normalized_dr23[m] - b * normalized_dr01[m]);
                forces.force23[m] =
                        pcrd.scalarForce * invdr23 * (a * normalized_dr01[m] - b * normalized_dr23[m]);
            }
        }
        else
        {
            /* No forces to apply for ill-defined cases*/
            clear_dvec(forces.force01);
            clear_dvec(forces.force23);
        }
    }
    else if (params.eGeom == epullgANGLEAXIS)
    {
        double cos_theta, cos_theta2;

        /* The angle-axis force is exactly the same as the angle force (above) except that in
           this case the second vector (dr23) is replaced by the pull vector. */
        cos_theta  = cos(spatialData.value);
        cos_theta2 = gmx::square(cos_theta);

        if (cos_theta2 < 1)
        {
            double a, b;
            double invdr01;
            dvec   normalized_dr01;

            invdr01 = 1. / dnorm(spatialData.dr01);
            dsvmul(invdr01, spatialData.dr01, normalized_dr01);
            a = -gmx::invsqrt(1 - cos_theta2); /* comes from d/dx acos(x) */
            b = a * cos_theta;

            for (int m = 0; m < DIM; m++)
            {
                forces.force01[m] =
                        pcrd.scalarForce * invdr01 * (a * spatialData.vec[m] - b * normalized_dr01[m]);
            }
        }
        else
        {
            clear_dvec(forces.force01);
        }
    }
    else if (params.eGeom == epullgDIHEDRAL)
    {
        double m2, n2, tol, sqrdist_32;
        dvec   dr32;
        /* Note: there is a small difference here compared to the
           dihedral force calculations in the bondeds (ref: Bekker 1994).
           There rij = ri - rj, while here dr01 = r1 - r0.
           However, all distance vectors occur in form of cross or inner products
           so that two signs cancel and we end up with the same expressions.
           Also, we treat the more general case of 6 groups (0..5) instead of 4 (i, j, k, l).
         */
        m2 = diprod(spatialData.planevec_m, spatialData.planevec_m);
        n2 = diprod(spatialData.planevec_n, spatialData.planevec_n);
        dsvmul(-1, spatialData.dr23, dr32);
        sqrdist_32 = diprod(dr32, dr32);
        tol        = sqrdist_32 * GMX_REAL_EPS; /* Avoid tiny angles */
        if ((m2 > tol) && (n2 > tol))
        {
            double a_01, a_23_01, a_23_45, a_45;
            double inv_dist_32, inv_sqrdist_32, dist_32;
            dvec   u, v;
            inv_dist_32    = gmx::invsqrt(sqrdist_32);
            inv_sqrdist_32 = inv_dist_32 * inv_dist_32;
            dist_32        = sqrdist_32 * inv_dist_32;

            /* Forces on groups 0, 1 */
            a_01 = pcrd.scalarForce * dist_32 / m2; /* scalarForce is -dV/dphi */
            dsvmul(-a_01, spatialData.planevec_m,
                   forces.force01); /* added sign to get force on group 1, not 0 */

            /* Forces on groups 4, 5 */
            a_45 = -pcrd.scalarForce * dist_32 / n2;
            dsvmul(a_45, spatialData.planevec_n, forces.force45); /* force on group 5 */

            /* Force on groups 2, 3 (defining the axis) */
            a_23_01 = -diprod(spatialData.dr01, dr32) * inv_sqrdist_32;
            a_23_45 = -diprod(spatialData.dr45, dr32) * inv_sqrdist_32;
            dsvmul(-a_23_01, forces.force01, u); /* added sign to get force from group 0, not 1 */
            dsvmul(a_23_45, forces.force45, v);
            dvec_sub(u, v, forces.force23); /* force on group 3 */
        }
        else
        {
            /* No force to apply for ill-defined cases */
            clear_dvec(forces.force01);
            clear_dvec(forces.force23);
            clear_dvec(forces.force45);
        }
    }
    else
    {
        for (int m = 0; m < DIM; m++)
        {
            forces.force01[m] = pcrd.scalarForce * spatialData.vec[m];
        }
    }

    return forces;
}


/* We use a global mutex for locking access to the pull data structure
 * during registration of external pull potential providers.
 * We could use a different, local mutex for each pull object, but the overhead
 * is extremely small here and registration is only done during initialization.
 */
static gmx::Mutex registrationMutex;

using Lock = gmx::lock_guard<gmx::Mutex>;

void register_external_pull_potential(struct pull_t* pull, int coord_index, const char* provider)
{
    GMX_RELEASE_ASSERT(pull != nullptr, "register_external_pull_potential called before init_pull");
    GMX_RELEASE_ASSERT(provider != nullptr,
                       "register_external_pull_potential called with NULL as provider name");

    if (coord_index < 0 || coord_index >= gmx::ssize(pull->coord))
    {
        gmx_fatal(FARGS,
                  "Module '%s' attempted to register an external potential for pull coordinate %d "
                  "which is out of the pull coordinate range %d - %zu\n",
                  provider, coord_index + 1, 1, pull->coord.size());
    }

    pull_coord_work_t* pcrd = &pull->coord[coord_index];

    if (pcrd->params.eType != epullEXTERNAL)
    {
        gmx_fatal(
                FARGS,
                "Module '%s' attempted to register an external potential for pull coordinate %d "
                "which of type '%s', whereas external potentials are only supported with type '%s'",
                provider, coord_index + 1, epull_names[pcrd->params.eType], epull_names[epullEXTERNAL]);
    }

    GMX_RELEASE_ASSERT(pcrd->params.externalPotentialProvider != nullptr,
                       "The external potential provider string for a pull coordinate is NULL");

    if (gmx_strcasecmp(provider, pcrd->params.externalPotentialProvider) != 0)
    {
        gmx_fatal(FARGS,
                  "Module '%s' attempted to register an external potential for pull coordinate %d "
                  "which expects the external potential to be provided by a module named '%s'",
                  provider, coord_index + 1, pcrd->params.externalPotentialProvider);
    }

    /* Lock to avoid (extremely unlikely) simultaneous reading and writing of
     * pcrd->bExternalPotentialProviderHasBeenRegistered and
     * pull->numUnregisteredExternalPotentials.
     */
    Lock registrationLock(registrationMutex);

    if (pcrd->bExternalPotentialProviderHasBeenRegistered)
    {
        gmx_fatal(FARGS,
                  "Module '%s' attempted to register an external potential for pull coordinate %d "
                  "more than once",
                  provider, coord_index + 1);
    }

    pcrd->bExternalPotentialProviderHasBeenRegistered = true;
    pull->numUnregisteredExternalPotentials--;

    GMX_RELEASE_ASSERT(pull->numUnregisteredExternalPotentials >= 0,
                       "Negative unregistered potentials, the pull code is inconsistent");
}


static void check_external_potential_registration(const struct pull_t* pull)
{
    if (pull->numUnregisteredExternalPotentials > 0)
    {
        size_t c;
        for (c = 0; c < pull->coord.size(); c++)
        {
            if (pull->coord[c].params.eType == epullEXTERNAL
                && !pull->coord[c].bExternalPotentialProviderHasBeenRegistered)
            {
                break;
            }
        }

        GMX_RELEASE_ASSERT(c < pull->coord.size(),
                           "Internal inconsistency in the pull potential provider counting");

        gmx_fatal(FARGS,
                  "No external provider for external pull potentials have been provided for %d "
                  "pull coordinates. The first coordinate without provider is number %zu, which "
                  "expects a module named '%s' to provide the external potential.",
                  pull->numUnregisteredExternalPotentials, c + 1,
                  pull->coord[c].params.externalPotentialProvider);
    }
}


/* Pull takes care of adding the forces of the external potential.
 * The external potential module  has to make sure that the corresponding
 * potential energy is added either to the pull term or to a term
 * specific to the external module.
 */
void apply_external_pull_coord_force(struct pull_t*        pull,
                                     int                   coord_index,
                                     double                coord_force,
                                     const t_mdatoms*      mdatoms,
                                     gmx::ForceWithVirial* forceWithVirial)
{
    pull_coord_work_t* pcrd;

    GMX_ASSERT(coord_index >= 0 && coord_index < gmx::ssize(pull->coord),
               "apply_external_pull_coord_force called with coord_index out of range");

    if (pull->comm.bParticipate)
    {
        pcrd = &pull->coord[coord_index];

        GMX_RELEASE_ASSERT(
                pcrd->params.eType == epullEXTERNAL,
                "The pull force can only be set externally on pull coordinates of external type");

        GMX_ASSERT(pcrd->bExternalPotentialProviderHasBeenRegistered,
                   "apply_external_pull_coord_force called for an unregistered pull coordinate");

        /* Set the force */
        pcrd->scalarForce = coord_force;

        /* Calculate the forces on the pull groups */
        PullCoordVectorForces pullCoordForces = calculateVectorForces(*pcrd);

        /* Add the forces for this coordinate to the total virial and force */
        if (forceWithVirial->computeVirial_ && pull->comm.isMasterRank)
        {
            matrix virial = { { 0 } };
            add_virial_coord(virial, *pcrd, pullCoordForces);
            forceWithVirial->addVirialContribution(virial);
        }

        apply_forces_coord(pull, coord_index, pullCoordForces, mdatoms,
                           as_rvec_array(forceWithVirial->force_.data()));
    }

    pull->numExternalPotentialsStillToBeAppliedThisStep--;
}

/* Calculate the pull potential and scalar force for a pull coordinate.
 * Returns the vector forces for the pull coordinate.
 */
static PullCoordVectorForces do_pull_pot_coord(struct pull_t* pull,
                                               int            coord_ind,
                                               t_pbc*         pbc,
                                               double         t,
                                               real           lambda,
                                               real*          V,
                                               tensor         vir,
                                               real*          dVdl)
{
    pull_coord_work_t& pcrd = pull->coord[coord_ind];

    assert(pcrd.params.eType != epullCONSTRAINT);

    double dev = get_pull_coord_deviation(pull, coord_ind, pbc, t);

    calc_pull_coord_scalar_force_and_potential(&pcrd, dev, lambda, V, dVdl);

    PullCoordVectorForces pullCoordForces = calculateVectorForces(pcrd);

    add_virial_coord(vir, pcrd, pullCoordForces);

    return pullCoordForces;
}

real pull_potential(struct pull_t*        pull,
                    const t_mdatoms*      md,
                    t_pbc*                pbc,
                    const t_commrec*      cr,
                    double                t,
                    real                  lambda,
                    const rvec*           x,
                    gmx::ForceWithVirial* force,
                    real*                 dvdlambda)
{
    real V = 0;

    assert(pull != nullptr);

    /* Ideally we should check external potential registration only during
     * the initialization phase, but that requires another function call
     * that should be done exactly in the right order. So we check here.
     */
    check_external_potential_registration(pull);

    if (pull->comm.bParticipate)
    {
        real dVdl = 0;

        pull_calc_coms(cr, pull, md, pbc, t, x, nullptr);

        rvec*      f             = as_rvec_array(force->force_.data());
        matrix     virial        = { { 0 } };
        const bool computeVirial = (force->computeVirial_ && MASTER(cr));
        for (size_t c = 0; c < pull->coord.size(); c++)
        {
            pull_coord_work_t* pcrd;
            pcrd = &pull->coord[c];

            /* For external potential the force is assumed to be given by an external module by a
               call to apply_pull_coord_external_force */
            if (pcrd->params.eType == epullCONSTRAINT || pcrd->params.eType == epullEXTERNAL)
            {
                continue;
            }

            PullCoordVectorForces pullCoordForces = do_pull_pot_coord(
                    pull, c, pbc, t, lambda, &V, computeVirial ? virial : nullptr, &dVdl);

            /* Distribute the force over the atoms in the pulled groups */
            apply_forces_coord(pull, c, pullCoordForces, md, f);
        }

        if (MASTER(cr))
        {
            force->addVirialContribution(virial);
            *dvdlambda += dVdl;
        }
    }

    GMX_ASSERT(pull->numExternalPotentialsStillToBeAppliedThisStep == 0,
               "Too few or too many external pull potentials have been applied the previous step");
    /* All external pull potentials still need to be applied */
    pull->numExternalPotentialsStillToBeAppliedThisStep = pull->numCoordinatesWithExternalPotential;

    return (MASTER(cr) ? V : 0.0);
}

void pull_constraint(struct pull_t*   pull,
                     const t_mdatoms* md,
                     t_pbc*           pbc,
                     const t_commrec* cr,
                     double           dt,
                     double           t,
                     rvec*            x,
                     rvec*            xp,
                     rvec*            v,
                     tensor           vir)
{
    assert(pull != nullptr);

    if (pull->comm.bParticipate)
    {
        pull_calc_coms(cr, pull, md, pbc, t, x, xp);

        do_constraint(pull, pbc, xp, v, MASTER(cr), vir, dt, t);
    }
}

void dd_make_local_pull_groups(const t_commrec* cr, struct pull_t* pull)
{
    gmx_domdec_t* dd;
    pull_comm_t*  comm;
    gmx_bool      bMustParticipate;

    dd = cr->dd;

    comm = &pull->comm;

    /* We always make the master node participate, such that it can do i/o,
     * add the virial and to simplify MC type extensions people might have.
     */
    bMustParticipate = (comm->bParticipateAll || comm->isMasterRank);

    for (pull_group_work_t& group : pull->group)
    {
        if (!group.globalWeights.empty())
        {
            group.localWeights.resize(group.atomSet.numAtomsLocal());
            for (size_t i = 0; i < group.atomSet.numAtomsLocal(); ++i)
            {
                group.localWeights[i] = group.globalWeights[group.atomSet.collectiveIndex()[i]];
            }
        }

        GMX_ASSERT(bMustParticipate || dd != nullptr,
                   "Either all ranks (including this rank) participate, or we use DD and need to "
                   "have access to dd here");

        /* We should participate if we have pull or pbc atoms */
        if (!bMustParticipate
            && (group.atomSet.numAtomsLocal() > 0
                || (group.epgrppbc == epgrppbcREFAT && group.pbcAtomSet->numAtomsLocal() > 0)))
        {
            bMustParticipate = TRUE;
        }
    }

    if (!comm->bParticipateAll)
    {
        /* Keep currently not required ranks in the communicator
         * if they needed to participate up to 20 decompositions ago.
         * This avoids frequent rebuilds due to atoms jumping back and forth.
         */
        const int64_t history_count = 20;
        gmx_bool      bWillParticipate;
        int           count[2];

        /* Increase the decomposition counter for the current call */
        comm->setup_count++;

        if (bMustParticipate)
        {
            comm->must_count = comm->setup_count;
        }

        bWillParticipate =
                bMustParticipate
                || (comm->bParticipate && comm->must_count >= comm->setup_count - history_count);

        if (debug && dd != nullptr)
        {
            fprintf(debug, "Our DD rank (%3d) pull #atoms>0 or master: %s, will be part %s\n", dd->rank,
                    gmx::boolToString(bMustParticipate), gmx::boolToString(bWillParticipate));
        }

        if (bWillParticipate)
        {
            /* Count the number of ranks that we want to have participating */
            count[0] = 1;
            /* Count the number of ranks that need to be added */
            count[1] = comm->bParticipate ? 0 : 1;
        }
        else
        {
            count[0] = 0;
            count[1] = 0;
        }

        /* The cost of this global operation will be less that the cost
         * of the extra MPI_Comm_split calls that we can avoid.
         */
        gmx_sumi(2, count, cr);

        /* If we are missing ranks or if we have 20% more ranks than needed
         * we make a new sub-communicator.
         */
        if (count[1] > 0 || 6 * count[0] < 5 * comm->nparticipate)
        {
            if (debug)
            {
                fprintf(debug, "Creating new pull subcommunicator of size %d\n", count[0]);
            }

#if GMX_MPI
            if (comm->mpi_comm_com != MPI_COMM_NULL)
            {
                MPI_Comm_free(&comm->mpi_comm_com);
            }

            /* This might be an extremely expensive operation, so we try
             * to avoid this splitting as much as possible.
             */
            assert(dd != nullptr);
            MPI_Comm_split(dd->mpi_comm_all, bWillParticipate ? 0 : 1, dd->rank, &comm->mpi_comm_com);
#endif

            comm->bParticipate = bWillParticipate;
            comm->nparticipate = count[0];

            /* When we use the previous COM for PBC, we need to broadcast
             * the previous COM to ranks that have joined the communicator.
             */
            for (pull_group_work_t& group : pull->group)
            {
                if (group.epgrppbc == epgrppbcPREVSTEPCOM)
                {
                    GMX_ASSERT(comm->bParticipate || !MASTER(cr),
                               "The master rank has to participate, as it should pass an up to "
                               "date prev. COM "
                               "to bcast here as well as to e.g. checkpointing");

                    gmx_bcast(sizeof(dvec), group.x_prev_step, cr->mpi_comm_mygroup);
                }
            }
        }
    }

    /* Since the PBC of atoms might have changed, we need to update the PBC */
    pull->bSetPBCatoms = TRUE;
}

static void init_pull_group_index(FILE*              fplog,
                                  const t_commrec*   cr,
                                  int                g,
                                  pull_group_work_t* pg,
                                  gmx_bool           bConstraint,
                                  const ivec         pulldim_con,
                                  const gmx_mtop_t*  mtop,
                                  const t_inputrec*  ir,
                                  real               lambda)
{
    /* With EM and BD there are no masses in the integrator.
     * But we still want to have the correct mass-weighted COMs.
     * So we store the real masses in the weights.
     */
    const bool setWeights = (pg->params.nweight > 0 || EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD);

    /* In parallel, store we need to extract localWeights from weights at DD time */
    std::vector<real>& weights = ((cr && PAR(cr)) ? pg->globalWeights : pg->localWeights);

    const SimulationGroups& groups = mtop->groups;

    /* Count frozen dimensions and (weighted) mass */
    int    nfrozen = 0;
    double tmass   = 0;
    double wmass   = 0;
    double wwmass  = 0;
    int    molb    = 0;
    for (int i = 0; i < pg->params.nat; i++)
    {
        int ii = pg->params.ind[i];
        if (bConstraint && ir->opts.nFreeze)
        {
            for (int d = 0; d < DIM; d++)
            {
                if (pulldim_con[d] == 1
                    && ir->opts.nFreeze[getGroupType(groups, SimulationAtomGroupType::Freeze, ii)][d])
                {
                    nfrozen++;
                }
            }
        }
        const t_atom& atom = mtopGetAtomParameters(mtop, ii, &molb);
        real          m;
        if (ir->efep == efepNO)
        {
            m = atom.m;
        }
        else
        {
            m = (1 - lambda) * atom.m + lambda * atom.mB;
        }
        real w;
        if (pg->params.nweight > 0)
        {
            w = pg->params.weight[i];
        }
        else
        {
            w = 1;
        }
        if (EI_ENERGY_MINIMIZATION(ir->eI))
        {
            /* Move the mass to the weight */
            w *= m;
            m = 1;
        }
        else if (ir->eI == eiBD)
        {
            real mbd;
            if (ir->bd_fric != 0.0F)
            {
                mbd = ir->bd_fric * ir->delta_t;
            }
            else
            {
                if (groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling].empty())
                {
                    mbd = ir->delta_t / ir->opts.tau_t[0];
                }
                else
                {
                    mbd = ir->delta_t
                          / ir->opts.tau_t[groups.groupNumbers[SimulationAtomGroupType::TemperatureCoupling][ii]];
                }
            }
            w *= m / mbd;
            m = mbd;
        }
        if (setWeights)
        {
            weights.push_back(w);
        }
        tmass += m;
        wmass += m * w;
        wwmass += m * w * w;
    }

    if (wmass == 0)
    {
        /* We can have single atom groups with zero mass with potential pulling
         * without cosine weighting.
         */
        if (pg->params.nat == 1 && !bConstraint && pg->epgrppbc != epgrppbcCOS)
        {
            /* With one atom the mass doesn't matter */
            wwmass = 1;
        }
        else
        {
            gmx_fatal(FARGS, "The total%s mass of pull group %d is zero",
                      pg->params.weight ? " weighted" : "", g);
        }
    }
    if (fplog)
    {
        fprintf(fplog, "Pull group %d: %5d atoms, mass %9.3f", g, pg->params.nat, tmass);
        if (pg->params.weight || EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
        {
            fprintf(fplog, ", weighted mass %9.3f", wmass * wmass / wwmass);
        }
        if (pg->epgrppbc == epgrppbcCOS)
        {
            fprintf(fplog, ", cosine weighting will be used");
        }
        fprintf(fplog, "\n");
    }

    if (nfrozen == 0)
    {
        /* A value != 0 signals not frozen, it is updated later */
        pg->invtm = -1.0;
    }
    else
    {
        int ndim = 0;
        for (int d = 0; d < DIM; d++)
        {
            ndim += pulldim_con[d] * pg->params.nat;
        }
        if (fplog && nfrozen > 0 && nfrozen < ndim)
        {
            fprintf(fplog,
                    "\nWARNING: In pull group %d some, but not all of the degrees of freedom\n"
                    "         that are subject to pulling are frozen.\n"
                    "         For constraint pulling the whole group will be frozen.\n\n",
                    g);
        }
        pg->invtm  = 0.0;
        pg->wscale = 1.0;
    }
}

struct pull_t* init_pull(FILE*                     fplog,
                         const pull_params_t*      pull_params,
                         const t_inputrec*         ir,
                         const gmx_mtop_t*         mtop,
                         const t_commrec*          cr,
                         gmx::LocalAtomSetManager* atomSets,
                         real                      lambda)
{
    struct pull_t* pull;
    pull_comm_t*   comm;

    pull = new pull_t;

    /* Copy the pull parameters */
    pull->params = *pull_params;
    /* Avoid pointer copies */
    pull->params.group = nullptr;
    pull->params.coord = nullptr;

    for (int i = 0; i < pull_params->ngroup; ++i)
    {
        pull->group.emplace_back(pull_params->group[i],
                                 atomSets->add({ pull_params->group[i].ind,
                                                 pull_params->group[i].ind + pull_params->group[i].nat }),
                                 pull_params->bSetPbcRefToPrevStepCOM);
    }

    if (cr != nullptr && DOMAINDECOMP(cr))
    {
        /* Set up the global to local atom mapping for PBC atoms */
        for (pull_group_work_t& group : pull->group)
        {
            if (group.epgrppbc == epgrppbcREFAT || group.epgrppbc == epgrppbcPREVSTEPCOM)
            {
                /* pbcAtomSet consists of a single atom */
                group.pbcAtomSet = std::make_unique<gmx::LocalAtomSet>(
                        atomSets->add({ &group.params.pbcatom, &group.params.pbcatom + 1 }));
            }
        }
    }

    pull->bPotential   = FALSE;
    pull->bConstraint  = FALSE;
    pull->bCylinder    = FALSE;
    pull->bAngle       = FALSE;
    pull->bXOutAverage = pull_params->bXOutAverage;
    pull->bFOutAverage = pull_params->bFOutAverage;

    GMX_RELEASE_ASSERT(pull->group[0].params.nat == 0,
                       "pull group 0 is an absolute reference group and should not contain atoms");

    pull->numCoordinatesWithExternalPotential = 0;

    for (int c = 0; c < pull->params.ncoord; c++)
    {
        /* Construct a pull coordinate, copying all coordinate parameters */
        pull->coord.emplace_back(pull_params->coord[c]);

        pull_coord_work_t* pcrd = &pull->coord.back();

        switch (pcrd->params.eGeom)
        {
            case epullgDIST:
            case epullgDIRRELATIVE: /* Direction vector is determined at each step */
            case epullgANGLE:
            case epullgDIHEDRAL: break;
            case epullgDIR:
            case epullgDIRPBC:
            case epullgCYL:
            case epullgANGLEAXIS:
                copy_rvec_to_dvec(pull_params->coord[c].vec, pcrd->spatialData.vec);
                break;
            default:
                /* We allow reading of newer tpx files with new pull geometries,
                 * but with the same tpx format, with old code. A new geometry
                 * only adds a new enum value, which will be out of range for
                 * old code. The first place we need to generate an error is
                 * here, since the pull code can't handle this.
                 * The enum can be used before arriving here only for printing
                 * the string corresponding to the geometry, which will result
                 * in printing "UNDEFINED".
                 */
                gmx_fatal(FARGS,
                          "Pull geometry not supported for pull coordinate %d. The geometry enum "
                          "%d in the input is larger than that supported by the code (up to %d). "
                          "You are probably reading a tpr file generated with a newer version of "
                          "Gromacs with an binary from an older version of Gromacs.",
                          c + 1, pcrd->params.eGeom, epullgNR - 1);
        }

        if (pcrd->params.eType == epullCONSTRAINT)
        {
            /* Check restrictions of the constraint pull code */
            if (pcrd->params.eGeom == epullgCYL || pcrd->params.eGeom == epullgDIRRELATIVE
                || pcrd->params.eGeom == epullgANGLE || pcrd->params.eGeom == epullgDIHEDRAL
                || pcrd->params.eGeom == epullgANGLEAXIS)
            {
                gmx_fatal(FARGS,
                          "Pulling of type %s can not be combined with geometry %s. Consider using "
                          "pull type %s.",
                          epull_names[pcrd->params.eType], epullg_names[pcrd->params.eGeom],
                          epull_names[epullUMBRELLA]);
            }

            pull->bConstraint = TRUE;
        }
        else
        {
            pull->bPotential = TRUE;
        }

        if (pcrd->params.eGeom == epullgCYL)
        {
            pull->bCylinder = TRUE;
        }
        else if (pcrd->params.eGeom == epullgANGLE || pcrd->params.eGeom == epullgDIHEDRAL
                 || pcrd->params.eGeom == epullgANGLEAXIS)
        {
            pull->bAngle = TRUE;
        }

        /* We only need to calculate the plain COM of a group
         * when it is not only used as a cylinder group.
         * Also the absolute reference group 0 needs no COM computation.
         */
        for (int i = 0; i < pcrd->params.ngroup; i++)
        {
            int groupIndex = pcrd->params.group[i];
            if (groupIndex > 0 && !(pcrd->params.eGeom == epullgCYL && i == 0))
            {
                pull->group[groupIndex].needToCalcCom = true;
            }
        }

        /* With non-zero rate the reference value is set at every step */
        if (pcrd->params.rate == 0)
        {
            /* Initialize the constant reference value */
            if (pcrd->params.eType != epullEXTERNAL)
            {
                low_set_pull_coord_reference_value(
                        pcrd, c,
                        pcrd->params.init * pull_conversion_factor_userinput2internal(&pcrd->params));
            }
            else
            {
                /* With an external pull potential, the reference value loses
                 * it's meaning and should not be used. Setting it to zero
                 * makes any terms dependent on it disappear.
                 * The only issue this causes is that with cylinder pulling
                 * no atoms of the cylinder group within the cylinder radius
                 * should be more than half a box length away from the COM of
                 * the pull group along the axial direction.
                 */
                pcrd->value_ref = 0.0;
            }
        }

        if (pcrd->params.eType == epullEXTERNAL)
        {
            GMX_RELEASE_ASSERT(
                    pcrd->params.rate == 0,
                    "With an external potential, a pull coordinate should have rate = 0");

            /* This external potential needs to be registered later */
            pull->numCoordinatesWithExternalPotential++;
        }
        pcrd->bExternalPotentialProviderHasBeenRegistered = false;
    }

    pull->numUnregisteredExternalPotentials             = pull->numCoordinatesWithExternalPotential;
    pull->numExternalPotentialsStillToBeAppliedThisStep = 0;

    pull->pbcType = ir->pbcType;
    switch (pull->pbcType)
    {
        case PbcType::No: pull->npbcdim = 0; break;
        case PbcType::XY: pull->npbcdim = 2; break;
        default: pull->npbcdim = 3; break;
    }

    if (fplog)
    {
        gmx_bool bAbs, bCos;

        bAbs = FALSE;
        for (const pull_coord_work_t& coord : pull->coord)
        {
            if (pull->group[coord.params.group[0]].params.nat == 0
                || pull->group[coord.params.group[1]].params.nat == 0)
            {
                bAbs = TRUE;
            }
        }

        fprintf(fplog, "\n");
        if (pull->bPotential)
        {
            fprintf(fplog, "Will apply potential COM pulling\n");
        }
        if (pull->bConstraint)
        {
            fprintf(fplog, "Will apply constraint COM pulling\n");
        }
        // Don't include the reference group 0 in output, so we report ngroup-1
        int numRealGroups = pull->group.size() - 1;
        GMX_RELEASE_ASSERT(numRealGroups > 0,
                           "The reference absolute position pull group should always be present");
        fprintf(fplog, "with %zu pull coordinate%s and %d group%s\n", pull->coord.size(),
                pull->coord.size() == 1 ? "" : "s", numRealGroups, numRealGroups == 1 ? "" : "s");
        if (bAbs)
        {
            fprintf(fplog, "with an absolute reference\n");
        }
        bCos = FALSE;
        // Don't include the reference group 0 in loop
        for (size_t g = 1; g < pull->group.size(); g++)
        {
            if (pull->group[g].params.nat > 1 && pull->group[g].params.pbcatom < 0)
            {
                /* We are using cosine weighting */
                fprintf(fplog, "Cosine weighting is used for group %zu\n", g);
                bCos = TRUE;
            }
        }
        if (bCos)
        {
            please_cite(fplog, "Engin2010");
        }
    }

    pull->bRefAt = FALSE;
    pull->cosdim = -1;
    for (size_t g = 0; g < pull->group.size(); g++)
    {
        pull_group_work_t* pgrp;

        pgrp = &pull->group[g];
        if (pgrp->params.nat > 0)
        {
            /* There is an issue when a group is used in multiple coordinates
             * and constraints are applied in different dimensions with atoms
             * frozen in some, but not all dimensions.
             * Since there is only one mass array per group, we can't have
             * frozen/non-frozen atoms for different coords at the same time.
             * But since this is a very exotic case, we don't check for this.
             * A warning is printed in init_pull_group_index.
             */

            gmx_bool bConstraint;
            ivec     pulldim, pulldim_con;

            /* Loop over all pull coordinates to see along which dimensions
             * this group is pulled and if it is involved in constraints.
             */
            bConstraint = FALSE;
            clear_ivec(pulldim);
            clear_ivec(pulldim_con);
            for (const pull_coord_work_t& coord : pull->coord)
            {
                gmx_bool bGroupUsed = FALSE;
                for (int gi = 0; gi < coord.params.ngroup; gi++)
                {
                    if (coord.params.group[gi] == static_cast<int>(g))
                    {
                        bGroupUsed = TRUE;
                    }
                }

                if (bGroupUsed)
                {
                    for (int m = 0; m < DIM; m++)
                    {
                        if (coord.params.dim[m] == 1)
                        {
                            pulldim[m] = 1;

                            if (coord.params.eType == epullCONSTRAINT)
                            {
                                bConstraint    = TRUE;
                                pulldim_con[m] = 1;
                            }
                        }
                    }
                }
            }

            /* Determine if we need to take PBC into account for calculating
             * the COM's of the pull groups.
             */
            switch (pgrp->epgrppbc)
            {
                case epgrppbcREFAT: pull->bRefAt = TRUE; break;
                case epgrppbcCOS:
                    if (pgrp->params.weight != nullptr)
                    {
                        gmx_fatal(FARGS,
                                  "Pull groups can not have relative weights and cosine weighting "
                                  "at same time");
                    }
                    for (int m = 0; m < DIM; m++)
                    {
                        if (m < pull->npbcdim && pulldim[m] == 1)
                        {
                            if (pull->cosdim >= 0 && pull->cosdim != m)
                            {
                                gmx_fatal(FARGS,
                                          "Can only use cosine weighting with pulling in one "
                                          "dimension (use mdp option pull-coord?-dim)");
                            }
                            pull->cosdim = m;
                        }
                    }
                    break;
                case epgrppbcNONE: break;
            }

            /* Set the indices */
            init_pull_group_index(fplog, cr, g, pgrp, bConstraint, pulldim_con, mtop, ir, lambda);
        }
        else
        {
            /* Absolute reference, set the inverse mass to zero.
             * This is only relevant (and used) with constraint pulling.
             */
            pgrp->invtm  = 0;
            pgrp->wscale = 1;
        }
    }

    /* If we use cylinder coordinates, do some initialising for them */
    if (pull->bCylinder)
    {
        for (const pull_coord_work_t& coord : pull->coord)
        {
            if (coord.params.eGeom == epullgCYL)
            {
                if (pull->group[coord.params.group[0]].params.nat == 0)
                {
                    gmx_fatal(FARGS,
                              "A cylinder pull group is not supported when using absolute "
                              "reference!\n");
                }
            }
            const auto& referenceGroup = pull->group[coord.params.group[0]];
            pull->dyna.emplace_back(referenceGroup.params, referenceGroup.atomSet,
                                    pull->params.bSetPbcRefToPrevStepCOM);
        }
    }

    /* The gmx_omp_nthreads module might not be initialized here, so max(1,) */
    pull->nthreads = std::max(1, gmx_omp_nthreads_get(emntDefault));
    pull->comSums.resize(pull->nthreads);

    comm = &pull->comm;

#if GMX_MPI
    /* Use a sub-communicator when we have more than 32 ranks, but not
     * when we have an external pull potential, since then the external
     * potential provider expects each rank to have the coordinate.
     */
    comm->bParticipateAll = (cr == nullptr || !DOMAINDECOMP(cr) || cr->dd->nnodes <= 32
                             || pull->numCoordinatesWithExternalPotential > 0
                             || getenv("GMX_PULL_PARTICIPATE_ALL") != nullptr);
    /* This sub-commicator is not used with comm->bParticipateAll,
     * so we can always initialize it to NULL.
     */
    comm->mpi_comm_com = MPI_COMM_NULL;
    comm->nparticipate = 0;
    comm->isMasterRank = (cr == nullptr || MASTER(cr));
#else
    /* No MPI: 1 rank: all ranks pull */
    comm->bParticipateAll = TRUE;
    comm->isMasterRank    = true;
#endif
    comm->bParticipate = comm->bParticipateAll;
    comm->setup_count  = 0;
    comm->must_count   = 0;

    if (!comm->bParticipateAll && fplog != nullptr)
    {
        fprintf(fplog, "Will use a sub-communicator for pull communication\n");
    }

    comm->pbcAtomBuffer.resize(pull->group.size());
    comm->comBuffer.resize(pull->group.size() * c_comBufferStride);
    if (pull->bCylinder)
    {
        comm->cylinderBuffer.resize(pull->coord.size() * c_cylinderBufferStride);
    }

    /* We still need to initialize the PBC reference coordinates */
    pull->bSetPBCatoms = TRUE;

    pull->out_x = nullptr;
    pull->out_f = nullptr;

    return pull;
}

static void destroy_pull(struct pull_t* pull)
{
#if GMX_MPI
    if (pull->comm.mpi_comm_com != MPI_COMM_NULL)
    {
        MPI_Comm_free(&pull->comm.mpi_comm_com);
    }
#endif

    delete pull;
}

void preparePrevStepPullCom(const t_inputrec* ir,
                            pull_t*           pull_work,
                            const t_mdatoms*  md,
                            t_state*          state,
                            const t_state*    state_global,
                            const t_commrec*  cr,
                            bool              startingFromCheckpoint)
{
    if (!ir->pull || !ir->pull->bSetPbcRefToPrevStepCOM)
    {
        return;
    }
    allocStatePrevStepPullCom(state, pull_work);
    if (startingFromCheckpoint)
    {
        if (MASTER(cr))
        {
            state->pull_com_prev_step = state_global->pull_com_prev_step;
        }
        if (PAR(cr))
        {
            /* Only the master rank has the checkpointed COM from the previous step */
            gmx_bcast(sizeof(double) * state->pull_com_prev_step.size(),
                      &state->pull_com_prev_step[0], cr->mpi_comm_mygroup);
        }
        setPrevStepPullComFromState(pull_work, state);
    }
    else
    {
        t_pbc pbc;
        set_pbc(&pbc, ir->pbcType, state->box);
        initPullComFromPrevStep(cr, pull_work, md, &pbc, state->x.rvec_array());
        updatePrevStepPullCom(pull_work, state);
    }
}

void finish_pull(struct pull_t* pull)
{
    check_external_potential_registration(pull);

    if (pull->out_x)
    {
        gmx_fio_fclose(pull->out_x);
    }
    if (pull->out_f)
    {
        gmx_fio_fclose(pull->out_f);
    }

    destroy_pull(pull);
}

gmx_bool pull_have_potential(const struct pull_t* pull)
{
    return pull->bPotential;
}

gmx_bool pull_have_constraint(const struct pull_t* pull)
{
    return pull->bConstraint;
}

bool pull_have_constraint(const pull_params_t* pullParameters)
{
    if (pullParameters == nullptr)
    {
        return false;
    }
    for (int c = 0; c < pullParameters->ncoord; c++)
    {
        if (pullParameters->coord[c].eType == epullCONSTRAINT)
        {
            return true;
        }
    }
    return false;
}
