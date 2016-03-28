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

#include "pull.h"

#include "config.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull_internal.h"
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

bool pull_coordinate_is_angletype(const t_pull_coord *pcrd)
{
    return (pcrd->eGeom == epullgANGLE ||
            pcrd->eGeom == epullgDIHEDRAL ||
            pcrd->eGeom == epullgANGLEAXIS);
}

static bool pull_coordinate_is_directional(const t_pull_coord *pcrd)
{
    return (pcrd->eGeom == epullgDIR ||
            pcrd->eGeom == epullgDIRPBC ||
            pcrd->eGeom == epullgDIRRELATIVE ||
            pcrd->eGeom == epullgCYL);
}

const char *pull_coordinate_units(const t_pull_coord *pcrd)
{
    return pull_coordinate_is_angletype(pcrd) ?  "deg" : "nm";
}

double pull_conversion_factor_userinput2internal(const t_pull_coord *pcrd)
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

double pull_conversion_factor_internal2userinput(const t_pull_coord *pcrd)
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

static std::string append_before_extension(std::string pathname,
                                           std::string to_append)
{
    /* Appends to_append before last '.' in pathname */
    size_t extPos = pathname.find_last_of('.');
    if (extPos == std::string::npos)
    {
        return pathname + to_append;
    }
    else
    {
        return pathname.substr(0, extPos) + to_append +
               pathname.substr(extPos, std::string::npos);
    }
}

static void pull_print_group_x(FILE *out, ivec dim,
                               const pull_group_work_t *pgrp)
{
    int m;

    for (m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", pgrp->x[m]);
        }
    }
}

static void pull_print_coord_dr_components(FILE *out, const ivec dim, const dvec dr)
{
    for (int m = 0; m < DIM; m++)
    {
        if (dim[m])
        {
            fprintf(out, "\t%g", dr[m]);
        }
    }
}

static void pull_print_coord_dr(FILE *out, const pull_coord_work_t *pcrd,
                                gmx_bool bPrintRefValue,
                                gmx_bool bPrintComponents)
{
    double unit_factor = pull_conversion_factor_internal2userinput(&pcrd->params);

    fprintf(out, "\t%g", pcrd->value*unit_factor);

    if (bPrintRefValue)
    {
        fprintf(out, "\t%g", pcrd->value_ref*unit_factor);
    }

    if (bPrintComponents)
    {
        pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->dr01);
        if (pcrd->params.ngroup >= 4)
        {
            pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->dr23);
        }
        if (pcrd->params.ngroup >= 6)
        {
            pull_print_coord_dr_components(out, pcrd->params.dim, pcrd->dr45);
        }
    }
}

static void pull_print_x(FILE *out, struct pull_t *pull, double t)
{
    int c;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        pull_print_coord_dr(out, pcrd,
                            pull->params.bPrintRefValue && pcrd->params.eType != epullEXTERNAL,
                            pull->params.bPrintComp);

        if (pull->params.bPrintCOM)
        {
            if (pcrd->params.eGeom == epullgCYL)
            {
                pull_print_group_x(out, pcrd->params.dim, &pull->dyna[c]);
            }
            else
            {
                pull_print_group_x(out, pcrd->params.dim,
                                   &pull->group[pcrd->params.group[0]]);
            }
            for (int g = 1; g < pcrd->params.ngroup; g++)
            {
                pull_print_group_x(out, pcrd->params.dim, &pull->group[pcrd->params.group[g]]);
            }
        }
    }
    fprintf(out, "\n");
}

static void pull_print_f(FILE *out, struct pull_t *pull, double t)
{
    int c;

    fprintf(out, "%.4f", t);

    for (c = 0; c < pull->ncoord; c++)
    {
        fprintf(out, "\t%g", pull->coord[c].f_scal);
    }
    fprintf(out, "\n");
}

void pull_print_output(struct pull_t *pull, gmx_int64_t step, double time)
{
    GMX_ASSERT(pull->numExternalPotentialsStillToBeAppliedThisStep == 0, "pull_print_output called before all external pull potentials have been applied");

    if ((pull->params.nstxout != 0) && (step % pull->params.nstxout == 0))
    {
        pull_print_x(pull->out_x, pull, time);
    }

    if ((pull->params.nstfout != 0) && (step % pull->params.nstfout == 0))
    {
        pull_print_f(pull->out_f, pull, time);
    }
}

static void set_legend_for_coord_components(const pull_coord_work_t *pcrd, int coord_index, char **setname, int *nsets_ptr)
{
    /*  Loop over the distance vectors and print their components. Each vector is made up of two consecutive groups. */
    for (int g = 0; g < pcrd->params.ngroup; g += 2)
    {
        /* Loop over the components */
        for (int m  = 0; m < DIM; m++)
        {
            if (pcrd->params.dim[m])
            {
                char legend[STRLEN];

                if (g == 0 && pcrd->params.ngroup <= 2)
                {
                    /*  For the simplest case we print a simplified legend without group indices, just the cooordinate index
                        and which dimensional component it is. */
                    sprintf(legend, "%d d%c", coord_index + 1, 'X' + m);
                }
                else
                {
                    /* Otherwise, print also the group indices (relative to the coordinate, not the global ones). */
                    sprintf(legend, "%d g %d-%d d%c", coord_index + 1, g + 1, g + 2, 'X' + m);
                }

                setname[*nsets_ptr] = gmx_strdup(legend);
                (*nsets_ptr)++;
            }
        }
    }
}

static FILE *open_pull_out(const char *fn, struct pull_t *pull,
                           const gmx_output_env_t *oenv,
                           gmx_bool bCoord,
                           const ContinuationOptions &continuationOptions)
{
    FILE  *fp;
    int    nsets, c, m;
    char **setname, buf[50];

    if (continuationOptions.appendFiles)
    {
        fp = gmx_fio_fopen(fn, "a+");
    }
    else
    {
        fp = gmx_fio_fopen(fn, "w+");
        if (bCoord)
        {
            sprintf(buf, "Position (nm%s)", pull->bAngle ? ", deg" : "");
            xvgr_header(fp, "Pull COM",  "Time (ps)", buf,
                        exvggtXNY, oenv);
        }
        else
        {
            sprintf(buf, "Force (kJ/mol/nm%s)", pull->bAngle ? ", kJ/mol/rad" : "");
            xvgr_header(fp, "Pull force", "Time (ps)", buf,
                        exvggtXNY, oenv);
        }

        /* With default mdp options only the actual coordinate value is printed (1),
         * but optionally the reference value (+ 1),
         * the group COMs for all the groups (+ ngroups_max*DIM)
         * and the components of the distance vectors can be printed (+ (ngroups_max/2)*DIM).
         */
        snew(setname, pull->ncoord*(1 + 1 + c_pullCoordNgroupMax*DIM + c_pullCoordNgroupMax/2*DIM));

        nsets = 0;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (bCoord)
            {
                /* The order of this legend should match the order of printing
                 * the data in print_pull_x above.
                 */

                /* The pull coord distance */
                sprintf(buf, "%d", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
                if (pull->params.bPrintRefValue &&
                    pull->coord[c].params.eType != epullEXTERNAL)
                {
                    sprintf(buf, "%d ref", c+1);
                    setname[nsets] = gmx_strdup(buf);
                    nsets++;
                }
                if (pull->params.bPrintComp)
                {
                    set_legend_for_coord_components(&pull->coord[c], c, setname, &nsets);
                }

                if (pull->params.bPrintCOM)
                {
                    for (int g = 0; g < pull->coord[c].params.ngroup; g++)
                    {
                        /* Legend for reference group position */
                        for (m = 0; m < DIM; m++)
                        {
                            if (pull->coord[c].params.dim[m])
                            {
                                sprintf(buf, "%d g %d %c", c+1, g + 1, 'X'+m);
                                setname[nsets] = gmx_strdup(buf);
                                nsets++;
                            }
                        }
                    }
                }
            }
            else
            {
                /* For the pull force we always only use one scalar */
                sprintf(buf, "%d", c+1);
                setname[nsets] = gmx_strdup(buf);
                nsets++;
            }
        }
        if (nsets > 1)
        {
            xvgr_legend(fp, nsets, (const char**)setname, oenv);
        }
        for (c = 0; c < nsets; c++)
        {
            sfree(setname[c]);
        }
        sfree(setname);
    }

    return fp;
}

/* Apply forces in a mass weighted fashion for part of the pull group */
static void apply_forces_grp_part(const pull_group_work_t *pgrp,
                                  int ind_start, int ind_end,
                                  const t_mdatoms *md,
                                  const dvec f_pull, int sign, rvec *f)
{
    double inv_wm = pgrp->mwscale;

    for (int i = ind_start; i < ind_end; i++)
    {
        int    ii    = pgrp->ind_loc[i];
        double wmass = md->massT[ii];
        if (pgrp->weight_loc)
        {
            wmass *= pgrp->weight_loc[i];
        }

        for (int d = 0; d < DIM; d++)
        {
            f[ii][d] += sign * wmass * f_pull[d] * inv_wm;
        }
    }
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_grp(const pull_group_work_t *pgrp,
                             const t_mdatoms *md,
                             const dvec f_pull, int sign, rvec *f,
                             int nthreads)
{
    if (pgrp->params.nat == 1 && pgrp->nat_loc == 1)
    {
        /* Only one atom and our rank has this atom: we can skip
         * the mass weighting, which means that this code also works
         * for mass=0, e.g. with a virtual site.
         */
        for (int d = 0; d < DIM; d++)
        {
            f[pgrp->ind_loc[0]][d] += sign*f_pull[d];
        }
    }
    else
    {
        if (pgrp->nat_loc <= c_pullMaxNumLocalAtomsSingleThreaded)
        {
            apply_forces_grp_part(pgrp, 0, pgrp->nat_loc, md, f_pull, sign, f);
        }
        else
        {
#pragma omp parallel for num_threads(nthreads) schedule(static)
            for (int th = 0; th < nthreads; th++)
            {
                int ind_start = (pgrp->nat_loc*(th + 0))/nthreads;
                int ind_end   = (pgrp->nat_loc*(th + 1))/nthreads;
                apply_forces_grp_part(pgrp, ind_start, ind_end,
                                      md, f_pull, sign, f);
            }
        }
    }
}

/* Apply forces in a mass weighted fashion to a cylinder group */
static void apply_forces_cyl_grp(const pull_group_work_t *pgrp,
                                 const double dv_corr,
                                 const t_mdatoms *md,
                                 const dvec f_pull, double f_scal,
                                 int sign, rvec *f,
                                 int gmx_unused nthreads)
{
    double inv_wm = pgrp->mwscale;

    /* The cylinder group is always a slab in the system, thus large.
     * Therefore we always thread-parallelize this group.
     */
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int i = 0; i < pgrp->nat_loc; i++)
    {
        int    ii     = pgrp->ind_loc[i];
        double mass   = md->massT[ii];
        double weight = pgrp->weight_loc[i];
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
            f[ii][m] += sign*inv_wm*(mass*weight*f_pull[m] +
                                     pgrp->mdw[i][m]*dv_com*f_scal);
        }
    }
}

/* Apply torque forces in a mass weighted fashion to the groups that define
 * the pull vector direction for pull coordinate pcrd.
 */
static void apply_forces_vec_torque(const struct pull_t     *pull,
                                    const pull_coord_work_t *pcrd,
                                    const t_mdatoms         *md,
                                    rvec                    *f)
{
    double inpr;
    int    m;
    dvec   f_perp;

    /* The component inpr along the pull vector is accounted for in the usual
     * way. Here we account for the component perpendicular to vec.
     */
    inpr = 0;
    for (m = 0; m < DIM; m++)
    {
        inpr += pcrd->dr01[m]*pcrd->vec[m];
    }
    /* The torque force works along the component of the distance vector
     * of between the two "usual" pull groups that is perpendicular to
     * the pull vector. The magnitude of this force is the "usual" scale force
     * multiplied with the ratio of the distance between the two "usual" pull
     * groups and the distance between the two groups that define the vector.
     */
    for (m = 0; m < DIM; m++)
    {
        f_perp[m] = (pcrd->dr01[m] - inpr*pcrd->vec[m])/pcrd->vec_len*pcrd->f_scal;
    }

    /* Apply the force to the groups defining the vector using opposite signs */
    apply_forces_grp(&pull->group[pcrd->params.group[2]], md,
                     f_perp, -1, f, pull->nthreads);
    apply_forces_grp(&pull->group[pcrd->params.group[3]], md,
                     f_perp,  1, f, pull->nthreads);
}

/* Apply forces in a mass weighted fashion */
static void apply_forces_coord(struct pull_t * pull, int coord,
                               const t_mdatoms * md,
                               rvec *f)
{
    /* Here it would be more efficient to use one large thread-parallel
     * region instead of potential parallel regions within apply_forces_grp.
     * But there could be overlap between pull groups and this would lead
     * to data races.
     */

    const pull_coord_work_t *pcrd = &pull->coord[coord];

    if (pcrd->params.eGeom == epullgCYL)
    {
        apply_forces_cyl_grp(&pull->dyna[coord], pcrd->cyl_dev, md,
                             pcrd->f01, pcrd->f_scal, -1, f,
                             pull->nthreads);

        /* Sum the force along the vector and the radial force */
        dvec f_tot;
        for (int m = 0; m < DIM; m++)
        {
            f_tot[m] = pcrd->f01[m] + pcrd->f_scal*pcrd->ffrad[m];
        }
        apply_forces_grp(&pull->group[pcrd->params.group[1]], md,
                         f_tot, 1, f, pull->nthreads);
    }
    else
    {
        if (pcrd->params.eGeom == epullgDIRRELATIVE)
        {
            /* We need to apply the torque forces to the pull groups
             * that define the pull vector.
             */
            apply_forces_vec_torque(pull, pcrd, md, f);
        }

        if (pull->group[pcrd->params.group[0]].params.nat > 0)
        {
            apply_forces_grp(&pull->group[pcrd->params.group[0]], md,
                             pcrd->f01, -1, f, pull->nthreads);
        }
        apply_forces_grp(&pull->group[pcrd->params.group[1]], md,
                         pcrd->f01, 1, f, pull->nthreads);

        if (pcrd->params.ngroup >= 4)
        {
            apply_forces_grp(&pull->group[pcrd->params.group[2]], md,
                             pcrd->f23, -1, f, pull->nthreads);
            apply_forces_grp(&pull->group[pcrd->params.group[3]], md,
                             pcrd->f23,  1, f, pull->nthreads);
        }
        if (pcrd->params.ngroup >= 6)
        {
            apply_forces_grp(&pull->group[pcrd->params.group[4]], md,
                             pcrd->f45, -1, f, pull->nthreads);
            apply_forces_grp(&pull->group[pcrd->params.group[5]], md,
                             pcrd->f45,  1, f, pull->nthreads);
        }
    }
}

real max_pull_distance2(const pull_coord_work_t *pcrd,
                        const t_pbc             *pbc)
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
            if (m < pbc->ndim_ePBC && pcrd->vec[m] != 0)
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

    return 0.25*max_d2;
}

/* This function returns the distance based on coordinates xg and xref.
 * Note that the pull coordinate struct pcrd is not modified.
 */
static void low_get_pull_coord_dr(const struct pull_t *pull,
                                  const pull_coord_work_t *pcrd,
                                  const t_pbc *pbc,
                                  dvec xg, dvec xref, double max_dist2,
                                  dvec dr)
{
    const pull_group_work_t *pgrp0 = &pull->group[pcrd->params.group[0]];

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

    dvec dref = {0, 0, 0};
    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        for (int m = 0; m < DIM; m++)
        {
            dref[m] = pcrd->value_ref*pcrd->vec[m];
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
        if (pcrd->params.dim[m] && !(directional && pcrd->vec[m] == 0))
        {
            dr2 += dr[m]*dr[m];
        }
    }
    /* Check if we are close to switching to another periodic image */
    if (max_dist2 > 0 && dr2 > 0.98*0.98*max_dist2)
    {
        /* Note that technically there is no issue with switching periodic
         * image, as pbc_dx_d returns the distance to the closest periodic
         * image. However in all cases where periodic image switches occur,
         * the pull results will be useless in practice.
         */
        gmx_fatal(FARGS, "Distance between pull groups %d and %d (%f nm) is larger than 0.49 times the box size (%f).\n%s",
                  pcrd->params.group[0], pcrd->params.group[1],
                  sqrt(dr2), sqrt(0.98*0.98*max_dist2),
                  pcrd->params.eGeom == epullgDIR ? "You might want to consider using \"pull-geometry = direction-periodic\" instead.\n" : "");
    }

    if (pcrd->params.eGeom == epullgDIRPBC)
    {
        dvec_inc(dr, dref);
    }
}

/* This function returns the distance based on the contents of the pull struct.
 * pull->coord[coord_ind].dr, and potentially vector, are updated.
 */
static void get_pull_coord_dr(struct pull_t *pull,
                              int            coord_ind,
                              const t_pbc   *pbc)
{
    double             md2;
    pull_coord_work_t *pcrd;
    pull_group_work_t *pgrp0, *pgrp1;

    pcrd = &pull->coord[coord_ind];

    if (pcrd->params.eGeom == epullgDIRPBC)
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
        const pull_group_work_t *pgrp2, *pgrp3;
        dvec                     vec;
        int                      m;

        pgrp2 = &pull->group[pcrd->params.group[2]];
        pgrp3 = &pull->group[pcrd->params.group[3]];

        pbc_dx_d(pbc, pgrp3->x, pgrp2->x, vec);

        for (m = 0; m < DIM; m++)
        {
            vec[m] *= pcrd->params.dim[m];
        }
        pcrd->vec_len = dnorm(vec);
        for (m = 0; m < DIM; m++)
        {
            pcrd->vec[m] = vec[m]/pcrd->vec_len;
        }
        if (debug)
        {
            fprintf(debug, "pull coord %d vector: %6.3f %6.3f %6.3f normalized: %6.3f %6.3f %6.3f\n",
                    coord_ind,
                    vec[XX], vec[YY], vec[ZZ],
                    pcrd->vec[XX], pcrd->vec[YY], pcrd->vec[ZZ]);
        }
    }

    pgrp0 = &pull->group[pcrd->params.group[0]];
    pgrp1 = &pull->group[pcrd->params.group[1]];

    low_get_pull_coord_dr(pull, pcrd, pbc,
                          pgrp1->x,
                          pcrd->params.eGeom == epullgCYL ? pull->dyna[coord_ind].x : pgrp0->x,
                          md2,
                          pcrd->dr01);

    if (pcrd->params.ngroup >= 4)
    {
        pull_group_work_t *pgrp2, *pgrp3;
        pgrp2 = &pull->group[pcrd->params.group[2]];
        pgrp3 = &pull->group[pcrd->params.group[3]];

        low_get_pull_coord_dr(pull, pcrd, pbc,
                              pgrp3->x,
                              pgrp2->x,
                              md2,
                              pcrd->dr23);
    }
    if (pcrd->params.ngroup >= 6)
    {
        pull_group_work_t *pgrp4, *pgrp5;
        pgrp4 = &pull->group[pcrd->params.group[4]];
        pgrp5 = &pull->group[pcrd->params.group[5]];

        low_get_pull_coord_dr(pull, pcrd, pbc,
                              pgrp5->x,
                              pgrp4->x,
                              md2,
                              pcrd->dr45);
    }
}

/* Modify x so that it is periodic in [-pi, pi)
 * It is assumed that x is in [-3pi, 3pi) so that x
 * needs to be shifted by at most one period.
 */
static void make_periodic_2pi(double *x)
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
static void low_set_pull_coord_reference_value(pull_coord_work_t *pcrd,
                                               int                coord_ind,
                                               double             value_ref)
{
    GMX_ASSERT(pcrd->params.eType != epullEXTERNAL, "The pull coord reference value should not be used with type external-potential");

    if (pcrd->params.eGeom == epullgDIST)
    {
        if (value_ref < 0)
        {
            gmx_fatal(FARGS, "Pull reference distance for coordinate %d (%f) needs to be non-negative", coord_ind + 1, value_ref);
        }
    }
    else if (pcrd->params.eGeom == epullgANGLE || pcrd->params.eGeom == epullgANGLEAXIS)
    {
        if (value_ref < 0 || value_ref > M_PI)
        {
            gmx_fatal(FARGS, "Pull reference angle for coordinate %d (%f) needs to be in the allowed interval [0,180] deg", coord_ind + 1, value_ref*pull_conversion_factor_internal2userinput(&pcrd->params));
        }
    }
    else if (pcrd->params.eGeom == epullgDIHEDRAL)
    {
        /* Allow pulling to be periodic for dihedral angles by remapping the reference value to the interval [-pi, pi). */
        make_periodic_2pi(&value_ref);
    }

    pcrd->value_ref = value_ref;
}

static void update_pull_coord_reference_value(pull_coord_work_t *pcrd, int coord_ind, double t)
{
    /* With zero rate the reference value is set initially and doesn't change */
    if (pcrd->params.rate != 0)
    {
        double value_ref = (pcrd->params.init + pcrd->params.rate*t)*pull_conversion_factor_userinput2internal(&pcrd->params);
        low_set_pull_coord_reference_value(pcrd, coord_ind, value_ref);
    }
}

/* Returns the dihedral angle. Updates the plane normal vectors m, n. */
static double get_dihedral_angle_coord(pull_coord_work_t *pcrd)
{
    double phi, sign;
    dvec   dr32; /* store instead of dr23? */

    dsvmul(-1, pcrd->dr23, dr32);
    dcprod(pcrd->dr01, dr32, pcrd->planevec_m);  /* Normal of first plane */
    dcprod(dr32, pcrd->dr45, pcrd->planevec_n);  /* Normal of second plane */
    phi = gmx_angle_between_dvecs(pcrd->planevec_m, pcrd->planevec_n);

    /* Note 1: the sign below is opposite of that in the bondeds or Bekker 1994
     * because there r_ij = ri - rj, while here dr01 = r_1 - r_0
     * Note 2: the angle between the plane normal vectors equals pi only when
     * both planes coincide. Thus, when phi = pi dr01 will lie in both planes and
     * we get a positive sign below. Thus, the range of the dihedral angle is (-180, 180].
     */
    sign = (diprod(pcrd->dr01, pcrd->planevec_n) < 0.0) ? 1.0 : -1.0;
    return sign*phi;
}

/* Calculates pull->coord[coord_ind].value.
 * This function also updates pull->coord[coord_ind].dr.
 */
static void get_pull_coord_distance(struct pull_t *pull,
                                    int            coord_ind,
                                    const t_pbc   *pbc)
{
    pull_coord_work_t *pcrd;
    int                m;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_dr(pull, coord_ind, pbc);

    switch (pcrd->params.eGeom)
    {
        case epullgDIST:
            /* Pull along the vector between the com's */
            pcrd->value = dnorm(pcrd->dr01);
            break;
        case epullgDIR:
        case epullgDIRPBC:
        case epullgDIRRELATIVE:
        case epullgCYL:
            /* Pull along vec */
            pcrd->value = 0;
            for (m = 0; m < DIM; m++)
            {
                pcrd->value += pcrd->vec[m]*pcrd->dr01[m];
            }
            break;
        case epullgANGLE:
            pcrd->value = gmx_angle_between_dvecs(pcrd->dr01, pcrd->dr23);
            break;
        case epullgDIHEDRAL:
            pcrd->value = get_dihedral_angle_coord(pcrd);
            break;
        case epullgANGLEAXIS:
            pcrd->value = gmx_angle_between_dvecs(pcrd->dr01, pcrd->vec);
            break;
        default:
            gmx_incons("Unsupported pull type in get_pull_coord_distance");
    }
}

/* Returns the deviation from the reference value.
 * Updates pull->coord[coord_ind].dr, .value and .value_ref.
 */
static double get_pull_coord_deviation(struct pull_t *pull,
                                       int            coord_ind,
                                       const t_pbc   *pbc,
                                       double         t)
{
    pull_coord_work_t *pcrd;
    double             dev = 0;

    pcrd = &pull->coord[coord_ind];

    get_pull_coord_distance(pull, coord_ind, pbc);

    update_pull_coord_reference_value(pcrd, coord_ind, t);

    /* Determine the deviation */
    dev = pcrd->value - pcrd->value_ref;

    /* Check that values are allowed */
    if (pcrd->params.eGeom == epullgDIST && pcrd->value == 0)
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

double get_pull_coord_value(struct pull_t *pull,
                            int            coord_ind,
                            const t_pbc   *pbc)
{
    get_pull_coord_distance(pull, coord_ind, pbc);

    return pull->coord[coord_ind].value;
}

static void clear_pull_forces_coord(pull_coord_work_t *pcrd)
{
    clear_dvec(pcrd->f01);
    clear_dvec(pcrd->f23);
    clear_dvec(pcrd->f45);
    pcrd->f_scal = 0;
}

void clear_pull_forces(struct pull_t *pull)
{
    int i;

    /* Zeroing the forces is only required for constraint pulling.
     * It can happen that multiple constraint steps need to be applied
     * and therefore the constraint forces need to be accumulated.
     */
    for (i = 0; i < pull->ncoord; i++)
    {
        clear_pull_forces_coord(&pull->coord[i]);
    }
}

/* Apply constraint using SHAKE */
static void do_constraint(struct pull_t *pull, t_pbc *pbc,
                          rvec *x, rvec *v,
                          gmx_bool bMaster, tensor vir,
                          double dt, double t)
{

    dvec         *r_ij;   /* x[i] com of i in prev. step. Obeys constr. -> r_ij[i] */
    dvec          unc_ij; /* xp[i] com of i this step, before constr.   -> unc_ij  */
    dvec         *rnew;   /* current 'new' positions of the groups */
    double       *dr_tot; /* the total update of the coords */
    dvec          vec;
    double        inpr;
    double        lambda, rm, invdt = 0;
    gmx_bool      bConverged_all, bConverged = FALSE;
    int           niter = 0, g, c, ii, j, m, max_iter = 100;
    double        a;
    dvec          tmp, tmp3;

    snew(r_ij,   pull->ncoord);
    snew(dr_tot, pull->ncoord);

    snew(rnew, pull->ngroup);

    /* copy the current unconstrained positions for use in iterations. We
       iterate until rinew[i] and rjnew[j] obey the constraints. Then
       rinew - pull.x_unc[i] is the correction dr to group i */
    for (g = 0; g < pull->ngroup; g++)
    {
        copy_dvec(pull->group[g].xp, rnew[g]);
    }

    /* Determine the constraint directions from the old positions */
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

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

        if (debug)
        {
            fprintf(debug, "Pull coord %d dr %f %f %f\n",
                    c, pcrd->dr01[XX], pcrd->dr01[YY], pcrd->dr01[ZZ]);
        }

        if (pcrd->params.eGeom == epullgDIR ||
            pcrd->params.eGeom == epullgDIRPBC)
        {
            /* Select the component along vec */
            a = 0;
            for (m = 0; m < DIM; m++)
            {
                a += pcrd->vec[m]*pcrd->dr01[m];
            }
            for (m = 0; m < DIM; m++)
            {
                r_ij[c][m] = a*pcrd->vec[m];
            }
        }
        else
        {
            copy_dvec(pcrd->dr01, r_ij[c]);
        }

        if (dnorm2(r_ij[c]) == 0)
        {
            gmx_fatal(FARGS, "Distance for pull coordinate %d is zero with constraint pulling, which is not allowed.", c + 1);
        }
    }

    bConverged_all = FALSE;
    while (!bConverged_all && niter < max_iter)
    {
        bConverged_all = TRUE;

        /* loop over all constraints */
        for (c = 0; c < pull->ncoord; c++)
        {
            pull_coord_work_t *pcrd;
            pull_group_work_t *pgrp0, *pgrp1;
            dvec               dr0, dr1;

            pcrd  = &pull->coord[c];

            if (pcrd->params.eType != epullCONSTRAINT)
            {
                continue;
            }

            update_pull_coord_reference_value(pcrd, c, t);

            pgrp0 = &pull->group[pcrd->params.group[0]];
            pgrp1 = &pull->group[pcrd->params.group[1]];

            /* Get the current difference vector */
            low_get_pull_coord_dr(pull, pcrd, pbc,
                                  rnew[pcrd->params.group[1]],
                                  rnew[pcrd->params.group[0]],
                                  -1, unc_ij);

            if (debug)
            {
                fprintf(debug, "Pull coord %d, iteration %d\n", c, niter);
            }

            rm = 1.0/(pgrp0->invtm + pgrp1->invtm);

            switch (pcrd->params.eGeom)
            {
                case epullgDIST:
                    if (pcrd->value_ref <= 0)
                    {
                        gmx_fatal(FARGS, "The pull constraint reference distance for group %d is <= 0 (%f)", c, pcrd->value_ref);
                    }

                    {
                        double q, c_a, c_b, c_c;

                        c_a = diprod(r_ij[c], r_ij[c]);
                        c_b = diprod(unc_ij, r_ij[c])*2;
                        c_c = diprod(unc_ij, unc_ij) - gmx::square(pcrd->value_ref);

                        if (c_b < 0)
                        {
                            q      = -0.5*(c_b - std::sqrt(c_b*c_b - 4*c_a*c_c));
                            lambda = -q/c_a;
                        }
                        else
                        {
                            q      = -0.5*(c_b + std::sqrt(c_b*c_b - 4*c_a*c_c));
                            lambda = -c_c/q;
                        }

                        if (debug)
                        {
                            fprintf(debug,
                                    "Pull ax^2+bx+c=0: a=%e b=%e c=%e lambda=%e\n",
                                    c_a, c_b, c_c, lambda);
                        }
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp1->invtm, r_ij[c], dr1);
                    dsvmul( lambda*rm*pgrp0->invtm, r_ij[c], dr0);
                    dr_tot[c] += -lambda*dnorm(r_ij[c]);
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    /* A 1-dimensional constraint along a vector */
                    a = 0;
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pcrd->vec[m];
                        a     += unc_ij[m]*vec[m];
                    }
                    /* Select only the component along the vector */
                    dsvmul(a, vec, unc_ij);
                    lambda = a - pcrd->value_ref;
                    if (debug)
                    {
                        fprintf(debug, "Pull inpr %e lambda: %e\n", a, lambda);
                    }

                    /* The position corrections dr due to the constraints */
                    dsvmul(-lambda*rm*pgrp1->invtm, vec, dr1);
                    dsvmul( lambda*rm*pgrp0->invtm, vec, dr0);
                    dr_tot[c] += -lambda;
                    break;
                default:
                    gmx_incons("Invalid enumeration value for eGeom");
                    /* Keep static analyzer happy */
                    clear_dvec(dr0);
                    clear_dvec(dr1);
            }

            /* DEBUG */
            if (debug)
            {
                int g0, g1;

                g0 = pcrd->params.group[0];
                g1 = pcrd->params.group[1];
                low_get_pull_coord_dr(pull, pcrd, pbc, rnew[g1], rnew[g0], -1, tmp);
                low_get_pull_coord_dr(pull, pcrd, pbc, dr1, dr0, -1, tmp3);
                fprintf(debug,
                        "Pull cur %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        rnew[g0][0], rnew[g0][1], rnew[g0][2],
                        rnew[g1][0], rnew[g1][1], rnew[g1][2], dnorm(tmp));
                fprintf(debug,
                        "Pull ref %8s %8s %8s   %8s %8s %8s d: %8.5f\n",
                        "", "", "", "", "", "", pcrd->value_ref);
                fprintf(debug,
                        "Pull cor %8.5f %8.5f %8.5f j:%8.5f %8.5f %8.5f d: %8.5f\n",
                        dr0[0], dr0[1], dr0[2],
                        dr1[0], dr1[1], dr1[2],
                        dnorm(tmp3));
            } /* END DEBUG */

            /* Update the COMs with dr */
            dvec_inc(rnew[pcrd->params.group[1]], dr1);
            dvec_inc(rnew[pcrd->params.group[0]], dr0);
        }

        /* Check if all constraints are fullfilled now */
        for (c = 0; c < pull->ncoord; c++)
        {
            pull_coord_work_t *pcrd;

            pcrd = &pull->coord[c];

            if (pcrd->params.eType != epullCONSTRAINT)
            {
                continue;
            }

            low_get_pull_coord_dr(pull, pcrd, pbc,
                                  rnew[pcrd->params.group[1]],
                                  rnew[pcrd->params.group[0]],
                                  -1, unc_ij);

            switch (pcrd->params.eGeom)
            {
                case epullgDIST:
                    bConverged =
                        fabs(dnorm(unc_ij) - pcrd->value_ref) < pull->params.constr_tol;
                    break;
                case epullgDIR:
                case epullgDIRPBC:
                case epullgCYL:
                    for (m = 0; m < DIM; m++)
                    {
                        vec[m] = pcrd->vec[m];
                    }
                    inpr = diprod(unc_ij, vec);
                    dsvmul(inpr, vec, unc_ij);
                    bConverged =
                        fabs(diprod(unc_ij, vec) - pcrd->value_ref) < pull->params.constr_tol;
                    break;
            }

            if (!bConverged)
            {
                if (debug)
                {
                    fprintf(debug, "NOT CONVERGED YET: Group %d:"
                            "d_ref = %f, current d = %f\n",
                            g, pcrd->value_ref, dnorm(unc_ij));
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
        invdt = 1/dt;
    }

    /* update atoms in the groups */
    for (g = 0; g < pull->ngroup; g++)
    {
        const pull_group_work_t *pgrp;
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
        copy_dvec(dr, tmp);
        for (j = 0; j < pgrp->nat_loc; j++)
        {
            ii = pgrp->ind_loc[j];
            if (pgrp->weight_loc)
            {
                dsvmul(pgrp->wscale*pgrp->weight_loc[j], dr, tmp);
            }
            for (m = 0; m < DIM; m++)
            {
                x[ii][m] += tmp[m];
            }
            if (v)
            {
                for (m = 0; m < DIM; m++)
                {
                    v[ii][m] += invdt*tmp[m];
                }
            }
        }
    }

    /* calculate the constraint forces, used for output and virial only */
    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;

        pcrd = &pull->coord[c];

        if (pcrd->params.eType != epullCONSTRAINT)
        {
            continue;
        }

        /* Accumulate the forces, in case we have multiple constraint steps */
        pcrd->f_scal += dr_tot[c]/((pull->group[pcrd->params.group[0]].invtm + pull->group[pcrd->params.group[1]].invtm)*dt*dt);

        if (vir != nullptr && pcrd->params.eGeom != epullgDIRPBC && bMaster)
        {
            double f_invr;

            /* Add the pull contribution to the virial */
            /* We have already checked above that r_ij[c] != 0 */
            f_invr = pcrd->f_scal/dnorm(r_ij[c]);

            for (j = 0; j < DIM; j++)
            {
                for (m = 0; m < DIM; m++)
                {
                    vir[j][m] -= 0.5*f_invr*r_ij[c][j]*r_ij[c][m];
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
            vir[j][m] -= 0.5*f[j]*dr[m];
        }
    }
}

/* Adds the pull contribution to the virial */
static void add_virial_coord(tensor vir, const pull_coord_work_t *pcrd)
{
    if (vir != nullptr && pcrd->params.eGeom != epullgDIRPBC)
    {
        /* Add the pull contribution for each distance vector to the virial. */
        add_virial_coord_dr(vir, pcrd->dr01, pcrd->f01);
        if (pcrd->params.ngroup >= 4)
        {
            add_virial_coord_dr(vir, pcrd->dr23, pcrd->f23);
        }
        if (pcrd->params.ngroup >= 6)
        {
            add_virial_coord_dr(vir, pcrd->dr45, pcrd->f45);
        }
    }
}

static void calc_pull_coord_scalar_force_and_potential(pull_coord_work_t *pcrd,
                                                       double dev, real lambda,
                                                       real *V, real *dVdl)
{
    real   k, dkdl;

    k    = (1.0 - lambda)*pcrd->params.k + lambda*pcrd->params.kB;
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
            if ((pcrd->params.eType == epullFLATBOTTOM && dev < 0) ||
                (pcrd->params.eType == epullFLATBOTTOMHIGH && dev > 0))
            {
                dev = 0;
            }

            pcrd->f_scal  =       -k*dev;
            *V           += 0.5*   k*gmx::square(dev);
            *dVdl        += 0.5*dkdl*gmx::square(dev);
            break;
        case epullCONST_F:
            pcrd->f_scal  =   -k;
            *V           +=    k*pcrd->value;
            *dVdl        += dkdl*pcrd->value;
            break;
        case epullEXTERNAL:
            gmx_incons("the scalar pull force should not be calculated internally for pull type external");
            break;
        default:
            gmx_incons("Unsupported pull type in do_pull_pot");
    }
}

static void calc_pull_coord_vector_force(pull_coord_work_t *pcrd)
{
    /* The geometry of the coordinate determines how the scalar force relates to the force on each group */

    if (pcrd->params.eGeom == epullgDIST)
    {
        double invdr01 = pcrd->value > 0 ? 1./pcrd->value : 0.;
        for (int m = 0; m < DIM; m++)
        {
            pcrd->f01[m] = pcrd->f_scal*pcrd->dr01[m]*invdr01;
        }
    }
    else if (pcrd->params.eGeom == epullgANGLE)
    {

        double cos_theta, cos_theta2;

        cos_theta  = cos(pcrd->value);
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
            double b       = a*cos_theta;
            double invdr01 = 1./dnorm(pcrd->dr01);
            double invdr23 = 1./dnorm(pcrd->dr23);
            dvec   normalized_dr01, normalized_dr23;
            dsvmul(invdr01, pcrd->dr01, normalized_dr01);
            dsvmul(invdr23, pcrd->dr23, normalized_dr23);

            for (int m = 0; m < DIM; m++)
            {
                /* Here, f_scal is -dV/dtheta */
                pcrd->f01[m] = pcrd->f_scal*invdr01*(a*normalized_dr23[m] - b*normalized_dr01[m]);
                pcrd->f23[m] = pcrd->f_scal*invdr23*(a*normalized_dr01[m] - b*normalized_dr23[m]);
            }
        }
        else
        {
            /* No forces to apply for ill-defined cases*/
            clear_pull_forces_coord(pcrd);
        }
    }
    else if (pcrd->params.eGeom == epullgANGLEAXIS)
    {
        double cos_theta, cos_theta2;

        /* The angle-axis force is exactly the same as the angle force (above) except that in
           this case the second vector (dr23) is replaced by the pull vector. */
        cos_theta  = cos(pcrd->value);
        cos_theta2 = gmx::square(cos_theta);

        if (cos_theta2 < 1)
        {
            double a, b;
            double invdr01;
            dvec   normalized_dr01;

            invdr01 = 1./dnorm(pcrd->dr01);
            dsvmul(invdr01, pcrd->dr01, normalized_dr01);
            a       = -gmx::invsqrt(1 - cos_theta2); /* comes from d/dx acos(x) */
            b       = a*cos_theta;

            for (int m = 0; m < DIM; m++)
            {
                pcrd->f01[m] = pcrd->f_scal*invdr01*(a*pcrd->vec[m] - b*normalized_dr01[m]);
            }
        }
        else
        {
            clear_pull_forces_coord(pcrd);
        }
    }
    else if (pcrd->params.eGeom == epullgDIHEDRAL)
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
        m2 = diprod(pcrd->planevec_m, pcrd->planevec_m);
        n2 = diprod(pcrd->planevec_n, pcrd->planevec_n);
        dsvmul(-1, pcrd->dr23, dr32);
        sqrdist_32 = diprod(dr32, dr32);
        tol        = sqrdist_32*GMX_REAL_EPS; /* Avoid tiny angles */
        if ((m2 > tol) && (n2 > tol))
        {
            double a_01, a_23_01, a_23_45, a_45;
            double inv_dist_32, inv_sqrdist_32, dist_32;
            dvec   u, v;
            inv_dist_32    = gmx::invsqrt(sqrdist_32);
            inv_sqrdist_32 = inv_dist_32*inv_dist_32;
            dist_32        = sqrdist_32*inv_dist_32;

            /* Forces on groups 0, 1 */
            a_01 = pcrd->f_scal*dist_32/m2;             /* f_scal is -dV/dphi */
            dsvmul(-a_01, pcrd->planevec_m, pcrd->f01); /* added sign to get force on group 1, not 0 */

            /* Forces on groups 4, 5 */
            a_45 = -pcrd->f_scal*dist_32/n2;
            dsvmul(a_45, pcrd->planevec_n, pcrd->f45); /* force on group 5 */

            /* Force on groups 2, 3 (defining the axis) */
            a_23_01 = -diprod(pcrd->dr01, dr32)*inv_sqrdist_32;
            a_23_45 = -diprod(pcrd->dr45, dr32)*inv_sqrdist_32;
            dsvmul(-a_23_01, pcrd->f01, u);                     /* added sign to get force from group 0, not 1 */
            dsvmul(a_23_45, pcrd->f45, v);
            dvec_sub(u, v, pcrd->f23);                          /* force on group 3 */
        }
        else
        {
            /* No force to apply for ill-defined cases */
            clear_pull_forces_coord(pcrd);
        }
    }
    else
    {
        for (int m = 0; m < DIM; m++)
        {
            pcrd->f01[m] = pcrd->f_scal*pcrd->vec[m];
        }
    }
}


/* We use a global mutex for locking access to the pull data structure
 * during registration of external pull potential providers.
 * We could use a different, local mutex for each pull object, but the overhead
 * is extremely small here and registration is only done during initialization.
 */
static gmx::Mutex registrationMutex;

using Lock = gmx::lock_guard<gmx::Mutex>;

void register_external_pull_potential(struct pull_t *pull,
                                      int            coord_index,
                                      const char    *provider)
{
    GMX_RELEASE_ASSERT(pull != nullptr, "register_external_pull_potential called before init_pull");
    GMX_RELEASE_ASSERT(provider != nullptr, "register_external_pull_potential called with NULL as provider name");

    if (coord_index < 0 || coord_index > pull->ncoord - 1)
    {
        gmx_fatal(FARGS, "Module '%s' attempted to register an external potential for pull coordinate %d which is out of the pull coordinate range %d - %d\n",
                  provider, coord_index + 1, 1, pull->ncoord);
    }

    pull_coord_work_t *pcrd = &pull->coord[coord_index];

    if (pcrd->params.eType != epullEXTERNAL)
    {
        gmx_fatal(FARGS, "Module '%s' attempted to register an external potential for pull coordinate %d which of type '%s', whereas external potentials are only supported with type '%s'",
                  provider, coord_index + 1, epull_names[pcrd->params.eType], epull_names[epullEXTERNAL]);
    }

    GMX_RELEASE_ASSERT(pcrd->params.externalPotentialProvider != nullptr, "The external potential provider string for a pull coordinate is NULL");

    if (gmx_strcasecmp(provider, pcrd->params.externalPotentialProvider) != 0)
    {
        gmx_fatal(FARGS, "Module '%s' attempted to register an external potential for pull coordinate %d which expects the external potential to be provided by a module named '%s'",
                  provider, coord_index + 1, pcrd->params.externalPotentialProvider);
    }

    /* Lock to avoid (extremely unlikely) simultaneous reading and writing of
     * pcrd->bExternalPotentialProviderHasBeenRegistered and
     * pull->numUnregisteredExternalPotentials.
     */
    Lock registrationLock(registrationMutex);

    if (pcrd->bExternalPotentialProviderHasBeenRegistered)
    {
        gmx_fatal(FARGS, "Module '%s' attempted to register an external potential for pull coordinate %d more than once",
                  provider, coord_index + 1);
    }

    pcrd->bExternalPotentialProviderHasBeenRegistered = true;
    pull->numUnregisteredExternalPotentials--;

    GMX_RELEASE_ASSERT(pull->numUnregisteredExternalPotentials >= 0, "Negative unregistered potentials, the pull code is inconsistent");
}


static void check_external_potential_registration(const struct pull_t *pull)
{
    if (pull->numUnregisteredExternalPotentials > 0)
    {
        int c;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (pull->coord[c].params.eType == epullEXTERNAL &&
                !pull->coord[c].bExternalPotentialProviderHasBeenRegistered)
            {
                break;
            }
        }

        GMX_RELEASE_ASSERT(c < pull->ncoord, "Internal inconsistency in the pull potential provider counting");

        gmx_fatal(FARGS, "No external provider for external pull potentials have been provided for %d pull coordinates. The first coordinate without provider is number %d, which expects a module named '%s' to provide the external potential.",
                  pull->numUnregisteredExternalPotentials,
                  c + 1, pull->coord[c].params.externalPotentialProvider);
    }
}


/* Pull takes care of adding the forces of the external potential.
 * The external potential module  has to make sure that the corresponding
 * potential energy is added either to the pull term or to a term
 * specific to the external module.
 */
void apply_external_pull_coord_force(struct pull_t        *pull,
                                     int                   coord_index,
                                     double                coord_force,
                                     const t_mdatoms      *mdatoms,
                                     gmx::ForceWithVirial *forceWithVirial)
{
    pull_coord_work_t *pcrd;

    GMX_ASSERT(coord_index >= 0 && coord_index < pull->ncoord, "apply_external_pull_coord_force called with coord_index out of range");

    if (pull->comm.bParticipate)
    {
        pcrd = &pull->coord[coord_index];

        GMX_RELEASE_ASSERT(pcrd->params.eType == epullEXTERNAL, "The pull force can only be set externally on pull coordinates of external type");

        GMX_ASSERT(pcrd->bExternalPotentialProviderHasBeenRegistered, "apply_external_pull_coord_force called for an unregistered pull coordinate");

        /* Set the force */
        pcrd->f_scal = coord_force;

        /* Calculate the forces on the pull groups */
        calc_pull_coord_vector_force(pcrd);

        /* Add the forces for this coordinate to the total virial and force */
        if (forceWithVirial->computeVirial_)
        {
            matrix virial = { { 0 } };
            add_virial_coord(virial, pcrd);
            forceWithVirial->addVirialContribution(virial);
        }

        apply_forces_coord(pull, coord_index, mdatoms,
                           as_rvec_array(forceWithVirial->force_.data()));
    }

    pull->numExternalPotentialsStillToBeAppliedThisStep--;
}

/* Calculate the pull potential and scalar force for a pull coordinate */
static void do_pull_pot_coord(struct pull_t *pull, int coord_ind, t_pbc *pbc,
                              double t, real lambda,
                              real *V, tensor vir, real *dVdl)
{
    pull_coord_work_t *pcrd;
    double             dev;

    pcrd = &pull->coord[coord_ind];

    assert(pcrd->params.eType != epullCONSTRAINT);

    dev = get_pull_coord_deviation(pull, coord_ind, pbc, t);

    calc_pull_coord_scalar_force_and_potential(pcrd, dev, lambda, V, dVdl);

    calc_pull_coord_vector_force(pcrd);

    add_virial_coord(vir, pcrd);
}

real pull_potential(struct pull_t *pull, t_mdatoms *md, t_pbc *pbc,
                    t_commrec *cr, double t, real lambda,
                    rvec *x, gmx::ForceWithVirial *force, real *dvdlambda)
{
    real V = 0;

    assert(pull != NULL);

    /* Ideally we should check external potential registration only during
     * the initialization phase, but that requires another function call
     * that should be done exactly in the right order. So we check here.
     */
    check_external_potential_registration(pull);

    if (pull->comm.bParticipate)
    {
        real dVdl = 0;

        pull_calc_coms(cr, pull, md, pbc, t, x, nullptr);

        rvec       *f             = as_rvec_array(force->force_.data());
        matrix      virial        = { { 0 } };
        const bool  computeVirial = (force->computeVirial_ && MASTER(cr));
        for (int c = 0; c < pull->ncoord; c++)
        {
            /* For external potential the force is assumed to be given by an external module by a call to
               apply_pull_coord_external_force */
            if (pull->coord[c].params.eType == epullCONSTRAINT || pull->coord[c].params.eType == epullEXTERNAL)
            {
                continue;
            }

            do_pull_pot_coord(pull, c, pbc, t, lambda,
                              &V,
                              computeVirial ? virial : nullptr,
                              &dVdl);

            /* Distribute the force over the atoms in the pulled groups */
            apply_forces_coord(pull, c, md, f);
        }

        if (MASTER(cr))
        {
            force->addVirialContribution(virial);
            *dvdlambda += dVdl;
        }
    }

    GMX_ASSERT(pull->numExternalPotentialsStillToBeAppliedThisStep == 0, "Too few or too many external pull potentials have been applied the previous step");
    /* All external pull potentials still need to be applied */
    pull->numExternalPotentialsStillToBeAppliedThisStep =
        pull->numCoordinatesWithExternalPotential;

    return (MASTER(cr) ? V : 0.0);
}

void pull_constraint(struct pull_t *pull, t_mdatoms *md, t_pbc *pbc,
                     t_commrec *cr, double dt, double t,
                     rvec *x, rvec *xp, rvec *v, tensor vir)
{
    assert(pull != NULL);

    if (pull->comm.bParticipate)
    {
        pull_calc_coms(cr, pull, md, pbc, t, x, xp);

        do_constraint(pull, pbc, xp, v, MASTER(cr), vir, dt, t);
    }
}

static void make_local_pull_group(gmx_ga2la_t *ga2la,
                                  pull_group_work_t *pg, int start, int end)
{
    int i, ii;

    pg->nat_loc = 0;
    for (i = 0; i < pg->params.nat; i++)
    {
        ii = pg->params.ind[i];
        if (ga2la)
        {
            if (!ga2la_get_home(ga2la, ii, &ii))
            {
                ii = -1;
            }
        }
        if (ii >= start && ii < end)
        {
            /* This is a home atom, add it to the local pull group */
            if (pg->nat_loc >= pg->nalloc_loc)
            {
                pg->nalloc_loc = over_alloc_dd(pg->nat_loc+1);
                srenew(pg->ind_loc, pg->nalloc_loc);
                if (pg->epgrppbc == epgrppbcCOS || pg->params.weight != nullptr)
                {
                    srenew(pg->weight_loc, pg->nalloc_loc);
                }
            }
            pg->ind_loc[pg->nat_loc] = ii;
            if (pg->params.weight != nullptr)
            {
                pg->weight_loc[pg->nat_loc] = pg->params.weight[i];
            }
            pg->nat_loc++;
        }
    }
}

void dd_make_local_pull_groups(t_commrec *cr, struct pull_t *pull, t_mdatoms *md)
{
    gmx_domdec_t   *dd;
    pull_comm_t    *comm;
    gmx_ga2la_t    *ga2la;
    gmx_bool        bMustParticipate;
    int             g;

    dd = cr->dd;

    comm = &pull->comm;

    if (dd)
    {
        ga2la = dd->ga2la;
    }
    else
    {
        ga2la = nullptr;
    }

    /* We always make the master node participate, such that it can do i/o
     * and to simplify MC type extensions people might have.
     */
    bMustParticipate = (comm->bParticipateAll || dd == nullptr || DDMASTER(dd));

    for (g = 0; g < pull->ngroup; g++)
    {
        int a;

        make_local_pull_group(ga2la, &pull->group[g],
                              0, md->homenr);

        /* We should participate if we have pull or pbc atoms */
        if (!bMustParticipate &&
            (pull->group[g].nat_loc > 0 ||
             (pull->group[g].params.pbcatom >= 0 &&
              ga2la_get_home(dd->ga2la, pull->group[g].params.pbcatom, &a))))
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
        const gmx_int64_t history_count = 20;
        gmx_bool          bWillParticipate;
        int               count[2];

        /* Increase the decomposition counter for the current call */
        comm->setup_count++;

        if (bMustParticipate)
        {
            comm->must_count = comm->setup_count;
        }

        bWillParticipate =
            bMustParticipate ||
            (comm->bParticipate &&
             comm->must_count >= comm->setup_count - history_count);

        if (debug && dd != nullptr)
        {
            fprintf(debug, "Our DD rank (%3d) pull #atoms>0 or master: %d, will be part %d\n",
                    dd->rank, bMustParticipate, bWillParticipate);
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
        if (count[1] > 0 || 6*count[0] < 5*comm->nparticipate)
        {
            if (debug)
            {
                fprintf(debug, "Creating new pull subcommunicator of size %d\n",
                        count[0]);
            }

#if GMX_MPI
            if (comm->mpi_comm_com != MPI_COMM_NULL)
            {
                MPI_Comm_free(&comm->mpi_comm_com);
            }

            /* This might be an extremely expensive operation, so we try
             * to avoid this splitting as much as possible.
             */
            assert(dd != NULL);
            MPI_Comm_split(dd->mpi_comm_all, bWillParticipate ? 0 : 1, dd->rank,
                           &comm->mpi_comm_com);
#endif

            comm->bParticipate = bWillParticipate;
            comm->nparticipate = count[0];
        }
    }

    /* Since the PBC of atoms might have changed, we need to update the PBC */
    pull->bSetPBCatoms = TRUE;
}

static void init_pull_group_index(FILE *fplog, t_commrec *cr,
                                  int g, pull_group_work_t *pg,
                                  gmx_bool bConstraint, ivec pulldim_con,
                                  const gmx_mtop_t *mtop,
                                  const t_inputrec *ir, real lambda)
{
    if (EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
    {
        /* There are no masses in the integrator.
         * But we still want to have the correct mass-weighted COMs.
         * So we store the real masses in the weights.
         * We do not set nweight, so these weights do not end up in the tpx file.
         */
        if (pg->params.nweight == 0)
        {
            snew(pg->params.weight, pg->params.nat);
        }
    }

    if (cr && PAR(cr))
    {
        pg->nat_loc    = 0;
        pg->nalloc_loc = 0;
        pg->ind_loc    = nullptr;
        pg->weight_loc = nullptr;
    }
    else
    {
        pg->nat_loc = pg->params.nat;
        pg->ind_loc = pg->params.ind;
        if (pg->epgrppbc == epgrppbcCOS)
        {
            snew(pg->weight_loc, pg->params.nat);
        }
        else
        {
            pg->weight_loc = pg->params.weight;
        }
    }

    const gmx_groups_t *groups = &mtop->groups;

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
                if (pulldim_con[d] == 1 &&
                    ir->opts.nFreeze[ggrpnr(groups, egcFREEZE, ii)][d])
                {
                    nfrozen++;
                }
            }
        }
        const t_atom &atom = mtopGetAtomParameters(mtop, ii, &molb);
        real          m;
        if (ir->efep == efepNO)
        {
            m = atom.m;
        }
        else
        {
            m = (1 - lambda)*atom.m + lambda*atom.mB;
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
            w                   *= m;
            m                    = 1;
            pg->params.weight[i] = w;
        }
        else if (ir->eI == eiBD)
        {
            real mbd;
            if (ir->bd_fric)
            {
                mbd = ir->bd_fric*ir->delta_t;
            }
            else
            {
                if (groups->grpnr[egcTC] == nullptr)
                {
                    mbd = ir->delta_t/ir->opts.tau_t[0];
                }
                else
                {
                    mbd = ir->delta_t/ir->opts.tau_t[groups->grpnr[egcTC][ii]];
                }
            }
            w                   *= m/mbd;
            m                    = mbd;
            pg->params.weight[i] = w;
        }
        tmass  += m;
        wmass  += m*w;
        wwmass += m*w*w;
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
        fprintf(fplog,
                "Pull group %d: %5d atoms, mass %9.3f",
                g, pg->params.nat, tmass);
        if (pg->params.weight ||
            EI_ENERGY_MINIMIZATION(ir->eI) || ir->eI == eiBD)
        {
            fprintf(fplog, ", weighted mass %9.3f", wmass*wmass/wwmass);
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
        pg->invtm  = -1.0;
    }
    else
    {
        int ndim = 0;
        for (int d = 0; d < DIM; d++)
        {
            ndim += pulldim_con[d]*pg->params.nat;
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

struct pull_t *
init_pull(FILE *fplog, const pull_params_t *pull_params, const t_inputrec *ir,
          int nfile, const t_filenm fnm[],
          const gmx_mtop_t *mtop, t_commrec *cr,
          const gmx_output_env_t *oenv, real lambda,
          gmx_bool bOutFile,
          const ContinuationOptions &continuationOptions)
{
    struct pull_t *pull;
    pull_comm_t   *comm;
    int            g, c, m;

    snew(pull, 1);

    /* Copy the pull parameters */
    pull->params       = *pull_params;
    /* Avoid pointer copies */
    pull->params.group = nullptr;
    pull->params.coord = nullptr;

    pull->ncoord       = pull_params->ncoord;
    pull->ngroup       = pull_params->ngroup;
    snew(pull->coord, pull->ncoord);
    snew(pull->group, pull->ngroup);

    pull->bPotential  = FALSE;
    pull->bConstraint = FALSE;
    pull->bCylinder   = FALSE;
    pull->bAngle      = FALSE;

    for (g = 0; g < pull->ngroup; g++)
    {
        pull_group_work_t *pgrp;
        int                i;

        pgrp = &pull->group[g];

        /* Copy the pull group parameters */
        pgrp->params = pull_params->group[g];

        /* Avoid pointer copies by allocating and copying arrays */
        snew(pgrp->params.ind, pgrp->params.nat);
        for (i = 0; i < pgrp->params.nat; i++)
        {
            pgrp->params.ind[i] = pull_params->group[g].ind[i];
        }
        if (pgrp->params.nweight > 0)
        {
            snew(pgrp->params.weight, pgrp->params.nweight);
            for (i = 0; i < pgrp->params.nweight; i++)
            {
                pgrp->params.weight[i] = pull_params->group[g].weight[i];
            }
        }
    }

    pull->numCoordinatesWithExternalPotential = 0;

    for (c = 0; c < pull->ncoord; c++)
    {
        pull_coord_work_t *pcrd;
        int                calc_com_start, calc_com_end, g;

        pcrd = &pull->coord[c];

        /* Copy all pull coordinate parameters */
        pcrd->params = pull_params->coord[c];

        switch (pcrd->params.eGeom)
        {
            case epullgDIST:
            case epullgDIRRELATIVE:  /* Direction vector is determined at each step */
            case epullgANGLE:
            case epullgDIHEDRAL:
                break;
            case epullgDIR:
            case epullgDIRPBC:
            case epullgCYL:
            case epullgANGLEAXIS:
                copy_rvec_to_dvec(pull_params->coord[c].vec, pcrd->vec);
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
                gmx_fatal(FARGS, "Pull geometry not supported for pull coordinate %d. The geometry enum %d in the input is larger than that supported by the code (up to %d). You are probably reading a tpr file generated with a newer version of Gromacs with an binary from an older version of Gromacs.",
                          c + 1, pcrd->params.eGeom, epullgNR - 1);
        }

        if (pcrd->params.eType == epullCONSTRAINT)
        {
            /* Check restrictions of the constraint pull code */
            if (pcrd->params.eGeom == epullgCYL ||
                pcrd->params.eGeom == epullgDIRRELATIVE ||
                pcrd->params.eGeom == epullgANGLE ||
                pcrd->params.eGeom == epullgDIHEDRAL ||
                pcrd->params.eGeom == epullgANGLEAXIS)
            {
                gmx_fatal(FARGS, "Pulling of type %s can not be combined with geometry %s. Consider using pull type %s.",
                          epull_names[pcrd->params.eType],
                          epullg_names[pcrd->params.eGeom],
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
        else if (pcrd->params.eGeom == epullgANGLE || pcrd->params.eGeom == epullgDIHEDRAL || pcrd->params.eGeom == epullgANGLEAXIS)
        {
            pull->bAngle = TRUE;
        }

        /* We only need to calculate the plain COM of a group
         * when it is not only used as a cylinder group.
         */
        calc_com_start = (pcrd->params.eGeom == epullgCYL         ? 1 : 0);
        calc_com_end   = pcrd->params.ngroup;

        for (g = calc_com_start; g <= calc_com_end; g++)
        {
            pull->group[pcrd->params.group[g]].bCalcCOM = TRUE;
        }

        /* With non-zero rate the reference value is set at every step */
        if (pcrd->params.rate == 0)
        {
            /* Initialize the constant reference value */
            if (pcrd->params.eType != epullEXTERNAL)
            {
                low_set_pull_coord_reference_value(pcrd, c, pcrd->params.init*pull_conversion_factor_userinput2internal(&pcrd->params));
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
            GMX_RELEASE_ASSERT(pcrd->params.rate == 0, "With an external potential, a pull coordinate should have rate = 0");

            /* This external potential needs to be registered later */
            pull->numCoordinatesWithExternalPotential++;
        }
        pcrd->bExternalPotentialProviderHasBeenRegistered = false;
    }

    pull->numUnregisteredExternalPotentials             =
        pull->numCoordinatesWithExternalPotential;
    pull->numExternalPotentialsStillToBeAppliedThisStep = 0;

    pull->ePBC = ir->ePBC;
    switch (pull->ePBC)
    {
        case epbcNONE: pull->npbcdim = 0; break;
        case epbcXY:   pull->npbcdim = 2; break;
        default:       pull->npbcdim = 3; break;
    }

    if (fplog)
    {
        gmx_bool bAbs, bCos;

        bAbs = FALSE;
        for (c = 0; c < pull->ncoord; c++)
        {
            if (pull->group[pull->coord[c].params.group[0]].params.nat == 0 ||
                pull->group[pull->coord[c].params.group[1]].params.nat == 0)
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
        fprintf(fplog, "with %d pull coordinate%s and %d group%s\n",
                pull->ncoord, pull->ncoord == 1 ? "" : "s",
                pull->ngroup, pull->ngroup == 1 ? "" : "s");
        if (bAbs)
        {
            fprintf(fplog, "with an absolute reference\n");
        }
        bCos = FALSE;
        for (g = 0; g < pull->ngroup; g++)
        {
            if (pull->group[g].params.nat > 1 &&
                pull->group[g].params.pbcatom < 0)
            {
                /* We are using cosine weighting */
                fprintf(fplog, "Cosine weighting is used for group %d\n", g);
                bCos = TRUE;
            }
        }
        if (bCos)
        {
            please_cite(fplog, "Engin2010");
        }
    }

    pull->bRefAt   = FALSE;
    pull->cosdim   = -1;
    for (g = 0; g < pull->ngroup; g++)
    {
        pull_group_work_t *pgrp;

        pgrp           = &pull->group[g];
        pgrp->epgrppbc = epgrppbcNONE;
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
            for (c = 0; c < pull->ncoord; c++)
            {
                pull_coord_work_t *pcrd;
                int                gi;
                gmx_bool           bGroupUsed;

                pcrd = &pull->coord[c];

                bGroupUsed = FALSE;
                for (gi = 0; gi < pcrd->params.ngroup; gi++)
                {
                    if (pcrd->params.group[gi] == g)
                    {
                        bGroupUsed = TRUE;
                    }
                }

                if (bGroupUsed)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (pcrd->params.dim[m] == 1)
                        {
                            pulldim[m] = 1;

                            if (pcrd->params.eType == epullCONSTRAINT)
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
            GMX_RELEASE_ASSERT(pull->npbcdim <= DIM, "npbcdim must be <= the number of dimensions");
            for (m = 0; m < pull->npbcdim; m++)
            {
                GMX_RELEASE_ASSERT(m <= DIM, "npbcdim must be <= the number of dimensions");
                if (pulldim[m] == 1 && pgrp->params.nat > 1)
                {
                    if (pgrp->params.pbcatom >= 0)
                    {
                        pgrp->epgrppbc = epgrppbcREFAT;
                        pull->bRefAt   = TRUE;
                    }
                    else
                    {
                        if (pgrp->params.weight != nullptr)
                        {
                            gmx_fatal(FARGS, "Pull groups can not have relative weights and cosine weighting at same time");
                        }
                        pgrp->epgrppbc = epgrppbcCOS;
                        if (pull->cosdim >= 0 && pull->cosdim != m)
                        {
                            gmx_fatal(FARGS, "Can only use cosine weighting with pulling in one dimension (use mdp option pull-coord?-dim)");
                        }
                        pull->cosdim = m;
                    }
                }
            }
            /* Set the indices */
            init_pull_group_index(fplog, cr, g, pgrp,
                                  bConstraint, pulldim_con,
                                  mtop, ir, lambda);
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
        snew(pull->dyna, pull->ncoord);

        for (c = 0; c < pull->ncoord; c++)
        {
            const pull_coord_work_t *pcrd;

            pcrd = &pull->coord[c];

            if (pcrd->params.eGeom == epullgCYL)
            {
                if (pull->group[pcrd->params.group[0]].params.nat == 0)
                {
                    gmx_fatal(FARGS, "A cylinder pull group is not supported when using absolute reference!\n");
                }
            }
        }
    }

    /* The gmx_omp_nthreads module might not be initialized here, so max(1,) */
    pull->nthreads = std::max(1, gmx_omp_nthreads_get(emntDefault));
    snew(pull->sum_com, pull->nthreads);

    comm = &pull->comm;

#if GMX_MPI
    /* Use a sub-communicator when we have more than 32 ranks */
    comm->bParticipateAll = (cr == nullptr || !DOMAINDECOMP(cr) ||
                             cr->dd->nnodes <= 32 ||
                             getenv("GMX_PULL_PARTICIPATE_ALL") != nullptr);
    /* This sub-commicator is not used with comm->bParticipateAll,
     * so we can always initialize it to NULL.
     */
    comm->mpi_comm_com    = MPI_COMM_NULL;
    comm->nparticipate    = 0;
#else
    /* No MPI: 1 rank: all ranks pull */
    comm->bParticipateAll = TRUE;
#endif
    comm->bParticipate    = comm->bParticipateAll;
    comm->setup_count     = 0;
    comm->must_count      = 0;

    if (!comm->bParticipateAll && fplog != nullptr)
    {
        fprintf(fplog, "Will use a sub-communicator for pull communication\n");
    }

    comm->rbuf     = nullptr;
    comm->dbuf     = nullptr;
    comm->dbuf_cyl = nullptr;

    /* We still need to initialize the PBC reference coordinates */
    pull->bSetPBCatoms = TRUE;

    /* Only do I/O when we are doing dynamics and if we are the MASTER */
    pull->out_x = nullptr;
    pull->out_f = nullptr;
    if (bOutFile)
    {
        /* Check for px and pf filename collision, if we are writing
           both files */
        std::string px_filename, pf_filename;
        std::string px_appended, pf_appended;
        try
        {
            px_filename  = std::string(opt2fn("-px", nfile, fnm));
            pf_filename  = std::string(opt2fn("-pf", nfile, fnm));
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;


        if ((pull->params.nstxout != 0) &&
            (pull->params.nstfout != 0) &&
            (px_filename == pf_filename))
        {
            if (!opt2bSet("-px", nfile, fnm) && !opt2bSet("-pf", nfile, fnm))
            {
                /* We are writing both pull files but neither set directly. */
                try
                {
                    px_appended   = append_before_extension(px_filename, "_pullx");
                    pf_appended   = append_before_extension(pf_filename, "_pullf");
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                pull->out_x   = open_pull_out(px_appended.c_str(), pull, oenv,
                                              TRUE, continuationOptions);
                pull->out_f = open_pull_out(pf_appended.c_str(), pull, oenv,
                                            FALSE, continuationOptions);
                return pull;
            }
            else
            {
                /* If one of -px and -pf is set but the filenames are identical: */
                gmx_fatal(FARGS, "Identical pull_x and pull_f output filenames %s",
                          px_filename.c_str());
            }
        }
        if (pull->params.nstxout != 0)
        {
            pull->out_x = open_pull_out(opt2fn("-px", nfile, fnm), pull, oenv,
                                        TRUE, continuationOptions);
        }
        if (pull->params.nstfout != 0)
        {
            pull->out_f = open_pull_out(opt2fn("-pf", nfile, fnm), pull, oenv,
                                        FALSE, continuationOptions);
        }
    }

    return pull;
}

static void destroy_pull_group(pull_group_work_t *pgrp)
{
    if (pgrp->ind_loc != pgrp->params.ind)
    {
        sfree(pgrp->ind_loc);
    }
    if (pgrp->weight_loc != pgrp->params.weight)
    {
        sfree(pgrp->weight_loc);
    }
    sfree(pgrp->mdw);
    sfree(pgrp->dv);

    sfree(pgrp->params.ind);
    sfree(pgrp->params.weight);
}

static void destroy_pull(struct pull_t *pull)
{
    int i;

    for (i = 0; i < pull->ngroup; i++)
    {
        destroy_pull_group(&pull->group[i]);
    }
    for (i = 0; i < pull->ncoord; i++)
    {
        if (pull->coord[i].params.eGeom == epullgCYL)
        {
            destroy_pull_group(&(pull->dyna[i]));
        }
    }
    sfree(pull->group);
    sfree(pull->dyna);
    sfree(pull->coord);

#if GMX_MPI
    if (pull->comm.mpi_comm_com != MPI_COMM_NULL)
    {
        MPI_Comm_free(&pull->comm.mpi_comm_com);
    }
#endif
    sfree(pull->comm.rbuf);
    sfree(pull->comm.dbuf);
    sfree(pull->comm.dbuf_cyl);

    sfree(pull);
}

void finish_pull(struct pull_t *pull)
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

gmx_bool pull_have_potential(const struct pull_t *pull)
{
    return pull->bPotential;
}

gmx_bool pull_have_constraint(const struct pull_t *pull)
{
    return pull->bConstraint;
}
