/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"

#include <sys/types.h>
#include <math.h>
#include "typedefs.h"
#include "physics.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "vec.h"
#include "coulomb.h"
#include "calc_verletbuf.h"
#include "../mdlib/nbnxn_consts.h"

#ifdef GMX_NBNXN_SIMD
/* The include below sets the SIMD instruction type (precision+width)
 * for all nbnxn SIMD search and non-bonded kernel code.
 */
#ifdef GMX_NBNXN_HALF_WIDTH_SIMD
#define GMX_USE_HALF_WIDTH_SIMD_HERE
#endif
#include "gmx_simd_macros.h"
#endif


/* The code in this file estimates a pairlist buffer length
 * given a target energy drift per atom per picosecond.
 * This is done by estimating the drift given a buffer length.
 * Ideally we would like to have a tight overestimate of the drift,
 * but that can be difficult to achieve.
 *
 * Significant approximations used:
 *
 * Uniform particle density. UNDERESTIMATES the drift by rho_global/rho_local.
 *
 * Interactions don't affect particle motion. OVERESTIMATES the drift on longer
 * time scales. This approximation probably introduces the largest errors.
 *
 * Only take one constraint per particle into account: OVERESTIMATES the drift.
 *
 * For rotating constraints assume the same functional shape for time scales
 * where the constraints rotate significantly as the exact expression for
 * short time scales. OVERESTIMATES the drift on long time scales.
 *
 * For non-linear virtual sites use the mass of the lightest constructing atom
 * to determine the displacement. OVER/UNDERESTIMATES the drift, depending on
 * the geometry and masses of constructing atoms.
 *
 * Note that the formulas for normal atoms and linear virtual sites are exact,
 * apart from the first two approximations.
 *
 * Note that apart from the effect of the above approximations, the actual
 * drift of the total energy of a system can be order of magnitude smaller
 * due to cancellation of positive and negative drift for different pairs.
 */


/* Struct for unique atom type for calculating the energy drift.
 * The atom displacement depends on mass and constraints.
 * The energy jump for given distance depend on LJ type and q.
 */
typedef struct
{
    real     mass;     /* mass */
    int      type;     /* type (used for LJ parameters) */
    real     q;        /* charge */
    gmx_bool bConstr;  /* constrained, if TRUE, use #DOF=2 iso 3 */
    real     con_mass; /* mass of heaviest atom connected by constraints */
    real     con_len;  /* constraint length to the heaviest atom */
} atom_nonbonded_kinetic_prop_t;

/* Struct for unique atom type for calculating the energy drift.
 * The atom displacement depends on mass and constraints.
 * The energy jump for given distance depend on LJ type and q.
 */
typedef struct
{
    atom_nonbonded_kinetic_prop_t prop; /* non-bonded and kinetic atom prop. */
    int                           n;    /* #atoms of this type in the system */
} verletbuf_atomtype_t;

void verletbuf_get_list_setup(gmx_bool                bGPU,
                              verletbuf_list_setup_t *list_setup)
{
    list_setup->cluster_size_i     = NBNXN_CPU_CLUSTER_I_SIZE;

    if (bGPU)
    {
        list_setup->cluster_size_j = NBNXN_GPU_CLUSTER_SIZE;
    }
    else
    {
#ifndef GMX_NBNXN_SIMD
        list_setup->cluster_size_j = NBNXN_CPU_CLUSTER_I_SIZE;
#else
        list_setup->cluster_size_j = GMX_SIMD_WIDTH_HERE;
#ifdef GMX_NBNXN_SIMD_2XNN
        /* We assume the smallest cluster size to be on the safe side */
        list_setup->cluster_size_j /= 2;
#endif
#endif
    }
}

static gmx_bool
atom_nonbonded_kinetic_prop_equal(const atom_nonbonded_kinetic_prop_t *prop1,
                                  const atom_nonbonded_kinetic_prop_t *prop2)
{
    return (prop1->mass     == prop2->mass &&
            prop1->type     == prop2->type &&
            prop1->q        == prop2->q &&
            prop1->bConstr  == prop2->bConstr &&
            prop1->con_mass == prop2->con_mass &&
            prop1->con_len  == prop2->con_len);
}

static void add_at(verletbuf_atomtype_t **att_p, int *natt_p,
                   const atom_nonbonded_kinetic_prop_t *prop,
                   int nmol)
{
    verletbuf_atomtype_t *att;
    int                     natt, i;

    if (prop->mass == 0)
    {
        /* Ignore massless particles */
        return;
    }

    att  = *att_p;
    natt = *natt_p;

    i = 0;
    while (i < natt && !atom_nonbonded_kinetic_prop_equal(prop, &att[i].prop))
    {
        i++;
    }

    if (i < natt)
    {
        att[i].n += nmol;
    }
    else
    {
        (*natt_p)++;
        srenew(*att_p, *natt_p);
        (*att_p)[i].prop = *prop;
        (*att_p)[i].n    = nmol;
    }
}

static void get_vsite_masses(const gmx_moltype_t *moltype,
                             const gmx_ffparams_t *ffparams,
                             real *vsite_m,
                             int *n_nonlin_vsite)
{
    int            ft, i;
    const t_ilist *il;

    *n_nonlin_vsite = 0;

    /* Check for virtual sites, determine mass from constructing atoms */
    for (ft = 0; ft < F_NRE; ft++)
    {
        if (IS_VSITE(ft))
        {
            il = &moltype->ilist[ft];

            for (i = 0; i < il->nr; i += 1+NRAL(ft))
            {
                const t_iparams *ip;
                real             cam[5], inv_mass, m_aj;
                int              a1, j, aj, coeff;

                ip = &ffparams->iparams[il->iatoms[i]];

                a1 = il->iatoms[i+1];

                if (ft != F_VSITEN)
                {
                    for (j = 1; j < NRAL(ft); j++)
                    {
                        cam[j] = moltype->atoms.atom[il->iatoms[i+1+j]].m;
                        if (cam[j] == 0)
                        {
                            cam[j] = vsite_m[il->iatoms[i+1+j]];
                        }
                        if (cam[j] == 0)
                        {
                            gmx_fatal(FARGS, "In molecule type '%s' %s construction involves atom %d, which is a virtual site of equal or high complexity. This is not supported.",
                                      *moltype->name,
                                      interaction_function[ft].longname,
                                      il->iatoms[i+1+j]+1);
                        }
                    }
                }

                switch (ft)
                {
                    case F_VSITE2:
                        /* Exact */
                        vsite_m[a1] = (cam[1]*cam[2])/(cam[2]*sqr(1-ip->vsite.a) + cam[1]*sqr(ip->vsite.a));
                        break;
                    case F_VSITE3:
                        /* Exact */
                        vsite_m[a1] = (cam[1]*cam[2]*cam[3])/(cam[2]*cam[3]*sqr(1-ip->vsite.a-ip->vsite.b) + cam[1]*cam[3]*sqr(ip->vsite.a) + cam[1]*cam[2]*sqr(ip->vsite.b));
                        break;
                    case F_VSITEN:
                        /* Exact */
                        inv_mass = 0;
                        for (j = 0; j < 3*ip->vsiten.n; j += 3)
                        {
                            aj    = il->iatoms[i+j+2];
                            coeff = ip[il->iatoms[i+j]].vsiten.a;
                            if (moltype->atoms.atom[aj].ptype == eptVSite)
                            {
                                m_aj = vsite_m[aj];
                            }
                            else
                            {
                                m_aj = moltype->atoms.atom[aj].m;
                            }
                            if (m_aj <= 0)
                            {
                                gmx_incons("The mass of a vsiten constructing atom is <= 0");
                            }
                            inv_mass += coeff*coeff/m_aj;
                        }
                        vsite_m[a1] = 1/inv_mass;
                        break;
                    default:
                        /* Use the mass of the lightest constructing atom.
                         * This is an approximation.
                         * If the distance of the virtual site to the
                         * constructing atom is less than all distances
                         * between constructing atoms, this is a safe
                         * over-estimate of the displacement of the vsite.
                         * This condition holds for all H mass replacement
                         * vsite constructions, except for SP2/3 groups.
                         * In SP3 groups one H will have a F_VSITE3
                         * construction, so even there the total drift
                         * estimate shouldn't be far off.
                         */
                        assert(j >= 1);
                        vsite_m[a1] = cam[1];
                        for (j = 2; j < NRAL(ft); j++)
                        {
                            vsite_m[a1] = min(vsite_m[a1], cam[j]);
                        }
                        (*n_nonlin_vsite)++;
                        break;
                }
                if (gmx_debug_at)
                {
                    fprintf(debug, "atom %4d %-20s mass %6.3f\n",
                            a1, interaction_function[ft].longname, vsite_m[a1]);
                }
            }
        }
    }
}

static void get_verlet_buffer_atomtypes(const gmx_mtop_t      *mtop,
                                        verletbuf_atomtype_t **att_p,
                                        int                   *natt_p,
                                        int                   *n_nonlin_vsite)
{
    verletbuf_atomtype_t          *att;
    int                            natt;
    int                            mb, nmol, ft, i, a1, a2, a3, a;
    const t_atoms                 *atoms;
    const t_ilist                 *il;
    const t_iparams               *ip;
    atom_nonbonded_kinetic_prop_t *prop;
    real                          *vsite_m;
    int                            n_nonlin_vsite_mol;

    att  = NULL;
    natt = 0;

    if (n_nonlin_vsite != NULL)
    {
        *n_nonlin_vsite = 0;
    }

    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        nmol = mtop->molblock[mb].nmol;

        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;

        /* Check for constraints, as they affect the kinetic energy.
         * For virtual sites we need the masses and geometry of
         * the constructing atoms to determine their velocity distribution.
         */
        snew(prop, atoms->nr);
        snew(vsite_m, atoms->nr);

        for (ft = F_CONSTR; ft <= F_CONSTRNC; ft++)
        {
            il = &mtop->moltype[mtop->molblock[mb].type].ilist[ft];

            for (i = 0; i < il->nr; i += 1+NRAL(ft))
            {
                ip         = &mtop->ffparams.iparams[il->iatoms[i]];
                a1         = il->iatoms[i+1];
                a2         = il->iatoms[i+2];
                if (atoms->atom[a2].m > prop[a1].con_mass)
                {
                    prop[a1].con_mass = atoms->atom[a2].m;
                    prop[a1].con_len  = ip->constr.dA;
                }
                if (atoms->atom[a1].m > prop[a2].con_mass)
                {
                    prop[a2].con_mass = atoms->atom[a1].m;
                    prop[a2].con_len  = ip->constr.dA;
                }
            }
        }

        il = &mtop->moltype[mtop->molblock[mb].type].ilist[F_SETTLE];

        for (i = 0; i < il->nr; i += 1+NRAL(F_SETTLE))
        {
            ip         = &mtop->ffparams.iparams[il->iatoms[i]];
            a1         = il->iatoms[i+1];
            a2         = il->iatoms[i+2];
            a3         = il->iatoms[i+3];
            /* Usually the mass of a1 (usually oxygen) is larger than a2/a3.
             * If this is not the case, we overestimate the displacement,
             * which leads to a larger buffer (ok since this is an exotic case).
             */
            prop[a1].con_mass = atoms->atom[a2].m;
            prop[a1].con_len  = ip->settle.doh;

            prop[a2].con_mass = atoms->atom[a1].m;
            prop[a2].con_len  = ip->settle.doh;

            prop[a3].con_mass = atoms->atom[a1].m;
            prop[a3].con_len  = ip->settle.doh;
        }

        get_vsite_masses(&mtop->moltype[mtop->molblock[mb].type],
                         &mtop->ffparams,
                         vsite_m,
                         &n_nonlin_vsite_mol);
        if (n_nonlin_vsite != NULL)
        {
            *n_nonlin_vsite += nmol*n_nonlin_vsite_mol;
        }

        for (a = 0; a < atoms->nr; a++)
        {
            if (atoms->atom[a].ptype == eptVSite)
            {
                prop[a].mass = vsite_m[a];
            }
            else
            {
                prop[a].mass = atoms->atom[a].m;
            }
            prop[a].type     = atoms->atom[a].type;
            prop[a].q        = atoms->atom[a].q;
             /* We consider an atom constrained, #DOF=2, when it is
             * connected with constraints to (at least one) atom with
             * a mass of more than 0.4x its own mass. This is not a critical
             * parameter, since with roughly equal masses the unconstrained
             * and constrained displacement will not differ much (and both
             * overestimate the displacement).
             */
            prop[a].bConstr = (prop[a].con_mass > 0.4*prop[a].mass);

            add_at(&att, &natt, &prop[a], nmol);
        }

        sfree(vsite_m);
        sfree(prop);
    }

    if (gmx_debug_at)
    {
        for (a = 0; a < natt; a++)
        {
            fprintf(debug, "type %d: m %5.2f t %d q %6.3f con %d con_m %5.3f con_l %5.3f n %d\n",
                    a, att[a].prop.mass, att[a].prop.type, att[a].prop.q,
                    att[a].prop.bConstr, att[a].prop.con_mass, att[a].prop.con_len, 
                    att[a].n);
        }
    }

    *att_p  = att;
    *natt_p = natt;
}

/* This function computes two components of the estimate of the variance
 * in the displacement of one atom in a system of two constrained atoms.
 * Returns in sigma2_2d the variance due to rotation of the constrained
 * atom around the atom to which it constrained.
 * Returns in sigma2_3d the variance due to displacement of the COM
 * of the whole system of the two constrained atoms.
 *
 * Note that we only take a single constraint (the one to the heaviest atom)
 * into account. If an atom has multiple constraints, this will result in
 * an overestimate of the displacement, which gives a larger drift and buffer.
 */
static void constrained_atom_sigma2(real kT_fac,
                                    const atom_nonbonded_kinetic_prop_t *prop,
                                    real *sigma2_2d,
                                    real *sigma2_3d)
{
    real sigma2_rot;
    real com_dist;
    real sigma2_rel;
    real scale;

    /* Here we decompose the motion of a constrained atom into two
     * components: rotation around the COM and translation of the COM.
     */

    /* Determine the variance for the displacement of the rotational mode */
    sigma2_rot = kT_fac/(prop->mass*(prop->mass + prop->con_mass)/prop->con_mass);

    /* The distance from the atom to the COM, i.e. the rotational arm */
    com_dist = prop->con_len*prop->con_mass/(prop->mass + prop->con_mass);

    /* The variance relative to the arm */
    sigma2_rel = sigma2_rot/(com_dist*com_dist);
    /* At 6 the scaling formula has slope 0,
     * so we keep sigma2_2d constant after that.
     */
    if (sigma2_rel < 6)
    {
        /* A constrained atom rotates around the atom it is constrained to.
         * This results in a smaller linear displacement than for a free atom.
         * For a perfectly circular displacement, this lowers the displacement
         * by: 1/arcsin(arc_length)
         * and arcsin(x) = 1 + x^2/6 + ...
         * For sigma2_rel<<1 the displacement distribution is erfc
         * (exact formula is provided below). For larger sigma, it is clear
         * that the displacement can't be larger than 2*com_dist.
         * It turns out that the distribution becomes nearly uniform.
         * For intermediate sigma2_rel, scaling down sigma with the third
         * order expansion of arcsin with argument sigma_rel turns out
         * to give a very good approximation of the distribution and variance.
         * Even for larger values, the variance is only slightly overestimated.
         * Note that the most relevant displacements are in the long tail.
         * This rotation approximation always overestimates the tail (which
         * runs to infinity, whereas it should be <= 2*com_dist).
         * Thus we always overestimate the drift and the buffer size.
         */
        scale      = 1/(1 + sigma2_rel/6);
        *sigma2_2d = sigma2_rot*scale*scale;
    }
    else
    {
        /* sigma_2d is set to the maximum given by the scaling above.
         * For large sigma2 the real displacement distribution is close
         * to uniform over -2*con_len to 2*com_dist.
         * Our erfc with sigma_2d=sqrt(1.5)*com_dist (which means the sigma
         * of the erfc output distribution is con_dist) overestimates
         * the variance and additionally has a long tail. This means
         * we have a (safe) overestimation of the drift.
         */
        *sigma2_2d = 1.5*com_dist*com_dist;
    }

    /* The constrained atom also moves (in 3D) with the COM of both atoms */
    *sigma2_3d = kT_fac/(prop->mass + prop->con_mass);
}

static void get_atom_sigma2(real kT_fac,
                            const atom_nonbonded_kinetic_prop_t *prop,
                            real *sigma2_2d,
                            real *sigma2_3d)
{
    if (prop->bConstr)
    {
        /* Complicated constraint calculation in a separate function */
        constrained_atom_sigma2(kT_fac, prop, sigma2_2d, sigma2_3d);
    }
    else
    {
        /* Unconstrained atom: trivial */
        *sigma2_2d = 0;
        *sigma2_3d = kT_fac/prop->mass;
    }
}

static void approx_2dof(real s2, real x, real *shift, real *scale)
{
    /* A particle with 1 DOF constrained has 2 DOFs instead of 3.
     * This code is also used for particles with multiple constraints,
     * in which case we overestimate the displacement.
     * The 2DOF distribution is sqrt(pi/2)*erfc(r/(sqrt(2)*s))/(2*s).
     * We approximate this with scale*Gaussian(s,r+shift),
     * by matching the distribution value and derivative at x.
     * This is a tight overestimate for all r>=0 at any s and x.
     */
    real ex, er;

    ex = exp(-x*x/(2*s2));
    er = gmx_erfc(x/sqrt(2*s2));

    *shift = -x + sqrt(2*s2/M_PI)*ex/er;
    *scale = 0.5*M_PI*exp(ex*ex/(M_PI*er*er))*er;
}

static real ener_drift(const verletbuf_atomtype_t *att, int natt,
                       const gmx_ffparams_t *ffp,
                       real kT_fac,
                       real md_ljd, real md_ljr, real md_el, real dd_el,
                       real r_buffer,
                       real rlist, real boxvol)
{
    double drift_tot, pot1, pot2, pot;
    int    i, j;
    real   s2i_2d, s2i_3d, s2j_2d, s2j_3d, s2, s;
    int    ti, tj;
    real   md, dd;
    real   sc_fac, rsh;
    double c_exp, c_erfc;

    drift_tot = 0;

    /* Loop over the different atom type pairs */
    for (i = 0; i < natt; i++)
    {
        get_atom_sigma2(kT_fac, &att[i].prop, &s2i_2d, &s2i_3d);
        ti = att[i].prop.type;

        for (j = i; j < natt; j++)
        {
            get_atom_sigma2(kT_fac, &att[j].prop, &s2j_2d, &s2j_3d);
            tj = att[j].prop.type;

            /* Add up the up to four independent variances */
            s2 = s2i_2d + s2i_3d + s2j_2d + s2j_3d; 

            /* Note that attractive and repulsive potentials for individual
             * pairs will partially cancel.
             */
            /* -dV/dr at the cut-off for LJ + Coulomb */
            md =
                md_ljd*ffp->iparams[ti*ffp->atnr+tj].lj.c6 +
                md_ljr*ffp->iparams[ti*ffp->atnr+tj].lj.c12 +
                md_el*att[i].prop.q*att[j].prop.q;

            /* d2V/dr2 at the cut-off for Coulomb, we neglect LJ */
            dd = dd_el*att[i].prop.q*att[j].prop.q;

            rsh    = r_buffer;
            sc_fac = 1.0;
            /* For constraints: adapt r and scaling for the Gaussian */
            if (att[i].prop.bConstr)
            {
                real sh, sc;

                approx_2dof(s2i_2d, r_buffer*s2i_2d/s2, &sh, &sc);
                rsh    += sh;
                sc_fac *= sc;
            }
            if (att[j].prop.bConstr)
            {
                real sh, sc;

                approx_2dof(s2j_2d, r_buffer*s2j_2d/s2, &sh, &sc);
                rsh    += sh;
                sc_fac *= sc;
            }

            /* Exact contribution of an atom pair with Gaussian displacement
             * with sigma s to the energy drift for a potential with
             * derivative -md and second derivative dd at the cut-off.
             * The only catch is that for potentials that change sign
             * near the cut-off there could be an unlucky compensation
             * of positive and negative energy drift.
             * Such potentials are extremely rare though.
             *
             * Note that pot has unit energy*length, as the linear
             * atom density still needs to be put in.
             */
            c_exp  = exp(-rsh*rsh/(2*s2))/sqrt(2*M_PI);
            c_erfc = 0.5*gmx_erfc(rsh/(sqrt(2*s2)));
            s      = sqrt(s2);

            pot1 = sc_fac*
                md/2*((rsh*rsh + s2)*c_erfc - rsh*s*c_exp);
            pot2 = sc_fac*
                dd/6*(s*(rsh*rsh + 2*s2)*c_exp - rsh*(rsh*rsh + 3*s2)*c_erfc);
            pot = pot1 + pot2;

            if (gmx_debug_at)
            {
                fprintf(debug, "n %d %d d s %.3f %.3f %.3f %.3f con %d md %8.1e dd %8.1e pot1 %8.1e pot2 %8.1e pot %8.1e\n",
                        att[i].n, att[j].n,
                        sqrt(s2i_2d), sqrt(s2i_3d),
                        sqrt(s2j_2d), sqrt(s2j_3d),
                        att[i].prop.bConstr+att[j].prop.bConstr,
                        md, dd, pot1, pot2, pot);
            }

            /* Multiply by the number of atom pairs */
            if (j == i)
            {
                pot *= (double)att[i].n*(att[i].n - 1)/2;
            }
            else
            {
                pot *= (double)att[i].n*att[j].n;
            }
            /* We need the line density to get the energy drift of the system.
             * The effective average r^2 is close to (rlist+sigma)^2.
             */
            pot *= 4*M_PI*sqr(rlist + s)/boxvol;

            /* Add the unsigned drift to avoid cancellation of errors */
            drift_tot += fabs(pot);
        }
    }

    return drift_tot;
}

static real surface_frac(int cluster_size, real particle_distance, real rlist)
{
    real d, area_rel;

    if (rlist < 0.5*particle_distance)
    {
        /* We have non overlapping spheres */
        return 1.0;
    }

    /* Half the inter-particle distance relative to rlist */
    d = 0.5*particle_distance/rlist;

    /* Determine the area of the surface at distance rlist to the closest
     * particle, relative to surface of a sphere of radius rlist.
     * The formulas below assume close to cubic cells for the pair search grid,
     * which the pair search code tries to achieve.
     * Note that in practice particle distances will not be delta distributed,
     * but have some spread, often involving shorter distances,
     * as e.g. O-H bonds in a water molecule. Thus the estimates below will
     * usually be slightly too high and thus conservative.
     */
    switch (cluster_size)
    {
        case 1:
            /* One particle: trivial */
            area_rel = 1.0;
            break;
        case 2:
            /* Two particles: two spheres at fractional distance 2*a */
            area_rel = 1.0 + d;
            break;
        case 4:
            /* We assume a perfect, symmetric tetrahedron geometry.
             * The surface around a tetrahedron is too complex for a full
             * analytical solution, so we use a Taylor expansion.
             */
            area_rel = (1.0 + 1/M_PI*(6*acos(1/sqrt(3))*d +
                                      sqrt(3)*d*d*(1.0 +
                                                   5.0/18.0*d*d +
                                                   7.0/45.0*d*d*d*d +
                                                   83.0/756.0*d*d*d*d*d*d)));
            break;
        default:
            gmx_incons("surface_frac called with unsupported cluster_size");
            area_rel = 1.0;
    }

    return area_rel/cluster_size;
}

void calc_verlet_buffer_size(const gmx_mtop_t *mtop, real boxvol,
                             const t_inputrec *ir, real drift_target,
                             const verletbuf_list_setup_t *list_setup,
                             int *n_nonlin_vsite,
                             real *rlist)
{
    double                resolution;
    char                 *env;

    real                  particle_distance;
    real                  nb_clust_frac_pairs_not_in_list_at_cutoff;

    verletbuf_atomtype_t *att  = NULL;
    int                   natt = -1, i;
    double                reppow;
    real                  md_ljd, md_ljr, md_el, dd_el;
    real                  elfac;
    real                  kT_fac, mass_min;
    int                   ib0, ib1, ib;
    real                  rb, rl;
    real                  drift;

    /* Resolution of the buffer size */
    resolution = 0.001;

    env = getenv("GMX_VERLET_BUFFER_RES");
    if (env != NULL)
    {
        sscanf(env, "%lf", &resolution);
    }

    /* In an atom wise pair-list there would be no pairs in the list
     * beyond the pair-list cut-off.
     * However, we use a pair-list of groups vs groups of atoms.
     * For groups of 4 atoms, the parallelism of SSE instructions, only
     * 10% of the atoms pairs are not in the list just beyond the cut-off.
     * As this percentage increases slowly compared to the decrease of the
     * Gaussian displacement distribution over this range, we can simply
     * reduce the drift by this fraction.
     * For larger groups, e.g. of 8 atoms, this fraction will be lower,
     * so then buffer size will be on the conservative (large) side.
     *
     * Note that the formulas used here do not take into account
     * cancellation of errors which could occur by missing both
     * attractive and repulsive interactions.
     *
     * The only major assumption is homogeneous particle distribution.
     * For an inhomogeneous system, such as a liquid-vapor system,
     * the buffer will be underestimated. The actual energy drift
     * will be higher by the factor: local/homogeneous particle density.
     *
     * The results of this estimate have been checked againt simulations.
     * In most cases the real drift differs by less than a factor 2.
     */

    /* Worst case assumption: HCP packing of particles gives largest distance */
    particle_distance = pow(boxvol*sqrt(2)/mtop->natoms, 1.0/3.0);

    get_verlet_buffer_atomtypes(mtop, &att, &natt, n_nonlin_vsite);
    assert(att != NULL && natt >= 0);

    if (debug)
    {
        fprintf(debug, "particle distance assuming HCP packing: %f nm\n",
                particle_distance);
        fprintf(debug, "energy drift atom types: %d\n", natt);
    }

    reppow = mtop->ffparams.reppow;
    md_ljd = 0;
    md_ljr = 0;
    if (ir->vdwtype == evdwCUT)
    {
        /* -dV/dr of -r^-6 and r^-repporw */
        md_ljd = -6*pow(ir->rvdw, -7.0);
        md_ljr = reppow*pow(ir->rvdw, -(reppow+1));
        /* The contribution of the second derivative is negligible */
    }
    else
    {
        gmx_fatal(FARGS, "Energy drift calculation is only implemented for plain cut-off Lennard-Jones interactions");
    }

    elfac = ONE_4PI_EPS0/ir->epsilon_r;

    /* Determine md=-dV/dr and dd=d^2V/dr^2 */
    md_el = 0;
    dd_el = 0;
    if (ir->coulombtype == eelCUT || EEL_RF(ir->coulombtype))
    {
        real eps_rf, k_rf;

        if (ir->coulombtype == eelCUT)
        {
            eps_rf = 1;
            k_rf   = 0;
        }
        else
        {
            eps_rf = ir->epsilon_rf/ir->epsilon_r;
            if (eps_rf != 0)
            {
                k_rf = pow(ir->rcoulomb, -3.0)*(eps_rf - ir->epsilon_r)/(2*eps_rf + ir->epsilon_r);
            }
            else
            {
                /* epsilon_rf = infinity */
                k_rf = 0.5*pow(ir->rcoulomb, -3.0);
            }
        }

        if (eps_rf > 0)
        {
            md_el = elfac*(pow(ir->rcoulomb, -2.0) - 2*k_rf*ir->rcoulomb);
        }
        dd_el = elfac*(2*pow(ir->rcoulomb, -3.0) + 2*k_rf);
    }
    else if (EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD)
    {
        real b, rc, br;

        b     = calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
        rc    = ir->rcoulomb;
        br    = b*rc;
        md_el = elfac*(b*exp(-br*br)*M_2_SQRTPI/rc + gmx_erfc(br)/(rc*rc));
        dd_el = elfac/(rc*rc)*(2*b*(1 + br*br)*exp(-br*br)*M_2_SQRTPI + 2*gmx_erfc(br)/rc);
    }
    else
    {
        gmx_fatal(FARGS, "Energy drift calculation is only implemented for Reaction-Field and Ewald electrostatics");
    }

    /* Determine the variance of the atomic displacement
     * over nstlist-1 steps: kT_fac
     * For inertial dynamics (not Brownian dynamics) the mass factor
     * is not included in kT_fac, it is added later.
     */
    if (ir->eI == eiBD)
    {
        /* Get the displacement distribution from the random component only.
         * With accurate integration the systematic (force) displacement
         * should be negligible (unless nstlist is extremely large, which
         * you wouldn't do anyhow).
         */
        kT_fac = 2*BOLTZ*ir->opts.ref_t[0]*(ir->nstlist-1)*ir->delta_t;
        if (ir->bd_fric > 0)
        {
            /* This is directly sigma^2 of the displacement */
            kT_fac /= ir->bd_fric;

            /* Set the masses to 1 as kT_fac is the full sigma^2,
             * but we divide by m in ener_drift().
             */
            for (i = 0; i < natt; i++)
            {
                att[i].prop.mass = 1;
            }
        }
        else
        {
            real tau_t;

            /* Per group tau_t is not implemented yet, use the maximum */
            tau_t = ir->opts.tau_t[0];
            for (i = 1; i < ir->opts.ngtc; i++)
            {
                tau_t = max(tau_t, ir->opts.tau_t[i]);
            }

            kT_fac *= tau_t;
            /* This kT_fac needs to be divided by the mass to get sigma^2 */
        }
    }
    else
    {
        kT_fac = BOLTZ*ir->opts.ref_t[0]*sqr((ir->nstlist-1)*ir->delta_t);
    }

    mass_min = att[0].prop.mass;
    for (i = 1; i < natt; i++)
    {
        mass_min = min(mass_min, att[i].prop.mass);
    }

    if (debug)
    {
        fprintf(debug, "md_ljd %e md_ljr %e\n", md_ljd, md_ljr);
        fprintf(debug, "md_el %e dd_el %e\n", md_el, dd_el);
        fprintf(debug, "sqrt(kT_fac) %f\n", sqrt(kT_fac));
        fprintf(debug, "mass_min %f\n", mass_min);
    }

    /* Search using bisection */
    ib0 = -1;
    /* The drift will be neglible at 5 times the max sigma */
    ib1 = (int)(5*2*sqrt(kT_fac/mass_min)/resolution) + 1;
    while (ib1 - ib0 > 1)
    {
        ib = (ib0 + ib1)/2;
        rb = ib*resolution;
        rl = max(ir->rvdw, ir->rcoulomb) + rb;

        /* Calculate the average energy drift at the last step
         * of the nstlist steps at which the pair-list is used.
         */
        drift = ener_drift(att, natt, &mtop->ffparams,
                           kT_fac,
                           md_ljd, md_ljr, md_el, dd_el, rb,
                           rl, boxvol);

        /* Correct for the fact that we are using a Ni x Nj particle pair list
         * and not a 1 x 1 particle pair list. This reduces the drift.
         */
        /* We don't have a formula for 8 (yet), use 4 which is conservative */
        nb_clust_frac_pairs_not_in_list_at_cutoff =
            surface_frac(min(list_setup->cluster_size_i, 4),
                         particle_distance, rl)*
            surface_frac(min(list_setup->cluster_size_j, 4),
                         particle_distance, rl);
        drift *= nb_clust_frac_pairs_not_in_list_at_cutoff;

        /* Convert the drift to drift per unit time per atom */
        drift /= ir->nstlist*ir->delta_t*mtop->natoms;

        if (debug)
        {
            fprintf(debug, "ib %3d %3d %3d rb %.3f %dx%d fac %.3f drift %f\n",
                    ib0, ib, ib1, rb,
                    list_setup->cluster_size_i, list_setup->cluster_size_j,
                    nb_clust_frac_pairs_not_in_list_at_cutoff,
                    drift);
        }

        if (fabs(drift) > drift_target)
        {
            ib0 = ib;
        }
        else
        {
            ib1 = ib;
        }
    }

    sfree(att);

    *rlist = max(ir->rvdw, ir->rcoulomb) + ib1*resolution;
}
