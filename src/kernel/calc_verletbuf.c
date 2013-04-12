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

/* Struct for unique atom type for calculating the energy drift.
 * The atom displacement depends on mass and constraints.
 * The energy jump for given distance depend on LJ type and q.
 */
typedef struct
{
    real     mass; /* mass */
    int      type; /* type (used for LJ parameters) */
    real     q;    /* charge */
    int      con;  /* constrained: 0, else 1, if 1, use #DOF=2 iso 3 */
    int      n;    /* total #atoms of this type in the system */
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

static void add_at(verletbuf_atomtype_t **att_p, int *natt_p,
                   real mass, int type, real q, int con, int nmol)
{
    verletbuf_atomtype_t *att;
    int                   natt, i;

    if (mass == 0)
    {
        /* Ignore massless particles */
        return;
    }

    att  = *att_p;
    natt = *natt_p;

    i = 0;
    while (i < natt &&
           !(mass == att[i].mass &&
             type == att[i].type &&
             q    == att[i].q &&
             con  == att[i].con))
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
        (*att_p)[i].mass = mass;
        (*att_p)[i].type = type;
        (*att_p)[i].q    = q;
        (*att_p)[i].con  = con;
        (*att_p)[i].n    = nmol;
    }
}

static void get_verlet_buffer_atomtypes(const gmx_mtop_t      *mtop,
                                        verletbuf_atomtype_t **att_p,
                                        int                   *natt_p,
                                        int                   *n_nonlin_vsite)
{
    verletbuf_atomtype_t *att;
    int                   natt;
    int                   mb, nmol, ft, i, j, a1, a2, a3, a;
    const t_atoms        *atoms;
    const t_ilist        *il;
    const t_atom         *at;
    const t_iparams      *ip;
    real                 *con_m, *vsite_m, cam[5];

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

        /* Check for constraints, as they affect the kinetic energy */
        snew(con_m, atoms->nr);
        snew(vsite_m, atoms->nr);

        for (ft = F_CONSTR; ft <= F_CONSTRNC; ft++)
        {
            il = &mtop->moltype[mtop->molblock[mb].type].ilist[ft];

            for (i = 0; i < il->nr; i += 1+NRAL(ft))
            {
                a1         = il->iatoms[i+1];
                a2         = il->iatoms[i+2];
                con_m[a1] += atoms->atom[a2].m;
                con_m[a2] += atoms->atom[a1].m;
            }
        }

        il = &mtop->moltype[mtop->molblock[mb].type].ilist[F_SETTLE];

        for (i = 0; i < il->nr; i += 1+NRAL(F_SETTLE))
        {
            a1         = il->iatoms[i+1];
            a2         = il->iatoms[i+2];
            a3         = il->iatoms[i+3];
            con_m[a1] += atoms->atom[a2].m + atoms->atom[a3].m;
            con_m[a2] += atoms->atom[a1].m + atoms->atom[a3].m;
            con_m[a3] += atoms->atom[a1].m + atoms->atom[a2].m;
        }

        /* Check for virtual sites, determine mass from constructing atoms */
        for (ft = 0; ft < F_NRE; ft++)
        {
            if (IS_VSITE(ft))
            {
                il = &mtop->moltype[mtop->molblock[mb].type].ilist[ft];

                for (i = 0; i < il->nr; i += 1+NRAL(ft))
                {
                    ip = &mtop->ffparams.iparams[il->iatoms[i]];

                    a1 = il->iatoms[i+1];

                    for (j = 1; j < NRAL(ft); j++)
                    {
                        cam[j] = atoms->atom[il->iatoms[i+1+j]].m;
                        if (cam[j] == 0)
                        {
                            cam[j] = vsite_m[il->iatoms[i+1+j]];
                        }
                        if (cam[j] == 0)
                        {
                            gmx_fatal(FARGS, "In molecule type '%s' %s construction involves atom %d, which is a virtual site of equal or high complexity. This is not supported.",
                                      *mtop->moltype[mtop->molblock[mb].type].name,
                                      interaction_function[ft].longname,
                                      il->iatoms[i+1+j]+1);
                        }
                    }

                    switch (ft)
                    {
                        case F_VSITE2:
                            /* Exact except for ignoring constraints */
                            vsite_m[a1] = (cam[2]*sqr(1-ip->vsite.a) + cam[1]*sqr(ip->vsite.a))/(cam[1]*cam[2]);
                            break;
                        case F_VSITE3:
                            /* Exact except for ignoring constraints */
                            vsite_m[a1] = (cam[2]*cam[3]*sqr(1-ip->vsite.a-ip->vsite.b) + cam[1]*cam[3]*sqr(ip->vsite.a) + cam[1]*cam[2]*sqr(ip->vsite.b))/(cam[1]*cam[2]*cam[3]);
                            break;
                        default:
                            /* Use the mass of the lightest constructing atom.
                             * This is an approximation.
                             * If the distance of the virtual site to the
                             * constructing atom is less than all distances
                             * between constructing atoms, this is a safe
                             * over-estimate of the displacement of the vsite.
                             * This condition holds for all H mass replacement
                             * replacement vsite constructions, except for SP2/3
                             * groups. In SP3 groups one H will have a F_VSITE3
                             * construction, so even there the total drift
                             * estimation shouldn't be far off.
                             */
                            assert(j >= 1);
                            vsite_m[a1] = cam[1];
                            for (j = 2; j < NRAL(ft); j++)
                            {
                                vsite_m[a1] = min(vsite_m[a1], cam[j]);
                            }
                            if (n_nonlin_vsite != NULL)
                            {
                                *n_nonlin_vsite += nmol;
                            }
                            break;
                    }
                }
            }
        }

        for (a = 0; a < atoms->nr; a++)
        {
            at = &atoms->atom[a];
            /* We consider an atom constrained, #DOF=2, when it is
             * connected with constraints to one or more atoms with
             * total mass larger than 1.5 that of the atom itself.
             */
            add_at(&att, &natt,
                   at->m, at->type, at->q, con_m[a] > 1.5*at->m, nmol);
        }

        sfree(vsite_m);
        sfree(con_m);
    }

    if (gmx_debug_at)
    {
        for (a = 0; a < natt; a++)
        {
            fprintf(debug, "type %d: m %5.2f t %d q %6.3f con %d n %d\n",
                    a, att[a].mass, att[a].type, att[a].q, att[a].con, att[a].n);
        }
    }

    *att_p  = att;
    *natt_p = natt;
}

static void approx_2dof(real s2, real x,
                        real *shift, real *scale)
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
    real   s2i, s2j, s2, s;
    int    ti, tj;
    real   md, dd;
    real   sc_fac, rsh;
    double c_exp, c_erfc;

    drift_tot = 0;

    /* Loop over the different atom type pairs */
    for (i = 0; i < natt; i++)
    {
        s2i = kT_fac/att[i].mass;
        ti  = att[i].type;

        for (j = i; j < natt; j++)
        {
            s2j = kT_fac/att[j].mass;
            tj  = att[j].type;

            /* Note that attractive and repulsive potentials for individual
             * pairs will partially cancel.
             */
            /* -dV/dr at the cut-off for LJ + Coulomb */
            md =
                md_ljd*ffp->iparams[ti*ffp->atnr+tj].lj.c6 +
                md_ljr*ffp->iparams[ti*ffp->atnr+tj].lj.c12 +
                md_el*att[i].q*att[j].q;

            /* d2V/dr2 at the cut-off for Coulomb, we neglect LJ */
            dd = dd_el*att[i].q*att[j].q;

            s2  = s2i + s2j;

            rsh    = r_buffer;
            sc_fac = 1.0;
            /* For constraints: adapt r and scaling for the Gaussian */
            if (att[i].con)
            {
                real sh, sc;
                approx_2dof(s2i, r_buffer*s2i/s2, &sh, &sc);
                rsh    += sh;
                sc_fac *= sc;
            }
            if (att[j].con)
            {
                real sh, sc;
                approx_2dof(s2j, r_buffer*s2j/s2, &sh, &sc);
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
                fprintf(debug, "n %d %d d s %.3f %.3f con %d md %8.1e dd %8.1e pot1 %8.1e pot2 %8.1e pot %8.1e\n",
                        att[i].n, att[j].n, sqrt(s2i), sqrt(s2j),
                        att[i].con+att[j].con,
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
                att[i].mass = 1;
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

    mass_min = att[0].mass;
    for (i = 1; i < natt; i++)
    {
        mass_min = min(mass_min, att[i].mass);
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
