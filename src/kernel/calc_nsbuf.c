/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.03
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <math.h>
#include "typedefs.h"
#include "physics.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "vec.h"
#include "coulomb.h"
#include "calc_nsbuf.h"

typedef struct
{
    real     mass;
    int      type;
    real     q;
    gmx_bool con;
    int      n;
} att_t;

static void add_at(att_t **att_p,int *natt_p,
                   real mass,int type,real q,gmx_bool con,int nmol)
{
    att_t *att;
    int natt,i;

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
        srenew(*att_p,*natt_p);
        (*att_p)[i].mass = mass;
        (*att_p)[i].type = type;
        (*att_p)[i].q    = q;
        (*att_p)[i].con  = con;
        (*att_p)[i].n    = nmol;

        if (gmx_debug_at)
        {
            fprintf(debug,"m %5.2f t %d q %6.3f con %d nmol %d\n",
                    mass,type,q,con,nmol);
        }
    }
}

static void get_ns_buffer_atomtypes(const gmx_mtop_t *mtop,
                                    att_t **att_p,int *natt_p,
                                    int *n_nonlin_vsite)
{
    att_t *att;
    int natt;
    int mb,nmol,ft,i,j,a1,a2,a3,a;
    const t_atoms *atoms;
    const t_ilist *il;
    const t_atom *at;
    const t_iparams *ip;
    real *con_m,*vsite_m,cam[5];

    att  = NULL;
    natt = 0;

    *n_nonlin_vsite = 0;

    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        nmol = mtop->molblock[mb].nmol;

        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;

        snew(con_m,atoms->nr);
        snew(vsite_m,atoms->nr);

        for(ft=F_CONSTR; ft<=F_CONSTRNC; ft++)
        {
            il = &mtop->moltype[mtop->molblock[mb].type].ilist[ft];

            for(i=0; i<il->nr; i+=3)
            {
                a1 = il->iatoms[i+1];
                a2 = il->iatoms[i+2];
                con_m[a1] += atoms->atom[a2].m;
                con_m[a2] += atoms->atom[a1].m;
            }
        }

        il = &mtop->moltype[mtop->molblock[mb].type].ilist[F_SETTLE];

        for(i=0; i<il->nr; i+=2)
        {
            a1 = il->iatoms[i+1];
            a2 = il->iatoms[i+1]+1;
            a3 = il->iatoms[i+1]+2;
            con_m[a1] += atoms->atom[a2].m + atoms->atom[a3].m;
            con_m[a2] += atoms->atom[a1].m + atoms->atom[a3].m;
            con_m[a3] += atoms->atom[a1].m + atoms->atom[a2].m;
        }

        for(ft=0; ft<F_NRE; ft++)
        {
            if (IS_VSITE(ft))
            {
                il = &mtop->moltype[mtop->molblock[mb].type].ilist[ft];

                for(i=0; i<il->nr; i+=1+NRAL(ft))
                {
                    ip = &mtop->ffparams.iparams[il->iatoms[i]];

                    a1 = il->iatoms[i+1];

                    for(j=1; j<NRAL(ft); j++)
                    {
                        cam[j] = atoms->atom[il->iatoms[i+1+j]].m;
                        if (cam[j] == 0)
                        {
                            cam[j] = vsite_m[il->iatoms[i+1+j]];
                        }
                        if (cam[j] == 0)
                        {
                            gmx_fatal(FARGS,"In molecule type '%s' %s construction involves atom %d, which is a virtual site of equal or high complexity. This is not supported.",
                                      *mtop->moltype[mtop->molblock[mb].type].name,
                                      interaction_function[ft].longname,
                                      il->iatoms[i+1+j]+1);
                        }
                    }

                    switch(ft)
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
                        vsite_m[a1] = cam[1];
                        for(j=2; j<NRAL(ft); j++)
                        {
                            vsite_m[a1] = min(vsite_m[a1],cam[j]);
                        }
                        *n_nonlin_vsite += nmol;
                        break;
                    }
                }
            }
        }

        for(a=0; a<atoms->nr; a++)
        {
            at = &atoms->atom[a];
            add_at(&att,&natt,
                   at->m,at->type,at->q,con_m[a] > 1.5*at->m,nmol);
        }

        sfree(vsite_m);
        sfree(con_m);
    }

    *att_p  = att;
    *natt_p = natt;
}

static real ener_drift(const att_t *att,int natt,
                       const gmx_ffparams_t *ffp,
                       real kT_fac,real vol_fac,
                       real d_ljd,real d_ljr,real d_el,real dd_el,
                       real rb)
{
    double drift_tot,drift;
    int  i,j;
    real s2i,s2j,s2,sfi,sfj,s;
    int  ti,tj;
    real d,dd;
    double c_exp,c_erfc;

    drift_tot = 0;

    for(i=0; i<natt; i++)
    {
        s2i = kT_fac/att[i].mass;
        ti  = att[i].type;

        for(j=i; j<natt; j++)
        {
            s2j = kT_fac/att[j].mass;
            tj = att[j].type;

            d =
                d_ljd*ffp->iparams[ti*ffp->atnr+tj].lj.c6 +
                d_ljr*ffp->iparams[ti*ffp->atnr+tj].lj.c12 +
                d_el*att[i].q*att[j].q;

            dd = dd_el*att[i].q*att[j].q;

            s2  = s2i + s2j;
            sfi = s2i/s2;
            sfj = s2j/s2;
            s   = sqrt(s2);

            c_exp  = exp(-rb*rb/(2*s2))/sqrt(2*M_PI);
            c_erfc = 0.5*erfc(rb/(sqrt(2*s2)));

            /* Exact contribution of a mass to the energy drift
             * for a potential with derivative -d and second derivative dd
             * at the cut-off. The only catch is that for potentials that
             * change sign near the cut-off there could be unlucky compensation
             * of positive and negative energy drift. Such potentials are
             * extremely rare.
             */
            drift =
                d*((rb*rb + s2)*c_erfc - rb*s*c_exp) +
                0.5*dd*(s*(rb*rb + 2*s2)*c_exp - rb*(rb*rb + 3*s2)*c_erfc);

            if (gmx_debug_at)
            {
                fprintf(debug,"n %d %d d s %.3f %.3f con %d %d %.3e dd %.3e drift %.3e\n",
                        att[i].n,att[j].n,sqrt(s2i),sqrt(s2j),att[i].con,att[j].con,
                        d,dd,drift);
            }

            /* Check if with constraints the distribution with 2 DOFs
             * is smaller than that with 3 DOFs at the buffer distance.
             * If so, approximate the 2 DOFs distribution with a Gaussian
             * with equal sigma and equal value at the buffer distance.
             * This is a tight overestimation.
             */
            /* Scale c_erfc to obtain the exact 2 DOFs distribution */
            c_erfc *= sqrt(0.5*M_PI);

            if (att[i].con && rb*sfi*c_exp > c_erfc*sqrt(s2i))
            {
                drift /= rb*sfi*c_exp/(c_erfc*sqrt(s2i));
            }
            if (att[j].con && rb*sfj*c_exp > c_erfc*sqrt(s2j))
            {
                drift /= rb*sfj*c_exp/(c_erfc*sqrt(s2j));
            }

            /* Add the density factor, without volume */
            if (j == i)
            {
                drift *= (double)att[i].n*(att[i].n - 1)/2;
            }
            else
            {
                drift *= (double)att[i].n*att[j].n;
            }

            drift_tot += fabs(drift);
        }
    }

    /* The volume term of the density */
    drift_tot *= vol_fac;

    return drift_tot;
}

void calc_ns_buffer_size(const gmx_mtop_t *mtop,real boxvol,
                         const t_inputrec *ir,real drift_target,
                         int *n_nonlin_vsite,
                         real *rlist)
{
    real resolution;

    real nb_cell_rel_pairs_not_in_list_at_cutoff;

    att_t *att;
    int  natt,i;
    real kT_fac,mass_min,vol_fac;
    real elfac,eps_rf,krf;
    real d_ljd,d_ljr,d_el,dd_el;
    real beta,beta_r;
    int  ib0,ib1,ib;
    real rb,rl;
    real drift;

    /* Resolution of the buffer size */
    resolution = 0.01;

    /* In an atom wise neighborlist there would be no pairs in the list
     * beyond the neighbor-list cut-off.
     * However, we use a neighborlist of groups vs groups of atoms.
     * For groups of 4 atoms, the parallelizm of SSE instructions,
     * only 10% of the atoms pairs are not in the list just after the cut-off.
     * As this percentage inceases slowly compared to the decrease of the
     * Gaussian displacement distribution over this range, we can simply
     * reduce the drift by this fraction.
     * For larger groups, e.g. of 8 atoms, this fraction will be lower,
     * but then buffer size will be on the conservative (large) side.
     */
    nb_cell_rel_pairs_not_in_list_at_cutoff = 0.1;

    get_ns_buffer_atomtypes(mtop,&att,&natt,n_nonlin_vsite);

    if (debug)
    {
        fprintf(debug,"energy drift atom types: %d\n",natt);
    }

    d_ljd = 0;
    d_ljr = 0;
    if (ir->vdwtype == evdwCUT)
    {
        d_ljd = 6*pow(ir->rvdw,-7.0);
        d_ljr = -12*pow(ir->rvdw,-13.0);
    }
    else
    {
        gmx_fatal(FARGS,"Energy drift calculation is only implemented for plain cut-off Lennard-Jones interactions");
    }

    elfac = ONE_4PI_EPS0/ir->epsilon_r;

    d_el  = 0;
    dd_el = 0;
    if (ir->coulombtype == eelCUT || EEL_RF(ir->coulombtype))
    {
        if (ir->coulombtype == eelCUT)
        {
            eps_rf = 1;
            krf = 0;
        }
        else
        {
            eps_rf = ir->epsilon_rf;
            if (eps_rf != 0)
            {
                krf = pow(ir->rcoulomb,-3.0)*(eps_rf - ir->epsilon_r)/(2*eps_rf + ir->epsilon_r);
            }
            else
            {
                /* epsilon_rf = infinity */
                krf = 0.5*pow(ir->rcoulomb,-3.0);
            }
        }

        if (eps_rf > 0)
        {
            d_el = elfac*(pow(ir->rcoulomb,-2.0) - 2*krf*ir->rcoulomb);
        }
        dd_el = elfac*2*krf;
    }
    else if (EEL_PME(ir->coulombtype) || ir->coulombtype == eelEWALD)
    {
        beta = calc_ewaldcoeff(ir->rcoulomb,ir->ewald_rtol);
        beta_r = beta*ir->rcoulomb;
        d_el = elfac*(2*beta*exp(-beta_r*beta_r)/(sqrt(M_PI)*ir->rcoulomb) + erfc(beta_r)/(ir->rcoulomb*ir->rcoulomb));
    }
    else
    {
        gmx_fatal(FARGS,"Energy drift calculation is only implemented for Reaction-Field and Ewald electrostatics");
    }

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
            for(i=0; i<natt; i++)
            {
                att[i].mass = 1;
            }
        }
        else
        {
            real tau_t;

            /* Per group tau_t is not implemented yet, use the maximum */
            tau_t = ir->opts.tau_t[0];
            for(i=1; i<ir->opts.ngtc; i++)
            {
                tau_t = max(tau_t,ir->opts.tau_t[i]);
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
    for(i=1; i<natt; i++)
    {
        mass_min = min(mass_min,att[i].mass);
    }

    if (debug)
    {
        fprintf(debug,"d_ljd %e d_ljr %e\n",d_ljd,d_ljr);
        fprintf(debug,"d_el %e dd_el %e\n",d_el,dd_el);
        fprintf(debug,"sqrt(kT_fac) %f\n",sqrt(kT_fac));
        fprintf(debug,"mass_min %f\n",mass_min);
    }

    /* Search using bisection */
    ib0 = -1;
    /* The drift should be neglible at 5 times the max sigma */
    ib1 = (int)(5*2*sqrt(kT_fac/mass_min)/resolution) + 1;
    while (ib1 - ib0 > 1)
    {
        ib = (ib0 + ib1)/2;
        rb = ib*resolution;
        rl = max(ir->rvdw,ir->rcoulomb) + rb;
        vol_fac = 4*M_PI*rl*rl/boxvol;
        
        /* Calculate the average energy drift at the last step
         * of the nstlist steps at which the neighbor list is used.
         */
        drift = ener_drift(att,natt,&mtop->ffparams,
                           kT_fac,
                           vol_fac,
                           d_ljd,d_ljr,d_el,dd_el,rb);

        /* Convert the drift to drift per unit time per atom */
        drift /= ir->nstlist*ir->delta_t*mtop->natoms;

        drift *= nb_cell_rel_pairs_not_in_list_at_cutoff;

        if (debug)
        {
            fprintf(debug,"ib %d %d %d rb %g drift %f\n",ib0,ib,ib1,rb,drift);
        }

        if (drift > drift_target)
        {
            ib0 = ib;
        }
        else
        {
            ib1 = ib;
        }
    }

    sfree(att);

    *rlist = max(ir->rvdw,ir->rcoulomb) + ib1*resolution;
}
