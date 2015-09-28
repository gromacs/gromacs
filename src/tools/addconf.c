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

#include <stdlib.h>
#include <string.h>
#include "vec.h"
#include "macros.h"
#include "smalloc.h"
#include "addconf.h"
#include "force.h"
#include "gstat.h"
#include "names.h"
#include "nsgrid.h"
#include "mdatoms.h"
#include "nrnb.h"
#include "ns.h"
#include "mtop_util.h"
#include "chargegroup.h"

static real box_margin;

static real max_dist(rvec *x, real *r, int start, int end)
{
    real maxd;
    int  i, j;

    maxd = 0;
    for (i = start; i < end; i++)
    {
        for (j = i+1; j < end; j++)
        {
            maxd = max(maxd, sqrt(distance2(x[i], x[j]))+0.5*(r[i]+r[j]));
        }
    }

    return 0.5*maxd;
}

static gmx_bool outside_box_minus_margin2(rvec x, matrix box)
{
    return ( (x[XX] < 2*box_margin) || (x[XX] > box[XX][XX]-2*box_margin) ||
             (x[YY] < 2*box_margin) || (x[YY] > box[YY][YY]-2*box_margin) ||
             (x[ZZ] < 2*box_margin) || (x[ZZ] > box[ZZ][ZZ]-2*box_margin) );
}

static gmx_bool outside_box_plus_margin(rvec x, matrix box)
{
    return ( (x[XX] < -box_margin) || (x[XX] > box[XX][XX]+box_margin) ||
             (x[YY] < -box_margin) || (x[YY] > box[YY][YY]+box_margin) ||
             (x[ZZ] < -box_margin) || (x[ZZ] > box[ZZ][ZZ]+box_margin) );
}

static int mark_res(int at, gmx_bool *mark, int natoms, t_atom *atom, int *nmark)
{
    int resind;

    resind = atom[at].resind;
    while ( (at > 0) && (resind == atom[at-1].resind) )
    {
        at--;
    }
    while ( (at < natoms) && (resind == atom[at].resind) )
    {
        if (!mark[at])
        {
            mark[at] = TRUE;
            (*nmark)++;
        }
        at++;
    }

    return at;
}

static real find_max_real(int n, real radius[])
{
    int  i;
    real rmax;

    rmax = 0;
    if (n > 0)
    {
        rmax = radius[0];
        for (i = 1; (i < n); i++)
        {
            rmax = max(rmax, radius[i]);
        }
    }
    return rmax;
}

static void combine_atoms(t_atoms *ap, t_atoms *as,
                          rvec xp[], rvec *vp, rvec xs[], rvec *vs,
                          t_atoms **a_comb, rvec **x_comb, rvec **v_comb)
{
    t_atoms *ac;
    rvec    *xc, *vc = NULL;
    int      i, j, natot, res0;

    /* Total number of atoms */
    natot = ap->nr+as->nr;

    snew(ac, 1);
    init_t_atoms(ac, natot, FALSE);

    snew(xc, natot);
    if (vp && vs)
    {
        snew(vc, natot);
    }

    /* Fill the new structures */
    for (i = j = 0; (i < ap->nr); i++, j++)
    {
        copy_rvec(xp[i], xc[j]);
        if (vc)
        {
            copy_rvec(vp[i], vc[j]);
        }
        memcpy(&(ac->atom[j]), &(ap->atom[i]), sizeof(ap->atom[i]));
        ac->atom[j].type = 0;
    }
    res0 = ap->nres;
    for (i = 0; (i < as->nr); i++, j++)
    {
        copy_rvec(xs[i], xc[j]);
        if (vc)
        {
            copy_rvec(vs[i], vc[j]);
        }
        memcpy(&(ac->atom[j]), &(as->atom[i]), sizeof(as->atom[i]));
        ac->atom[j].type    = 0;
        ac->atom[j].resind += res0;
    }
    ac->nr   = j;
    ac->nres = ac->atom[j-1].resind+1;
    /* Fill all elements to prevent uninitialized memory */
    for (i = 0; i < ac->nr; i++)
    {
        ac->atom[i].m     = 1;
        ac->atom[i].q     = 0;
        ac->atom[i].mB    = 1;
        ac->atom[i].qB    = 0;
        ac->atom[i].type  = 0;
        ac->atom[i].typeB = 0;
        ac->atom[i].ptype = eptAtom;
    }

    /* Return values */
    *a_comb = ac;
    *x_comb = xc;
    *v_comb = vc;
}

static t_forcerec *fr = NULL;

void do_nsgrid(FILE *fp, gmx_bool bVerbose,
               matrix box, rvec x[], t_atoms *atoms, real rlong,
               const output_env_t oenv)
{
    gmx_mtop_t     *mtop;
    gmx_localtop_t *top;
    t_mdatoms      *md;
    t_block        *cgs;
    t_inputrec     *ir;
    t_nrnb          nrnb;
    t_commrec      *cr;
    int            *cg_index;
    gmx_moltype_t  *molt;
    gmx_ffparams_t *ffp;
    ivec           *nFreeze;
    int             i, m, natoms;
    rvec            box_size;
    real           *lambda, *dvdl;

    natoms = atoms->nr;

    /* Charge group index */
    snew(cg_index, natoms);
    for (i = 0; (i < natoms); i++)
    {
        cg_index[i] = i;
    }

    /* Topology needs charge groups and exclusions */
    snew(mtop, 1);
    init_mtop(mtop);
    mtop->natoms = natoms;
    /* Make one moltype that contains the whol system */
    mtop->nmoltype = 1;
    snew(mtop->moltype, mtop->nmoltype);
    molt        = &mtop->moltype[0];
    molt->name  = mtop->name;
    molt->atoms = *atoms;
    stupid_fill_block(&molt->cgs, mtop->natoms, FALSE);
    stupid_fill_blocka(&molt->excls, natoms);
    /* Make one molblock for the whole system */
    mtop->nmolblock = 1;
    snew(mtop->molblock, mtop->nmolblock);
    mtop->molblock[0].type       = 0;
    mtop->molblock[0].nmol       = 1;
    mtop->molblock[0].natoms_mol = natoms;
    /* Initialize a single energy group */
    mtop->groups.grps[egcENER].nr = 1;
    mtop->groups.ngrpnr[egcENER]  = 0;
    mtop->groups.grpnr[egcENER]   = NULL;

    ffp = &mtop->ffparams;

    ffp->ntypes = 1;
    ffp->atnr   = 1;
    ffp->reppow = 12;
    snew(ffp->functype, 1);
    snew(ffp->iparams, 1);
    ffp->iparams[0].lj.c6  = 1;
    ffp->iparams[0].lj.c12 = 1;

    /* inputrec structure */
    snew(ir, 1);
    ir->coulombtype = eelCUT;
    ir->vdwtype     = evdwCUT;
    ir->ndelta      = 2;
    ir->ns_type     = ensGRID;
    snew(ir->opts.egp_flags, 1);

    top = gmx_mtop_generate_local_top(mtop, ir);

    /* Some nasty shortcuts */
    cgs  = &(top->cgs);

    /* mdatoms structure */
    snew(nFreeze, 2);
    snew(md, 1);
    md = init_mdatoms(fp, mtop, FALSE);
    atoms2md(mtop, ir, 0, NULL, 0, mtop->natoms, md);
    sfree(nFreeze);

    /* forcerec structure */
    if (fr == NULL)
    {
        fr = mk_forcerec();
    }
    snew(cr, 1);
    cr->nnodes   = 1;
    /* cr->nthreads = 1; */

    /*    ir->rlist       = ir->rcoulomb = ir->rvdw = rlong;
       printf("Neighborsearching with a cut-off of %g\n",rlong);
       init_forcerec(stdout,fr,ir,top,cr,md,box,FALSE,NULL,NULL,NULL,TRUE);*/
    fr->cg0     = 0;
    fr->hcg     = top->cgs.nr;
    fr->nWatMol = 0;

    /* Prepare for neighboursearching */
    init_nrnb(&nrnb);

    /* Init things dependent on parameters */
    ir->rlistlong = ir->rlist = ir->rcoulomb = ir->rvdw = rlong;
    /* create free energy data to avoid NULLs */
    snew(ir->fepvals, 1);
    printf("Neighborsearching with a cut-off of %g\n", rlong);
    init_forcerec(stdout, oenv, fr, NULL, ir, mtop, cr, box, FALSE,
                  NULL, NULL, NULL, NULL, NULL, TRUE, -1);
    if (debug)
    {
        pr_forcerec(debug, fr, cr);
    }

    /* Calculate new stuff dependent on coords and box */
    for (m = 0; (m < DIM); m++)
    {
        box_size[m] = box[m][m];
    }
    calc_shifts(box, fr->shift_vec);
    put_charge_groups_in_box(fp, 0, cgs->nr, fr->ePBC, box, cgs, x, fr->cg_cm);

    /* Do the actual neighboursearching */
    snew(lambda, efptNR);
    snew(dvdl, efptNR);
    init_neighbor_list(fp, fr, md->homenr);
    search_neighbours(fp, fr, x, box, top,
                      &mtop->groups, cr, &nrnb, md, lambda, dvdl, NULL, TRUE, FALSE, FALSE);

    if (debug)
    {
        dump_nblist(debug, cr, fr, 0);
    }

    if (bVerbose)
    {
        fprintf(stderr, "Successfully made neighbourlist\n");
    }
}

gmx_bool bXor(gmx_bool b1, gmx_bool b2)
{
    return (b1 && !b2) || (b2 && !b1);
}

void add_conf(t_atoms *atoms, rvec **x, rvec **v, real **r, gmx_bool bSrenew,
              int ePBC, matrix box, gmx_bool bInsert,
              t_atoms *atoms_solvt, rvec *x_solvt, rvec *v_solvt, real *r_solvt,
              gmx_bool bVerbose, real rshell, int max_sol, const output_env_t oenv)
{
    t_nblist       *nlist;
    t_atoms        *atoms_all;
    real            max_vdw, *r_prot, *r_all, n2, r2, ib1, ib2;
    int             natoms_prot, natoms_solvt;
    int             i, j, jj, m, j0, j1, jjj, jnres, jnr, inr, iprot, is1, is2;
    int             prev, resnr, nresadd, d, k, ncells, maxincell;
    int             dx0, dx1, dy0, dy1, dz0, dz1;
    int             ntest, nremove, nkeep;
    rvec            dx, xi, xj, xpp, *x_all, *v_all;
    gmx_bool       *remove, *keep;
    int             bSolSol;

    natoms_prot  = atoms->nr;
    natoms_solvt = atoms_solvt->nr;
    if (natoms_solvt <= 0)
    {
        fprintf(stderr, "WARNING: Nothing to add\n");
        return;
    }

    if (ePBC == epbcSCREW)
    {
        gmx_fatal(FARGS, "Sorry, %s pbc is not yet supported", epbc_names[ePBC]);
    }

    if (bVerbose)
    {
        fprintf(stderr, "Calculating Overlap...\n");
    }

    /* Set margin around box edges to largest solvent dimension.
     * The maximum distance between atoms in a solvent molecule should
     * be calculated. At the moment a fudge factor of 3 is used.
     */
    r_prot     = *r;
    box_margin = 3*find_max_real(natoms_solvt, r_solvt);
    max_vdw    = max(3*find_max_real(natoms_prot, r_prot), box_margin);
    fprintf(stderr, "box_margin = %g\n", box_margin);

    snew(remove, natoms_solvt);

    nremove = 0;
    if (!bInsert)
    {
        for (i = 0; i < atoms_solvt->nr; i++)
        {
            if (outside_box_plus_margin(x_solvt[i], box) )
            {
                i = mark_res(i, remove, atoms_solvt->nr, atoms_solvt->atom, &nremove);
            }
        }
        fprintf(stderr, "Removed %d atoms that were outside the box\n", nremove);
    }

    /* Define grid stuff for genbox */
    /* Largest VDW radius */
    snew(r_all, natoms_prot+natoms_solvt);
    for (i = j = 0; i < natoms_prot; i++, j++)
    {
        r_all[j] = r_prot[i];
    }
    for (i = 0; i < natoms_solvt; i++, j++)
    {
        r_all[j] = r_solvt[i];
    }

    /* Combine arrays */
    combine_atoms(atoms, atoms_solvt, *x, v ? *v : NULL, x_solvt, v_solvt,
                  &atoms_all, &x_all, &v_all);

    /* Do neighboursearching step */
    do_nsgrid(stdout, bVerbose, box, x_all, atoms_all, max_vdw, oenv);

    /* check solvent with solute */
    nlist = &(fr->nblists[0].nlist_sr[eNL_VDW]);
    fprintf(stderr, "nri = %d, nrj = %d\n", nlist->nri, nlist->nrj);
    for (bSolSol = 0; (bSolSol <= (bInsert ? 0 : 1)); bSolSol++)
    {
        ntest = nremove = 0;
        fprintf(stderr, "Checking %s-Solvent overlap:",
                bSolSol ? "Solvent" : "Protein");
        for (i = 0; (i < nlist->nri && nremove < natoms_solvt); i++)
        {
            inr = nlist->iinr[i];
            j0  = nlist->jindex[i];
            j1  = nlist->jindex[i+1];
            rvec_add(x_all[inr], fr->shift_vec[nlist->shift[i]], xi);

            for (j = j0; (j < j1 && nremove < natoms_solvt); j++)
            {
                jnr = nlist->jjnr[j];
                if (jnr < 0)
                {
                    /* skip padding */
                    continue;
                }
                copy_rvec(x_all[jnr], xj);

                /* Check solvent-protein and solvent-solvent */
                is1 = inr-natoms_prot;
                is2 = jnr-natoms_prot;

                /* Check if at least one of the atoms is a solvent that is not yet
                 * listed for removal, and if both are solvent, that they are not in the
                 * same residue.
                 */
                if ((!bSolSol &&
                     bXor((is1 >= 0), (is2 >= 0)) && /* One atom is protein */
                     ((is1 < 0) || ((is1 >= 0) && !remove[is1])) &&
                     ((is2 < 0) || ((is2 >= 0) && !remove[is2]))) ||

                    (bSolSol  &&
                     (is1 >= 0) && (!remove[is1]) &&                   /* is1 is solvent */
                     (is2 >= 0) && (!remove[is2]) &&                   /* is2 is solvent */
                     (bInsert ||                                       /* when inserting also check inside the box */
                      (outside_box_minus_margin2(x_solvt[is1], box) && /* is1 on edge */
                       outside_box_minus_margin2(x_solvt[is2], box))   /* is2 on edge */
                     ) &&
                     (atoms_solvt->atom[is1].resind !=                 /* Not the same residue */
                      atoms_solvt->atom[is2].resind)))
                {

                    ntest++;
                    rvec_sub(xi, xj, dx);
                    n2 = norm2(dx);
                    r2 = sqr(r_all[inr]+r_all[jnr]);
                    if (n2 < r2)
                    {
                        if (bInsert)
                        {
                            nremove = natoms_solvt;
                            for (k = 0; k < nremove; k++)
                            {
                                remove[k] = TRUE;
                            }
                        }
                        /* Need only remove one of the solvents... */
                        if (is2 >= 0)
                        {
                            (void) mark_res(is2, remove, natoms_solvt, atoms_solvt->atom,
                                            &nremove);
                        }
                        else if (is1 >= 0)
                        {
                            (void) mark_res(is1, remove, natoms_solvt, atoms_solvt->atom,
                                            &nremove);
                        }
                        else
                        {
                            fprintf(stderr, "Neither atom is solvent%d %d\n", is1, is2);
                        }
                    }
                }
            }
        }
        if (!bInsert)
        {
            fprintf(stderr, " tested %d pairs, removed %d atoms.\n", ntest, nremove);
        }
    }
    if (debug)
    {
        for (i = 0; i < natoms_solvt; i++)
        {
            fprintf(debug, "remove[%5d] = %s\n", i, bool_names[remove[i]]);
        }
    }

    /* Search again, now with another cut-off */
    if (rshell > 0)
    {
        do_nsgrid(stdout, bVerbose, box, x_all, atoms_all, rshell, oenv);
        nlist = &(fr->nblists[0].nlist_sr[eNL_VDW]);
        fprintf(stderr, "nri = %d, nrj = %d\n", nlist->nri, nlist->nrj);
        nkeep = 0;
        snew(keep, natoms_solvt);
        for (i = 0; i < nlist->nri; i++)
        {
            inr = nlist->iinr[i];
            j0  = nlist->jindex[i];
            j1  = nlist->jindex[i+1];

            for (j = j0; j < j1; j++)
            {
                jnr = nlist->jjnr[j];

                /* Check solvent-protein and solvent-solvent */
                is1 = inr-natoms_prot;
                is2 = jnr-natoms_prot;

                /* Check if at least one of the atoms is a solvent that is not yet
                 * listed for removal, and if both are solvent, that they are not in the
                 * same residue.
                 */
                if (is1 >= 0 && is2 < 0)
                {
                    mark_res(is1, keep, natoms_solvt, atoms_solvt->atom, &nkeep);
                }
                else if (is1 < 0 && is2 >= 0)
                {
                    mark_res(is2, keep, natoms_solvt, atoms_solvt->atom, &nkeep);
                }
            }
        }
        fprintf(stderr, "Keeping %d solvent atoms after proximity check\n",
                nkeep);
        for (i = 0; i < natoms_solvt; i++)
        {
            remove[i] = remove[i] || !keep[i];
        }
        sfree(keep);
    }
    /* count how many atoms and residues will be added and make space */
    if (bInsert)
    {
        j     = atoms_solvt->nr;
        jnres = atoms_solvt->nres;
    }
    else
    {
        j     = 0;
        jnres = 0;
        for (i = 0; ((i < atoms_solvt->nr) &&
                     ((max_sol == 0) || (jnres < max_sol))); i++)
        {
            if (!remove[i])
            {
                j++;
                if ((i == 0) ||
                    (atoms_solvt->atom[i].resind != atoms_solvt->atom[i-1].resind))
                {
                    jnres++;
                }
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "Will add %d atoms in %d residues\n", j, jnres);
    }
    if (!bInsert)
    {
        /* Flag the remaing solvent atoms to be removed */
        jjj = atoms_solvt->atom[i-1].resind;
        for (; (i < atoms_solvt->nr); i++)
        {
            if (atoms_solvt->atom[i].resind > jjj)
            {
                remove[i] = TRUE;
            }
            else
            {
                j++;
            }
        }
    }

    if (bSrenew)
    {
        srenew(atoms->resinfo,  atoms->nres+jnres);
        srenew(atoms->atomname, atoms->nr+j);
        srenew(atoms->atom,     atoms->nr+j);
        srenew(*x,              atoms->nr+j);
        if (v)
        {
            srenew(*v,       atoms->nr+j);
        }
        srenew(*r,              atoms->nr+j);
    }

    /* add the selected atoms_solvt to atoms */
    if (atoms->nr > 0)
    {
        resnr = atoms->resinfo[atoms->atom[atoms->nr-1].resind].nr;
    }
    else
    {
        resnr = 0;
    }
    prev    = -1;
    nresadd = 0;
    for (i = 0; i < atoms_solvt->nr; i++)
    {
        if (!remove[i])
        {
            if (prev == -1 ||
                atoms_solvt->atom[i].resind != atoms_solvt->atom[prev].resind)
            {
                nresadd++;
                atoms->nres++;
                resnr++;
                atoms->resinfo[atoms->nres-1] =
                    atoms_solvt->resinfo[atoms_solvt->atom[i].resind];
                atoms->resinfo[atoms->nres-1].nr = resnr;
                /* calculate shift of the solvent molecule using the first atom */
                copy_rvec(x_solvt[i], dx);
                put_atoms_in_box(ePBC, box, 1, &dx);
                rvec_dec(dx, x_solvt[i]);
            }
            atoms->atom[atoms->nr]     = atoms_solvt->atom[i];
            atoms->atomname[atoms->nr] = atoms_solvt->atomname[i];
            rvec_add(x_solvt[i], dx, (*x)[atoms->nr]);
            if (v)
            {
                copy_rvec(v_solvt[i], (*v)[atoms->nr]);
            }
            (*r)[atoms->nr]               = r_solvt[i];
            atoms->atom[atoms->nr].resind = atoms->nres-1;
            atoms->nr++;
            prev = i;
        }
    }
    if (bSrenew)
    {
        srenew(atoms->resinfo,  atoms->nres+nresadd);
    }

    if (bVerbose)
    {
        fprintf(stderr, "Added %d molecules\n", nresadd);
    }

    sfree(remove);
    done_atom(atoms_all);
    sfree(x_all);
    sfree(v_all);
}
