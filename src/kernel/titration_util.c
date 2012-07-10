/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: poldata.c,v 1.20 2009/05/17 13:56:55 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#include <string.h>
#include "smalloc.h"
#include "typedefs.h"
#include "hackblock.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "constr.h"
#include "vec.h"
#include "txtdump.h"
#include "futil.h"
#include "titration_util.h"
#include "titration_parm.h"

const char *qhopregimes[etQhopNR] = { "NONE", "SE", "Intermediate", "TST" };

static void swap_rvec(rvec a,rvec b)
{
    rvec tmp;
    
    copy_rvec(a,tmp);
    copy_rvec(b,a);
    copy_rvec(tmp,b);
}

void qhop_tautomer_swap(const titration_t T,
                        rvec x[], rvec v[], rvec fbefore[], rvec fafter[],
                        int prim, int sec)
{
    int i, xv;
    t_qhop_atom *qprim, *qsec;

    qprim = &T->qhop_atoms[prim];
    qsec  = &T->qhop_atoms[sec];

    if (prim == sec)
    {
        return;
    }

    /* Swap titrating atom */
    swap_rvec(x[qprim->atom_id],x[qsec->atom_id]);

    if (v != NULL)
        swap_rvec(v[qprim->atom_id],v[qsec->atom_id]);
    if (fbefore != NULL)
        swap_rvec(fbefore[qprim->atom_id],fbefore[qsec->atom_id]);
    if (fafter != NULL)
        swap_rvec(fafter[qprim->atom_id],fafter[qsec->atom_id]);

    /* swap hydrogens */
    for (i=0; i < qsec->nr_protons; i++)
        /* qsec->nr_protons <= qprim->nr_protons */
    {
        swap_rvec(x[qprim->protons[i]],x[qsec->protons[i]]);
      
        if (v != NULL)
            swap_rvec(v[qprim->protons[i]],v[qsec->protons[i]]);
        if (fbefore != NULL)
            swap_rvec(fbefore[qprim->protons[i]],fbefore[qsec->protons[i]]);
        if (fafter != NULL)
            swap_rvec(fafter[qprim->protons[i]],fafter[qsec->protons[i]]);
    }
}

int find_inert_atomtype(const gmx_mtop_t *mtop, const t_forcerec *fr)
{
    int n, i, nvdwparam, nti, ti;
    real *vdwparam, c[3];

    /* Since we ran grompp with -norenum we know that
     * the last atomtype is the inert hydrogen.
     * This may change in future versions, however. */

    n = mtop->atomtypes.nr-1;


    /* Validate that this is indeed an inert hydrogen. */

    nvdwparam = fr->bBHAM ? 3 : 2;

    nti = nvdwparam * fr->ntype * n; /* The n-row */

    for (i=0; i<=n; i++)
    {
        vdwparam = fr->nbfp;
        ti = nti + nvdwparam * i; /* The offset, i-column */

        if (fr->bBHAM)
        {
            /* Buckingham */
            c[0]               = vdwparam[ti];
            c[1]               = vdwparam[ti+1];
            c[2]               = vdwparam[ti+2];
        }
        else
        {
            /* LJ */
            c[0]               = vdwparam[ti];
            c[1]               = vdwparam[ti+1];
            c[2]               = 0.0;
        }

        if ((c[0] != 0.0) ||
            (c[1] != 0.0) ||
            (c[2] != 0.0))
        {
            gmx_fatal(FARGS, "Atomtype %d is not inert!", n+1);
        }
    }

    return n;
}

static void qhop_attach_ilib(gmx_localtop_t *top, const qhop_db *db)
{
    int nr, i;
    qhop_subres *res;

    nr = top->idef.ntypes + db->rb.ni;

    srenew(top->idef.iparams,  nr);
    srenew(top->idef.functype, nr);

    for(i=0; (i<db->rb.ni); i++) {
        top->idef.iparams[top->idef.ntypes + i] = db->rb.ilib[i];
        top->idef.functype[top->idef.ntypes + i] = db->rb.ftlib[i];
    }

    top->idef.ntypes = nr;
}

void set_proton_presence(qhop_H_exist *Hext, const int atomid, const gmx_bool present)
{
    Hext->H[Hext->atomid2H[atomid]] = present ? (char)1: (char)0;
}

gmx_bool get_proton_presence(const qhop_H_exist *Hext, const int atomid)
{
    return Hext->H[Hext->atomid2H[atomid]] == (char)1;
}

void unindex_ilists(t_qhop_residue *qres)
{
    int i, j;
  
    /*Loop over bonded*/
    for (i=0; i<F_NRE; i++)
    {
        {
            for (j=0; j < qres->bindex.nr[i]; j++)
            {
                qres->bindex.indexed[i][j] = FALSE;
            }
        }
    }
    qres->nr_indexed = 0;
}

void index_ilists(t_qhop_residue *qres,
                  const qhop_db *db,
                  const gmx_localtop_t *top,
                  const t_commrec *cr)
{
    int rt, b, i, prev, ni, *ra, ia, ft, itype, nia, *l2g, bt, j;
    gmx_bool bMatch, bDD;
    char b2[STRLEN],buf[STRLEN];
  
    const t_ilist *ilist = top->idef.il;

    if (db->rb.qrt[qres->rtype].bWater)
    {
        return;
    }

    if (top == NULL)
    {
        gmx_fatal(FARGS, "NULL topology passed to index_ilists().");
    }

    /* bDD = DOMAINDECOMP(cr); */

    rt = qres->rtype;

    /* l2g = bDD ? cr->dd->gatindex : NULL; /\* local to global atom translation *\/ */

    unindex_ilists(qres);

    /* Do the local top */
    for (bt=0; bt < ebtsNR; bt++)
    {

        itype = db->rb.btype[bt];
        if (itype < 0)
        {
            /* No bondeds of this type */
            continue;
        }

        nia = interaction_function[itype].nratoms;

        for (b=0; b < qres->bindex.nr[itype]; b++)
        {
            /* get residue-local atom numbers */
            ra = db->rb.qrt[rt].ba[bt][b];
	  
            /* Loop the ilist */

            i=0;

            if (itype == F_PDIHS)
            {
                /* If we have dihedrals of type 9 we may
                 * have several instances of the same
                 * set of atoms. We can avoid indexing the
                 * same ilist entry by starting our search
                 * at the previous instace if there is one,
                 * or at 0 if there is none.
                 */
                bMatch = FALSE;

                for (prev=b-1; prev>=0; prev--)
                {
                    bMatch = TRUE;

                    for (j=0; j<nia; j++)
                    {
                        if (ra[j] != db->rb.qrt[rt].ba[bt][prev][j])
                        {
                            bMatch = FALSE;
                            break;
                        }
                    }

                    if (!bMatch)
                    {
                        /* No? Try the atoms in reverse order. */
                        bMatch = TRUE;
                        for (j=0; j<nia; j++)
                        {
                            if (ra[j] != db->rb.qrt[rt].ba[bt][prev][nia-j-1])
                            {
                                bMatch = FALSE;
                                break;
                            }
                        }
                    }

                    if (bMatch)
                    {
                        i = qres->bindex.ilist_pos[itype][prev] + nia + 1;
                        break;
                    }

                } /* i-loop */

                if (!bMatch)
                {
                    /* No, we didn't find that set of atoms.
                     * Start search from the beginning. */
                    i = 0;
                }

            }
            else
            {
                i = 0; /* Start from the beginning for this next interaction. */
            }

            bMatch = FALSE;

            while((i < ilist[itype].nr) && !bMatch)
            {
                ft = i++; /* functype/parameter index */
	      
                bMatch = TRUE;
	      
                for (j=0; j<nia; j++)
                {
                    /* Translate atoms in ilist to residue local numbers */
                    ia = /* bDD ? */
/* 		    l2g[ilist[itype].iatoms[i+j]] - qres->atoms[0] : /\* local != global*\/ */
                        ilist[itype].iatoms[i+j] - qres->atoms[0];       /* local != global*/

                    if (ia != ra[j])
                    {
                        /* No match. On to the next interaction; */
                        bMatch = FALSE;
                        break;
                    }
                }

                if (!bMatch)
                {
                    /* Test the atoms in reverse order */

                    bMatch = TRUE;

                    for (j=0; j<nia; j++)
                    {
                        /* Translate atoms in ilist to residue local numbers */
                        ia = /* bDD ? */
/* 			l2g[ilist[itype].iatoms[i+j]] - qres->atoms[0] : /\* local != global*\/ */
                            ilist[itype].iatoms[i+j] - qres->atoms[0];       /* local != global*/

                        if (ia != ra[nia-j-1])
                        {
                            /* No match. On to the next interaction; */
                            bMatch = FALSE;
                            break;
                        }
                    }
                }

                if (bMatch)
                { 
                    /* Ok, we've found it */
                    break;
                }

                i += (nia); /* proceed with next interaction */

            } /* while-loop over the ilist */

            if (bMatch)
            {
                /* interaction found. */
                qres->bindex.ilist_pos[itype][b] = (t_iatom)(i-1); /* Must be -1 since we want the parameter
                                                                    * index that preceeds the atoms too. */
                /* Should one store ft somewhere? No it never changes. */
                qres->bindex.indexed[itype][b] = TRUE;
                qres->nr_indexed++;
            }
            else
            {
                /* Constraints may be absent from the ilists */
#ifdef DEBUG
                if (itype != F_CONSTR)
                {
                    buf[0] = '\0';
                    for (j=0; j<nia; j++)
                    {
                        sprintf(b2," %d%c", ra[j], (j==nia-1) ? '.' : ',');
                        strncat(buf,b2,STRLEN-1-strlen(buf));
                    }

                    gmx_fatal(FARGS, "Interaction of type %s not found in ilist.\nResidue-local atom indices are %s", interaction_function[itype].name, buf);
                }
#endif
            }

        } /* bonded (b) loop */

    } /* bonded type (bt) loop */

}

/* returns the index in the ilib where bond to this proton is found. */
int qhop_get_proton_bond_params(const qhop_db *db, const titration_t T,
                                t_qhop_atom *qatom, gmx_localtop_t *top, 
                                int proton_id, const t_commrec *cr)
{
    int b, bi, i, ai;
    const int niatoms = 3;
    gmx_bool bMatch;
    t_restp *rtp;
    t_qhop_residue *qres;

    qres = &(T->qhop_residues[qatom->qres_id]);

    rtp  = &(db->rtp[db->rb.qrt[qres->rtype].subres[qres->subres].irtp]);

    /* find the parameters for this constraint.
     * We assume that there's only one constraint
     * involving this hydrogen. This is NOT the case with
     * angle constraints. */
    for (b=0; b < rtp->rb[ebtsBONDS].nb; b++)
    {
        /* bi will index the list of bonded interactions of qres->rtype*/
        bi = db->rb.qrt[qres->rtype].subres[qres->subres].biMap[ebtsBONDS][b];

        for (i=0; i < (niatoms-1); i++)
        {
            /* Assume global atom numbering. CHANGE!*/
            ai = proton_id - qres->atoms[0]; /* residue local atom id of the proton. */
            if (db->rb.qrt[qres->rtype].ba[ebtsBONDS][bi][i] == ai)
            {
                return db->rb.qrt[qres->rtype].subres[qres->subres].findex[ebtsBONDS][bi];
            }
        }
    }
    gmx_fatal(FARGS, "Constraint not found!");
  
    return 0;
}

void qhop_constrain(t_qhop_residue *qres, titration_t T, const qhop_db *db, gmx_localtop_t *top, t_mdatoms *md, int proton_id, gmx_constr_t constr, const t_inputrec *ir, const t_commrec *cr)
{
    int nr, nalloc, heavy, h, rt, r, b, bi, i, params, ai, *g2l, *iatoms;
    const int niatoms = 3;
    gmx_bool bMatch;
    char *heavyname, *hname;
    t_qhop_atom *qatom;
    t_restp *rtp;

    bMatch = FALSE;

    if (qres->nr_indexed == 0)
    {
        index_ilists(qres, db, top, cr);
    }
#if 0
    for (heavy=0; (heavy < qres->nr_titrating_sites) && (!bMatch); heavy++)
    {
        qatom = &(T->qhop_atoms[qres->titrating_sites[heavy]]);
        for (h=0; h < qatom->nr_protons; h++)
        {
            if (qatom->protons[h] == proton_id)
            {
                heavyname = qatom->atomname;
                heavy = qatom->atom_id;
                bMatch = TRUE;
                break;
            }
        }
        if (bMatch)
        {
            break;
        }
    }

    if (!bMatch)
    {
        /* doing nothing. The proton wasn't found on this residue. */
        return;
    }

    params =
        qhop_get_proton_bond_params(db, qr, qatom, top, proton_id, cr)  /* position in ilib */
        + top->idef.ntypes - db->rb.ni;                                 /* plus the offset */
  
    /* if (DOMAINDECOMP(cr)) */
/*     { */
/*       g2l = cr->dd->gatindex; /\* use this to map global to local atomnumbers below *\/ */
/*     } */

    /* Make a new entry at the end of the constraint ilist. */
    nr = top->idef.il[F_CONSTR].nr;
    nalloc = top->idef.il[F_CONSTR].nalloc;
  
    if (nr+niatoms > nalloc)
    {
        srenew(top->idef.il[F_CONSTR].iatoms, nalloc+niatoms);
        top->idef.il[F_CONSTR].nalloc += niatoms;
    }

    top->idef.il[F_CONSTR].iatoms[nr]   = params;
    top->idef.il[F_CONSTR].iatoms[nr+1] = /* DOMAINDECOMP(cr) ? g2l[heavy] : */ heavy;
    top->idef.il[F_CONSTR].iatoms[nr+2] = /* DOMAINDECOMP(cr) ? g2l[h]     : */ proton_id;
    top->idef.il[F_CONSTR].nr += 3;

    /* Update bonded index */

    h = proton_id - qres->atoms[0]; /* Make residue local index */

    /* Now find the bonded term we're adding */
    for (b=0; b < qres->bindex.nr[F_CONSTR]; b++)
    {
        bi = db->rb.res[qres->rtype][qres->res].biMap[ebtsBONDS][b]; /* Map to restype bondeds */

        rtp = &(db->rtp[db->rb.res[qres->rtype][qres->res].irtp]);

        iatoms = db->rb.ba[qres->rtype][ebtsBONDS][bi];

        if (h == iatoms[0] || h == iatoms[1])
        {
            qres->bindex.ilist_pos[F_CONSTR][b] = nr;
            break;
        }
    }
#endif
    /* Now hack the t_gmx_constr */
    set_constraints(constr, top,(t_inputrec *)ir, md,(t_commrec *)cr);

    /* gmx_fatal(FARGS, "Could not find the constraint in the rtp data."); */
}

void qhop_deconstrain(t_qhop_residue *qres, const qhop_db *db, gmx_localtop_t *top, t_mdatoms *md, int proton_id, gmx_constr_t constr, const t_inputrec *ir, const t_commrec *cr)
{
    int b, ip, *iatoms;
    const int niatoms = 3;
#if 0
    if (qres->nr_indexed == 0)
    {
        index_ilists(qres, db, top, cr);
    }

    for (b=0; b < qres->bindex.nr[F_CONSTR]; b++)
    {
        ip = qres->bindex.ilist_pos[F_CONSTR][b];

        if (ip == -1)
        {
            /* Can't remove this constraint, it's not really indexed
             * which probably means it's not there. */
            return;
        }

        iatoms = &(top->idef.il[F_CONSTR].iatoms[ip]);

        if (iatoms[1] == proton_id || iatoms[2] == proton_id)
        {
            /* Here it is. Remove it from the list. */

            /* Here the ilist is shifted by niatoms to erase the constraint in question. */
            memmove((void*)&(iatoms[0]),
                    (void*)&(iatoms[niatoms]),
                    (top->idef.il[F_CONSTR].nr-ip-niatoms) * sizeof(int));

            top->idef.il[F_CONSTR].nr -= niatoms;

            qres->bindex.ilist_pos[F_CONSTR][b] = -1;

            break;
        }
    }
#endif
    /* Now hack the t_gmx_constr */
    set_constraints(constr, top, (t_inputrec *)ir, md,(t_commrec *)cr);
}


/* Reads the info in the db->qhop_H_exist and finds the corresponding
 * residue subtype in the qhop_db.
 * Returns the index in db->res[rt][]
 * This function is used to set the correct protonation states
 * at the start of the simulation.
 * Right now matching is based on counting the protons on
 * titrating sites in the qhop_res and in the rtp. Thus,
 * the existence map may need slight reshuffling to really
 * represent the global protonation state. */
int which_subRes(const titration_t T,
                 qhop_db *db, const int resnr)
{
    int
        h, i, j, k, t, nH, *nHres, a, r, b, aa,
        atomid, nsubres, natoms, qrtp;
    char
        *Hname, *DAname;
    t_atom
        atom;
    t_qhop_residue
        *qres;
    gmx_bool
        bSameRes, bNew, bMatch;
    char
        **DAlist, **DAlist2;
    int
        *n_H, *n_H2, nHtot, nDA, nDA2, DA;
    qhop_reactant
        *reac;
    t_restp
        *rtp;

    range_check(resnr,0,T->nr_qhop_residues);
    qres    = &(T->qhop_residues[resnr]);

    range_check(qres->rtype,0,db->rb.nrestypes);
    nsubres = db->rb.qrt[qres->rtype].nsubres;


    /* To better handle the tautomerism, maybe one should
       count hydrogens on the donors/acceptors instead
       of matching the names. Here we match the names.
    */

    nDA = qres->nr_titrating_sites;
    snew(n_H, nDA);
    snew(DAlist, nDA);

    /* ************************************************
       Count the protons present on each titrating site
       (qhop_atom) for this qhop_residue.
    */

    /* loop over donors/acceptors in the residue */
    nHtot = 0;
    for (i=0; i < nDA; i++)
    {
        k = qres->titrating_sites[i];

        /* Stash the donor/acceptor name */
        DAlist[i] = T->qhop_atoms[k].atomname;

        /* Proton roll call */
        for (j=0; j < T->qhop_atoms[k].nr_protons; j++)
        {
            /* get the atom id of the protons */
            h = T->qhop_atoms[k].protons[j];

            /* is it present? */
            if (get_proton_presence(&(db->H_map), h))
            {
                n_H[i]++;
                nHtot++;
            }
        }
    }

    /* Ok. n_H now stores the number of active protons on all donors/acceptors.
       Donor/acceptor names are stored in DAlist. */


    /* ********************************************************
       Count the protons on each titrating site for different
       residue subtypes and see if any matches the qhop_residue
    */

    /* Loop over residue subtypes */
    bSameRes = FALSE;
    for (r=0; (r<nsubres) && !bSameRes; r++)
    {
        n_H2 = NULL;
        DAlist2 = NULL;
        nDA2 = 0;

        bMatch = TRUE;

        /* DA==0 => acceptors, DA==1 => donors */
        for (DA=0; DA<2 && bMatch; DA++)
        {
            /* Go through the donors and acceptors and make a list */
            for (j=0;
                 (j < (DA==0 ?
                       db->rb.qrt[qres->rtype].subres[r].na :
                       db->rb.qrt[qres->rtype].subres[r].nd))
                     && bMatch;
                 j++)
            {
                reac = (DA==0 ?
                        &(db->rb.qrt[qres->rtype].subres[r].acc[j]) :
                        &(db->rb.qrt[qres->rtype].subres[r].don[j]));

                /* Tautomer loop. There may be chemically equivalent
                 * sites, e.g. the oxygens in a carboxylate. */
                for (t=0; t<reac->nname && bMatch; t++)
                {
                    /* Have we already seen this donor/acceptor? */
                    bNew = TRUE;
                    for (k=0; (k < nDA2) && bNew; k++)
                    {
                        if (strcmp(DAlist2[k], reac->name[t]) == 0)
                        {
                            /* We already had it in the list */
                            bNew = FALSE;
                        }
                    }
		  
                    if (bNew)
                    {
                        srenew(DAlist2, nDA2+1);
                        srenew(n_H2, nDA2+1);

                        n_H2[nDA2] = 0;
                        DAlist2[nDA2] = reac->name[t];
		      
                        rtp = &(db->rtp[db->rb.qrt[qres->rtype].subres[r].irtp]);
                        for (a=0; a<rtp->natom; a++)
                        {
                            if (strcmp(reac->name[t], *(rtp->atomname[a])) == 0)
                            {
                                /* Here's the titrating site. Find hydrogens.
                                 * Assume a contiguous strtch of hydrogens followed
                                 * by the next heavy atom. */
                                for (aa=a+1; aa<rtp->natom; aa++)
                                {
                                    if (*(rtp->atomname[aa]) && *(rtp->atomname[aa][0]) == 'H')
                                    {
                                        /* A(nother) hydrogen. */
                                        n_H2[nDA2]++;
                                    }
                                    else
                                    {
                                        /* Next heavy atom. */
                                        break;
                                    }
                                }
                                break;
                            }
                        }

                        nDA2++;
                    }
                }
            }
        }

        /* Now compare the hydrogen counts to see if this is the flavor we want. */
        /* First check if the total number of hydrogens are the same for both lists */
        nH = 0;
        for (j=0; j<nDA; j++)
        {
            nH += n_H[j];
        }

        for (j=0; j<nDA2; j++)
        {
            nH -= n_H2[j];
        }
      
        bSameRes = (nH==0); /* Same total number of hydrogens? */

        /* Sum up hydrogens on tautomeric sites and compare */

        /* DA==0 => acceptors, DA==1 => donors */
        for (DA=0; DA<2 && bSameRes; DA++)
        {
            for (j=0;
                 (j < (DA==0 ?
                       db->rb.qrt[qres->rtype].subres[r].na :
                       db->rb.qrt[qres->rtype].subres[r].nd))
                     && bSameRes;
                 j++)
            {
                reac = (DA==0 ?
                        &(db->rb.qrt[qres->rtype].subres[r].acc[j]) :
                        &(db->rb.qrt[qres->rtype].subres[r].don[j]));

                nH = 0;
	      
                for (t=0; t<reac->nname && bSameRes; t++)
                {
                    for (i=0; i<nDA; i++)
                    {
                        if (strcmp(DAlist[i], reac->name[t]) == 0)
                        {
                            nH += n_H[i];
                            break;
                        }
                    }

                    for (i=0; i<nDA2; i++)
                    {
                        if (strcmp(DAlist2[i], reac->name[t]) == 0)
                        {
                            nH -= n_H2[i];
                            break;
                        }
                    }

                    /* If they are the same, then nH should be zero here. */
                    bSameRes = (nH==0);
                }
            }
        }

        if (bSameRes)
        {
            /* Adjust the existence map */
	  
            /* Zero all hydrogens */
            for (i=0; i < qres->nr_titrating_sites; i++)
            {
                for (j=0; j < T->qhop_atoms[qres->titrating_sites[i]].nr_protons; j++)
                {
                    set_proton_presence(&(db->H_map),
                                        T->qhop_atoms[qres->titrating_sites[i]].protons[j],
                                        FALSE);
                }
            }
	  
            /* Go through the rtp */
            rtp = &(db->rtp[db->rb.qrt[qres->rtype].subres[r].irtp]);
	  
            for (i=0; i < rtp->natom; i++)
            {
                /* Is it a hydrogen atom? */
                if (rtp->atom[i].atomnumber == 1)
                {
                    /* Match with atoms in qres */
                    for (j=0; j < qres->nr_atoms; j++)
                    {
                        if (strcmp(*(rtp->atomname[i]), qres->atomnames[j]) == 0)
                        {
                            set_proton_presence(&(db->H_map),
                                                qres->atoms[j],
                                                TRUE);
                        }
                    }
                }
            }
        }

        /* Clean up */
        if (n_H2 != NULL)
        {
            sfree(n_H2);
        }
      
        if (DAlist2 != NULL)
        {
            sfree(DAlist2);
        }
      
        if (bSameRes)
        {
            break;
        }
    }

    if (n_H != NULL)
    {
        sfree(n_H);
    }
      
    if (r >= nsubres)
    {
        gmx_fatal(FARGS, "Didn't find the subres of %s with %d H, most likely due to inconsistencies in\nforce field files or incorrect hydrogen existence map.\nr = %d, nsubres = %d, resnr = %d",
                  db->rb.qrt[qres->rtype].canonical,nHtot,r,nsubres,resnr);
    }
    /* printf("subres %d\n",r);*/
    return r;
}

void qhop_swap_bondeds(t_qhop_residue *swapres,
                       qhop_subres *prod,
                       qhop_db *db,
                       gmx_localtop_t *top,
                       const t_commrec *cr)
{
    int i, bt, b, bp, p, it, offset;
    t_restp *rtp, *rtpr;
    qhop_subres *qres;


    if  (swapres->nr_indexed == 0)
    {
        index_ilists(swapres, db, top, cr);
    }

    qres = &(db->rb.qrt[swapres->rtype].subres[swapres->subres]);

    rtp  = &(db->rtp[db->rb.qrt[swapres->rtype].irtp]);
    rtpr = &(db->rtp[qres->irtp]);

    /* rb->ilib and top->idef.functype are offset by this much: */
    offset = top->idef.ntypes - db->rb.ni;
  
    /* printf("ntypes = %d, rb.ni = %d, offset = %d\n",top->idef.ntypes,db->rb.ni,offset); */
  
    for (bt=0; bt < ebtsNR; bt++)
    {
        it = db->rb.btype[bt]; /* Which ilist */

        for (b=0; b < rtp->rb[bt].nb; b++)
        {
            /* Which bonded is this in the subtype? */

            p = db->rb.inull[bt]; /* Point to the dummy parameter.
                                   * If this bond isn't present in the
                                   * biMap, then it should be a dummy interaction. */

            for (bp=0; bp < rtpr->rb[bt].nb; bp++)
            {
                if (qres->biMap[bt][bp] == b)
                {
                    p = prod->findex[bt][bp]; /* Parameter index */
                    break;
                }
            }

            /* Change bonded parameters */

            /* We need this check, since all constraints may not be
             * present in the ilists and may therefore be unindexed. */
            if (swapres->bindex.indexed[it][b])
            {
                /* printf("Changing functype from %d to %d\n",
                   top->idef.il[it].iatoms[swapres->bindex.ilist_pos[it][b]],
                   p+offset); */
                top->idef.il[it].iatoms[swapres->bindex.ilist_pos[it][b]] =
                    p + offset;
            }
        }
    }
  
}

/* We change vdv by changing atomtype. */
void qhop_swap_vdws(const t_qhop_residue *swapres,
                    const qhop_subres *prod,
                    t_mdatoms *md,
                    const qhop_db *db)
{

    int i, j;
    t_restp *rtp;

    rtp = &(db->rtp[prod->irtp]);

    /* For now, assume that all atomtypes are present, e.g. that grompp was run with -norenum */

    /* Just make it all inert... */
    for (i=0; i < swapres->nr_atoms; i++)
    {
        md->typeA[swapres->atoms[i]] = db->inertH;
    }

    /* ...then set the atomtypes for the atoms in the residue */
    for (i=0; i < prod->niatom; i++)
    {
        md->typeA[swapres->atoms[prod->iatomMap[i]]] = rtp->atom[i].type;
    }
}

static void low_level_swap_m_and_q(t_mdatoms *md, const t_atom *atom, const int atomid)
{
    real m, q;
  
    if (atom == NULL)
    {
        m = 0;
        q = 0;
    }
    else
    {
        m = atom->m;
        q = atom->q;
    }

    md->chargeA[atomid] = q;

    if (md->massA != NULL)
    {
        md->massA[atomid] = m;
    }
    md->massT[atomid] = m;
    md->invmass[atomid] = (m==0 ? 0 : 1/m);
  
}

/* Hacks mdatoms such that the masses and
 * charges of swapres match those of prod.
 * Works on an atomname basis
 */
void qhop_swap_m_and_q(const t_qhop_residue *swapres,
                       const qhop_subres *prod,
                       t_mdatoms *md,
                       const qhop_db *db, titration_t T)
{
    int i, j, k, nH;
    t_restp *rtp;
    real m;
    gmx_bool bWater;
    t_qhop_atom *qa;

    bWater = db->rb.qrt[swapres->rtype].bWater;
  
    rtp = &(db->rtp[prod->irtp]);

    /* Compare atomnames with atomnames in rtp.
     * Change m and q upon match. */

    for (i=0; i < swapres->nr_atoms; i++)
    {
        /* Zero m and q for all atoms
         * Otherwise any atoms that was
         * just turned off will retain
         * their m and q. */

        low_level_swap_m_and_q(md, NULL, swapres->atoms[i]);

        for (j=0; j < rtp->natom; j++)
        {
            if (strcmp(*(rtp->atomname[j]), swapres->atomnames[i]) == 0)
            {
                /* Now set the q and m for real */
                low_level_swap_m_and_q(md, &(rtp->atom[j]), swapres->atoms[i]);
                break;
            }
        }
    }

    if (bWater)
    {
        /* Masses are different for water. */
        qa = &(T->qhop_atoms[swapres->titrating_sites[0]]);
        nH = 0;

        for (i=0; i<qa->nr_protons; i++)
        {
            if (get_proton_presence(&db->H_map, qa->protons[i]))
            {
                nH++;
            }
        }

        /* Hardcoded = ugly. Sorry. This will have to go. */
        md->massT[qa->protons[0]] = 1.008;
        md->massT[qa->protons[1]] = 1.008;
        md->massT[qa->protons[2]] = 0; /* vsite */
        md->massT[qa->protons[3]] = 0; /* vsite */
        md->massT[qa->protons[4]] = 0; /* vsite */
        md->massT[qa->atom_id]    = 15.9994 /*+ (nH-2)*1.008*/;
        for(k=0; (k<=4); k++) 
        {
            if (md->massA != NULL)
                md->massA[qa->protons[k]] = md->massT[qa->protons[k]];
            if (md->massT[qa->protons[k]] != 0)
                md->invmass[qa->protons[k]] = 1.0/md->massT[qa->protons[k]];
        }
        if (md->massA != NULL)
            md->massA[qa->atom_id] = md->massT[qa->atom_id];
        md->invmass[qa->atom_id]   = 1.0/md->massT[qa->atom_id];
    }
}

void set_interactions(titration_t T,
                      qhop_db *qdb,
                      gmx_localtop_t *top,
                      t_mdatoms *md,
                      t_qhop_residue *qres,
                      t_commrec *cr)
{
    int i, j, k, nri, bt, b, Gr, RB, R;
    qhop_subres *reac;
    t_restp *res;
    t_qhop_atom *QA;
  
    /* point pointers at the right residue */
    Gr   = qres->res_nr; /* global res id. */
    RB   = qres->rtype;
    R    = qres->subres;
    reac = &(qdb->rb.qrt[RB].subres[R]);
    QA = T->qhop_atoms;

    for (i=0; i<qres->nr_titrating_sites; i++)
    {
        /* find among qdb->res[][].acc/.don */
        for (j=0; j < reac->na; j++)
        {
            /* loop over tautomeric sites */
            for (k=0; k < reac->acc[j].nname; k++)
            {
                if (strcmp(reac->acc[j].name[k], QA[qres->titrating_sites[i]].atomname) == 0)
                {
                    QA[qres->titrating_sites[i]].state |= eQACC;
                    break;
                }
            }
        }

        /* find among qdb->res[][].acc/.don */
        for (j=0; j<reac->nd; j++)
        {
            /* loop over tautomeric sites */
            for (k=0; k<reac->don[j].nname; k++)
            {
                if (strcmp(reac->don[j].name[k], QA[qres->titrating_sites[i]].atomname) == 0)
                {
                    QA[qres->titrating_sites[i]].state |= eQDON;
                    break;
                }
            }
        }
    }
    qhop_swap_bondeds(qres, reac, qdb, top, cr);
    qhop_swap_vdws(qres, reac, md, qdb);
    qhop_swap_m_and_q(qres, reac, md, qdb, T);
    /* Non-bonded */
}

/* Sets the bqhopdonor[] and bqhopacceptor[] arrays in a t_mdatoms. */
static void qhop_atoms2md(t_mdatoms *md, const titration_t T)
{
    int i, j, k;
    t_qhop_atom *a;
    qhop_db *db;
    gmx_bool bWater;

    db = T->db;

    /* Should probably set the massT array too. */
    srenew(md->bqhopdonor,md->nalloc);
    srenew(md->bqhopacceptor,md->nalloc);
    for(k=0; (k<md->nalloc); k++) 
    {
        md->bqhopdonor[k] = FALSE;
        md->bqhopacceptor[k] = FALSE;
    }

    for (i=0; i < T->nr_qhop_atoms; i++)
    {
        a = &(T->qhop_atoms[i]);

        md->bqhopacceptor[a->atom_id] = (a->state & eQACC) != 0;
        md->bqhopdonor[a->atom_id]    = (a->state & eQDON) != 0;

        bWater = T->db->rb.qrt[T->qhop_residues[a->qres_id].rtype].bWater;

        for (j=0; (!bWater && (j < a->nr_protons)); j++)
        {
            if (db->H_map.H[db->H_map.atomid2H[a->protons[j]]] == 0)
            {
                md->massT[a->protons[j]] = 0;
            }
        }
    }
}

/* Sets the interactions according to the hydrogen existence map.
 * This requires a finalized t_mdatoms. */
static void finalize_titration_t(titration_t T, gmx_localtop_t *top, 
                                 t_mdatoms *md, t_commrec *cr)
{
    int i, j, nb, ft;
    qhop_db *db;
    t_qhop_residue *qres;
    /* #define DUMPTOP 1 */
#ifdef DUMPTOP 
    FILE *fp;
    if (NULL != (fp=ffopen("before.txt","w"))) 
    {
        pr_ltop(fp,0,"before",top,TRUE);
        fclose(fp);
    }
#endif
    db = T->db;

    for (i=0; i < T->nr_qhop_residues; i++)
    {
        qres = &(T->qhop_residues[i]);

        /* What flavour of the residue shall we set it to? */
        qres->subres = which_subRes(T, db, i);
      
        /* Alocate memory for the bonded index. */
        /* This should probably go to qhop_toputil.c */
        snew(qres->bindex.ilist_pos, F_NRE);
        snew(qres->bindex.indexed,   F_NRE);
        for (j=0; j < ebtsNR; j++)
        {
            ft = db->rb.btype[j];
	  
            if (ft >= 0)
            {
                nb = db->rtp[db->rb.qrt[qres->rtype].irtp].rb[j].nb;
                qres->bindex.nr[ft] = nb;

                snew((qres->bindex.ilist_pos[ft]), nb);

                snew((qres->bindex.indexed[ft]), nb);
            }
        }
        set_interactions(T, db, top, md, qres, cr);
    }
#ifdef DUMPTOP 
    if (NULL != (fp=ffopen("after.txt","w"))) 
    {
        pr_ltop(fp,0,"after",top,TRUE);
        fclose(fp);
    }
#endif

    qhop_atoms2md(md, T);
}


void finalize_titration(titration_t T,gmx_localtop_t *top,t_mdatoms *mdatoms,
                        t_commrec *cr)
{
    /* Complete the t_mdatoms and qhoprec, extend the topology. */
    make_ilib(T->db);
    qhop_attach_ilib(top, T->db);
    finalize_titration_t(T, top, mdatoms, cr);
}
