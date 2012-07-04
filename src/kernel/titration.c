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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "typedefs.h"
#include "mtop_util.h"
#include "random.h"
#include "gmx_random.h"
#include "constr.h"
#include "mdrun.h"
#include "typedefs.h"
#include "hackblock.h"
#include "lambertw.h"
#include "titration_db.h"
#include "titration.h"
#include "titration_util.h"
//#include "forcerec.h"

#define DOO   0.35
#define BOND  0.15
#define HBOND  0.2
#define MAXHOPANG 60
#define MAX_MC_TRANSFERS 4

typedef enum { eHOP_FORWARD, eHOP_BACKWARD } eHOP_t;

static void swap_rvec(rvec a,rvec b)
{
    rvec tmp;
    
    copy_rvec(a,tmp);
    copy_rvec(b,a);
    copy_rvec(tmp,b);
}

void fold_inactive_protons(const titration_t T, rvec x[], rvec v[])
{
    int r, t, h;
    t_qhop_residue *qres;
    t_qhop_atom *qatom;
    qhop_db *db;

    db = T->db;

    for (r=0; r < T->nr_qhop_residues; r++)
    {
        qres = &(T->qhop_residues[r]);
        if (!db->rb.qrt[qres->rtype].bWater)
        {
            for (t=0; t < qres->nr_titrating_sites; t++)
            {
                qatom = &(T->qhop_atoms[qres->titrating_sites[t]]);

                for (h=0; h < qatom->nr_protons; h++)
                {
                    if (!get_proton_presence(&(db->H_map), qatom->protons[h]))
                    {
                        copy_rvec(x[qatom->atom_id], x[qatom->protons[h]]);
                        copy_rvec(v[qatom->atom_id], v[qatom->protons[h]]);
                    }
                }
            }
        }
    }
}

/* Makes a list of tautomeric sites and returns its length.
 * ai is the qhop_atom id of the primary site. */
static int get_tautomeric_sites(const qhop_db *db, const titration_t T, int ai, int **tautomers, int DA)
{
    int i, nt, *t, rt, t_nr, r, reac, a;
    t_qhop_residue *qres;
    t_qhop_atom *qatom;
    qhop_reactant *R;
    /* t_restp *rtp; */

    t = NULL;
    nt = 0;
    qatom = &(T->qhop_atoms[ai]);
    qres  = &(T->qhop_residues[qatom->qres_id]);
    rt    = qres->rtype;
    r     = qres->subres;
  
    snew(t, ++nt);
    t[0] = ai;

    if (qres->nr_titrating_sites > 1)
    {
/*       rtp = &(db->rtp[db->rb.res[rt][r].rtp]); */

        for (reac=0;
             reac < (DA==eQACC) ? db->rb.qrt[rt].subres[r].na : 
                 db->rb.qrt[rt].subres[r].nd;
             reac++)
        {
            R = (DA==eQACC) ?
                &(db->rb.qrt[rt].subres[r].acc[reac]) :
                &(db->rb.qrt[rt].subres[r].don[reac]);
            if (strcmp(qatom->atomname, R->name[0]) == 0)
            {
                /* This is the primary. Collect the tautomers. */
                for (a=1; a < R->nname; a++)
                {
                    srenew(t, nt+1);
                    memset(&(t[nt]),0,sizeof(t[nt]));
                    for (t_nr=0; t_nr < qres->nr_titrating_sites; t_nr++)
                    {
                        if (strcmp(T->qhop_atoms[qres->titrating_sites[t_nr]].atomname, R->name[a]) == 0)
                        {
                            t[nt++] = qres->titrating_sites[t_nr];
                            break;
                        }
                    }
                }
	      
                break;
            }
        }
    }
  
    *tautomers = t;
    return nt;
}

/* *secname holds the name of the acceptor/donor. If this is indeed
 * the primary titratable site, then this function returns the
 * qhop_atom id of that atom. If not, it returns the qhop_atom id
 * of the primary site. */
static int qhop_get_primary(titration_t T, int hop_id, int DA)
{
    qhop_db *db;
    char *name;
    int rt, r, a, reac, t;
    qhop_subres *qres;
    qhop_reactant *qreac;
    t_qhop_residue *qhopres = &(T->qhop_residues[T->qhop_atoms[hop_id].qres_id]);
    t_qhop_atom *qatom = &(T->qhop_atoms[hop_id]);
    char *secname = T->qhop_atoms[hop_id].atomname;
    
    if (DA != eQDON && DA != eQACC)
    {
        gmx_fatal(FARGS, "Called qhop_get_primary() with DA != eQDON || eQACC");
    }

    db = T->db;
    rt = qhopres->rtype;
    r  = qhopres->subres;
    qres = &(db->rb.qrt[rt].subres[r]);

    name = qatom->atomname;
  
    for (reac=0;
         reac < (DA == eQACC) ? qres->na : qres->nd;
         reac++)
    {
        qreac = (DA == eQACC) ?
            &(qres->acc[reac]) :
            &(qres->don[reac]);

        for (t=0; t < qreac->nname; t++)
        {
            if (strcmp(qreac->name[t], name) == 0)
            {
                for (a=0; a < qhopres->nr_titrating_sites; a++)
                {
                    if (strcmp(qreac->name[0], T->qhop_atoms[qhopres->titrating_sites[a]].atomname) == 0)
                    {
                        return qhopres->titrating_sites[a];
                    }
                }
            }
        }
    }

    /* If we've gotten this far, then something is wrong. */
    gmx_fatal(FARGS,
              "Primary tautomeric site could not be determined:\n"
              "  %s atom id      %i (%s).",
              DA == eQDON ? "donor" : "acceptor",
              T->qhop_atoms[hop_id].atom_id,
              T->qhop_atoms[hop_id].atomname);
  
    return -1;
}

static qhop_parameters *get_qhop_params (titration_t T,t_hop *hop, qhop_db_t db)
{/* (char *donor_name,char *acceptor_name, qhop_db_t db){ */

    /* parameters of the current donor acceptor combination are based on 
     * the residue name.
     */

    char *donor_name, *acceptor_name,
        *donor_atom, *acceptor_atom;
  
    qhop_parameters
        *q;

    snew(q,1);

    if (db == NULL)
    {
        gmx_fatal(FARGS, "Qhop database not initialized.");
    }
    else
    {
        donor_name    = db->rb.qrt
            [T->qhop_residues[T->qhop_atoms[hop->donor_id].qres_id].rtype].subres
            [T->qhop_residues[T->qhop_atoms[hop->donor_id].qres_id].subres].name;

        acceptor_name = db->rb.qrt
            [T->qhop_residues[T->qhop_atoms[hop->acceptor_id].qres_id].rtype].subres
            [T->qhop_residues[T->qhop_atoms[hop->acceptor_id].qres_id].subres].name;

        donor_atom    = T->qhop_atoms[hop->primary_d].atomname;
        acceptor_atom = T->qhop_atoms[hop->primary_a].atomname;

        if (NULL != debug)
            fprintf(debug, "--- PARAMETERS FOR DONOR %s (%s) - ACCEPTOR %s (%s) ---\n",
                    donor_name, donor_atom, acceptor_name, acceptor_atom);

        if (qhop_db_get_parameters(db, donor_name, acceptor_name,
                                   donor_atom, acceptor_atom, q) != 1)
            gmx_fatal(FARGS, "Parameters not found in qhop database for donor %s-%s acceptor %s-%s.",donor_name,donor_atom,
                      acceptor_name,acceptor_atom);
    }

    qhop_db_print(q);
  
    return (q);
}


static int ilist_iterator_start(const t_ilist *ilist, const int itype, int* index, t_iatom *atomlist)
{
    if (ilist->nr == 0)
    {
        *index = 0;
        atomlist = NULL;
        return 0;
    }

    *index = 0;
    atomlist = &(ilist->iatoms[0]);
    return interaction_function[itype].nratoms;
}

/* Returns the length of the interaction in atoms, including the type. */
static int ilist_iterator(const t_ilist *ilist, const int itype, int *index, t_iatom *atomlist)
{
    int ilen;
    if (*index < ilist->nr)
    {
        return 0;
    }
  
    /* Go to next interaction */
    ilen = interaction_function[itype].nratoms;
    *index += ilen + 1;         /* +1 for the type. */
    atomlist = &(ilist->iatoms[*index]);

    return ilen;
}

/* Is atom id among the qhop_atoms?
 * Returns the matching qhop_atom,
 * or -1 if not found.*/
static int is_qhopatom(const titration_t T, const int id)
{
    int i;
    for (i=0; i<T->nr_qhop_atoms; i++)
    {
        if (T->qhop_atoms[i].atom_id == id)
        {
            return i;
        }
    }

    return -1;
}

/* returns the restype for residue with name rname */
static int is_restype(const qhop_db *db, const titration_t T, const char *rname)
{
    int i;

    if (rname == NULL || rname[0]=='\0')
    {
        return -1;
    }

    for (i=0; i<T->nr_qhop_residues; i++)
    {
        if (strcmp(rname, db->rb.qrt[T->qhop_residues[i].rtype].canonical) == 0)
        {
            return i;
        }
    }

    return -2;
}

/* Attaches proton atom_id to t_qhop_atom qatom if it's not already attached. */
static void attach_proton(t_qhop_atom *qatom, const int id)
{
    int i;

    for (i=0; i<qatom->nr_protons; i++)
    {
        if (qatom->protons[i] == id)
        {
            return;
        }
    }

    srenew(qatom->protons, qatom->nr_protons+1);
    qatom->protons[qatom->nr_protons++] = id;
}

/* Finds the hydrogens connected to the t_qhop_atoms in qr. */
static void scan_for_H(titration_t T, const qhop_db *db)
{

#ifdef COMPLICATED_SCAN
    /*REDUNDANT (I hope)*/
    int i, j, a,  natoms, tsite, a1, a2, a3, qa;
    t_iatom *atom;

    natoms = ilist_iterator_start(ilist, itype, &i, atom);
   
    if (natoms == 0)
    {
        return;
    }

    /* Uses global atom numbers. Make this work with node-local numbering! */

    do
    {
        /* Remember that the interaction's first element is the interaction type,
         * so the first atom is ilist[1] */

        switch(itype)	
        {
        case F_BONDS:
        case F_G96BONDS:
        case F_MORSE:
        case F_CUBICBONDS:
        case F_CONSTR:
        case F_HARMONIC:
        case F_FENEBONDS:
        case F_TABBONDS:
            for(a=1; a<=2; a++)
            {
                a1 = atom[a];
                if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1) /* Make sure atomnumber is set! */
                {
                    /* is it connected to a t_qhop_atom? */
                    for (i=0; i<T->nr_qhop_residues; i++)
                    {
                        if (T->qhop_residues[i].rtype = rt)
                        {
                            for (tsite = 0;
                                 tsite<T->qhop_residues[i].nr_titrating_sites;
                                 tsite++)
                            {
                                a2 = atom[a%2 + 1];
                                qa = T->qhop_residues[i].titrating_sites[tsite];
                                if (strcmp(T->qhop_atoms[qa].atomname,
                                           *(mtop->moltype[mt].atoms.atomname[a2])) == 0)
                                {
                                    attach_proton(&(T->qhop_atoms[qa]), a2);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            break;
	  
        case F_SETTLE:
            /* This implies an oxygen and two attached hydrogens */
            a1 = atom[1]; /* O */
            a2 = atom[2]; /* H */
            a3 = atom[3]; /* H */

            for (i=0; i<T->nr_qhop_residues; i++)
            {
                if (T->qhop_residues[i].rtype = rt)
                {
                    for (tsite = 0;
                         tsite<T->qhop_residues[i].nr_titrating_sites;
                         tsite++)
                    {
                        qa = T->qhop_residues[i].titrating_sites[tsite];
                        if (strcmp(T->qhop_atoms[qa].atomname,
                                   *(mtop->moltype[mt].atoms.atomname[a1])) == 0)
                        {
                            attach_proton(&(T->qhop_atoms[qa]), a2);
                            attach_proton(&(T->qhop_atoms[qa]), a3);
                        }
                    }
                }
            }
            break;

        case F_VSITE2:
            /* Can't see when this vsite type would be applicable,
             * but let's implement it anyway. */
            a1 = atom[1]; /* H */
            if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1)
            {
                /* is it connected to a t_qhop_atom? */
                for (i=0; i<T->nr_qhop_residues; i++)
                {
                    if (T->qhop_residues[i].rtype = rt)
                    {
                        for (tsite = 0;
                             tsite<T->qhop_residues[i].nr_titrating_sites;
                             tsite++)
                        {
                            for (a=2; a<=3; a++)
                            {
                                a2 = atom[a];
                                qa = T->qhop_residues[i].titrating_sites[tsite];
                                if (strcmp(T->qhop_atoms[qa].atomname,
                                           *(mtop->moltype[mt].atoms.atomname[a2])) == 0)
                                {
                                    attach_proton(&(T->qhop_atoms[qa]), a2);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            break;

        case F_VSITE3:
        case F_VSITE3FD:
        case F_VSITE3FAD:
        case F_VSITE3OUT: /* <- Not entirely sure this works with an amine */
        case F_VSITE4FD:
            a1 = atom[1]; /* H */
            if (mtop->moltype[mt].atoms.atom[a1].atomnumber = 1)
            {
                /* is it connected to a t_qhop_atom? */
                for (i=0; i<T->nr_qhop_residues; i++)
                {
                    if (T->qhop_residues[i].rtype = rt)
                    {
                        for (tsite = 0;
                             tsite<T->qhop_residues[i].nr_titrating_sites;
                             tsite++)
                        {
                            a2 = atom[2];
                            qa = T->qhop_residues[i].titrating_sites[tsite];
                        }
                    }
                }
            }

            if (db->H_map.atomid2H[atom[1]] >= 0)
            {
                /* Second atom is attached to the H. */
                tsite = is_qhopatom(qr, atom[2]);
                if (tsite != -1)
                {
                    attach_proton(&(T->qhop_atoms[tsite]), atom[1]);
                    break;
                }
            }
            break;

        default:
            fprintf(stderr, "Unsupported connection type for qhop.\n");
        }

    }
    while(natoms = ilist_iterator(ilist, itype, &i, atom) != 0);

#else
    /* Here's the simple scan */

    gmx_bool bNext;
    int  a, aa, r, qa;
    t_qhop_atom *qatoms;
    t_qhop_residue *qresidues;

    qatoms = T->qhop_atoms;
    qresidues = T->qhop_residues;

    for (qa=0; qa < T->nr_qhop_atoms; qa++)
    {
        bNext = FALSE;
        r = qatoms[qa].qres_id;

        for (a=0; a < qresidues[r].nr_atoms && !bNext; a++)
        {
            if (qresidues[r].atoms[a] == qatoms[qa].atom_id)
            {
                /* Atart scanning for hydrogens. Attach all H found
                 * to the qhop_atom qa and break as soon as a heavy
                 * atom is encountered.
                 * This assumes that all hydrogens, and the hydrogens only,
                 * have names starting with 'H', and that hydrogens come
                 * immediately after the heavy atoms they're attached to. */

                for (aa = a+1; aa < qresidues[r].nr_atoms; aa++)
                {
                    if (qresidues[r].atomnames[aa] != NULL &&
                        qresidues[r].atomnames[aa][0] == 'H')
                    {
                        attach_proton(&(qatoms[qa]), qresidues[r].atoms[aa]);
                    }
                    else
                    {
                        bNext = TRUE;
                        break;
                    }
                }
            }
        }
    }
    
#endif 

}

/* This function now makes get_qhop_residue obsolete. */
/* remove commrec? */
static int get_qhop_atoms(FILE *fplog,
                          t_commrec *cr,
                          gmx_mtop_t *mtop,
                          titration_t T,
                          t_inputrec *ir,
                          qhop_db *qdb){

    /* reads in the qhop donor and acceptor atoms, finds their protons
       (the ones that can move), identifies their residue numbers and
       stores the information in the qhop_atoms struct.
    */

    /* Uses the qhop_resblocks in qdb as a basis for finding acceptors and donors. */
    int
        rb, r, a,
        reac, nreac[2], AD,
        **H, nH,
        q_atoms_nr,   q_atoms_max,
        q_residue_nr, q_residue_max,
        i, j, k, o, ngrps, resnr, mt;

    gmx_bool
        bMatchRB, bMatch;

    char
        *resname,
        *atomname;

    t_qhop_atom
        *q_atoms,
        *qa;

    t_qhop_residue          *q_residue;
    t_restp                 *rtp;
    qhop_reactant           *qreac[2];
    gmx_groups_t            *groups;
    gmx_mtop_atomloop_all_t aloop;
    t_atom                  *atom;
    t_ilist                 *ilist;
    
    q_atoms_nr = 0;
    q_residue_nr = 0;
    q_atoms_max = 1000;
    q_residue_max = 1000;

    snew(q_atoms, q_atoms_max);
    snew(q_residue, q_residue_max);
    
    groups = &mtop->groups;

    ngrps = ir->opts.ngTitrationH;

    if (NULL != fplog)
        fprintf(fplog, "ir->opts.ngTitrationH = %d\n", ir->opts.ngTitrationH);

    aloop = gmx_mtop_atomloop_all_init(mtop);

    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom)) 
    {
        /* ***************************
           Set the H-existence array */

        /* group loop */
        for(j=0; j<ngrps; j++)
        {
            if (ggrpnr(groups, egcTitrationH ,i) == j)
            {
                if (NOTSET == qdb->H_map.atomid2H[i])
                    gmx_fatal(FARGS,"Inconsistency in h_exist: atom %d is not a hydrogen",i+1);
                qdb->H_map.H[qdb->H_map.atomid2H[i]] = 1;
            }
        }

        /* *********************
           Find the Titration atoms */

        /* Do the resname and atomname show up in the qhop_resblocks? */

        if(q_atoms_nr >= q_atoms_max)
        {
            q_atoms_max += 1000;
            srenew(q_atoms,q_atoms_max);
            memset(&(q_atoms[q_atoms_max-1000]),0,1000*sizeof(q_atoms[0]));
        }

        if(q_residue_nr >= q_residue_max)
        {
            q_residue_max += 1000;
            srenew(q_residue,q_residue_max);
            memset(&(q_residue[q_residue_max-1000]),0,1000*sizeof(q_residue[0]));
        }

        gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

        bMatchRB = FALSE;
        bMatch = FALSE;
        T->global_atom_to_qhop_atom[i] = NOTSET;

        for (rb=0; (rb < qdb->rb.nrestypes) && !bMatchRB; rb++)
        {
            if (strcmp(qdb->rb.qrt[rb].canonical, resname) == 0)
            {
                /* We have a matching res(block) */
                bMatchRB = TRUE;

                for (r=0; (r < qdb->rb.qrt[rb].nsubres) && !bMatch; r++)
                {
                    /* Is the atom found among the donors/acceptors ?*/
                    nreac[0] = qdb->rb.qrt[rb].subres[r].na;
                    nreac[1] = qdb->rb.qrt[rb].subres[r].nd;
                    qreac[0] = qdb->rb.qrt[rb].subres[r].acc;
                    qreac[1] = qdb->rb.qrt[rb].subres[r].don;

                    /* AD==0 => acceptors, 1 => donors. This enables looping over both. */
                    for (AD=0;
                         !bMatch && AD<2;
                         AD++)		  
                    {
	      
                        for (reac = 0;
                             !bMatch && (reac < nreac[AD]);
                             reac++)
                        {
                            for (a=0;
                                 !bMatch && a<qreac[AD][reac].nname;
                                 a++)
                            {
                                if (strcmp(qreac[AD][reac].name[a], atomname) == 0)
                                {
                                    bMatch = TRUE;
                                    q_atoms[q_atoms_nr].resname      = qdb->rb.qrt[rb].canonical;
                                    q_atoms[q_atoms_nr].atomname     = qreac[AD][reac].name[a];
                                    q_atoms[q_atoms_nr].res_id       = resnr-1;
                                    q_atoms[q_atoms_nr].atom_id      = i;
                                    q_atoms[q_atoms_nr].nr_protons   = 0; /* To be decided */
                                    q_atoms[q_atoms_nr].protons      = NULL;
                                    q_atoms[q_atoms_nr].nr_acceptors = NOTSET;
                                    q_atoms[q_atoms_nr].acceptors    = NULL;
                                    T->global_atom_to_qhop_atom[i]   = q_atoms_nr;

                                    /* Short circuiting makes this conditional valid */
                                    if ((q_atoms_nr==0)
                                        || (q_atoms[q_atoms_nr].res_id !=
                                            q_atoms[q_atoms_nr-1].res_id))
                                    {
                                        /* New residue. */
                                        q_atoms[q_atoms_nr].qres_id                = q_residue_nr;
                                        q_residue[q_residue_nr].rtype              = rb;
                                        q_residue[q_residue_nr].subres             = NOTSET;     /* To be decided. */
                                        q_residue[q_residue_nr].atoms              = NULL;     /* To be set. */
                                        q_residue[q_residue_nr].atomnames          = NULL;     /* To be set. */
                                        q_residue[q_residue_nr].nr_atoms           = 0;     /* To be decided. */
                                        q_residue[q_residue_nr].nr_titrating_sites = 0; /* To be decided. */
                                        q_residue[q_residue_nr].res_nr             = NOTSET;  /* To be set */
                                        q_residue[q_residue_nr].res_nr             = resnr-1;

                                        q_residue_nr++;
                                    }
                                    else 
                                    {
                                        q_atoms[q_atoms_nr].qres_id = q_atoms[q_atoms_nr-1].qres_id;
                                    }

                                    if (q_residue_nr > 0)
                                    {
                                        srenew(q_residue[q_residue_nr-1].titrating_sites,
                                               q_residue[q_residue_nr-1].nr_titrating_sites + 1);

                                        q_residue[q_residue_nr-1].titrating_sites[q_residue[q_residue_nr-1].nr_titrating_sites] = q_atoms_nr;
                                        q_residue[q_residue_nr-1].nr_titrating_sites++;
                                    }
                                    else
                                    {
                                        gmx_fatal(FARGS, "q_residue_nr = %i, which shouldn't be able to happen!\n", q_residue_nr);
                                    }

                                    q_atoms_nr++;		    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if (NULL != fplog)
        fprintf(fplog,"There are %d Titration donors\n", q_atoms_nr);
    if (q_atoms_nr == 0)
        gmx_fatal(FARGS,"Titration MD calculation requested but there are no proton donors.");
    /* Assign atoms to each qhop_residue */
    aloop = gmx_mtop_atomloop_all_init(mtop);

    while (gmx_mtop_atomloop_all_next(aloop, &i, &atom))
    {
        gmx_mtop_atomloop_all_names(aloop, &atomname, &resnr, &resname);

        for (r=0; r < q_residue_nr; r++)
        {
            if (resnr-1 == q_residue[r].res_nr)
            {
                rtp = &(qdb->rtp[q_residue[r].rtype]);
	      
                srenew(q_residue[r].atoms,     q_residue[r].nr_atoms + 1);
                srenew(q_residue[r].atomnames, q_residue[r].nr_atoms + 1);

                q_residue[r].atoms    [q_residue[r].nr_atoms] = i;
                q_residue[r].atomnames[q_residue[r].nr_atoms] = strdup(atomname);
                q_residue[r].nr_atoms++;

                /* The following block was made redundant by the call to scan_for_H in qhop_init() */

/* 	      /\* Is it a hydrogen? *\/ */
/* 	      if (atom->atomnumber == 1) */
/* 		{ */
/* 		  /\* Find the right titrating site and add the H to the proton array. *\/ */
/* 		  bMatch = FALSE; */

/* 		  for (j=0; j < q_residue[r].nr_titrating_sites && !bMatch; j++) */
/* 		    { */
/* 		      qa = &(q_atoms[q_residue[r].titrating_sites[j]]); */

/* 		      /\* Loop over the t_restp to see which donor/acceptor it is bound to *\/ */
/* 		      for (k=0; k < rtp->rb[ebtsBONDS].nb && !bMatch; k++) */
/* 			{ */
/* 			  for (a=0; a<2; a++) */
/* 			    { */
/* 			      if (strcmp(rtp->rb[ebtsBONDS].b[k].a[a%2], qa->atomname) == 0 && */
/* 				  strcmp(rtp->rb[ebtsBONDS].b[k].a[(a+1)%2], atomname) == 0) */
/* 				{ */
/* 				  srenew(qa->protons, qa->nr_protons+1); */
/* 				  qa->protons[qa->nr_protons] = i; */
/* 				  qa->nr_protons++; */

/* 				  bMatch = TRUE; */
/* 				  break; */
/* 				} */
/* 			    } */
/* 			} */
/* 		    } */
/* 		} */
            }
        }
    }  

    /* Clean up */

    srenew(q_atoms, q_atoms_nr);     /* size it down, as it is most likely oversized. */
    srenew(q_residue, q_residue_nr); /* size it down, as it is most likely oversized. */
    T->qhop_atoms       = q_atoms;
    T->qhop_residues    = q_residue;
    T->nr_qhop_atoms    = q_atoms_nr;
    T->nr_qhop_residues = q_residue_nr;


    /* Scan the idef for the hydrogens we may have missed. */

    /* Tried gmx_mtop_ilistloop_*, but I couldn't get it to work.*/


    /* loop over moltypes */
/*   for (mt=0; mt < mtop->nmoltype; mt++) */
/*     { */
/*       /\* loop over residues in molecule *\/ */
/*       for (r=0; r < mtop->moltype[mt].atoms.nres; r++) */
/* 	{ */
/* 	  rt = is_restype(qr, *(mtop->moltype[mt].atoms.resinfo[r].name)); */

/* 	  if (rt >= 0) */
/* 	    { */
/* 	      for(a=0; a < mtop->moltype[mt].atoms.nr; a++) */
/* 		{ */
/* 		  if (mtop->moltype[mt].atoms.atom[a].res == r) */
/* 		    { */
		     
/* 		    } */
/* 		} */
    /* Try the g_hbond way. */

/* 	      for (i=0; i<F_NRE; i++) */
/* 		{ */
/* 		  if (i == F_POSRES) */
/* 		    { */
/* 		      continue; */
/* 		    } */
		  
/* 		  for (j=0; */
/* 		       j < ilist[i]->nr; */
/* 		       j+= interaction_function[i].nratoms+1) */
/* 		    { */
/* 		      if (i == F_SETTLE) */
/* 			{ */
/* 			  a1 = ilist[i]->iatoms[j+1]; */
/* 			  a2 = ilist[i]->iatoms[j+2]; */
/* 			  a3 = ilist[i]->iatoms[j+3]; */
/* 			  for (k=0; k < T->nr_qhop_residues; k++) */
/* 			    { */
/* 			      if (T->qhop_residues[k] == rt) */
/* 				{ */
/* 				  for (a=0; a < T->qhop_residues[k].nr_titratable_sites; a++) */
/* 				    { */
/* 				      qa = T->qhop_residues[k].titratable_sites[a]; */
/* 				      if (strcmp(T->qhop_atoms[qa], */
/* 						 *(mtop->moltype[mt].t_atoms.atomname[a1])) == 0) */
/* 					{ */
/* 					  attach_proton(&(T->qhop_atoms[qa]), a2); */
/* 					  attach_proton(&(T->qhop_atoms[qa]), a3); */
/* 					} */
/* 				    } */
/* 				} */
/* 			    } */
/* 			  continue; */
/* 			} */
		      
/* 		      if (IS_CHEMBOND(i)) */
/* 			{ */
/* 			  for (o=0; o<2; o++) */
/* 			    { */
/* 			      a1 = ilist[i]->iatoms[j+1+o]; */
/* 			      a2 = ilist[i]->iatoms[j+1+((o+1)%2)]; */
/* 			      if (mtop->moltype[mt].t_atoms.atom[a1].atomnumber = 1) */
/* 				{ */
/* 				  for (k=0; k < T->nr_qhop_residues; k++) */
/* 				    { */
/* 				      if (T->qhop_residues[k] == rt) */
/* 					{ */
/* 					  for (a=0; a < T->qhop_residues[k].nr_titratable_sites; a++) */
/* 					    { */
/* 					      qa = T->qhop_residues[k].titratable_sites[a]; */
/* 					      if (strcmp(T->qhop_atoms[qa], */
/* 							 *(mtop->moltype[mt].t_atoms.atomname[a2])) == 0) */
/* 						{ */
/* 						  attach_proton(&(T->qhop_atoms[qa]), a1); */
/* 						} */
/* 					    } */
/* 					} */
/* 				    } */
/* 				  break; */
/* 				} */
/* 			    } */
/* 			  continue; */
/* 			} */

/* 		      if (IS_VSITE(i)) */
/* 			{ */
/* 			  a1 = ilist[i]->iatoms[j+1]; */
/* 			  /\* Now step backwards until a non-hydrogen is found.*\/ */
/* 			  a2 = a1-1; */
			  
/* 			} */
/* 		    } */
/* 		} */

/* 	      for (i=0; i<F_NRE; i++) */
/* 		{ */
/* 		  if (IS_CHEMBOND(i) */
/* 		      || i == F_SETTLE */
/* 		      || (i >= F_VSITE2 && i <= F_VSITE4FD)) */
/* 		    { */
/* 		      scan_ilist_for_hydrogens(mtop, qhoprec, &(qdb->H_map), rt, ilist[i], i); */
/* 		    } */
/* 		} */
/* 	    } */
/* 	} */
/*     } */
    return(q_atoms_nr);
} /* get_qhop_atoms */




/* obsolete */
static int int_comp(const void *a, const void *b){
    return (*(int *)a) - (*(int *)b);
}

/* obsolete */
static int qhop_atom_comp(const void *a, const void *b){
  
    return (int)(((t_qhop_atom *)a)->atom_id)-(int)(((t_qhop_atom *)b)->atom_id);
  
}

/* OBSOLETE! */
/* remove commrec? */
/* static void get_qhop_residue(t_commrec *cr, gmx_mtop_t *mtop,  */
/* 			     t_qhop_atom    *qhop_atom, */
/* 			     t_qhop_residue *qhop_residue){ */
  
/*   /\* from the qhopatoms struct, we now retreive the corresponding */
/*      residues. Of these residues we store the charges for each of the */
/*      states, the residue can be in. The latter are picked up from a */
/*      file somewhere. */
/*   *\/   */

/*  /\* We'll do this in an entirely new way. We no longer have a t_qhop_residue per t_qhop_atom, */
/*   * but one per qhopable residue. As some residues may have several titrating sites the numbers may differ. *\/ */

/* /\*   gmx_mtop_atomloop_all_t aloop; *\/ */
/* /\*   t_atom *\/ */
/* /\*     *atom; *\/ */
/*   int */
/*     maxatom = 100, nr_atoms=0,j,resnr; */
/*   char *resname, *atomname; */
//  snew(resname,6);
//snew(atomname,6);
/*   snew(qhop_residue->atoms,maxatom); */
/*   qhop_residue->res_nr    = qhop_atom->res_id; */
/*   qhop_residue->atoms[0] = qhop_atom->atom_id; /\* we quicksort later... *\/ */

/*   aloop = gmx_mtop_atomloop_all_init(mtop); */
/*   while (gmx_mtop_atomloop_all_next(aloop,&j,&atom)) { */
/*     gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
/*     if(resnr == qhop_atom->res_id){ */
/*       if(nr_atoms >= maxatom){ */
/* 	maxatom += 100; */
/* 	srenew(qhop_residue->atoms,maxatom); */
/*       } */
/*       qhop_residue->atoms[nr_atoms++] = j; */
/*     } */
/*   }     */
/*   qhop_residue->nr_atoms = nr_atoms; */
/*   srenew(qhop_residue->atoms,nr_atoms); */
/*   qsort(qhop_residue->atoms,nr_atoms, */
/* 	(size_t)sizeof(qhop_residue->atoms[0]),int_comp); */
//  free(atomname);
//free(resname);
/* } /\* get_qhop_residue *\/ */

/* Call this one to decide wether e.g. LYS or LYSH should be used,
 * based upon the hydrogen existence function.
 * This should only be required at t=0, because at other times we know
 * from and to what residue subtype a transfer is going. */

/* OBSOLETE! */
/* static void determine_protonation_state(t_qhoprec *qhoprec, qhop_db *qdb, gmx_mtop_t *top,  int resnr) */
/* { */
/*   int a, r=0, nres=0, i=0, nH, nHref; */
/*   int *H; */
/*   gmx_bool match; */
/*   t_qhop_residue *res=NULL; */
/*   qhop_res *qres; */
/*   gmx_mtop_atomloop_all_t aloop; */
/*   t_atom    *atom; */

/*   for (r=0; r<T->nr_qhop_residues; r++) */
/*     if (resnr = T->qhop_residues[r].res_nr) */
/*       { */
/* 	res = &(T->qhop_residues[r]); */
/* 	break; */
/*       } */
  
/*   if (res == NULL) */
/*     gmx_fatal(FARGS, "Couldn't find residue with global resnr %i in the qhoprec."); */

/*   /\* extract the hydrogens in res to a short list *\/ */

/*   snew(H, 10); */
/*   aloop = gmx_mtop_atomloop_all_init(mtop); */
/*   while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) */
/*     { */
/*       gmx_mtop_atomloop_all_names(aloop,&atomname,&resnr,&resname); */
/*       if (nH % 10 == 0) */
/* 	{ */
/* 	  nH+=10; */
/* 	  srenew(H, nH); */
/* 	} */
      
/*       for (a=0; a<res->nr_atoms; a++) */
/* 	{ */
	  
/* 	} */
/*     } */
/*   /\* Which restypes are available? *\/ */
/*   nres = qdb->rb.nres[res->rtype]; /\* How many to choose from? *\/ */
/*   match = FALSE; */
/*   for (r=0; r<nres && !match; r++) */
/*     { */
/*       /\* compare the names of the existing and non-existing hydrogens */
/*        * with those in the residue subtypes *\/ */
      
/*     } */
  

/*   /\* specify residue type from the H-existence array. *\/ */
/*   sfree(H); */
/* } */

static void get_protons(t_commrec *cr, gmx_mtop_t *mtop, 
                        t_qhop_atom    *qhop_atom,
                        t_qhop_residue *qhop_residue,
                        rvec *x, t_pbc *pbc){
    /* searching through the residue atoms to find the protons. By using
       a distance criteriumm we avoid the graph (don't know if there is
       such thing with v-sites, which we eventually want to use for
       qhop).
    */
    int
        i,nr_protons,max_protons;
    t_atom
        *atom;
    rvec
        vec;
    real
        dist;

    nr_protons = 0;
    max_protons = 10;

    snew(qhop_atom->protons,max_protons);

    for(i=0; i < qhop_residue->nr_atoms; i++)
    {
        gmx_mtop_atomnr_to_atom(mtop,qhop_residue->atoms[i],&atom);

        if (atom->atomnumber <= 1)
        {
            /* we need  v-sites as well */
            pbc_dx(pbc,
                   x[qhop_atom->atom_id],
                   x[qhop_residue->atoms[i]],
                   vec);

            dist = sqrt(iprod(vec,vec));

            if(dist < BOND )
            {
                if(nr_protons >= max_protons)
                {
                    max_protons += 10;
                    srenew(qhop_atom->protons,max_protons);
                }
                qhop_atom->protons[nr_protons++] = qhop_residue->atoms[i];
            }
        }
    }

    qhop_atom->nr_protons = nr_protons;
    //srenew(qhop_atom->protons, nr_protons);
} /* get_protons */ 


/* Related functions wrapped up in a new function
 * to de-clutter init_titration(). */
static void qhop_connect_rtp_library(FILE *fplog,qhop_db *db)
{
    /* Make tables of interacting atoms for the bonded rtp data
       so we don't have to go via atomnames all the time */    
    qhop_db_names2nrs(fplog,db);
    /* So we can superimpose the subresidues on the residue type */
    qhop_db_map_subres_atoms(db);   
    /* ditto. */
    qhop_db_map_subres_bondeds(db); 
}

static gmx_bool has_constraints(gmx_mtop_t *mtop)
{
    int i;
    
    for(i=0; (i<mtop->ffparams.ntypes); i++)
        if ((mtop->ffparams.functype[i] == F_CONSTR) || (mtop->ffparams.functype[i] == F_SETTLE))
            return TRUE;
    return FALSE;
}

void init_titration(FILE *fplog,const char *ff,
                    t_commrec *cr, gmx_mtop_t *mtop, t_inputrec *ir, 
                    t_forcerec *fr, matrix box, t_mdatoms *md)
{
    int 
        nr_qhop_atoms, nr_qhop_residues, i, j, target, ft, nb;
    t_qhop_residue
        *qres;
    qhop_db *db;
    titration_t T;
    
    fprintf(stderr, "Titration MD computation requested. May the force be with you....\n");

    snew(T,1);
    fr->titration = T;
    
    nr_qhop_atoms    = 0;
    nr_qhop_residues = 0;

    if ((db = qhop_db_read((char *)ff, mtop)) == NULL)
    {
        gmx_fatal(FARGS,"Can not read qhop database information");
    }
    
    snew(T->global_atom_to_qhop_atom, mtop->natoms);

    nr_qhop_atoms = get_qhop_atoms(fplog,cr, mtop, T ,ir, db);

    /* Find hydrogens and fill out the T->qhop_atoms[].protons */
    scan_for_H(T, db);

    qhop_connect_rtp_library(fplog,db);

    db->inertH = find_inert_atomtype(mtop, fr);
    if (NULL != fplog)
        fprintf(fplog,"Inert hydrogen for titration calculations has type %d\n",db->inertH);
    db->bConstraints = has_constraints(mtop);
    if (NULL != fplog)
        fprintf(fplog,"Constraints in the topology: %s\n",yesno_names[db->bConstraints]);
        
    T->hop = NULL;
    T->db = db;
} /* init_titration */


static void find_acceptors(FILE *fplog,t_commrec *cr, t_forcerec *fr, 
                           titration_t T,rvec *x, 
                           t_pbc *pbc, t_mdatoms *md, qhop_H_exist *Hmap)
{
    /* using the qhop nblists with i particles the donor atoms and
       j-particles the acceptor atoms, we search which of the acceptor
       atoms are available to accept a proton from the donor atoms. We
       make use of simple geometric criteria, such as distance and DHA
       angles. The acceptor atoms are stored in acceptor arrays of each
       of the donor atoms.
    */
    int         i,j,jmax=100,k,acc,don;
    t_nblist    qhoplist;
    t_qhop_atom donor,acceptor;
    rvec        dx,veca,vecb,vecc;
    int         nr_hop=0,closest;
    gmx_bool    found_proton;
    real        ang, ang2,r_close;
    
    qhoplist = fr->titration_nblist;
    
    for(i=0; i < qhoplist.nri; i++)
    {
        /* every qhop atom should be a unique i
           particle I think  */
        if(qhoplist.iinr[i] < 0 || !md->bqhopdonor[qhoplist.iinr[i]])
        {
            continue;
        }

        don      = T->global_atom_to_qhop_atom[qhoplist.iinr[i]];
        range_check(don,0,T->nr_qhop_atoms);
        donor    = T->qhop_atoms[don];

        for(j = qhoplist.jindex[i];
            j < qhoplist.jindex[i+1];
            j++)
        {
            if(!md->bqhopacceptor[qhoplist.jjnr[j]])
            {
                continue;
            }

            acc =  T->global_atom_to_qhop_atom[qhoplist.jjnr[j]];
            range_check(acc,0,T->nr_qhop_atoms);
            acceptor = T->qhop_atoms[acc];

            /* check whether the acceptor can take one of the donor's protons
             */
            pbc_dx(pbc, x[acceptor.atom_id], x[donor.atom_id], dx);
            if (norm(dx) <= DOO)
            {
                /* find the proton that can go. I assume that there can only be one!
                 */
                found_proton = FALSE;
                r_close = 100000;
                closest = -1;
                for (k=0; k < donor.nr_protons; k++)
                {
                    if(get_proton_presence(Hmap, donor.protons[k]) /* md->chargeA[donor.protons[k]] > 0.000001 */)
                    {
                        pbc_dx(pbc, x[donor.protons[k]], x[acceptor.atom_id], dx);
                        /* if(norm(dx)< HBOND) */
                        {
                            found_proton = TRUE;
                            /* norm2 is faster */
                            if (norm2(dx) < r_close)
                            {
                                r_close =  norm2(dx);
                                closest = k;
                            }
                            /* break; */
                        }
                    }
                }

                k = closest;

                if(found_proton)
                {
                    /* compute DHA angle */
                    /* use gmx_angle() instead? */
                    pbc_dx(pbc, x[donor.atom_id],    x[donor.protons[k]], veca);
                    pbc_dx(pbc, x[donor.protons[k]], x[acceptor.atom_id], vecb);
                    pbc_dx(pbc, x[acceptor.atom_id], x[donor.atom_id],    vecc);
		
                    /*ang2 = RAD2DEG*gmx_angle(vecc, veca);*/
                    ang = RAD2DEG*acos( (norm2(vecb) + norm2(veca) - norm2(vecc) )/
                                        (2 * norm(veca) * norm(vecb)));
                    /*fprintf(fplog,"ang = %g, ang2 = %g\n",ang,ang2);*/
                    /* ((M_PI-ang)*180.0/(M_PI) >= (HOPANG)) */
                    /*		  if ((ang)*180.0/(M_PI) >= (HOPANG))*/
                    if (ang > 180-MAXHOPANG)
                    {
                        /* the geometric conditions are right for a hop to take place
                         */
                        if (nr_hop >= T->max_nr_hop)
                        {
                            T->max_nr_hop++;
                            srenew(T->hop,T->max_nr_hop);
                            memset(&(T->hop[nr_hop]),0,sizeof(T->hop[0]));
                            snew(T->hop[nr_hop].f_before,md->nr);
                            snew(T->hop[nr_hop].f_after,md->nr);
                        }
                        T->hop[nr_hop].donor_id    = don;
                        T->hop[nr_hop].acceptor_id = acc;
                        T->hop[nr_hop].proton_id   = k;
                        T->hop[nr_hop].regime      = NOTSET;
                        T->hop[nr_hop].primary_d   = NOTSET;
                        T->hop[nr_hop].primary_a   = NOTSET;
                        T->hop[nr_hop].rda         = norm(vecc);
                        T->hop[nr_hop].ang         = ang;
                        T->hop[nr_hop].prob        = 0;
                        T->hop[nr_hop].bDonated    = FALSE;
                        nr_hop++;
                    }
                }
            }
        }
    }

    T->nr_hop = nr_hop;

} /* find_acceptors */


static void build_rotation_matrix(const real angle, const rvec axis, matrix *r)
{
    real theta;

    theta = angle*M_PI/180.0;

    (*r)[0][0] = axis[0]*axis[0]+(1.0-axis[0]*axis[0])*cos(theta);
    (*r)[0][1] = axis[0]*axis[1]*(1.0-cos(theta))-axis[2]*sin(theta);
    (*r)[0][2] = axis[0]*axis[2]*(1.0-cos(theta))+axis[1]*sin(theta);
  
    (*r)[1][0] = axis[0]*axis[1]*(1.0-cos(theta))+axis[2]*sin(theta);
    (*r)[1][1] = axis[1]*axis[1]+(1.0-axis[1]*axis[1])*cos(theta);
    (*r)[1][2] = axis[1]*axis[2]*(1.0-cos(theta))-axis[0]*sin(theta);
  
    (*r)[2][0] = axis[0]*axis[2]*(1.0-cos(theta))-axis[1]*sin(theta);
    (*r)[2][1] = axis[1]*axis[2]*(1.0-cos(theta))+axis[0]*sin(theta);
    (*r)[2][2] = axis[2]*axis[2]+(1.0-axis[2]*axis[2])*cos(theta);  
}

/* If bFlip==TRUE, then the rotation axis is O-H1 + O-H2 */
/* Rotate the group around the donor/acceptor. */
static void _rotate_group(t_pbc *pbc, rvec *x, t_qhop_atom *res,
                          real angle, rvec axis)
{
    int
        i,j;
    matrix
        r;
    rvec
        xnew ,temp;
    real
        n;

    clear_rvec(xnew);
    clear_rvec(temp);

    build_rotation_matrix(angle, axis, &r);

    for(j=0; j < res->nr_protons; j++)
    {
        rvec_sub(x[res->protons[j]], x[res->atom_id], temp); /* shift origin */
        mvmul(r, temp, xnew);                                /* rotate */
        rvec_add(xnew, x[res->atom_id], x[res->protons[j]]); /* re-shift origin */
    }
} /* _rotate_group */

static void rotate_water(t_pbc *pbc, rvec *x, t_qhop_atom *res, real angle)
{
    rvec axis, temp;

    /* The rotation axis is defined by the O-<H> vector. */

    pbc_dx(pbc, x[res->protons[2]], x[res->atom_id], axis);
    pbc_dx(pbc, x[res->protons[3]], x[res->atom_id], temp);  rvec_inc(axis, temp);
    pbc_dx(pbc, x[res->protons[4]], x[res->atom_id], temp);  rvec_inc(axis, temp);

    unitv(axis, axis);

    _rotate_group(pbc, x, res, angle, axis);
}

static void flip_water(t_pbc *pbc, rvec *x, t_qhop_atom *res)
{
    rvec axis, temp;
  
    /* The rotation axis is defined by the bisecor of the water */

    pbc_dx(pbc, x[res->protons[0]], x[res->atom_id], axis);
    pbc_dx(pbc, x[res->protons[1]], x[res->atom_id], temp);
    rvec_inc(axis, temp);

    unitv(axis, axis);

    _rotate_group(pbc, x, res, 180, axis);
}

static gmx_bool is_DA(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod, int DA)
{
    int i, j, rt;
    qhop_reactant *qres;

    if (DA != eQDON && DA != eQACC)
    {
        gmx_fatal(FARGS, "Called is_DA() with DA != eQDON || eQACC");
    }

    rt = qr[qa->qres_id].rtype;

    if ((DA == eQACC ?
         db->rb.qrt[rt].subres[prod].na :
         db->rb.qrt[rt].subres[prod].nd)
        > 1)
    {
        /* Loop over acceptors on this residue */
        for (i=0;
             i < (DA == eQACC ?
                  db->rb.qrt[rt].subres[prod].na :
                  db->rb.qrt[rt].subres[prod].nd);
             i++)
        {
            /* Loop over tautomers */
            qres =  (DA == eQACC) ?
                &(db->rb.qrt[rt].subres[prod].acc[i]) :
                &(db->rb.qrt[rt].subres[prod].don[i]);

            /* See if an atom matches */
            for (j=0; j < qres->nname; j++)
            {
                if (strcmp(qres->name[j], qa->atomname) == 0)
                {
                    return TRUE;
                }
            }
        }
    }

    return FALSE;
}

/* is_donor() and is_acceptor() are wrappers for is_DA() */
static gmx_bool is_donor(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod)
{
    return is_DA(db, qr, qa, prod, eQDON);
}

static gmx_bool is_acceptor(qhop_db *db, t_qhop_residue *qr, t_qhop_atom *qa, int prod)
{
    return is_DA(db, qr, qa, prod, eQACC);
}

/* Don't use qhop_titrate directly. Use qhop_protonate()
 * or qhop_deprotonate().
 *
 * Arguments are like qhop_(de)protonate(),
 * with the addition of DA which,
 * if DA==eQDON, makes a deprotonation and,
 * if DA==eQACC, a protonation.
 * Returns the subres being the product. */
static int qhop_titrate(qhop_db *db, titration_t T,
                        const t_inputrec *ir,
                        const t_commrec *cr, gmx_localtop_t *top,
                        gmx_constr_t constr,
                        t_qhop_atom *qatom, t_mdatoms *md,
                        const gmx_bool bWater, const gmx_bool bSwapBondeds, const int DA)
{
    int
        i, rt, r, reactant, t;
    qhop_subres
        *prod_res;
    t_qhop_residue
        *qres;
    gmx_bool bDonor;

    if (DA != eQDON && DA != eQACC)
    {
        gmx_fatal(FARGS, "Called qhop_titrate() with DA != eQDON || eQACC");
    }
  
    bDonor = (DA & eQDON);

    qres = &(T->qhop_residues[qatom->qres_id]);
    rt = qres->rtype;
    r = qres->subres;

    prod_res = NULL;


    /* Find out what the product is. */

    /* Loop over acceptors or donors */
    for (reactant = 0;
         (reactant < (bDonor ? (db->rb.qrt[rt].subres[r].nd) : (db->rb.qrt[rt].subres[r].na)))
             && prod_res==NULL;
         reactant++)
    {
        /* Loop over tautomers */
        for (t = 0;
             t < (bDonor ? (db->rb.qrt[rt].subres[r].don[reactant].nname) : (db->rb.qrt[rt].subres[r].acc[reactant].nname));
             t++)
        {
            /* Is this the reactant? */
            if (!strcasecmp(bDonor ?
                            db->rb.qrt[rt].subres[r].don[reactant].name[t] :
                            db->rb.qrt[rt].subres[r].acc[reactant].name[t],
                            qatom->atomname))
            {
                /* This is the reactant! */
                prod_res = bDonor ?
                    db->rb.qrt[rt].subres[r].don[reactant].productdata :
                    db->rb.qrt[rt].subres[r].acc[reactant].productdata;
                break;
            }
        }
    }

    if (prod_res == NULL)
    {
        gmx_fatal(FARGS, "Could not find the product for this reaction.");
    }

    /* Find the index in db->rb.res[rt][] */
    for (i=0; i<db->rb.qrt[rt].nsubres; i++)
    {
        if (strcmp(prod_res->name,
                   db->rb.qrt[rt].subres[i].name) == 0)
        {
            break;
        }
    }

    if (i == db->rb.qrt[rt].nsubres)
    {
        gmx_fatal(FARGS, "Could not find the product for this reaction.");
    }

    /* Set the residue subtype */
    qres->subres = i;

    /* Change interactions 
       For water we don't change any bonded or VdW. For now. 
       WARNING 
    */
    if (bWater)
    {
        /* Ugly hardcoding 
           WARNING 
        */
        set_proton_presence(&(db->H_map), qatom->protons[0], bDonor);
        set_proton_presence(&(db->H_map), qatom->protons[1], bDonor);
        set_proton_presence(&(db->H_map), qatom->protons[2], !bDonor);
        set_proton_presence(&(db->H_map), qatom->protons[3], !bDonor);
        set_proton_presence(&(db->H_map), qatom->protons[4], !bDonor);
    }
    else
    {
        if (bSwapBondeds)
        {
            if (db->bConstraints)
            {
                if (bDonor)
                {
                    qhop_deconstrain(qres, db, top, md, qatom->protons[qatom->nr_protons-1], constr, ir, cr);
                }	    
                else
                {
                    qhop_constrain(qres, T, db, top, md, qatom->protons[qatom->nr_protons-1], constr, ir, cr);
                }
            }
            unindex_ilists(qres); /* Why? For serial code this would be redundant */
            qhop_swap_bondeds(qres, prod_res, db, top, cr);
        }
        set_proton_presence(&(db->H_map),
                            qatom->protons[qatom->nr_protons-1],
                            !bDonor);
    }

    qhop_swap_vdws(qres, prod_res, md, db);

    qhop_swap_m_and_q(qres, prod_res, md, db, T);

    return i;
}

static void wipe_titrating_sites(titration_t T, t_qhop_residue *reactant, t_mdatoms *md)
{
    int i;
    t_qhop_atom *qa;

    qa = T->qhop_atoms;

    /* Wipe donor/acceptor info of qhop_atoms for this residue. */
    for (i=0; i<reactant->nr_titrating_sites; i++)
    {
        qa[i].state = eQNONE;
        md->bqhopdonor[qa[i].atom_id]=FALSE;
        md->bqhopacceptor[qa[i].atom_id]=FALSE;
    }
}


/* Wrappers for qhop_titrate() */
static void qhop_protonate(qhop_db *db, titration_t T,
                           const t_inputrec *ir,
                           const t_commrec *cr, gmx_localtop_t *top,
                           gmx_constr_t constr,
                           t_qhop_atom *qatom,
                           t_mdatoms *md, gmx_bool bWater,
                           gmx_bool bSwapBondeds, gmx_bool bRealHop)
{
    int prod, i, nt, *t;
    gmx_bool isacc;
    t_qhop_atom *qa;

    if (bRealHop)
    {
        wipe_titrating_sites(T, &(T->qhop_residues[qatom->qres_id]), md);
    }

    prod = qhop_titrate(db, T, ir, cr, top, constr, qatom, md, bWater, bSwapBondeds, eQACC);
  
    if (!bRealHop)
    {
        return;
    }

    /* Since it just got a proton it must become a donor. */

    /* Also make tautomeric sites donors */
    nt = get_tautomeric_sites(db, T, T->global_atom_to_qhop_atom[qatom->atom_id], &t, eQDON);

    for (i=0; i<nt; i++)
    {
        qa = &(T->qhop_atoms[t[i]]);

        md->bqhopdonor[qa->atom_id] = TRUE;
        qa->state |= eQDON;

        /* Is is still an acceptor? */
        isacc =  is_acceptor(db, T->qhop_residues, qa, prod);
        md->bqhopacceptor[qa->atom_id] = isacc;

        if (isacc)
        {
            qa->state |= eQACC;
        }
        else
        {
            qa->state &= ~eQACC;
        }
    }
}

static void qhop_deprotonate(qhop_db *db, titration_t T,
                      const t_inputrec *ir,
                      const t_commrec *cr, gmx_localtop_t *top,
                      gmx_constr_t constr,
                      t_qhop_atom *qatom,
                      t_mdatoms *md, gmx_bool bWater,
                      gmx_bool bSwapBondeds, gmx_bool bRealHop)
{
    int prod, i, nt, *t;
    gmx_bool isdon;
    t_qhop_atom *qa;

    if (bRealHop)
    {
        wipe_titrating_sites(T, &(T->qhop_residues[qatom->qres_id]), md);
    }

    prod = qhop_titrate(db, T, ir, cr, top, constr, qatom, md, bWater, bSwapBondeds, eQDON);

    if (!bRealHop)
    {
        return;
    }

    /* Since it just lost a proton it must become an acceptor. */

    /* Also make tautomeric sites acceptors */
    nt = get_tautomeric_sites(db, T, T->global_atom_to_qhop_atom[qatom->atom_id], &t, eQACC);

    for (i=0; i<nt; i++)
    {
        qa = &(T->qhop_atoms[t[i]]);

        md->bqhopacceptor[qa->atom_id] = TRUE;
        qa->state |= eQACC;

        /* Is it still a donor? */
        isdon = is_donor(db, T->qhop_residues, qa, prod);
        /* md->bqhopdonor[qatom->atom_id] = isdon; 
           Changed to below DvdS 2011-03-05 */
        md->bqhopdonor[qa->atom_id] = isdon;

        if (isdon)
        {
            qa->state |= eQDON;
        }
        else
        {
            qa->state &= ~eQDON;
        }
    }

}

static gmx_bool change_protonation(FILE *fplog,
                                   t_commrec *cr, const t_inputrec *ir, 
                                   titration_t T, 
                                   t_mdatoms *md, t_hop *hop, rvec *x,rvec *v, 
                                   rvec *f_before, rvec *f_after,
                                   eHOP_t ehop,gmx_mtop_t *mtop,
                                   gmx_localtop_t *top, gmx_constr_t constr,
                                   t_pbc *pbc,gmx_bool bRealHop)
{
    /* alters the topology in order to evaluate the new energy. In case
       of hydronium donating a proton, we might have to change the
       position of the resulting HW1, and HW2 in order to create a
       correct water configuration after the hop.  This is easiest done
       by a reflection in the plane of Hw1, OW, HW2 in the line
       OW-HW. We do not need to put this all back, if the proposal is
       not accepted. Water-Hydronium hops will be topology
       swaps. hydronium-residue hops will lead to annihilation of
       hydronium, and its index. residue-water hops lead to creation of
       a new hydronium.
    */

    /* Note that we still have special treatment of water and its symmetry. */

    int
        i, p, ft, ap;
    t_qhop_residue
        *donor,*acceptor;
    t_qhop_atom
        *donor_atom,*acceptor_atom;
    gmx_bool
        bAccWater, bDonWater;
    double
        ang;
    rvec
        x_tmp,
        w_normal,
        w_don,
        OH1,
        OH2,
        w_third,
        *vvv,*fbefore,*fafter;
    real
        planedist, d;
    qhop_db 
        *db;
    const int
        /* 0 and 1 are ordinary water protons, 2 is the "out of plane" one
         * in hydronium, while 3 and 4 are in the water plane.
         * See if we can avoid this hardcoding. */
        out_of_plane_proton = 2;

    if (bRealHop) 
    {
        vvv = v;
        fbefore = f_before;
        fafter = f_after;
    }
    else
    {
        vvv = NULL;
        fbefore = NULL;
        fafter = NULL;
    }
    db = T->db;
    
    if (NULL != debug)
    {
        fprintf(debug,"There are %d atoms. bRealHop = %d\n",md->nr,bRealHop);
        fprintf(debug,"Donors: %d - %d, Acceptors: %d -  %d\n",
                hop->primary_d, hop->donor_id,hop->primary_a, hop->acceptor_id);
        pr_rvecs(debug,0,"before change_protonation",x,md->nr);
    }   
    if (eHOP_FORWARD == ehop)
    {
        qhop_tautomer_swap(T, x, vvv, fbefore, fafter, hop->primary_d, hop->donor_id);
        qhop_tautomer_swap(T, x, vvv, fbefore, fafter, hop->primary_a, hop->acceptor_id);
    }

    donor_atom    = &T->qhop_atoms[hop->primary_d];
    acceptor_atom = &T->qhop_atoms[hop->primary_a];
    donor         = &T->qhop_residues[donor_atom->qres_id];
    acceptor      = &T->qhop_residues[acceptor_atom->qres_id];

    bAccWater = db->rb.qrt[acceptor->rtype].bWater;
    bDonWater = db->rb.qrt[donor->rtype].bWater;

    /* in case of our special hydronium/waters, where the protons are
       simply a dummy charge, we might have to alter the positions of
       the real atoms.
    */
  
    /* check if we need to do rotation */ 

    if (eHOP_BACKWARD == ehop)
        ang = -120;
    else
        ang = 120;

    if (bDonWater)
    {
        /* assumuming the user uses the correct water topology..... */
        /* WARNING */
        fprintf(fplog,"WARNING: rotating water, proton_id = %d\n",hop->proton_id);
        switch (hop->proton_id)
        {
        case 3:
            rotate_water(pbc, x,donor_atom,-ang);
            break;
	  
        case 4:
            rotate_water(pbc, x,donor_atom, ang);
            break;
	  
        default:
            break;
        }
    }
    else
    {
        /* We have similar problems when there are many protons on a titratable site.
           Flipping (or similar) may be needed. */

        /* Take away the LAST proton. That makes sense from an atomname point of view.
         * Ammonium would have H1 H2 H3 H4 while ammonia would have H1 H2 H3. */

        /* So, we substitute the proton at the hopping position for
         * the last proton on the site. */
        if (hop->proton_id != donor_atom->nr_protons-1)
        {
            swap_rvec(x[donor_atom->protons[hop->proton_id]], 
                      x[donor_atom->protons[donor_atom->nr_protons-1]]);
            /* WARNING velocity ? */
            if (bRealHop)
                swap_rvec(v[donor_atom->protons[hop->proton_id]], 
                          v[donor_atom->protons[donor_atom->nr_protons-1]]);
        }
    }

    /* Check if the proton ends up on the wrong side of the acceptor */
    if (bAccWater)
    {
        if (eHOP_FORWARD == ehop)
            /* Do we need to flip? */
        {
            pbc_dx(pbc,
                   x[donor_atom->atom_id],
                   x[acceptor_atom->atom_id],
                   w_don);
            pbc_dx(pbc,
                   x[acceptor_atom->protons[0]],
                   x[acceptor_atom->atom_id],
                   OH1);
            pbc_dx(pbc,
                   x[acceptor_atom->protons[1]],
                   x[acceptor_atom->atom_id],
                   OH2);
            pbc_dx(pbc,
                   x[acceptor_atom->protons[out_of_plane_proton]],
                   x[acceptor_atom->atom_id],
                   w_third);

            cprod(OH1, OH2, w_normal);
	  
            hop->bFlip = (iprod(w_don, w_normal) * iprod(w_third, w_normal) < 0);
	  
            /* So, if <AD, (DH1 x DH2)> is negative,
               then the vsite construction needs flipping. */

            /* ERRATA: above line should be correct, but the c-parameter
             * is *negative* in the water/H3O model I'm using.
             * Therefore I check if the projection of the out-of-plane-proton
             * on the water normal points in the same
             * direction as the projection of the AD-vector.
             * <AD, w_normal> * <w_third, w_normal> > 0  <==>  no flipping needed.*/
        }
      
        /* If (eHOP_BACKWARD == ehop), we already tested for flipping */
      
        if (hop->bFlip)
        {
            flip_water(pbc, x, acceptor_atom);
        }
    }
    if (bRealHop && bDonWater && bAccWater)
    {
        /* Just swap atom positions. Don't use qhop_titrate(). */
        for (i=0; i<donor->nr_atoms; i++)
        {
            /* coords */
            swap_rvec(x[acceptor->atoms[i]], x[donor->atoms[i]]);
            /* vel */
            swap_rvec(vvv[acceptor->atoms[i]], vvv[donor->atoms[i]]);
            /* force */
            swap_rvec(fbefore[acceptor->atoms[i]], fbefore[donor->atoms[i]]);
            swap_rvec(fafter[acceptor->atoms[i]], fafter[donor->atoms[i]]);
        }
    }
    else
    {
        /* (de)protonate */
        if( ((donor_atom->state & eQDON) && (acceptor_atom->state & eQACC))
            || ((acceptor_atom->state & eQDON) && (donor_atom->state & eQACC)) )
        {
            if (eHOP_BACKWARD == ehop)
            {
                qhop_protonate  (db, T, ir, cr, top, constr, donor_atom,    md, bDonWater, bRealHop, bRealHop);
                qhop_deprotonate(db, T, ir, cr, top, constr, acceptor_atom, md, bAccWater, bRealHop, bRealHop);
            }
            else
            {
                qhop_protonate  (db, T, ir, cr, top, constr, acceptor_atom, md, bAccWater, bRealHop, bRealHop);
                qhop_deprotonate(db, T, ir, cr, top, constr, donor_atom,    md, bDonWater, bRealHop, bRealHop);
            }
        }
        
        if (!bAccWater)
        {
            ap = acceptor_atom->nr_protons-1; /* This is the proton to add. */
            if (eHOP_BACKWARD == ehop)
            {
                /* Put proton back. It's the last proton on the acdeptor/donor that 
                   goes in/out of existence. Diprotic acids may be the exceptions, 
                   but let's not worry about that now. 
                   WARNING.
                */
                copy_rvec(hop->xold, x[acceptor_atom->protons[ap]]);
                if (bRealHop)
                    copy_rvec(hop->vold, v[acceptor_atom->protons[ap]]);
            }
            else
            {
                /* We need to reposition the hydrogen correctly.
                   Save old position and velocity. Needed for undo */
                copy_rvec(x[acceptor_atom->protons[ap]], hop->xold); 
                if (bRealHop)
                    copy_rvec(v[acceptor_atom->protons[ap]], hop->vold); 
                
                /* Reposition the proton on the DA-line */
                pbc_dx(pbc,x[donor_atom->atom_id],x[acceptor_atom->atom_id],x_tmp);
                unitv(x_tmp, x_tmp);
                
                /* What type of bond? */
                /* Constraints WARNING flag */
                ft = db->rb.btype[ebtsBONDS];
                
                /* We can safely assume that all protons on the
                   acceptor are equal with respect to bonded params,
                   so we just pick the first one. 
                   WARNING. 
                */
                p = (qhop_get_proton_bond_params(db, T, acceptor_atom, top,
                                                 acceptor_atom->protons[ap], cr)  
                     + top->idef.ntypes - db->rb.ni);
                /*   position in ilib + offset           */
                
                d = top->idef.iparams[p].harmonic.rA;
                
                svmul(d, x_tmp, x_tmp);
                rvec_add(x[acceptor_atom->atom_id], x_tmp, x[acceptor_atom->protons[ap]]);
                copy_rvec(x[acceptor_atom->protons[ap]], hop->xnew);
            }
        }
    }
    if (eHOP_BACKWARD == ehop)
    {
        qhop_tautomer_swap(T, x, vvv, fbefore, fafter, hop->primary_d, hop->donor_id );
        qhop_tautomer_swap(T, x, vvv, fbefore, fafter, hop->primary_a, hop->acceptor_id );
    }
    if (NULL != debug)
        pr_rvecs(debug,0,"after change_protonation",x,md->nr);

    return TRUE; /* Ok, this is not very good for error checking. */
} /* change_protonation */

static void evaluate_energy(FILE *fplog,
                            t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
                            gmx_wallcycle_t wcycle, 
                            gmx_localtop_t *top,gmx_mtop_t *mtop, 
                            gmx_groups_t *groups,t_state *state,
                            t_mdatoms *md,t_fcdata *fcd,
                            t_graph *graph, t_forcerec *fr, 
                            gmx_vsite_t *vsite,rvec mu_tot,
                            gmx_large_int_t step,
                            gmx_enerdata_t *enerd,
                            rvec *f)
{
    real t,etot,terminate;
    tensor force_vir;
    int forceflags;
    titration_t T;
    
    T = fr->titration;
    
    forceflags = GMX_FORCE_STATECHANGED;
    if (!T->bFreshNlists)
    {
        forceflags |= GMX_FORCE_NS;
    }

    terminate = 0;

    if (vsite)
    {
        construct_vsites(NULL,vsite,state->x,nrnb,1,NULL,
                         top->idef.iparams,top->idef.il,
                         fr->ePBC,fr->bMolPBC,graph,cr,state->box);
    }
  
    do_force(fplog,cr,ir,step,nrnb,wcycle,top,mtop,groups,
             state->box,state->x,&state->hist,
             f,force_vir,md,enerd,fcd,
             state->lambda,graph,
             fr,vsite,mu_tot,step*ir->delta_t,NULL,NULL,FALSE,
             GMX_FORCE_ALLFORCES|GMX_FORCE_VIRIAL|GMX_FORCE_STATECHANGED|GMX_FORCE_NS);

    T->bFreshNlists = TRUE;

    /* Communicate stuff when parallel */
#ifdef GMX_MPI

    if (PAR(cr))
    {
        //    wallcycle_start(wcycle,ewcMoveE);
        gmx_impl("qhop in parallel");
        /* global_stat(NULL,cr,enerd,NULL,NULL,mu_tot,
           ir,NULL,FALSE,NULL,NULL,NULL,NULL,&terminate); */

        //    wallcycle_stop(wcycle,ewcMoveE);
    }
#endif
} /* evaluate_energy */

static void print_qhop_info(titration_t T)
{
    int i,j,k;
    
    /* print qhop atoms */
    for(i=0; i < T->nr_qhop_atoms; i++)
    {
        fprintf(stderr,"qhop atom [%d], atom_id [%d], res [%d] bDonor = %d, %d protons\n",i,
                T->qhop_atoms[i].atom_id,
                T->qhop_atoms[i].res_id,
                /*T->qhop_atoms[i].bdonor*/
                0,T->qhop_atoms[i].nr_protons);
    }
    for(i=0; i < T->nr_qhop_residues; i++)
    {
        for(k=0; k < T->qhop_residues[i].nr_atoms; k++)
        {
            fprintf(stderr,"%3d ",T->qhop_residues[i].atoms[k]);
            fprintf(stderr,"\n");
        }
    }
} /* print_qhop_info */

static void print_hop(t_hop hop)
{
    fprintf(stderr,"donor atom %d, acceptor atom %d, proton %d\n",
            hop.donor_id,hop.acceptor_id,hop.proton_id);
} /* print_hop */

/* return the residue index (in tres) if i in an atom index in a qhop_residue of interest */
static int is_res_nr(int i, int *tres, int nres, titration *t_rec)
{
    int r, ri, a;

    for (r=0; r<nres; r++)
    {
        ri = tres[r];
        for (a=0; a < t_rec->qhop_residues[ri].nr_atoms; a++)
        {
            if (i == t_rec->qhop_residues[ri].atoms[a])
            {
                return r;
            }
        }
    }

    return -1;
}

/* Redundant */
static gmx_bool is_in_reactant(int a, t_qhop_residue *reac, t_atoms *atoms)
{
    /* Is this OK? It's fast :-) */
    return reac->res_nr == atoms->atom[a].resind;

    /* int i; */

    /* for (i=0; i < reac->nr_atoms; i++) */
    /* { */
    /*     if (a == reac->atoms[i]) */
    /*         return TRUE; */
    /* } */
    /* return FALSE; */
}

static t_nblist* reduce_nblist(titration *t_rec,
                               int nhops,
                               t_hop *hops,
                               t_nblist *nlist)
{
    /* The idea here is to make a reduced nblist from an full one.
       It will contain all neighbors that are part of
       {donors} U {acceptors}.
       Later, further reduction will produce nblists that connect
       individual pairs of reactants by calling this function again. */

    int i, j, h, min_id, max_id, npairs, nres, maxres, ii, jj, nj0, nj1, restemp, r, r1, r2, rr1, rr2, ri1, ri2, jtot, nri;
    int *pairs=NULL, *res=NULL, *tres=NULL;
    gmx_bool bNew, bSort;
    t_nblist *nlist_short = NULL;
    t_qhop_residue *d, *a;
    t_hop *hop;

    /* min_id = INT_MAX; */
    /* max_id = 0; */

    maxres = t_rec->nr_hop*2; /* decent estimate. Not too big, not too small. */

    snew(res,  maxres);
    snew(tres, maxres);

    snew(nlist_short,1);

    /* nlist_short->enlist = nlist->enlist; */
    /* nlist_short->il_code = nlist->il_code; */
    *nlist_short = *nlist;
    nlist_short->nri = 0;
    nlist_short->maxnri = 0;
    nlist_short->nrj = 0;
    nlist_short->maxnrj = 0;
    nlist_short->iinr = NULL;
    nlist_short->iinr_end = NULL;
    nlist_short->jindex = NULL;
    nlist_short->jjnr = NULL;
    nlist_short->jjnr_end = NULL;
    /* Just copy the mutex and counter?
       It will not be run at the same
       time as the ordinary nonbondeds anyway */
    nlist_short->mtx = nlist->mtx;
    nlist_short->count = nlist->count;
    nres = 0;
    /* Find what pairs of residues are connected by looking at the hops */
    for (h=0; h < nhops; h++)
    {
        hop = &(hops[h]);
        d = &(t_rec->qhop_residues[hop->donor_id]);
        a = &(t_rec->qhop_residues[hop->acceptor_id]);

        bNew = TRUE;

        for (i=0; i<nres; i++)
        {
            if (res[i] == d->res_nr)
            {
                bNew = FALSE;             
                break;
            }
        }
        if (bNew)
        {
            while (maxres <= nres)
            {
                maxres += 4; /* Other chunk size? */
                srenew(res,  maxres);
                srenew(tres, maxres);
            }

            res[nres]    = d->res_nr;
            tres[nres++] = hop->donor_id;
        }
        
        bNew = TRUE;
        for (i=0; i<nres; i++)
        {
            if (res[i] == a->res_nr)
            {
                bNew = FALSE;
                break;
            }
        }
        if (bNew)
        {
            while (maxres <= nres)
            {
                maxres += 4; /* Other chunk size? */
                srenew(res,  maxres);
                srenew(tres, maxres);
            }

            res[nres]    = a->res_nr;
            tres[nres++] = hop->acceptor_id;
        }
    }

    /* Sort the residue list */
    bSort = FALSE;
    while (!bSort)
    {
        bSort = TRUE;
        for (r=0; r<nres-1; r++)
        {
            if (res[r]>res[r+1])
            {
                restemp = res[r];
                res[r] = res[r+1];
                res[r+1] = restemp;

                restemp = tres[r];
                tres[r] = tres[r+1];
                tres[r+1] = restemp;

                bSort = FALSE;
            }
        }
    }

    snew(pairs, nres*nres); /* Make it square. Who cares? It'll be small anyway. */

    /* Make a bolean matrix showing what residues interact */
    for (h=0; h < nhops; h++)
    {
        hop = &(hops[h]);

        r1 = (t_rec->qhop_residues[hop->donor_id]).res_nr;
        r2 = (t_rec->qhop_residues[hop->acceptor_id]).res_nr;
        
        /* Get positions in res[] */
        for (rr1=0; rr1<nres; rr1++)
        {
            if (r1 == res[rr1])
            {
                break;
            }
        }

        for (rr2=0; rr2<nres; rr2++)
        {
            if (r2 == res[rr2])
            {
                break;
            }
        }

        if (rr1 >= nres || rr2 >= nres)
        {
            gmx_fatal(FARGS, "Couldn't find residue in H-transfer list");
        }


        if (r1<r2)
        {
            pairs[rr1*nres + rr2] = 1;
        }
        else
        {
            pairs[rr2*nres + rr1] = 1;
        }
    }    


    /* scan nblist */

    /* The nblist is not very intuitive at first glance.
     * ---- i-particles ----
     * loop n from 0 to nlist->nri
     *   ii = nlist->iinr[n] = atom id of the i-particle.
     *   x = shX+ x[ii*3 + 0]
     *   y = shY+ x[ii*3 + 1]
     *   z = shZ+ x[ii*3 + 2],
     *     where sh? are shift vectors in XYZ:
     *   shX = shiftvec[is3 + 0]
     *   shX = shiftvec[is3 + 1]
     *   shX = shiftvec[is3 + 2],
     *     where is3 = nlist->shift[n] * 3
     *
     * ---- j-particles ----
     * loop k from nlist->jindex[n] to nlist->jindex[n+1]
     *   jnr = nlist->jjnr[k] = atom id of the j-particle
     *   x = x[j3 + 0]
     *   y = x[j3 + 1]
     *   z = x[j3 + 2]
     *
     * Here, we're only interested in the indices and shift vectors.
     */

    jtot = 0; /* corroesponds to the last element of nlist_short->jjnr */

    snew(nlist_short->jjnr, 1);
    snew(nlist_short->jindex, 1);
    nlist_short->jindex[0] = 0;


    for (i=0; i < nlist->nri; i++)
    {        
        bNew = TRUE; /* This is a new i-particle */

        ii = nlist->iinr[i];

        ri1 = is_res_nr(ii, tres, nres, t_rec);

        if (ri1 < 0)
        {
            continue;
        }

        r1 = res[ri1];

        /* Find pairs that link reacting groups */
        nj0 = nlist->jindex[i];
        nj1 = nlist->jindex[i+1];

        for (j=nj0; j<nj1; j++)
        {
            jj = nlist->jjnr[j];
            
            ri2 = is_res_nr(jj, tres, nres, t_rec);

            if (ri2 < 0)
            {
                /* These are not the particles you're looking for */
                continue;
            }

            r2 = res[ri2];

            /* If atoms are within the same reacting residue, or if they are in different residues that react, go! */
            if (pairs[ri1*nres + ri2] || pairs[ri2*nres + ri1] || ri1==ri2)
            {
                /* Add to new nblist. */
                if (bNew)
                {
                    bNew = FALSE;

                    nri = nlist_short->nri;

                    srenew(nlist_short->iinr, nri+1);
                    nlist_short->iinr[nri] = ii;

                    srenew(nlist_short->shift, nri+1);
                    nlist_short->shift[nri] = nlist->shift[i];

                    srenew(nlist_short->jindex, nri+2);

                    srenew(nlist_short->gid, nri+1);
                    nlist_short->gid[nri] = nlist->gid[i];

                    nlist_short->nri++;
                }

                srenew(nlist_short->jjnr, jtot+1);
                nlist_short->jjnr[jtot] = jj;
                jtot++;
                nlist_short->jindex[nri+1] = jtot; /* This will keep on changing */
            }
        }
    }

    nlist_short->nrj = jtot;

    sfree(res);
    sfree(pairs);

    return nlist_short;
}

static void get_self_energy(FILE *fplog, t_inputrec *ir, t_commrec *cr, t_qhop_residue donor, 
                            t_qhop_residue acceptor, t_mdatoms *md, 
                            gmx_localtop_t *top,
                            gmx_mtop_t *mtop,
                            gmx_wallcycle_t wcycle,
                            gmx_groups_t *groups,
                            t_fcdata *fcd, t_graph *graph, gmx_vsite_t *vsite,
                            rvec mu_tot,
                            rvec *x, t_pbc *pbc, t_forcerec *fr, matrix box, history_t *hist,
                            real *Ecoul_self,real *Evdw_self, t_nblists *nlists_pair,
                            t_nrnb *nrnb, int step)
{
    /* the internal energy of the hoppers together Is difficult because
       the internal energy of the individal residues, eg. if part of
       protein, will also change upon hopping. For water does not
       matter.
    */
    int
        iat, ia, jat, ja, i, j, nti, nvdwparam,  vdwparam, ivdw, DA, e, efirst, elast, ni;
    gmx_bool
        bSkip;
    real
        qa, qq, Ecoul, Evdw, ek, ec, eps, drsq,
        rinv, rinvsq, r6, r12, c6, c12, ce1, ce2, Ec;
    qhop_db
        *db;
    rvec
        vec;
    t_qhop_residue res;
    titration_t T;
    t_nblists _nblists_tmp;
    t_ilist ilist[F_NRE];
    t_nblists *nblists_tmp;
    t_idef idef_tmp;
    tensor force_vir;
    gmx_enerdata_t *enerd;

    T = fr->titration;                        
    db = T->db;

    Ecoul = 0.0;
    Evdw  = 0.0;
    ek    = 0.0;
    ec    = 0.0;

    if (EEL_RF(fr->eeltype))
    {
        ek = fr->k_rf;
        ec = fr->c_rf;
    }

    eps = fr->epsfac;

    if (fr->bvdwtab)
    {
        ivdw = 3;
    }
    else if (fr->bBHAM)
    {
        ivdw = 2;
    }
    else 
    {
        ivdw = 1;
    }

    nvdwparam = (ivdw==2) ? 3 : 2;

    /* Intermolecular part */
    for(i=0; i < donor.nr_atoms; i++)
    {
        ia = donor.atoms[i];
        iat = md->typeA[ia];
        qa = md->chargeA[ia];
        nti = nvdwparam * fr->ntype * iat;

        for(j=0; j < acceptor.nr_atoms; j++)
        {
            ja = acceptor.atoms[j];
            jat = md->typeA[ja];
            vdwparam = nti + nvdwparam * jat;

            /* Coulomb */
            qq = qa * md->chargeA[ja];

            pbc_dx(pbc, x[ia], x[ja], vec);
            drsq = norm2(vec);
            rinv = gmx_invsqrt(drsq);
            rinvsq = rinv*rinv;

            Ec = qq*(rinv+ek*drsq-ec); /* = qq*rinv if not RF. */
            
            /* Hack to remove the erfc call. */
            if (0 && fr->bEwald)
            {
                Ecoul += Ec*erfc(fr->ewaldcoeff/rinv);
            }
            else
            {
                Ecoul +=  Ec;
            }

            /* Add vdw! */
            switch (ivdw)
            {
            case 1:
            case 3:
                /* LJ */
                r6  = rinvsq*rinvsq*rinvsq;
                r12 = r6*r6;
                c6  = fr->nbfp[vdwparam];
                c12 = fr->nbfp[vdwparam+1];
                Evdw += c12*r12 - c6*r6;

                break;

            case 2:
                /* Buckingham */
                r6  = rinvsq*rinvsq*rinvsq;
                c6  = fr->nbfp[vdwparam];
                ce1 = fr->nbfp[vdwparam+1];
                ce2 = fr->nbfp[vdwparam+2];
                Evdw += ce1*exp(-ce2*drsq*rinv) - c6*r6;
                break;

                /*case 3: */
                /* Tabulated */
                /* gmx_fatal(FARGS, "Tabulated vdw not yet supported with qhop. Stay tuned.");
                break;
                */
            default:
                gmx_fatal(FARGS, "Unknown vdw type!");
                break;
            }
        }
    }

#if 0 /* Removing 14 for now */
    /* Now remove intramolecular non-bonded stuff */
    for (DA=0; DA<2; DA++)
    {
        /*  DA==0 => DONOR
         *  DA==1 => ACCEPTOR  */
      
        res = DA==0 ? donor : acceptor;

        for(i=0; i < res.nr_atoms; i++)
        {
            ia = res.atoms[i];
            iat = md->typeA[ia];
            qa = md->chargeA[ia];
            nti = nvdwparam * fr->ntype * iat;
	  
            efirst = top->excls.index[ia];
            elast  = top->excls.index[ia+1];

            for(j=i; j < res.nr_atoms; j++)
            {
                ja = res.atoms[j];
                jat = md->typeA[ja];
                vdwparam = nti + nvdwparam * jat;
	      
                /* Is this an exclusion? if so, skip it! */
                bSkip = FALSE;
                for (e=efirst; e<elast; e++)
                {
                    if (ja == top->excls.a[e])
                    {
                        bSkip = TRUE;
                        break;
                    }
                }

                if (bSkip)
                {
                    continue;
                }

                /* Coulomb */
                qq = qa * md->chargeA[ja];

                pbc_dx(pbc, x[ia], x[ja], vec);
                drsq = norm2(vec);
                rinv = gmx_invsqrt(drsq);
                rinvsq = rinv*rinv;

                Ec = qq*(rinv+ek*drsq-ec);
                if (fr-bEwald)
                {
                    Ecoul += Ec*erfc(fr->ewaldcoeff/rinv);
                }
                else
                {
                    Ecoul +=  Ec;
                }

                /* Add vdw! */
                switch (ivdw)
                {
                case 1:
                    /* LJ */
                    r6  = rinvsq*rinvsq*rinvsq;
                    r12 = r6*r6;
                    c6  = fr->nbfp[vdwparam];
                    c12 = fr->nbfp[vdwparam+1];
                    Evdw += c12*r12 - c6*r6;

                    break;

                case 2:
                    /* Buckingham */
                    r6  = rinvsq*rinvsq*rinvsq;
                    c6  = fr->nbfp[vdwparam];
                    ce1 = fr->nbfp[vdwparam+1];
                    ce2 = fr->nbfp[vdwparam+2];
                    Evdw += ce1*exp(-ce2*drsq*rinv) - c6*r6;
                    break;

                case 3:
                    /* Tabulated */
                    gmx_fatal(FARGS, "Tabulated vdw not yet supported with qhop. Stay tuned.");
                    break;

                default:
                    gmx_fatal(FARGS, "Unknown vdw type!");
                    break;
                }
            }
        }
    }
#endif /* 1-4 interactions */

    Ecoul *= eps;

    *Ecoul_self = Ecoul;
    *Evdw_self  = Evdw;

#ifdef TITRATION_NB_KERNELS

    /* reaction field correction STILL needed!!!! nu even geen zin in*/
    /* return; */

    /* Instead: use the non-bonded kernels to get the numbers right.
     * Strip the nblist and feed it to the kernel. */
    int flags;
    /* int min_id, max_id; */
    /* gmx_bool bReac; */
    /* t_nblist *nlist_short; */
    
    /* /\* Can't do this with CG *\/ */
    /* if (nlist->enlist == enlistCG_CG) */
    /* { */
    /*     gmx_fatal(FARGS, "Can't do titration with CG."); */
    /* } */
    
    /* /\* === Strip nblist === *\/ */
    /* nlist_short = reduce_nblist(t_rec, nlist) */

    /* /\* Find highest and lowest atom id in the titrating residues to speed up search. *\/ */
    
    /* min_id = INT_MAX; */
    /* max_id = 0; */

    /* for (i=0; i<donor->nr_atoms; i++) */
    /* { */
    /*     min_id = min(min_id, donor->atoms[i]); */
    /*     max_id = max(max_id, donor->atoms[i]); */
    /* } */

    /* for (i=0; i<acceptor->nr_atoms; i++) */
    /* { */
    /*     min_id = min(min_id, acceptor->atoms[i]); */
    /*     max_id = max(max_id, acceptor->atoms[i]); */
    /* } */

    /* === Extract the 1-4 interactions that needs undoing === */
    
#define N_TERMS 3
    const int term[N_TERMS] = {F_LJ14, F_COUL14, F_LJC14_Q};

    memset(ilist, 0, sizeof(t_ilist)*F_NRE); /* Just make sure it's all NULL */
    for (i=0; i < N_TERMS; i++)
    {
        ni = 0;
        for (j=0; j < donor.bindex.nr[term[i]] ; j++)
        {
            {
                /* Always ft + 2 atoms = 3 elements for these interactions */
                ni = j*3;
                srenew(ilist[term[i]].iatoms, ni+3);
                memcpy(&(ilist[term[i]].iatoms[ni]), &(top->idef.il[term[i]].iatoms[donor.bindex.ilist_pos[i][ni]]), 3);
            }
        }
        if (ilist[term[i]].iatoms)
        {
            ilist[term[i]].nr     = ni+3;
            ilist[term[i]].nalloc = ni+3;
        }
    }

#undef N_TERMS

    /* Add the LJ/Coul */

    /* === Call non-boded kernels === */

    snew(enerd, 1);

    init_enerdata(groups->grps[egcENER].nr, ir->n_flambda, enerd);

    /* Temporarily switch the nblists and t_idef in ft, top and mtop . */

    nblists_tmp = fr->nblists;
    fr->nblists = nlists_pair;
    idef_tmp = top->idef;
    for (i=0; i<F_NRE; i++)
    {
        top->idef.il[i] = ilist[i];
    }

    do_force(fplog, cr, ir, step, nrnb, wcycle, top, mtop, groups, box, x, hist,
             T->f,
             force_vir,
             md,
             enerd,
             fcd,
             0 /* lambda */,
             graph,
             fr,
             vsite,
             mu_tot,
             step*ir->delta_t,
             NULL,
             NULL,
             FALSE,
             GMX_FORCE_BONDED | /* GMX_FORCE_NONBONDED | */ GMX_FORCE_SUBSYSTEM);

    /* Swap back. */
    fr->nblists = nblists_tmp;
    top->idef = idef_tmp;

    *Ecoul_self = enerd->term[F_COUL_SR]
        + enerd->term[F_COUL_LR]
        + enerd->term[F_RF_EXCL]
        + enerd->term[F_COUL_RECIP]
        /* + enerd->term[F_COUL14] */
        + enerd->term[F_RF_EXCL]
        + enerd->term[F_POLARIZATION];

    *Evdw_self = enerd->term[F_LJ]
        /* + enerd->term[F_LJ14] */
        + enerd->term[F_BHAM]
        + enerd->term[F_LJ_LR]
        + enerd->term[F_BHAM_LR]
        + enerd->term[F_DISPCORR];

        /* And these? */
  /* + enerd->term[F_LJC14_Q] */
  /* + enerd->term[F_LJC_PAIRS_NB] */

  /* + enerd->term[F_COUL_SR] */
  /* + enerd->term[F_COUL_LR] */
  /* + enerd->term[F_COUL_RECIP] */
  /* + enerd->term[F_DPD] */
  /* + enerd->term[F_POLARIZATION] */
  /* + enerd->term[F_WATER_POL] */
  /* + enerd->term[F_THOLE_POL] */
#endif /* TITRATION_NB_KERNELS */

    return;

} /* get_self_energy */

static real calc_S(const qhop_parameters *p,
                   const t_hop *hop)
{
    real d = hop->rda - (p->t_A);

    return (p->s_A) * (d * d) + (p->v_A);
}

static real calc_V(const qhop_parameters *p,
                   const t_hop *hop)
{
    return (p->s_C) * exp( -(p->t_C) * (hop->rda - 0.20)) + (p->v_C); 
}


/**
 * Return the transfer probability for an arbitrary time window
 *
 * Models the transfer as a poisson process,
 * p(t) =  1-exp(-k*t) = 1 - ( 1-p(0.01 ps) )^t
 * 
 */
static real poisson_prob(real p,            ///< Hopping probability over a 0.01 ps window
                         real qhop_interval ///< The interval between qhop attempts in picoseconds
    )
{
    return 1-pow((1-p), (qhop_interval/0.01));
}

static real calc_K(qhop_parameters *p, t_hop *hop)
{
    return p->k_1 * exp(-p->k_2 * (hop->rda-0.23)) + p->k_3;
}

static real calc_M(qhop_parameters *p, t_hop *hop)
{
    return p->m_1 * exp(-p->m_2 * (hop->rda-0.23)) + p->m_3;
}
/**
 * \brief Calculates the SE-regime validity limit and stores it in hop->El.
 *
 * \param p    Pointer to hopping parameters
 * \param hop  Pointer to a t_hop */
static void compute_E12_left(qhop_parameters *p,t_hop *hop)
{
  
    real
        K, M;

    /* compute the value of E12 that is required for 
     * p_SE(rda,10fs) = 0.1. 
     */
    K = calc_K(p, hop);
    M = calc_M(p, hop);

    hop->El = -(atanh( 2*(0.1-0.5) )-M) / K;
}


/**
 * Calculates the TST-regime validity limit and stores it in hop->Er.
 */
static void compute_E12_right(qhop_parameters *p, ///< Pointer to hopping parameters
                              t_hop *hop         ///< The hop
    )
{
  
    real
        S, T, V, Eb, f, g, A;

    /* compute the value of E12 at which Exp(Eb/kT)=100
     * Eb = S +T*E12 + V*E12*E12
     *
     *      -T + Sqrt(T*T - 4*V*(S-Eb))
     * E12 = --------------------------
     *                 2*V
     *
     * where Eb = 1/beta*ln(100)
     */
  
    /*  Eb = (BOLTZ*Temp)*log(100);
     *  No, we must use RGAS, since the units are in kJ/mol.
     */


    /* It gets worse. The quantity which is supposed to be
     * a factor log(100) larger than the thermal energy is
     * not Eb, but Eb-hbo. This leads to the following equation:
     *  Eb = f * exp(-g *Eb) + A  (1),
     * where A := h+log(100). Basically one needs to take the
     * ground state vibration into account. Eq. (1) can be
     * solved with the Labert W function (a.k.a. the product
     * logarithm). Using that Eb, we can solve the quadratic
     * equation described at the top of this function definition.
     */

    A = p->h + (log(100) * BOLTZ * hop->T);
    g = p->g;
    f = p->f;
    Eb = (A * g + LambertW(exp(-A * g) * f * g)) / g;

    /* Eb = RGAS * Temp * log(100); */
    /* d = hop->rda - (p->t_A); */
    S = calc_S(p, hop);
    T = p->s_B;
    V = calc_V(p, hop);

    hop->Er = (-T + sqrt(T*T - 4*V*(S-Eb))) / (2*V);
  
    /*   E12_right = ((-99.50-0.0760*Temp)*(rda*10.-(2.3450+0.000410*Temp))*(rda*10.-(2.3450+0.000410*Temp))+10.3)*CAL2JOULE; */
} /* compute_E12_right */


/**
 * Calculates the barrier height.
 */
static real compute_Eb(qhop_parameters *p, ///< Pointer to qhop_parameters			   
                       t_hop *hop
    )
{
    real
        Eb, S, T, V, temp;
  
    /* temp = rda - p->t_A; */
    S = calc_S(p, hop); /* p->s_A * (temp*temp) + p->v_A; */
    T = p->s_B;
    V = calc_V(p, hop); /* p->s_C*exp(-p->t_C*(rda-0.20))+p->v_C; */
    Eb = S + T*hop->E12 + V*hop->E12*hop->E12;
    return(Eb);
}

/**
 * \brief Calculates the TST transfer probability for a 10 fs time window and hbar omega/2.
 */
static real compute_rate_TST(qhop_parameters *p, ///< Pointer to qhop_parameters
                             t_hop *hop         ///< The hop
    )
{
    real
        Q, R, kappa, half_hbar_omega, ETST, pTST, E_M, rda, T;

    rda = hop->rda;
    T = hop->T;

    if (hop->E12 > 0.0)
    {
        E_M  = hop->Eb - hop->E12;
    }
    else
    {
        E_M = hop->Eb;
    }

    Q =
        p->q_1
        + p->q_2 * T
        + p->q_3 * T*T;

    R =
        p->r_1
        + p->r_2 * T
        + p->r_3 * T*T;

    kappa =
        exp(p->p_1 + Q*E_M + R*E_M*E_M);

    half_hbar_omega = p->f * exp(-p->g * hop->Eb) + p->h;
    ETST =-(hop->Eb-half_hbar_omega)/(BOLTZ * T);
    pTST = (kappa*BOLTZMANN*T/(PLANCK1 * 1e12)) * exp(ETST) * 0.01;

    hop->hbo = half_hbar_omega;

    return (pTST);
} /* compute_prob_TST */

/**
 * Calculates the barrierless transfer probability for a 10 fs time window.
 */
static real compute_rate_SE(qhop_parameters *p, ///< Pointer to qhop_parameters
                            t_hop *hop          ///< The hop
    )
{
    real
        K, M, pSE;
  
    K = calc_K(p, hop);
    M = calc_M(p, hop);
    pSE = (0.5 * tanh(-K * hop->E12 + M) + 0.5);
    return(pSE);
} /* compute_prob_SE */

static void dump_enerd(FILE *fplog,const char *s,gmx_enerdata_t *e,
                       const char *alg)
{
    int k;
    
    if (NULL != fplog)
    {
        for(k=0; (k<F_NRE); k++)
            if (e->term[k] != 0)
                fprintf(fplog,"%s: %s %s %g\n",alg,s,interaction_function[k].name,
                        e->term[k]);
    }
}

/**
 * Calculates the transfer probability for the intermediate regime over a 10 fs time window.
 */
static real compute_rate_log(qhop_parameters *p, ///< Pointer to qhop_parameters
                             t_hop *hop         ///< The hop 
    )
{
  
    /* we now compute the probability over a 10fs interval. We return the
     * probability per timestep.
     */
    real 
        rSE,rTST,rate;
  
    rSE  = compute_rate_SE(p, hop);
    rTST = compute_rate_TST(p, hop);
    rate = rSE * pow((rTST/rSE), (hop->E12 - hop->El)/(hop->Er - hop->El));
    /* See, the log-space interpolation translates to the expression in the first argument above.
     * No need to log/exp back and forth. */

    return(rate);
} /* compute_prob_log */

/* dE is the energy difference according to the force field, 'self'-interactions excluded. */
static void compute_E12(const qhop_parameters *p, t_hop *hop, real dE)
{
    hop->E12_0 =
        p->alpha
        + p->beta  * hop->rda
        + p->gamma * hop->rda * hop->rda;

    hop->E12 = hop->E12_0 + dE;
}

static void get_hop_prob(FILE *fplog,
                         t_commrec *cr, t_inputrec *ir, t_nrnb *nrnb,
                         gmx_wallcycle_t wcycle, 
                         gmx_localtop_t *top,gmx_mtop_t *mtop,
                         gmx_constr_t constr,
                         gmx_groups_t *groups,t_state *state,
                         t_mdatoms *md,t_fcdata *fcd,
                         t_graph *graph, t_forcerec *fr, 
                         gmx_vsite_t *vsite,rvec mu_tot,
                         t_hop *hop, t_pbc *pbc,gmx_large_int_t step, qhop_db *db,
                         gmx_bool *bHaveEbefore,
                         gmx_enerdata_t *Ebefore, gmx_enerdata_t *Eafter, t_nblists *nlists_reduced)
{
    /* compute the hopping probability based on the Q-hop criteria */
    int i,j,nnlists;
    real Ebefore_mm, Ebefore_self,
        Eafter_mm, Eafter_self,
        Ebefore_self_coul, Ebefore_self_vdw,
        Eafter_self_coul, Eafter_self_vdw,
        Edelta, Edelta_coul,
        r_TST, r_SE, r_log, Thop;
    qhop_parameters *p;
    titration_t T; 
    t_nblists *nlists_pair;
   
    T = fr->titration;
    
    /* Is this the primary titrating site, or do we need swapping? */
    hop->primary_d = qhop_get_primary(T,hop->donor_id,eQDON);
    hop->primary_a = qhop_get_primary(T,hop->acceptor_id,eQACC);

    p = get_qhop_params(T, hop, db);

    if (!(*bHaveEbefore))
    {
        evaluate_energy(fplog,cr, ir, nrnb, wcycle,top, mtop,
                        groups, state, md, fcd, graph,
                        fr,vsite,mu_tot,step,Ebefore,hop->f_before);
        dump_enerd(debug,"Before",Ebefore,eTitrationAlg_names[ir->titration_alg]);
        *bHaveEbefore = TRUE;
    }
    /* This needs to be evaluated for each possible hop */

#ifdef TITRATION_NB_KERNELS
    /* Make a nblist containing only the atoms form the residues involved in this hop.
       This is a slice of the full nblist(s), so it can be used for a good energy evaluation. */
    nnlists = fr->nnblists;
    nlists_pair = NULL;
    snew(nlists_pair, nnlists);
    
    for (i=0; i<nnlists; i++)
    {
        nlists_pair[i] = nlists_reduced[i];
        for (j=0; j<eNL_NR; j++)
        {
            nlists_pair[i].nlist_sr[j] = *(reduce_nblist(T, 1, hop, &(nlists_reduced[i].nlist_sr[j])));
            nlists_pair[i].nlist_lr[j] = *(reduce_nblist(T, 1, hop, &(nlists_reduced[i].nlist_lr[j])));
        }
    }
#endif

    get_self_energy(fplog, ir, cr,T->qhop_residues[T->qhop_atoms[hop->donor_id].qres_id],
                    T->qhop_residues[T->qhop_atoms[hop->acceptor_id].qres_id],
                    md, top, mtop, wcycle, groups, fcd, graph, vsite, mu_tot, state->x, pbc, fr,
                    state->box, &(state->hist),
                    &Ebefore_self_coul, &Ebefore_self_vdw, nlists_pair, nrnb, step);
    Ebefore_mm = (Ebefore->term[F_EPOT] - 
                  (Ebefore->term[F_COUL14] + Ebefore->term[F_LJ14]));
    Ebefore_self = (Ebefore_self_coul + Ebefore_self_vdw); 

    if (change_protonation(fplog, 
                           cr, ir, T, md, hop, state->x, state->v, hop->f_before, hop->f_after,
                           eHOP_FORWARD, mtop, top, constr, pbc, FALSE))
    {
        evaluate_energy(fplog,cr,ir,nrnb,wcycle,top,mtop,groups,state,md,
                        fcd,graph,fr,vsite,mu_tot,step, Eafter,hop->f_after);
        dump_enerd(debug,"After ",Eafter,eTitrationAlg_names[ir->titration_alg]);
        get_self_energy(fplog, ir, cr, T->qhop_residues[T->qhop_atoms[hop->donor_id].qres_id],
                        T->qhop_residues[T->qhop_atoms[hop->acceptor_id].qres_id],
                        md, top, mtop, wcycle, groups, fcd, graph, vsite, mu_tot, state->x, pbc, fr,
                        state->box, &(state->hist),
                        &Eafter_self_coul, &Eafter_self_vdw, nlists_pair, nrnb, step);
        Eafter_mm = (Eafter->term[F_EPOT] - 
                     (Eafter->term[F_COUL14] + Eafter->term[F_LJ14]));
        Eafter_self = (Eafter_self_coul + Eafter_self_vdw); 
        
        /* Correction for electronic polarization. If epsilon_r = 1, this term is 0 */
        Edelta_coul = (((Eafter->term[F_COUL_SR] + Eafter->term[F_COUL_LR] + 
                         Eafter->term[F_RF_EXCL] + Eafter->term[F_COUL_RECIP] -
                         Eafter_self_coul) -
                        (Ebefore->term[F_COUL_SR] + Ebefore->term[F_COUL_LR] + 
                         Ebefore->term[F_RF_EXCL] + Ebefore->term[F_COUL_RECIP] - 
                         Ebefore_self_coul)) * 
                       (1.0/ir->titration_epsilon_r - 1.0));
                       
        Edelta = ((Eafter_mm - Eafter_self) - (Ebefore_mm - Ebefore_self) +
                  Edelta_coul);
        hop->DE_Environment = Eafter_mm - Ebefore_mm;
        
        if (NULL != fplog)
            fprintf(fplog,"%s: Edelta %g, Eafter_mm %g, Eafter_self %g, Ebefore_mm %g, Ebefore_self %g, Edelta_coul %g, DE_Environment %g, DE_Self %g don: %d (%s) acc: %d (%s)\n%s: Eafter_self_coul %g, Ebefore_self_coul %g, Eafter_self_vdw %g, Ebefore_self_vdw %g\n",
                    eTitrationAlg_names[ir->titration_alg],
                    Edelta,Eafter_mm,Eafter_self,Ebefore_mm,Ebefore_self,Edelta_coul,
                    hop->DE_Environment,
                    Eafter_self-Ebefore_self,
                    T->qhop_atoms[hop->donor_id].res_id,
                    T->qhop_atoms[hop->donor_id].resname,
                    T->qhop_atoms[hop->acceptor_id].res_id,
                    T->qhop_atoms[hop->acceptor_id].resname,
                    eTitrationAlg_names[ir->titration_alg],
                    Eafter_self_coul,Ebefore_self_coul,Eafter_self_vdw, Ebefore_self_vdw);

        compute_E12(p, hop, Edelta);

        if (eTitrationAlgMC != ir->titration_alg)
        {
            hop->Eb        = compute_Eb(p,hop);
        }

        Thop = ir->delta_t * ir->titration_freq;

        if (eTitrationAlgQhop == ir->titration_alg)
        {
            compute_E12_right(p, hop);
            compute_E12_left (p, hop);
            
            if(hop->E12 > hop->Er)
            {
                /* Classical TST regime */
                r_TST = compute_rate_TST(p, hop);
                hop->prob = poisson_prob(r_TST, Thop);
                hop->regime = etQhopTST;
            }
            else if (hop->E12 > hop->El)
            {
                /* Intermediate regime */
                r_log = compute_rate_log(p, hop);
                hop->prob = poisson_prob(r_log, Thop);
                hop->regime = etQhopI;
            }
        }
        if ((eTitrationAlgICE == ir->titration_alg) || (hop->E12 <= hop->El))
        {
            /* Schroedinger regime */
            r_SE = compute_rate_SE (p, hop);
            hop->prob = poisson_prob(r_SE, (ir->delta_t * ir->titration_freq));
            hop->regime = etQhopSE;
        }
        if (eTitrationAlgMC == ir->titration_alg)
        {
            if (hop->E12 <= 0)
            {
                hop->prob = 1;
            }
            else
            {
                hop->prob = exp(-hop->E12 / (BOLTZ * hop->T));
            }
        }
        /* now undo the move */
        if (!change_protonation(fplog,
                                cr, ir, T, md, hop, state->x, state->v, hop->f_before,
                                hop->f_after,
                                eHOP_BACKWARD, mtop, top, constr, pbc, FALSE))
        {
            gmx_fatal(FARGS,"Oops, cannot undo the change in protonation.");
        }
    }
    else
    {
        fprintf(stderr, "change_protonation() returned FALSE!");
    }
  
    sfree(p);
#ifdef TITRATION_NB_KERNELS
    sfree(nlists_pair)
#endif
}
	    
real distance_dependence(rvec x, rvec c, real rc, t_pbc *pbc)
{
    real r;
    rvec dx;
  
    pbc_dx(pbc, x, c, dx);
    r = norm(dx);

    if (r < rc)
        return exp(-r/rc);
        /* return 1.0 - exp(-r/rc);
           This has gotta be wrong. The function should have
           a maximum at r=0 and decline with increasing r.
           Replace by exp(-r/rc) */
    else
        return 0.0;
}  /* distance_dependence */

static real check_ekin(rvec *v, t_mdatoms *md,char *title)
{
    int    i,j,k;
    real   m_2,ekin,mtot,mi;
    rvec   vcm,dv;
    tensor ekt;
    
    clear_rvec(vcm);
    mtot = 0;
    for(i=0;i<md->nr;i++)
    {
        if (md->ptype[i] == eptAtom) {
            mi    = md->massT[i];
            mtot += mi;
            for(j=0;j<DIM;j++)
            {
                vcm[j]+=v[i][j]*mi;
            }
        }
    }
    for(j=0;j<DIM;j++)
    {
        vcm[j] /= mtot;
    }
    
    clear_mat(ekt);
    /* WARNING this routine is probably not needed */
    for(i=0;i<md->nr;i++)
    {
        if (md->ptype[i] == eptAtom) {
            m_2 = 0.5*md->massT[i];
            rvec_sub(v[i],vcm,dv);
            for(j=0;j<DIM;j++)
            {
                for(k=0; (k<DIM); k++)
                {
                    ekt[j][k] += m_2*dv[j]*dv[k];
                }
            }   
        }
    }
    ekin = trace(ekt);
    if (debug) {
        pr_rvec(debug,0,"vcm",vcm,DIM,FALSE);
        pr_rvecs(debug,0,title,ekt,DIM);
        fprintf(debug,"%s: scalar ekin: %g\n",title,ekin);
    }
    return ekin;
} /* check_ekin */

/* Scales the velocities within radius rc. */
static gmx_bool scale_v(FILE *fplog, t_inputrec *ir,
                        rvec *x, rvec *v, t_mdatoms *md, 
                        int donor_atom, int acceptor_atom,
                        real DE,t_pbc *pbc)
{
    /* computes new velocities, to get ridf of the DE. Returns the new
       kinetic energy;
    */
    real Cinv, C, a0, a1, a2, mfi2, v2, C2, p2, q2, pq, lambda, mass;
    rvec center, p, q, dv;
    real *fi;
    int  i, j;


    /* compute the center of the reactants. 
     */
    for(j=0; j<DIM ; j++)
    {
        center[j] = 0.5*(x[donor_atom][j] + x[acceptor_atom][j]);
    }

    /* first compute f, C, p and q */  
    snew(fi, md->nr);
    Cinv = 0;
    clear_rvec(p);
    clear_rvec(q);
    for(i=0; i < md->nr; i++)
    {
        fi[i] = distance_dependence(x[i], center, ir->titration_vscale_radius, pbc);
        if (fi[i] > 0) 
        {
            mass = md->massT[i];
            Cinv += mass*fi[i];
            for(j=0; j < DIM; j++)
            {
                p[j] += mass * v[i][j];
                q[j] += mass * v[i][j] * fi[i];
            }	
        }
    }

    C  = 1/Cinv;
    C2 = C*C;

    /* with C, p and q, we can compute a, b and c that we need to find
       the solution for the quadratic equation of David and Erik M.
    */

    a0 = 2*DE;
    a1 = 0;
    a2 = 0;
    p2 = iprod(p,p);
    q2 = iprod(q,q);
    pq = iprod(p,q);
    for(i=0; i < md->nr; i++) 
    {
        mass = md->massT[i];
        mfi2 = mass * fi[i] * fi[i];
        if ((fi[i] > 0) && (mass > 0))
        {
            v2 = iprod(v[i],v[i]);
            a0 += mass * ( fi[i]*fi[i]*C2 * p2 - v2);
            a1 += 2 * mfi2 * ( C * iprod(v[i], p) - C2 * pq);
            a2 += mfi2 * ( v2 - 2 * C * iprod(v[i], q) + C2 * q2);
        }
    }
  
    /* With a2,a1,a0 we can find lambda. We need only the positive solution : */
    if (a1*a1-4*a0*a2 > 0) 
    {
        lambda = (-a1 + sqrt(a1*a1 - 4*a2*a0))/(2*a2);
    }
    else 
    {
        lambda = -1;
        if (NULL != fplog)
            fprintf(fplog,"%s: Not enough energy to perform a proton transfer. Better luck next time.\n",eTitrationAlg_names[ir->titration_alg]);
    }
    if (NULL != fplog)
        fprintf(fplog,"%s: lambda = %g. a2 = %g a1 = %g, a0 = %g\n",
                eTitrationAlg_names[ir->titration_alg],lambda,a2,a1,a0);
    if (lambda > 0.0 )
    {
        /* with lambda, we compute dv: */
        for(j=0; j<DIM; j++)
        {
            dv[j] = C * (p[j] - lambda*q[j]) / lambda;
        }
        
        /* with lambda and dv we scale the velocities of all particles */
        for(i=0; i < md->nr; i++)
        {
            if (fi[i] > 0) 
            {
                for(j=0; j < DIM; j++)
                {
                    v[i][j] = lambda * fi[i] * (v[i][j] + dv[j]);
                }
            }
        }
    }
    sfree(fi);
    
    return (lambda > 0);
} /* scale_v */

static gmx_bool scale_velocities(FILE *fplog,
                                 t_commrec *cr,t_inputrec *ir, t_nrnb *nrnb,
                                 gmx_wallcycle_t wcycle,  gmx_localtop_t *top,
                                 gmx_mtop_t *mtop, gmx_groups_t *groups,
                                 t_state *state,  t_mdatoms *md, 
                                 titration_t T,gmx_large_int_t step,
                                 gmx_constr_t constr,t_pbc *pbc,  
                                 t_hop *hop,gmx_ekindata_t *ekindata,
                                 real veta,real vetanew,rvec *f_hop,
                                 rvec *f_old)
{
    /* takes as input the total MM potential energy change for the hop
       and alters the kinetic energy so that the total energy remains
       constant.
    */
    gmx_bool bPscal, bConverged, bSufficientEkin;
    int      donor_atom, acceptor_atom, iter, maxiter, i, m;
    real     ekin_old, ekin_new, dvdl, ekin_before, DE, mI, dt_2m;
    char     stepstr[STEPSTRSIZE];
    rvec     *v;
    gmx_bool TEST = TRUE;
    
    if (ir->titration_vscale_radius <= 0)
        return TRUE;
    
    /* Compute new velocities */
    if ( TEST ) {
        snew(v,md->nr);
        for(i=0; (i<md->nr); i++) {
            mI = md->massT[i];
            if (mI > 0) {
                dt_2m = 0.5*ir->delta_t/mI;
                for(m=0; (m<DIM); m++) {
                    v[i][m] = state->v[i][m] + (f_old[i][m]+f_hop[i][m])*dt_2m;
                }
            }
        }
    }
    else
        v = state->v;
    bPscal = (ir->epc != epcNO) && ((step % ir->nstpcouple) == 0);
    
    /* iterate until the velocities satisfy both constraints and
     * energy/momentum conservation
     */

    donor_atom    = T->qhop_atoms[hop->donor_id].atom_id;
    acceptor_atom = T->qhop_atoms[hop->acceptor_id].atom_id;

    /* Compute ekin based on old velocities */
    ekin_before = check_ekin(v,md,"before");
    ekin_old = ekin_before;
    ekin_new = 0;
    /*  DE = hop->E12;*/
    DE = hop->DE_Environment;
   
    if (T->db->bConstraints)
        maxiter = 10;
    else
        maxiter = 10;
        
    bSufficientEkin = TRUE;
    bConverged = FALSE;
    for(iter=0; (iter<maxiter) && !bConverged && bSufficientEkin; iter++)
    {
        bSufficientEkin = scale_v(fplog,ir, state->x, v,
                                  md, donor_atom, acceptor_atom,
                                  DE, pbc);
    
        if (bSufficientEkin)
        {
            /* constrain the velocities */ 
            if (T->db->bConstraints)
            {
                if (NULL != fplog)
                    fprintf(fplog,"%s: Before constr. ekin_old = %f\n",
                            eTitrationAlg_names[ir->titration_alg],ekin_old);

                constrain(NULL,FALSE,FALSE,constr,&top->idef,ir,ekindata,cr,
                          step,1,md,state->x,v,v,state->box,
                          state->lambda,&dvdl,NULL,NULL,nrnb,
                          econqVeloc,bPscal,veta,vetanew);
                /* Correct the velocities for constraints */
                
                
                /* And send the velocities around again, bit expensive, 
                   so there must be an easier way
                   of doing this...
                */
            }
            
            ekin_new = check_ekin(v ,md, "new");   
            if (NULL != fplog)
            {
                fprintf(fplog,"%s: iteration %d, ekin_new = %f, ekin_old = %f. DE = %f\n",
                        eTitrationAlg_names[ir->titration_alg],iter,ekin_new,ekin_old,DE);
                pr_rvecs(fplog,0,"Original Ekin Tensor",
                         ekindata->ekin,DIM);
                fprintf(fplog,"Scalar ekin: %g\n",trace(ekindata->ekin));
            }
            DE = DE + (ekin_new - ekin_old);
            ekin_old = ekin_new;
            bConverged = (DE*DE < 1); /* hard coded treshold.... */
        }
    }
    if (NULL != fplog)
    {
        gmx_step_str(step, stepstr);
        fprintf(fplog,"%s: Energy correction at step %s: %f, DE_Environment: %f\n",
                eTitrationAlg_names[ir->titration_alg],stepstr,ekin_new-ekin_before,hop->DE_Environment);
    }
    if (bSufficientEkin && bConverged && TEST) {
        /* Copy new velocities */
        for(i=0; (i<md->nr); i++) {
            mI = md->massT[i];
            if (mI > 0) {
                dt_2m = 0.5*ir->delta_t/mI;
                for(m=0; (m<DIM); m++) {
                    state->v[i][m] = v[i][m] - (f_old[i][m]+f_hop[i][m])*dt_2m;
                }
            }
        }
    }
    sfree(v);
    
    return (bSufficientEkin && bConverged);
} /* scale_velocities */

static int hopcomp(const void *a,const void *b)
{
    t_hop *ha = (t_hop *)a;
    t_hop *hb = (t_hop *)b;
    real dd;
    if ((dd = (ha->prob-hb->prob)) < 0)
        return -1;
    else if (dd > 0)
        return 1;
    return 0;
}
  
/**
 * Changes the order of the hops
 */
static void scramble_hops(titration_t T, int mode)      ///< etQhopMode????
{
    int   i,r;
    t_hop tmp_hop;

    switch (mode)
    {
    case eTitrationModeOne:
    case eTitrationModeGillespie:
        break;

    case eTitrationModeRandlist:
        /* Shuffle the hops! */
        for (i=0; i<T->nr_hop-1; i++)
        {
            r = i+gmx_rng_uniform_uint32(T->rng) % (T->nr_hop - i);
            tmp_hop = T->hop[i];
            T->hop[i] = T->hop[r];
            T->hop[r] = tmp_hop;
        }
        break;

    case eTitrationModeList:
        /* Sort according to probability */
        qsort(T->hop,T->nr_hop,sizeof(T->hop[0]),hopcomp);
        break;

    default:
        gmx_fatal(FARGS, "That titration_mode is unsupported: %i", mode);
    }
}

int do_titration(FILE *fplog,
                 t_commrec *cr,
                 t_inputrec *ir,
                 t_nrnb *nrnb,
                 gmx_wallcycle_t wcycle, 
                 gmx_localtop_t *top,
                 gmx_mtop_t *mtop,
                 gmx_groups_t *groups,
                 t_state *state,
                 t_mdatoms *md,
                 t_fcdata *fcd,
                 t_graph *graph,
                 t_forcerec *fr,
                 gmx_constr_t constr,
                 gmx_vsite_t *vsite,
                 rvec mu_tot,
                 gmx_bool bBornRadii,
                 real Temperature,
                 gmx_large_int_t step,
                 gmx_ekindata_t *ekindata,
                 tensor force_vir,
                 rvec *ftitration,
                 real *DE_Titration)
{
    char
        stepstr[STEPSTRSIZE];
    int
        i,j,a,nTransfers, nnlists;
    gmx_bool
        bAgain,bHop;
    real 
        rnr;
    int 
        start_seed=0,iHop=0;
    t_pbc pbc;
    qhop_db
        *db;
    gmx_bool bHaveEbefore;
    gmx_enerdata_t *Ebefore,*Eafter;
    real
        veta, vetanew, DE_Env;
    titration_t T;
    
    t_nblists *nlists_reduced;
    
    T = fr->titration;
    db = T->db;

    set_pbc_dd(&pbc,fr->ePBC,DOMAINDECOMP(cr) ? cr->dd : NULL, FALSE, state->box);
    
    snew(Ebefore,1);
    snew(Eafter,1);
    init_enerdata(groups->grps[egcENER].nr, ir->n_flambda, Ebefore);
    init_enerdata(groups->grps[egcENER].nr, ir->n_flambda, Eafter);

    /* We need to regenerate the hop nlists only when needed. 
       Doesn't work yet though. */
    T->bFreshNlists = FALSE; 
    T->bHaveEbefore = FALSE;

    if(T->rng == NULL)
    {
        if(MASTER(cr))
        {
            start_seed = ir->titration_random_seed;
        }
#ifdef GMX_MPI
        if(PAR(cr))
        {
            MPI_Bcast(&start_seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
#endif
        T->rng     = gmx_rng_init(start_seed);
        T->rng_int = gmx_rng_init(start_seed);
    } 

    find_acceptors(fplog,cr, fr, T, state->x, &pbc, md, &(db->H_map));

    DE_Env = 0;
    
    if(T->nr_hop > 0)
    {
#ifdef TITRATION_NB_KERNELS
        /* make first reduction of neighborlists */
        nnlists = fr->nnblists;
        nlists_reduced = NULL;
        snew(nlists_reduced, nnlists);

        for (i=0; i<nnlists; i++)
        {
            /* if (!fr->nblists[i]) */
            /* { */
            /*     continue; */
            /* } */

            nlists_reduced[i] = fr->nblists[i];

            for (j=0; j<eNL_NR; j++)
            {
                nlists_reduced[i].nlist_sr[j] = *(reduce_nblist(T, T->nr_hop, T->hop, &(fr->nblists[i].nlist_sr[j])));
                nlists_reduced[i].nlist_lr[j] = *(reduce_nblist(T, T->nr_hop, T->hop, &(fr->nblists[i].nlist_lr[j])));
            }                
        }
#endif

        bHaveEbefore = FALSE;
        if (ir->titration_mode == eTitrationModeList) 
        {
            /* For this Q-hop style mode we have to compute all probabilities 
               first. That could mean unnecessary work! */
            for(i=0; i<T->nr_hop; i++)
            {
                T->hop[i].T = Temperature;
                get_hop_prob(fplog,cr,ir,nrnb,wcycle,top,mtop,constr,groups,state,md,
                             fcd,graph,fr,vsite,mu_tot,
                             &(T->hop[i]),&pbc,step, db,
                             &bHaveEbefore,Ebefore,Eafter, nlists_reduced);
            }
        }
        scramble_hops(T, ir->titration_mode);

        if (ir->titration_alg == eTitrationAlgMC)
        {
            /* Correct the probability for the different number of possible hops from different donors:
             *   p = p' * NTRANSFERS/MAXTRANSFERS,
             * where NTRANSFERS are the actual number of possible transfers from a site, and MAXTRANSFERS
             * is a constant>=NTRANSFERS. */

            for (i=0; i<T->nr_hop; i++)
            {
                nTransfers = 0;
                for (j=0; j<T->nr_hop; j++)
                {
                    if (T->hop[i].donor_id == T->hop[j].donor_id)
                    {
                        nTransfers++;
                    }
                }
                T->hop[i].prob *= nTransfers/MAX_MC_TRANSFERS;
            }
        }
        
        for(i=0; (i <T->nr_hop); i++)
        {
            if (!T->hop[i].bDonated)
            {
                if (ir->titration_mode != eTitrationModeList) 
                {
                    T->hop[i].T = Temperature;
                    get_hop_prob(fplog,cr,ir,nrnb,wcycle,top,mtop,constr,groups,state,md,
                                 fcd,graph,fr,vsite,mu_tot,
                                 &(T->hop[i]),&pbc,step, db,
                                 &bHaveEbefore,Ebefore,Eafter, nlists_reduced);
                }
            
                rnr = gmx_rng_uniform_real(T->rng); 
                if (NULL != debug)
                {
                    /* some printing to the log file */
                    fprintf(debug,
                            "\n%d. don %d acc %d. E12 = %4.4f, DE_Environment = %4.4f, Eb = %4.4f, hbo = %4.4f, rda = %2.4f, ang = %2.4f, El = %4.4f, Er = %4.4f prob. = %f, ran. = %f (%s)\n", i,
                            T->qhop_atoms[T->hop[i].donor_id].res_id,
                            T->qhop_atoms[T->hop[i].acceptor_id].res_id,
                            T->hop[i].E12, T->hop[i].DE_Environment, T->hop[i].Eb,
                            T->hop[i].hbo, T->hop[i].rda, T->hop[i].ang,
                            T->hop[i].El, T->hop[i].Er,
                            T->hop[i].prob, rnr, qhopregimes[T->hop[i].regime]);
                }
 
                /* Attempt hop */
                if (T->hop[i].prob > rnr)
                {
                    /* Hoppen maar! */
                    int k;
                    
                    if (0)
                        for(k=0; (k<md->nr); k++) 
                        {
                            fprintf(fplog,"B_ATOM %5d  M %10.5f V %10.5f  %10.5f  %10.5f\n",
                                    k,md->massT[k],state->v[k][XX],
                                    state->v[k][YY],state->v[k][ZZ]);
                        }
                    fprintf(fplog,"EKIN before CP %g\n",
                            check_ekin(state->v,md,"before"));
                    
                    /* Note different indices to T->hop */
                    bHop = change_protonation(fplog,
                                              cr, ir, T, md, &(T->hop[i]), 
                                              state->x, state->v, T->hop[0].f_before, 
                                              T->hop[i].f_after, 
                                              eHOP_FORWARD, mtop, top, constr, &pbc, TRUE);
                    if (0)
                        for(k=0; (k<md->nr); k++) 
                        {
                            fprintf(fplog,"A_ATOM %5d  M %10.5f V %10.5f  %10.5f  %10.5f\n",
                                    k,md->massT[k],state->v[k][XX],
                                    state->v[k][YY],state->v[k][ZZ]);
                        }
                    fprintf(fplog,"EKIN after CP %g\n",
                            check_ekin(state->v,md,"after"));
                    DE_Env += T->hop[i].DE_Environment;
                        
                    if (bHop)
                        iHop = i;
                        
                    if (bHop && 0)
                    {
                        veta = 0;
                        vetanew = 0;
                        bHop = scale_velocities(fplog,cr,ir,nrnb,wcycle,top,mtop,groups,
                                                state,md,T,step,constr,
                                                &pbc,&(T->hop[i]),ekindata,veta,vetanew,
                                                T->hop[i].f_after,ftitration);
                    }
                        
                    if (bHop && (ir->titration_mode != eTitrationModeOne))
                    {
                        /* Zap all hops whose reactants were just consumed. */
                        for (j = i+1; j < T->nr_hop; j++)
                        {
                            if (T->qhop_atoms[T->hop[j].donor_id].res_id    == T->qhop_atoms[T->hop[i].donor_id].res_id    ||
                                T->qhop_atoms[T->hop[j].acceptor_id].res_id == T->qhop_atoms[T->hop[i].acceptor_id].res_id ||
                                T->qhop_atoms[T->hop[j].donor_id].res_id    == T->qhop_atoms[T->hop[i].acceptor_id].res_id ||
                                T->qhop_atoms[T->hop[j].acceptor_id].res_id == T->qhop_atoms[T->hop[i].donor_id].res_id)
                            {
                                T->hop[j].regime = etQhopNONE;
                                T->hop[j].bDonated = TRUE;
                                if (NULL != debug)
                                    fprintf(debug," * Zapping hop number %i *\n", j);
                            }
                        }
                    }
                    if (!bHop)
                    {

                        /* Note different indices to T->hop */
                        if (!change_protonation(fplog,
                                                cr, ir, T, md, &(T->hop[i]), 
                                                state->x, state->v, T->hop[0].f_before, 
                                                T->hop[i].f_after, 
                                                eHOP_BACKWARD, mtop, top, 
                                                constr, &pbc, TRUE))
                            gmx_fatal(FARGS,"Oops could not hop back. WTH?");
                    }
                }
                else
                {
                    bHop = FALSE;
                }

                if ((NULL != fplog))
                {
                    gmx_step_str(step, stepstr);
                    fprintf(fplog,"%s: Ekin = %g\n",
                            eTitrationAlg_names[ir->titration_alg],
                            check_ekin(state->v,md,"Check"));
                    fprintf(fplog,"%s: %5s at step %s. E12 = %8.3f, hopper: %d don: %d (%s) acc: %d (%s)",
                            eTitrationAlg_names[ir->titration_alg],
                            (bHop ? "P-hop" : "no hop"),
                            stepstr,T->hop[i].E12, i,
                            T->qhop_atoms[T->hop[i].donor_id].res_id,
                            T->qhop_atoms[T->hop[i].donor_id].resname,
                            T->qhop_atoms[T->hop[i].acceptor_id].res_id,
                            T->qhop_atoms[T->hop[i].acceptor_id].resname);
                    fprintf(fplog,
                            " DE_Environment = %4.4f, Eb = %4.4f, rda = %2.4f, ang = %2.4f, El = %4.4f, Er = %4.4f prob. = %f, ran. = %f (%s)\n", 
                            T->hop[i].DE_Environment, T->hop[i].Eb,
                            T->hop[i].rda, T->hop[i].ang,
                            T->hop[i].El, T->hop[i].Er,
                            T->hop[i].prob, rnr, qhopregimes[T->hop[i].regime]);
                }

                if (ir->titration_mode == eTitrationModeOne)
                {
                    break;
                }
            }
        }
    }  
#ifdef TITRATION_NB_KERNELS
    sfree(nlist_reduced);
#endif
    if (!bHop)
        return eTitration_No;
    else {
        if (ir->titration_alg == eTitrationAlgQhop)
            return eTitration_NoEcorr;
        else
        {
            real hdt,ddek2,ekbefore,ekafter,ekafter2,dek2 = 0;
            rvec df,dv,vnew,vnew2,dv_after,dv_before,vbefore,vafter;
            
            /* How should we treat multiple hops? Check indices in T->hop
             */
            if (getenv("EK2") != NULL)
            {
                hdt = 0.5*ir->delta_t;
                ekbefore = ekafter = ekafter2 = 0;
                for(i=0; (i<md->nr); i++) 
                {
                    copy_rvec(T->hop[iHop].f_after[i],ftitration[i]);
                    /* Do we need to swap forces as well? */
                    rvec_sub(T->hop[iHop].f_after[i],T->hop[0].f_before[i],df);
                    svmul(md->invmass[i]*hdt,T->hop[0].f_before[i],dv_before);
                    svmul(md->invmass[i]*hdt,T->hop[iHop].f_after[i],dv_after);
                    
                    rvec_add(state->v[i],dv_before,vbefore);
                    ekbefore += 0.5*md->massT[i]*iprod(vbefore,vbefore);
                    
                    svmul(md->invmass[i]*hdt,df,dv);
                    rvec_add(vbefore,dv,vafter);
                    ekafter += 0.5*md->massT[i]*iprod(vafter,vafter);
                    
                    rvec_add(state->v[i],dv_after,vnew2);
                    ekafter2 += 0.5*md->massT[i]*iprod(vnew2,vnew2);
                    
                    /*rvec_add(state->v[i],dv,vnew);*/
                    ddek2 = md->massT[i]*(2*iprod(vbefore,dv)+iprod(dv,dv));
                    dek2 += ddek2;
                    fprintf(fplog,"df[%5d] = %8.3f  %8.3f  %8.3f  dek2 = %8.3f\n",
                            i,df[XX],df[YY],df[ZZ],ddek2);
                }
                dek2 /= 2;
                fprintf(fplog,"dek2 = %g ekbefore = %g ekafter = %g ekafter2 = %g\n",
                        dek2,ekbefore,ekafter,ekafter2);
            }
            *DE_Titration = DE_Env+dek2;
            return eTitration_Ecorr;
        }
    }
}

/* paramStr is the string from the rtp file.
 * bts = ebtsBONDS, ..., ebtsCMAP
 * bt is the type found in t_rbonded.type
 */

static void str2bonded(const char *paramStr,
                       const int bt, const int bts,
                       t_iparams *params, int *iatoms)
{

#ifndef _ALL_BPARAMS
#define _ALL_BPARAMS &p[0],&p[1],&p[2],&p[3],&p[4],&p[5],&p[6],&p[7],&p[8],&p[9],&p[10],&p[11] /* p[0:MAXFORCEPARAM] */
#endif

    int ft, i, niatoms=btsNiatoms[bts];
    real p[MAXFORCEPARAM];

    int np, strOffset;
    char format[1+2+MAXATOMLIST*2+MAXFORCEPARAM*2];
    /* 2 for the functype,
     * MAXATOMS*2 for the %i belonging to each atom,
     * MAXFORCEPARAMS*2 for the %f belonging to each forceparameter,
     * 1 for the '\0'
     */
    memset(p, 0, sizeof(real)*MAXFORCEPARAM);
    memset(params, 0, sizeof(params[0])); /* zap the params to zero */
    strOffset = btsNiatoms[bts]*2+1; /* +1 for the initial whitespace, *2 for the "%i" */
    sprintf(format," %s","%i%i%i%i%i%i%i"); /* the atoms */
    for (i=0; i<MAXFORCEPARAM; i++)
        sprintf(&(format[i*2+2+strOffset]),"%%f");
  
    np = sscanf(paramStr,"%d%d%d%d%d%d",&iatoms[0],&iatoms[1],&iatoms[2],
                &iatoms[3],&iatoms[4],&iatoms[5]);
    if (np < btsNiatoms[bts])
        gmx_fatal(FARGS,"Not enough arguments on line '%s'. Expected %d.",
                  paramStr,btsNiatoms[bts]);
    
    /*
    case 2:
        np = sscanf(paramStr, format,
                    iatoms[0],
                    iatoms[1],
                    _ALL_BPARAMS);
        break;
    case 3:
        np = sscanf(paramStr, format,
                    iatoms[0],
                    iatoms[1],
                    iatoms[2],
                    _ALL_BPARAMS);
        break;
    case 4:
        np = sscanf(paramStr, format,
                    iatoms[0],
                    iatoms[1],
                    iatoms[2],
                    iatoms[3],
                    _ALL_BPARAMS);
        break;
    case 5:
        np = sscanf(paramStr, format,
                    iatoms[0],
                    iatoms[1],
                    iatoms[2],
                    iatoms[3],
                    iatoms[4],
                    _ALL_BPARAMS);
        break;
    case 6:
        np = sscanf(paramStr, format,
                    iatoms[0],
                    iatoms[1],
                    iatoms[2],
                    iatoms[3],
                    iatoms[4],
                    iatoms[5],
                    _ALL_BPARAMS);
        break;
        default:*/
        /* This reeeeally shouldn't happen, unless a new bonded type is added with more than 5 atoms. */
/*        gmx_fatal(FARGS, "Wrong number of atoms requested: %d", btsNiatoms[bts]);
    }
*/
    /*   ====  Now set the params  ====
     * Note that tabulated potentials, vsiten, orires and cmap
     * are not supported at this stage, simply because such information
     * is not available in the rtp file. Other bonded types are striktly
     * speaking also not available, e.g. vsite, but since their parameters
     * all reals it's trivial to just pass the parameters to the t_iparams.
     */
    if (bts==ebtsPDIHS && bt==1) /* Ordinary proper dihedral. */
    {
        params->pdihs.phiA = p[0];
        params->pdihs.cpA  = p[1];
        params->pdihs.mult = (int)p[2]; /* <- This is why this one needs special treatment */
        params->pdihs.phiA = p[3];
        params->pdihs.cpA  = p[4];
    }
    else
    {
        for (i=0; i<MAXFORCEPARAM; i++)
            (*params).generic.buf[i] = p[i];
    }
}

static void qhop_stash_bonded(qhop_db_t db, gmx_mtop_t *mtop)
{
    int rt, r, nrt, nres, bt, bi,
        natoms, /* nqatoms,  */i, j, ft,
        iatoms[MAXATOMLIST],
        iatomsMapped[MAXATOMLIST],
        *atomMap=NULL,
        ilistEntry[MAXATOMLIST+1];

    t_restp *rtpr, *rtprt;

    /* WE NEED TO DETERMINE THE FUNCTYPE ft !!!!*/

    t_iparams params;
  
    ft = -1;
    nrt = db->rb.nrestypes;
    for (rt=0; rt < nrt; rt++)
    {
        nres = db->rb.qrt[rt].nsubres;
        rtprt = &(db->rtp[db->rb.qrt[rt].irtp]);
        /* nqatoms = rtp->natom; */

        for (r=0; r < nres; r++)
        {
            rtpr = &(db->rtp[db->rb.qrt[rt].subres[r].irtp]);
            /* rtpIndex = db->rb.res[rt][res].rtp; */

            /* Now map the atoms in the residue subtype (XXX)
             * to atoms in the residue type (qXXX) */
/* 	  natoms = rtpr->natom; */
/* 	  srenew(atomMap, natoms); */

/* 	  for (i=0; i < natoms; i++) */
/* 	    { */

            /* No need to map atoms here anymore. We have the iatomMap[]. */

/* 	      atomMap[i] = -1; */
/* 	      for (j=0; j<nqatoms; j++) */
/* 		{ */
/* 		  /\* Match them on a name basis *\/ */
/* 		  if (gmx_strcasecmp(*(rtpr->atomname[i]), */
/* 				     *(db->rtp[db->rb.rtp[rt]].atomname[j])) == 0) */
/* 		    { */
/* 		      /\* It's a match! *\/ */
/* 		      atomMap[i] = j; */
/* 		      break; */
/* 		    } */
/* 		  if (atomMap[i] == -1) */
/* 		    gmx_fatal(FARGS, "Could not map atom %s in %s onto the atoms in %s.", */
/* 			      *(rtpr->atomname[i]), */
/* 			      rtpr->resname, */
/* 			      db->rtp[db->rb.rtp[rt]].resname); */
/* 		} */
/* 	    } */

            /* Loop over all bonded interactions for this residue */
            for (bt=0; bt<ebtsNR; bt++)
            {
                for (bi=0; bi < rtpr->rb[bt].nb; bi++)
                {
                    /* Parse the params and atoms. */
                    str2bonded(rtpr->rb[bt].b[bi].s,
                               bt,
                               bt<=ebtsIDIHS ? db->bts[bt] : 0,
                               &params,
                               iatoms);

                    /* Store the functype/iparams index */
                    ilistEntry[0] = gmx_mtop_append_itype(mtop, ft, params);

                    /* find corresponding atoms in qXXX via iatomMap.
                     * note that the atoms in ilistEntry[1,...] are relative
                     * to the start of the residue, so they must be offset when
                     * actually adding this stuff to the ilist. */
                    for (i=0; i < btsNiatoms[bt]; i++)
                        ilistEntry[i+1] = db->rb.qrt[rt].subres[r].iatomMap[iatoms[i]];
                }
            } /* bt */
        } /* res */
    } /* res */
}

