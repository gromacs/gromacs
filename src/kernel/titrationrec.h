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
#ifndef _TITRATIONREC_H
#define _TITRATIONREC_H

#include "resall.h"
#include "gpp_atomtype.h"
#include "gmx_random.h"
#include "constr.h"

/* -------      Nomenclature      ------- *
 * All types with the suffix _t are       *
 * pointers to types with the same        *
 * base name. Types without the _t suffix *
 * are structures or basic types.         *
 * -------------------------------------- */

enum { etQhopNONE, etQhopSE, etQhopI, etQhopTST, etQhopNR };

extern const char *qhopregimes[etQhopNR];

enum {eQNONE=0, eQACC=1, eQDON=2, eQACCDON=3, eQNR=4};
/* Trivial. But note that eQACCDON = eQACC | eQDON */

enum {eQIDEF, eQBONLIB}; /* use the first to get the pointer to the idef,
			    the other for the bonded stuff stored for a
			    certain residue. */

/* Symbolic name for qhop_db. We can't include
 * gmx_qhop_types.h since that breaks compilation. */
typedef struct qhop_db *qhop_db_s; 

typedef struct {
  int donor_id, acceptor_id, proton_id,
    regime, primary_d, primary_a;     
  /* These is used for proton tautomerism
     if primary != donor/acceptor then
     we swap the coordinates of the primary
     and donor/acceptor when titrating. */
  real E12_0, E12, DE_Environment, Eb,
    rda, ang, prob, hbo, kappa,
    Er, El, T;
  /* We put temperature T here, partially because one might
     want local temperatures in the future, but mainly to
     reduce the number of function arguments passed around
     in this more "object oriented" approach. The latter
     goes for most of the real data members here. */
  gmx_bool bFlip,bDonated;
    rvec xold,vold,   /* Where the acceptor proton used to be */
        xnew;    /* Where the acceptor proton may end up */
    rvec *f; /* We need to store the forces from each trial hop in case we do more than one */
} t_hop;

/* Keeps track of where bonded interactions are in an ilist.
   We use it to quickly change bonded interactions for a residue
   upon (de)protonation. */
typedef struct {
  t_iatom **ilist_pos; /* dimensions: [F_NRE][#bonded_of this type] */

  int nr[F_NRE]; /* How many interactions in ilist_pos[] arrays? */

  gmx_bool **indexed; /* Are the bonded interactions indexed.
		       * Dimensions: [F_NRE][#bonded_of this type] */

} qhop_bonded_index;

typedef struct {
  int atom_id; /* global atom index */
  int res_id;  /* global residue index */
  int qres_id; /* qhop_residue index */

  char *atomname, *resname; /* will point to names in the depths of qhop_resblocks */

  int state; /* None, acceptor, donor or both? This can perhaps speed up the search for potential hops.
	      * Takes on values eQNONE, eQACC, eQDON, eQACCDON. */
  
  int nr_protons; /* the number or protons that can hop away */
  int *protons;  /* the atom_id of the protons that can hop. These
		    are simply searched upon initialization. It is
	            assumed that we use complete topologies,
	            i.e. with all protons present. Protons are on/off
	            depending on their charges.

		    No they are not anymore. There is an exisence array for all hydrogens in the system
		    in qhop_db (maybe it should move to t_qhoprec?).
		  */

  int nr_acceptors; /* known after nbsearching. */
  int *acceptors;   /* j particles that fulfil additional geometric
		      criteria */
  /* upon accepting a proton, state becomes eQDON or eQACCDON; */
} t_qhop_atom;

typedef struct {
  /* rtype and res are needed to switch between interaction parameters. */
  int    rtype; /* index in qhop_db.restype[]  and qhop_db.res[]. Will not change upon (de)protonation */
  int    subres;   /* index in qhop_db.res[rtype][]. Changes upon (de)protonation. */

  int    *atoms; /* global atom numbers belonging to the residue */
  int    nr_atoms;
  char   **atomnames;
  int    nr_titrating_sites;
  int    *titrating_sites; /* points back to acceptor donors in the qhop_atoms structure. */

  int    state; /* 0 is competely deprotonated, 1 is fully protonated.
		   For histidines, we have to decide the state from
		   the atom name.... */
  /* Let's change to eQACC for deprotonated, eQACCDON for ampholytic and eQDON for fully protonated. */
  
  int    res_nr; /* Global resnr. */

  qhop_bonded_index bindex; /* Makes ilistPosG and ilistPosL redundant. For now anyway. */
  int nr_indexed;      /* are the ilists indexed? */
} t_qhop_residue;

typedef struct titration {
  gmx_constr_t   constr; 
  int            nr_hop,max_nr_hop;
  t_hop          *hop;
  int            nr_qhop_atoms;
  t_qhop_atom    *qhop_atoms;
  int            nr_qhop_residues;
  t_qhop_residue *qhop_residues;
  qhop_db_s      db;
  gmx_rng_t      rng, rng_int;
  gmx_bool       bFreshNlists, bHaveEbefore;
  int *global_atom_to_qhop_atom;
} titration;

/* I tried to strip the gmx_ suffix from the type names 
 * whenever possible, to keep the names a bit shorter. */

typedef struct qhop_db *qhop_db_t;
typedef struct qhop *qhop_t;
typedef struct qhop_resblocks *qhop_resblocks_t;
typedef struct qhop_subres *qhop_subres_t;
typedef struct qhop_interactions *qhop_interactions_t;
typedef struct qhop_H_exist *qhop_H_exist_t;

typedef struct qhop_H_exist {
  int     nH;   /* Number of hydrogens. Length of H. */
  char    *H;   /* existence function for all hydrogens.
                   I chose char since it's compact. */
  atom_id *H2atomid;   /* maps the elements in H to atom_id */
  int     *atomid2H;   /* the inverse function */
} qhop_H_exist;

/********************
 * From gmx_qhop_db *
 *                  */

typedef struct qhop_parameters {
  real  alpha, beta, gamma;
  real  k_1, k_2, k_3, m_1, m_2, m_3;
  real  s_A, t_A, v_A, s_B, s_C, t_C, v_C;
  real  f, g, h;
  real  p_1, q_1, q_2, q_3, r_1, r_2, r_3;
} qhop_parameters;

typedef struct qhop_resinfo_t {
  int id,charge,natom;
  int ndonor,*donor,nacceptor,*acceptor;
} qhop_resinfo_t;

/**********************
 * From gmx_qhop_parm *
 *                    */
	
typedef struct qhop {
  char *donor,*acceptor, *don_atom, *acc_atom;
  int nparam,nparam_c;
  char **value,**unit,**name;
} qhop;

typedef struct qhop_reactant {
  int nname;      /* Length of **name */
  char **name;    /* A list of acceptor/donor atoms, due to proton 
		     tautomerism, eg. the two oxygens in a carbonyl. */
  char *product;  /* What will the res turn into when this 
		     donor/aceptor reacts? */
  qhop_subres_t productdata; /* Pointer to the product qhop_subres */
  int  nH;        /* Number of protons */
  char **H    ;   /* Proton names */
} qhop_reactant;

typedef struct qhop_subres {
  char *name;    /* Name of the residue */
  int na, nd;    /* Number of acceptors and donors */
  qhop_reactant *acc, *don;
  int *iatomMap; /* maps the atoms to the atoms in restype qXXX */
  int **biMap;   /* maps the bonded interactions to the bondeds in restype qXXX
		  * Matches the rtp data. */
  int niatom;    /* = db->rtp[rtp].natom */
  /*   int nft; */
  int irtp;      /* indexes the t_restp-array rtp in qhop_db. */

  /* Here are arrays that indexes functypes/parameters.
   * it must match the rtp data.
   * Since the ilib is attached to the end of the existing
   * t_iparams* in the local topology, findex is indirectly
   * an index to the parameters in the local topology:
   * ilib[i] = top->idef.iparams[i + top->idef.ntypes - db->rb.ni]

   * findex[bt][b] */
  int       **findex; /* indexes the parameters in the ilib...  */
  int      **mfindex; /* ... and the molecular topology.        */
} qhop_subres;

typedef struct {
    char *canonical;     /* Name of the "residue family", e.g. ASP */
    char *protonated;    /* Name of the fully protonated variant, e.g. ASH */
  char *description;   /* Descriptive comment */
  int  nsubres;        /* Number of related residues, e.g. 2 for 
			  qLYS: {LYS, LYSH} */
  qhop_subres *subres; /* has size [nsubres[i] */
  gmx_bool bWater;     /* Is this a water? */
  gmx_bool bInTop;     /* Is this present in the topology? */
  int irtp;            /* The root of all evil:
			  indexes the t_restp-array rtp in qhop_db. One element 
			  for each restype. Note that this is for the "residue 
			  families" only, e.g. qLYS. Every related residue 
			  has its own index to the t_restp-array in 
			  qhop_reactant. */
  int **ba[ebtsNR];  /* Residue local atom numbers for the bonded interactions
		      * first dimension == ebtsNR
		      * second dimension  == bond/angle/dihedral
		      * Matches the bonded interactions in the rtp-data */
} qhop_restype;

typedef struct qhop_resblocks {
  int nrestypes;     /* Number of restypes */
  qhop_restype *qrt; /* The actual restypes */

  /* The following is for the interaction library *
   * It stores the interaction parameters for a residue.
   * They are needed to switch between protonation states. */

  int  nf;           /* number of extra files */
  char **files;      /* extra files containg additional parameters. */


  int btype[ebtsNR]; /* which forcetype correpond to the bonded interactions
		      * in the rtp data? E.g. F_ANGLE ... */
  int ftype[ebtsNR]; /* Which functype correpond to the bonded interactions
		      * in the rtp data? E.g. dihedral type 9 */

  int ni;            /* Size of ilib below. */

  t_iparams *ilib;   /* The actual interaction parameters.
		      * ilib will be appended to the iparams
		      * in a gmx_localtop_t.idef */
  t_functype *ftlib; /* Functypes for all parameters in ilib. */
  int inull[ebtsNR]; /* Here are the the dummy interactions in ilib. 
		      * We need one per bonded type since the functypes
		      * differ for the dummy interactions although their
		      * parameters are the same. */
} qhop_resblocks;

typedef struct currentRes {
  int Acc;  /* Is it an acceptor? If Acc==exmlACCEPTOR then yes. */
  int rt;   /* Finds the current restype in qhop_resblocks */
  int r;    /* Finds the current res in in qhop_resblocks[rt] */
  int da;   /* Finds the current don or acc in qhop_resblocks[rt][r] */
} currentRes;

typedef struct qhop_db {
  int              inertH;    /* The atomtype for the inert hydrogen. */
  gmx_bool         bConstraints; /* True if constraints are present */
  int              nrtp;
  t_restp          *rtp;
  /* Replacing resinfo with more elaborate structures */
  /*qhop_resinfo_t    *resinfo; */
  int              bts[4];
  int              nrexcl;
  gmx_bool         bAllDih,bHH14,bRemoveDih;
  gpp_atomtype_t   atype;
  t_symtab         tab;
  int              ngqh;
  qhop_t           *gqh;
  qhop_resblocks   rb;
  qhop_parameters  *qhop_param;
  int              nres;
  qhop_H_exist     H_map;
  t_idef           *idef; /* Must point to the idef used by the integrator. */
} qhop_db;

#endif
