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
#ifndef _GMX_QHOP_TYPES_H
#define _GMX_QHOP_TYPES_H

#include "resall.h"
#include "gpp_atomtype.h"

/* -------      Nomenclature      ------- *
 * All types with the suffix _t are       *
 * pointers to types with the same        *
 * base name. Types without the _t suffix *
 * are structures or basic types.         *
 * -------------------------------------- */

/* I tried to strip the gmx_ suffix from the type names 
 * whenever possible, to keep the names a bit shorter. */

typedef struct qhop_db *qhop_db_t;

typedef struct qhop *qhop_t;
typedef struct qhop_resblocks *qhop_resblocks_t;
typedef struct qhop_subres *qhop_subres_t;
typedef struct qhop_interactions *qhop_interactions_t;
typedef struct qhop_H_exist *qhop_H_exist_t;

enum {etQhopNONE, etQhopSE, etQhopI, etQhopTST};

static const char *qhopregimes[] = {"NONE", "SE", "Intermediate", "TST"};

typedef struct qhop_H_exist {
  int     nH;   /* Number of hydrogens. Length of H. */
  char    *H;   /* existence function for all hydrogens.
	          * I chose char since it's compact. */
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
  char *restype;       /* Name of the "residue family", e.g. qLYS. */
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
  int              constrain; /* What do we constrain? none=0, h-bonds=1, 
				 all-bonds=2, */
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
