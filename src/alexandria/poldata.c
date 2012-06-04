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
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include "network.h"
#include "string2.h"
#include "smalloc.h"
#include "names.h"
#include "gmx_fatal.h"
#include "poldata.h"
#include "gmx_simple_comm.h"

typedef struct {
    char   *desc,*type,*elem,*miller_equiv,*vdwparams,*charge;
    double polarizability,sig_pol;
} t_ffatype;

typedef struct {
    char   *elem,*name,*type,*neighbors,*geometry;
    char   **nb;
    int    numbonds,numnb;
} t_brule;

typedef struct {
    char   *atom1,*atom2,*params;
    char   elem1[4],elem2[4];
    double length,sigma,bondorder; 
} t_gt_bond;

typedef struct {
    char   *atom1,*atom2,*atom3,*params;
    double angle,sigma; 
} t_gt_angle;

typedef struct {
    char   *atom1,*atom2,*atom3,*atom4,*params;
    double dihedral,sigma; 
} t_gt_dihedral;

typedef struct {
    char   *elem;
    double polarizability;
} t_bosque;

typedef struct {
    char   *name,*alexandria_equiv;
    int    atomnumber;
    double tau_ahc,alpha_ahp;
} t_miller;

typedef struct {
    char *central;
    char *attached;
    int  numattach;
} t_symcharges;

#define EEMBUFSIZE 256
#define MAXZETA    12
typedef struct {
    int    eqg_model,nzeta,row[MAXZETA];
    char   name[EEMBUFSIZE],zetastr[EEMBUFSIZE],qstr[EEMBUFSIZE],rowstr[EEMBUFSIZE];
    double J0,chi0,q[MAXZETA],zeta[MAXZETA]; 
} t_eemprops;

typedef struct {
    int  eqg_model;
    char *epref;
} t_epref;

typedef struct gmx_poldata {
    char *filename;
    int nalexandria,nalexandria_c;
    t_ffatype *alexandria;
    int nbrule,nbrule_c;
    t_brule *brule;
    char *alexandria_polar_unit;
    char *alexandria_forcefield;
    int nexcl;
    double fudgeQQ,fudgeLJ;
    char *gt_vdw_function,*gt_combination_rule;
    int gt_vdw_ftype,gt_comb_rule;
    char *gt_bond_function;
    int ngt_bond,ngt_bond_c,gt_bond_ftype;
    char *gt_length_unit;
    t_gt_bond *gt_bond;
    char *gt_angle_function;
    int ngt_angle,ngt_angle_c,gt_angle_ftype;
    char *gt_angle_unit;
    t_gt_angle *gt_angle;
    char *gt_dihedral_function[egdNR];
    int ngt_dihedral[egdNR],ngt_dihedral_c[egdNR],gt_dihedral_ftype[egdNR];
    t_gt_dihedral *gt_dihedral[egdNR];
    int nmiller,nmiller_c;
    t_miller *miller;
    char *miller_tau_unit,*miller_ahp_unit;
    int nbosque,nbosque_c;
    t_bosque *bosque;
    char *bosque_polar_unit;
    int nsymcharges,nsymcharges_c;
    t_symcharges *symcharges;
    int nep,nep_c;
    t_eemprops *eep;
    int ner,ner_c;
    t_epref *epr;
} gmx_poldata;

#define assign_str(dst,src)  if (dst) { if (src) *dst = strdup(src); else *dst = NULL; }
#define assign_scal(dst,src) if (dst) *dst = src

gmx_poldata_t gmx_poldata_init()
{
    gmx_poldata_t pd;
    int i;
      
    snew(pd,1);
    
    /* Initiate some crucial variables */
    pd->nexcl = NOTSET;
    pd->fudgeQQ = NOTSET;
    pd->fudgeLJ = NOTSET;
    pd->gt_comb_rule = NOTSET;
    pd->gt_vdw_ftype = NOTSET;
    pd->gt_bond_ftype = NOTSET;
    pd->gt_angle_ftype = NOTSET;
    for(i=0; (i<egdNR); i++)
        pd->gt_dihedral_ftype[i] = NOTSET;
    
    return pd;
}

void gmx_poldata_set_filename(gmx_poldata_t pd,char *fn2)
{
    if (NULL == fn2) {
        fprintf(stderr,"Trying to set poldata filename to NULL\n");
        return;
    }
    if (NULL != pd->filename) {
        fprintf(stderr,"Overwriting poldata filename from %s to %s\n",
                pd->filename,fn2);
        sfree(pd->filename);
    }
    pd->filename = strdup(fn2);
}
    
void gmx_poldata_set_vdw_function(gmx_poldata_t pd,char *func)
{
    int i;
    
    if (NULL != pd->gt_vdw_function)
        sfree(pd->gt_vdw_function);
    for(i=0; (i<F_NRE); i++)
        if (strcasecmp(interaction_function[i].name,func) == 0)
            break;
    if (i == F_NRE)
        gmx_fatal(FARGS,"Van der Waals function '%s' does not exist in gromacs",func);
    pd->gt_vdw_ftype = i;
    pd->gt_vdw_function = strdup(func);
}

char *gmx_poldata_get_vdw_function(gmx_poldata_t pd)
{
    return pd->gt_vdw_function;
}

int gmx_poldata_get_vdw_ftype(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_vdw_ftype)
        gmx_fatal(FARGS,"No valid function type for Van der Waals function in %s",
                  pd->filename);
    return pd->gt_vdw_ftype;
}

void gmx_poldata_set_nexcl(gmx_poldata_t pd,int nexcl)
{
    pd->nexcl = nexcl;
}

int gmx_poldata_get_nexcl(gmx_poldata_t pd)
{
    if (NOTSET == pd->nexcl) 
        gmx_fatal(FARGS,"Nexclusions not set in %s",pd->filename);
    return pd->nexcl;
}

void gmx_poldata_set_fudgeQQ(gmx_poldata_t pd,double fudgeQQ)
{
    pd->fudgeQQ = fudgeQQ;
}

double gmx_poldata_get_fudgeQQ(gmx_poldata_t pd)
{
    if (NOTSET == pd->fudgeQQ) 
        gmx_fatal(FARGS,"fudgeQQ not set in %s",pd->filename);
    return pd->fudgeQQ;
}

void gmx_poldata_set_fudgeLJ(gmx_poldata_t pd,double fudgeLJ)
{
    pd->fudgeLJ = fudgeLJ;
}

double gmx_poldata_get_fudgeLJ(gmx_poldata_t pd)
{
    if (NOTSET == pd->fudgeLJ) 
        gmx_fatal(FARGS,"fudgeLJ not set in %s",pd->filename);
    return pd->fudgeLJ;
}

void gmx_poldata_set_combination_rule(gmx_poldata_t pd,char *func)
{
    int i;
    
    if (NULL != pd->gt_combination_rule)
        sfree(pd->gt_combination_rule);
    for(i=0; (i<eCOMB_NR); i++)
        if (strcasecmp(ecomb_names[i],func) == 0)
            break;
    if (i == eCOMB_NR)
        gmx_fatal(FARGS,"Combination rule '%s' does not exist in gromacs",func);
    pd->gt_comb_rule = i;
    pd->gt_combination_rule = strdup(func);
}

char *gmx_poldata_get_combination_rule(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_vdw_ftype)
        gmx_fatal(FARGS,"No valid function type for Van der Waals function in %s",
                  pd->filename);
    return pd->gt_combination_rule;
}

int gmx_poldata_get_comb_rule(gmx_poldata_t pd)
{
    if (NOTSET == pd->gt_comb_rule)
        gmx_fatal(FARGS,"No valid combination rule in %s",pd->filename);
    return pd->gt_comb_rule;
}

void gmx_poldata_add_atype(gmx_poldata_t pd,char *elem,char *desc,
                           char *gt_type,
                           char *miller_equiv,char *charge,
                           double polarizability,double sig_pol,
                           char *vdwparams)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    char    buf[EEMBUFSIZE],**ptr;
    int     i,j;
  
    for(i=0; (i<pold->nalexandria); i++) 
        if (strcasecmp(pold->alexandria[i].type,gt_type) == 0)
            break;
    if (i == pold->nalexandria) {
        pold->nalexandria++;
        srenew(pold->alexandria,pold->nalexandria);
  
        sp = &(pold->alexandria[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->type           = strdup(gt_type);
        sp->miller_equiv   = strdup(miller_equiv);
        sp->charge         = strdup(charge);
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
        sp->vdwparams      = strdup(vdwparams);
    }
    else 
        fprintf(stderr,"Atom type %s was already added to poldata record\n",gt_type);
}

void gmx_poldata_add_bonding_rule(gmx_poldata_t pd,
                                  char *gt_name,char *gt_type,
                                  char *geometry,int numbonds,
                                  char *neighbors)
{
    t_brule *sp;
    char    buf[EEMBUFSIZE],**ptr;
    int     i,j;
  
    for(j=0; (j<pd->nalexandria); j++) 
        if (strcasecmp(pd->alexandria[j].type,gt_type) == 0)
            break;
    if (j < pd->nalexandria) {
        for(i=0; (i<pd->nbrule); i++) 
            if (strcasecmp(pd->brule[i].name,gt_name) == 0)
                break;
        if (i == pd->nbrule) {
            pd->nbrule++;
            srenew(pd->brule,pd->nbrule);
            
            sp = &(pd->brule[i]);
            sp->elem           = strdup(pd->alexandria[j].elem);
            sp->name           = strdup(gt_name);
            sp->type           = strdup(gt_type);
            sp->neighbors      = strdup(neighbors);
            sp->nb             = split(' ',neighbors);
            for(sp->numnb=0; (sp->nb[sp->numnb] != NULL); sp->numnb++)
                ;
            sp->geometry       = strdup(geometry);
            sp->numbonds       = numbonds;
        }
        else 
            fprintf(stderr,"Atom %s was already added to pdata record\n",gt_name);
    }
    else if (NULL != debug) {
        fprintf(debug,"Ignoring bonding rule involving unknown atom type %s\n",
                gt_type);
    }
}

int gmx_poldata_get_bonding_rule(gmx_poldata_t pd,
                                 char **gt_name,char **gt_type,
                                 char **geometry,int *numbonds,
                                 char **neighbors)
{
    if (pd->nbrule_c < pd->nbrule) {
        assign_str(gt_name,pd->brule[pd->nbrule_c].name);
        assign_str(gt_type,pd->brule[pd->nbrule_c].type);
        assign_str(geometry,pd->brule[pd->nbrule_c].geometry);
        assign_scal(numbonds,pd->brule[pd->nbrule_c].numbonds);
        assign_str(neighbors,pd->brule[pd->nbrule_c].neighbors);
        
        pd->nbrule_c++;
        
        return 1;
    }
    return 0;
}

int gmx_poldata_get_natypes(gmx_poldata_t pd)
{
    return pd->nalexandria;
}

int gmx_poldata_get_ngt_bond(gmx_poldata_t pd)
{
    return pd->ngt_bond;
}

int gmx_poldata_get_ngt_angle(gmx_poldata_t pd)
{
    return pd->ngt_angle;
}

int gmx_poldata_get_ngt_dihedral(gmx_poldata_t pd,int egd)
{
    return pd->ngt_dihedral[egd];
}

void gmx_poldata_set_bond_function(gmx_poldata_t pd,char *fn)
{
    int i;
    
    if (NULL != pd->gt_bond_function)
        sfree(pd->gt_bond_function);
    for(i=0; (i<F_NRE); i++)
        if (strcasecmp(interaction_function[i].name,fn) == 0)
            break;
    if (i == F_NRE)
        gmx_fatal(FARGS,"Bond function '%s' does not exist in gromacs",fn);
    pd->gt_bond_ftype = i;
    pd->gt_bond_function = strdup(fn);
}

char *gmx_poldata_get_bond_function(gmx_poldata_t pd)
{
    return pd->gt_bond_function;
}

void gmx_poldata_set_angle_function(gmx_poldata_t pd,char *fn)
{
    int i;
    
    if (NULL != pd->gt_angle_function)
        sfree(pd->gt_angle_function);
    for(i=0; (i<F_NRE); i++)
        if (strcasecmp(interaction_function[i].name,fn) == 0)
            break;
    if (i == F_NRE)
        gmx_fatal(FARGS,"Angle function '%s' does not exist in gromacs",fn);
    pd->gt_angle_ftype = i;
    pd->gt_angle_function = strdup(fn);
}

char *gmx_poldata_get_angle_function(gmx_poldata_t pd)
{
    return pd->gt_angle_function;
}

int gmx_poldata_get_bond_ftype(gmx_poldata_t pd)
{
    return pd->gt_bond_ftype;
}

int gmx_poldata_get_angle_ftype(gmx_poldata_t pd)
{
    return pd->gt_angle_ftype;
}

int gmx_poldata_get_dihedral_ftype(gmx_poldata_t pd,int egd)
{
    return pd->gt_dihedral_ftype[egd];
}

void gmx_poldata_set_dihedral_function(gmx_poldata_t pd,int egd,char *func)
{
    int i;

    if (NULL != pd->gt_dihedral_function[egd])
        sfree(pd->gt_dihedral_function[egd]);
    for(i=0; (i<F_NRE); i++)
        if (strcasecmp(interaction_function[i].name,func) == 0)
            break;
    if (i == F_NRE)
        gmx_fatal(FARGS,"Dihedral function '%s' does not exist in gromacs",func);
    pd->gt_dihedral_ftype[egd] = i;
    pd->gt_dihedral_function[egd] = strdup(func);
}

char *gmx_poldata_get_dihedral_function(gmx_poldata_t pd,int egd)
{
    return pd->gt_dihedral_function[egd];
}

void gmx_poldata_set_atype_polarizability(gmx_poldata_t pd,char *gt_type,
                                          double polarizability,double sig_pol)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    int i;
  
    for(i=0; (i<pold->nalexandria); i++) {
        if (strcasecmp(gt_type,pold->alexandria[i].type) == 0) {
            sp = &(pold->alexandria[i]);
            sp->polarizability = polarizability;
            sp->sig_pol = sig_pol;
        }
    }
}		

char *gmx_poldata_get_polar_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->alexandria_polar_unit;
}

char *gmx_poldata_get_length_unit(gmx_poldata_t pd)
{
    if (NULL == pd->gt_length_unit)
        gmx_fatal(FARGS,"No length unit in %s",
                  (NULL != pd->filename) ? pd->filename : "unknown");

    return pd->gt_length_unit;
}

void gmx_poldata_set_polar_unit(gmx_poldata_t pd,char *polar_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->alexandria_polar_unit   = strdup(polar_unit);
}

char *gmx_poldata_get_force_field(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->alexandria_forcefield;
}

void gmx_poldata_set_force_field(gmx_poldata_t pd,char *forcefield)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->alexandria_forcefield   = strdup(forcefield);
}

void gmx_poldata_set_length_unit(gmx_poldata_t pd,char *length_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->gt_length_unit   = strdup(length_unit);
}

char *gmx_poldata_get_geometry(gmx_poldata_t pd,char *gt_atom)
{
    int i;
  
    if (gt_atom)
        for(i=0; (i<pd->nbrule); i++) 
            if (strcasecmp(pd->brule[i].name,gt_atom) == 0) 
                return pd->brule[i].geometry;
      
    return NULL;
}

char *gmx_poldata_get_desc(gmx_poldata_t pd,char *gt_type)
{
    int i;
  
    if (gt_type)
        for(i=0; (i<pd->nalexandria); i++) 
            if (strcasecmp(pd->alexandria[i].type,gt_type) == 0) 
                return pd->alexandria[i].desc;
      
    return NULL;
}

char *gmx_poldata_get_type(gmx_poldata_t pd,char *gt_atom)
{
    int i;
  
    if (gt_atom)
        for(i=0; (i<pd->nbrule); i++) 
            if (strcasecmp(pd->brule[i].name,gt_atom) == 0) 
                return pd->brule[i].type;
  
    return NULL;
}

static gmx_bool strcasestr_start(char *needle,char *haystack)
{
    char *ptr;
    
    ptr = strcasestr(haystack,needle);
    
    return (ptr == haystack);
}

static int count_neighbors(t_brule *brule,int nbond,char *nbhybrid[])
{
    int i,j,ni=0,*jj,i_found;
        
    snew(jj,nbond+1);
    for(i=0; (i<brule->numnb); i++)
    {
        i_found = 0;
        for(j=0; (j<nbond); j++) 
        {
            if (
                (NULL != nbhybrid[j]) &&
                (jj[j] == 0) &&
                (i_found == 0) &&
                strcasestr_start(brule->nb[i],nbhybrid[j])
                )
            {
                i_found = 1;
                jj[j] = j+1;
                ni++;
            }
        }
    }
    sfree(jj);

    return ni;
}

char *gmx_poldata_get_charge(gmx_poldata_t pd,char *gt_type)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (gt_type)
        for(i=0; (i<pold->nalexandria); i++) 
            if (strcasecmp(pold->alexandria[i].type,gt_type) == 0) 
                return pold->alexandria[i].charge;
      
    return NULL;
}

char *gmx_poldata_get_atom(gmx_poldata_t pd,char *elem,
                           int nbond,char *neighbors[],
                           const char *geometry)
{
    int i,nnb,pbest=-1,best=-1;
  
    for(i=0; (i<pd->nbrule); i++) 
    {
        nnb = count_neighbors(&(pd->brule[i]),nbond,neighbors);
        if ((strcasecmp(pd->brule[i].elem,elem) == 0) &&
            (strcasecmp(pd->brule[i].geometry,geometry) == 0) &&
            (nbond == pd->brule[i].numbonds) &&
            (nnb == pd->brule[i].numnb))
        {
            if ((NULL != debug) && (best != -1) && (nnb == pbest))
            {
                fprintf(debug,"poldata clash: %s and %s are equally well suited, nnb = %d. Using %s.\n",
                        pd->brule[best].name,pd->brule[i].name,nnb,
                        pd->brule[best].name);
            }
            else if (nnb > pbest)
            {
                pbest = nnb;
                best = i;
            }
        }
    }
    
    if (best != -1)
        return pd->brule[best].name;
    else
        return NULL;
}

int gmx_poldata_type_polarizability(gmx_poldata_t pd,char *gt_type,double *polar,
                                    double *sig_pol)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;

    for(i=0; (i<pold->nalexandria); i++) 
        if (strcasecmp(gt_type,pold->alexandria[i].type) == 0)
        {
            *polar   = pold->alexandria[i].polarizability;
            *sig_pol = pold->alexandria[i].sig_pol;
            return 1;
        }
    return 0;
}

char *gmx_poldata_get_miller_equiv(gmx_poldata_t pd,char *gt_type)
{
    t_ffatype *sp;
    int i;

    for(i=0; (i<pd->nalexandria); i++) 
        if (strcasecmp(gt_type,pd->alexandria[i].type) == 0) 
            return pd->alexandria[i].miller_equiv;
    return NULL;
}

int gmx_poldata_get_atype(gmx_poldata_t pd,char **elem,char **desc,
                          char **gt_type,char **miller_equiv,
                          char **charge,double *polarizability,double *sig_pol,
                          char **vdwparams)
{
    t_ffatype *sp;
    int i;
    
    if (pd->nalexandria_c < pd->nalexandria) {
        sp = &(pd->alexandria[pd->nalexandria_c]);
        assign_str(charge,sp->charge);
        assign_scal(polarizability,sp->polarizability);
        assign_scal(sig_pol,sp->sig_pol);
        assign_str(elem,sp->elem);
        assign_str(desc,sp->desc);
        assign_str(gt_type,sp->type);
        assign_str(miller_equiv,sp->miller_equiv);
        assign_str(vdwparams,sp->vdwparams);
        pd->nalexandria_c++;
        return 1;
    }
    else
        pd->nalexandria_c = 0;
    
    return 0;
}

int gmx_poldata_search_atype(gmx_poldata_t pd,char *key,
                             char **elem,char **desc,
                             char **gt_type,char **miller_equiv,
                             char **charge,double *polarizability,double *sig_pol,
                             char **vdwparams)
{
    t_ffatype *sp;
    int i;
  
    for(i=0; (i<pd->nalexandria); i++) 
        if (strcasecmp(key,pd->alexandria[i].type) == 0) 
            break;
      
    if (i<pd->nalexandria) {
        sp = &(pd->alexandria[i]);
        assign_str(charge,sp->charge);
        assign_scal(polarizability,sp->polarizability);
        assign_scal(sig_pol,sp->sig_pol);
        assign_str(elem,sp->elem);
        assign_str(desc,sp->desc);
        assign_str(gt_type,sp->type);
        assign_str(miller_equiv,sp->miller_equiv);
        assign_str(vdwparams,sp->vdwparams);
    
        return 1;
    }
    else
        return 0;
}
				  
double gmx_poldata_get_bondorder(gmx_poldata_t pd,char *elem1,char *elem2,
                                 double distance,double toler)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    char *ba1,*ba2;
    double dev,dev_best = 100000;
    int j,i,i_best=-1;
  
    if ((NULL == elem1) || (NULL == elem2))
        return 0.0;
    for(i=0; (i<pold->ngt_bond); i++) {
        if (0 == strlen(pold->gt_bond[i].elem1)) {
            for(j=0; (j<pold->nalexandria); j++) 
                if (strcasecmp(pold->alexandria[j].type,pold->gt_bond[i].atom2) == 0)
                    strcmp(pold->gt_bond[i].elem1,pold->alexandria[j].elem);
        }
        if (0 == strlen(pold->gt_bond[i].elem2)) {
            for(j=0; (j<pold->nalexandria); j++) 
                if (strcasecmp(pold->alexandria[j].type,pold->gt_bond[i].atom2) == 0)
                    strcmp(pold->gt_bond[i].elem2,pold->alexandria[j].elem);
        }
        ba1 = pold->gt_bond[i].elem1;
        ba2 = pold->gt_bond[i].elem2;
        if (((strcasecmp(ba1,elem1) == 0) && (strcasecmp(ba2,elem2) == 0)) ||
            ((strcasecmp(ba1,elem2) == 0) && (strcasecmp(ba2,elem1) == 0))) {
            dev = fabs((pold->gt_bond[i].length - distance)/pold->gt_bond[i].length);
            if (dev < dev_best) {
                dev_best = dev;
                i_best = i;
            }
        }
    }
    if (dev_best < toler)
        return pold->gt_bond[i_best].bondorder;
  
    return 0.0;
}
				  
void gmx_poldata_add_miller(gmx_poldata_t pd,char *name,
                            int atomnumber,
                            double tau_ahc,double alpha_ahp,
                            char *alexandria_equiv)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_miller *mil;
  
    pold->nmiller++;
    srenew(pold->miller,pold->nmiller);
    mil = &(pold->miller[pold->nmiller-1]);
    mil->name       = strdup(name);
    mil->atomnumber = atomnumber;
    mil->tau_ahc    = tau_ahc;
    mil->alpha_ahp  = alpha_ahp;
    if (alexandria_equiv)
        mil->alexandria_equiv = strdup(alexandria_equiv);
    else
        mil->alexandria_equiv = NULL;
}
				  
void gmx_poldata_set_miller_units(gmx_poldata_t pd,char *tau_unit,char *ahp_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->miller_tau_unit = strdup(tau_unit);
    pold->miller_ahp_unit = strdup(ahp_unit);
}

void gmx_poldata_get_miller_units(gmx_poldata_t pd,char **tau_unit,
                                  char **ahp_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    assign_str(tau_unit,pold->miller_tau_unit);
    assign_str(ahp_unit,pold->miller_ahp_unit);
}

char *gmx_poldata_get_miller(gmx_poldata_t pd,char *name,
                             int *atomnumber,double *tau_ahc,
                             double *alpha_ahp,char **alexandria_equiv)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_miller *mil;
    int i;
  
    if (name) {
        for(i=0; (i<pold->nmiller); i++) 
            if (strcasecmp(name,pold->miller[i].name) == 0) 
                break;
    }
    else
        i = pold->nmiller_c;
    
    if (i < pold->nmiller) {
        mil = &(pold->miller[i]);
        assign_scal(atomnumber,mil->atomnumber);
        assign_scal(tau_ahc,mil->tau_ahc);
        assign_scal(alpha_ahp,mil->alpha_ahp);
        assign_str(alexandria_equiv,mil->alexandria_equiv);
        if (!name)
            pold->nmiller_c++;
    
        return mil->name;
    }
    else
        pold->nmiller_c = 0;
    
    return NULL;
}

void gmx_poldata_add_bosque(gmx_poldata_t pd,char *elem,
                            double polarizability)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_bosque *bs;
  
    pold->nbosque++;
    srenew(pold->bosque,pold->nbosque);
    bs = &(pold->bosque[pold->nbosque-1]);
    bs->elem           = strdup(elem);
    bs->polarizability = polarizability;
}

char *gmx_poldata_get_bosque(gmx_poldata_t pd,char *name,char *elem,
                             double *polarizability)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if ((NULL != name) && 
        ((strchr(name,'+') != NULL) || (strchr(name,'-') != NULL)))
        return NULL;
    if (elem) {
        for(i=0; (i<pold->nbosque); i++) 
            if (strcasecmp(elem,pold->bosque[i].elem) == 0) 
                break;
    }
    else 
        i = pold->nbosque_c;
    
    if (i < pold->nbosque) {
        assign_scal(polarizability,pold->bosque[i].polarizability);
        if (!elem)
            pold->nbosque_c++;
    
        return pold->bosque[i].elem;
    }
    else
        pold->nbosque_c = 0;
    
    return NULL;
}
				  
void gmx_poldata_set_bosque_unit(gmx_poldata_t pd,char *polar_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->bosque_polar_unit   = strdup(polar_unit);
}
				  
char *gmx_poldata_get_bosque_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->bosque_polar_unit;
}

static int search_atomtype(gmx_poldata_t pd,char *atom)
{
    int j;
    
    for(j=0; (j<pd->nalexandria); j++) 
        if (strcasecmp(pd->alexandria[j].type,atom) == 0) 
            return j;
    return -1;
}

/* 
 * gt_bond stuff 
 */
int gmx_poldata_set_bond_params(gmx_poldata_t pd,char *atom1,char *atom2,
                                char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_bond *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_bond); i++) {
        gt_b = &(pold->gt_bond[i]);
        if ((((strcasecmp(gt_b->atom1,atom1) == 0) && 
              (strcasecmp(gt_b->atom2,atom2) == 0)) ||
             ((strcasecmp(gt_b->atom1,atom2) == 0) && 
              (strcasecmp(gt_b->atom2,atom1) == 0))))
            break;
    }
    if (i < pold->ngt_bond) {
        if (NULL != gt_b->params)
            sfree(gt_b->params);
        gt_b->params = strdup(params);
        return 1;
    }
    return 0;
}

int gmx_poldata_add_bond(gmx_poldata_t pd,char *atom1,char *atom2,
                         double length,double sigma,double bondorder,char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_bond *gt_b;
    int j,a1,a2;

    if (-1 == (a1 = search_atomtype(pd,atom1)))
        return 0;
    if (-1 == (a2 = search_atomtype(pd,atom2)))
        return 0;
    if (gmx_poldata_set_bond_params(pd,atom1,atom2,params) == 0) {
        pold->ngt_bond++;
        srenew(pold->gt_bond,pold->ngt_bond);
        gt_b = &(pold->gt_bond[pold->ngt_bond-1]);
        gt_b->atom1     = strdup(atom1);
        strncpy(gt_b->elem1,pold->alexandria[a1].elem,sizeof(gt_b->elem1)-1);
        gt_b->atom2     = strdup(atom2);
        strncpy(gt_b->elem2,pold->alexandria[a2].elem,sizeof(gt_b->elem2)-1);
        gt_b->length    = length;
        gt_b->sigma     = sigma;
        gt_b->bondorder = bondorder;
        gt_b->params    = strdup(params);
    }
    return 1;
}

int gmx_poldata_get_bond(gmx_poldata_t pd,char **atom1,char **atom2,
                         double *length,double *sigma,double *bondorder,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_bond *gt_b;
    int i;
  
    if (pold->ngt_bond_c < pold->ngt_bond) {
        gt_b = &(pold->gt_bond[pold->ngt_bond_c]);
        assign_str(atom1,gt_b->atom1);
        assign_str(atom2,gt_b->atom2);
        assign_scal(length,gt_b->length);
        assign_scal(sigma,gt_b->sigma);
        assign_scal(bondorder,gt_b->bondorder);
        assign_str(params,gt_b->params);
        pold->ngt_bond_c++;
    
        return pold->ngt_bond_c;
    }
    pold->ngt_bond_c = 0;
  
    return 0;
}

int gmx_poldata_search_bond(gmx_poldata_t pd,char *atom1,char *atom2,
                            double *length,double *sigma,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_bond *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_bond); i++) {
        gt_b = &(pold->gt_bond[i]);
        if (((strcasecmp(gt_b->atom1,atom1) == 0) && 
             (strcasecmp(gt_b->atom2,atom2) == 0)) ||
            ((strcasecmp(gt_b->atom1,atom2) == 0) && 
             (strcasecmp(gt_b->atom2,atom1) == 0))) {
            assign_scal(length,gt_b->length);
            assign_scal(sigma,gt_b->sigma);
            assign_str(params,gt_b->params);
        
            return i+1;
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
int gmx_poldata_set_angle_params(gmx_poldata_t pd,char *atom1,char *atom2,
                                 char *atom3,char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_angle *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_angle); i++) {
        gt_b = &(pold->gt_angle[i]);
        if ((strcasecmp(gt_b->atom2,atom2) == 0) && 
            (((strcasecmp(gt_b->atom1,atom1) == 0) && 
              (strcasecmp(gt_b->atom3,atom3) == 0)) ||
             ((strcasecmp(gt_b->atom1,atom3) == 0) && 
              (strcasecmp(gt_b->atom3,atom1) == 0))))
            break;
    }
    if (i < pold->ngt_angle) {
        if (NULL != gt_b->params)
            sfree(gt_b->params);
        gt_b->params = strdup(params);
        return 1;
    }
    return 0;
}
				   
int gmx_poldata_add_angle(gmx_poldata_t pd,char *atom1,char *atom2,
                          char *atom3,double angle,double sigma,char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_angle *gt_b;
    int i;
  
    if ((-1 == search_atomtype(pd,atom1)) ||
        (-1 == search_atomtype(pd,atom2)) ||
        (-1 == search_atomtype(pd,atom3)))
        return 0;

    if (0 == gmx_poldata_set_angle_params(pd,atom1,atom2,atom3,params)) {
        pold->ngt_angle++;
        srenew(pold->gt_angle,pold->ngt_angle);
        gt_b = &(pold->gt_angle[pold->ngt_angle-1]);
        gt_b->atom1   = strdup(atom1);
        gt_b->atom2   = strdup(atom2);
        gt_b->atom3   = strdup(atom3);
        gt_b->angle   = angle;
        gt_b->sigma   = sigma;
        gt_b->params  = strdup(params);
    }
    return 1;
}

int gmx_poldata_get_angle(gmx_poldata_t pd,char **atom1,char **atom2,
                          char **atom3,double *angle,double *sigma,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_angle *gt_b;
    int i;
  
    if (pold->ngt_angle_c < pold->ngt_angle) {
        gt_b = &(pold->gt_angle[pold->ngt_angle_c]);
        assign_str(atom1,gt_b->atom1);
        assign_str(atom2,gt_b->atom2);
        assign_str(atom3,gt_b->atom3);
        assign_scal(angle,gt_b->angle);
        assign_scal(sigma,gt_b->sigma);
        assign_str(params,gt_b->params);
        pold->ngt_angle_c++;
    
        return pold->ngt_angle_c;
    }
    pold->ngt_angle_c = 0;
  
    return 0;
}

int gmx_poldata_search_angle(gmx_poldata_t pd,char *atom1,char *atom2,
                             char *atom3,double *angle,double *sigma,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_angle *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_angle); i++) {
        gt_b = &(pold->gt_angle[i]);
        if ((strcasecmp(gt_b->atom2,atom2) == 0) &&
            (((strcasecmp(gt_b->atom1,atom1) == 0) && 
              (strcasecmp(gt_b->atom3,atom3) == 0)) ||
             ((strcasecmp(gt_b->atom1,atom3) == 0) && 
              (strcasecmp(gt_b->atom3,atom1) == 0)))) {
            assign_scal(angle,gt_b->angle);
            assign_scal(sigma,gt_b->sigma);
            assign_str(params,gt_b->params);
        
            return i+1;
        }
    }
    return 0;
}

void gmx_poldata_set_angle_unit(gmx_poldata_t pd,char *angle_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->gt_angle_unit   = strdup(angle_unit);
}

char *gmx_poldata_get_angle_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->gt_angle_unit;
}

/*
 * gt_dihedral stuff
 */
static int gtd_comp(const void *a,const void *b) 
{
    t_gt_dihedral *gt_a = (t_gt_dihedral *)a;
    t_gt_dihedral *gt_b = (t_gt_dihedral *)b;
    int n;
    
    if (0 == (n = strcasecmp(gt_a->atom1,gt_b->atom1)))
        if (0 == (n = strcasecmp(gt_a->atom2,gt_b->atom2)))
            if (0 == (n = strcasecmp(gt_a->atom3,gt_b->atom3)))
                n = strcasecmp(gt_a->atom4,gt_b->atom4);
                
    return n;
}   

static t_gt_dihedral *search_dihedral(gmx_poldata_t pd,int egd,
                                      char *atom1,char *atom2,
                                      char *atom3,char *atom4)
{
    t_gt_dihedral gt_a,*gt_res,*gt_dptr;
    int nd;
    
    gt_dptr = pd->gt_dihedral[egd];
    nd = pd->ngt_dihedral[egd];
    gt_a.atom1 = atom1;
    gt_a.atom2 = atom2;
    gt_a.atom3 = atom3;
    gt_a.atom4 = atom4;
    gt_res = bsearch(&gt_a,gt_dptr,nd,sizeof(gt_a),&gtd_comp);
    if (NULL == gt_res) {
        gt_a.atom1 = atom4;
        gt_a.atom2 = atom3;
        gt_a.atom3 = atom2;
        gt_a.atom4 = atom1;
        gt_res = bsearch(&gt_a,gt_dptr,nd,sizeof(gt_a),gtd_comp);
    }
    return gt_res;
}

int gmx_poldata_set_dihedral_params(gmx_poldata_t pd,int egd,
                                    char *atom1,char *atom2,
                                    char *atom3,char *atom4,char *params)
{
    t_gt_dihedral *gt_res;
    
    gt_res = search_dihedral(pd,egd,atom1,atom2,atom3,atom4);
    if (NULL != gt_res) {
        if (NULL != gt_res->params)
            sfree(gt_res->params);
        gt_res->params = strdup(params);
        return 1;
    }
    return 0;
}

int gmx_poldata_add_dihedral(gmx_poldata_t pd,int egd,
                             char *atom1,char *atom2,
                             char *atom3,char *atom4,double dihedral,
                             double sigma,char *params)
{
    t_gt_dihedral *gt_b;
    int i;

    if ((-1 == search_atomtype(pd,atom1)) ||
        (-1 == search_atomtype(pd,atom2)) ||
        (-1 == search_atomtype(pd,atom3)) ||
        (-1 == search_atomtype(pd,atom4)))
        return 0;

    if (0 == gmx_poldata_set_dihedral_params(pd,egd,atom1,atom2,
                                             atom3,atom4,params)) {  
        pd->ngt_dihedral[egd]++;
        srenew(pd->gt_dihedral[egd],pd->ngt_dihedral[egd]);
        gt_b = &(pd->gt_dihedral[egd][pd->ngt_dihedral[egd]-1]);
        gt_b->atom1   = strdup(atom1);
        gt_b->atom2   = strdup(atom2);
        gt_b->atom3   = strdup(atom3);
        gt_b->atom4   = strdup(atom4);
        gt_b->dihedral   = dihedral;
        gt_b->sigma   = sigma;
        gt_b->params = strdup(params);
        qsort(pd->gt_dihedral[egd],pd->ngt_dihedral[egd],sizeof(pd->gt_dihedral[egd][0]),
              gtd_comp);
    }
    return 1;
}

int gmx_poldata_get_dihedral(gmx_poldata_t pd,int egd,
                             char **atom1,char **atom2,
                             char **atom3,char **atom4,double *dihedral,
                             double *sigma,char **params)
{
    t_gt_dihedral *gt_b;
    int i;
  
    if (pd->ngt_dihedral_c[egd] < pd->ngt_dihedral[egd]) {
        gt_b = &(pd->gt_dihedral[egd][pd->ngt_dihedral_c[egd]]);
        assign_str(atom1,gt_b->atom1);
        assign_str(atom2,gt_b->atom2);
        assign_str(atom3,gt_b->atom3);
        assign_str(atom4,gt_b->atom4);
        assign_scal(dihedral,gt_b->dihedral);
        assign_scal(sigma,gt_b->sigma);
        assign_str(params,gt_b->params);
        pd->ngt_dihedral_c[egd]++;
    
        return pd->ngt_dihedral_c[egd];
    }
    pd->ngt_dihedral_c[egd] = 0;
  
    return 0;
}

int gmx_poldata_search_dihedral(gmx_poldata_t pd,int egd,
                                char *atom1,char *atom2,
                                char *atom3,char *atom4,
                                double *dihedral,double *sigma,char **params)
{
    t_gt_dihedral *gt_res;
    
    gt_res = search_dihedral(pd,egd,atom1,atom2,atom3,atom4);
    if (NULL != gt_res) {
        assign_scal(dihedral,gt_res->dihedral);
        assign_scal(sigma,gt_res->sigma);
        assign_str(params,gt_res->params);
        
        return 1 + (int) (gt_res - pd->gt_dihedral[egd]);
    }
    return 0;
}

void gmx_poldata_add_symcharges(gmx_poldata_t pd,char *central,
                                char *attached,int numattach)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_symcharges *sc;
    int i;
  
    for(i=0; (i<pold->nsymcharges); i++) {
        sc = &(pold->symcharges[i]);
        if ((strcasecmp(sc->central,central) == 0) &&
            (strcasecmp(sc->attached,attached) == 0) &&
            (sc->numattach == numattach))
            break;
    }
    if (i == pold->nsymcharges) {
        pold->nsymcharges++;
        srenew(pold->symcharges,pold->nsymcharges);
        sc = &(pold->symcharges[i]);
        sc->central  = strdup(central);
        sc->attached = strdup(attached);
        sc->numattach   = numattach;
    }
}

int gmx_poldata_get_symcharges(gmx_poldata_t pd,char **central,
                               char **attached,int *numattach)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_symcharges *sc;
    int i;
  
    if (pold->nsymcharges_c < pold->nsymcharges) {
        sc = &(pold->symcharges[pold->nsymcharges_c]);
        assign_str(central,sc->central);
        assign_str(attached,sc->attached);
        assign_scal(numattach,sc->numattach);
        pold->nsymcharges_c++;
    
        return 1;
    }
    pold->nsymcharges_c = 0;
  
    return 0;
}

int gmx_poldata_search_symcharges(gmx_poldata_t pd,char *central,
                                  char *attached,int numattach)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_symcharges *sc;
    int i;
  
    for(i=0; (i<pold->nsymcharges); i++) {
        sc = &(pold->symcharges[i]);
        if ((strcasecmp(sc->central,central) == 0) &&
            (strcasecmp(sc->attached,attached) == 0) &&
            (sc->numattach == numattach))
            return 1;
    }
  
    return 0;
}

/* Electrostatics properties */
static t_eemprops *get_eep(gmx_poldata *pold,int eqg_model,char *name) 
{
    int i,j;
    t_eemprops *eep;
  
    for(i=0; (i<pold->nep); i++) 
        if ((strcasecmp(pold->eep[i].name,name) == 0) &&
            (pold->eep[i].eqg_model == eqg_model))
            return &(pold->eep[i]);
    return NULL;
}

void gmx_poldata_set_eemprops(gmx_poldata_t pd,
                              int eqg_model,char *name,
                              double J0,double chi0,char *zeta,char *q,char *row)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_eemprops *eep;
    char **sz,**sq,**sr;
    int  n;
    
    eep = get_eep(pold,eqg_model,name);
    if (NULL == eep) {
        srenew(pold->eep,++pold->nep);
        eep = &(pold->eep[pold->nep-1]);
    }
    eep->eqg_model = eqg_model;
    strncpy(eep->name,name,EEMBUFSIZE-1);
    eep->name[EEMBUFSIZE-1] = '\0';
    eep->J0    = J0;
    sz = split(' ',zeta);
    sq = split(' ',q);
    sr = split(' ',row);
    strncpy(eep->zetastr,zeta,EEMBUFSIZE-1);
    strncpy(eep->qstr,q,EEMBUFSIZE-1);
    strncpy(eep->rowstr,row,EEMBUFSIZE-1);
    n = 0;
    while ((NULL != sz[n]) && (NULL != sq[n]) && (NULL != sr[n]))
    {
        if (n < MAXZETA) 
        {
            eep->zeta[n] = atof(sz[n]);
            eep->q[n]    = atof(sq[n]);
            eep->row[n]  = atoi(sr[n]);
        }
        sfree(sz[n]);
        sfree(sq[n]);
        sfree(sr[n]);
        n++;
    }
    if ((NULL != sz[n]) || (NULL != sq[n]) || (NULL != sr[n])) 
    {
        if (NULL != sz[n])
            fprintf(stderr,"Warning: more zeta values than q/row values for %s n = %d\n",name,n);
        if (NULL != sq[n])
            fprintf(stderr,"Warning: more q values than zeta/row values for %s n = %d\n",name,n);
        if (NULL != sr[n])
            fprintf(stderr,"Warning: more row values than q/zeta values for %s n = %d\n",name,n);
    }
    eep->nzeta = n;
    if (n >= MAXZETA)
    {
        fprintf(stderr,"More than %d zeta and/or q values for %s\n",MAXZETA,eep->name);
        eep->nzeta = MAXZETA;
    }
    for(; (n < MAXZETA); n++) 
    {
        eep->zeta[n] = 0;
        eep->q[n]    = 0;
        eep->row[n]  = 0;
    }
    eep->chi0  = chi0;
}

int gmx_poldata_get_eemprops(gmx_poldata_t pd,
                             int *eqg_model,char **name,
                             double *J0,double *chi0,char **zeta,char **q,char **row)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    if (pold->nep_c < pold->nep) {
        assign_scal(eqg_model,pold->eep[pold->nep_c].eqg_model);
        assign_str(name,pold->eep[pold->nep_c].name);
        assign_scal(J0,pold->eep[pold->nep_c].J0);
        assign_str(zeta,pold->eep[pold->nep_c].zetastr);
        assign_str(q,pold->eep[pold->nep_c].qstr);
        assign_str(row,pold->eep[pold->nep_c].rowstr);
        assign_scal(chi0,pold->eep[pold->nep_c].chi0);
        pold->nep_c++;
        return 1;
    }
    else {
        pold->nep_c = 0;
        return 0;
    }
}

int gmx_poldata_get_numprops(gmx_poldata_t pd,int eqg_model)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i,n=0;
  
    for(i=0; (i<pold->nep); i++) 
        if (pold->eep[i].eqg_model == eqg_model)
            n++;
      
    return n;
}

int gmx_poldata_have_pol_support(gmx_poldata_t pd,char *gt_type)
{
    int i;
    
    for(i=0; (i<pd->nalexandria); i++)
        if (strcasecmp(gt_type,pd->alexandria[i].type) == 0)
            return 1;
    return 0;
}

int gmx_poldata_have_eem_support(gmx_poldata_t pd,int eqg_model,char *name,
                                 gmx_bool bAllowZeroParameters)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_eemprops *eep = get_eep(pold,eqg_model,name);
  
    return (eep && (bAllowZeroParameters || ((eep->J0 > 0) && (eep->chi0 > 0))));
}

double gmx_poldata_get_j00(gmx_poldata_t pd,int eqg_model,char *name)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_eemprops *eer;
  
    if ((eer = get_eep(pold,eqg_model,name)) != NULL)
        return eer->J0;
    else
        gmx_fatal(FARGS,"No J0 data for eqg_model %d and name %s",
                  eqg_model,name);
    return -1;
}

char *gmx_poldata_get_qstr(gmx_poldata_t pd,int eqg_model,char *name)
{
    t_eemprops *eer;
    
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        return eer->qstr;
    }
    return NULL;
}

char *gmx_poldata_get_rowstr(gmx_poldata_t pd,int eqg_model,char *name)
{
    t_eemprops *eer;
    
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        return eer->rowstr;
    }
    return NULL;
}

int gmx_poldata_get_row(gmx_poldata_t pd,int eqg_model,char *name,int zz)
{
    t_eemprops *eer;
  
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        range_check(zz,0,eer->nzeta);
        return eer->row[zz];
    }
    return -1;
}

double gmx_poldata_get_zeta(gmx_poldata_t pd,int eqg_model,char *name,int zz)
{
    t_eemprops *eer;
  
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        if ((zz < 0) || (zz >= eer->nzeta))
            printf("Bleh\n");
        range_check(zz,0,eer->nzeta);
        return eer->zeta[zz];
    }
    return -1;
}

int gmx_poldata_get_nzeta(gmx_poldata_t pd,int eqg_model,char *name)
{
    t_eemprops *eer;
  
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        return eer->nzeta;
    }
    return 0;
}

double gmx_poldata_get_q(gmx_poldata_t pd,int eqg_model,char *name,int zz)
{
    t_eemprops *eer;
  
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
    {
        range_check(zz,0,eer->nzeta);
        return eer->q[zz];
    }
    return -1;
}

double gmx_poldata_get_chi0(gmx_poldata_t pd,int eqg_model,char *name)
{
    t_eemprops *eer;
    
    if ((eer = get_eep(pd,eqg_model,name)) != NULL)
        return eer->chi0;
    else
        gmx_fatal(FARGS,"No chi0 data for eqg_model %d and name %s",eqg_model,name);
    return -1;
}

void gmx_poldata_set_epref(gmx_poldata_t pd,int eqg_model,char *epref)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    for(i=0; (i<pold->ner); i++) {
        if (pold->epr[i].eqg_model == eqg_model)
        {
            if (pold->epr[i].epref)
                sfree(pold->epr[i].epref);
            pold->epr[i].epref = strdup(epref);
            break;
        }
    }
    if (i == pold->ner) 
    {
        srenew(pold->epr,++pold->ner);
        pold->epr[i].eqg_model = eqg_model;
        pold->epr[i].epref = strdup(epref);
    }
}

char *gmx_poldata_get_epref(gmx_poldata_t pd,int eqg_model)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    for(i=0; (i<pold->ner); i++) 
        if (pold->epr[i].eqg_model == eqg_model)
            return pold->epr[i].epref;
    return NULL;
}

int gmx_poldata_list_epref(gmx_poldata_t pd,int *eqg_model,char **epref)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (pold->ner_c < pold->ner)
    {
        assign_scal(eqg_model,pold->epr[pold->ner_c].eqg_model);
        assign_str(epref,pold->epr[pold->ner_c].epref);
        pold->ner_c++;
        return 1;
    }
    pold->ner_c = 0;

    return 0;
}

void gmx_poldata_comm_eemprops(gmx_poldata_t pd,t_commrec *cr)
{
    int i,j,nep;
    t_eemprops *ep;
    
    if (NULL != debug)
        fprintf(debug,"Going to update eemprops on node %d\n",cr->nodeid);
    if (MASTER(cr))
    {
        for(i=1; (i<cr->nnodes); i++) 
        {
            gmx_send_int(cr,i,pd->nep);
            gmx_send(cr,i,pd->eep,pd->nep*sizeof(pd->eep[0]));
        }
    }
    else 
    {
        nep = gmx_recv_int(cr,0);
        if (nep != pd->nep)
            gmx_fatal(FARGS,"Inconsistency in number of EEM parameters");
        snew(ep,pd->nep);
        gmx_recv(cr,0,ep,pd->nep*sizeof(ep[0]));
        for(i=0; (i<pd->nep); i++) 
            pd->eep[i] = ep[i];
        sfree(ep);
    }
    if (NULL != debug) 
    {
        fprintf(debug,"  EEP  Atom      Chi      J00     Zeta\n");
        for(i=0; (i<pd->nep); i++)
        {
            fprintf(debug,"%5s %5s %8.3f %8.3f",
                    get_eemtype_name(pd->eep[i].eqg_model),
                    pd->eep[i].name,pd->eep[i].chi0,pd->eep[i].J0);
            for(j=0; (j<pd->eep[i].nzeta); j++)
                fprintf(debug," %8.3f",pd->eep[i].zeta[j]);
            fprintf(debug,"\n");
        }
    }
}

void gmx_poldata_comm_force_parameters(gmx_poldata_t pd,t_commrec *cr)
{
    int i,j,nep;
    t_eemprops *ep;
    
    if (NULL != debug)
        fprintf(debug,"Going to update force parameters on node %d\n",cr->nodeid);
    if (MASTER(cr))
    {
        for(i=1; (i<cr->nnodes); i++) 
        {
            gmx_send_int(cr,i,pd->nep);
            gmx_send(cr,i,pd->eep,pd->nep*sizeof(pd->eep[0]));
        }
    }
    else 
    {
        nep = gmx_recv_int(cr,0);
        if (nep != pd->nep)
            gmx_fatal(FARGS,"Inconsistency in number of EEM parameters");
        snew(ep,pd->nep);
        gmx_recv(cr,0,ep,pd->nep*sizeof(ep[0]));
        for(i=0; (i<pd->nep); i++) 
            pd->eep[i] = ep[i];
        sfree(ep);
    }
    if (NULL != debug) 
    {
        fprintf(debug,"  EEP  Atom      Chi      J00     Zeta\n");
        for(i=0; (i<pd->nep); i++)
        {
            fprintf(debug,"%5s %5s %8.3f %8.3f",
                    get_eemtype_name(pd->eep[i].eqg_model),
                    pd->eep[i].name,pd->eep[i].chi0,pd->eep[i].J0);
            for(j=0; (j<pd->eep[i].nzeta); j++)
                fprintf(debug," %8.3f",pd->eep[i].zeta[j]);
            fprintf(debug,"\n");
        }
    }
}

void gmx_poldata_check_consistency(FILE *fp,gmx_poldata pd)
{
    
}

typedef struct {
    const char *name,*ref;
    gmx_bool bWeight;
} t_eemtype_props;

static t_eemtype_props eemtype_props[eqgNR] = { 
    { "None",     "None",	         FALSE },   
    { "Yang", 	"Yang2006b",     TRUE },    
    { "Bultinck", "Bultinck2002a", FALSE },
    { "Rappe",    "Rappe1991a",	 TRUE },  
    { "AXp",    	"Maaren2010a",	 FALSE },    
    { "AXs", 	    "Maaren2010a",   TRUE },
    { "AXg",      "Maaren2010a",   TRUE },
    { "ESP",      "Kollman1991a",  FALSE },
    { "RESP",     "Kollman1991a",  FALSE }
};

int name2eemtype(const char *name)
{
    int i;
  
    for(i=0; (i<eqgNR); i++) {
        if (strcasecmp(name,eemtype_props[i].name) == 0)
            return i;
    }
    return -1;
}

const char *get_eemtype_name(int eem)
{
    int i;
  
    if ((eem >= 0) && (eem < eqgNR))
        return eemtype_props[eem].name;
    
    return NULL;
}

