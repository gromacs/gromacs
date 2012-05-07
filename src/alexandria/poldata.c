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
#include "gmx_fatal.h"
#include "poldata.h"
#include "gmx_simple_comm.h"

typedef struct {
    char   *desc,*name,*type,*elem,*miller_equiv,*neighbors,*geometry,*vdwparams,*charge;
    char   **nb;
    int    numbonds,numnb;
    double polarizability,sig_pol;
} t_ffatype;

typedef struct {
    char   *atom1,*atom2,*params;
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
    int nspoel,nalexandria_c;
    t_ffatype *spoel;
    char *alexandria_polar_unit;
    char *alexandria_forcefield;
    int ngt_bond,ngt_bond_c;
    char *gt_length_unit;
    t_gt_bond *gt_bond;
    int ngt_angle,ngt_angle_c;
    char *gt_angle_unit;
    t_gt_angle *gt_angle;
    int ngt_dihedral,ngt_dihedral_c;
    t_gt_dihedral *gt_dihedral;
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
    gmx_poldata *pold;
  
    snew(pold,1);
  
    return (gmx_poldata_t) pold;
}

void gmx_poldata_add_ffatype(gmx_poldata_t pd,char *elem,char *desc,
                             char *gt_name,char *gt_type,
                             char *miller_equiv,char *charge,
                             char *geometry,int numbonds,
                             char *neighbors,
                             double polarizability,double sig_pol,
                             char *vdwparams)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    char    buf[EEMBUFSIZE],**ptr;
    int     i,j;
  
    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(pold->spoel[i].name,gt_name) == 0)
            break;
    if (i == pold->nspoel) {
        pold->nspoel++;
        srenew(pold->spoel,pold->nspoel);
  
        sp = &(pold->spoel[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->name           = strdup(gt_name);
        sp->type           = strdup(gt_type);
        sp->miller_equiv   = strdup(miller_equiv);
        sp->neighbors      = strdup(neighbors);
        sp->nb             = split(' ',neighbors);
        for(sp->numnb=0; (sp->nb[sp->numnb] != NULL); sp->numnb++)
            ;
        sp->charge         = strdup(charge);
        sp->geometry       = strdup(geometry);
        sp->numbonds       = numbonds;
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
        sp->vdwparams      = strdup(vdwparams);
    }
    else 
        fprintf(stderr,"Atom %s was already added to poldata record\n",gt_name);
}

int gmx_poldata_get_natypes(gmx_poldata_t pd)
{
    return pd->nspoel;
}

void gmx_poldata_set_ffatype(gmx_poldata_t pd,char *gt_type,
                           double polarizability,double sig_pol)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    int i;
  
    for(i=0; (i<pold->nspoel); i++) {
        if (strcasecmp(gt_type,pold->spoel[i].type) == 0) {
            sp = &(pold->spoel[i]);
            sp->polarizability = polarizability;
            sp->sig_pol = sig_pol;
        }
    }
}		

char *gmx_poldata_get_ffatype_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->alexandria_polar_unit;
}

char *gmx_poldata_get_length_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->gt_length_unit;
}

void gmx_poldata_set_ffatype_unit(gmx_poldata_t pd,char *polar_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->alexandria_polar_unit   = strdup(polar_unit);
}

char *gmx_poldata_get_ffatype_ff(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->alexandria_forcefield;
}

void gmx_poldata_set_ffatype_ff(gmx_poldata_t pd,char *forcefield)
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
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (gt_atom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,gt_atom) == 0) 
                return pold->spoel[i].geometry;
      
    return NULL;
}

char *gmx_poldata_get_desc(gmx_poldata_t pd,char *gt_atom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (gt_atom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,gt_atom) == 0) 
                return pold->spoel[i].desc;
      
    return NULL;
}

char *gmx_poldata_get_gt_type(gmx_poldata_t pd,char *gt_atom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (gt_atom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,gt_atom) == 0) 
                return pold->spoel[i].type;
  
    return NULL;
}

static gmx_bool strcasestr_start(char *needle,char *haystack)
{
    char *ptr;
    
    ptr = strcasestr(haystack,needle);
    
    return (ptr == haystack);
}

static int count_neighbors(t_ffatype *spoel,int nbond,char *nbhybrid[])
{
    int i,j,ni=0,*jj,i_found;
        
    snew(jj,nbond+1);
    for(i=0; (i<spoel->numnb); i++)
    {
        i_found = 0;
        for(j=0; (j<nbond); j++) 
        {
            if (
                (NULL != nbhybrid[j]) &&
                (jj[j] == 0) &&
                (i_found == 0) &&
                strcasestr_start(spoel->nb[i],nbhybrid[j])
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

char *gmx_poldata_get_charge(gmx_poldata_t pd,char *gt_atom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (gt_atom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,gt_atom) == 0) 
                return pold->spoel[i].charge;
      
    return NULL;
}

char *gmx_poldata_get_gt_atom(gmx_poldata_t pd,char *elem,
                             int nbond,char *neighbors[],
                             const char *geometry)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i,nnb,pbest=-1,best=-1;
  
    for(i=0; (i<pold->nspoel); i++) 
    {
        nnb = count_neighbors(&(pold->spoel[i]),nbond,neighbors);
        if ((strcasecmp(pold->spoel[i].elem,elem) == 0) &&
            (strcasecmp(pold->spoel[i].geometry,geometry) == 0) &&
            (nbond == pold->spoel[i].numbonds) &&
            (nnb == pold->spoel[i].numnb))
        {
            if ((NULL != debug) && (best != -1) && (nnb == pbest))
            {
                fprintf(debug,"poldata clash: %s and %s are equally well suited, nnb = %d. Using %s.\n",
                        pold->spoel[best].name,pold->spoel[i].name,nnb,
                        pold->spoel[best].name);
            }
            else if (nnb > pbest)
            {
                pbest = nnb;
                best = i;
            }
        }
    }
    
    if (best != -1)
        return pold->spoel[best].name;
    else
        return NULL;
}

int gmx_poldata_gt_type_polarizability(gmx_poldata_t pd,char *gt_type,double *polar,
                                       double *sig_pol)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;

    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(gt_type,pold->spoel[i].type) == 0)
        {
            *polar   = pold->spoel[i].polarizability;
            *sig_pol = pold->spoel[i].sig_pol;
            return 1;
        }
    return 0;
}

char *gmx_poldata_get_miller_equiv(gmx_poldata_t pd,char *gt_name)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    int i;

    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(gt_name,pold->spoel[i].name) == 0) 
            return pold->spoel[i].miller_equiv;
    return NULL;
}

char *gmx_poldata_get_ffatype(gmx_poldata_t pd,char *gt_name,char **elem,char **desc,
                              char **gt_type,char **miller_equiv,
                              char **charge,char **geometry,
                              int *numbonds,char **neighbors,
                              double *polarizability,double *sig_pol,
                              char **vdwparams)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_ffatype *sp;
    int i;
  
    if (gt_name) {
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(gt_name,pold->spoel[i].type) == 0) 
                break;
    }
    else
        i = pold->nalexandria_c;
      
    if (i<pold->nspoel) {
        sp = &(pold->spoel[i]);
        assign_str(charge,sp->charge);
        assign_str(geometry,sp->geometry);
        assign_scal(polarizability,sp->polarizability);
        assign_scal(sig_pol,sp->sig_pol);
        assign_scal(numbonds,sp->numbonds);
        assign_str(elem,sp->elem);
        assign_str(desc,sp->desc);
        assign_str(gt_type,sp->type);
        assign_str(miller_equiv,sp->miller_equiv);
        assign_str(vdwparams,sp->vdwparams);
        assign_str(neighbors,sp->neighbors);
        if (!gt_name)
            pold->nalexandria_c++;
    
        return sp->name;
    }
    else
        pold->nalexandria_c = 0;
    
    return NULL;
}
				  
double gmx_poldata_get_bondorder(gmx_poldata_t pd,char *elem1,char *elem2,
                                 double distance,double toler)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    double dev,dev_best = 100;
    char *ba1,*ba2;
    int j,i,i_best=-1;
  
    if ((NULL == elem1) || (NULL == elem2))
        return 0.0;
    for(i=0; (i<pold->ngt_bond); i++) {
        ba1 = pold->gt_bond[i].atom1;
        ba2 = pold->gt_bond[i].atom2;
        for(j=0; (j<pold->nspoel); j++) {
            if (strcasecmp(pold->spoel[j].name,pold->gt_bond[i].atom1) == 0)
                ba1 = pold->spoel[j].elem;
            if (strcasecmp(pold->spoel[j].name,pold->gt_bond[i].atom2) == 0)
                ba2 = pold->spoel[j].elem;
        }
        if (((NULL != ba1) && (NULL != ba2)) &&
            (((strcasecmp(ba1,elem1) == 0) && (strcasecmp(ba2,elem2) == 0)) ||
             ((strcasecmp(ba1,elem2) == 0) && (strcasecmp(ba2,elem1) == 0)))) {
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

/* 
 * gt_bond stuff 
 */
void gmx_poldata_add_gt_bond(gmx_poldata_t pd,char *atom1,char *atom2,
                             double length,double sigma,double bondorder,char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_bond *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_bond); i++) {
        gt_b = &(pold->gt_bond[i]);
        if ((((strcasecmp(gt_b->atom1,atom1) == 0) && 
              (strcasecmp(gt_b->atom2,atom2) == 0)) ||
             ((strcasecmp(gt_b->atom1,atom2) == 0) && 
              (strcasecmp(gt_b->atom2,atom1) == 0))) &&
            (gt_b->bondorder == bondorder))
            break;
    }
    if (i < pold->ngt_bond) {
        if (gt_b->length != length) {
            fprintf(stderr,"Updating bond %s-%s, bondorder %g, length %g with different bondlength %g.\n",
                    atom1,atom2,bondorder,gt_b->length,length);
            gt_b->length = length;
            gt_b->sigma = sigma;
            gt_b->params = strdup(params);
        }
    }
    else {
        pold->ngt_bond++;
        srenew(pold->gt_bond,pold->ngt_bond);
        gt_b = &(pold->gt_bond[pold->ngt_bond-1]);
        gt_b->atom1   = strdup(atom1);
        gt_b->atom2   = strdup(atom2);
        gt_b->length = length;
        gt_b->sigma = sigma;
        gt_b->bondorder = bondorder;
        gt_b->params = strdup(params);
    }
}

int gmx_poldata_get_gt_bond(gmx_poldata_t pd,char **atom1,char **atom2,
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
    
        return 1;
    }
    pold->ngt_bond_c = 0;
  
    return 0;
}

int gmx_poldata_search_gt_bond(gmx_poldata_t pd,char *atom1,char *atom2,
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
        
            return 1;
        }
    }
    return 0;
}

/*
 * gt_angle stuff
 */
void gmx_poldata_add_gt_angle(gmx_poldata_t pd,char *atom1,char *atom2,
                              char *atom3,double angle,double sigma,char *params)
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
        if (gt_b->angle != angle) {
            fprintf(stderr,"Updating angle %s-%s-%s with different angle %g (was %g).\n",
                    atom1,atom2,atom3,angle,gt_b->angle);
            gt_b->angle = angle;
            gt_b->sigma = sigma;
            gt_b->params = strdup(params);
        }
    }
    else {
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
}
				   
int gmx_poldata_get_gt_angle(gmx_poldata_t pd,char **atom1,char **atom2,
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
    
        return 1;
    }
    pold->ngt_angle_c = 0;
  
    return 0;
}

int gmx_poldata_search_gt_angle(gmx_poldata_t pd,char *atom1,char *atom2,
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
        
            return 1;
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
void gmx_poldata_add_gt_dihedral(gmx_poldata_t pd,char *atom1,char *atom2,
                                 char *atom3,char *atom4,double dihedral,
                                 double sigma,char *params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_dihedral *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_dihedral); i++) {
        gt_b = &(pold->gt_dihedral[i]);
        if (((strcasecmp(gt_b->atom1,atom1) == 0) && 
             (strcasecmp(gt_b->atom2,atom2) == 0) &&
             (strcasecmp(gt_b->atom3,atom3) == 0) && 
             (strcasecmp(gt_b->atom4,atom4) == 0)) ||
            ((strcasecmp(gt_b->atom1,atom4) == 0) && 
             (strcasecmp(gt_b->atom2,atom3) == 0) &&
             (strcasecmp(gt_b->atom3,atom2) == 0) && 
             (strcasecmp(gt_b->atom4,atom1) == 0))) 
            break;
    }
    if (i < pold->ngt_dihedral) {
        if (gt_b->dihedral != dihedral) {
            fprintf(stderr,"Updating dihedral %s-%s-%s-%s with different dihedral %g (was %g).\n",
                    atom1,atom2,atom3,atom4,dihedral,gt_b->dihedral);
            gt_b->dihedral = dihedral;
            gt_b->sigma = sigma;
            gt_b->params = strdup(params);
        }
    }
    else {
        pold->ngt_dihedral++;
        srenew(pold->gt_dihedral,pold->ngt_dihedral);
        gt_b = &(pold->gt_dihedral[pold->ngt_dihedral-1]);
        gt_b->atom1   = strdup(atom1);
        gt_b->atom2   = strdup(atom2);
        gt_b->atom3   = strdup(atom3);
        gt_b->atom4   = strdup(atom4);
        gt_b->dihedral   = dihedral;
        gt_b->sigma   = sigma;
        gt_b->params = strdup(params);
    }
}
				   
int gmx_poldata_get_gt_dihedral(gmx_poldata_t pd,char **atom1,char **atom2,
                                char **atom3,char **atom4,double *dihedral,
                                double *sigma,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_dihedral *gt_b;
    int i;
  
    if (pold->ngt_dihedral_c < pold->ngt_dihedral) {
        gt_b = &(pold->gt_dihedral[pold->ngt_dihedral_c]);
        assign_str(atom1,gt_b->atom1);
        assign_str(atom2,gt_b->atom2);
        assign_str(atom3,gt_b->atom3);
        assign_str(atom4,gt_b->atom4);
        assign_scal(dihedral,gt_b->dihedral);
        assign_scal(sigma,gt_b->sigma);
        assign_str(params,gt_b->params);
        pold->ngt_dihedral_c++;
    
        return 1;
    }
    pold->ngt_dihedral_c = 0;
  
    return 0;
}

int gmx_poldata_search_gt_dihedral(gmx_poldata_t pd,char *atom1,char *atom2,
                                   char *atom3,char *atom4,
                                   double *dihedral,double *sigma,char **params)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_gt_dihedral *gt_b;
    int i;
  
    for(i=0; (i<pold->ngt_dihedral); i++) {
        gt_b = &(pold->gt_dihedral[i]);
        if (((strcasecmp(gt_b->atom1,atom1) == 0) && 
             (strcasecmp(gt_b->atom2,atom2) == 0) &&
             (strcasecmp(gt_b->atom3,atom3) == 0) && 
             (strcasecmp(gt_b->atom4,atom4) == 0)) ||
            ((strcasecmp(gt_b->atom1,atom4) == 0) && 
             (strcasecmp(gt_b->atom2,atom3) == 0) &&
             (strcasecmp(gt_b->atom3,atom2) == 0) && 
             (strcasecmp(gt_b->atom4,atom1) == 0))) {
            assign_scal(dihedral,gt_b->dihedral);
            assign_scal(sigma,gt_b->sigma);
            assign_str(params,gt_b->params);
        
            return 1;
        }
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
    
    for(i=0; (i<pd->nspoel); i++)
        if (strcasecmp(gt_type,pd->spoel[i].type) == 0)
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

