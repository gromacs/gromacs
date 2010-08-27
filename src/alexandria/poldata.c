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
    char   *desc,*name,*type,*elem,*miller_equiv,*spref,*neighbors,*geometry;
    char   **nb;
    int    charge,numbonds,numnb;
    double polarizability,sig_pol;
} t_alexandria;

typedef struct {
    char   *atom1,*atom2;
    double length,bondorder; 
} t_smbond;

typedef struct {
    char   *atom1,*atom2,*atom3;
    double angle; 
} t_smangle;

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
    double J0,q[MAXZETA],zeta[MAXZETA],chi0; 
} t_eemprops;

typedef struct {
    int  eqg_model;
    char *epref;
} t_epref;

typedef struct gmx_poldata {
    int nspoel,nalexandria_c;
    t_alexandria *spoel;
    char *alexandria_polar_unit;
    int nsmbond,nsmbond_c;
    char *sm_length_unit;
    t_smbond *smbond;
    int nsmangle,nsmangle_c;
    char *sm_angle_unit;
    t_smangle *smangle;
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

void gmx_poldata_add_alexandria(gmx_poldata_t pd,char *elem,char *desc,
                           char *smname,char *smtype,
                           char *miller_equiv,
                           int charge,
                           char *geometry,int numbonds,
                           char *neighbors,
                           double polarizability,double sig_pol,
                           char *spref)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_alexandria *sp;
    char    buf[EEMBUFSIZE];
    int     i;
  
    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(pold->spoel[i].name,smname) == 0)
            break;
    if (i == pold->nspoel) {
        pold->nspoel++;
        srenew(pold->spoel,pold->nspoel);
  
        sp = &(pold->spoel[i]);
        sp->elem           = strdup(elem);
        sp->desc           = strdup(desc);
        sp->name           = strdup(smname);
        sp->type           = strdup(smtype);
        sp->miller_equiv   = strdup(miller_equiv);
        sp->neighbors      = strdup(neighbors);
        sp->nb             = split(' ',neighbors);
        for(sp->numnb=0; (sp->nb[sp->numnb] != NULL); sp->numnb++)
            ;
        sp->charge         = charge;
        sp->geometry       = strdup(geometry);
        sp->numbonds       = numbonds;
        sp->polarizability = polarizability;
        sp->sig_pol        = sig_pol;
        sp->spref          = strdup(spref);
    }
    else 
        fprintf(stderr,"Atom %s was already added to poldata record\n",smname);
}

void gmx_poldata_set_alexandria(gmx_poldata_t pd,char *smtype,
                           double polarizability,double sig_pol,char *spref)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_alexandria *sp;
    int i;
  
    for(i=0; (i<pold->nspoel); i++) {
        if (strcasecmp(smtype,pold->spoel[i].type) == 0) {
            sp = &(pold->spoel[i]);
            sp->polarizability = polarizability;
            sp->sig_pol = sig_pol;
            if (sp->spref)
                sfree(sp->spref);
            sp->spref = strdup(spref);
        }
    }
}		

char *gmx_poldata_get_alexandria_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->alexandria_polar_unit;
}

char *gmx_poldata_get_length_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->sm_length_unit;
}

void gmx_poldata_set_alexandria_unit(gmx_poldata_t pd,char *polar_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->alexandria_polar_unit   = strdup(polar_unit);
}

void gmx_poldata_set_length_unit(gmx_poldata_t pd,char *length_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->sm_length_unit   = strdup(length_unit);
}

double gmx_poldata_get_alpha(gmx_poldata_t pd,char *smatom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (smatom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,smatom) == 0) 
                return pold->spoel[i].polarizability;
      
    return -1;
}

char *gmx_poldata_get_geometry(gmx_poldata_t pd,char *smatom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (smatom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,smatom) == 0) 
                return pold->spoel[i].geometry;
      
    return NULL;
}

char *gmx_poldata_get_desc(gmx_poldata_t pd,char *smatom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (smatom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,smatom) == 0) 
                return pold->spoel[i].desc;
      
    return NULL;
}

char *gmx_poldata_get_smtype(gmx_poldata_t pd,char *smatom)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;
  
    if (smatom)
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(pold->spoel[i].name,smatom) == 0) 
                return pold->spoel[i].type;
  
    return NULL;
}

static int count_neighbors(t_alexandria *spoel,int nbond,char *nbhybrid[])
{
    int i,j,ni=0,*jj,i_found;
        
    snew(jj,nbond+1);
    for(i=0; (i<spoel->numnb); i++)
    {
        i_found = 0;
        for(j=0; (j<nbond); j++) 
        {
            if ((NULL != nbhybrid[j]) &&
                (jj[j] == 0) &&
                (i_found == 0) &&
                (NULL != strcasestr(nbhybrid[j],spoel->nb[i])))
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

char *gmx_poldata_get_smatom(gmx_poldata_t pd,char *elem,
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

int gmx_poldata_smtype_polarizability(gmx_poldata_t pd,char *smtype,double *polar,
                                      double *sig_pol)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    int i;

    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(smtype,pold->spoel[i].type) == 0)
        {
            *polar   = pold->spoel[i].polarizability;
            *sig_pol = pold->spoel[i].sig_pol;
            return 1;
        }
    return 0;
}

char *gmx_poldata_get_miller_equiv(gmx_poldata_t pd,char *smname)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_alexandria *sp;
    int i;

    for(i=0; (i<pold->nspoel); i++) 
        if (strcasecmp(smname,pold->spoel[i].name) == 0) 
            return pold->spoel[i].miller_equiv;
    return NULL;
}

char *gmx_poldata_get_alexandria(gmx_poldata_t pd,char *smname,char **elem,char **desc,
                            char **smtype,char **miller_equiv,
                            int *charge,char **geometry,
                            int *numbonds,char **neighbors,
                            double *polarizability,double *sig_pol,char **spref)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_alexandria *sp;
    int i;
  
    if (smname) {
        for(i=0; (i<pold->nspoel); i++) 
            if (strcasecmp(smname,pold->spoel[i].type) == 0) 
                break;
    }
    else
        i = pold->nalexandria_c;
      
    if (i<pold->nspoel) {
        sp = &(pold->spoel[i]);
        assign_scal(charge,sp->charge);
        assign_str(geometry,sp->geometry);
        assign_scal(polarizability,sp->polarizability);
        assign_scal(sig_pol,sp->sig_pol);
        assign_scal(numbonds,sp->numbonds);
        assign_str(elem,sp->elem);
        assign_str(desc,sp->desc);
        assign_str(smtype,sp->type);
        assign_str(miller_equiv,sp->miller_equiv);
        assign_str(spref,sp->spref);
        assign_str(neighbors,sp->neighbors);
        if (!smname)
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
    int i,i_best=-1;
  
    if ((NULL == elem1) || (NULL == elem2))
        return 0.0;
    for(i=0; (i<pold->nsmbond); i++) {
        if (((strcasecmp(pold->smbond[i].atom1,elem1) == 0) &&
             (strcasecmp(pold->smbond[i].atom2,elem2) == 0)) ||
            ((strcasecmp(pold->smbond[i].atom1,elem2) == 0) &&
             (strcasecmp(pold->smbond[i].atom2,elem1) == 0))) {
            dev = fabs((pold->smbond[i].length - distance)/pold->smbond[i].length);
            if (dev < dev_best) {
                dev_best = dev;
                i_best = i;
            }
        }
    }
    if (dev_best < toler)
        return pold->smbond[i_best].bondorder;
  
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
				  
void gmx_poldata_add_smbond(gmx_poldata_t pd,char *atom1,char *atom2,
                            double length,double bondorder)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smbond *smb;
    int i;
  
    for(i=0; (i<pold->nsmbond); i++) {
        smb = &(pold->smbond[i]);
        if ((((strcasecmp(smb->atom1,atom1) == 0) && 
              (strcasecmp(smb->atom2,atom2) == 0)) ||
             ((strcasecmp(smb->atom1,atom2) == 0) && 
              (strcasecmp(smb->atom2,atom1) == 0))) &&
            (smb->bondorder == bondorder))
            break;
    }
    if (i < pold->nsmbond) {
        if (smb->length != length)
            gmx_fatal(FARGS,"Trying to add bond %s-%s, bondorder %g, with different bondlength %g",
                      atom1,atom2,bondorder,length);
    }
    else {
        pold->nsmbond++;
        srenew(pold->smbond,pold->nsmbond);
        smb = &(pold->smbond[pold->nsmbond-1]);
        smb->atom1   = strdup(atom1);
        smb->atom2   = strdup(atom2);
        smb->length = length;
        smb->bondorder = bondorder;
    }
}

int gmx_poldata_get_smbond(gmx_poldata_t pd,char **atom1,char **atom2,
                           double *length,double *bondorder)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smbond *smb;
    int i;
  
    if (pold->nsmbond_c < pold->nsmbond) {
        smb = &(pold->smbond[pold->nsmbond_c]);
        assign_str(atom1,smb->atom1);
        assign_str(atom2,smb->atom2);
        assign_scal(length,smb->length);
        assign_scal(bondorder,smb->bondorder);
        pold->nsmbond_c++;
    
        return 1;
    }
    pold->nsmbond_c = 0;
  
    return 0;
}

int gmx_poldata_search_smbond(gmx_poldata_t pd,char *atom1,char *atom2,
                              double *length)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smbond *smb;
    int i;
  
    for(i=0; (i<pold->nsmbond); i++) {
        smb = &(pold->smbond[i]);
        if (((strcasecmp(smb->atom1,atom1) == 0) && 
             (strcasecmp(smb->atom2,atom2) == 0)) ||
            ((strcasecmp(smb->atom1,atom2) == 0) && 
             (strcasecmp(smb->atom2,atom1) == 0))) {
            assign_scal(length,smb->length);

            return 1;
        }
    }
    return 0;
}

void gmx_poldata_add_smangle(gmx_poldata_t pd,char *atom1,char *atom2,
                             char *atom3,double angle)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smangle *smb;
    int i;
  
    for(i=0; (i<pold->nsmangle); i++) {
        smb = &(pold->smangle[i]);
        if ((strcasecmp(smb->atom2,atom2) == 0) && 
            (((strcasecmp(smb->atom1,atom1) == 0) && 
              (strcasecmp(smb->atom3,atom3) == 0)) ||
             ((strcasecmp(smb->atom1,atom3) == 0) && 
              (strcasecmp(smb->atom3,atom1) == 0))))
            break;
    }
    if (i < pold->nsmangle) {
        if (smb->angle != angle)
            gmx_fatal(FARGS,"Trying to add angle %s-%s-%s with different angle %g",
                      atom1,atom2,atom3,angle);
    }
    else {
        pold->nsmangle++;
        srenew(pold->smangle,pold->nsmangle);
        smb = &(pold->smangle[pold->nsmangle-1]);
        smb->atom1   = strdup(atom1);
        smb->atom2   = strdup(atom2);
        smb->atom3   = strdup(atom3);
        smb->angle   = angle;
    }
}
				   
int gmx_poldata_get_smangle(gmx_poldata_t pd,char **atom1,char **atom2,
                            char **atom3,double *angle)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smangle *smb;
    int i;
  
    if (pold->nsmangle_c < pold->nsmangle) {
        smb = &(pold->smangle[pold->nsmangle_c]);
        assign_str(atom1,smb->atom1);
        assign_str(atom2,smb->atom2);
        assign_str(atom3,smb->atom3);
        assign_scal(angle,smb->angle);
        pold->nsmangle_c++;
    
        return 1;
    }
    pold->nsmangle_c = 0;
  
    return 0;
}

int gmx_poldata_search_smangle(gmx_poldata_t pd,char *atom1,char *atom2,
                               char *atom3,double *angle)
{
    gmx_poldata *pold = (gmx_poldata *) pd;
    t_smangle *smb;
    int i;
  
    for(i=0; (i<pold->nsmangle); i++) {
        smb = &(pold->smangle[i]);
        if ((strcasecmp(smb->atom2,atom2) == 0) &&
            (((strcasecmp(smb->atom1,atom1) == 0) && 
              (strcasecmp(smb->atom3,atom3) == 0)) ||
             ((strcasecmp(smb->atom1,atom3) == 0) && 
              (strcasecmp(smb->atom3,atom1) == 0)))) {
            assign_scal(angle,smb->angle);
      
            return 1;
        }
    }
    return 0;
}

void gmx_poldata_set_angle_unit(gmx_poldata_t pd,char *angle_unit)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    pold->sm_angle_unit   = strdup(angle_unit);
}

char *gmx_poldata_get_angle_unit(gmx_poldata_t pd)
{
    gmx_poldata *pold = (gmx_poldata *) pd;

    return pold->sm_angle_unit;
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

int gmx_poldata_have_pol_support(gmx_poldata_t pd,char *smtype)
{
    int i;
    
    for(i=0; (i<pd->nspoel); i++)
        if (strcasecmp(smtype,pd->spoel[i].type) == 0)
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
    int i,nep;
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
        for(i=0; (i<pd->nep); i++)
            fprintf(debug,"EEP: %5d Atom: %5s  Chi: %8.3f  J00:  %8.3f  Z1: %8.3f  Z2: %8.3f\n",
                    pd->eep[i].eqg_model,pd->eep[i].name,
                    pd->eep[i].chi0,pd->eep[i].J0,
                    pd->eep[i].zeta[0],pd->eep[i].zeta[1]);
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

