/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop.c,v 1.32 2009/06/01 06:13:18 spoel Exp $
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "molprop.h"
	
const char *emp_name[empNR] = 
{ 
    "potential", "dipole", "quadrupole", "polarizability", "energy" 
};

typedef struct 
{
    char *cname;
    int  cnumber;
} t_catom;

typedef struct 
{
    char    *compname;
    int     ncatom,ncatom_c;
    t_catom *catom;
} t_composition;

typedef struct 
{
    int    ref;
    char   *type,*unit;
    double xx,yy,zz,xy,xz,yz;
} t_quadrupole;

typedef struct 
{
    int    ref;
    char   *type,*unit;
    double value,error;
} t_mpenergy;

typedef struct 
{
    int    ref;
    char   *type,*unit;
    double x,y,z;
    double aver,error;
} t_dipole;

typedef t_dipole t_polarizability;

typedef struct { 
    char *xyz_unit,*V_unit;
    int espid;
    double x,y,z,V;
} t_potential;

typedef struct 
{
    char   *type,*unit;
    double q;
} t_charge;

typedef struct 
{
    char     *name,*unit;
    double   x,y,z;
    int      atomid,nq,nq_c;
    t_charge *q;
} t_calc_atom;

typedef struct 
{
    int              eMP;
    char             *reference,*conformation;
    int              nprop[empNR],nprop_c[empNR];
    t_polarizability *polar;
    t_dipole         *dipole;
    t_mpenergy       *energy;
    t_quadrupole     *quadrupole;
    t_potential      *potential;
    /* The items below are used only in a calculation */
    char             *program,*method,*basisset;
    int              natom,natom_c;
    t_calc_atom      *catom;
} t_calc_exp;

typedef struct gmx_molprop
{
    char          *formula,*molname,*iupac,*cas,*cid,*inchi;
    double        mass;
    int           charge,multiplicity;
    int           nexper,nexper_c;
    t_calc_exp    *exper;
    int           ncalc,ncalc_c;
    t_calc_exp    *calc;
    int           ncomposition,ncomp_c;
    t_composition *composition;
    int           ncategory,ncateg_c;
    char          **category;
} gmx_molprop;

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

gmx_molprop_t gmx_molprop_init()
{
    gmx_molprop *mp;
  
    snew(mp,1);
  
    return (gmx_molprop_t) mp;
}

void gmx_molprop_check_consistency(gmx_molprop_t mpt)
{
    int i,j;
    double dx,dy,dz,dd;
        
    for(i=0; (i<mpt->ncalc); i++) 
    {
        if (mpt->calc[i].nprop[empPOTENTIAL] >= mpt->calc[i].natom)
        {
            for(j=0; (j<mpt->calc[i].natom); j++)
            {
                dx = mpt->calc[i].catom[j].x - mpt->calc[i].potential[j].x;
                dy = mpt->calc[i].catom[j].y - mpt->calc[i].potential[j].y;
                dz = mpt->calc[i].catom[j].z - mpt->calc[i].potential[j].z;
                dd = sqrt(dx*dx+dy*dy+dz*dz);
                if (dd > 1e-3)
                {
                    fprintf(stderr,"Mol. %s atom %d different coordinate ( %g %g %g )\n",
                            mpt->molname,j+1,dx,dy,dz);
                }
            }
        }
    }
}

void gmx_molprop_delete_composition(gmx_molprop_t mp,const char *compname)
{
    int i,j;
  
    for(i=0; (i<mp->ncomposition); i++)
        if (strcasecmp(compname,mp->composition[i].compname) == 0)
            break;
    if (i < mp->ncomposition) 
    {
        for(j=0; (j<mp->composition[i].ncatom); j++) 
        {
            sfree(mp->composition[i].catom[j].cname);
            mp->composition[i].catom[j].cname = NULL;
            mp->composition[i].catom[j].cnumber = 0;
        }
        sfree(mp->composition[i].catom);
        mp->composition[i].catom = NULL;
        mp->composition[i].ncatom = 0;
    }  
}

static void t_calc_exp_delete(t_calc_exp *ce)
{
    int j;
    
    sfree(ce->reference);
    sfree(ce->method);
    sfree(ce->program);
    sfree(ce->basisset);
    for(j=0; (j<ce->natom); j++) 
    {
        sfree(ce->catom[j].name);
        sfree(ce->catom[j].unit);
    }
}

void gmx_molprop_delete(gmx_molprop_t mp)
{
    int i,j;
  
    for(i=0; (i<mp->ncomposition); i++)
        gmx_molprop_delete_composition(mp,mp->composition[i].compname);
    sfree(mp->composition);
    mp->ncomposition = mp->ncomp_c = 0;
    for(i=0; (i<mp->ncategory); i++)
        sfree(mp->category[i]);
    sfree(mp->category);
    mp->ncategory = mp->ncateg_c = 0;
    mp->nexper = mp->nexper_c = 0;
    mp->ncalc = mp->ncalc_c = 0;
    sfree(mp->molname);
    sfree(mp->formula);
    sfree(mp->iupac);
    sfree(mp->cas);
    sfree(mp->cid);
    sfree(mp->inchi);
    for(i=0; (i<mp->nexper); i++) 
        t_calc_exp_delete(&mp->exper[i]);
    for(i=0; (i<mp->ncalc); i++) 
        t_calc_exp_delete(&mp->calc[i]);
  
    if (0)
        printf("gmx_molprop_delete not completely implemented yet\n");
}

void gmx_molprop_set_mass(gmx_molprop_t mp,double mass)
{
    if (mass > 1) 
    {
        if ((NULL != debug) && (mp->mass > 1) && (fabs(mp->mass-mass) > 1e-5))
            fprintf(debug,"Changing mass for %s from %g to %g\n",
                    mp->molname,mp->mass,mass);
        mp->mass=mass;
    }
}

double gmx_molprop_get_mass(gmx_molprop_t mp)
{
    return mp->mass;
}

void gmx_molprop_set_charge(gmx_molprop_t mp,int charge)
{
    if ((NULL != debug) && (mp->charge != 0) && (mp->charge != charge))
        fprintf(debug,"Changing charge for %s from %d to %d\n",
                mp->molname,mp->charge,charge);
    mp->charge = charge;
}

int gmx_molprop_get_charge(gmx_molprop_t mp)
{
    return mp->charge;
}

void gmx_molprop_set_multiplicity(gmx_molprop_t mp,int multiplicity)
{
    if (multiplicity > 0) 
    {
        if ((NULL != debug) && (mp->multiplicity != 0) && (mp->multiplicity != multiplicity))
            fprintf(debug,"Changing multiplicity for %s from %d to %d\n",
                    mp->molname,mp->multiplicity,multiplicity);
        mp->multiplicity=multiplicity;
    }
}

int gmx_molprop_get_multiplicity(gmx_molprop_t mp)
{
    return mp->multiplicity;
}

void gmx_molprop_set_formula(gmx_molprop_t mp,char *formula)
{
    if (formula) 
    {
        if (strlen(formula) > 0)
        {
            if (mp->formula != NULL) 
            {
                if ((NULL != debug) && (strcasecmp(formula,mp->formula) != 0))
                    fprintf(debug,"Replacing molecule %s formula '%s' by '%s'\n",
                            mp->molname,mp->formula,formula);
                sfree(mp->formula);
            }
            mp->formula = strdup(formula);
        }
        else if (NULL == mp->formula)
            mp->formula = strdup(formula);
    }
}

char *gmx_molprop_get_formula(gmx_molprop_t mp)
{
    return mp->formula;
}

void gmx_molprop_set_molname(gmx_molprop_t mp,char *molname)
{
    if (molname) 
    {
        if (mp->molname != NULL) {
            if (NULL != debug)
                fprintf(debug,"Replacing molname %s by %s\n",mp->molname,molname);
            sfree(mp->molname);
        }
        mp->molname = strdup(molname);
    }
    else 
        gmx_incons("Trying to set molname without actually supplying a valid argument");
}

char *gmx_molprop_get_molname(gmx_molprop_t mp)
{
    return mp->molname;
}

void gmx_molprop_set_iupac(gmx_molprop_t mp,char *iupac)
{
    assign_str(&mp->iupac,iupac);
}			   

char *gmx_molprop_get_iupac(gmx_molprop_t mp)
{
    return mp->iupac;
}

void gmx_molprop_set_inchi(gmx_molprop_t mp,char *inchi)
{
    assign_str(&mp->inchi,inchi);
}			   

char *gmx_molprop_get_inchi(gmx_molprop_t mp)
{
    return mp->inchi;
}

void gmx_molprop_set_cas(gmx_molprop_t mp,char *cas)
{
    assign_str(&mp->cas,cas);
}			   

char *gmx_molprop_get_cas(gmx_molprop_t mp)
{
    return mp->cas;
}

void gmx_molprop_set_cid(gmx_molprop_t mp,char *cid)
{
    assign_str(&mp->cid,cid);
}			   

char *gmx_molprop_get_cid(gmx_molprop_t mp)
{
    return mp->cid;
}

void gmx_molprop_add_composition(gmx_molprop_t mp,const char *compname)
{
    int i,j;
  
    for(i=0; (i<mp->ncomposition); i++) 
    {
        if(strcasecmp(compname,mp->composition[i].compname) == 0)
            return;
    }
    mp->ncomposition++;
    srenew(mp->composition,mp->ncomposition);
    mp->composition[i].compname = strdup(compname);
    mp->composition[i].ncatom = 0;
    mp->composition[i].ncatom_c = 0;
    mp->composition[i].catom = NULL;
}

char *gmx_molprop_get_composition(gmx_molprop_t mp)
{
    if (mp->ncomp_c < mp->ncomposition) 
    {
        return mp->composition[mp->ncomp_c++].compname;
    }
    mp->ncomp_c = 0;
    return NULL;
}

int gmx_molprop_get_natom(gmx_molprop_t mpt)
{
    char *comp,*atom;
    int  n=0,natom;
    
    while (0 !=  gmx_molprop_get_composition_atom(mpt,"bosque",
                                                  &atom,&natom))
    {
        n += natom;
        sfree(atom);
    }
    return n;
}

int gmx_molprop_count_composition_atoms(gmx_molprop_t mp,
                                        char *compname,char *atom)
{
    int i,j,n;
  
    for(i=0; (i<mp->ncomposition); i++)
        if (strcasecmp(compname,mp->composition[i].compname) == 0)
            break;
    if (i < mp->ncomposition) 
    {
        if (NULL == atom) 
        {
            n = 0;
            for(j=0; (j<mp->composition[i].ncatom); j++)
                n = n + mp->composition[i].catom[j].cnumber;
            return n;
        }
        else 
        {
            for(j=0; (j<mp->composition[i].ncatom); j++)
                if (strcasecmp(mp->composition[i].catom[j].cname,atom) == 0)
                    return mp->composition[i].catom[j].cnumber;
        }  
    }
    return 0;
}
	
static int catom_comp(const void *a,const void *b)
{
    t_catom *ca = (t_catom *)a;
    t_catom *cb = (t_catom *)b;
    
    return strcasecmp(ca->cname,cb->cname);
}
				       
void gmx_molprop_add_composition_atom(gmx_molprop_t mp,const char *compname,
                                      char *atomname,int natom)
{
    int i=-1,j;
  
    if (compname) 
    {
        for(i=0; (i<mp->ncomposition); i++)
            if(strcasecmp(compname,mp->composition[i].compname) == 0)
                break;
        if (i == mp->ncomposition) 
            gmx_molprop_add_composition(mp,compname);
    }
    else if (mp->ncomposition > 0) 
        i = mp->ncomposition-1;
    
    if (i >= 0) 
    {
        for(j=0; (j<mp->composition[i].ncatom); j++) 
        {
            if (strcasecmp(mp->composition[i].catom[j].cname,atomname) == 0) 
            {
                break;
            }
        }
        if (j == mp->composition[i].ncatom) 
        {
            mp->composition[i].ncatom++;
            srenew(mp->composition[i].catom,mp->composition[i].ncatom);
            mp->composition[i].catom[j].cnumber = 0;
            mp->composition[i].catom[j].cname = strdup(atomname);
            mp->composition[i].catom[j].cnumber += natom;
            qsort(mp->composition[i].catom,mp->composition[i].ncatom,
                  sizeof(mp->composition[i].catom[0]),catom_comp);
        }
        else
            mp->composition[i].catom[j].cnumber += natom;
    }
    else
        gmx_incons("gmx_molprop_add_composition_atom called with invalid arguments");
}

void gmx_molprop_replace_composition_atom(gmx_molprop_t mp,char *compname,
                                          char *oldatom,char *newatom)
{
    int i=-1,j,k;
  
    if (compname) 
    {
        for(i=0; (i<mp->ncomposition); i++)
            if(strcasecmp(compname,mp->composition[i].compname) == 0)
                break;
    
        if (i >= 0) 
        {
            for(j=0; (j<mp->composition[i].ncatom); j++) 
            {
                if (strcasecmp(mp->composition[i].catom[j].cname,oldatom) == 0) 
                {
                    sfree(mp->composition[i].catom[j].cname);
                    mp->composition[i].catom[j].cname = strdup(newatom);
                    break;
                }
            }
        }
    }
}

int gmx_molprop_get_composition_atom(gmx_molprop_t mp,char *compname,
                                     char **atomname,int *natom)
{
    t_composition *cc;
    int i,nc;
  
    for(i=0; (i<mp->ncomposition); i++)
        if (strcasecmp(mp->composition[i].compname,compname) == 0)
            break;
    if (i < mp->ncomposition) 
    {
        cc = &(mp->composition[i]);
        nc = cc->ncatom_c;
        if (nc < cc->ncatom) 
        {
            assign_str(atomname,cc->catom[nc].cname);
            assign_scal(natom,cc->catom[nc].cnumber);
            cc->ncatom_c++;
            return 1;
        }
        cc->ncatom_c = 0;
    }
  
    return 0;
}

static int ref2index(int ref,int *eMPtype)
{
    if (ref < 0) 
    {
        *eMPtype = eMOLPROP_Calc;
        return -ref-1;
    }
    else if (ref > 0) 
    {
        *eMPtype = eMOLPROP_Exp;
        return ref-1;
    }
    else 
        gmx_fatal(FARGS,"Internal error: ref = %d",ref);
        
    *eMPtype = 0;
    return 0; /* To make compiler happy */
}

static int index2ref(int index,int eMPtype)
{
    switch (eMPtype) 
    {
    case eMOLPROP_Exp:
        return index+1;
    case eMOLPROP_Calc:
        return -(index+1);
    default:
        gmx_fatal(FARGS,"Internal error: eMPtype = %d",eMPtype);
    }
    return 0; /* To make compiler happy */
}

static void init_ce(t_calc_exp *ce,const char *program,const char *method,
                    const char *basisset,const char *reference,const char *conformation)
{
    int j;
    
    assign_str(&ce->program,program);
    assign_str(&ce->method,method);
    assign_str(&ce->basisset,basisset);
    assign_str(&ce->reference,reference);
    assign_str(&ce->conformation,conformation);
    ce->natom = 0;
    ce->natom_c = 0;
    ce->catom = NULL;
    ce->polar = NULL;
    ce->dipole = NULL;
    ce->energy = NULL;
    ce->quadrupole = NULL;
    ce->potential = NULL;
    for(j=0; (j<empNR); j++) {
        ce->nprop[j]   = 0;
        ce->nprop_c[j] = 0;
    }
}

void gmx_molprop_add_experiment(gmx_molprop_t mp,const char *reference,
                                const char *conformation,int *expref)
{
    int i,j;
  
    for(i=0; (i<mp->nexper); i++) 
    {
        if ((strcasecmp(mp->exper[i].reference,reference) == 0) &&
            (strcasecmp(mp->exper[i].conformation,conformation) == 0))
            break;
    }
    *expref = index2ref(i,eMOLPROP_Exp);
    if (i == mp->nexper) 
    {						    
        srenew(mp->exper,++mp->nexper);
        init_ce(&(mp->exper[mp->nexper-1]),NULL,NULL,NULL,
                reference,conformation);
    }
}

int gmx_molprop_get_experiment(gmx_molprop_t mp,char **reference,
                               char **conformation,int *expref)
{
    if (mp->nexper_c < mp->nexper) 
    {
        assign_str(reference,mp->exper[mp->nexper_c].reference);
        assign_str(conformation,mp->exper[mp->nexper_c].conformation);
        *expref = index2ref(mp->nexper_c,eMOLPROP_Exp);
        mp->nexper_c++;
    
        return 1;
    }
    *expref      = index2ref(0,eMOLPROP_Exp);
    mp->nexper_c = 0;
  
    return 0;
}

static void ce_add_energy(t_calc_exp *ce,
                          char *type,char *unit,
                          double value,double error)
{
    int i = ce->nprop[empENERGY]++;
    
    srenew(ce->energy,i+1);
    ce->energy[i].type = strdup(type);
    ce->energy[i].unit = strdup(unit);
    ce->energy[i].value = value;
    ce->energy[i].error = error;
}

void gmx_molprop_add_energy(gmx_molprop_t mp,int ref,
                            char *type,char *unit,
                            double value,double error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMP == eMOLPROP_Calc)
        ce_add_energy(&mp->calc[i],type,unit,value,error);
    else
        ce_add_energy(&mp->exper[i],type,unit,value,error);
}
	
static int ce_get_energy(t_calc_exp *ce,
                          char **type,char **unit,
                          double *value,double *error)
{
    int i = ce->nprop_c[empENERGY]++;

    if (i < ce->nprop[empENERGY]) 
    {
        assign_str(type,ce->energy[i].type);
        assign_str(unit,ce->energy[i].unit);
        assign_scal(value,ce->energy[i].value);
        assign_scal(error,ce->energy[i].error);
        return 1;
    }
    else 
    {
        ce->nprop_c[empENERGY] = 0;
        return 0;
    }
}

int gmx_molprop_get_energy(gmx_molprop_t mp,int ref,
                           char **type,char **unit,
                           double *value,double *error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMOLPROP_Calc == eMP) 
        return ce_get_energy(&(mp->calc[i]),type,unit,value,error);
    else
        return ce_get_energy(&(mp->exper[i]),type,unit,value,error);
}

static void ce_add_polar(t_calc_exp *ce,
                         const char *type,const char *unit,
                         double x,double y,double z,
                         double aver,double error)
{
    int i = ce->nprop[empPOLARIZABILITY]++;

    srenew(ce->polar,i+1);
    ce->polar[i].type  = strdup(type);
    ce->polar[i].unit  = strdup(unit);
    ce->polar[i].x     = x;
    ce->polar[i].y     = y;
    ce->polar[i].z     = z;
    ce->polar[i].aver  = aver;
    ce->polar[i].error = error;
}

void gmx_molprop_add_polar(gmx_molprop_t mp,int ref,
                           const char *type,const char *unit,
                           double x,double y,double z,
                           double aver,double error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMP == eMOLPROP_Calc)
        ce_add_polar(&mp->calc[i],type,unit,x,y,z,aver,error);
    else 
        ce_add_polar(&mp->exper[i],type,unit,x,y,z,aver,error);
}

static int ce_get_polar(t_calc_exp *ce,
                         char **type,char **unit,
                         double *x,double *y,double *z,
                         double *aver,double *error)
{
    int i = ce->nprop_c[empPOLARIZABILITY]++;
    
    if (i < ce->nprop[empPOLARIZABILITY])
    {
        assign_str(type,ce->polar[i].type);
        assign_str(unit,ce->polar[i].unit);
        assign_scal(x,ce->polar[i].x);
        assign_scal(y,ce->polar[i].y);
        assign_scal(z,ce->polar[i].z);
        assign_scal(aver,ce->polar[i].aver);
        assign_scal(error,ce->polar[i].error);
        return 1;
    }
    else 
    {
        ce->nprop_c[empPOLARIZABILITY] = 0;
        return 0;
    }
}

int gmx_molprop_get_polar(gmx_molprop_t mp,int ref,
                          char **type,char **unit,
                          double *x,double *y,double *z,
                          double *aver,double *error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMP == eMOLPROP_Calc) 
        return ce_get_polar(&(mp->calc[i]),type,unit,x,y,z,aver,error);
    else
        return ce_get_polar(&(mp->exper[i]),type,unit,x,y,z,aver,error);
}
	
void gmx_molprop_add_potential(gmx_molprop_t mpt,int ref,
                               char *xyz_unit,char *V_unit,
                               int espid,
                               double x,double y,double z,
                               double V)
{
    int        i,np,eMP;
    t_calc_exp *ce;
        
    i = ref2index(ref,&eMP);

    if (eMOLPROP_Calc == eMP)
    {
        ce = &(mpt->calc[i]);
        np = ce->nprop[empPOTENTIAL]++;
        srenew(ce->potential,ce->nprop[empPOTENTIAL]);
        ce->potential[np].xyz_unit  = strdup(xyz_unit);
        ce->potential[np].V_unit    = strdup(V_unit);
        ce->potential[np].espid = espid;
        ce->potential[np].x = x;
        ce->potential[np].y = y;
        ce->potential[np].z = z;
        ce->potential[np].V = V;
    }
}

int gmx_molprop_get_potential(gmx_molprop_t mpt,int ref,
                              char **xyz_unit,char **V_unit,
                              int *espid,
                              double *x,double *y,double *z,
                              double *V)
{
    int i,eMP;
    t_calc_exp *ce;
    
    i = ref2index(ref,&eMP);
    if (eMOLPROP_Calc == eMP)
    {
        ce = &(mpt->calc[i]);

        i = ce->nprop_c[empPOTENTIAL]++;
    
        if (i < ce->nprop[empPOTENTIAL]) 
        {
            assign_str(xyz_unit,ce->potential[i].xyz_unit);
            assign_str(V_unit,ce->potential[i].V_unit);
            assign_scal(espid,ce->potential[i].espid);
            assign_scal(x,ce->potential[i].x);
            assign_scal(y,ce->potential[i].y);
            assign_scal(z,ce->potential[i].z);
            assign_scal(V,ce->potential[i].V);
            return 1;
        }
        else 
        {
            ce->nprop_c[empPOTENTIAL] = 0;
            return 0;
        }
    }
    return 0;
}

static void ce_add_dipole(t_calc_exp *ce,		     
                          char *type,char *unit,
                          double x,double y,double z,
                          double aver,double error)  
{
    int i;
    
    i = ce->nprop[empDIPOLE]++;
    srenew(ce->dipole,i+1);
    ce->dipole[i].type  = strdup(type);
    ce->dipole[i].unit  = strdup(unit);
    ce->dipole[i].x     = x;
    ce->dipole[i].y     = y;
    ce->dipole[i].z     = z;
    ce->dipole[i].aver  = aver;
    ce->dipole[i].error = error;
}

void gmx_molprop_add_dipole(gmx_molprop_t mp,int ref,
                            char *type,char *unit,
                            double x,double y,double z,
                            double aver,double error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);

    if (eMP == eMOLPROP_Calc)
        ce_add_dipole(&(mp->calc[i]),type,unit,x,y,z,aver,error);
    else
        ce_add_dipole(&(mp->exper[i]),type,unit,x,y,z,aver,error);
}

static int ce_get_dipole(t_calc_exp *ce,
                         char **type,char **unit,
                         double *x,double *y,double *z,
                         double *aver,double *error)
{
    int i = ce->nprop_c[empDIPOLE]++;
    
    if (i < ce->nprop[empDIPOLE]) 
    {
        assign_str(type,ce->dipole[i].type);
        assign_str(unit,ce->dipole[i].unit);
        assign_scal(x,ce->dipole[i].x);
        assign_scal(y,ce->dipole[i].y);
        assign_scal(z,ce->dipole[i].z);
        assign_scal(aver,ce->dipole[i].aver);
        assign_scal(error,ce->dipole[i].error);
        return 1;
    }
    else 
    {
        ce->nprop_c[empDIPOLE] = 0;
        return 0;
    }
}

int gmx_molprop_get_dipole(gmx_molprop_t mp,int ref,
                           char **type,char **unit,
                           double *x,double *y,double *z,
                           double *aver,double *error)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMOLPROP_Calc == eMP)
        return ce_get_dipole(&(mp->calc[i]),type,unit,x,y,z,aver,error);
    else
        return ce_get_dipole(&(mp->exper[i]),type,unit,x,y,z,aver,error);
}

static void ce_add_quadrupole(t_calc_exp *ce,
                              char *type,char *unit,
                              double xx,double yy,double zz,
                              double xy,double xz,double yz)
{
    int i,j;
    
    i = ce->nprop[empQUADRUPOLE]++;
    srenew(ce->quadrupole,i+1);
    ce->quadrupole[i].type  = strdup(type);
    ce->quadrupole[i].unit  = strdup(unit);
    ce->quadrupole[i].xx    = xx;
    ce->quadrupole[i].yy    = yy;
    ce->quadrupole[i].zz    = zz;
    ce->quadrupole[i].xy    = xy;
    ce->quadrupole[i].xz    = xz;
    ce->quadrupole[i].yz    = yz;
}

void gmx_molprop_add_quadrupole(gmx_molprop_t mp,int ref,
                                char *type,char *unit,
                                double xx,double yy,double zz,
                                double xy,double xz,double yz)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);

    if (eMP == eMOLPROP_Calc)
        ce_add_quadrupole(&(mp->calc[i]),type,unit,xx,yy,zz,xy,xz,yz);
    else
        ce_add_quadrupole(&(mp->exper[i]),type,unit,xx,yy,zz,xy,xz,yz);
}

static int ce_get_quadrupole(t_calc_exp *ce,
                             char **type,char **unit,
                             double *xx,double *yy,double *zz,
                             double *xy,double *xz,double *yz)
{
    int i = ce->nprop_c[empQUADRUPOLE]++;
    
    if (i < ce->nprop[empQUADRUPOLE]) 
    {
        assign_str(type,ce->quadrupole[i].type);
        assign_str(unit,ce->quadrupole[i].unit);
        assign_scal(xx,ce->quadrupole[i].xx);
        assign_scal(yy,ce->quadrupole[i].yy);
        assign_scal(zz,ce->quadrupole[i].zz);
        assign_scal(xy,ce->quadrupole[i].xy);
        assign_scal(xz,ce->quadrupole[i].xz);
        assign_scal(yz,ce->quadrupole[i].yz);
        
        return 1;
    }
    else 
    {
        ce->nprop_c[empQUADRUPOLE] = 0;
        return 0;
    }
}

int gmx_molprop_get_quadrupole(gmx_molprop_t mp,int ref,
                               char **type,char **unit,
                               double *xx,double *yy,double *zz,
                               double *xy,double *xz,double *yz)
{
    int i,eMP;
    
    i = ref2index(ref,&eMP);
    if (eMOLPROP_Calc == eMP)
        return ce_get_quadrupole(&(mp->calc[i]),type,unit,xx,yy,zz,xy,xz,yz);
    else
        return ce_get_quadrupole(&(mp->exper[i]),type,unit,xx,yy,zz,xy,xz,yz);
}

void gmx_molprop_reset_experiment(gmx_molprop_t mp)
{
    int i,j;
    
    mp->nexper_c = 0;
    for(i=0; (i<mp->nexper); i++)
        for(j=0; (j<empNR); j++)
            mp->exper[i].nprop_c[j] = 0;
}

void gmx_molprop_reset_calculation(gmx_molprop_t mp)
{
    int i,j;
    
    for(i=0; (i<mp->ncalc); i++)
    {
        for(j=0; (j<empNR); j++)
            mp->calc[i].nprop_c[j] = 0;
        mp->calc[i].natom_c = 0;
    }
    mp->ncalc_c = 0;
}

void gmx_molprop_reset(gmx_molprop_t mp)
{
    gmx_molprop_reset_experiment(mp);
    gmx_molprop_reset_calculation(mp);
    gmx_molprop_reset_category(mp);
}

void gmx_molprop_add_category(gmx_molprop_t mp,char *category)
{
    int i;
  
    for(i=0; (i<mp->ncategory); i++) 
        if (strcasecmp(mp->category[i],category) == 0)
            break;
    if (i == mp->ncategory) 
    {
        mp->ncategory++;
        srenew(mp->category,mp->ncategory);
        mp->category[i] = strdup(category);
    }
} 

void gmx_molprop_reset_category(gmx_molprop_t mp)
{
    mp->ncateg_c = 0;
}

void gmx_molprop_reset_composition(gmx_molprop_t mp)
{
    mp->ncomp_c = 0;
}

char *gmx_molprop_get_category(gmx_molprop_t mp)
{
    int i;
  
    if (mp->ncateg_c < mp->ncategory) 
    {
        mp->ncateg_c++;
        return mp->category[mp->ncateg_c-1];
    }
    gmx_molprop_reset_category(mp);
  
    return NULL;
}

int gmx_molprop_search_category(gmx_molprop_t mp,char *catname)
{
    int i;
  
    for(i=0; (i<mp->ncategory); i++)
        if (strcasecmp(mp->category[i],catname) == 0)
            return 1;
  
    return 0;
}

static void merge_props(gmx_molprop_t dst,gmx_molprop_t src,int dstref,int srcref)
{
    double value,error,x,y,z,xx,yy,zz,xy,xz,yz,V;
    int    espid;
    char   *type,*unit;
   
    while(gmx_molprop_get_energy(src,srcref,&type,&unit,&value,&error) == 1) 
    {
        gmx_molprop_add_energy(dst,dstref,type,unit,value,error);
        sfree(type);
        sfree(unit);
    }
    while(gmx_molprop_get_dipole(src,srcref,&type,&unit,&x,&y,&z,&value,&error) == 1) 
    {
        gmx_molprop_add_dipole(dst,dstref,type,unit,x,y,z,value,error);
        sfree(type);
        sfree(unit);
    }
    while(gmx_molprop_get_quadrupole(src,srcref,&type,&unit,&xx,&yy,&zz,&xy,&xz,&yz) == 1) 
    {
        gmx_molprop_add_quadrupole(dst,dstref,type,unit,xx,yy,zz,xy,xz,yz);
        sfree(type);
        sfree(unit);
    }
    while(gmx_molprop_get_polar(src,srcref,&type,&unit,&x,&y,&z,&value,&error) == 1) 
    {
        gmx_molprop_add_polar(dst,dstref,type,unit,x,y,z,value,error);
        sfree(type);
        sfree(unit);
    }
    while(gmx_molprop_get_potential(src,srcref,&type,&unit,&espid,
                                    &x,&y,&z,&V) == 1) 
    {
        gmx_molprop_add_potential(dst,dstref,type,unit,espid,x,y,z,V);
        sfree(type);
        sfree(unit);
    }
}

void gmx_molprop_merge(gmx_molprop_t dst,gmx_molprop_t src)
{
    gmx_molprop *ddd = (gmx_molprop *)dst;
    gmx_molprop *sss = (gmx_molprop *)src;
    int i,nd,ns,atomid,ndcomp,dstref,srcref,datomref,satomref,mult,smult;
    char *tmp,*name,*type,*program,*method,*basisset,*catom,*formula;
    int cnumber,charge,multiplicity;
    char *value,*error,*unit,*reference,*atomname,*coords,*conformation;
    char *stmp,*dtmp;
    double mass,x,y,z,q,sq;
    
    gmx_molprop_reset(src);
    gmx_molprop_reset(dst);
    
    while ((tmp = gmx_molprop_get_category(src)) != NULL) 
    {
        gmx_molprop_add_category(dst,tmp);
    }
    gmx_molprop_set_formula(dst,gmx_molprop_get_formula(src));
    gmx_molprop_set_mass(dst,gmx_molprop_get_mass(src));
    mult = gmx_molprop_get_multiplicity(dst);
    if (mult <= 1)
        gmx_molprop_set_multiplicity(dst,gmx_molprop_get_multiplicity(src));
    else {
        smult = gmx_molprop_get_multiplicity(src);
        if ((NULL != debug) && (smult != mult))
            fprintf(debug,"Not overriding multiplicity to %d when merging since it is %d (%s)\n",
                    smult,mult,gmx_molprop_get_molname(src));
    }
    q = gmx_molprop_get_charge(dst);
    if (q == 0)
        gmx_molprop_set_charge(dst,gmx_molprop_get_charge(src));
    else {
        sq = gmx_molprop_get_charge(src);
        if ((NULL != debug) && (sq != q))
            fprintf(debug,"Not overriding charge to %g when merging since it is %g (%s)\n",
                    sq,q,gmx_molprop_get_molname(src));
    }
            
    stmp = gmx_molprop_get_iupac(src);
    dtmp = gmx_molprop_get_iupac(dst);
    if (((NULL == dtmp) || ((NULL != dtmp) && (strlen(dtmp) == 0))) &&
        (NULL != stmp))
        gmx_molprop_set_iupac(dst,stmp);
    stmp = gmx_molprop_get_cas(src);
    dtmp = gmx_molprop_get_cas(dst);
    if (((NULL == dtmp) || ((NULL != dtmp) && (strlen(dtmp) == 0))) &&
        (NULL != stmp))
        gmx_molprop_set_cas(dst,stmp);
    stmp   = gmx_molprop_get_cid(src);
    dtmp   = gmx_molprop_get_cid(dst);
    if (((NULL == dtmp) || ((NULL != dtmp) && (strlen(dtmp) == 0))) &&
        (NULL != stmp))
        gmx_molprop_set_cid(dst,stmp);
    stmp = gmx_molprop_get_inchi(src);
    dtmp = gmx_molprop_get_inchi(dst);
    if (((NULL == dtmp) || ((NULL != dtmp) && (strlen(dtmp) == 0))) &&
        (NULL != stmp))
        gmx_molprop_set_inchi(dst,stmp);
  
    while (gmx_molprop_get_experiment(src,&reference,&conformation,&srcref) == 1) 
    {
        gmx_molprop_add_experiment(dst,reference,conformation,&dstref);
        merge_props(dst,src,dstref,srcref);
        sfree(reference);
        sfree(conformation);
    }

    while (gmx_molprop_get_calculation(src,&program,&method,
                                       &basisset,&reference,&conformation,&srcref) == 1) 
    {
        gmx_molprop_add_calculation(dst,program,method,basisset,
                                    reference,conformation,&dstref);
        merge_props(dst,src,dstref,srcref);
        while (gmx_molprop_calc_get_atom(src,srcref,&atomname,&atomid,&satomref) == 1) 
        {
            gmx_molprop_calc_add_atom(dst,dstref,atomname,atomid,&datomref);
            if (gmx_molprop_calc_get_atomcoords(src,srcref,satomref,&unit,&x,&y,&z) == 1)
            {
                gmx_molprop_calc_set_atomcoords(dst,dstref,datomref,unit,x,y,z);
                sfree(unit);
            }
            while (gmx_molprop_calc_get_atomcharge(src,srcref,satomref,&type,&unit,&q) == 1)
            {
                gmx_molprop_calc_set_atomcharge(dst,dstref,datomref,type,unit,q);
                sfree(type);
                sfree(unit);
            }
            
            sfree(atomname);
        }
        sfree(program);
        sfree(method);
        sfree(basisset);
        sfree(reference);
    }

  
    ndcomp = 0;
    for(i=0; (i<ddd->ncomposition); i++)
        ndcomp += ddd->composition[i].ncatom;
    if (ndcomp == 0) {
        while ((tmp = gmx_molprop_get_composition(src)) != NULL) 
        {
            gmx_molprop_add_composition(dst,tmp);
            while ((gmx_molprop_get_composition_atom(src,tmp,&catom,&cnumber)) == 1) 
            {
                gmx_molprop_add_composition_atom(dst,tmp,catom,cnumber);
                sfree(catom);
            }
        }
    }
    else 
    {
        if (NULL != debug)
            fprintf(debug,"Dst for %s (%s) already contains %d composition entries (source has %d).\n"
                    "Not adding source composition elements\n",
                    sss->molname,ddd->formula,ddd->ncomposition,sss->ncomposition);
    }
}

gmx_molprop_t gmx_molprop_copy(gmx_molprop_t mp)
{
    gmx_molprop_t dst = gmx_molprop_init();
  
    gmx_molprop_set_molname(dst,gmx_molprop_get_molname(mp));
  
    gmx_molprop_merge(dst,mp);
  
    return dst;
}

void gmx_molprop_add_calculation(gmx_molprop_t mp,const char *program,char *method,
                                 char *basisset,char *reference,
                                 char *conformation,int *calcref)
{
    int j;
    
    srenew(mp->calc,++mp->ncalc);
    init_ce(&(mp->calc[mp->ncalc-1]),program,method,basisset,
            reference,conformation);
    *calcref = index2ref(mp->ncalc-1,eMOLPROP_Calc);
}

void gmx_molprop_calc_add_atom(gmx_molprop_t mp,int calcref,
                               char *atomname,int atomid,int *atomref)
{
    t_calc_exp *calc;
    int index,emp;

    index = ref2index(calcref,&emp);
    calc = &(mp->calc[index]);

    srenew(calc->catom,++calc->natom);
    calc->catom[calc->natom-1].name   = strdup(atomname);
    calc->catom[calc->natom-1].atomid = atomid;
    calc->catom[calc->natom-1].x    = 0;
    calc->catom[calc->natom-1].y    = 0;
    calc->catom[calc->natom-1].z    = 0;
    calc->catom[calc->natom-1].unit = NULL;
    calc->catom[calc->natom-1].nq   = 0;
    calc->catom[calc->natom-1].nq_c = 0;
    calc->catom[calc->natom-1].q    = NULL;
    *atomref = index2ref(calc->natom-1,eMOLPROP_Exp);
}

int gmx_molprop_get_calculation(gmx_molprop_t mp,char **program,
                                char **method,char **basisset,
                                char **reference,char **conformation,int *calcref)
{
    if (mp->ncalc_c < mp->ncalc) 
    {
        assign_str(program,mp->calc[mp->ncalc_c].program);
        assign_str(method,mp->calc[mp->ncalc_c].method);
        assign_str(basisset,mp->calc[mp->ncalc_c].basisset);
        assign_str(reference,mp->calc[mp->ncalc_c].reference);
        assign_str(conformation,mp->calc[mp->ncalc_c].conformation);
        *calcref = index2ref(mp->ncalc_c,eMOLPROP_Calc);
        mp->ncalc_c++;
    
        return 1;
    }
    *calcref = 0;
    mp->ncalc_c = 0;
  
    return 0;
}

int gmx_molprop_calc_get_atom(gmx_molprop_t mp,int calcref,char **atomname,int *atomid,
                              int *atomref)
{
    t_calc_exp *calc;
    int index,emp;
  
    index = ref2index(calcref,&emp);
    calc  = &(mp->calc[index]);
    if (calc->natom_c < calc->natom) 
    {
        assign_str(atomname,calc->catom[calc->natom_c].name);
        assign_scal(atomid,calc->catom[calc->natom_c].atomid);
        *atomref = index2ref(calc->natom_c,eMOLPROP_Exp);
        calc->natom_c++;
      
        return 1;
    }
    calc->natom_c = 0;
    *atomref = 0;
  
    return 0;
}

void gmx_molprop_calc_set_atomcoords(gmx_molprop_t mp,int calcref,int atomref,
                                     char *unit,double x,double y,double z)
{
    t_calc_exp *calc;
    int cindex,aindex,emp;
  
    cindex = ref2index(calcref,&emp);
    range_check(cindex,0,mp->ncalc);
    aindex = ref2index(atomref,&emp);
    calc   = &mp->calc[cindex];
    range_check(aindex,0,calc->natom);
  
    calc->catom[aindex].x = x;
    calc->catom[aindex].y = y;
    calc->catom[aindex].z = z;
    if (NULL == unit)
        gmx_fatal(FARGS,"Can not set atomic coordinates without unit");
    calc->catom[aindex].unit = strdup(unit);
}
					   
int gmx_molprop_calc_get_atomcoords(gmx_molprop_t mp,int calcref,int atomref,
                                    char **unit,double *x,double *y,double *z)
{
    t_calc_exp *calc;
    int cindex,aindex,emp;
  
    cindex = ref2index(calcref,&emp);
    range_check(cindex,0,mp->ncalc);
    aindex = ref2index(atomref,&emp);
    calc   = &mp->calc[cindex];
    range_check(aindex,0,calc->natom);

    assign_scal(x,calc->catom[aindex].x);
    assign_scal(y,calc->catom[aindex].y);
    assign_scal(z,calc->catom[aindex].z);
    assign_str(unit,calc->catom[aindex].unit);
    
    return 1;
}
					   
void gmx_molprop_calc_set_atomcharge(gmx_molprop_t mp,int calcref,int atomref,
                                     char *type,char *unit,double q)
{
    t_calc_exp *calc;
    t_calc_atom *ca;
    int cindex,aindex,emp;
  
    cindex = ref2index(calcref,&emp);
    range_check(cindex,0,mp->ncalc);
    aindex = ref2index(atomref,&emp);
    calc   = &mp->calc[cindex];
    range_check(aindex,0,calc->natom);
    ca = &(calc->catom[aindex]);
    srenew(ca->q,++ca->nq);
    ca->q[ca->nq-1].type = strdup(type);
    ca->q[ca->nq-1].unit = strdup(unit);
    ca->q[ca->nq-1].q = q;  
}

int gmx_molprop_calc_get_atomcharge(gmx_molprop_t mp,int calcref,int atomref,
                                    char **type,char **unit,double *q)
{
    t_calc_exp *calc;
    t_calc_atom *ca;
    int cindex,aindex,emp;
  
    cindex = ref2index(calcref,&emp);
    range_check(cindex,0,mp->ncalc);
    aindex = ref2index(atomref,&emp);
    calc   = &mp->calc[cindex];
    range_check(aindex,0,calc->natom);
    ca = &(calc->catom[aindex]);
  
    if (ca->nq_c < ca->nq) 
    {
        assign_str(type,ca->q[ca->nq_c].type);
        assign_str(unit,ca->q[ca->nq_c].unit);
        assign_scal(q,ca->q[ca->nq_c].q);
        ca->nq_c++;
    
        return 1;
    }
    ca->nq_c = 0;
  
    return 0;
}

