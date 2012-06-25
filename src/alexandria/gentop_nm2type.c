/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gentop_nm2type.c,v 1.14 2009/05/17 13:56:55 spoel Exp $
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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "assert.h"
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "smalloc.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "symtab.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "bondf.h"
#include "names.h"
#include "atomprop.h"
#include "poldata.h"
#include "gentop_nm2type.h"
#include "gentop_vsite.h"

enum { egmLinear, egmPlanar, egmRingPlanar, egmTetrahedral, egmNR };
static const char *geoms[egmNR] = { "linear", "planar", "ring-planar", "tetrahedral" };

static gmx_bool is_planar(rvec xi,rvec xj,rvec xk,rvec xl,t_pbc *pbc,
                      real phi_toler)
{
    int t1,t2,t3;
    rvec r_ij,r_kj,r_kl,m,n;
    real sign,phi;
  
    phi = RAD2DEG*dih_angle(xi,xj,xk,xl,pbc,r_ij,r_kj,r_kl,m,n,&sign,&t1,&t2,&t3);
		  
    return (fabs(phi) < phi_toler);
}

static gmx_bool is_linear(rvec xi,rvec xj,rvec xk,t_pbc *pbc,
                      real th_toler)
{
    int t1,t2,t3;
    rvec r_ij,r_kj;
    real costh,sign,th;
  
    th = fabs(RAD2DEG*bond_angle(xi,xj,xk,pbc,r_ij,r_kj,&costh,&t1,&t2));
		  
    return (th > th_toler) || (th < 180-th_toler);
}

typedef struct {
    int    nb,geom,atomnr,nat,myat,nn;
    int    *bbb;
    double nbo;
    double *valence;
    char   *elem,*aname;
    char   **nb_hybrid,**at_ptr,**tp;
} t_atsel;

static void calc_bondorder(FILE *fp,
                           gmx_poldata_t pd,int natoms,t_atsel ats[],t_params *bonds,
                           double bondorder[])
{
    int    ai,aj,j,lu;
    double toler,vtol,tol,dist,bo;
    char   *unit;
    
    unit = gmx_poldata_get_length_unit(pd);
    lu = string2unit(unit);
    /* Tolerance is 0.002 nm */
    toler = gmx2convert(0.002,lu);
    vtol = 1e-4;

    for(j=0; (j<bonds->nr); j++) 
    {
        ai = bonds->param[j].AI;
        aj = bonds->param[j].AJ;

        if ((NULL != ats[ai].tp[ats[ai].myat]) && (NULL != ats[aj].tp[ats[aj].myat]))
        {
            dist = gmx2convert(bonds->param[j].c[0],lu);
            tol = toler;
            bo = 0;
            while ((0 == bo) && (tol < 10*toler)) 
            {
                bo = gmx_poldata_atype_bondorder(pd,
                                                 ats[ai].tp[ats[ai].myat],
                                                 ats[aj].tp[ats[aj].myat],dist,tol);
                tol *= 2;
            }
            bondorder[j] = bo;
        }
        if (NULL != fp)
        {
            fprintf(fp,"Bond %s%d-%s%d order %g dist %g %s\n",
                    ats[ai].tp[ats[ai].myat],ai+1,
                    ats[aj].tp[ats[aj].myat],aj+1,
                    bondorder[j],dist,unit);
        }
    }
}

static double minimize_valence(FILE *fp,
                               int natoms,t_atsel ats[],const char *molname,
                               t_params *bonds,gmx_poldata_t pd,
                               double bondorder[])
{
    int    iter,iter2,maxiter=10;
    int    i,j,ai,aj,ntot,j_best,jjj,nnull;
    double valence,ddval,dval,dval_best;
    char   buf[1024],b2[64];
    
    dval_best = natoms*4;
    for(iter=0; (iter<maxiter) && (dval_best > 0); iter++) 
    {
        nnull = 1;
        for(iter2=0; (iter2<maxiter) && (nnull > 0); iter2++)
        {
            nnull = 0;
            for(i=0; (i<natoms); i++) 
            {
                for(j=0; (j<ats[i].nb); j++) 
                {
                    aj = ats[i].bbb[j];
                    if (NULL != ats[aj].at_ptr)
                        ats[i].nb_hybrid[j] = ats[aj].at_ptr[ats[aj].myat];
                    else
                        ats[i].nb_hybrid[j] = ats[aj].elem;
                }
                
                ats[i].at_ptr = gmx_poldata_get_atoms(pd,ats[i].elem,ats[i].nb,
                                                      ats[i].nb_hybrid,
                                                      geoms[ats[i].geom]);
                if ((NULL == ats[i].at_ptr) || (NULL == ats[i].at_ptr[0]))
                    nnull++;
            }
        }
        
        if (0 < nnull)
        {
            fprintf(stderr,"Could not find bonding rules to describe all atoms in %s\n",
                    molname);
            return 1;
        }
        /* Select the first possible type */
        ntot = 1;
        for(i=0; (i<natoms); i++) 
        {
            ats[i].nat=0;
            while (NULL != ats[i].at_ptr[ats[i].nat])
            {
                srenew(ats[i].tp,ats[i].nat+1);
                ats[i].tp[ats[i].nat] = gmx_poldata_get_type(pd,ats[i].at_ptr[ats[i].nat]);
                if (0 == gmx_poldata_type_valence(pd,ats[i].tp[ats[i].nat],&valence))
                    gmx_fatal(FARGS,"Can not find valence for %s",
                              ats[i].tp[ats[i].nat]);
                srenew(ats[i].valence,ats[i].nat+1);
                ats[i].valence[ats[i].nat] = valence;
                ats[i].nat++;
            }
            if (0 == ats[i].nat)
                return 1;
            ats[i].myat = 0;
            ntot *= ats[i].nat;
            ats[i].nn = ntot;
        }
        if (NULL != fp)
        {
            fprintf(fp,"There are %d possible combinations of atomtypes in %s\n",
                    ntot,molname);
        }
        if (ntot < 10000) 
        {
            /* Loop over all possibilities */
            j_best = NOTSET;
            for(j=0; (j<ntot) && (dval_best > 0); j++)
            {
                /* Set the atomtypes */
                buf[0] = '\0';
                for(i=0; (i<natoms); i++) 
                {
                    if (i > 0)
                        jjj = (j/ats[i-1].nn);
                    else
                        jjj = j;
                    ats[i].myat = jjj % ats[i].nat;
                    sprintf(b2," %s",ats[i].at_ptr[ats[i].myat]);
                    strcat(buf,b2);
                    ats[i].nbo = 0;
                }
                /* Compute the bond orders */
                calc_bondorder(fp,pd,natoms,ats,bonds,bondorder);
                
                /* Compute the valences based on bond order */
                for(i=0; (i<bonds->nr); i++)
                {
                    ai   = bonds->param[i].a[0];
                    aj   = bonds->param[i].a[1];
                    ats[ai].nbo += bondorder[i];
                    ats[aj].nbo += bondorder[i];
                }
                /* Now check the total deviation from valence */
                dval = 0;
                for(i=0; (i<natoms); i++) 
                {
                    ddval = fabs(ats[i].nbo-ats[i].valence[ats[i].myat]);
                    dval += ddval;
                    if ((ddval > 0) && (NULL != fp))
                    {
                        fprintf(fp,"val for %s%d is %g (should be %g)\n",
                                ats[i].tp[ats[i].myat],i+1,ats[i].nbo,
                                ats[i].valence[ats[i].myat]);
                    }
                }
                if (NULL != fp)
                {
                    fprintf(fp,"j = %d myat = %s dval = %g\n",j,buf,dval);
                }
                if (dval < dval_best)
                {
                    dval_best = dval;
                    j_best = j;
                }
            }
        }
    }
    if (0 < dval_best)
    {
        fprintf(stderr,"Could not find suitable atomtypes for %s. Valence deviation is %g\n",
                molname,dval_best);
    }
    return dval_best;
}

int nm2type(FILE *fp,char *molname,gmx_poldata_t pd,gmx_atomprop_t aps,
            t_symtab *tab,t_atoms *atoms,gmx_bool bRing[],double bondorder[],
            gpp_atomtype_t atype,int *nbonds,t_params *bonds,
            char **gt_atoms,rvec x[],t_pbc *pbc,real th_toler,real phi_toler,
            gentop_vsite_t gvt)
{
    int     cur = 0;
#define   prev (1-cur)
    int     i,j,k,m,n,nonebond,nresolved;
    int     ai,aj,nnull,nra,nrp,aa[6];
    int     type;
    char    *aname_m,*aname_n,*gt_type;
    char    *envptr;
    char    **ptr;
    char    *empty_string = "";
    t_atsel *ats;
    double  mm,dval;
    int     bLinear,bPlanar,imm,iter;
    real    value;
    t_atom  *atom;
    t_param *param;

    envptr = getenv("MOLNAME");
    if ((NULL != envptr) && (0 == strcasecmp(molname,envptr)))
        printf("BOE! Molname = %s\n",molname);
    
    snew(atom,1);
    snew(param,1);
    nonebond = 0;
    snew(ats,atoms->nr);
    for(i=0; (i<atoms->nr); i++) 
    {
        /* Atom information */
        atoms->atom[i].type = NOTSET;
        ats[i].aname = *atoms->atomname[i];
        if ((ats[i].elem = gmx_atomprop_element(aps,atoms->atom[i].atomnumber)) == NULL)
            ats[i].elem = empty_string;
        ats[i].atomnr = atoms->atom[i].atomnumber;
        if ((ats[i].atomnr < 0) && (gmx_atomprop_query(aps,epropElement,"???",
                                                       ats[i].aname,&value)))
        {
            ats[i].atomnr = gmx_nint(value);
            atoms->atom[i].atomnumber = ats[i].atomnr;
        }
        /* Bond information */
        if (nbonds[i] == 1)
            nonebond++;
        snew(ats[i].bbb,nbonds[i]);
        snew(ats[i].nb_hybrid,nbonds[i]);
        ats[i].nb  = 0;
    }
    /* Fill local bonding arrays */
    for(j=0; (j<bonds->nr); j++) {
        ai = bonds->param[j].AI;
        aj = bonds->param[j].AJ;
        ats[ai].bbb[ats[ai].nb++] = aj;
        ats[aj].bbb[ats[aj].nb++] = ai;
    }
    /* Check for consistency in number of bonds */
    for(i=0; (i<atoms->nr); i++) 
    {
        if (ats[i].nb != nbonds[i]) 
        {
            gmx_fatal(FARGS,"Inconsistency between number of nbonds[%d] = %d and bonds structure (%d)",i,nbonds[i],ats[i].nb); 
        }
    }
    for(i=0; (i<atoms->nr); i++) 
    {
        /* Now test initial geometry */
        ats[i].geom = egmTetrahedral;
        if ((ats[i].nb == 2) && is_linear(x[i],x[ats[i].bbb[0]],x[ats[i].bbb[1]],pbc,th_toler)) 
        {
            ats[i].geom = egmLinear;
            if (NULL != gvt)
                gentop_vsite_add_linear(gvt,ats[i].bbb[0],i,ats[i].bbb[1]);
        }
        else if ((ats[i].nb == 3) && is_planar(x[i],x[ats[i].bbb[0]],x[ats[i].bbb[1]],x[ats[i].bbb[2]],
                                               pbc,phi_toler))
        {
            if (bRing[i]) 
                ats[i].geom = egmRingPlanar;
            else
                ats[i].geom = egmPlanar;
            if (NULL != gvt)
                gentop_vsite_add_planar(gvt,i,ats[i].bbb[0],ats[i].bbb[1],ats[i].bbb[2],
                                        nbonds);
        }
    }
    /* Iterate geometry test to look for rings. In a six-ring the information
     * spreads in at most three steps around the ring, so four iterations is sufficient.
     */
    for(iter=0; (iter<4); iter++)
    {
        nra = 0;
        for(i=0; (i<atoms->nr); i++) {
            if ((ats[i].nb == 3) && (ats[i].geom == egmRingPlanar))
            {
                /* If two of three neighbors are in a ring, I am too */
                nrp = 0;
                for(k=0; (k<3); k++) 
                {
                    if ((ats[ats[i].bbb[k]].geom == egmRingPlanar) ||
                        (1 &&
                         (ats[ats[i].bbb[k]].geom == egmTetrahedral) &&
                         (ats[ats[i].bbb[k]].nb == 2)))
                        nrp++;
                }
                if (nrp < 2)
                    ats[i].geom = egmPlanar;
            }
            if (ats[i].geom == egmRingPlanar)
                nra++;
        }
    }
    if (fp)
    {
        fprintf(fp,"There are %d atoms with one bond in %s.\n",
                nonebond,molname);
        fprintf(fp,"There are %d atoms in ring structures in %s.\n",
                nra,molname);
        for(i=0; (i<atoms->nr); i++)
        {
            fprintf(fp,"Geometry = %s (atom %s with %d bonds)\n",
                    geoms[ats[i].geom],ats[i].aname,ats[i].nb);
        }
    }

    /* Now we will determine for each atom the possible atom types in an iterative
     * fashion until convergence (deviation from valence criteria dval == 0).
     */
    dval = minimize_valence(fp,atoms->nr,ats,molname,bonds,pd,bondorder);
    
    nresolved = 0;
    if (0 == dval)
    {
        /* Now set the contents of the atoms structure */
        for(i=0; (i<atoms->nr); i++) 
        {
            gt_atoms[i] = strdup(ats[i].at_ptr[ats[i].myat]);
            gt_type = gmx_poldata_get_type(pd,gt_atoms[i]);
            mm = 0;
            if (gmx_atomprop_query(aps,epropMass,"???",ats[i].elem,&value)) 
                mm = value;
            else 
                fprintf(stderr,"Can not find a mass for element %s",ats[i].elem);
            
            type = get_atomtype_type(gt_type,atype);
            if (type == NOTSET)
            {
                /* Store mass and polarizability in atomtype as well */
                atom->qB = 0;
                atom->m  = mm;
                type    = add_atomtype(atype,tab,atom,gt_type,param,
                                       atom->type,0,0,0,ats[i].atomnr,0,0);
            }
            
            /* atoms->atom[i].q  = 0;*/
            atoms->atom[i].qB = 0;
            atoms->atom[i].m  = atoms->atom[i].mB = mm;
            atoms->atom[i].type  = type;
            atoms->atom[i].typeB = type;
            atoms->atomtype[i] = put_symtab(tab,gt_type);
            
            nresolved++;
        }
    }
    
    /* fprintf(stderr,"\n");*/
    if (nresolved < atoms->nr) {
        char *missing,atom[32];
        snew(missing,12*(atoms->nr-nresolved));
        for(i=0; (i<atoms->nr); i++) {
            if (gt_atoms[i] == NULL) {
                sprintf(atom," %s-%d",*atoms->atomname[i],i+1);
                strcat(missing,atom);
            }
        }
        if (fp)
            fprintf(fp,"Could not resolve atomtypes %s for: %s.\n",
                    missing,molname);
        sfree(missing);
    }
    sfree(atom);
    sfree(param);
            
    return nresolved;
}
