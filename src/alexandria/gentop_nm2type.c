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

int nm2type(FILE *fp,char *molname,gmx_poldata_t pd,gmx_atomprop_t aps,
            t_symtab *tab,t_atoms *atoms,gmx_bool bRing[],
            gpp_atomtype_t atype,int *nbonds,t_params *bonds,
            char **gt_atoms,rvec x[],t_pbc *pbc,real th_toler,real phi_toler,
            gentop_vsite_t gvt)
{
    int     cur = 0;
#define   prev (1-cur)
    int     i,j,k,m,n,nresolved,nbh,maxbond,nb,nonebond;
    int     ai,aj,iter=0,maxiter=10,aa[6];
    int     *bbb,type,atomnr_i;
    char    *aname_i,*aname_m,*aname_n,*gt_atom,*gt_type,*elem_i;
    char    **nbhybrid;
    char    *empty_string = "";
    double  mm;
    int     bLinear,bPlanar,geom;
    real    value;
    t_atom  *atom;
    t_param *param;
    enum { egmLinear, egmPlanar, egmRingPlanar, egmTetrahedral, egmNR };
    const char *geoms[egmNR] = { "linear", "planar", "ring-planar", "tetrahedral" };

    snew(atom,1);
    snew(param,1);
    maxbond = 0;
    nonebond = 0;
    for(i=0; (i<atoms->nr); i++) {
        maxbond = max(maxbond,nbonds[i]);
        if (nbonds[i] == 1)
            nonebond++;
    }
    if (debug)
        fprintf(debug,"Max number of bonds per atom in %s is %d.\nThere are %d atoms with one bond.\n",
                molname,maxbond,nonebond);
    snew(bbb,maxbond);
    snew(nbhybrid,maxbond);
    
    /*if (strcasecmp(molname,"2-pentanol") == 0)
      printf("BOE\n");*/
    for(i=0; (i<atoms->nr); i++) {
        atoms->atom[i].type = -NOTSET;
    }
    do {
        nresolved = 0;
  
        for(i=0; (i<atoms->nr); i++) {
            /* In the first iteration do only the atoms with one bond, e.g.
             * hydrogen, or halides, but also carbonyl oxygens or nitrile nitrogen.
             */
            if (1 /*(((iter == 0) && (nbonds[i] == 1)) || (iter > 0))*/ /*&&
                                                                    (atoms->atom[i].type == -NOTSET)*/) {
                aname_i = *atoms->atomname[i];
                if ((elem_i = gmx_atomprop_element(aps,atoms->atom[i].atomnumber)) == NULL)
                    elem_i = empty_string;
                atomnr_i = atoms->atom[i].atomnumber;
                if ((atomnr_i < 0) && (gmx_atomprop_query(aps,epropElement,"???",
                                                          aname_i,&value)))
                    atomnr_i = gmx_nint(value);
	
                if (debug) 
                    fprintf(debug,"%s%d has bonds to: ",aname_i,i+1);
                nb  = 0;
                nbh = 0;
                for(j=0; (j<bonds->nr); j++) {
                    ai = bonds->param[j].AI;
                    aj = bonds->param[j].AJ;
                    if (ai == i) {
                        bbb[nb++] = aj;
                        if (atoms->atom[aj].atomnumber == 1)
                            nbh++;
                    }
                    else if (aj == i) {
                        bbb[nb++] = ai;
                        if (atoms->atom[ai].atomnumber == 1)
                            nbh++;
                    }
                }
                for(j=0; (j<nb); j++) {
                    nbhybrid[j] = gt_atoms[bbb[j]];
                    if (debug)
                        fprintf(debug," %s%d (%s)",*atoms->atomname[bbb[j]],
                                bbb[j]+1,nbhybrid[j]);
                }
                if (debug)
                    fprintf(debug,"\n");
	
                geom = egmTetrahedral;
                if ((nb == 2) && is_linear(x[i],x[bbb[0]],x[bbb[1]],pbc,th_toler)) 
                {
                    geom = egmLinear;
                    if ((iter == 0) && (NULL != gvt))
                        gentop_vsite_add_linear(gvt,bbb[0],i,bbb[1]);
                }
                else if ((nb == 3) && is_planar(x[i],x[bbb[0]],x[bbb[1]],x[bbb[2]],
                                                pbc,phi_toler))
                {
                    if (bRing[i]) 
                        geom = egmRingPlanar;
                    else
                        geom = egmPlanar;
                    if ((iter == 0) && (NULL != gvt))
                        gentop_vsite_add_planar(gvt,i,bbb[0],bbb[1],bbb[2],nbonds);
                }
                if (debug)
                    fprintf(debug,"Geometry = %s (atom %s with %d bonds)\n",
                            geoms[geom],aname_i,nb);
                if ((gt_atom = gmx_poldata_get_atom(pd,elem_i,nb,nbhybrid,
                                                    geoms[geom])) != NULL)
                {
                    if (debug)
                        fprintf(debug,"Selected %s\n",gt_atom);
                    gt_type = gmx_poldata_get_type(pd,gt_atom);
                    gt_atoms[i] = strdup(gt_atom);
                    mm = 0;
                    if (gmx_atomprop_query(aps,epropMass,"???",elem_i,&value)) 
                        mm = value;
                    else 
                        fprintf(stderr,"Can not find a mass for element %s",elem_i);
                        
                    type = get_atomtype_type(gt_type,atype);
                    if (type == NOTSET)
                    {
                        /* Store mass and polarizability in atomtype as well */
                        atom->qB = 0;
                        atom->m  = mm;
                        type    = add_atomtype(atype,tab,atom,gt_type,param,
                                               atom->type,0,0,0,atomnr_i,0,0);
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
        }
        if (debug)
            fprintf(debug,"Iter %d nresolved %d/%d\n",iter,nresolved,atoms->nr);
        iter++;
    } while (iter < maxiter);
    
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
        if (debug)
            fprintf(debug,"Could not resolve atomtypes %s for: %s.\n",
                    missing,molname);
        if (fp)
            fprintf(fp,"Could not resolve atomtypes %s for: %s.\n",
                    missing,molname);
        sfree(missing);
    }
    sfree(bbb);
    sfree(atom);
    sfree(param);
    sfree(nbhybrid);
    
    return nresolved;
}
