/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_util.c,v 1.51 2009/06/01 06:13:18 spoel Exp $
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
#include <stdlib.h>
#include "smalloc.h"
#include "maths.h"
#include "typedefs.h"
#include "string2.h"
#include "futil.h"
#include "physics.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "atomprop.h"
#include "gpp_atomtype.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
//#include "molprop_tables.h"
#include "gentop_vsite.h"
#include "gentop_core.h"
#include "gentop_comm.h"
//#include "molselect.h"

int get_val(gmx_molprop_t mp,int expref,char *type,int emp,
            double *value,double *error,double vec[3],
            tensor quadrupole) 
{
    int done = 0;
    char *unit,*mytype;
    double vv,ee,x,y,z,xx,yy,zz,xy,xz,yz;
          
    switch (emp) {
    case empENERGY:
        while ((done == 0) && (gmx_molprop_get_energy(mp,expref,&mytype,&unit,&vv,&ee) == 1))
        {
            if ((NULL == type) || (strcasecmp(type,mytype) == 0)) 
            {
                *value = vv;
                *error = ee;
                done = 1;
            }
            sfree(unit);
            sfree(mytype);
        }
        break;
    case empDIPOLE:
        while ((done == 0) && (gmx_molprop_get_dipole(mp,expref,&mytype,&unit,&x,&y,&z,&vv,&ee) == 1))
        {
            if ((NULL == type) || (strcasecmp(type,mytype) == 0)) 
            {
                *value = vv;
                *error = ee;
                vec[XX] = x;
                vec[YY] = y;
                vec[ZZ] = z;
                done = 1;
            }
            sfree(unit);
            sfree(mytype);
        }
        break;
    case empPOLARIZABILITY:
        while ((done == 0) && (gmx_molprop_get_polar(mp,expref,&mytype,&unit,&xx,&yy,&zz,&vv,&ee) == 1))
        {
            if ((NULL == type) || (strcasecmp(type,mytype) == 0)) 
            {
                *value = vv;
                *error = ee;
                vec[XX] = xx;
                vec[YY] = yy;
                vec[ZZ] = zz;
                done = 1;
            }
            sfree(unit);
            sfree(mytype);
        }
        break;
    case empQUADRUPOLE:
        while ((done == 0) && (gmx_molprop_get_quadrupole(mp,expref,&mytype,&unit,&xx,&yy,&zz,&xy,&xz,&yz) == 1))
        {
            if ((NULL == type) || (strcasecmp(type,mytype) == 0)) 
            {
                quadrupole[XX][XX] = xx;
                quadrupole[XX][YY] = xy;
                quadrupole[XX][ZZ] = xz;
                quadrupole[YY][XX] = 0;
                quadrupole[YY][YY] = yy;
                quadrupole[YY][ZZ] = yz;
                quadrupole[ZZ][XX] = 0;
                quadrupole[ZZ][YY] = 0;
                quadrupole[ZZ][ZZ] = zz;
                done = 1;
            }
            sfree(unit);
            sfree(mytype);
        }
        break;
    default:
        break;
    }
    return done;
}

int mp_get_prop_ref(gmx_molprop_t mp,int emp,int iQM,char *lot,
                    char *conf,char *type,double *value,double *error,
                    char **ref,char **mylot,
                    double vec[3],tensor quadrupole)
{
    char *buf,*method,*reference,*name,*basisset,*program,*conformation,*expconf;
    int expref,n,k,done=0;
    char **qm,**ll;
    
    if ((iQM == iqmExp) || (iQM == iqmBoth)) 
    {
        gmx_molprop_reset_experiment(mp);
        reference = NULL;
        while ((done == 0) &&
               (gmx_molprop_get_experiment(mp,&reference,
                                           &expconf,&expref) == 1)) 
        {
            if ((NULL == conf) || (strcasecmp(conf,expconf) == 0))
                done = get_val(mp,expref,type,emp,value,error,vec,quadrupole);
            sfree(expconf);
        }
        gmx_molprop_reset_experiment(mp);
        if (done != 0) {
            if (NULL != ref)
                *ref   = strdup(reference);
            else
                sfree(reference);
            if (NULL != mylot)
                *mylot = strdup("Experiment");
        }
        else {
            if (NULL != reference)
                sfree(reference);
        }
    }
    
    if ((done == 0) && ((iQM == iqmBoth) || (iQM == iqmQM)))
    {
        gmx_molprop_reset_calculation(mp);
        if (NULL != lot) 
        {
            qm = split(':',lot);
            n = 0;
            while ((done == 0) && (NULL != qm[n]))
            {
                ll = split('/',qm[n]);
                if ((NULL != ll[0]) && (NULL != ll[1])) 
                {
                    while ((done == 0) &&
                           (gmx_molprop_get_calculation(mp,&program,&method,&basisset,
                                                        &reference,&conformation,&expref) == 1)) 
                    {
                        if (((strcasecmp(method,ll[0]) == 0) &&
                             (strcasecmp(basisset,ll[1]) == 0)) &&
                            ((NULL == conf) || (strcasecmp(conf, conformation) == 0)))
                        {
                            done = get_val(mp,expref,type,emp,value,error,vec,quadrupole);
                        }
                        sfree(program);
                        sfree(basisset);
                        sfree(method);
                        if ((done != 0) && (NULL != ref))
                            *ref = strdup(reference);
                        sfree(reference);
                        sfree(conformation);
                    }
                }
                k = 0;
                while (ll[k] != NULL)
                {
                    sfree(ll[k]);
                    k++;
                }
                if (done != 0) 
                {
                    if (NULL != mylot)
                        *mylot = strdup(qm[n]);
                }
                n++;
            }
            k = 0;
            while (qm[k] != NULL)
            {
                sfree(qm[k]);
                k++;
            }
        }
        gmx_molprop_reset_calculation(mp);
    }
    return done;
}

int mp_get_prop(gmx_molprop_t mp,int emp,int iQM,char *lot,
                char *conf,char *type,double *value)
{
    double error,vec[3];
    tensor quad;

    return mp_get_prop_ref(mp,emp,iQM,lot,conf,type,value,&error,
                           NULL,NULL,vec,quad);
}                

int gmx_molprop_get_calc_lot(gmx_molprop_t mpt,char *lot)
{
    int  k,expref = 0,done = 0;
    char **ll,*method,*program,*basisset;
    
    ll = split('/',lot);
    if ((NULL != ll[0]) && (NULL != ll[1])) 
    {
        while ((done == 0) &&
               (gmx_molprop_get_calculation(mpt,NULL,&method,&basisset,
                                            NULL,NULL,&expref) == 1)) 
        {
            if (((strcasecmp(method,ll[0]) == 0) &&
                 (strcasecmp(basisset,ll[1]) == 0)))
                done = 1;
            sfree(method);
            sfree(basisset);
        }
    }
    k = 0;
    while (ll[k] != NULL)
    {
        sfree(ll[k]);
        k++;
    }
    return expref;
}

enum { elgcOK, elgcNOATOMS, elgcNOTOP, elgcNR };

static int lo_gen_composition(gmx_molprop_t mp,gmx_poldata_t pd,gmx_atomprop_t aps,
                              gmx_bool *bSpoel,gmx_bool *bMiller,gmx_bool *bBosque,
                              double th_toler,double phi_toler)
{
    t_atoms    *atoms;
    rvec       *x=NULL;
    int        i,j,natom=-1,bOK=FALSE;
    int        nexcl=1;
    char       **anames=NULL,**smnames;
    int        *nbonds,bDone=0,calcref,atomref,atomid,ftb;
    t_params   *plist=NULL;
    t_excls    *excls;
    double     btol=0.2;
    gpp_atomtype_t atype;
    t_symtab   symtab;
    matrix     box;
    t_pbc      pbc;
    char       *program,*method,*basisset,*reference,*conformation,*miller,*elem;
    char       *atomname,*unit,*molname,*type;
    double     xx,yy,zz;
    gentop_vsite_t gvt;
    
    clear_mat(box);
    set_pbc(&pbc,epbcNONE,box);
    molname = gmx_molprop_get_molname(mp);
    while (!bDone && (gmx_molprop_get_calculation(mp,&program,&method,
                                                  &basisset,&reference,
                                                  &conformation,&calcref) == 1)) 
    {
        /* This assumes we have either all atoms or none. A consistency check could be
         * to compare the number of atoms to the formula */
        natom = 0;
        while (gmx_molprop_calc_get_atom(mp,calcref,&atomname,&atomid,&atomref) == 1) 
        {
            srenew(anames,natom+1);
            srenew(x,natom+1);
            anames[natom] = strdup(atomname);
            if (gmx_molprop_calc_get_atomcoords(mp,calcref,atomref,&unit,&xx,&yy,&zz) == 1) 
            {
                int myunit;
                if (NULL == unit)
                    myunit = eg2cAngstrom;
                else
                {
                    myunit = string2unit(unit);
                    if (-1 == myunit)
                        myunit = eg2cAngstrom; 
                }
                x[natom][XX] = convert2gmx(xx,myunit);
                x[natom][YY] = convert2gmx(yy,myunit);
                x[natom][ZZ] = convert2gmx(zz,myunit);
            }
            else
            {
                fprintf(stderr,"Can not find coordinates for atom %s%d in %s",
                        atomname,natom,molname);
            }
            natom++;
            sfree(atomname);
            if (NULL != unit)
                sfree(unit);
        }
        bDone = (natom > 0);
        if (debug && bDone) 
            fprintf(debug,"LO_COMP: %s %s/%s\n",molname,method,basisset);
        sfree(program);
        sfree(basisset);
        sfree(method);
        sfree(reference);
        sfree(conformation);
    }
    if (bDone) 
    {
        snew(atoms,1);
        open_symtab(&symtab);
        init_t_atoms(atoms,natom,FALSE);
        for(i=0; (i<atoms->nr); i++) 
        {
            atoms->atomname[i] = put_symtab(&symtab,anames[i]);
            atoms->atom[i].atomnumber = gmx_atomprop_atomnumber(aps,anames[i]);
        }
        snew(nbonds,atoms->nr);
        snew(smnames,atoms->nr);
        snew(plist,F_NRE);
        mk_bonds(pd,atoms,x,NULL,plist,nbonds,TRUE,TRUE,TRUE,
                 nexcl,&excls,FALSE,NULL,aps,btol,TRUE);
        ftb = gmx_poldata_get_bond_ftype(pd);
        /* Setting the atom types: this depends on the bonding */
        gvt = gentop_vsite_init(egvtLINEAR);
        if ((atype = set_atom_type(NULL,molname,&symtab,atoms,&(plist[ftb]),
                                   nbonds,smnames,
                                   pd,aps,x,&pbc,th_toler,phi_toler,gvt)) != NULL) 
            bOK = TRUE;
        gentop_vsite_done(&gvt);
        if (bOK) 
        {
            /* Gather all the atomtypes */
            gmx_molprop_delete_composition(mp,"spoel");
            gmx_molprop_add_composition(mp,"spoel");

            for(i=0; (i<natom); i++) 
            {
                type = get_atomtype_name(atoms->atom[i].type,atype);
                if (NULL == type)
                    return elgcNOTOP;
                gmx_molprop_add_composition_atom(mp,"spoel",type,1);
                if ((miller = gmx_poldata_get_miller_equiv(pd,smnames[i])) != NULL)
                    gmx_molprop_add_composition_atom(mp,"miller",miller,1);
                if ((elem  = gmx_atomprop_element(aps,atoms->atom[i].atomnumber)) != NULL)
                    gmx_molprop_add_composition_atom(mp,"bosque",elem,1);
                sfree(smnames[i]);
            }
            if (debug)
                fprintf(debug,"\n");
            done_atomtype(atype);
        }
        sfree(smnames);
        close_symtab(&symtab);
        free_symtab(&symtab);
    }
    gmx_molprop_reset_calculation(mp);
    
    *bSpoel  = gmx_molprop_count_composition_atoms(mp,(char *)"spoel",NULL) == natom;
    *bMiller = gmx_molprop_count_composition_atoms(mp,(char *)"miller",NULL) == natom;
    *bBosque = gmx_molprop_count_composition_atoms(mp,(char *)"bosque",NULL) == natom;
    
    if (!bDone)
        return elgcNOATOMS;
    else if (bOK)
        return elgcOK;
    else
        return elgcNOTOP;
}

void generate_composition(int nmol,gmx_molprop_t mp[],gmx_poldata_t pd,
                          gmx_atomprop_t ap,gmx_bool bForceGenComp,
                          double th_toler,double phi_toler)
{
    int    i,j,cnumber,elgc,atomnumber,nok=0,ntest=0,natom;
    char   *compname,*catom,*miller_equiv=NULL,*elem,*alexandria_equiv,*atomname,*molname;
    gmx_bool   bSpoel,bMiller,bBosque;
    double mass;
    real   mm,mh;
    
    (void) gmx_atomprop_query(ap,epropMass,"???","H",&mh);
    for(j=0; (j<nmol); j++) 
    {
        mass    = 0;
        bSpoel  = TRUE;
        bMiller = TRUE;
        bBosque = TRUE;
        
        molname = gmx_molprop_get_molname(mp[j]);
        /*if (strcasecmp("bromomethane",molname) == 0)
          fprintf(stderr,"Psst: %s\n",molname);*/
        gmx_molprop_delete_composition(mp[j],"miller");
        gmx_molprop_delete_composition(mp[j],"bosque");
        gmx_molprop_add_composition(mp[j],"miller");
        gmx_molprop_add_composition(mp[j],"bosque");
        
        elgc = elgcOK;
        if (bForceGenComp || 
            (gmx_molprop_count_composition_atoms(mp[j],(char *)"spoel",NULL) == 0)) {
            elgc = lo_gen_composition(mp[j],pd,ap,&bSpoel,&bMiller,&bBosque,th_toler,phi_toler);
            if (elgc == elgcOK) 
            {
                if (debug) 
                {
                    fprintf(debug,"LO_COMP: ");
                    while(gmx_molprop_get_composition_atom(mp[j],(char *)"spoel",&atomname,&natom) == 1)
                    {
                        fprintf(debug," %s:%d",atomname,natom);
                    }
                    fprintf(debug,"\n");
                }
                nok++;
            }
            else if (debug && (elgc == elgcNOTOP))
                fprintf(debug,"Failed to make composition for %s\n",
                        gmx_molprop_get_molname(mp[j]));
            if (elgc != elgcNOATOMS)
                ntest++;
        }
        else
        {
            /* Let's use an existing definition of spoel atoms */
            while (gmx_molprop_get_composition_atom(mp[j],(char *)"spoel",
                                                    &catom,&cnumber) == 1) {
                if (miller_equiv)
                    miller_equiv = NULL;
                if (1 == gmx_poldata_search_atype(pd,catom,&elem,NULL,NULL,&miller_equiv,
                                                  NULL,NULL,NULL,NULL)) 
                {
                    (void) gmx_atomprop_query(ap,epropMass,"???",elem,&mm);
                    mass += mm*cnumber;
                    /* See if this miller atom exists */
                    if (gmx_poldata_get_miller(pd,miller_equiv,&atomnumber,
                                               NULL,NULL,&alexandria_equiv) != NULL) 
                    {
                        gmx_molprop_add_composition_atom(mp[j],"miller",
                                                         miller_equiv,cnumber);
                    }
                    else 
                        bMiller = FALSE;
                    /* See if this bosque atom exists */
                    if (gmx_poldata_get_bosque(pd,NULL,elem,NULL) != NULL) 
                    {
                        gmx_molprop_add_composition_atom(mp[j],"bosque",
                                                         elem,cnumber);
                    }
                    else 
                        bBosque = FALSE;
                    sfree(elem);
                    sfree(miller_equiv);
                }
                else if (gmx_poldata_get_miller(pd,catom,&atomnumber,
                                                NULL,NULL,&alexandria_equiv) != NULL) 
                {
                    elem = gmx_atomprop_element(ap,atomnumber);
                    (void) gmx_atomprop_query(ap,epropMass,"???",elem,&mm);
                    mass += mm*cnumber;
                    if (alexandria_equiv) 
                    {
                        gmx_molprop_replace_composition_atom(mp[j],(char *)"spoel",
                                                             catom,alexandria_equiv);
                        gmx_molprop_add_composition_atom(mp[j],(char *)"miller",
                                                         catom,cnumber);
                        if (gmx_poldata_get_bosque(pd,NULL,elem,NULL) != NULL) 
                        {
                            gmx_molprop_add_composition_atom(mp[j],(char *)"bosque",
                                                             elem,cnumber);
                        }
                        else 
                            bBosque = FALSE;
                    }
                    else
                        bSpoel = FALSE;
                    sfree(alexandria_equiv);
                }
                else 
                {
                    bSpoel = FALSE;
                    bMiller = FALSE;
                    bBosque = FALSE;
                    fprintf(stderr,"Atom %s not found in poldata for molecule %s\n",
                            catom,gmx_molprop_get_molname(mp[j]));
                }
            }
            gmx_molprop_set_mass(mp[j],mass);
        }
        if (!bMiller)
            gmx_molprop_delete_composition(mp[j],"miller");
        if (!bBosque)
            gmx_molprop_delete_composition(mp[j],"bosque");
        if (!bSpoel) 
            gmx_molprop_delete_composition(mp[j],"spoel");
    }
    printf("Generated composition for %d out of %d molecules.\n",nok,ntest);
}

void generate_formula(int nmol,gmx_molprop_t mp[],gmx_atomprop_t ap)
{
    int  i,j,jj,cnumber,an;
    char formula[1280],number[32];
    char *mform,*compname,*catom;
    int  *ncomp;
    real value;

    for(i=0; (i<nmol); i++) {
        snew(ncomp,110);  
        formula[0] = '\0';
        while (gmx_molprop_get_composition_atom(mp[i],(char *)"bosque",&catom,&cnumber) == 1) {
            if (gmx_atomprop_query(ap,epropElement,"???",catom,&value)) {
                an = gmx_nint(value);
                range_check(an,0,110);
                if (an > 0)
                    ncomp[an] += cnumber;
            }
        }
        if (ncomp[6] > 0) 
        {
            strcat(formula,"C");
            if (ncomp[6] > 1) 
            {
                sprintf(number,"%d",ncomp[6]);
                strcat(formula,number);
            }
            ncomp[6] = 0;
            if (ncomp[1] > 0) 
            {
                strcat(formula,"H");
                if (ncomp[1] > 1) 
                {
                    sprintf(number,"%d",ncomp[1]);
                    strcat(formula,number);
                }
                ncomp[1] = 0;
            }
        }
        for(j=1; (j<110); j++) 
        {
            if (ncomp[j] > 0)
            {
                strcat(formula,gmx_atomprop_element(ap,j));
                if (ncomp[j] > 1) 
                {
                    sprintf(number,"%d",ncomp[j]);
                    strcat(formula,number);
                }
            }
        }
        mform = gmx_molprop_get_formula(mp[i]);
        if (strlen(formula) > 0) 
        {
            if (debug) 
            {
                if (mform && (strcasecmp(formula,mform) != 0))
                    fprintf(debug,"Formula '%s' does match '%s' based on composition for %s.\n",
                            mform,formula,gmx_molprop_get_molname(mp[i]));
            }
            gmx_molprop_set_formula(mp[i],formula);
        }
        else if ((!mform || (strlen(mform) == 0)) && debug)
            fprintf(debug,"Empty composition and formula for %s\n",
                    gmx_molprop_get_molname(mp[i]));
        sfree(ncomp);
    }
}
  
gmx_molprop_t atoms_2_molprop(char *molname,t_atoms *atoms,char **smnames,
                              gmx_atomprop_t ap,gmx_poldata_t pd,gmx_bool bForceGenComp,
                              double th_toler,double ph_toler)
{
    gmx_molprop_t mp;
    int ai,aj,i,j;
    char buf[32],*ptr;
  
    mp = gmx_molprop_init();
    gmx_molprop_set_molname(mp,molname);
    if (!bForceGenComp)
    {
        gmx_molprop_add_composition(mp,"spoel");
        
        for(i=0; (i<atoms->nr); i++) 
            gmx_molprop_add_composition_atom(mp,"spoel",smnames[i],1);
    }
    generate_composition(1,&mp,pd,ap,bForceGenComp,th_toler,ph_toler);
    generate_formula(1,&mp,ap);
    gmx_molprop_set_molname(mp,molname);
  
    return mp;  
}

int molprop_2_atoms(gmx_molprop_t mp,gmx_atomprop_t ap,
                    t_symtab *tab,const char *lot,
                    t_atoms *atoms,const char *q_algorithm,
                    rvec **x)
{
    char   *program,*method,*basisset,*name,*type,*value;
    char   *atomname,*unit,*coords,**aa=NULL,*molnm;
    char   **ptr,**ptr1,buf[1024],format[1024];
    double *q=NULL,qq,xx,yy,zz;
    int    i,nn,done,natom,calcref,atomref,atomid;
  
    (*x)  = NULL;
    ptr1  = split(':',lot);
    nn    = 0;
    done  = 0;
    natom = 0;
    molnm = gmx_molprop_get_molname(mp);
    while ((done == 0) && (ptr1[nn] != NULL))
    {
        ptr = split('/',ptr1[nn]);
        nn++;
        while ((done == 0) &&
               (gmx_molprop_get_calculation(mp,&program,&method,&basisset,
                                            NULL,NULL,&calcref) == 1)) 
        {
            if (NULL != debug)
                fprintf(debug,"LOT = %s/%s molnm = %s\n",method,basisset,molnm);
            if ((strcasecmp(method,ptr[0]) == 0) &&
                (strcasecmp(basisset,ptr[1]) == 0)) 
            {
                natom = 0;
                while (gmx_molprop_calc_get_atom(mp,calcref,&atomname,&atomid,&atomref) == 1) 
                {
                    int myunit;
                    
                    srenew(aa,natom+1);
                    srenew((*x),natom+1);
                    srenew(q,natom+1);
                    q[natom] = 0;
                    gmx_molprop_calc_get_atomcoords(mp,calcref,atomref,&unit,&xx,&yy,&zz);
                    myunit = string2unit(unit);
                    if (-1 == myunit)
                    {
                        myunit = eg2cAngstrom;
                    }
                    (*x)[natom][XX] = convert2gmx(xx,myunit);
                    (*x)[natom][YY] = convert2gmx(yy,myunit);
                    (*x)[natom][ZZ] = convert2gmx(zz,myunit);
                    aa[natom] = strdup(atomname);
                    while (gmx_molprop_calc_get_atomcharge(mp,calcref,atomref,&type,&unit,&qq) == 1) 
                    {
                        int myunit;
                        
                        myunit = string2unit(unit);
                        if (-1 == myunit)
                            myunit = eg2cElectron;
                        if (strcasecmp(type,q_algorithm) == 0) 
                            q[natom] = convert2gmx(qq,myunit);
                        sfree(type);
                        sfree(unit);
                    }
                    natom++;
                    sfree(atomname);
                }
                init_t_atoms(atoms,natom,FALSE);
                for(i=0; (i<natom); i++) {
                    t_atoms_set_resinfo(atoms,i,tab,molnm,1,' ',1,' ');
                    atoms->atomname[i] = put_symtab(tab,aa[i]);
                    atoms->atom[i].atomnumber = gmx_atomprop_atomnumber(ap,aa[i]);
                    atoms->atom[i].q = q[i]; 
                    atoms->atom[i].resind = 0;
                    sfree(aa[i]);
                }
                atoms->nres = 1;
                sfree(q);
                sfree(aa);
                done = 1;
            }
        }
        gmx_molprop_reset_calculation(mp);
    }
    if (NULL != debug)
        fprintf(debug,"Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molnm,lot,natom);
    return natom;
}

static void gmx_molprop_dump(gmx_molprop_t mp,int i,FILE *fp)
{
    char *cat;
  
    if (fp) {
        fprintf(fp,"\nmolprop %d\n",i);
        fprintf(fp,"molname:      %s\n",gmx_molprop_get_molname(mp));
        fprintf(fp,"formula:      %s\n",gmx_molprop_get_formula(mp));
        fprintf(fp,"mass:         %g\n",gmx_molprop_get_mass(mp));
        fprintf(fp,"charge:       %d\n",gmx_molprop_get_charge(mp));
        fprintf(fp,"multiplicity: %d\n",gmx_molprop_get_multiplicity(mp));
        fprintf(fp,"category:    ");
        while ((cat = gmx_molprop_get_category(mp)) != NULL)
            fprintf(fp," %s",cat);
        fprintf(fp,"\n");
    }
}

int my_strcmp(char *a,char *b)
{
    if ((NULL != a) && (NULL != b))
        return strcasecmp(a,b);
    else
        return 1;
}

static void merge_doubles(int *np,gmx_molprop_t mp[],char *doubles,
                          gmx_bool bForceMerge)
{
    FILE *fp;
    int  i,j,nporig=*np,ndouble=0;
    gmx_bool bForm,bName,bDouble;
    char *molname[2];
    char *form[2];
    gmx_molprop_t **mymp;
    int  cur = 0;
#define prev (1-cur)
  
    if (NULL != doubles)
        fp = fopen(doubles,"w");
    else
        fp = NULL;
    snew(mymp,*np);
    for(i=0; (i<*np); i++) 
        mymp[i] = &(mp[i]);
    
    for(i=0; (i<*np); ) {
        bDouble = FALSE;
        gmx_molprop_dump(*(mymp[i]),i,debug);
        molname[cur] = gmx_molprop_get_molname(*(mymp[i]));
        if (strcmp(molname[cur],"benzene") == 0)
            printf("KOEKOEK\n");
        form[cur]    = gmx_molprop_get_formula(*(mymp[i]));
        if (i > 0) {
            bForm = (my_strcmp(form[prev],form[cur]) == 0);
            bName = (my_strcmp(molname[prev],molname[cur]) == 0);
            if (bName) {
                if (!bForm && debug)
                    fprintf(debug,"%s %s with formulae %s - %s\n",
                            bForceMerge ? "Merging molecules" : "Found molecule",
                            molname[prev],form[prev],form[cur]);
                if (bForceMerge || bForm) {
                    if (fp)
                        fprintf(fp,"%5d  %s\n",ndouble+1,molname[prev]);
                    gmx_molprop_merge(*(mymp[i-1]),*(mymp[i]));
                    gmx_molprop_delete(*(mymp[i]));
                    for(j=i+1; (j<*np); j++) {
                        mymp[j-1] = mymp[j];
                    }
                    bDouble = TRUE;
                    ndouble++;
                    (*np)--;
                    molname[cur]  = gmx_molprop_get_molname(*(mymp[i]));
                    form[cur]     = gmx_molprop_get_formula(*(mymp[i]));
                }
            }
        }
        if (!bDouble) {
            i++;
            cur = prev;
        }
    }
    for(i=0; (i<*np); i++) {
        mp[i] = *(mymp[i]);
        gmx_molprop_dump(mp[i],10000+i,debug);
    }
    sfree(mymp);
    if (fp)
        fclose(fp);
    printf("There were %d double entries, leaving %d out of %d after merging.\n",
           ndouble,*np,nporig);
}

static void dump_mp(int np,gmx_molprop_t mp[])
{
    FILE *fp;
    int  i,j,k;
  
    fp = fopen("dump_mp.dat","w");
  
    for(i=0; (i<np); ) {
        for(j=i; (j<np-1) && (strcasecmp(gmx_molprop_get_formula(mp[i]),
                                         gmx_molprop_get_formula(mp[j+1])) == 0); j++)
            ;
        if (j > i) {
            for(k=i; (k<=j); k++)
                fprintf(fp,"%-20s  %s\n",
                        gmx_molprop_get_formula(mp[k]),
                        gmx_molprop_get_molname(mp[k]));
            fprintf(fp,"\n");
        }
        i=j+1;
    }
  
    fclose(fp);
}

gmx_molprop_t *merge_xml(int nfile,char **filens,char *outf,
                         char *sorted,char *doubles,int *nmolprop,
                         gmx_atomprop_t ap,gmx_poldata_t pd,
                         gmx_bool bForceMerge,gmx_bool bForceGenComp,
                         double th_toler,double ph_toler)
{
    gmx_molprop_t *mp=NULL,*mpout=NULL;
    int       i,j,np,npout=0,tmp;
    char      buf[100];
  
    for(i=0; (i<nfile); i++) 
    {
        mp = gmx_molprops_read(filens[i],&np);
        generate_composition(np,mp,pd,ap,bForceGenComp,th_toler,ph_toler);
        generate_formula(np,mp,ap);
        srenew(mpout,npout+np);
        for(j=0; (j<np); j++) 
        {
            mpout[npout+j] = gmx_molprop_copy(mp[j]);
            gmx_molprop_delete(mp[j]);
        }
        sfree(mp);
        npout += np;
    }
    tmp = npout;
    gmx_molprop_sort(npout,mpout,empSORT_Molname,NULL,NULL);
    merge_doubles(&npout,mpout,doubles,bForceMerge);
    printf("There are %d total molecules after merging %d doubles.\n",
           npout,tmp);
      
    if (outf) 
    {
        printf("There are %d entries to store in output file %s\n",npout,outf);
        gmx_molprops_write(outf,npout,mpout,0);
    }
    if (sorted) 
    {
        gmx_molprop_sort(npout,mpout,empSORT_Formula,NULL,NULL);
        gmx_molprops_write(sorted,npout,mpout,0);
        dump_mp(npout,mpout);
    }
  
    *nmolprop = npout; 
  
    return mpout;
}

static int comp_mp_molname(const void *a,const void *b)
{
    gmx_molprop_t *ma = (gmx_molprop_t *) a;
    gmx_molprop_t *mb = (gmx_molprop_t *) b;
    char *mma = gmx_molprop_get_molname(*ma);
    char *mmb = gmx_molprop_get_molname(*mb);
  
    if (mma && mmb)
        return strcasecmp(mma,mmb);
    else
        return 0;
}

static int comp_mp_formula(const void *a,const void *b)
{
    int r;
    gmx_molprop_t *ma = (gmx_molprop_t *)a;
    gmx_molprop_t *mb = (gmx_molprop_t *)b;
    char *fma = gmx_molprop_get_formula(*ma);
    char *fmb = gmx_molprop_get_formula(*mb);
  
    r = strcasecmp(fma,fmb);
  
    if (r == 0) 
        return comp_mp_molname(a,b);
    else 
        return r;
}

gmx_atomprop_t my_aps;

static int comp_mp_elem(const void *a,const void *b)
{
    int i,r;
    char *elem;
    gmx_molprop_t *ma = (gmx_molprop_t *)a;
    gmx_molprop_t *mb = (gmx_molprop_t *)b;
  
    r = gmx_molprop_count_composition_atoms(*ma,(char *)"bosque",(char *)"C")-
        gmx_molprop_count_composition_atoms(*mb,(char *)"bosque",(char *)"C");
    if (r != 0)
        return r;
  
    for(i=1; (i<=109); i++) {
        if (i != 6) {
            elem = gmx_atomprop_element(my_aps,i);
            r = gmx_molprop_count_composition_atoms(*ma,(char *)"bosque",elem)-
                gmx_molprop_count_composition_atoms(*mb,(char *)"bosque",elem);
            if (r != 0)
                return r;
        }
    }
    return comp_mp_molname(a,b);
}

gmx_molselect_t my_gms;

static int comp_mp_selection(const void *a,const void *b)
{
    int ia,ib;
    gmx_molprop_t *ma = (gmx_molprop_t *)a;
    gmx_molprop_t *mb = (gmx_molprop_t *)b;

    ia = gmx_molselect_index(my_gms,gmx_molprop_get_iupac(*ma));
    ib = gmx_molselect_index(my_gms,gmx_molprop_get_iupac(*mb));
    
    return (ia - ib);
}

void gmx_molprop_sort(int np,gmx_molprop_t mp[],int alg,gmx_atomprop_t apt,
                      gmx_molselect_t gms)
{
    int size = sizeof(mp[0]);
  
    /* gmx_molprops_write("debug.xml",np,mp); */
    switch(alg) {
    case empSORT_Molname:
        qsort(mp,np,size,comp_mp_molname);
        break;
    case empSORT_Formula:
        qsort(mp,np,size,comp_mp_formula);
        break;
    case empSORT_Composition:
        my_aps = apt;
        qsort(mp,np,size,comp_mp_elem);
        my_aps = NULL;
        break;
    case empSORT_Selection:
        if (NULL != gms) 
        {
            my_gms = gms;
            qsort(mp,np,size,comp_mp_selection);
        }
        else
            gmx_fatal(FARGS,"Need molecule selection to sort on");
        break;
    default:
        gmx_incons("Invalid algorithm for gmx_molprop_sort");
    }
}

static void add_qmc_conf(t_qmcount *qmc,char *conformation)
{
    int j;
    
    for(j=0; (j<qmc->nconf); j++) 
    {
        if (strcasecmp(qmc->conf[j],conformation) == 0)
            break;
    }
    if (j == qmc->nconf) 
    {
        qmc->nconf++;
        srenew(qmc->conf,qmc->nconf);
        qmc->conf[j] = strdup(conformation);
    }
}

static int get_qmc_count(t_qmcount *qmc,char *method,char *basis,char *type)
{
    int j;
    
    for(j=0; (j<qmc->n); j++) 
    {
        if ((strcasecmp(qmc->method[j],method) == 0) &&
            (strcasecmp(qmc->basis[j],basis) == 0) &&
            (strcasecmp(qmc->type[j],type) == 0))
            return qmc->count[j];
    }
    return 0;
}

static void add_qmc_calc(t_qmcount *qmc,char *method,char *basis,char *type)
{
    int j;
    
    for(j=0; (j<qmc->n); j++) 
    {
        if ((strcasecmp(qmc->method[j],method) == 0) &&
            (strcasecmp(qmc->basis[j],basis) == 0) &&
            (strcasecmp(qmc->type[j],type) == 0))
            break;
    }
    if (j == qmc->n) 
    {
        srenew(qmc->method,qmc->n+1);
        srenew(qmc->basis,qmc->n+1);
        srenew(qmc->type,qmc->n+1);
        srenew(qmc->count,qmc->n+1);
        srenew(qmc->lot,qmc->n+1);
        qmc->method[qmc->n] = strdup(method);
        qmc->basis[qmc->n]  = strdup(basis);
        qmc->type[qmc->n]   = strdup(type);
        snew(qmc->lot[qmc->n],strlen(method)+strlen(basis)+strlen(type)+16);
        sprintf(qmc->lot[qmc->n],"%s/%s/%s",method,basis,type);
        qmc->count[qmc->n]  = 1;
        qmc->n++;
    }
    else
    {
        qmc->count[j]++;
    }
}

t_qmcount *find_calculations(int np,gmx_molprop_t mp[],int emp,char *fc_str)
{
    char *program,*method,*basis,*reference,*conformation;
    char **qm,**ll;
    int  i,j,k,n,calcref,expref;
    double value,error,vec[3];
    tensor quadrupole;
    t_qmcount *qmc;
    char *types[] = { (char *)"elec", (char *)"ESP", (char *)"RESP", 
                      (char *)"empirical", (char *)"refractive index", 
                      (char *)"HF", (char *)"MP2", (char *)"Molecular Beam" };
    int ntypes;
    
    ntypes = (sizeof(types)/sizeof(types[0]));
    snew(qmc,1);
    for(i=0; (i<np); i++) 
    {
        while(gmx_molprop_get_experiment(mp[i],&reference,&conformation,&expref) == 1) 
        {
            add_qmc_conf(qmc,conformation);
        }
        while(gmx_molprop_get_calculation(mp[i],&program,&method,&basis,
                                          &reference,&conformation,&calcref) == 1) 
        {
            add_qmc_conf(qmc,conformation);
        }
    }
    if (NULL != fc_str) 
    {
        qm = split(':',fc_str);
        n = 0;
        while (NULL != qm[n])
        {
            ll = split('/',qm[n]);
            if ((NULL != ll[0]) && (NULL != ll[1]) && (NULL != ll[2])) 
            {
                add_qmc_calc(qmc,ll[0],ll[1],ll[2]);
            }
            n++;
        }
    }
    
    for(i=0; (i<np); i++) 
    {
        while(gmx_molprop_get_calculation(mp[i],&program,&method,&basis,
                                          &reference,&conformation,&calcref) == 1) 
        {
            for(k=0; (k<ntypes); k++) 
            {
                if ((NULL == fc_str) || (get_qmc_count(qmc,method,basis,types[k]) > 0))
                {
                    if (get_val(mp[i],calcref,types[k],emp,&value,&error,vec,quadrupole) == 1) 
                    {
                        add_qmc_calc(qmc,method,basis,types[k]);
                    }
                }
            }
            sfree(program);
            sfree(method);
            sfree(basis);
            sfree(reference);
            sfree(conformation);
        }
    }
    for(i=0; (i<qmc->n); i++)
    {
        /* Since we initialized these we have counted one extra */
        if (NULL != fc_str) 
        {
            qmc->count[i]--;
        }
        
    }
        
    return qmc;
}

static void gmx_send_polar(t_commrec *cr,int dest,char *type,char *unit,
                           double x,double y,double z,double aver,double error)
{
    gmx_send_int(cr,dest,1);
    gmx_send_str(cr,dest,type);
    gmx_send_str(cr,dest,unit);
    gmx_send_double(cr,dest,x);
    gmx_send_double(cr,dest,y);
    gmx_send_double(cr,dest,z);
    gmx_send_double(cr,dest,aver);
    gmx_send_double(cr,dest,error);
}

static int gmx_recv_polar(t_commrec *cr,int src,char **type,char **unit,
                          double *x,double *y,double *z,double *aver,double *error)
{
    int kk = gmx_recv_int(cr,src);
    if (kk == 1) 
    {
        *type  = gmx_recv_str(cr,src);
        *unit  = gmx_recv_str(cr,src);
        *x     = gmx_recv_double(cr,src);
        *y     = gmx_recv_double(cr,src);
        *z     = gmx_recv_double(cr,src);
        *aver  = gmx_recv_double(cr,src);
        *error = gmx_recv_double(cr,src);
    }
    return kk;
}

static void gmx_send_quad(t_commrec *cr,int dest,char *type,char *unit,
                          double xx,double yy,double zz,
                          double xy,double xz,double yz)
{
    gmx_send_int(cr,dest,1);
    gmx_send_str(cr,dest,type);
    gmx_send_str(cr,dest,unit);
    gmx_send_double(cr,dest,xx);
    gmx_send_double(cr,dest,yy);
    gmx_send_double(cr,dest,zz);
    gmx_send_double(cr,dest,xy);
    gmx_send_double(cr,dest,xz);
    gmx_send_double(cr,dest,yz);
}

static int gmx_recv_quad(t_commrec *cr,int src,char **type,char **unit,
                         double *xx,double *yy,double *zz,
                         double *xy,double *xz,double *yz)
{
    int kk = gmx_recv_int(cr,src);
    if (kk == 1) 
    {
        *type  = gmx_recv_str(cr,src);
        *unit  = gmx_recv_str(cr,src);
        *xx    = gmx_recv_double(cr,src);
        *yy    = gmx_recv_double(cr,src);
        *zz    = gmx_recv_double(cr,src);
        *xy    = gmx_recv_double(cr,src);
        *xz    = gmx_recv_double(cr,src);
        *yz    = gmx_recv_double(cr,src);
    }
    return kk;
}

static void gmx_send_energy(t_commrec *cr,int dest,char *type,char *unit,
                            double aver,double error)
{
    gmx_send_int(cr,dest,1);
    gmx_send_str(cr,dest,type);
    gmx_send_str(cr,dest,unit);
    gmx_send_double(cr,dest,aver);
    gmx_send_double(cr,dest,error);
}

static int gmx_recv_energy(t_commrec *cr,int src,char **type,char **unit,
                           double *aver,double *error)
{
    int kk = gmx_recv_int(cr,src);
    if (kk == 1) 
    {
        *type  = gmx_recv_str(cr,src);
        *unit  = gmx_recv_str(cr,src);
        *aver  = gmx_recv_double(cr,src);
        *error = gmx_recv_double(cr,src);
    }
    return kk;
}

void gmx_molprop_send(t_commrec *cr,int dest,gmx_molprop_t mp)
{
    char *ptr,*atomname,*reference,*conformation,*type,*unit;
    char *program,*method,*basisset,*xyz_unit,*V_unit;
    int  natom,espid,atomid,expref,atomref;
    int  ccc=0;
    int  ncalc;
    double x,y,z,V,xx,yy,zz,xy,xz,yz,aver,value,error,q;
    
    gmx_molprop_reset(mp);
#ifdef GMX_MPI
    /* Generic stuff */
    gmx_send_double(cr,dest,gmx_molprop_get_mass(mp));
    gmx_send_int(cr,dest,gmx_molprop_get_charge(mp));
    gmx_send_int(cr,dest,gmx_molprop_get_multiplicity(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_formula(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_molname(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_iupac(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_cas(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_cid(mp));
    gmx_send_str(cr,dest,gmx_molprop_get_inchi(mp));

    /* Composition */
    while ((ptr = gmx_molprop_get_composition(mp)) != NULL) 
    {
        gmx_send_str(cr,dest,ptr);
        while (gmx_molprop_get_composition_atom(mp,ptr,&atomname,&natom) == 1) {
            gmx_send_int(cr,dest,natom);
            gmx_send_str(cr,dest,atomname);
        }
        gmx_send_int(cr,dest,0);
    }
    gmx_send_str(cr,dest,NULL);

    /* Categories */
    while ((ptr = gmx_molprop_get_category(mp)) != NULL)
    {
        gmx_send_str(cr,dest,ptr);
    }
    gmx_send_str(cr,dest,NULL);

    /* Experiments */
    while (gmx_molprop_get_experiment(mp,&reference,&conformation,&expref) == 1)
    {
        gmx_send_int(cr,dest,1);
        gmx_send_str(cr,dest,reference);
        gmx_send_str(cr,dest,conformation);
        while (gmx_molprop_get_polar(mp,expref,&type,&unit,
                                     &x,&y,&z,&aver,&error) == 1) 
        {
            gmx_send_polar(cr,dest,type,unit,x,y,z,aver,error);
        }
        gmx_send_int(cr,dest,0);
        while (gmx_molprop_get_dipole(mp,expref,&type,&unit,
                                      &x,&y,&z,&aver,&error) == 1) 
        {
            gmx_send_polar(cr,dest,type,unit,x,y,z,aver,error);
        }
        gmx_send_int(cr,dest,0);
        while (gmx_molprop_get_energy(mp,expref,&type,&unit,&value,&error) == 1)
        {
            gmx_send_energy(cr,dest,type,unit,value,error);
        }
        gmx_send_int(cr,dest,0);
    }
    gmx_send_int(cr,dest,0);

    /* Calculations */
    ncalc = 0;
    while (gmx_molprop_get_calculation(mp,&program,&method,&basisset,
                                       &reference,&conformation,&expref) == 1)
    {
        gmx_send_int(cr,dest,1);
        gmx_send_str(cr,dest,program);
        gmx_send_str(cr,dest,method);
        gmx_send_str(cr,dest,basisset);
        gmx_send_str(cr,dest,reference);
        gmx_send_str(cr,dest,conformation);
        while (gmx_molprop_get_polar(mp,expref,&type,&unit,
                                     &x,&y,&z,&aver,&error) == 1) 
        {
            gmx_send_polar(cr,dest,type,unit,x,y,z,aver,error);
        }
        gmx_send_int(cr,dest,0);
        while (gmx_molprop_get_dipole(mp,expref,&type,&unit,
                                      &x,&y,&z,&aver,&error) == 1) 
        {
            gmx_send_polar(cr,dest,type,unit,x,y,z,aver,error);
        }
        gmx_send_int(cr,dest,0);
        while (gmx_molprop_get_energy(mp,expref,&type,&unit,&value,&error) == 1)
        {
            gmx_send_energy(cr,dest,type,unit,value,error);
        }
        gmx_send_int(cr,dest,0);
        while (gmx_molprop_get_quadrupole(mp,expref,&type,&unit,
                                          &xx,&yy,&zz,&xy,&xz,&yz) == 1)
        {
            gmx_send_quad(cr,dest,type,unit,xx,yy,zz,xy,xz,yz);
        }
        gmx_send_int(cr,dest,0);

        while (gmx_molprop_get_potential(mp,expref,&xyz_unit,&V_unit,&espid,&x,&y,&z,&V) == 1)
        {
            gmx_send_int(cr,dest,1);
            gmx_send_str(cr,dest,xyz_unit);
            gmx_send_str(cr,dest,V_unit);
            gmx_send_int(cr,dest,espid);
            gmx_send_double(cr,dest,x);
            gmx_send_double(cr,dest,y);
            gmx_send_double(cr,dest,z);
            gmx_send_double(cr,dest,V);
        }
        gmx_send_int(cr,dest,0);
        
        while (gmx_molprop_calc_get_atom(mp,expref,&atomname,&atomid,&atomref) == 1)
        {
            gmx_send_str(cr,dest,atomname);
            gmx_send_int(cr,dest,atomid);
            if (gmx_molprop_calc_get_atomcoords(mp,expref,atomref,&unit,
                                                &x,&y,&z) == 1)
            {
                gmx_send_str(cr,dest,unit);
                gmx_send_double(cr,dest,x);
                gmx_send_double(cr,dest,y);
                gmx_send_double(cr,dest,z);
                while (gmx_molprop_calc_get_atomcharge(mp,expref,atomref,
                                                       &type,&unit,&q) == 1)
                {
                    gmx_send_str(cr,dest,type);
                    gmx_send_str(cr,dest,unit);
                    gmx_send_double(cr,dest,q);
                }
                gmx_send_str(cr,dest,NULL);
            }
        }
        gmx_send_str(cr,dest,NULL);
        ncalc++;
    }
    if (NULL != debug) 
        fprintf(debug,"Sent %d calculations to %d for mol %s\n",ncalc,dest,
                gmx_molprop_get_molname(mp));
    gmx_send_int(cr,dest,0);

#endif
}

gmx_molprop_t gmx_molprop_receive(t_commrec *cr,int src)
{
    gmx_molprop_t mp = gmx_molprop_init();
    char   *ptr,*atomname,*type,*unit,*reference,*conformation,*comp;
    char   *program,*method,*basisset,*xyz_unit,*V_unit;
    int    nexp,kk,natom,expref,atomref,ccc=0,espid,atomid;
    int    ncalc;
    double x,y,z,V,xx,yy,zz,xy,xz,yz,aver,error,q;
    
#ifdef GMX_MPI
    /* Generic stuff */
    gmx_molprop_set_mass(mp,gmx_recv_double(cr,src));
    gmx_molprop_set_charge(mp,gmx_recv_int(cr,src));
    gmx_molprop_set_multiplicity(mp,gmx_recv_int(cr,src));
    gmx_molprop_set_formula(mp,gmx_recv_str(cr,src));
    gmx_molprop_set_molname(mp,gmx_recv_str(cr,src));
    gmx_molprop_set_iupac(mp,gmx_recv_str(cr,src));
    gmx_molprop_set_cas(mp,gmx_recv_str(cr,src));
    gmx_molprop_set_cid(mp,gmx_recv_str(cr,src));
    gmx_molprop_set_inchi(mp,gmx_recv_str(cr,src));

    /* Composition */
    while ((ptr = gmx_recv_str(cr,src)) != NULL) 
    {
        gmx_molprop_add_composition(mp,ptr);
        while ((natom = gmx_recv_int(cr,src)) > 0)
        {
            comp = gmx_recv_str(cr,src);
            gmx_molprop_add_composition_atom(mp,ptr,comp,natom);
            sfree(comp);
        }
        sfree(ptr);
    }

    /* Categories */
    while ((ptr = gmx_recv_str(cr,src)) != NULL)
    {
        gmx_molprop_add_category(mp,ptr);
        sfree(ptr);
    }

    /* Experiments */
    while (gmx_recv_int(cr,src) == 1)
    {
        reference    = gmx_recv_str(cr,src);
        conformation = gmx_recv_str(cr,src);
        gmx_molprop_add_experiment(mp,reference,conformation,&expref);
        while (gmx_recv_polar(cr,src,&type,&unit,&x,&y,&z,&aver,&error) == 1)
        {
            gmx_molprop_add_polar(mp,expref,type,unit,x,y,z,aver,error);
            sfree(type);
            sfree(unit);
        }
        while (gmx_recv_polar(cr,src,&type,&unit,&x,&y,&z,&aver,&error) == 1)
        {
            gmx_molprop_add_dipole(mp,expref,type,unit,x,y,z,aver,error);
            sfree(type);
            sfree(unit);
        }
        while (gmx_recv_energy(cr,src,&type,&unit,&aver,&error) == 1)
        {
            gmx_molprop_add_energy(mp,expref,type,unit,aver,error);
            sfree(type);
            sfree(unit);
        }
    }
    /* Calculations */
    ncalc = 0;
    while (gmx_recv_int(cr,src) == 1)
    {
        program      = gmx_recv_str(cr,src);
        method       = gmx_recv_str(cr,src);
        basisset     = gmx_recv_str(cr,src);
        reference    = gmx_recv_str(cr,src);
        conformation = gmx_recv_str(cr,src);
        gmx_molprop_add_calculation(mp,program,method,basisset,
                                    reference,conformation,&expref);
        while (gmx_recv_polar(cr,src,&type,&unit,&x,&y,&z,&aver,&error) == 1)
        {
            gmx_molprop_add_polar(mp,expref,type,unit,x,y,z,aver,error);
            sfree(type);
            sfree(unit);
        }
        while (gmx_recv_polar(cr,src,&type,&unit,&x,&y,&z,&aver,&error) == 1)
        {
            gmx_molprop_add_dipole(mp,expref,type,unit,x,y,z,aver,error);
            sfree(type);
            sfree(unit);
        }
        while (gmx_recv_energy(cr,src,&type,&unit,&aver,&error) == 1)
        {
            gmx_molprop_add_energy(mp,expref,type,unit,aver,error);
            sfree(type);
            sfree(unit);
        }
        while (gmx_recv_quad(cr,src,&type,&unit,&xx,&yy,&zz,&xy,&xz,&yz) == 1)
        {
            gmx_molprop_add_quadrupole(mp,expref,type,unit,xx,yy,zz,xy,xz,yz);
            sfree(type);
            sfree(unit);
        }

        while (gmx_recv_int(cr,src) == 1) 
        {
            xyz_unit = gmx_recv_str(cr,src);
            V_unit = gmx_recv_str(cr,src);
            espid = gmx_recv_int(cr,src);
            x = gmx_recv_double(cr,src);
            y = gmx_recv_double(cr,src);
            z = gmx_recv_double(cr,src);
            V = gmx_recv_double(cr,src);
            gmx_molprop_add_potential(mp,expref,xyz_unit,V_unit,espid,x,y,z,V);
            sfree(xyz_unit);
            sfree(V_unit);
        }
        while ((atomname = gmx_recv_str(cr,src)) != NULL)
        {
            atomid = gmx_recv_int(cr,src);
            gmx_molprop_calc_add_atom(mp,expref,atomname,atomid,&atomref);
            if (NULL != (unit = gmx_recv_str(cr,src)))
            {
                x = gmx_recv_double(cr,src);
                y = gmx_recv_double(cr,src);
                z = gmx_recv_double(cr,src);
                gmx_molprop_calc_set_atomcoords(mp,expref,atomref,unit,x,y,z);
                sfree(unit);
                while ((type = gmx_recv_str(cr,src)) != NULL)
                {
                    unit = gmx_recv_str(cr,src);
                    q    = gmx_recv_double(cr,src);
                    gmx_molprop_calc_set_atomcharge(mp,expref,atomref,type,unit,q);
                    sfree(type);
                    sfree(unit);
                }
                sfree(atomname);
            }
        }
        ncalc++;
    }
    if (NULL != debug) 
        fprintf(debug,"Reveived %d calculations from %d for mol %s\n",ncalc,src,
                gmx_molprop_get_molname(mp));

#endif
    return mp;
}

gmx_bool gmx_molprop_support_composition(gmx_molprop_t mp,char *composition)
{
    gmx_bool bComp  = FALSE;
    char *ptr;
    
    if (NULL == composition)
        return TRUE;
        
    while ((ptr = gmx_molprop_get_composition(mp)) != NULL) {
        if (strcmp(composition,ptr) == 0)
            bComp = TRUE;
    }
    if (debug && !bComp)
        fprintf(debug,"No %s composition for molecule %s\n",composition,
                gmx_molprop_get_molname(mp));
    return bComp;   
}
