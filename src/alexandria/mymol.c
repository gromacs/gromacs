/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: molprop_xml.c,v 1.23 2009/06/01 06:13:18 spoel Exp $
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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "toputil.h"
#include "gmx_fatal.h"
#include "force.h"
#include "vec.h"
#include "smalloc.h"
#include "convparm.h"
#include "mtop_util.h"
#include "shellfc.h"
#include "mdatoms.h"
#include "symtab.h"
#include "mymol.h"
#include "gpp_nextnb.h"
#include "statutil.h"
#include "string2.h"
#include "gpp_atomtype.h"
#include "gentop_vsite.h"
#include "molprop_util.h"
#include "gentop_core.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "gmx_simple_comm.h"

const char *immsg(int imm)  
{ 
    static const char *msg[immNR] = {
        "OK", "Zero Dipole", "No Quadrupole", "Charged", "Error", 
        "Atom type problem", "Atom number problem", "Converting from molprop",
        "Determining bond order",
        "QM Inconsistency (ESP dipole does not match Elec)", 
        "Not in training set", "No experimental data"
    };
    assert(imm >= 0);
    assert(imm < immNR);
  
    return msg[imm];
}

static void generate_nbparam(int ftype,int comb,double ci[],double cj[],
                             t_iparams *ip)
{
    double sig,eps;
    
    switch(ftype) {
    case F_LJ:
        switch (comb) {
        case eCOMB_GEOMETRIC:
            /* Gromos rules */
            ip->lj.c6  = sqrt(ci[0] * cj[0]);
            ip->lj.c12 = sqrt(ci[1] * cj[1]);
            break;
            
        case eCOMB_ARITHMETIC:
            /* c0 and c1 are epsilon and sigma */
            sig = (ci[0]+cj[0])*0.5;
            eps = sqrt(ci[1]*cj[1]);
            ip->lj.c6  = 4*eps*pow(sig,6);
            ip->lj.c12 = 4*eps*pow(sig,12);
            
            break;
        case eCOMB_GEOM_SIG_EPS:
            /* c0 and c1 are epsilon and sigma */
            sig = sqrt(ci[0]*cj[0]);
            eps = sqrt(ci[1]*cj[1]);
            ip->lj.c6  = 4*eps*pow(sig,6);
            ip->lj.c12 = 4*eps*pow(sig,12);
            
            break;
        default:
            gmx_fatal(FARGS,"No such combination rule %d",comb);
        }
        break;
    default:
        gmx_fatal(FARGS,"No such function type supported %s",
                  interaction_function[ftype].name);
    }
}

static void do_init_mtop(gmx_mtop_t *mtop,int nmoltype,char **molname,
                         int natoms,t_atoms **atoms,gmx_poldata_t pd)
{
    init_mtop(mtop);
    if (nmoltype <= 0)
    {
        gmx_incons("Number of moltypes less than 1 in do_init_mtop");
    } 
    mtop->name = molname;
    mtop->nmoltype = nmoltype;
    snew(mtop->moltype,mtop->nmoltype);
    mtop->moltype[0].name = molname;
    mtop->nmolblock = nmoltype;
    snew(mtop->molblock,mtop->nmolblock);
    mtop->molblock[0].nmol = 1;
    mtop->molblock[0].type = 0;
    mtop->molblock[0].natoms_mol = natoms;
    
    /* Create a charge group block */
    stupid_fill_block(&(mtop->moltype[0].cgs),natoms,FALSE);
    
    mtop->natoms = natoms;
    init_t_atoms(&(mtop->moltype[0].atoms),natoms,FALSE);
    *atoms = &(mtop->moltype[0].atoms);
}

static void mk_ff_mtop(gmx_mtop_t *mtop,gmx_poldata_t pd)
{
    int i,j,tp,ftype,nrfp,comb_rule,ntype;
    char **ptr;
    char *params;
    double **c;

    ntype = gmx_poldata_get_natypes(pd);
    
    if (ntype <= 0)
    {
        gmx_incons("Number of atomtypes less than 1 in mk_ff_mtop");
    } 
    mtop->ffparams.reppow=12;
    mtop->ffparams.atnr=ntype;
    mtop->ffparams.fudgeQQ=gmx_poldata_get_fudgeQQ(pd);
    mtop->ffparams.ntypes=ntype*ntype;
    snew(mtop->ffparams.functype,ntype*ntype);
    snew(mtop->ffparams.iparams,ntype*ntype);
    
    comb_rule = gmx_poldata_get_comb_rule(pd);
    ftype = gmx_poldata_get_vdw_ftype(pd);
    nrfp  = NRFPA(ftype);
    
    /* Derive table of Van der Waals force field parameters */
    snew(c,ntype);
    i = 0;
    while (1 == gmx_poldata_get_atype(pd,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&params)) {
        ptr = split(' ',params);
        snew(c[i],MAXFORCEPARAM);
        j=0;
        while ((j<MAXFORCEPARAM) && (NULL != ptr[j])) {
            c[i][j] = atof(ptr[j]);
            sfree(ptr[j]);
            j++;
        }
        sfree(ptr);
        sfree(params);
        i++;
        if (i > ntype)
            gmx_fatal(FARGS,"Number of atom types apparently larger than %d",ntype);
    }
    /* Generate nonbonded matrix */
    for(i=tp=0; (i<ntype); i++) {
        for(j=0; (j<ntype); j++) {
            mtop->ffparams.functype[tp] = ftype;
            generate_nbparam(ftype,comb_rule,c[i],c[j],&mtop->ffparams.iparams[tp]);
            tp++;
        }
    }
    for(i=0; (i<ntype); i++) 
        sfree(c[i]);
    sfree(c);
}

static void excls_to_blocka(int natom,t_excls excls[],t_blocka *blocka)
{
    int i,j,k,nra;
    
    if (blocka->nr < natom)
    {
        srenew(blocka->index,natom+1);
    }
    nra = 0;
    for(i=0; (i<natom); i++)
        nra += excls[i].nr;
    snew(blocka->a,nra+1);
    nra = 0;
    for(i=j=0; (i<natom); i++)
    {
        blocka->index[i] = nra;
        for(k=0; (k<excls[i].nr); k++)
            blocka->a[j++] = excls[i].e[k];
        nra += excls[i].nr;
    }
    blocka->index[natom] = nra;
    blocka->nr = natom;
    blocka->nra = nra;
}

static void plist_to_mtop(gmx_poldata_t pd,t_params plist[],gmx_mtop_t *mtop)
{
    int i,j,k,l,nra,n=0,nrfp,ati,atj,tp;
    real c[MAXFORCEPARAM];
    double fudgeLJ;
    
    /* Generate pairs */
    fudgeLJ=gmx_poldata_get_fudgeLJ(pd);

    for(i=0; (i<F_NRE); i++)
    {
        nra = NRAL(i);
        nrfp = NRFPA(i);
        snew(mtop->moltype[0].ilist[i].iatoms,plist[i].nr*(nra+1));
        srenew(mtop->ffparams.functype,mtop->ffparams.ntypes+plist[i].nr);
        srenew(mtop->ffparams.iparams,mtop->ffparams.ntypes+plist[i].nr);
        k = 0;
        for(j=0; (j<plist[i].nr); j++)
        {
            if (i == F_LJ14) {
                ati = mtop->moltype[0].atoms.atom[plist[i].param[j].a[0]].type;
                atj = mtop->moltype[0].atoms.atom[plist[i].param[j].a[1]].type;
                tp = ati*mtop->ffparams.atnr+atj;
                c[0] = mtop->ffparams.iparams[tp].lj.c6*fudgeLJ;
                c[1] = mtop->ffparams.iparams[tp].lj.c12*fudgeLJ;
            }
            else {
                for(l=0; (l<nrfp); l++) {
                    c[l] = plist[i].param[j].c[l];
                    if (NOTSET == c[l])
                        c[l] = 0;
                }
            }
            n = enter_params(&mtop->ffparams,i,c,0,12,n,TRUE);
            mtop->moltype[0].ilist[i].iatoms[k++] = n;
            for(l=0; (l<nra); l++)
                mtop->moltype[0].ilist[i].iatoms[k++] = plist[i].param[j].a[l];
        }
        mtop->moltype[0].ilist[i].nr = k;
    }
}

void mtop_update_cgs(gmx_mtop_t *mtop)
{
    int i,j;
    
    for(i=0; (i<mtop->nmoltype); i++) 
    {
        if (mtop->moltype[i].atoms.nr > mtop->moltype[i].cgs.nr)
        {
            mtop->moltype[i].cgs.nr = mtop->moltype[i].atoms.nr;
            mtop->moltype[i].cgs.nalloc_index = mtop->moltype[i].atoms.nr+1;
            srenew(mtop->moltype[i].cgs.index,mtop->moltype[i].cgs.nr+1);
            for(j=0; (j<=mtop->moltype[i].cgs.nr); j++)
                mtop->moltype[i].cgs.index[j] = j;
        }
    }
}
            
static gmx_bool is_symmetric(t_mymol *mymol,real toler)
{
    int  i,j,m;
    real mm,tm;
    rvec com,test;
    gmx_bool *bSymm,bSymmAll;
  
    clear_rvec(com);
    tm = 0;
    for(i=0; (i<mymol->atoms->nr); i++) {
        mm  = mymol->atoms->atom[i].m;
        tm += mm; 
        for(m=0; (m<DIM); m++) {
            com[m] += mm*mymol->x[i][m];
        }
    }
    if (tm > 0) {
        for(m=0; (m<DIM); m++) 
            com[m] /= tm;
    }
    for(i=0; (i<mymol->atoms->nr); i++) 
        rvec_dec(mymol->x[i],com);
    
    snew(bSymm,mymol->atoms->nr);
    for(i=0; (i<mymol->atoms->nr); i++) {
        bSymm[i] = (norm(mymol->x[i]) < toler);
        for(j=i+1; (j<mymol->atoms->nr) && !bSymm[i]; j++) {
            rvec_add(mymol->x[i],mymol->x[j],test);
            if (norm(test) < toler) {
                bSymm[i] = TRUE;
                bSymm[j] = TRUE;
            }
        }
    }
    bSymmAll = TRUE;
    for(i=0; (i<mymol->atoms->nr); i++) {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for(i=0; (i<mymol->atoms->nr); i++) 
        rvec_inc(mymol->x[i],com);
  
    return bSymmAll;
}

static int init_mymol(FILE *fp,t_mymol *mymol,gmx_molprop_t mp,
                      gmx_bool bQM,char *lot,gmx_bool bZero,
                      gmx_poldata_t pd,gmx_atomprop_t aps,
                      int  iModel,t_commrec *cr,int *nwarn,
                      gmx_bool bCharged,const output_env_t oenv,
                      real th_toler,real ph_toler,
                      real dip_toler,real hfac,gmx_bool bH14,
                      gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                      int nexcl,gmx_bool bESP,
                      real watoms,real rDecrZeta,gmx_bool bPol,gmx_bool bFitZeta)
{
    int      ftb,i,j,k,ia,m,version,generation,step,*nbonds,tatomnumber,imm=immOK;
    char     *mylot=NULL,*myref=NULL;
    char     **molnameptr;
    rvec     xmin,xmax;
    tensor   quadrupole;
    double   value,error,vec[3];
    gentop_vsite_t gvt;
    real     btol = 0.2,Ha;
    t_pbc    pbc;
    t_excls  *newexcls=NULL;
    t_params plist[F_NRE];
    
    init_plist(plist);
    mymol->qtotal  = gmx_molprop_get_charge(mp);
    mymol->mult    = gmx_molprop_get_multiplicity(mp);
    mymol->natom   = gmx_molprop_get_natom(mp);
    init_enerdata(1,0,&mymol->enerd);
    ftb = gmx_poldata_get_bond_ftype(pd);
        
    /* Inputrec parameters */
    mymol->ir.tabext      = 2; /* nm */
    mymol->ir.ePBC        = epbcNONE;
    mymol->ir.epsilon_r   = 1;
    mymol->ir.vdwtype     = evdwCUT;
    mymol->ir.coulombtype = eelCUT;
    
    if (mymol->natom <= 0)
        imm = immAtomTypes;
    else
    {
        if (bPol)
            mymol->nalloc = mymol->natom*2+2;
        else
            mymol->nalloc = mymol->natom;
        
        if ((mymol->qtotal != 0) && !bCharged)
            imm = immCharged;
    }
    if (immOK == imm) 
    {
        mymol->molname  = strdup(gmx_molprop_get_molname(mp));
        /* Read coordinates */
        open_symtab(&(mymol->symtab));
        molnameptr = put_symtab(&(mymol->symtab),mymol->molname);
        if (strcmp(mymol->molname,"benzene") == 0)
            printf("BOE\n");
        do_init_mtop(&mymol->mtop,1,molnameptr,
                     mymol->natom,&(mymol->atoms),pd);
        if (molprop_2_atoms(mp,aps,&(mymol->symtab),lot,mymol->atoms,
                            (const char *)"ESP",&(mymol->x)) == 0)
            imm = immMolpropConv;
    }
    if (immOK == imm)
    {
        clear_rvec(xmin);
        clear_rvec(xmax);
        clear_rvec(mymol->coq);
        tatomnumber = 0;
        snew(mymol->qESP,mymol->nalloc);
        for(i=0; (i<mymol->atoms->nr); i++) 
        {
            mymol->qESP[i] = mymol->atoms->atom[i].q;
            for(m=0; (m<DIM); m++)
            {
                if (mymol->x[i][m] < xmin[m])
                    xmin[m] = mymol->x[i][m];
                else if (mymol->x[i][m] > xmax[m])
                    xmax[m] = mymol->x[i][m];
                mymol->coq[m] += mymol->x[i][m]*mymol->atoms->atom[i].atomnumber;
            }
            tatomnumber += mymol->atoms->atom[i].atomnumber;
        }
        if (tatomnumber > 0)
        {
            for(m=0; (m<DIM); m++)
                mymol->coq[m] /= tatomnumber;
            for(i=0; (i<mymol->atoms->nr); i++) 
                rvec_dec(mymol->x[i],mymol->coq);
        }
        else
        {
            imm = immAtomNumber;
        }
    }
    if (immOK == imm) 
    {
        if (!bZero && is_symmetric(mymol,0.01)) 
        {
            imm = immZeroDip;
        }
    }
    if (immOK == imm)
    {   
        clear_mat(mymol->box);
        for(m=0; (m<DIM); m++)
            mymol->box[m][m] = 2*(xmax[m]-xmin[m]) + 1;
        
        snew(nbonds,mymol->nalloc);
        snew(mymol->smnames,mymol->nalloc);
        snew(mymol->bRing,mymol->nalloc);
        mymol->nbond = mk_bonds(pd,mymol->atoms,mymol->x,NULL,plist,nbonds,
                                mymol->bRing,bH14,bAllDihedrals,bRemoveDoubleDihedrals,
                                nexcl,&mymol->excls,
                                TRUE,mymol->box,aps,btol,TRUE);
        snew(mymol->bondorder,mymol->nbond);
        mk_ff_mtop(&mymol->mtop,pd);
        
        /* Setting the atom types: this depends on the bonding */
        set_pbc(&pbc,epbcNONE,mymol->box);
        gvt = gentop_vsite_init(egvtALL);
        if ((mymol->atype = set_atom_type(NULL,mymol->molname,&(mymol->symtab),mymol->atoms,
                                          &(plist[ftb]),nbonds,mymol->bRing,mymol->bondorder,
                                          mymol->smnames,pd,aps,
                                          mymol->x,&pbc,th_toler,ph_toler,gvt)) == NULL) 
        {
            imm = immAtomTypes;
        }
        sfree(nbonds);
    }
    if (immOK == imm) 
    {
        gentop_vsite_generate_special(gvt,FALSE,mymol->atoms,&mymol->x,plist,
                                      &(mymol->symtab),mymol->atype,&mymol->excls,pd);

        gentop_vsite_done(&gvt);
        close_symtab(&(mymol->symtab));

        mymol->eSupport = eSupportLocal;
        if (mp_get_prop_ref(mp,empDIPOLE,(bQM ? iqmQM : iqmBoth),
                            lot,NULL,(char *)"elec",
                            &value,&error,&myref,&mylot,
                            vec,quadrupole) == 0)
        {
            if (!bZero)
                imm = immZeroDip;
            if (NULL != myref)
                sfree(myref);
            if (NULL != mylot)
                sfree(mylot);
        }
        else
        {
            mymol->dip_exp  = value;
            mymol->dip_err  = error;
            mymol->lot      = mylot;
            mymol->ref      = myref;
            for(m=0; (m<DIM); m++)
            {
                mymol->mu_exp[m] = vec[m];
            }
            mymol->mu_exp2 = sqr(value);
            if (error <= 0) {
                if (debug)
                    fprintf(debug,"WARNING: Error for %s is %g, assuming it is 10%%.\n",
                            gmx_molprop_get_molname(mp),error);
                (*nwarn)++;
                error = 0.1*value;
            }
            mymol->dip_weight = sqr(1.0/error);
        }
        if (mp_get_prop_ref(mp,empDIPOLE,iqmQM,
                            lot,NULL,(char *)"ESP",&value,&error,NULL,NULL,vec,quadrupole) != 0)
        {
            for(m=0; (m<DIM); m++)
            {
                mymol->mu_esp[m] = vec[m];
            }
        }
        if (mp_get_prop(mp,empENERGY,(bQM ? iqmQM : iqmBoth),
                        lot,NULL,"Hform",&value) != 0)
        {
            mymol->Hform = value;
            mymol->Emol = value;
            for(ia=0; (ia<mymol->atoms->nr); ia++) {
                if (gmx_atomprop_query(aps,epropHatomization,NULL,
                                       *mymol->atoms->atomname[ia],&Ha)) {
                    mymol->Emol -= Ha;
                }
                else {
                    mymol->Emol = 0;
                    break;
                }
            }
        }
        else {
            imm = immNoData;
        }
    }
    if (immOK == imm)
    {
        rvec dx;
        rvec_sub(mymol->mu_esp,mymol->mu_exp,dx);
        if (norm(dx) > dip_toler)
        {
            imm = immQMInconsistency;
        }
    }
    if (immOK == imm)
    {
        if (bQM)
        {
            if (mp_get_prop_ref(mp,empQUADRUPOLE,iqmQM,
                                lot,NULL,(char *)"elec",&value,&error,
                                NULL,NULL,vec,quadrupole) == 0)
                imm = immNoQuad;
            else
                copy_mat(quadrupole,mymol->Q_exp);
            if (mp_get_prop_ref(mp,empQUADRUPOLE,iqmQM,
                                lot,NULL,(char *)"ESP",&value,&error,
                                NULL,NULL,vec,quadrupole) != 0)
                copy_mat(quadrupole,mymol->Q_esp);
        }
    }
    if ((immOK == imm) && bESP)
    {
        char *xyz_unit,*V_unit;
        double x,y,z,V;
        int  cref,xu,vu,espnr,ftb;
        
        ftb = gmx_poldata_get_bond_ftype(pd);
        snew(mymol->symmetric_charges,mymol->nalloc);
        mymol->symmetric_charges = symmetrize_charges(TRUE,mymol->atoms,
                                                      &(plist[ftb]),pd,aps,(char *)"");
        
        mymol->gr = gmx_resp_init(pd,iModel,FALSE,0,0,mymol->qtotal,
                                  0,0,-1,TRUE,watoms,rDecrZeta,FALSE,1,
                                  bFitZeta,FALSE,NULL);
        if (NULL != mymol->gr)
        {
            gmx_resp_add_atom_info(mymol->gr,mymol->atoms,pd);
            gmx_resp_update_atomtypes(mymol->gr,mymol->atoms);
            gmx_resp_add_atom_symmetry(mymol->gr,pd,
                                       mymol->symmetric_charges);
            gmx_resp_add_atom_coords(mymol->gr,mymol->x);
            if (0 != (cref = gmx_molprop_get_calc_lot(mp,lot))) 
            {
                while (gmx_molprop_get_potential(mp,cref,&xyz_unit,
                                                 &V_unit,&espnr,&x,&y,&z,&V) == 1)
                {
                    /* Maybe not convert to gmx ? */
                    xu = string2unit(xyz_unit);
                    vu = string2unit(V_unit);
                    if (-1 == xu)
                        xu = eg2cAngstrom;
                    if (-1 == vu)
                        vu = eg2cHartree_e;
                    gmx_resp_add_point(mymol->gr,
                                       convert2gmx(x,xu),
                                       convert2gmx(y,xu),
                                       convert2gmx(z,xu),
                                       convert2gmx(V,vu));
                    sfree(xyz_unit);
                    sfree(V_unit);
                }
            }
        }
        else
            imm = immError;                          
    }
    else 
    {
        mymol->gr = NULL;
    }
    if (immOK == imm)
    {
        if (bPol) {
            mymol->vs = init_vsite(&mymol->mtop,cr);
            snew(mymol->buf,mymol->nalloc);
            snew(newexcls,mymol->nalloc);
            srenew(mymol->x,mymol->nalloc);
            add_shells(pd,mymol->nalloc,&(mymol->mtop.moltype[0].atoms),
                       mymol->atype,plist,mymol->x,&mymol->symtab,
                       &newexcls,mymol->smnames);
            mymol->mtop.natoms = mymol->mtop.moltype[0].atoms.nr;
            mymol->mtop.molblock[0].natoms_mol = mymol->mtop.natoms;
            excls_to_blocka(mymol->nalloc,newexcls,
                            &(mymol->mtop.moltype[0].excls));
            mtop_update_cgs(&mymol->mtop);
            reset_q(mymol->atoms);
            mymol->shell = init_shell_flexcon(debug,&mymol->mtop,0,mymol->x);
        }
        else {
            mymol->shell = NULL;
            excls_to_blocka(mymol->nalloc,mymol->excls,&(mymol->mtop.moltype[0].excls));
        }
        plist_to_mtop(pd,plist,&mymol->mtop);
        snew(mymol->f,mymol->nalloc);
        mymol->fr = mk_forcerec();
        init_forcerec(NULL,oenv,mymol->fr,NULL,&mymol->ir,&mymol->mtop,cr,
                      mymol->box,FALSE,NULL,NULL,NULL,NULL,NULL, TRUE,-1);
        init_state(&mymol->state,mymol->atoms->nr,1,1,1,0);
        mymol->ltop = gmx_mtop_generate_local_top(&(mymol->mtop),&(mymol->ir));
        mymol->md = init_mdatoms(NULL,&mymol->mtop,FALSE);
    }
    if (immOK != imm) 
    {
        mymol->qgen = NULL; 
    }
    if (immOK != imm) 
    {
        /* Remove temporary data */
        if (NULL != fp)
            fprintf(fp,"Could not add %s because of %s\n",mymol->molname,
                    immsg(imm));
    }
    else if (NULL != debug)
    {
        fprintf(debug,"Succesfully added %s\n",mymol->molname);
    }
    return imm;
}

static void add_index_count(t_index_count *ic,char *name,gmx_bool bConst)
{
    int i;
  
    for(i=0; (i<ic->n); i++) 
        if (strcasecmp(ic->name[i],name) == 0) {
            if (ic->bConst[i] != bConst)
                gmx_fatal(FARGS,"Trying to add atom %s as both constant and optimized",
                          name);
            else
                fprintf(stderr,"Trying to add %s twice\n",name);
        }
    if (i == ic->n) {
        ic->n++;
        srenew(ic->name,ic->n);
        srenew(ic->tot_count,ic->n);
        srenew(ic->count,ic->n);
        srenew(ic->bConst,ic->n);
        ic->name[i] = strdup(name);
        ic->tot_count[i] = 0;
        ic->count[i] = 0;
        ic->bConst[i] = bConst;
        if (bConst)
            ic->nconst++;
        else
            ic->nopt++;
    }
}

static void sum_index_count(t_index_count *ic,t_commrec *cr)
{
    int i;
    
    for(i=0; (i<ic->n); i++)
        ic->tot_count[i] = ic->count[i];
    if (cr->nnodes > 1)
        gmx_sumi(ic->n,ic->tot_count,cr);
}

static void inc_index_count(t_index_count *ic,char *name)
{
    int i;
  
    for(i=0; (i<ic->n); i++) 
        if (strcasecmp(ic->name[i],name) == 0) {
            ic->count[i]++;
            break;
        }
    if (i == ic->n)
        gmx_fatal(FARGS,"No such atom %s",name);
}

static void dec_index_count(t_index_count *ic,char *name)
{
    int i;
  
    for(i=0; (i<ic->n); i++) 
        if (strcasecmp(ic->name[i],name) == 0) {
            if (ic->count[i] > 0) {
                ic->count[i]--;
                break;
            }
            else
                gmx_fatal(FARGS,"Trying to decrease number of atoms %s below zero",
                          name);
        }
    if (i == ic->n)
        gmx_fatal(FARGS,"No such atom %s",name);
}

static int n_index_count(t_index_count *ic,char *name) 
{
    int i;
  
    for(i=0; (i<ic->n); i++) {
        if (strcasecmp(ic->name[i],name) == 0)
            return i;
    }
    return -1;
}

static int c_index_count(t_index_count *ic,char *name) 
{
    int i;
  
    for(i=0; (i<ic->n); i++) {
        if (strcasecmp(ic->name[i],name) == 0)
            return ic->tot_count[i];
    }
    return 0;
}

char *opt_index_count(t_index_count *ic)
{
    for(; (ic->nopt_c < ic->n); ic->nopt_c++)
        if (!ic->bConst[ic->nopt_c]) {
            return ic->name[ic->nopt_c++];
        }
    ic->nopt_c = 0;
  
    return NULL;
}

static gmx_bool const_index_count(t_index_count *ic,char *name) 
{
    int i;
  
    for(i=0; (i<ic->n); i++) {
        if (strcasecmp(ic->name[i],name) == 0)
            return ic->bConst[i];
    }
    return FALSE;
}

static void dump_index_count(t_index_count *ic,FILE *fp,
                             int iModel,gmx_poldata_t pd,
                             gmx_bool bFitZeta)
{
    int i,j,nZeta,nZopt;
    double zz;
    if (fp) {
        fprintf(fp,"Atom index for this optimization.\n");
        fprintf(fp,"Name  Number  Action   #Zeta\n");
        for(i=0; (i<ic->n); i++) {
            nZeta = gmx_poldata_get_nzeta(pd,iModel,ic->name[i]);
            nZopt = 0;
            for(j=0; (j<nZeta); j++)
            {
                zz = gmx_poldata_get_zeta(pd,iModel,ic->name[i],j);
                if (zz > 0)
                    nZopt++;
            }
            if (ic->bConst[i])
                fprintf(fp,"%-4s  %6d  Constant\n",ic->name[i],ic->count[i]);
            else
                fprintf(fp,"%-4s  %6d  Optimized %4d%s\n",
                        ic->name[i],ic->count[i],nZopt,
                        bFitZeta ? " optimized" : " constant");
        }
        fprintf(fp,"\n");
        fflush(fp);
    }
}

static int clean_index_count(t_index_count *ic,int minimum_data,FILE *fp)
{
    int i,j,nremove=0;
  
    for(i=0; (i<ic->n); ) {
        if (!ic->bConst[i] && (ic->tot_count[i] < minimum_data)) {
            if (fp)
                fprintf(fp,"Not enough support in data set for optimizing %s\n",
                        ic->name[i]);
            sfree(ic->name[i]);
            for(j=i; (j<ic->n-1); j++) {
                ic->name[j]   = ic->name[j+1];
                ic->count[j]  = ic->count[j+1];
                ic->tot_count[j]  = ic->tot_count[j+1];
                ic->bConst[j] = ic->bConst[j+1];
            }
            nremove++;
            ic->n--;
            ic->nopt--;
        }
        else 
            i++;
    }
    return nremove;
}

static int check_data_sufficiency(FILE *fp,int nmol,t_mymol mol[],
                                  int minimum_data,gmx_poldata_t pd,
                                  t_index_count *ic,gmx_atomprop_t aps,
                                  int iModel,char *opt_elem,char *const_elem,
                                  t_commrec *cr,gmx_bool bPol,
                                  gmx_bool bFitZeta)
{
    int i,j,nremove,nsupported;
    gmx_mtop_atomloop_all_t aloop;
    t_atom *atom;
    char   **ptr,*myname;
    int    k,at_global,mymodel,resnr;
  
    /* Parse opt_elem list to test which elements to optimize */
    if (const_elem) {
        ptr = split(' ',const_elem);
        for(k=0; (ptr[k]); k++) {
            if (gmx_poldata_have_eem_support(pd,iModel,ptr[k],FALSE)) 
                add_index_count(ic,ptr[k],TRUE);
        }
    }
    if (opt_elem) {
        ptr = split(' ',opt_elem);
        for(k=0; (ptr[k]); k++) {
            if (gmx_poldata_have_eem_support(pd,iModel,ptr[k],TRUE)) 
                add_index_count(ic,ptr[k],FALSE);
        } 
    }
    else {
        while (gmx_poldata_get_eemprops(pd,&mymodel,&myname,NULL,NULL,NULL,NULL,NULL) != 0) {
            if ((mymodel == iModel) &&
                !const_index_count(ic,myname) &&
                gmx_poldata_have_eem_support(pd,iModel,myname,FALSE))
                add_index_count(ic,myname,FALSE);
        }
    }
    sum_index_count(ic,cr);
    dump_index_count(ic,debug,iModel,pd,bFitZeta);
    for(i=0; (i<nmol); i++) {
        if (mol[i].eSupport != eSupportNo) 
        {
            aloop = gmx_mtop_atomloop_all_init(&mol[i].mtop);
            k = 0;
            while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom) &&
                   (mol[i].eSupport != eSupportNo)) 
            {
                if ((atom->atomnumber > 0) || !bPol)
                {
                    if (n_index_count(ic,*(mol[i].atoms->atomtype[k])) == -1) 
                    {
                        if (debug)
                            fprintf(debug,"Removing %s because of lacking support for atom %s\n",
                                    mol[i].molname,*(mol[i].atoms->atomtype[k]));
                        mol[i].eSupport = eSupportNo;
                    }
                }
                k++;
            }
            if (mol[i].eSupport != eSupportNo) {
                gmx_assert(k,mol[i].atoms->nr);
                aloop = gmx_mtop_atomloop_all_init(&mol[i].mtop);
                k = 0;
                while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) 
                {
                    if ((atom->atomnumber > 0) || !bPol)
                        inc_index_count(ic,*(mol[i].atoms->atomtype[k]));
                    k++;
                }
                gmx_assert(k,mol[i].atoms->nr);
            }
        }
    }
    do {
        sum_index_count(ic,cr);
        dump_index_count(ic,debug,iModel,pd,bFitZeta);
        nremove = 0;
        for(i=0; (i<nmol); i++) {
            if (mol[i].eSupport != eSupportNo) {
                j = 0;
                aloop = gmx_mtop_atomloop_all_init(&mol[i].mtop);
                while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) {
                    if (c_index_count(ic,*(mol[i].atoms->atomtype[j])) < minimum_data) {
                        if (debug)
                            fprintf(debug,"Removing %s because of no support for name %s\n",
                                    mol[i].molname,*(mol[i].atoms->atomtype[j]));
                        break;
                    }
                    j++;
                }
                if (j < mol[i].mtop.natoms) {
                    aloop = gmx_mtop_atomloop_all_init(&mol[i].mtop);
                    k = 0;
                    while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) 
                        dec_index_count(ic,*(mol[i].atoms->atomtype[k++]));
                    mol[i].eSupport = eSupportNo;
                    nremove++;
                }
            }
        }
        if (cr->nnodes > 1) 
        {
            /* Sum nremove */
            gmx_sumi(1,&nremove,cr);
        }
    } while (nremove > 0);
    nremove = clean_index_count(ic,minimum_data,debug);
    sum_index_count(ic,cr);
    dump_index_count(ic,fp,iModel,pd,bFitZeta);
    
    nsupported = 0;
    for(i=0; (i<nmol); i++) 
    {
        if (mol[i].eSupport == eSupportLocal)
        {
            if (NULL != debug)
                fprintf(debug,"Supported molecule %s on CPU %d\n",
                        mol[i].molname,cr->nodeid);
            nsupported++;
        }
    }
    if (cr->nnodes > 1) 
    {
        gmx_sumi(1,&nsupported,cr);
    }
    if (fp)
    {
        fprintf(fp,"Removed %d atomtypes\n",nremove);
        fprintf(fp,"There are %d supported molecules left.\n\n",nsupported);
    }
    return nsupported;
}

t_moldip *init_moldip(t_commrec *cr,gmx_bool bQM,gmx_bool bGaussianBug,
                      int  iModel,real rDecrZeta,real epsr,
                      real J0_0,real Chi0_0,real w_0,
                      real J0_1,real Chi0_1,real w_1,
                      real fc_bound,real fc_mu,real fc_quad,real fc_charge,
                      real fc_esp,char *fixchi,
                      gmx_bool bOptHfac,real hfac,
                      gmx_bool bPol,gmx_bool bFitZeta)
{
    t_moldip *md;
    
    snew(md,1);
    md->cr       = cr;
    md->bQM      = bQM;
    md->bDone    = FALSE;
    md->bFinal   = FALSE;
    md->bGaussianBug = bGaussianBug;
    md->bFitZeta = bFitZeta;
    md->iModel   = iModel;
    md->decrzeta = rDecrZeta;
    md->epsr     = epsr;
    md->J0_0       = J0_0;
    md->Chi0_0     = Chi0_0;
    md->w_0        = w_0;
    md->J0_1       = J0_1;
    md->Chi0_1     = Chi0_1;
    md->w_1        = w_1;
    md->fc[ermsMU]     = fc_mu;
    md->fc[ermsBOUNDS] = fc_bound;
    md->fc[ermsQUAD]   = fc_quad;
    md->fc[ermsCHARGE] = fc_charge;
    md->fc[ermsESP]    = fc_esp;
    md->fc[ermsForce2] = 1e-4;
    md->fc[ermsEPOT]   = 1;
    md->fixchi     = strdup(fixchi);
    md->hfac       = hfac;	  
    md->hfac0      = hfac;	  
    md->bOptHfac   = bOptHfac;
    md->bPol       = bPol;
  
    return md;
}

void read_moldip(t_moldip *md,
                 FILE *fp,const char *fn,const char *pd_fn,
                 int minimum_data,
                 gmx_bool bZero,gmx_bool bWeighted,
                 char *opt_elem,char *const_elem,
                 char *lot,gmx_bool bCharged,
                 output_env_t oenv,gmx_molselect_t gms,
                 real th_toler,real ph_toler,real dip_toler,
                 gmx_bool bH14,gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                 real watoms,gmx_bool bCheckSupport)
{
    char     **strings,buf[STRLEN];
    int      i,j,n,kk,nstrings,nwarn=0,nzero=0,nmol_cpu;
    double   dip,dip_err,dx,dy,dz;
    rvec     mu;
    int      nmol,nexcl,imm,imm_count[immNR];
    gmx_molprop_t *mp=NULL;
    char     *molname;
#ifdef GMX_MPI
    MPI_Status status;
#endif
    
    for(imm = 0; (imm<immNR); imm++)
        imm_count[imm] = 0;
        
    /* Read the EEM parameters */
    md->atomprop   = gmx_atomprop_init();
    
    /* Force field data */
    if ((md->pd = gmx_poldata_read(pd_fn,md->atomprop)) == NULL)
        gmx_fatal(FARGS,"Can not read the force field information. File %s missing or incorrect.",pd_fn);
    
    if ((n = gmx_poldata_get_numprops(md->pd,md->iModel)) == 0)
        gmx_fatal(FARGS,"File %s does not contain the requested parameters for model %d",pd_fn,md->iModel);
    nexcl = gmx_poldata_get_nexcl(md->pd);
    
    if (NULL != fp)
    {  
        fprintf(fp,"There are %d atom types in the input file %s:\n---\n",
                n,pd_fn);
        fprintf(fp,"---\n\n");
    }
    
    /* Read other stuff */
    if (MASTER(md->cr))
    {
        /* Now read the molecules */
        mp = gmx_molprops_read(fn,&nmol);
        for(i=0; (i<nmol); i++) 
            gmx_molprop_check_consistency(mp[i]);
        nmol_cpu = nmol/md->cr->nnodes + 1;
    }
    else {
        nmol_cpu = 0;
    }
    if (PAR(md->cr))
        gmx_sumi(1,&nmol_cpu,md->cr);
    if (MASTER(md->cr))
        snew(md->mymol,nmol);
    else
        snew(md->mymol,nmol_cpu);
    
    if (MASTER(md->cr)) 
    {
        for(i=n=0; (i<nmol); i++)
        {
            if (imsTrain == gmx_molselect_status(gms,gmx_molprop_get_iupac(mp[i])))
            {
                int dest = (n % md->cr->nnodes);
                
                imm = init_mymol(fp,&(md->mymol[n]),mp[i],md->bQM,lot,bZero,
                                 md->pd,md->atomprop,
                                 md->iModel,md->cr,&nwarn,bCharged,oenv,
                                 th_toler,ph_toler,dip_toler,md->hfac,bH14,
                                 bAllDihedrals,bRemoveDoubleDihedrals,nexcl,
                                 (md->fc[ermsESP] > 0),watoms,md->decrzeta,
                                 md->bPol,md->bFitZeta);
                if (immOK == imm)
                {
                    if (dest > 0)
                    {
                        md->mymol[n].eSupport = eSupportRemote;
                        /* Send another molecule */
                        gmx_send_int(md->cr,dest,1);
                        gmx_molprop_send(md->cr,dest,mp[i]);
                        imm = gmx_recv_int(md->cr,dest);
                        if (imm != immOK) 
                            fprintf(stderr,"Molecule %s was not accepted on node %d - error %s\n",
                                    md->mymol[n].molname,dest,immsg(imm));
                    }
                    else
                        md->mymol[n].eSupport = eSupportLocal;
                    if (immOK == imm)
                        n++;
                }
                if ((immOK != imm) && (NULL != debug))
                {
                    fprintf(debug,"IMM: Dest: %d %s - %s\n",
                            dest,gmx_molprop_get_molname(mp[i]),immsg(imm));
                }
            }
            else
                imm = immTest;
            imm_count[imm]++;
            /* Dispose of the molecule */
            gmx_molprop_delete(mp[i]);
        }
        /* Send signal done with transferring molecules */
        for(i=1; (i<md->cr->nnodes); i++) 
        {
            gmx_send_int(md->cr,i,0);
        }
    }
    else 
    {
        n = 0;
        while (gmx_recv_int(md->cr,0) == 1) 
        {
            /* Receive another molecule */
            gmx_molprop_t mpnew = gmx_molprop_receive(md->cr,0);
            imm = init_mymol(fp,&(md->mymol[n]),mpnew,md->bQM,lot,bZero,
                             md->pd,md->atomprop,
                             md->iModel,md->cr,&nwarn,bCharged,oenv,
                             th_toler,ph_toler,dip_toler,md->hfac,
                             bH14,bAllDihedrals,bRemoveDoubleDihedrals,
                             nexcl,(md->fc[ermsESP] > 0),
                             watoms,md->decrzeta,md->bPol,md->bFitZeta);
            md->mymol[n].eSupport = eSupportLocal;
            imm_count[imm]++;
            if (immOK == imm)
            {
                n++;
            }
            gmx_send_int(md->cr,0,imm);
            /* Dipose of the molecules */
            gmx_molprop_delete(mpnew);
        }
    }
    md->nmol = n;
    
    if (fp)
    {
        fprintf(fp,"There were %d warnings because of zero error bars.\n",nwarn);
        fprintf(fp,"Made topologies for %d out of %d molecules.\n",n,
                (MASTER(md->cr)) ? nmol : nmol_cpu);
        for(i=0; (i<immNR); i++)
            if (imm_count[i] > 0)
                fprintf(fp,"%d molecules - %s.\n",imm_count[i],immsg(i));
        if (imm_count[immOK] != nmol)
        {
            fprintf(fp,"Check %s.debug for more information.\nYou may have to use the -debug 1 flag.\n\n",ShortProgram());
        }
    }
    snew(md->ic,1);
    if (bCheckSupport) 
    {
        md->nmol_support =
            check_data_sufficiency(MASTER(md->cr) ? fp : NULL,md->nmol,md->mymol,
                                   minimum_data,md->pd,md->ic,md->atomprop,
                                   md->iModel,opt_elem,const_elem,md->cr,
                                   md->bPol,md->bFitZeta);
        if (md->nmol_support == 0)
            gmx_fatal(FARGS,"No support for any molecule!");
    }
    else {
        md->nmol_support = md->nmol;
    }
}


