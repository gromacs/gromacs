/* 
 * $Id: mymol2.cpp,v 1.23 2009/06/01 06:13:18 spoel Exp $
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
#include "pdb2top.h"
#include "confio.h"
#include "toputil.h"
#include "gmx_fatal.h"
#include "force.h"
#include "vec.h"
#include "smalloc.h"
#include "convparm.h"
#include "pdbio.h"
#include "macros.h"
#include "mtop_util.h"
#include "shellfc.h"
#include "mdatoms.h"
#include "symtab.h"
#include "gpp_nextnb.h"
#include "statutil.h"
#include "string2.h"
#include "gen_ad.h"
#include "gpp_atomtype.h"
#include "topdirs.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "gmx_simple_comm.h"
#include "gentop_vsite.hpp"
#include "gentop_core.hpp"
#include "molprop.hpp"
#include "molprop_util.hpp"
#include "molprop_xml.hpp"
#include "gauss_io.hpp"
#include "mymol.hpp"

static const char *gentop_version = "gentop 0.98";

static void get_force_constants(gmx_poldata_t pd,t_params plist[],t_atoms *atoms)
{
    int j,n,ft,k;
    double xx,sx,bo;
    char *params,**ptr;
    
    ft = gmx_poldata_get_bond_ftype(pd);
    for(j=0; (j<plist[ft].nr); j++) {
        if (1 == gmx_poldata_search_bond(pd,
                                         *atoms->atomtype[plist[ft].param[j].a[0]],
                                         *atoms->atomtype[plist[ft].param[j].a[1]],
                                         &xx,&sx,NULL,&bo,&params)) {
            ptr = split(' ',params);
            n = 0;
            while ((n<MAXFORCEPARAM) && (NULL != ptr[n])) {
                plist[ft].param[j].c[n] = atof(ptr[n]);
                sfree(ptr[n]);
                n++;
            }
            sfree(ptr);
        }
    }
    ft = gmx_poldata_get_angle_ftype(pd);
    for(j=0; (j<plist[ft].nr); j++) {
        if (1 == gmx_poldata_search_angle(pd,
                                          *atoms->atomtype[plist[ft].param[j].a[0]],
                                          *atoms->atomtype[plist[ft].param[j].a[1]],
                                          *atoms->atomtype[plist[ft].param[j].a[2]],
                                          &xx,&sx,NULL,&params)) {
            ptr = split(' ',params);
            n = 0;
            while ((n<MAXFORCEPARAM) && (NULL != ptr[n])) {
                plist[ft].param[j].c[n] = atof(ptr[n]);
                sfree(ptr[n]);
                n++;
            }
            sfree(ptr);
        }
    }
    for(k=0; (k<egdNR); k++) {
        ft = gmx_poldata_get_dihedral_ftype(pd,k);
        for(j=0; (j<plist[ft].nr); j++) {
            if (1 == gmx_poldata_search_dihedral(pd,k,
                                                 *atoms->atomtype[plist[ft].param[j].a[0]],
                                                 *atoms->atomtype[plist[ft].param[j].a[1]],
                                                 *atoms->atomtype[plist[ft].param[j].a[2]],
                                                 *atoms->atomtype[plist[ft].param[j].a[3]],
                                                 &xx,&sx,NULL,&params)) {
                ptr = split(' ',params);
                n = 0;
                while ((n<MAXFORCEPARAM) && (NULL != ptr[n])) {
                    plist[ft].param[j].c[n] = atof(ptr[n]);
                    sfree(ptr[n]);
                    n++;
                }
                sfree(ptr);
            }
        }
    }
}

namespace alexandria {

const char *immsg(immStatus imm)  
{ 
    static const char *msg[immNR] = {
        "Unknown status",
        "OK", "Zero Dipole", "No Quadrupole", "Charged", 
        "Atom type problem", "Atom number problem", "Converting from molprop",
        "Determining bond order", "RESP Initialization",
        "Charge generation",
        "QM Inconsistency (ESP dipole does not match Elec)", 
        "Not in training set", "No experimental data"
    };
    
    return msg[imm];
}

static void mv_plist(t_params *dst,t_params *src)
{
    int i;
    
    if (dst->maxnr < src->nr) {
        srenew(dst->param,src->nr);
        dst->maxnr = src->nr;
    }
    for(i=0; (i<src->nr); i++) {
        cp_param(&dst->param[i],&src->param[i]);
    }
    dst->nr=src->nr;
    
    src->nr = 0;
}

static void mv_plists(gmx_poldata_t pd,t_params plist[],bool bForward)
{
    int ft;
    
    /* Now move over to the appropriate function types */
    if (NOTSET == (ft = gmx_poldata_get_bond_ftype(pd)))
        gmx_fatal(FARGS,"Bond function type not set in force field file");
    if (F_BONDS != ft) {
        if (bForward)
            mv_plist(&plist[ft],&plist[F_BONDS]);
        else
            mv_plist(&plist[F_BONDS],&plist[ft]);
    }
    if (NOTSET == (ft = gmx_poldata_get_angle_ftype(pd)))
        gmx_fatal(FARGS,"Angle function type not set in force field file");
    if (F_ANGLES != ft) {
        if (bForward) 
            mv_plist(&plist[ft],&plist[F_ANGLES]);
        else
            mv_plist(&plist[F_ANGLES],&plist[ft]);
    }
    if (NOTSET == (ft = gmx_poldata_get_dihedral_ftype(pd,egdPDIHS)))
        gmx_fatal(FARGS,"Dihedral function type not set in force field file");
    if (F_PDIHS != ft) {
        if (bForward)
            mv_plist(&plist[ft],&plist[F_PDIHS]);
        else
            mv_plist(&plist[F_PDIHS],&plist[ft]);
    }
    if (NOTSET == (ft = gmx_poldata_get_dihedral_ftype(pd,egdIDIHS)))
        gmx_fatal(FARGS,"Improper function type not set in force field file");
    if (F_IDIHS != ft) {
        if (bForward) 
            mv_plist(&plist[ft],&plist[F_IDIHS]);
        else
            mv_plist(&plist[F_IDIHS],&plist[ft]);
    }
}

static void detect_rings(t_params *bonds,int natom,gmx_bool bRing[])
{
    /* Check for 4,5,6,7,8 rings.
     */
    int j,k,l,m,n,o,p,q,a1,a2,a3,a4,a5,a6,a7,a8,a9;
    
    for(j=0; (j<natom); j++) 
        bRing[j] = FALSE;
        
    for(a1=0; (a1<natom); a1++) 
    {
        for(j=0; (j<bonds->nr); j++) 
        {
            a2 = NOTSET;
            if (bonds->param[j].a[0] == a1)
                a2 = bonds->param[j].a[1];
            else if (bonds->param[j].a[1] == a1)
                a2 = bonds->param[j].a[0];
            if (a2 != NOTSET) {
                for(k=0; (k<bonds->nr); k++) 
                {
                    a3 = NOTSET;
                    if (bonds->param[k].a[0] == a2)
                        a3 = bonds->param[k].a[1];
                    else if (bonds->param[k].a[1] == a2)
                        a3 = bonds->param[k].a[0];
                    if ((a3 != NOTSET) && (a3 != a1)) 
                    {
                        for(l=0; (l<bonds->nr); l++) 
                        {
                            a4 = NOTSET;
                            if (bonds->param[l].a[0] == a3)
                                a4 = bonds->param[l].a[1];
                            else if (bonds->param[l].a[1] == a3)
                                a4 = bonds->param[l].a[0];
                            if ((a4 != NOTSET) && (a4 != a2)) 
                            {
                                for(m=0; (m<bonds->nr); m++) 
                                {
                                    a5 = NOTSET;
                                    if (bonds->param[m].a[0] == a4)
                                        a5 = bonds->param[m].a[1];
                                    else if (bonds->param[m].a[1] == a4)
                                        a5 = bonds->param[m].a[0];
                                    if ((a5 != NOTSET) && (a5 != a3)) 
                                    {
                                        if (a5 == a1) 
                                        {
                                            /* 4-ring */
                                            bRing[a1] = bRing[a2] = bRing[a3] = bRing[a4] = TRUE;
                                        }
                                        else if (a3 != a1) 
                                        {
                                            for(n=0; (n<bonds->nr); n++) 
                                            {
                                                a6 = NOTSET;
                                                if (bonds->param[n].a[0] == a5)
                                                    a6 = bonds->param[n].a[1];
                                                else if (bonds->param[n].a[1] == a5)
                                                    a6 = bonds->param[n].a[0];
                                                if ((a6 != NOTSET) && (a6 != a4)) 
                                                {
                                                    if (a6 == a1) 
                                                    {
                                                        /* 5-ring */
                                                        bRing[a1] = bRing[a2] = bRing[a3] = bRing[a4] = bRing[a5] = TRUE;
                                                    }
                                                    else 
                                                    {
                                                        for(o=0; (o<bonds->nr); o++) 
                                                        {
                                                            a7 = NOTSET;
                                                            if (bonds->param[o].a[0] == a6)
                                                                a7 = bonds->param[o].a[1];
                                                            else if (bonds->param[o].a[1] == a6)
                                                                a7 = bonds->param[o].a[0];
                                                            if ((a7 != NOTSET) && (a7 != a5)) 
                                                            {
                                                                if (a7 == a1) 
                                                                {
                                                                    /* 6-ring */
                                                                    bRing[a1] = bRing[a2] = bRing[a3] = 
                                                                        bRing[a4] = bRing[a5] = bRing[a6] = TRUE;
                                                                }
                                                                else  
                                                                {
                                                                    for(p=0; (p<bonds->nr); p++) 
                                                                    {
                                                                        a8 = NOTSET;
                                                                        if (bonds->param[p].a[0] == a7)
                                                                            a8 = bonds->param[p].a[1];
                                                                        else if (bonds->param[p].a[1] == a7)
                                                                            a8 = bonds->param[p].a[0];
                                                                        if ((a8 != NOTSET) && (a8 != a6)) 
                                                                        {
                                                                            if (a8 == a1) 
                                                                            {
                                                                                /* 7-ring */
                                                                                bRing[a1] = bRing[a2] = bRing[a3] = 
                                                                                    bRing[a4] = bRing[a5] = bRing[a6] = bRing[a7] = TRUE;
                                                                            }
                                                                            else  
                                                                            {
                                                                                for(q=0; (q<bonds->nr); q++) 
                                                                                {
                                                                                    a9 = NOTSET;
                                                                                    if (bonds->param[q].a[0] == a8)
                                                                                        a9 = bonds->param[q].a[1];
                                                                                    else if (bonds->param[q].a[1] == a8)
                                                                                        a9 = bonds->param[q].a[0];
                                                                                    if (a9 == a1) 
                                                                                    {
                                                                                        /* 8-ring */
                                                                                        bRing[a1] = bRing[a2] = bRing[a3] = 
                                                                                            bRing[a4] = bRing[a5] = bRing[a6] = bRing[a7] = bRing[a8] = TRUE;
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

int MyMol::MakeBonds(gmx_poldata_t pd,
                     gmx_conect gc,
                     t_params plist[],
                     int nbond[],
                     gmx_bool bH14,
                     gmx_bool bAllDihedrals,
                     gmx_bool bRemoveDoubleDihedrals,
                     gmx_bool bPBC,matrix box,gmx_atomprop_t aps,real tol,
                     gmx_bool bMovePlists)
{
    t_param  b;
    int      i,j,nbonds;
    t_pbc    pbc;
    rvec     dx;
    real     dx1,dx11;
    gmx_bool bBond;
    char     *elem_i,*elem_j;
    char     *length_unit;
    t_nextnb nnb;
    t_restp rtp;

    for(i=0; (i<MAXATOMLIST); i++)
        b.a[i] = -1;
    for(i=0; (i<MAXFORCEPARAM); i++)
        b.c[i] = 0.0;
    
    length_unit = gmx_poldata_get_length_unit(pd);
    if (bPBC)
        set_pbc(&pbc,-1,box);
    for(i=0; (i<topology->atoms.nr); i++) 
    {
        elem_i = gmx_atomprop_element(aps,topology->atoms.atom[i].atomnumber);
        for(j=i+1; (j<topology->atoms.nr); j++) 
        {
            elem_j = gmx_atomprop_element(aps,topology->atoms.atom[j].atomnumber);
            if (NULL != debug)
                fprintf(debug,"Testing %s-%d vs. %s-%d\n",elem_i,i,elem_j,j); 
            if (bPBC)
                pbc_dx(&pbc,x[i],x[j],dx);
            else
                rvec_sub(x[i],x[j],dx);
      
            dx11 = norm(dx);
            dx1 = gmx2convert(dx11,string2unit(length_unit));
      
            if (gc) 
            {
                bBond = gmx_conect_exist(gc,i,j);
            }
            else 
            {
                bBond = gmx_poldata_elem_is_bond(pd,elem_i,elem_j,dx1,tol);
            }
            if (bBond) 
            {
                b.a[0] = i;
                b.a[1] = j;
                b.c[0] = dx11;
                add_param_to_list (&(plist[F_BONDS]),&b);
                nbond[i]++;
                nbond[j]++;
            }
            if (debug) 
            {
                if (bBond)
                    fprintf(debug,"Bonding atoms %s-%d and %s-%d at distance %g\n",
                            *topology->atoms.atomname[i],i+1,*topology->atoms.atomname[j],j+1,dx1);
                else 
                    fprintf(debug,"Not bonding  %s-%d and %s-%d at distance %g\n",
                            *topology->atoms.atomname[i],i+1,*topology->atoms.atomname[j],j+1,dx1);
            }
        }
    }
    /* Make Angles and Dihedrals */
    snew(excls,topology->atoms.nr);
    init_nnb(&nnb,topology->atoms.nr,nexcl_+2);
    gen_nnb(&nnb,plist);
    detect_rings(&plist[F_BONDS],topology->atoms.nr,bRing);
    nbonds = plist[F_BONDS].nr;
    print_nnb(&nnb,"NNB");
    rtp.bKeepAllGeneratedDihedrals = TRUE;
    rtp.bRemoveDihedralIfWithImproper = TRUE;
    rtp.bGenerateHH14Interactions = TRUE;
    rtp.nrexcl = nexcl_;
    gen_pad(&nnb,&(topology->atoms),&rtp,plist,excls,NULL,FALSE);
    generate_excls(&nnb,nexcl_,excls);
    done_nnb(&nnb);
    
    /* Move the plist to the correct function */
    if (bMovePlists)
        mv_plists(pd,plist,true);
        
    return nbonds;
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
                         int natoms,gmx_poldata_t pd)
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
}

static void mk_ff_mtop(gmx_mtop_t *mtop,gmx_poldata_t pd)
{
    int i,j,tp,ftype,comb_rule,ntype;
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
            
bool MyMol::IsSymmetric(real toler)
{
    int  i,j,m;
    real mm,tm;
    rvec com,test;
    gmx_bool *bSymm,bSymmAll;
  
    clear_rvec(com);
    tm = 0;
    for(i=0; (i<topology->atoms.nr); i++) 
    {
        mm  = topology->atoms.atom[i].m;
        tm += mm; 
        for(m=0; (m<DIM); m++) 
        {
            com[m] += mm*x[i][m];
        }
    }
    if (tm > 0) 
    {
        for(m=0; (m<DIM); m++) 
            com[m] /= tm;
    }
    for(i=0; (i<topology->atoms.nr); i++)
        rvec_dec(x[i],com);
    
    snew(bSymm,topology->atoms.nr);
    for(i=0; (i<topology->atoms.nr); i++) 
    {
        bSymm[i] = (norm(x[i]) < toler);
        for(j=i+1; (j<topology->atoms.nr) && !bSymm[i]; j++) 
        {
            rvec_add(x[i],x[j],test);
            if (norm(test) < toler) 
            {
                bSymm[i] = TRUE;
                bSymm[j] = TRUE;
            }
        }
    }
    bSymmAll = TRUE;
    for(i=0; (i<topology->atoms.nr); i++) {
        bSymmAll = bSymmAll && bSymm[i];
    }
    sfree(bSymm);
    for(i=0; (i<topology->atoms.nr); i++) 
        rvec_inc(x[i],com);
  
    return bSymmAll;
}

static void fill_inputrec(t_inputrec *ir)
{
    ir->tabext      = 2; /* nm */
    ir->ePBC        = epbcNONE;
    ir->epsilon_r   = 1;
    ir->vdwtype     = evdwCUT;
    ir->coulombtype = eelCUT;
}

MyMol::MyMol(): MolProp()
{
    bHaveShells_ = false;
    bHaveVSites_ = false;
    cgnr_ = NULL;
    symmetric_charges_ = NULL;
    gr = NULL;
    immAtoms = immOK;
    immTopology = immOK;
    immCharges = immOK;
    open_symtab(&symtab);
    atype = init_atomtype();
    for(int i = 0; (i<ebtsNR); i++)
    {
        bts[i] = NOTSET;
    }
    clear_mat(box);
}

MyMol::~MyMol()
{
    sfree(cgnr_);
    sfree(symmetric_charges_);
    done_symtab(&symtab);
    done_atomtype(atype);
}

immStatus MyMol::GenerateAtoms(gmx_atomprop_t ap,
                               gmx_poldata_t pd,
                               const char *lot,
                               const char *q_algorithm)
{
    std::string molnm;
    
    int       myunit;
    double    q,xx,yy,zz;
    int       natom;
    gmx_bool  bDone = FALSE;
    immStatus imm = immOK;
        
    CalculationIterator ci = GetLot(lot);
    if (ci < EndCalculation()) 
    {
        t_param nb;
        
        memset(&nb,0,sizeof(nb));
        natom = 0;
        init_t_atoms(&(topology->atoms),NAtom(),FALSE);
        snew(x,NAtom());
        snew(topology->atoms.atomtype,NAtom());
        snew(topology->atoms.atomtypeB,NAtom());
        
        for(CalcAtomIterator cai=ci->BeginAtom(); (cai<ci->EndAtom()); cai++)
        {
            myunit = string2unit((char *)cai->GetUnit().c_str());
            if (myunit == -1)
                gmx_fatal(FARGS,"Unknown unit '%s' for atom coords",
                          cai->GetUnit().c_str());
            cai->GetCoords(&xx,&yy,&zz);
            x[natom][XX] = convert2gmx(xx,myunit);
            x[natom][YY] = convert2gmx(yy,myunit);
            x[natom][ZZ] = convert2gmx(zz,myunit);
            
            q = 0;
            for(AtomicChargeIterator qi=cai->BeginQ(); (qi<cai->EndQ()); qi++) 
            {
                if (strcmp(qi->GetType().c_str(),q_algorithm) == 0)
                {
                    myunit = string2unit((char *)qi->GetUnit().c_str());
                    q = convert2gmx(qi->GetQ(),myunit);
                    break;
                }
            }
            
            t_atoms_set_resinfo(&(topology->atoms),natom,&symtab,molnm.c_str(),1,' ',1,' ');
            topology->atoms.atomname[natom] = put_symtab(&symtab,cai->GetName().c_str());
            topology->atoms.atom[natom].atomnumber = gmx_atomprop_atomnumber(ap,cai->GetName().c_str());

            real mass;            
            if (gmx_atomprop_query(ap,epropMass,"???",cai->GetName().c_str(),&mass))
            {
                topology->atoms.atom[natom].m  = mass;
                topology->atoms.atom[natom].mB = mass;
            }
            strcpy(topology->atoms.atom[natom].elem,gmx_atomprop_element(ap,topology->atoms.atom[natom].atomnumber));
            
            topology->atoms.atom[natom].q  = 
                topology->atoms.atom[natom].qB = topology->atoms.atom[natom].q;
            topology->atoms.atom[natom].resind = 0;
            // First set the atomtype
            topology->atoms.atomtype[natom] = 
                topology->atoms.atomtypeB[natom] = put_symtab(&symtab,cai->GetObtype().c_str());
            
            natom++;
        }
        /* Change their atomtype from the OpenBabel internal type to the
         * one specified in our force field file (gentop.dat).
         */
        //translate_atomtypes(&topology->atoms,&(symtab),
        //                  gmx_poldata_get_force_field(pd));
        for(int i=0; (i<natom); i++) 
        {
            topology->atoms.atom[i].type = 
                topology->atoms.atom[i].typeB = add_atomtype(atype,&symtab,
                                                             &(topology->atoms.atom[i]),
                                                             *topology->atoms.atomtype[i],
                                                             &nb,
                                                             0,0.0,0.0,0.0,
                                                             topology->atoms.atom[i].atomnumber,
                                                             0.0,0.0);
        }
        topology->atoms.nr = natom;
        topology->atoms.nres = 1;
        bDone = ((natom == NAtom()) && (natom > 0));
    }
    else
    {
        imm = immAtomTypes;
    }
    if (NULL != debug)
        fprintf(debug,"Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molnm.c_str(),lot,natom);
    
    return imm;
}

immStatus MyMol::GenerateTopology(gmx_atomprop_t ap,
                                  gmx_poldata_t pd,
                                  const char *lot,
                                  const char *q_algorithm,
                                  bool bPol,
                                  int nexcl)
{
    alexandria::BondIterator bi;
    immStatus imm = immOK;
    int natom,nalloc,ftb;
    t_param b;
    t_nextnb nnb;
    t_restp *rtp;
   
    /* Set bts for topology output */
    bts[ebtsBONDS]  = gmx_poldata_get_bond_ftype(pd);
    bts[ebtsANGLES] = gmx_poldata_get_angle_ftype(pd);
    bts[ebtsIDIHS]  = gmx_poldata_get_dihedral_ftype(pd,egdIDIHS);
    bts[ebtsPDIHS]  = gmx_poldata_get_dihedral_ftype(pd,egdPDIHS);
 
    ftb = bts[ebtsBONDS];
    nexcl_ = nexcl;
    GenerateComposition(pd);
    if (NAtom() <= 0)
        imm = immAtomTypes;
    else
    {
        if (bPol)
            nalloc = NAtom()*2+2;
        else
            nalloc = NAtom();
        
        //if ((GetCharge() != 0) && !bCharged)
        //   imm = immCharged;
    }
    snew(topology,1);
    init_top(topology);
    
    /* Get atoms */
    imm = GenerateAtoms(ap,pd,lot,q_algorithm);

    if (immOK == imm)
    {    
        init_plist(plist);
        
        /* Store bonds in harmonic potential list first, update type later */
        ftb = F_BONDS;
        memset(&b,0,sizeof(b));
        for(bi=BeginBond(); (bi<EndBond()); bi++)
        {
            b.a[0] = bi->GetAi() - 1;
            b.a[1] = bi->GetAj() - 1;
            add_param_to_list(&(plist[ftb]),&b);
        }
        natom = NAtom();
        
        /* Make Angles and Dihedrals */
        snew(excls,natom);
        init_nnb(&nnb,natom,nexcl_+2);
        gen_nnb(&nnb,plist);
        //detect_rings(&plist[F_BONDS],topology->atoms.nr,bRing);
        //nbonds = plist[F_BONDS].nr;
        print_nnb(&nnb,"NNB");
        snew(rtp,1);
        rtp->bKeepAllGeneratedDihedrals = TRUE;
        rtp->bRemoveDihedralIfWithImproper = TRUE;
        rtp->bGenerateHH14Interactions = TRUE;
        rtp->nrexcl = nexcl_;
        generate_excls(&nnb,nexcl_,excls);
        gen_pad(&nnb,&(topology->atoms),rtp,plist,excls,NULL,FALSE);
        sfree(rtp);
        done_nnb(&nnb);
    }
    
    
    if (immOK == imm)
    {    
        get_force_constants(pd,plist,&topology->atoms);
    }
    
    return imm;
}

void MyMol::CalcMultipoles()
{
    int i,m;
    rvec mu,mm;
    real r2,dfac,q;
    gmx_mtop_atomloop_all_t aloop;
    t_atom *atom;     
    int    at_global;
    
    clear_rvec(mu);
    aloop = gmx_mtop_atomloop_all_init(&mtop);
    i = 0;
    clear_mat(Q_calc);
    clear_rvec(coq);
    while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) {
        q = atom->q;
        svmul(ENM2DEBYE*q,x[i],mm);
        rvec_inc(mu,mm);
        
        dfac = q*0.5*10*ENM2DEBYE;
        r2   = iprod(x[i],x[i]);
        for(m=0; (m<DIM); m++)
            Q_calc[m][m] += dfac*(3*sqr(x[i][m]) - r2);
        Q_calc[XX][YY] += dfac*3*(x[i][XX]+coq[XX])*(x[i][YY]+coq[YY]);
        Q_calc[XX][ZZ] += dfac*3*(x[i][XX]+coq[XX])*(x[i][ZZ]+coq[ZZ]);
        Q_calc[YY][ZZ] += dfac*3*(x[i][YY]+coq[YY])*(x[i][ZZ]+coq[ZZ]);
        
        i++;
    }
    gmx_assert(i,topology->atoms.nr);
    copy_rvec(mu,mu_calc);
    dip_calc = norm(mu);
}

immStatus MyMol::GenerateCharges(gmx_poldata_t pd,
                                 gmx_atomprop_t ap,
                                 int iModel,
                                 real hfac,real epsr,
                                 real qtol,int maxiter,int maxcycle,
                                 const char *lot,
                                 bool bSymmetricCharges,
                                 const char *symm_string)
{
    int i,ePBC,eQGEN;
    char qgen_msg[STRLEN];
    immStatus imm = immOK;
    
    if (immOK == imm)
    {
        if (bSymmetricCharges)
        {
            symmetric_charges_ = 
                symmetrize_charges(bSymmetricCharges,
                                   &topology->atoms,&(plist[F_BONDS]),
                                   pd,ap,symm_string);
        }
    }

    if (immOK == imm)
    {    
        //gmx_resp_get_atom_info(gr,&topology->atoms,&symtab,&x);
        gmx_resp_add_atom_info(gr,&topology->atoms,pd);
        gmx_resp_add_atom_symmetry(gr,symmetric_charges_);
        gmx_resp_update_atomtypes(gr,&(topology->atoms));
        gmx_resp_summary(stdout,gr,symmetric_charges_);
        gmx_resp_add_atom_coords(gr,x);
        
        CalculationIterator ci = GetLot(lot);
        if (ci != EndCalculation())
        {
            printf("There are %d potential points\n",ci->NPotential());
            for(ElectrostaticPotentialIterator epi=ci->BeginPotential(); (epi<ci->EndPotential()); epi++)
            {
                /* Maybe not convert to gmx ? */
                int xu = string2unit(epi->GetXYZunit().c_str());
                int vu = string2unit(epi->GetVunit().c_str());
                if (-1 == xu)
                    xu = eg2cAngstrom;
                if (-1 == vu)
                    vu = eg2cHartree_e;
                gmx_resp_add_point(gr,
                                   convert2gmx(epi->GetX(),xu),
                                   convert2gmx(epi->GetY(),xu),
                                   convert2gmx(epi->GetZ(),xu),
                                   convert2gmx(epi->GetV(),vu));
            }
        }

        /* Check which algorithm to use for charge generation */
        strcpy(qgen_msg,"");
        if (iModel == eqgNone)
        {
            printf("Using charges from gentop.dat\n");
            for(i=0; (i<topology->atoms.nr); i++)
            {
                char *qq = gmx_poldata_get_charge(pd,smnames[i]);
                double q;
                if (NULL != qq)
                    sscanf(qq,"%lf",&q);
                else
                    q = 0;
                topology->atoms.atom[i].q  = topology->atoms.atom[i].qB = q;
            }
            eQGEN = eQGEN_OK;
        }
        else 
        {
            qgen = gentop_qgen_init(pd,&topology->atoms,ap,x,iModel,hfac,GetCharge(),epsr);
            
            if (NULL == qgen)
                gmx_fatal(FARGS,"Can not generate charges for %s. Probably due to issues with atomtype detection or support.\n",GetMolname().c_str());
            eQGEN = generate_charges(NULL,
                                     qgen,gr,GetMolname().c_str(),pd,&topology->atoms,x,qtol,
                                     maxiter,maxcycle,ap,hfac);
            qgen_message(qgen,sizeof(qgen_msg),qgen_msg,gr);
            if (eQGEN_OK != eQGEN)
                imm = immChargeGeneration;
        }
    }
    
    return imm;
}

void MyMol::PrintConformation(const char *fn)
{
    char title[STRLEN];
    
    sprintf(title,"%s processed by %s",GetMolname().c_str(),ShortProgram());
    write_sto_conf(fn,title,&topology->atoms,x,NULL,epbcNONE,box);
}

static void write_zeta_q(FILE *fp,gentop_qgen_t qgen,
                         t_atoms *atoms,gmx_poldata_t pd,int iModel)
{
    int    i,j,k,nz,row;
    double zeta,q,qtot;
    gmx_bool   bAtom;
    
    if (NULL == qgen)
        return;
        
    fprintf(fp,"[ zeta_q ]\n");
    fprintf(fp,"; i     type    nq  row       zeta          q\n");
    k = -1;
    for(i=0; (i<atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
            k++;
        if (k == -1)
            gmx_fatal(FARGS,"The first atom must be a real atom, not a shell");
        nz = gentop_qgen_get_nzeta(qgen,k);
        if (nz != NOTSET)
        {
            fprintf(fp,"%5d %6s %5d",i+1,get_eemtype_name(iModel),(bAtom) ? nz-1 : 1);
            qtot = 0;
            for(j=(bAtom ? 0 : nz-1); (j<(bAtom ? nz-1 : nz)); j++)
            {
                row   = gentop_qgen_get_row(qgen,k,j);
                q     = gentop_qgen_get_q(qgen,k,j);
                zeta  = gentop_qgen_get_zeta(qgen,k,j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET)) 
                {
                    qtot += q;
                    fprintf(fp,"%5d %10g %10g",row,zeta,q);
                }
            }
            atoms->atom[i].q = qtot;
            fprintf(fp,"\n");
        }
    }
    fprintf(fp,"\n");
}

static void write_zeta_q2(gentop_qgen_t qgen,gpp_atomtype_t atype,
                          t_atoms *atoms,gmx_poldata_t pd,int iModel)
{
    FILE   *fp;
    int    i,j,k,nz,row;
    double zeta,q,qtot;
    gmx_bool   bAtom;
    
    if (NULL == qgen)
        return;

    fp = fopen("zeta_q.txt","w");        
    k = -1;
    for(i=0; (i<atoms->nr); i++)
    {
        bAtom = (atoms->atom[i].ptype == eptAtom);
        if (bAtom)
            k++;
        if (k == -1)
            gmx_fatal(FARGS,"The first atom must be a real atom, not a shell");
        nz = gentop_qgen_get_nzeta(qgen,k);
        if (nz != NOTSET)
        {
            fprintf(fp,"%6s  %5s  %5d",get_eemtype_name(iModel),
                    get_atomtype_name(atoms->atom[i].type,atype),
                    (bAtom) ? nz-1 : 1);
            qtot = 0;
            for(j=(bAtom ? 0 : nz-1); (j<(bAtom ? nz-1 : nz)); j++)
            {
                row   = gentop_qgen_get_row(qgen,k,j);
                q     = gentop_qgen_get_q(qgen,k,j);
                zeta  = gentop_qgen_get_zeta(qgen,k,j);
                if ((row != NOTSET) && (q != NOTSET) && (zeta != NOTSET)) 
                {
                    qtot += q;
                    fprintf(fp,"%5d %10g %10g",row,zeta,q);
                }
            }
            atoms->atom[i].q = qtot;
            fprintf(fp,"\n");
        }
    }
    fprintf(fp,"\n");
    fclose(fp);
}

void MyMol::PrintTopology(const char *fn,gmx_poldata_t pd,int iModel,
                          const char *forcefield,bool bVerbose)
{
    FILE   *fp;
    t_mols printmol;
    bool   bITP;
    int    i,bts2[ebtsNR];

    if (GetMolname().size() > 0)
    {
        printmol.name = strdup(GetMolname().c_str());
    }
    else if (GetFormula().size() > 0)
    {
        printmol.name = strdup(GetFormula().c_str());
    }
    else
    {
        printmol.name = strdup("Onbekend");
    }
    printmol.nr   = 1;
    
    /* Write topology file */
    bITP = (fn2ftp(fn) == efITP);
    fp = ffopen(fn,"w");
    if (!bITP)
        print_top_header(fp,fn,gentop_version,bITP,
                         GetForceField().c_str(),
                         1.0,"No Comment");
    
    /* Make pdb2gmx compatible bts array */
    for(i=0; (i<ebtsNR); i++)
        bts2[i] = NOTSET;
    for(i=1; (i<20) && (bts2[ebtsBONDS] == NOTSET); i++)
        if (ifunc_index(d_bonds,i) == bts[ebtsBONDS])
            bts2[ebtsBONDS] = i;
    for(i=1; (i<20) && (bts2[ebtsANGLES] == NOTSET); i++)
        if (ifunc_index(d_angles,i) == bts[ebtsANGLES])
            bts2[ebtsANGLES] = i;
    for(i=1; (i<20) && (bts2[ebtsPDIHS] == NOTSET); i++)
        if (ifunc_index(d_dihedrals,i) == bts[ebtsPDIHS])
            bts2[ebtsPDIHS] = i;
    for(i=1; (i<20) && (bts2[ebtsIDIHS] == NOTSET); i++)
        if (ifunc_index(d_dihedrals,i) == bts[ebtsIDIHS])
            bts2[ebtsIDIHS] = i;
    bts2[ebtsEXCLS] = 0;
    bts2[ebtsCMAP] = 0;
    for(i=0; (i<ebtsNR); i++)
        if (NOTSET == bts2[i])
            gmx_fatal(FARGS,"Could not find ftype for bts[%d]",i);
    
    
    if (bHaveShells_)
    {
        /* write_zeta_q(fp,qqgen,&topology->atoms,pd,iModel);*/
        write_zeta_q2(qgen,atype,&topology->atoms,pd,iModel);
    }
    write_top(fp,NULL,printmol.name,&topology->atoms,FALSE,bts2,plist,excls,atype,cgnr_,nexcl_);
    if (!bITP)
        print_top_mols(fp,printmol.name,forcefield,NULL,0,NULL,1,&printmol);
    
    if (bVerbose)
    {
        printf("There are %4d proper dihedrals, %4d impropers\n"
               "          %4d angles, %4d linear angles\n"
               "          %4d pairs, %4d bonds, %4d atoms\n"
               "          %4d polarizations\n",
               plist[bts[ebtsPDIHS]].nr,  plist[bts[ebtsIDIHS]].nr, 
               plist[bts[ebtsANGLES]].nr, plist[F_LINEAR_ANGLES].nr,
               plist[F_LJ14].nr,   plist[bts[ebtsBONDS]].nr,topology->atoms.nr,
               plist[F_POLARIZATION].nr);
    }

    fclose(fp);
}

void MyMol::PrintRTPEntry(const char *fn)
{
    print_rtp(fn,gentop_version,
              &topology->atoms,plist,cgnr_,asize(bts),bts);
}

void MyMol::GenerateVsitesShells(gmx_poldata_t pd,bool bGenVSites,bool bAddShells,
                                 bool bPairs,eDih edih)
{
    if (bGenVSites)
    {
        int i,anr = topology->atoms.nr;
        gentop_vsite_generate_special(gvt,bGenVSites,&topology->atoms,&x,plist,
                                      &symtab,atype,&excls,pd);
        if (topology->atoms.nr > anr) 
        {
            srenew(smnames,topology->atoms.nr);
            for(i=anr; (i<topology->atoms.nr); i++) 
            {
                smnames[i] = strdup("ML");
            }
            bHaveVSites_ = true;
        }
    }
    if (!bPairs)
        plist[F_LJ14].nr = 0;
    if (edih == edihNo)
        plist[bts[ebtsPDIHS]].nr = 0;
    
    if (bAddShells)
    {
        int nalloc = topology->atoms.nr*2+2;
        srenew(x,nalloc);
        srenew(smnames,nalloc);
        srenew(excls,nalloc);
        add_shells(pd,nalloc,&topology->atoms,atype,plist,x,&symtab,&excls,smnames);
        bHaveShells_ = true;
    }
    real mu = calc_dip(&topology->atoms,x);
}

void MyMol::GenerateChargeGroups(eChargeGroup ecg,bool bUsePDBcharge,
                                 const char *ndxfn,int nmol)
{
    real qtot,mtot;
    
    if ((cgnr_ = generate_charge_groups(ecg,&topology->atoms,
                                        &plist[bts[ebtsBONDS]],&plist[F_POLARIZATION],
                                        bUsePDBcharge,
                                        &qtot,&mtot)) == NULL)
        gmx_fatal(FARGS,"Error generating charge groups");
        
    if (ecg != ecgAtom)
        sort_on_charge_groups(cgnr_,&topology->atoms,plist,x,excls,smnames,
                              ndxfn,nmol);
}

void MyMol::GenerateCube(int iModel,
                         gmx_poldata_t pd,
                         real spacing,
                         const char *reffn,
                         const char *pcfn,
                         const char *pdbdifffn,
                         const char *potfn,
                         const char *rhofn,
                         const char *hisfn,
                         const char *difffn,
                         const char *diffhistfn,
                         output_env_t oenv)
{
    char       *gentop_version = (char *)"v0.99b";
    gmx_resp_t grref;
                        
    if (NULL != gr) 
    {
        /* This has to be done before the grid is f*cked up by 
           writing a cube file */
        grref = gmx_resp_copy(gr);
        gmx_resp_potcomp(gr,pcfn,pdbdifffn,oenv);
        if ((NULL != potfn) || (NULL != hisfn) || (NULL != rhofn) || 
            ((NULL != difffn) && (NULL != reffn)))
        {
            char buf[256];
            
            sprintf(buf,"Potential generated by %s based on %s charges",
                    gentop_version,
                    get_eemtype_name(iModel));
                    
            if (NULL != difffn)
            {
                gmx_resp_add_atom_info(grref,&topology->atoms,pd);
                gmx_resp_add_atom_symmetry(grref,symmetric_charges_);
                gmx_resp_read_cube(grref,reffn,FALSE);
                gmx_resp_copy_grid(gr,grref);
            }
            else 
            {
                gmx_resp_make_grid(gr,spacing,box,x);
            }
            if (NULL != rhofn)
            {
                sprintf(buf,"Electron density generated by %s based on %s charges",
                        gentop_version,get_eemtype_name(iModel));
                gmx_resp_calc_rho(gr);
                gmx_resp_write_rho(gr,rhofn,buf);
            }
            sprintf(buf,"Potential generated by %s based on %s charges",
                    gentop_version,get_eemtype_name(iModel));
            if (NULL != potfn)
            {
                gmx_resp_calc_pot(gr);
                gmx_resp_write_cube(gr,potfn,buf);
            }
            if (NULL != hisfn)
            {
                gmx_resp_write_histo(gr,hisfn,buf,oenv);
            }
            if ((NULL != difffn) || (NULL != diffhistfn))
            {
                sprintf(buf,"Potential difference generated by %s based on %s charges",
                        gentop_version,
                        get_eemtype_name(iModel));
                gmx_resp_write_diff_cube(grref,gr,difffn,diffhistfn,buf,oenv,0);
                gmx_resp_destroy(grref);
            }
        }
        gmx_resp_destroy(grref);
    }
}
   
immStatus MyMol::Initxx(FILE *fp,
                      alexandria::GaussAtomProp &gap,
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
    int      nalloc=0;
    int      ftb,i,ia,m,nbond,*nbonds,tatomnumber;
    immStatus imm=immOK;
    char     *mylot=NULL,*myref=NULL;
    char     **molnameptr;
    rvec     xmin,xmax;
    tensor   quadrupole;
    double   value,dv0,dv298,error,vec[3];
    gentop_vsite_t gvt;
    real     btol = 0.2;
    t_pbc    pbc;
    t_excls  *newexcls=NULL;
    t_param  b;
    // This is not good
    gmx_resp_t gr;
    
    ftb = gmx_poldata_get_bond_ftype(pd);
    if (immTopology != immOK)
    {
        printf("Sorry, no topology!\n");
    }
    init_enerdata(1,0,&enerd);
    /* Inputrec parameters */
    fill_inputrec(&ir);
    
    if (immOK == imm) 
    {
        /* Read coordinates */
        open_symtab(&symtab);
        molnameptr = put_symtab(&(symtab),GetMolname().c_str());
        do_init_mtop(&mtop,1,molnameptr,NAtom(),pd);
    }
    if (immOK == imm) 
    {
        /* Change their atomtype from the OpenBabel internal type to the
         * one specified in our force field file (gentop.dat).
         */
        //translate_atomtypes(&topology->atoms,&(symtab),
        //                  gmx_poldata_get_force_field(pd));
    }
    if (immOK == imm)
    {
        clear_rvec(xmin);
        clear_rvec(xmax);
        clear_rvec(coq);
        tatomnumber = 0;
        snew(qESP,nalloc);
        for(i=0; (i<topology->atoms.nr); i++) 
        {
            qESP[i] = topology->atoms.atom[i].q;
            for(m=0; (m<DIM); m++)
            {
                if (x[i][m] < xmin[m])
                    xmin[m] = x[i][m];
                else if (x[i][m] > xmax[m])
                    xmax[m] = x[i][m];
                coq[m] += x[i][m]*topology->atoms.atom[i].atomnumber;
            }
            tatomnumber += topology->atoms.atom[i].atomnumber;
        }
        if (tatomnumber > 0)
        {
            for(m=0; (m<DIM); m++)
                coq[m] /= tatomnumber;
            for(i=0; (i<topology->atoms.nr); i++) 
                rvec_dec(x[i],coq);
        }
        else
        {
            imm = immAtomNumber;
        }
    }
    if (immOK == imm) 
    {
        if (!bZero && IsSymmetric(0.01)) 
        {
            imm = immZeroDip;
        }
    }
    if (immOK == imm)
    {   
        clear_mat(box);
        for(m=0; (m<DIM); m++)
            box[m][m] = 2*(xmax[m]-xmin[m]) + 1;
        nbond = NBond();
        gvt = gentop_vsite_init(egvtALL);
        if (nbond > 0) 
        {
            alexandria::BondIterator bi;
            
            /* Use the bonds stored in molprop */
            i = 0;
            for(bi=BeginBond(); (bi<EndBond()); bi++)
            {
                b.a[0] = bi->GetAi();
                b.a[1] = bi->GetAj();
                add_param_to_list(&(plist[ftb]),&b),
                i++;
            }
            if (i != nbond)
                imm = immAtomTypes;
            else
            {
                /* Get exclusions and other stuff */
            }
        }
        else 
        {
            double *bondorder;
            snew(nbonds,nalloc);
            snew(smnames,nalloc);
            MakeBonds(pd,(gmx_conect)NULL,plist,nbonds,
                      bH14,bAllDihedrals,bRemoveDoubleDihedrals,
                      TRUE,box,aps,btol,TRUE);
            mk_ff_mtop(&mtop,pd);
        
            /* Setting the atom types: this depends on the bonding */
            set_pbc(&pbc,epbcNONE,box);
            snew(bondorder,nbond);
            
            if ((atype = set_atom_type(NULL,GetMolname().c_str(),&(symtab),&(topology->atoms),
                                       &(plist[ftb]),nbonds,bRing,bondorder,
                                       smnames,pd,aps,
                                       x,&pbc,th_toler,ph_toler,gvt)) == NULL) 
            {
                imm = immAtomTypes;
            }
            sfree(nbonds);
            sfree(bondorder);
            sfree(smnames);
        }
    }
    if (immOK == imm) 
    {
        gentop_vsite_generate_special(gvt,FALSE,&(topology->atoms),&x,plist,
                                      &(symtab),atype,&excls,pd);

        gentop_vsite_done(&gvt);
        close_symtab(&(symtab));

        eSupp = eSupportLocal;
        if (GetPropRef(MPO_DIPOLE,(bQM ? iqmQM : iqmBoth),
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
            dip_exp  = value;
            dip_err  = error;
            lot      = mylot;
            //ref      = myref;
            for(m=0; (m<DIM); m++)
            {
                mu_exp[m] = vec[m];
            }
            mu_exp2 = sqr(value);
            if (error <= 0) {
                if (debug)
                    fprintf(debug,"WARNING: Error for %s is %g, assuming it is 10%%.\n",
                            GetMolname().c_str(),error);
                (*nwarn)++;
                error = 0.1*value;
            }
            dip_weight = sqr(1.0/error);
        }
        if (GetPropRef(MPO_DIPOLE,iqmQM,
                       lot,NULL,(char *)"ESP",&value,&error,NULL,NULL,vec,quadrupole) != 0)
        {
            for(m=0; (m<DIM); m++)
            {
                mu_esp[m] = vec[m];
            }
        }
        if (GetProp(MPO_ENERGY,(bQM ? iqmQM : iqmBoth),
                    lot,NULL,(char *)"DHf(298.15K)",&value) != 0)
        {
            Hform = value;
            Emol = value;
            for(ia=0; (ia<topology->atoms.nr); ia++) {
                if (gap.GetValue(*topology->atoms.atomname[ia],
                                  (char *)"exp",(char *)"DHf(0K)",0,&dv0) &&
                    gap.GetValue(*topology->atoms.atomname[ia],
                                  (char *)"exp",(char *)"H(0K)-H(298.15K)",
                                  298.15,&dv298))
                {
                    Emol -= convert2gmx(dv0+dv298,eg2cHartree);
                }
                else {
                    Emol = 0;
                    break;
                }
            }
            if (ia < topology->atoms.nr)
            {
                imm = immNoData;
            }
        }
        else {
            imm = immNoData;
        }
    }
    if (immOK == imm)
    {
        rvec dx;
        rvec_sub(mu_esp,mu_exp,dx);
        if (norm(dx) > dip_toler)
        {
            imm = immQMInconsistency;
        }
    }
    if (immOK == imm)
    {
        if (bQM)
        {
            if (GetPropRef(MPO_QUADRUPOLE,iqmQM,
                           lot,NULL,(char *)"elec",&value,&error,
                           NULL,NULL,vec,quadrupole) == 0)
                imm = immNoQuad;
            else
                copy_mat(quadrupole,Q_exp);
            if (GetPropRef(MPO_QUADRUPOLE,iqmQM,
                           lot,NULL,(char *)"ESP",&value,&error,
                           NULL,NULL,vec,quadrupole) != 0)
                copy_mat(quadrupole,Q_esp);
        }
    }
    if ((immOK == imm) && bESP)
    {
        int  xu,vu;
        
        //snew(symmetric_charges,nalloc);
        symmetric_charges_ = symmetrize_charges(TRUE,&topology->atoms,
                                                &(plist[ftb]),pd,aps,(char *)"");
        
        gr = gmx_resp_init(iModel,FALSE,0,0,GetCharge(),
                           0,0,-1,TRUE,watoms,rDecrZeta,FALSE,1,
                           bFitZeta,FALSE,NULL);
        if (NULL != gr)
        {
            alexandria::CalculationIterator ci;
            alexandria::ElectrostaticPotentialIterator epi;
            
            gmx_resp_add_atom_info(gr,&(topology->atoms),pd);
            gmx_resp_update_atomtypes(gr,&(topology->atoms));
            gmx_resp_add_atom_symmetry(gr,
                                       symmetric_charges_);
            gmx_resp_add_atom_coords(gr,x);
            ci = GetLot(lot);
            if (ci != EndCalculation())
            {
                for(epi=ci->BeginPotential(); (epi<ci->EndPotential()); epi++)
                {
                    /* Maybe not convert to gmx ? */
                    xu = string2unit(epi->GetXYZunit().c_str());
                    vu = string2unit(epi->GetVunit().c_str());
                    if (-1 == xu)
                        xu = eg2cAngstrom;
                    if (-1 == vu)
                        vu = eg2cHartree_e;
                    gmx_resp_add_point(gr,
                                       convert2gmx(epi->GetX(),xu),
                                       convert2gmx(epi->GetY(),xu),
                                       convert2gmx(epi->GetZ(),xu),
                                       convert2gmx(epi->GetV(),vu));
                }
            }
        }
        else
            imm = immRespInit;                          
    }
    else 
    {
        gr = NULL;
    }
    if (immOK == imm)
    {
        if (bPol) {
            gmx_vsite_t *vs = init_vsite(&mtop,cr,TRUE);
            snew(buf,nalloc);
            snew(newexcls,nalloc);
            srenew(x,nalloc);
            add_shells(pd,nalloc,&(mtop.moltype[0].atoms),
                       atype,plist,x,&symtab,
                       &newexcls,smnames);
            mtop.natoms = mtop.moltype[0].atoms.nr;
            mtop.molblock[0].natoms_mol = mtop.natoms;
            excls_to_blocka(nalloc,newexcls,
                            &(mtop.moltype[0].excls));
            mtop_update_cgs(&mtop);
            reset_q(&(topology->atoms));
            shell = init_shell_flexcon(debug,&mtop,0,x);
        }
        else {
            shell = NULL;
            excls_to_blocka(nalloc,excls,&(mtop.moltype[0].excls));
        }
        plist_to_mtop(pd,plist,&mtop);
        snew(f,nalloc);
        fr = mk_forcerec();
        init_forcerec(NULL,oenv,fr,NULL,&ir,&mtop,cr,
                      box,FALSE,NULL,NULL,NULL,NULL,NULL,NULL,TRUE,-1);
        init_state(&state,topology->atoms.nr,1,1,1,0);
        ltop = gmx_mtop_generate_local_top(&(mtop),&(ir));
        md = init_mdatoms(NULL,&mtop,FALSE);
    }
    if (immOK != imm) 
    {
        qgen = NULL; 
    }
    if (immOK != imm) 
    {
        /* Remove temporary data */
        if (NULL != fp)
            fprintf(fp,"Could not add %s because of %s\n",GetMolname().c_str(),
                    immsg(imm));
    }
    else if (NULL != debug)
    {
        fprintf(debug,"Succesfully added %s\n",GetMolname().c_str());
    }
    return imm;
}

}
