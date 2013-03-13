/* 
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
#include <algorithm>
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
#include "gpp_nextnb.h"
#include "toputil.h"
#include "gen_ad.h"
#include "gmx_simple_comm.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "molselect.h"
#include "molprop.hpp"
#include "molprop_util.hpp"
#include "molprop_xml.hpp"
#include "gentop_vsite.hpp"
#include "gentop_core.hpp"

int mp_get_prop_ref(alexandria::MolProp mp,MolPropObservable mpo,int iQM,char *lot,
                    const char *conf,const char *type,double *value,double *error,
                    char **ref,char **mylot,
                    double vec[3],tensor quadrupole)
{
    alexandria::ExperimentIterator ei;
    alexandria::CalculationIterator ci;
    std::string reference,method,name,basisset,program,conformation,expconf;
    
    int n,k,done=0;
    char **qm,**ll;
    
    if ((iQM == iqmExp) || (iQM == iqmBoth)) 
    {
        for(ei=mp.BeginExperiment(); (ei<mp.EndExperiment()); ei++)
        {
            reference = ei->GetReference();
            expconf   = ei->GetConformation();
            
            if ((NULL == conf) || (strcasecmp(conf,expconf.c_str()) == 0))
                done = ei->GetVal(type,mpo,value,error,vec,quadrupole);
        }
        if (done != 0) {
            if (NULL != ref)
                *ref = strdup(reference.c_str());
            if (NULL != mylot)
                *mylot = strdup("Experiment");
        }
    }
    
    if ((done == 0) && ((iQM == iqmBoth) || (iQM == iqmQM)))
    {
        if (NULL != lot) 
        {
            qm = split(':',lot);
            n = 0;
            while ((done == 0) && (NULL != qm[n]))
            {
                ll = split('/',qm[n]);
                if ((NULL != ll[0]) && (NULL != ll[1])) 
                {
                    for(ci=mp.BeginCalculation(); (ci<mp.EndCalculation()); ci++)
                    {
                        basisset     = ci->GetBasisset();
                        method       = ci->GetMethod();
                        reference    = ci->GetReference();
                        conformation = ci->GetConformation();
                      
                        if (((strcasecmp(method.c_str(),ll[0]) == 0) &&
                             (strcasecmp(basisset.c_str(),ll[1]) == 0)) &&
                            ((NULL == conf) || (strcasecmp(conf,conformation.c_str()) == 0)))
                        {
                            done = ci->GetVal(type,mpo,value,error,vec,quadrupole);
                        }
                        if ((done != 0) && (NULL != ref))
                            *ref = strdup(reference.c_str());
                    }
                }
                k = 0;
                while (ll[k] != NULL)
                {
                    sfree(ll[k]);
                    k++;
                }
                sfree(ll);
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
    }
    return done;
}

int mp_get_prop(alexandria::MolProp mp,MolPropObservable mpo,int iQM,char *lot,
                char *conf,char *type,double *value)
{
    double error,vec[3];
    tensor quad;

    return mp_get_prop_ref(mp,mpo,iQM,lot,conf,type,value,&error,
                           NULL,NULL,vec,quad);
}                

enum { elgcOK, elgcNOATOMS, elgcNOTOP, elgcNR };

static int lo_gen_composition(alexandria::MolPropIterator mpi,
                              gmx_poldata_t pd,gmx_atomprop_t aps,
                              gmx_bool *bSpoel,gmx_bool *bMiller,gmx_bool *bBosque,
                              double th_toler,double phi_toler)
{
    std::string miller;
    alexandria::CalculationIterator ci;
    alexandria::CalcAtomIterator cai;
    alexandria::MolecularComposition mci_bosque("bosque"),mci_spoel("spoel"),mci_miller("miller");
    alexandria::AtomNumIterator ani;
    
    for(ci=mpi->BeginCalculation(); (ci<mpi->EndCalculation()); ci++) 
    {
        /* This assumes we have either all atoms or none. 
         * A consistency check could be
         * to compare the number of atoms to the formula */
        for(cai=ci->BeginAtom(); (cai<ci->EndAtom()); cai++)
        {
            alexandria::AtomNum anb(cai->GetName(),1),ans(cai->GetObtype(),1);
            mci_bosque.AddAtom(anb);
            mci_spoel.AddAtom(ans);
            
            miller = gmx_poldata_get_miller_equiv(pd,cai->GetObtype().c_str());
            if (miller.size() > 0)
            {
                alexandria::AtomNum anm(miller,1);
                mci_miller.AddAtom(anm);
            }
        }
    }
    
    *bSpoel  = mci_spoel.CountAtoms() > 0;
    *bMiller = mci_miller.CountAtoms() > 0;
    *bBosque = mci_bosque.CountAtoms() > 0;

    mpi->AddComposition(mci_bosque);
    mpi->AddComposition(mci_spoel);
    mpi->AddComposition(mci_miller);
    
    
    if (debug && *bSpoel) 
    {
        fprintf(debug,"LO_COMP: ");
        for(ani = mci_spoel.BeginAtomNum(); (ani < mci_spoel.EndAtomNum()); ani++)
        {
            fprintf(debug," %s:%d",ani->GetAtom().c_str(),ani->GetNumber());
        }
        fprintf(debug,"\n");
    }
    
    return elgcOK;
}

void generate_composition(std::vector<alexandria::MolProp> mp,gmx_poldata_t pd)
{
    int nOK = 0;
    alexandria::MolPropIterator mpi;
    
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
    {
        mpi->DeleteComposition("spoel");
        mpi->DeleteComposition("miller");
        mpi->DeleteComposition("bosque");
        if (TRUE == mpi->GenerateComposition(pd))
            nOK++;
        else if (debug)
        {
            fprintf(debug,"Failed to make composition for %s\n",
                    mpi->GetMolname().c_str());
        }
    }
    if (mp.size() > 1)
        printf("Generated composition for %d out of %d molecules.\n",nOK,(int)mp.size());
}

void generate_formula(std::vector<alexandria::MolProp> mp,gmx_atomprop_t ap)
{
    int  j,cnumber,an;
    char formula[1280],number[32];
    int  *ncomp;
    real value;
    std::string compname,catom,mform;
    alexandria::MolPropIterator mpi;
    alexandria::MolecularCompositionIterator mci;
    alexandria::AtomNumIterator ani;
    
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) {
        snew(ncomp,110);  
        formula[0] = '\0';
        mci=mpi->SearchMolecularComposition("bosque");
        if (mci != mpi->EndMolecularComposition()) 
        {
            for(ani=mci->BeginAtomNum(); (ani<mci->EndAtomNum()); ani++)
            {
                catom = ani->GetAtom();
                cnumber = ani->GetNumber();
                if (gmx_atomprop_query(ap,epropElement,"???",catom.c_str(),&value))
                {
                    an = gmx_nint(value);
                    range_check(an,0,110);
                    if (an > 0)
                        ncomp[an] += cnumber;
                }
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
        for(j=109; (j>=1); j--) 
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
        mform = mpi->GetFormula();
        if (strlen(formula) > 0) 
        {
            if (debug) 
            {
                if ((mform.size() > 0) && (strcasecmp(formula,mform.c_str()) != 0))
                    fprintf(debug,"Formula '%s' does match '%s' based on composition for %s.\n",
                            mform.c_str(),formula,mpi->GetMolname().c_str());
            }
            mpi->SetFormula(formula);
        }
        else if ((mform.size() == 0) && debug)
            fprintf(debug,"Empty composition and formula for %s\n",
                    mpi->GetMolname().c_str());
        sfree(ncomp);
    }
}
  
alexandria::MolProp atoms_2_molprop(char *molname,int natoms,char **smnames,
                                    gmx_atomprop_t ap,gmx_poldata_t pd)
{
    alexandria::MolProp mp;
    alexandria::MolecularComposition mci("spoel");
    int i;
        
    mp.SetMolname(molname);
        
    for(i=0; (i<natoms); i++) 
    {
        alexandria::AtomNum an(smnames[i],1);
        mci.AddAtom(an);
    }
    mp.AddComposition(mci);
    (void) mp.GenerateFormula(ap);
        
    return mp;  
}

static gmx_bool molprop_2_atoms(alexandria::MolProp mp,gmx_atomprop_t ap,
                                t_symtab *tab,const char *lot,
                                t_atoms *atoms,const char *q_algorithm,
                                rvec **x)
{
    alexandria::CalculationIterator ci;
    alexandria::CalcAtomIterator cai;
    alexandria::AtomicChargeIterator qi;
    std::string molnm;
    
    int myunit;
    double q,xx,yy,zz;
    int    i,natom;
    gmx_bool bDone = FALSE;
    
    (*x)  = NULL;
    molnm = mp.GetMolname();
        
    ci = mp.GetLot(lot);
    if (ci < mp.EndCalculation()) 
    {
        natom = 0;
        init_t_atoms(atoms,mp.NAtom(),FALSE);
        snew(*x,mp.NAtom());
        
        i = 0;
        for(cai=ci->BeginAtom(); (cai<ci->EndAtom()); cai++,i++)
        {
            myunit = string2unit((char *)cai->GetUnit().c_str());
            cai->GetCoords(&xx,&yy,&zz);
            (*x)[natom][XX] = convert2gmx(xx,myunit);
            (*x)[natom][YY] = convert2gmx(yy,myunit);
            (*x)[natom][ZZ] = convert2gmx(zz,myunit);
            
            q = 0;
            for(qi=cai->BeginQ(); (qi<cai->EndQ()); qi++) 
            {
                if (strcmp(qi->GetType().c_str(),q_algorithm) == 0)
                {
                    myunit = string2unit((char *)qi->GetUnit().c_str());
                    q = convert2gmx(qi->GetQ(),myunit);
                }
            }
            
            t_atoms_set_resinfo(atoms,i,tab,molnm.c_str(),1,' ',1,' ');
            atoms->atomname[i] = put_symtab(tab,cai->GetName().c_str());
            atoms->atom[i].atomnumber = gmx_atomprop_atomnumber(ap,cai->GetName().c_str());
            strcpy(atoms->atom[i].elem,gmx_atomprop_element(ap,atoms->atom[i].atomnumber));
            atoms->atom[i].q = q; 
            atoms->atom[i].resind = 0;
        }
        atoms->nres = 1;
        bDone = TRUE;
    }
    if (NULL != debug)
        fprintf(debug,"Tried to convert %s to gromacs. LOT is %s. Natoms is %d\n",
                molnm.c_str(),lot,natom);
    return bDone;
}

gmx_bool molprop_2_topology(alexandria::MolProp mp,gmx_atomprop_t ap,
                            gmx_poldata_t pd,
                            t_symtab *tab,const char *lot,
                            t_topology *top,const char *q_algorithm,
                            rvec **x,t_params plist[F_NRE],
                            int nexcl,t_excls **excls)
{
    alexandria::BondIterator bi;
    
    int natom,ftb;
    t_param b;
    t_nextnb nnb;
        
    /* Get atoms */
    if (FALSE == molprop_2_atoms(mp,ap,tab,lot,&(top->atoms),q_algorithm,x))
        return FALSE;
    
    /* Store bonds in harmonic potential list first, update type later */
    ftb = F_BONDS;
    memset(&b,0,sizeof(b));
    for(bi=mp.BeginBond(); (bi<mp.EndBond()); bi++)
    {
        b.a[0] = bi->GetAi() - 1;
        b.a[1] = bi->GetAj() - 1;
        add_param_to_list(&(plist[ftb]),&b);
    }
    natom = mp.NAtom();
    
    /* Make Angles and Dihedrals */
    snew((*excls),natom);
    init_nnb(&nnb,natom,nexcl+2);
    gen_nnb(&nnb,plist);
    //detect_rings(&plist[F_BONDS],atoms->nr,bRing);
    //nbonds = plist[F_BONDS].nr;
    print_nnb(&nnb,"NNB");
    gen_pad(&nnb,&top->atoms,TRUE,TRUE,TRUE,nexcl,plist,*excls,NULL,FALSE);
    generate_excls(&nnb,nexcl,*excls);
    done_nnb(&nnb);

    return TRUE;
}

int my_strcmp(const char *a,const char *b)
{
    if ((NULL != a) && (NULL != b))
        return strcasecmp(a,b);
    else
        return 1;
}

void merge_doubles(std::vector<alexandria::MolProp> &mp,char *doubles,
                   gmx_bool bForceMerge)
{
    alexandria::MolPropIterator mpi,mmm[2];
    std::string molname[2];
    std::string form[2];
    
    FILE *fp;
    int  i,ndouble=0;
    gmx_bool bForm,bName,bDouble;
    int  cur = 0;
#define prev (1-cur)
  
    if (NULL != doubles)
        fp = fopen(doubles,"w");
    else
        fp = NULL;
    i = 0;
    for(mpi=mp.begin(); (mpi<mp.end()); )
    {
        bDouble = FALSE;
        mpi->Dump(debug);
        mmm[cur] = mpi;
        molname[cur] = mpi->GetMolname();
        form[cur]    = mpi->GetFormula();
        if (i > 0) {
            bForm = (form[prev] == form[cur]);
            bName = (molname[prev] == molname[cur]);
            if (bName) {
                if (!bForm && debug)
                    fprintf(debug,"%s %s with formulae %s - %s\n",
                            bForceMerge ? "Merging molecules" : "Found molecule",
                            molname[prev].c_str(),form[prev].c_str(),form[cur].c_str());
                if (bForceMerge || bForm) {
                    if (fp)
                        fprintf(fp,"%5d  %s\n",ndouble+1,molname[prev].c_str());
                    mmm[prev]->Merge(*(mmm[cur]));
                    mpi = mp.erase(mmm[cur]);
                    
                    bDouble = TRUE;
                    ndouble++;
                    if (mpi != mp.end())
                    {
                        mmm[cur]      = mpi;
                        molname[cur]  = mmm[cur]->GetMolname();
                        form[cur]     = mmm[cur]->GetFormula();
                    }
                }
            }
        }
        if (!bDouble) {
            cur = prev;
            mpi++;
            i++;
        }
    }
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++)
    {
        mpi->Dump(debug);
    }
    if (fp)
        fclose(fp);
    printf("There were %d double entries, leaving %d after merging.\n",
           ndouble,(int)mp.size());
}

static void dump_mp(std::vector<alexandria::MolProp> mp)
{
    alexandria::MolPropIterator mpi;
    FILE *fp;
  
    fp = fopen("dump_mp.dat","w");
  
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
    {
        fprintf(fp,"%-20s  %s\n",mpi->GetFormula().c_str(),
                mpi->GetMolname().c_str());
    }
  
    fclose(fp);
}

void merge_xml(int nfile,char **filens,
               std::vector<alexandria::MolProp> &mpout, 
               char *outf,char *sorted,char *doubles,
               gmx_atomprop_t ap,gmx_poldata_t pd,
               gmx_bool bForceMerge,gmx_bool bForceGenComp,
               double th_toler,double ph_toler)
{
    std::vector<alexandria::MolProp> mp;
    alexandria::MolPropIterator mpi;
    int       i,npout=0,tmp;
  
    for(i=0; (i<nfile); i++)
    {
        if (!gmx_fexist(filens[i]))
            continue;
        MolPropRead(filens[i],mp);
        generate_composition(mp,pd);
        generate_formula(mp,ap);
        for(mpi=mp.begin(); (mpi<mp.end());)
        {
            mpout.push_back(*mpi);
            mpi = mp.erase(mpi);
        }
    }
    tmp = mpout.size();
    printf("mpout.size() = %u mpout.max_size() = %u\n",
           (unsigned int)mpout.size(),(unsigned int)mpout.max_size());
    for(mpi=mpout.begin(); (mpi<mpout.end()); mpi++)
    {
        mpi->Dump(debug);
    }
    MolPropSort(mpout,MPSA_MOLNAME,NULL,NULL);
    merge_doubles(mpout,doubles,bForceMerge);
    printf("There are %d total molecules before merging, %d after.\n",
           tmp,(int)mpout.size());
      
    if (outf) 
    {
        printf("There are %d entries to store in output file %s\n",npout,outf);
        MolPropWrite(outf,mpout,FALSE);
    }
    if (sorted) 
    {
        MolPropSort(mpout,MPSA_FORMULA,NULL,NULL);
        MolPropWrite(sorted,mpout,FALSE);
        dump_mp(mpout);
    }
}

static bool comp_mp_molname(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string mma = ma.GetMolname();
    std::string mmb = mb.GetMolname();
  
    //return 0;
    return (mma.compare(mmb) < 0); 
        //(strcmp(mma.c_str(),mmb.c_str()) <= 0);
}

static bool comp_mp_formula(alexandria::MolProp ma,
                            alexandria::MolProp mb)
{
    std::string fma = ma.GetFormula();
    std::string fmb = mb.GetFormula();
  
    if (fma.compare(fmb) < 0)
        return true;
    else
        return comp_mp_molname(ma,mb);
}

gmx_atomprop_t my_aps;

static bool comp_mp_elem(alexandria::MolProp ma,
                         alexandria::MolProp mb)
{
    int i;
    alexandria::MolecularCompositionIterator mcia,mcib;
    std::string bosque("bosque"),C("C");
 
    mcia = ma.SearchMolecularComposition(bosque);
    mcib = mb.SearchMolecularComposition(bosque);

    if ((mcia != ma.EndMolecularComposition()) &&
        (mcib != mb.EndMolecularComposition()))
    {
        if (mcia->CountAtoms(C) < mcib->CountAtoms(C))
            return true;
  
        for(i=1; (i<=109); i++) {
            if (i != 6) {
                std::string elem(gmx_atomprop_element(my_aps,i));
                
                if (mcia->CountAtoms(elem) < mcib->CountAtoms(elem))
                    return true;
            }
        }
    }
    return comp_mp_molname(ma,mb);
}

gmx_molselect_t my_gms;

static bool comp_mp_selection(alexandria::MolProp ma,alexandria::MolProp mb)
{
    int ia,ib;

    ia = gmx_molselect_index(my_gms,ma.GetIupac().c_str());
    ib = gmx_molselect_index(my_gms,mb.GetIupac().c_str());
    
    if (ia < ib)
        return true;
    else
        return false;
}

void MolPropSort(std::vector<alexandria::MolProp> &mp,
                 MolPropSortAlgorithm mpsa,gmx_atomprop_t apt,
                 gmx_molselect_t gms)
{
    switch(mpsa) {
    case MPSA_MOLNAME:
        std::sort(mp.begin(),mp.end(),comp_mp_molname);
        break;
    case MPSA_FORMULA:
        std::sort(mp.begin(),mp.end(),comp_mp_formula);
        break;
    case MPSA_COMPOSITION:
        my_aps = apt;
        std::sort(mp.begin(),mp.end(),comp_mp_elem);
        my_aps = NULL;
        break;
    case MPSA_SELECTION:
        if (NULL != gms) 
        {
            my_gms = gms;
            std::sort(mp.begin(),mp.end(),comp_mp_selection);
        }
        else
            gmx_fatal(FARGS,"Need molecule selection to sort on");
        break;
    default:
        gmx_incons("Invalid algorithm for MolPropSort");
    }
}

static void add_qmc_conf(t_qmcount *qmc,const char *conformation)
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

static int get_qmc_count(t_qmcount *qmc,const char *method,const char *basis,const char *type)
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

static void add_qmc_calc(t_qmcount *qmc,const char *method,const char *basis,const char *type)
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

t_qmcount *find_calculations(std::vector<alexandria::MolProp> mp,
                             MolPropObservable mpo,
                             const char *fc_str)
{
    alexandria::MolPropIterator     mpi;
    alexandria::ExperimentIterator  ei;
    alexandria::CalculationIterator ci;
    
    const char *method,*basis;
    char **qm,**ll;
    int  i,n;
    double value,error,vec[3];
    tensor quadrupole;
    t_qmcount *qmc;
    std::vector<std::string> types;
    std::vector<std::string>::iterator ti;
    
    snew(qmc,1);
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
    {
        for(ei=mpi->BeginExperiment(); (ei<mpi->EndExperiment()); ei++) 
        {
            add_qmc_conf(qmc,ei->GetConformation().c_str());
        }
        for(ci=mpi->BeginCalculation(); (ci<mpi->EndCalculation()); ci++) 
        {
            add_qmc_conf(qmc,ci->GetConformation().c_str());
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
                for(ti=types.begin(); (ti<types.end()); ti++)
                {
                    if (0 == strcasecmp(ti->c_str(),ll[2]))
                        break;
                }
                if (ti == types.end()) 
                {
                    types.push_back(ll[2]);
                }
            }
            n++;
        }
    }
    
    for(mpi=mp.begin(); (mpi<mp.end()); mpi++) 
    {
        for(ci=mpi->BeginCalculation(); (ci<mpi->EndCalculation()); ci++) 
        {
            method = ci->GetMethod().c_str();
            basis  = ci->GetBasisset().c_str();
            for(ti=types.begin(); (ti<types.end()); ti++) 
            {
                if ((NULL == fc_str) || (get_qmc_count(qmc,method,basis,ti->c_str()) > 0))
                {
                    if (ci->GetVal(ti->c_str(),mpo,&value,&error,vec,quadrupole) == 1) 
                    {
                        add_qmc_calc(qmc,method,basis,ti->c_str());
                    }
                }
            }
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

