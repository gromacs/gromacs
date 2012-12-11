/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gentop.c,v 1.26 2009/05/20 10:48:03 spoel Exp $
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
#include "typedefs.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "strdb.h" 
#include "futil.h"
#include "poldata.h"
#include "atomprop.h"
#include "molprop.h"
#include "molprop_util.h"
#include "string2.h"
#include "gmx_gauss_io.h"
#include "vec.h"

// Include Open Babel classes for OBMol and OBConversion
#ifdef HAVE_LIBOPENBABEL2
#include <iostream>
#include <fstream>
#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/residue.h>
#include <openbabel/obiter.h>
#include <openbabel/obconversion.h>
#include <openbabel/math/vector3.h>

using namespace std;

typedef struct gap {
    char *element, *method,*desc;
    real temp;
    real value;
} gap_t;

typedef struct gau_atomprop
{
    int   ngap;
    gap_t *gap;
} gau_atomprop;


typedef struct {
    int  id;
    real x,y,z;
    real V;
} t_espv;

static int esp_comp(const void *a,const void *b)
{
    t_espv *ea = (t_espv *)a;
    t_espv *eb = (t_espv *)b;
  
    if (ea->V < eb->V)
        return -1;
    else if (ea->V > ea->V)
        return 1;
    return 0;
}
 
static OpenBabel::OBConversion *read_babel(const char *g98,OpenBabel::OBMol *mol)
{
    ifstream         g98f;
    char             *g98z,*ptr;
  
    if (FALSE == gmx_fexist(g98)) 
    {
        snew(g98z,strlen(g98)+4);
        strcpy(g98z,g98);
        strcat(g98z,".gz");
        g98f.open(g98z,ios::in);
        sfree(g98z);
    }
    else 
    {
        g98f.open(g98,ios::in);
    }
    if (!g98f.is_open()) 
    {
        gmx_fatal(FARGS,"Can not open file %s for reading",g98);
    }
  
    // Read from g98f
    OpenBabel::OBConversion *conv = new OpenBabel::OBConversion(&g98f,&cout);
  
    // Try to set input format to G98 file if it is not clear from the extension,
    // that means, this code will equally well read sdf, pdb etc.
    ptr = (char *)strrchr(g98,'.');
    if ((NULL != ptr) && (strlen(ptr) >= 2))
        ptr++;
    else
        ptr = (char *)"g98";
    if (conv->SetInFormat(ptr))
    {
        if (conv->Read(mol))
        {
            g98f.close();
            
            return conv; // exit with success
        }
        else 
        {
            cerr << "Could not read input file " << g98 << " with OpenBabel2." << endl;
        }
    }
    else 
    {
        cerr << "Input file " << g98 << " has incomprehensible format." << endl;
    }
    g98f.close();
    
    return NULL;
}

static gmx_molprop_t gmx_molprop_read_babel(const char *g98,
                                            gmx_atomprop_t aps,gmx_poldata_t pd,
                                            char *molnm,char *iupac,char *conformation,
                                            char *basisset,real th_toler,real ph_toler,
                                            int maxpot,gmx_bool bVerbose)
{
    /* Read a gaussian log file */
    OpenBabel::OBMol mol;
    OpenBabel::OBAtomIterator OBai;
    OpenBabel::OBBondIterator OBbi;
    OpenBabel::OBConversion *conv;
    OpenBabel::OBAtom *OBa,*OBa1,*OBa2;
    OpenBabel::OBBond *OBb;
    OpenBabel::OBPairData *OBpd,*OBdhf;
    OpenBabel::OBVectorData *dipole;
    OpenBabel::OBMatrixData *quadrupole,*pol_tensor;
    OpenBabel::OBFreeGrid *esp;
    OpenBabel::OBElementTable *OBet;
    std::string formula,attr,value,inchi;
  
    t_espv *espv = NULL;
    const char *reference="Spoel2013a",*unknown="unknown";
    char *program,*method,*basis,*charge_model,*ptr,*g98ptr;
    int i,ii0,calcref,atomref,atomid,bondid,natom;
    gmx_molprop_t mpt;  
    double ii,deltai,dval;
    int k;
    const char *etypes[] = { "DHf(0K)", "DHf(298.15K)" };
  
    conv = read_babel(g98,&mol);
    if (NULL == conv)
    {
        fprintf(stderr,"Failed reading %s\n",g98);
        return NULL;
    }
    delete conv;
    
    /* Create new molprop calculation */ 
    mpt = gmx_molprop_init();

    conv = new OpenBabel::OBConversion(&cin,&cout);
    // Now extract classification info.
    if (conv->SetOutFormat("fpt"))
    {
        vector<string> vs;
        string ss;
        const char *exclude[] = { ">", "C_ONS_bond", "Rotatable_bond", "Conjugated_double_bond", "Conjugated_triple_bond", "Chiral_center_specified", "Cis_double_bond", "Bridged_rings", "Conjugated_tripple_bond", "Trans_double_bond" };
#define nexclude (sizeof(exclude)/sizeof(exclude[0]))
        char *dup,*ptr;
        unsigned int i,j;
        
        conv->AddOption("f", OpenBabel::OBConversion::OUTOPTIONS, "FP4");
        conv->AddOption("s");
        conv->Convert();
        ss = conv->WriteString(&mol,false);
        if (OpenBabel::tokenize(vs,ss))
        {
            for(i=0; (i<vs.size()); i++)
            {
                for(j=0; (j<nexclude); j++)
                {
                    if(strcasecmp(exclude[j],vs[i].c_str()) == 0)
                        break;
                }
                if (j == nexclude)
                {
                    dup = strdup(vs[i].c_str());
                    while (NULL != (ptr = strchr(dup,'_')))
                    {
                        *ptr = ' ';
                    }
                    gmx_molprop_add_category(mpt,dup);
                }
            }
        }
    }
    
    // Get bondorders.
    mol.PerceiveBondOrders();
    
    OBpd = (OpenBabel::OBPairData *)mol.GetData("basis");
    if ((NULL != basisset) && (strlen(basisset) > 0))
        basis = basisset;
    else if (NULL != OBpd)
    {
        basis = strdup(OBpd->GetValue().c_str());
        if (NULL != (ptr = strstr(basis," (5D, 7F)")))
            *ptr = '\0';
    }
    else
        basis = strdup(unknown);
  
    OBpd = (OpenBabel::OBPairData *)mol.GetData("program");
    if (NULL != OBpd)
        program = strdup(OBpd->GetValue().c_str());
    else
        program = strdup(unknown);
  
    OBpd = (OpenBabel::OBPairData *)mol.GetData("method");
    if (NULL != OBpd)
        method = strdup(OBpd->GetValue().c_str());
    else 
        method = strdup(unknown);
  
    g98ptr = (char *) strrchr(g98,'/');
    if (NULL == g98ptr) 
        g98ptr = (char *)g98;
    else {
        g98ptr++;
        if (strlen(g98ptr) == 0)
            g98ptr = (char *)g98;
    }
    gmx_molprop_add_calculation(mpt,program,method,basis,reference,
                                conformation,g98ptr,&calcref);
  
    gmx_molprop_set_charge(mpt,mol.GetTotalCharge());
    gmx_molprop_set_mass(mpt,mol.GetMolWt());
    gmx_molprop_set_multiplicity(mpt,mol.GetTotalSpinMultiplicity());
  
    formula = mol.GetFormula();
    if (formula.size() > 0)
        gmx_molprop_set_formula(mpt,formula.c_str());
    else
        gmx_molprop_set_formula(mpt,unknown);
  
    conv->SetOutFormat("inchi");
    inchi = conv->WriteString(&mol,true);

    delete conv;

    if (inchi.size() > 0)
        gmx_molprop_set_inchi(mpt,inchi.c_str());
    else
        gmx_molprop_set_inchi(mpt,unknown);
  
    if (NULL != molnm)
        gmx_molprop_set_molname(mpt,molnm);
    else
        gmx_molprop_set_molname(mpt,unknown);
    
    if (NULL != iupac)
        gmx_molprop_set_iupac(mpt,iupac);
    else
        gmx_molprop_set_iupac(mpt,unknown);
    
    for(k=0; (k<2); k++) 
    {
        OBdhf = (OpenBabel::OBPairData *)mol.GetData(etypes[k]);
        if (NULL != OBdhf)
        {
            value = OBdhf->GetValue();
            dval = convert2gmx(atof(value.c_str()),eg2c_kcal_mole);
            gmx_molprop_add_energy(mpt,calcref,etypes[k],
                                   unit2string(eg2c_kJ_mole),dval,0);
        }
    }
  
    /* Now add properties by extracting them from the OpenBabel structure */
    OBpd = (OpenBabel::OBPairData *) mol.GetData("PartialCharges");
    if (NULL != OBpd)
        charge_model = strdup(OBpd->GetValue().c_str());
    else
        charge_model = strdup(unknown);
  
    OBet = new OpenBabel::OBElementTable();
  
    OBai = mol.BeginAtoms();
    atomid = 1;
    for (OBa = mol.BeginAtom(OBai); (NULL != OBa); OBa = mol.NextAtom(OBai)) {
        gmx_molprop_calc_add_atom(mpt,calcref,
                                  OBet->GetSymbol(OBa->GetAtomicNum()),
                                  OBa->GetType(),atomid,&atomref);
        gmx_molprop_calc_set_atomcoords(mpt,calcref,atomref,unit2string(eg2c_pm),
                                        100*OBa->x(),100*OBa->y(),100*OBa->z());
        gmx_molprop_calc_set_atomcharge(mpt,calcref,atomref,charge_model,
                                        "e",OBa->GetPartialCharge());
        atomid++;
    }
    delete OBet;
  
    OBbi = mol.BeginBonds();
    bondid = 1;
    for (OBb = mol.BeginBond(OBbi); (NULL != OBb); OBb = mol.NextBond(OBbi)) {
        OBa1 = OBb->GetBeginAtom();
        OBa2 = OBb->GetEndAtom(); 
        gmx_molprop_add_bond(mpt,1+OBa1->GetIndex(),1+OBa2->GetIndex(),OBb->GetBondOrder());
        bondid++;
    }
  
    // Dipole
    dipole = (OpenBabel::OBVectorData *) mol.GetData("Dipole Moment");
    if (NULL != dipole)
    {
        OpenBabel::vector3 v3 = dipole->GetData();
        gmx_molprop_add_dipole(mpt,calcref,"electronic",
                               unit2string(eg2c_Angstrom),
                               v3.GetX(),v3.GetY(),v3.GetZ(),
                               v3.length(),0.0);
    }
  
    // Quadrupole
    quadrupole = (OpenBabel::OBMatrixData *) mol.GetData("Traceless Quadrupole Moment");
    if (NULL != quadrupole)
    {
        OpenBabel::matrix3x3 m3 = quadrupole->GetData();
        double mm[9];
        m3.GetArray(mm);
        gmx_molprop_add_quadrupole(mpt,calcref,"electronic",
                                   unit2string(eg2c_Buckingham),
                                   mm[0],mm[4],mm[8],mm[1],mm[2],mm[5]);
    }
  
    // Polarizability
    pol_tensor = (OpenBabel::OBMatrixData *) mol.GetData("Exact polarizability");
    if (NULL != pol_tensor)
    {
        OpenBabel::matrix3x3 m3 = pol_tensor->GetData();
        double mm[9],alpha,fac;
        int i;
        m3.GetArray(mm);
        fac = 1000*pow(convert2gmx(1,eg2c_Bohr),3);
        for(i=0; (i<9); i++)
            mm[i] *= fac; 
        //cout << "fac = " << fac << "\n";
        alpha = (mm[0]+mm[4]+mm[8])/3.0;
        gmx_molprop_add_polar(mpt,calcref,"electronic",
                              unit2string(eg2c_Angstrom3),
                              mm[0],mm[4],mm[8],alpha,0);
    }
  
    // Electrostatic potential
    esp = (OpenBabel::OBFreeGrid *) mol.GetData("Electrostatic Potential");
    if (NULL != esp)
    {
        OpenBabel::OBFreeGridPoint *fgp;
        OpenBabel::OBFreeGridPointIterator fgpi;
        int espid=0;
      
        fgpi = esp->BeginPoints();
        snew(espv,esp->NumPoints());
        for(fgp = esp->BeginPoint(fgpi); (NULL != fgp); fgp = esp->NextPoint(fgpi))
        {
            espv[espid].id = espid+1;
            espv[espid].x = 100*fgp->GetX();
            espv[espid].y = 100*fgp->GetY();
            espv[espid].z = 100*fgp->GetZ();
            espv[espid].V = fgp->GetV();
            espid++;
        }
        natom = mol.NumAtoms();
      
        if ((maxpot > 0) && (maxpot < esp->NumPoints())) 
        {
            qsort(espv+natom,espid-natom,sizeof(espv[0]),esp_comp);
            deltai = (esp->NumPoints()-natom)/maxpot;
        }
        else 
        {
            maxpot = esp->NumPoints();
            deltai = 1;
        }
        for(i=0; (i<natom); i++) 
        {
            gmx_molprop_add_potential(mpt,calcref,
                                      unit2string(eg2c_pm),
                                      unit2string(eg2c_Hartree_e),
                                      espv[i].id,
                                      espv[i].x,espv[i].y,espv[i].z,
                                      espv[i].V);
        }
        ii = natom;
        for(; (i<maxpot); i++) 
        {
            /* Convert to integer */
            ii0 = ii;
            gmx_molprop_add_potential(mpt,calcref,
                                      unit2string(eg2c_pm),
                                      unit2string(eg2c_Hartree_e),
                                      i,
                                      espv[ii0].x,espv[ii0].y,espv[ii0].z,
                                      espv[ii0].V);
            ii+=deltai;
        }
        sfree(espv);
    }
  
    return mpt;
}

#endif

/******************************************************************* 
 ******************************************************************* 
 * Code for in case we do not have openbabel 
 ******************************************************************* 
 *******************************************************************/

static int get_lib_file(const char *db,char ***strings)
{
    FILE *in;
    char **ptr=NULL;
    char buf[STRLEN];
    int  i,nstr,maxi;

    in=libopen(db);
  
    i=maxi=0;
    while (fgets2(buf,STRLEN-1,in)) {
        if (i>=maxi) {
            maxi+=50;
            srenew(ptr,maxi);
        }
        ptr[i] = strdup(buf);
        i++;
    }
    nstr=i;
    ffclose(in);
    srenew(ptr,nstr);
    *strings=ptr;
  
    return nstr;
}

/* read composite method atom data */
gau_atomprop_t read_gauss_data()
{
    char **strings=NULL,**ptr;
    int nstrings,i,j=0;
    gau_atomprop_t gaps;

    snew(gaps,1);
    nstrings = get_lib_file("atomization_energies.dat",&strings);
    
    /* First loop over strings to deduct basic stuff */
    for(i=0; (i<nstrings); i++) 
    {
        if ( strings[i][0] == '#') {
            continue;
        } 

        ptr = split('|', strings[i]);
        if ((NULL != ptr) && 
            (NULL != ptr[0]) && (NULL != ptr[1]) &&
            (NULL != ptr[2]) && (NULL != ptr[3]) &&
            (NULL != ptr[4]))
        {
            srenew(gaps->gap,j+1); 
            gaps->gap[j].element = ptr[0];
            gaps->gap[j].method = ptr[1];
            gaps->gap[j].desc = ptr[2];
            gaps->gap[j].temp = atof(ptr[3]);
            gaps->gap[j].value = atof(ptr[4]);
            j++;
        }
        sfree(strings[i]);
    }
    sfree(strings);
    
    gaps->ngap = j;
    return gaps;
}

void done_gauss_data(gau_atomprop_t gaps)
{
    fprintf(stderr,"Please clean up your gau_atomprop_t\n");
}

static int gau_atomprop_get_value(gau_atomprop_t gaps,char *element,char *method,
                                  char *desc,double temp,double *value)
{
    int i,found;
    double ttol = 0.01; /* K */
    
    found = 0;
    for(i=0; (i<gaps->ngap); i++) 
    {
        if ((0 == strcasecmp(gaps->gap[i].element,element)) &&
            (0 == strcasecmp(gaps->gap[i].method,method)) &&
            (0 == strcasecmp(gaps->gap[i].desc,desc)) &&
            (fabs(temp - gaps->gap[i].temp) < ttol)) 
        {
            *value = gaps->gap[i].value;
            found = 1;
            break;
        }
    }
    return found;
}

static int gmx_molprop_add_dhform(gmx_molprop_t mpt,int calcref,
                                  gau_atomprop_t gaps,
                                  char *method,double temp,
                                  double eg34,double ezpe,double etherm,
                                  gmx_bool bVerbose)
{
    char   *atomname;
    char   desc[128];
    int    natom,ntot;
    double vm,ve,vdh,vvm,vve,vvdh,ee,e0Hartree,eTHartree;
    
    if (gmx_molprop_get_natom(mpt) == 0)
        return 0;

    ntot = 0;
    vm = ve = vdh = 0;
    sprintf(desc,"%s(0K)",method);
    while (0  < gmx_molprop_get_composition_atom(mpt,"bosque",
                                                 &atomname,&natom))
    {
        /* There are natom atoms with name atomname */
        if ((1 == gau_atomprop_get_value(gaps,atomname,method,desc,0,&vvm)) &&
            (1 == gau_atomprop_get_value(gaps,atomname,(char *)"exp",
                                         (char *)"DHf(0K)",0,&vve)) &&
            (1 == gau_atomprop_get_value(gaps,atomname,(char *)"exp",
                                         (char *)"H(0K)-H(298.15K)",298.15,&vvdh)))
        {
         
            vm  += natom*vvm;
            ve  += natom*vve;
            vdh += natom*vvdh;
            ntot += natom;
        }
        else
        {
            return 0;
        }
    }
    /* Make sure units match! */
    e0Hartree = eg34+ve-vm;
    eTHartree = eg34-ezpe+etherm+ve-vm-vdh;
    
    ee = convert2gmx(e0Hartree,eg2c_Hartree);
    gmx_molprop_add_energy(mpt,calcref,"DHf(0K)",unit2string(eg2c_kJ_mole),ee,0);
    if (bVerbose)
        printf("natom %d e0Hartree %f eTHartree %f eg34 %f ve %f vm %f vdh %f ezpe %f etherm %f \n",
               ntot,e0Hartree,eTHartree,eg34,ve,vm,vdh,ezpe,etherm);
           
    ee = convert2gmx(eTHartree,eg2c_Hartree);
    gmx_molprop_add_energy(mpt,calcref,"DHf(298.15K)",unit2string(eg2c_kJ_mole),ee,0);
    
    return 1;
}

/* Read a line from a G03/G09 composite method (G3, G4, etc) record */
static gmx_bool gau_comp_meth_read_line(char *line,real *temp,real *pres)
{
    char **ptr1,**ptr2;
    int i;
    
    ptr1 = split('=', line);
    if ((NULL != ptr1) && (NULL != ptr1[1]))
    {
        ptr2 = split(' ',ptr1[1]);
        if ((NULL != ptr2) && (NULL != ptr2[0]) && (NULL != ptr1[2])) {
            *temp = atof(ptr2[0]);
            *pres = atof(ptr1[2]);
            for(i=0; (ptr2[i] != NULL); i++)
                sfree(ptr2[i]);
            sfree(ptr2);
        }
        for(i=0; (ptr1[i] != NULL); i++)
            sfree(ptr1[i]);
        sfree(ptr1);
        
        return TRUE;
    }
    return FALSE;
}

static gmx_bool read_polar(char *str,tensor T)
{
    char **ptr;
    int  k;
    gmx_bool bRes;
    
    bRes = TRUE;
    ptr = split(' ',str);
    if (NULL == ptr)
        return FALSE;
        
    if (NULL != ptr[2])
        T[XX][XX] = atof(ptr[2]);
    else 
        bRes = FALSE;
    if (NULL != ptr[3])
        T[YY][XX] = T[XX][YY] = atof(ptr[3]);
    else 
        bRes = FALSE;
    if (NULL != ptr[4])
        T[YY][YY] = atof(ptr[4]);
    else 
        bRes = FALSE;
    if (NULL != ptr[5])
        T[XX][ZZ] = T[ZZ][XX] = atof(ptr[5]);
    else 
        bRes = FALSE;
    if (NULL != ptr[6])
        T[YY][ZZ] = T[ZZ][YY] = atof(ptr[6]);
    else 
        bRes = FALSE;
    if (NULL != ptr[7])
        T[ZZ][ZZ] = atof(ptr[7]);
    else 
        bRes = FALSE;
    for(k=0; (k<=7); k++)
        if (NULL != ptr[k])
            sfree(ptr[k]);
    sfree(ptr);
    
    return bRes;
}

static gmx_bool read_dipole(char *str,rvec mu)
{
    char **ptr;
    int  k;
    gmx_bool bRes;
    
    bRes = TRUE;
    ptr = split(' ',str);
    if (NULL == ptr)
        return FALSE;
        
    if (NULL != ptr[1])
        mu[XX] = atof(ptr[1]);
    else 
        bRes = FALSE;
    if (NULL != ptr[3])
        mu[YY] = atof(ptr[3]);
    else 
        bRes = FALSE;
    if (NULL != ptr[5])
        mu[ZZ] = atof(ptr[5]);
    else 
        bRes = FALSE;
    for(k=0; (k<=7); k++)
        if (NULL != ptr[k])
            sfree(ptr[k]);
    sfree(ptr);
    
    return bRes;
}

static gmx_bool read_quad(char *str1,char *str2,tensor Q)
{
    gmx_bool bRes;
    rvec q1,q2;
    
    bRes = read_dipole(str1,q1);
    bRes = bRes && read_dipole(str2,q2);
    Q[XX][XX] = q1[XX];
    Q[YY][YY] = q1[YY];
    Q[ZZ][ZZ] = q1[ZZ];
    Q[XX][YY] = Q[YY][XX] = q2[XX];
    Q[XX][ZZ] = Q[ZZ][XX] = q2[YY];
    Q[YY][ZZ] = Q[ZZ][YY] = q2[ZZ];
    
    return bRes;
}

static gmx_molprop_t gmx_molprop_read_log(const char *fn,
                                          gmx_atomprop_t aps,gmx_poldata_t pd,
                                          char *molnm,char *iupac,char *conformation,
                                          char *basisset,gau_atomprop_t gaps,
                                          real th_toler,real ph_toler,
                                          int maxpot,gmx_bool bVerbose)
{
    /* Read a gaussian log file */
    char **strings=NULL;
    char sbuf[STRLEN];
    int  nstrings,atomicnumber=-1;
    double ee,x,y,z,V,myener,mass=0;
    int i,j,k,kk,anumber,natom,nesp,nelprop,charge,nfitpoints=-1;
    int calcref=NOTSET,atomref=NOTSET,len;
    gmx_bool bWarnESP=FALSE;
    gmx_bool bAtomicCenter;
    gmx_molprop_t mpt;
    real mm;
    char *atomname,*ginc,*hfener,*mp2ener,*g2ener,*g3ener,*g4ener,*cbsener;
    char *reference = (char *)"This Work";
    char *program=NULL,*method=NULL,*basis=NULL;
    char **ptr,**qtr,*mymeth;
    t_espv *espv=NULL;
    real temp,pres,ezpe,ezpe2,etherm,etherm2,comp_0K,comp_energy,comp_enthalpy,comp_free_energy,ii,deltai;
    gmx_bool bEtherm = FALSE, bEzpe = FALSE, bTemp = FALSE;
    gmx_bool bPolar, bQuad, bDipole;
    tensor polar,quad;
    rvec dipole;
    int status,ii0;

    nstrings = get_file(fn,&strings);
    
    /* Create new calculation */ 
    mpt = gmx_molprop_init();

    natom   = 0;
    nesp    = 0;
    nelprop = NOTSET;
    charge  = NOTSET;
    hfener  = NULL;
    mp2ener = NULL;
    g2ener  = NULL;
    g3ener  = NULL;
    g4ener  = NULL;
    cbsener = NULL;
    bPolar  = bQuad = bDipole = FALSE;
  
    /* First loop over strings to deduct basic stuff */
    for(i=0; (i<nstrings); i++) 
    {
        if (NULL != strstr(strings[i],"Revision")) 
        {
            program = strdup(strings[i]);
            trim(program);
            len = strlen(program);
            if ((len > 0) && (program[len-1] == ',')) 
            {
                program[len-1] = '\0';
            }
        }
        else if ((NULL != strstr(strings[i],"Standard basis:")) &&
                 (NULL == basisset))
        {
            ptr = split(' ',strings[i]);
            if (NULL != ptr[2]) 
            {
                basis = strdup(ptr[2]);
            }
        }
        else if (NULL != strstr(strings[i],"Temperature=")) 
        {
            status = gau_comp_meth_read_line(strings[i],&temp,&pres);
            if (bVerbose)
                printf("na gau_(): temp %f pres %f\n", temp, pres);
            bTemp = TRUE;
        }
        else if (NULL != strstr(strings[i],"Exact polarizability")) 
        {
            bPolar = read_polar(strings[i],polar);
        }
        else if (NULL != strstr(strings[i],"Dipole moment")) 
        {
            bDipole = read_dipole(strings[i+1],dipole);
            i += 1;
        }
        else if (NULL != strstr(strings[i],"Traceless Quadrupole moment")) 
        {
            bQuad = read_quad(strings[i+1],strings[i+2],quad);
            i += 2;
        }
        else if (NULL != strstr(strings[i],"E(ZPE)="))
        {
            status = gau_comp_meth_read_line(strings[i],&ezpe2,&etherm2);
            if (bVerbose)
                printf("na gau_(): ezpe2 %f etherm2 %f \n", ezpe2,etherm2);
        }
        else if (NULL != strstr(strings[i],"Zero-point correction=")) 
        {
            ptr = split(' ',strings[i]);
            bEzpe = TRUE;
            if (NULL != ptr[2]) 
            {
                ezpe = atof(strdup(ptr[2]));
            }
            if (bVerbose)
                printf("na gau_(): ezpe %f \n", ezpe);
        }
        else if (NULL != strstr(strings[i],"Thermal correction to Enthalpy=")) 
        {
          
            ptr = split(' ',strings[i]);
            if (NULL != ptr[4]) 
            {
                etherm = atof(strdup(ptr[4]));
            }
            bEtherm = TRUE;
            if (bVerbose)
                printf("na gau_(): etherm %f \n", etherm);
        }
        else if ((NULL != strstr(strings[i],"G2(0 K)=")) ||
                 (NULL != strstr(strings[i],"G3(0 K)=")) ||
                 (NULL != strstr(strings[i],"G4(0 K)=")) ||
                 (NULL != strstr(strings[i],"CBS-QB3 (0 K)=")))
        {
            status = gau_comp_meth_read_line(strings[i],&comp_0K,&comp_energy);
            if (bVerbose)
                printf("na gau_(): comp_0K %f comp_energy %f\n", comp_0K,comp_energy);
        }
        else if ((NULL != strstr(strings[i],"G2 Enthalpy=")) ||
                 (NULL != strstr(strings[i],"G3 Enthalpy=")) ||
                 (NULL != strstr(strings[i],"G4 Enthalpy=")) ||
                 (NULL != strstr(strings[i],"CBS-QB3 Enthalpy=")))
        {
            status = gau_comp_meth_read_line(strings[i],&comp_enthalpy,&comp_free_energy);
            if (bVerbose)
                printf("na gau_(): comp_enthalpy %f comp_free_energy %f\n", 
                       comp_enthalpy, comp_free_energy);
        }
        else if (NULL != strstr(strings[i],"GINC")) 
        {
            j = 0;
            while ((i+j<nstrings) && (strlen(strings[i+j]) > 0))
                j++;
            snew(ginc,80*(j+1));
            for(k=i; (k<i+j); k++) 
            {
                trim(strings[k]);
                strncat(ginc,strings[k],80);
            }
            i = k-1;
            ptr = split('\\',ginc);
            for(j=0; (NULL != ptr[j]); j++) 
            {
                if (NULL != strstr(ptr[j],"HF=")) 
                {
                    hfener = strdup(ptr[j]+3);
                }
                if (NULL != strstr(ptr[j],"MP2=")) 
                {
                    mp2ener = strdup(ptr[j]+4);
                }
                if (NULL != strstr(ptr[j],"G2=")) 
                {
                    g2ener = strdup(ptr[j]+3);
                }
                if (NULL != strstr(ptr[j],"G3=")) 
                {
                    g3ener = strdup(ptr[j]+3);
                }
                if (NULL != strstr(ptr[j],"G4=")) 
                {
                    g4ener = strdup(ptr[j]+3);
                }
                if (NULL != strstr(ptr[j],"CBSQB3=")) 
                {
                    cbsener = strdup(ptr[j]+7);
                }
            }
        }
        else if (NULL != strstr(strings[i],"#P")) 
        {
            if (NULL == method) 
            {
                /* The first method is the true method */
                ptr = split(' ',strings[i]);
                if (NULL != ptr[1]) {
                    qtr = split('/',ptr[1]);
                    if (NULL != qtr[0]) 
                    {
                        method = strdup(qtr[0]);
                    }
                }
            }
        }
    }
  
    /* Add the calculation */
    if ((NULL == basis) && (NULL != basisset))
        basis = basisset;
    if ((calcref == NOTSET) && (NULL != program) && (NULL != method) && (NULL != basis))
    {
        gmx_molprop_add_calculation(mpt,program,method,basis,reference,conformation,fn,&calcref);
        if (bPolar)
            gmx_molprop_add_polar(mpt,calcref,"elec","A^3",polar[XX][XX],polar[YY][YY],polar[ZZ][ZZ],
                                  trace(polar)/3.0,0);
        if (bDipole)
            gmx_molprop_add_dipole(mpt,calcref,"elec","Debye",
                                   dipole[XX],dipole[YY],dipole[ZZ],
                                   norm(dipole),0);
        if (bQuad)
            gmx_molprop_add_quadrupole(mpt,calcref,"elec","Buckingham",
                                       quad[XX][XX],quad[YY][YY],quad[ZZ][ZZ],
                                       quad[XX][YY],quad[XX][ZZ],quad[YY][ZZ]);
	
	
        for(i=0; (i<nstrings); i++) 
        {
            bAtomicCenter = (NULL != strstr(strings[i],"Atomic Center"));
          
            if (NULL != strstr(strings[i],"fitting atomic charges"))
            {
                if (1 == sscanf(strings[i],"%d",&kk))
                {
                    if (nfitpoints == -1)
                    {
                        nfitpoints = kk;
                    }
                    else
                    {
                        if (!bWarnESP)
                        {
                            /*   gmx_resp_warning(fn,i+1);*/
                            bWarnESP = TRUE;
                        }
                        if (kk != nfitpoints)
                            fprintf(stderr,"nfitpoints was %d, now is %d\n",
                                    nfitpoints,kk);
                        nfitpoints = kk;
                    }
                }
            }
            else if (NULL != strstr(strings[i],"Stoichiometry"))
            {
                if (1 == sscanf(strings[i],"%*s%s",sbuf))
                {
                    gmx_molprop_set_formula(mpt,sbuf);
                    if (NULL == molnm)
                        gmx_molprop_set_molname(mpt,sbuf);
                    else
                        gmx_molprop_set_molname(mpt,molnm);
                    if (NULL == iupac)
                        gmx_molprop_set_iupac(mpt,sbuf);
                    else
                        gmx_molprop_set_iupac(mpt,iupac);
                }
            }
            else if (NULL != strstr(strings[i],"Coordinates (Angstroms)"))
            {
                atomicnumber = i;
            }
            else if ((atomicnumber >= 0) && (i >= atomicnumber+3))
            {
                if (6 == sscanf(strings[i],"%d%d%d%lf%lf%lf",
                                &k,&anumber,&kk,&x,&y,&z))
                {
                    if (natom+1 != k)
                    {
                        if (!bWarnESP)
                        {
                            /* gmx_resp_warning(fn,i+1);*/
                            bWarnESP = TRUE;
                        }
                    }
                    else
                    {
                        natom++;
                        atomname = gmx_atomprop_element(aps,anumber);
                        gmx_molprop_calc_add_atom(mpt,calcref,atomname,atomname,natom,&atomref);
                        if (TRUE == gmx_atomprop_query(aps,epropMass,"",atomname,&mm)) {
                            mass += mm;
                        }
                        gmx_molprop_calc_set_atomcoords(mpt,calcref,atomref,"pm",
                                                        100*x,100*y,100*z);
                      
                    }
                }
                else
                    atomicnumber = -1;
            }
            else if (NULL != strstr(strings[i],"Charge ="))
            {
                if (1 != sscanf(strings[i],"%*s%*s%d",&charge))
                    fprintf(stderr,"Can not read the charge on line %d of file %s\n",i+1,fn);
            }
            else if (bAtomicCenter || (NULL != strstr(strings[i],"ESP Fit Center")))
            {
                if ((bAtomicCenter  && (4 != sscanf(strings[i],"%*s%*s%d%*s%*s%lf%lf%lf",&k,&x,&y,&z))) ||
                    (!bAtomicCenter && (4 != sscanf(strings[i],"%*s%*s%*s%d%*s%*s%lf%lf%lf",&k,&x,&y,&z))))
                    fprintf(stderr,"Warning something fishy on line %d of file %s\n",i+1,fn);
                else
                {
                    if (bAtomicCenter)
                    {
                        if (NULL != debug)
                            fprintf(debug,"Coordinates for atom %d found on line %d %8.3f  %8.3f  %8.3f\n",k,i,x,y,z);
                    }
                    if (k > nesp)
                    {
                        nesp++;
                        srenew(espv,nesp);
                    }
                    espv[k-1].x = 100*x;
                    espv[k-1].y = 100*y;
                    espv[k-1].z = 100*z;
                    if (NULL != debug)
                        fprintf(debug,"Coordinates for fit %d found on line %d\n",k,i);
                }
            }
            else if (NULL != strstr(strings[i],"Electrostatic Properties (Atomic Units)"))
            {
                if ((NOTSET != nelprop) && (!bWarnESP))
                {
                    /*gmx_resp_warning(fn,i+1);*/
                    bWarnESP = TRUE;
                }
                nelprop = i;
            }
            else if ((NOTSET != nelprop) && (i >= nelprop+6) && (i < nelprop+6+natom+nfitpoints))
            {
                if (2 != sscanf(strings[i],"%d%*s%lf",&k,&V))
                    fprintf(stderr,"Warning something fishy on line %d of file %s\n",i+1,fn);
                if ((k-1) != (i-nelprop-6))
                    fprintf(stderr,"Warning something fishy with fit center number on line %d of file %s\n",i+1,fn);
                if (k-1 >= nesp)
                    fprintf(stderr,"More potential points (%d) than fit centers (%d). Check line %d of file %s\n",
                            k,nesp,i+1,fn);
                espv[k-1].V = V;
                if (NULL != debug)
                    fprintf(debug,"Potential %d found on line %d\n",k,i);
            }
            sfree(strings[i]);
        }
        if (debug)
            fprintf(debug,"Found %d atoms, %d esp points, and %d fitpoints\n",natom,nesp,nfitpoints);
        sfree(strings);
      
        if ((charge == NOTSET) || (natom == 0)) 
            gmx_fatal(FARGS,"Error reading Gaussian log file. Charge = %d, natom = %d",
                      charge,natom);
        gmx_molprop_set_charge(mpt,charge);
        gmx_molprop_set_mass(mpt,mass);
      
        if ((maxpot > 0) && (maxpot < nesp)) {
            qsort(espv,nesp,sizeof(espv[0]),esp_comp);
            deltai = nesp/maxpot;
        }
        else {
            maxpot = nesp;
            deltai = 1;
        }
        ii = 0;
        for(i=0; (i<maxpot); i++) 
        {
            /* Convert to integer */
            ii0 = ii;
            gmx_molprop_add_potential(mpt,calcref,"pm","Hartree/e",i,
                                      espv[ii0].x,espv[ii0].y,espv[ii0].z,
                                      espv[ii0].V);
            ii+=deltai;
        }
        sfree(espv);
      
        /* Generate atomic composition, needed for energetics */
        generate_composition(1,&mpt,pd,aps,TRUE,th_toler,ph_toler);

        /* Add energies */
        if ((NULL != g4ener) || (NULL != g2ener) || (NULL != g3ener) || (NULL != cbsener))
        {
            if (NULL != g4ener) 
            {
                mymeth = (char *)"G4";
                myener = atof(g4ener);
            }
            else if (NULL != g3ener)
            {
                mymeth = (char *)"G3";
                myener = atof(g3ener);
            }
            else if (NULL != g2ener)
            {
                mymeth = (char *)"G2";
                myener = atof(g2ener);
            }
            else 
            {
                mymeth = (char *)"CBS-QB3";
                myener = atof(cbsener);
            }
            if (bVerbose)
                printf("myener %f\n",myener);
            if (bEtherm && bEzpe && bTemp)
            {
                if (0 == gmx_molprop_add_dhform(mpt,calcref,gaps,mymeth,temp,myener,ezpe,etherm,bVerbose))
                {
                    fprintf(stderr,"No support for atomic energies in %s, method %s\n",
                            gmx_molprop_get_molname(mpt),mymeth);
                }
            }
        }		
        else if (NULL != mp2ener)
        {
            ee = convert2gmx(atof(mp2ener),eg2c_Hartree);
            gmx_molprop_add_energy(mpt,calcref,"MP2",unit2string(eg2c_kJ_mole),ee,0);
        }
        else if (NULL != hfener)
        {
            ee = convert2gmx(atof(hfener),eg2c_Hartree);
            gmx_molprop_add_energy(mpt,calcref,"HF",unit2string(eg2c_kJ_mole),ee,0);
        }
            
        return mpt;
    }
    else {
        fprintf(stderr,"Error reading %s, program = '%s' basis = '%s' method = '%s'\n",
                fn,program,basis,method);
        for(i=0; (i<nstrings); i++) {
            sfree(strings[i]);
        }
        sfree(strings);
        return NULL;
    }
}

gmx_molprop_t gmx_molprop_read_gauss(const char *g98,gmx_bool bBabel,
                                     gmx_atomprop_t aps,gmx_poldata_t pd,
                                     char *molnm,char *iupac,char *conf,
                                     char *basis,gau_atomprop_t gaps,
                                     real th_toler,real ph_toler,
                                     int maxpot,gmx_bool bVerbose)
{
    gmx_molprop_t mp = NULL;
  
#ifdef HAVE_LIBOPENBABEL2
    if (bBabel)
        mp = gmx_molprop_read_babel(g98,aps,pd,molnm,iupac,conf,basis,
                                    th_toler,ph_toler,maxpot,bVerbose);
    else
        mp = gmx_molprop_read_log(g98,aps,pd,molnm,iupac,conf,basis,gaps,
                                  th_toler,ph_toler,maxpot,bVerbose);
#else
    mp = gmx_molprop_read_log(g98,aps,pd,molnm,iupac,conf,basis,gaps,
                              th_toler,ph_toler,maxpot,bVerbose);
#endif

    return mp;
}
