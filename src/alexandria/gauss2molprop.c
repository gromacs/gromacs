/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * $Id: gauss2molprop.c,v 1.26 2009/05/20 10:48:03 spoel Exp $
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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "typedefs.h"
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "filenm.h"
#include "pbc.h"
#include "pdbio.h"
#include "gpp_atomtype.h"
#include "atomprop.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "poldata.h"
#include "poldata_xml.h"

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

typedef struct gau_atomprop *gau_atomprop_t;

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
gau_atomprop_t read_gauss_data(const char *fn)
{
    FILE *fp;
    char **strings=NULL,**ptr;
    int nstrings,i,j=0;
    gau_atomprop_t gaps;

    snew(gaps,1);
    nstrings = get_lib_file(fn,&strings);
    
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

static int gau_atomprop_get_value(gau_atomprop_t gaps,char *element,char *method,char *desc,double temp,
                                  double *value)
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
            (1 == gau_atomprop_get_value(gaps,atomname,"exp","DHf(0K)",0,&vve)) &&
            (1 == gau_atomprop_get_value(gaps,atomname,"exp","H(0K)-H(298.15K)",298.15,&vvdh)))
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
    
    ee = convert2gmx(e0Hartree,eg2cHartree);
    gmx_molprop_add_energy(mpt,calcref,"DHf(0K)","kJ/mol",ee,0);
    if (bVerbose)
        printf("natom %d e0Hartree %f eTHartree %f eg34 %f ve %f vm %f vdh %f ezpe %f etherm %f \n",
               ntot,e0Hartree,eTHartree,eg34,ve,vm,vdh,ezpe,etherm);
           
    ee = convert2gmx(eTHartree,eg2cHartree);
    gmx_molprop_add_energy(mpt,calcref,"DHf(298.15K)","kJ/mol",ee,0);
    
    return 1;
}

/* Read a line from a G03/G09 composite method (G3, G4, etc) record */
int gau_comp_meth_read_line(char *line,real *temp,real *pres)
{
    char **ptr1,**ptr2;
    gmx_bool bThermResults=FALSE;
    
    ptr1 = split('=', line);
    if (NULL != ptr1[1]) 
    {
        ptr2 = split(' ',ptr1[1]);
        if (NULL != ptr2[0] && NULL != ptr1[2]) {
            *temp = atof(ptr2[0]);
            *pres = atof(ptr1[2]);
        }
    }

    return TRUE;
}

typedef struct {
    rvec esp;
    real V;
} t_espv;

static int espv_comp(const void *a,const void *b)
{
    t_espv *va = (t_espv *)a;
    t_espv *vb = (t_espv *)b;
    real dv = va->V - vb->V;
    
    if (dv < 0)
        return -1;
    else if (dv > 0)
        return 1;
    return 0;
}

gmx_molprop_t gmx_molprop_read_log(gmx_atomprop_t aps,gmx_poldata_t pd,
                                   const char *fn,char *molnm,char *iupac,
                                   gau_atomprop_t gaps,
                                   real th_toler,real ph_toler,
                                   int maxpot,gmx_bool bVerbose)
{
  /* Read a gaussian log file */
  char **strings=NULL;
  char sbuf[STRLEN];
  int  nstrings,atomiccenter=-1,atomicnumber=-1;
  double ee,x,y,z,V,myener,mass=0;
  int i,j,k,kk,zz,anumber,natom,nesp,nelprop,charge,nfitpoints=-1;
  int calcref=NOTSET,atomref=NOTSET,len,iMP;
  gmx_bool bWarnESP=FALSE;
  gmx_bool bAtomicCenter;
  gmx_molprop_t mpt;
  real mm;
  char *atomname,*ginc,*hfener,*mp2ener,*g2ener,*g3ener,*g4ener,*cbsener;
  char *conformation = "unknown";
  char *reference = "This Work";
  char *program=NULL,*method=NULL,*basis=NULL;
  char **ptr,**qtr,*mymeth;
  rvec xatom;
  t_espv *espv=NULL;
  char **ptr2;
  gmx_bool bThermResults=FALSE;
  real temp,pres,ezpe,ezpe2,etherm,etherm2,comp_0K,comp_energy,comp_enthalpy,comp_free_energy,atom_ener,temp_corr,ii,deltai;
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
      else if ((NULL != strstr(strings[i],"Standard basis:")) || 
               (NULL != strstr(strings[i],"Basis read from chk"))) 
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

          bThermResults = TRUE;
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
          if (NULL != ptr[2]) 
          {
              ezpe = atof(strdup(ptr[2]));
          }
          if (bVerbose)
              printf("na gau_(): ezpe %f \n", ezpe);

          bThermResults = TRUE;
      }
      else if (NULL != strstr(strings[i],"Thermal correction to Enthalpy=")) 
      {
          
          ptr = split(' ',strings[i]);
          if (NULL != ptr[4]) 
          {
              etherm = atof(strdup(ptr[4]));
          }

          if (bVerbose)
              printf("na gau_(): etherm %f \n", etherm);

          bThermResults = TRUE;
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
  if ((calcref == NOTSET) && (NULL != program) && (NULL != method) && (NULL != basis))
  {
      gmx_molprop_add_calculation(mpt,program,method,basis,reference,conformation,&calcref);
    
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
                      gmx_molprop_calc_add_atom(mpt,calcref,atomname,natom,&atomref);
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
                      xatom[XX] = 100*x;
                      xatom[YY] = 100*y;
                      xatom[ZZ] = 100*z;
                      if (NULL != debug)
                          fprintf(debug,"Coordinates for atom %d found on line %d %8.3f  %8.3f  %8.3f\n",k,i,x,y,z);
                  }
                  if (k > nesp)
                  {
                      nesp++;
                      srenew(espv,nesp);
                  }
                  espv[k-1].esp[XX] = 100*x;
                  espv[k-1].esp[YY] = 100*y;
                  espv[k-1].esp[ZZ] = 100*z;
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
          gmx_fatal(FARGS,"Error reading Gaussian log file.");
      gmx_molprop_set_charge(mpt,charge);
      gmx_molprop_set_mass(mpt,mass);
      
      if ((maxpot > 0) && (maxpot < nesp)) {
          qsort(espv,nesp,sizeof(espv[0]),espv_comp);
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
                                    espv[ii0].esp[XX],espv[ii0].esp[YY],espv[ii0].esp[ZZ],
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
              mymeth = "G4";
              myener = atof(g4ener);
          }
          else if (NULL != g3ener)
          {
              mymeth = "G3";
              myener = atof(g3ener);
          }
          else if (NULL != g2ener)
          {
              mymeth = "G2";
              myener = atof(g2ener);
          }
          else 
          {
              mymeth = "CBS-QB3";
              myener = atof(cbsener);
          }
          if (bVerbose)
              printf("myener %f\n",myener);
          if (0 == gmx_molprop_add_dhform(mpt,calcref,gaps,mymeth,temp,myener,ezpe,etherm,bVerbose))
          {
              fprintf(stderr,"No support for atomic energies in %s, method %s\n",
                      gmx_molprop_get_molname(mpt),mymeth);
          }
      }		
      else if (NULL != mp2ener)
      {
          ee = convert2gmx(atof(mp2ener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"MP2","kJ/mol",ee,0);
      }
      else if (NULL != hfener)
      {
          ee = convert2gmx(atof(hfener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"HF","kJ/mol",ee,0);
      }
            
      return mpt;
  }
  else {
    fprintf(stderr,"Error reading %s\n",fn);
    for(i=0; (i<nstrings); i++) {
      sfree(strings[i]);
    }
    sfree(strings);
    return NULL;
  }
}

int main(int argc, char *argv[])
{
  static const char *desc[] = {
    "gauss2molprop reads a series of Gaussian output files, and collects",
    "useful information, and saves it to molprop file."
  };
    
  t_filenm fnm[] = {
    { efLOG, "-g03",  "gauss",  ffRDMULT },
    { efDAT, "-f",    "atomization_energies.dat",  ffREAD },
    { efDAT, "-o",    "molprop", ffWRITE }
  };
#define NFILE asize(fnm)
  static gmx_bool bVerbose = FALSE;
  static char *molnm=NULL,*iupac=NULL;
  static real th_toler=170,ph_toler=5;
  static int  maxpot=0;
  static gmx_bool compress=FALSE;
  t_pargs pa[] = {
    { "-v",      FALSE, etBOOL, {&bVerbose},
      "Generate verbose terminal output." },
    { "-th_toler", FALSE, etREAL, {&th_toler},
      "HIDDENMinimum angle to be considered a linear A-B-C bond" },
    { "-ph_toler", FALSE, etREAL, {&ph_toler},
      "HIDDENMaximum angle to be considered a planar A-B-C/B-C-D torsion" },
    { "-compress", FALSE, etBOOL, {&compress},
      "Compress output XML files" },
    { "-molnm", FALSE, etSTR, {&molnm},
      "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
    { "-iupac", FALSE, etSTR, {&iupac},
      "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
    { "-maxpot", FALSE, etINT, {&maxpot},
      "Max number of potential points to add to the molprop file. If 0 all points are registered, else a selection of points evenly spread over the range of values is taken" }
  };
  output_env_t oenv;
  gmx_atomprop_t aps;
  gmx_poldata_t  pd;
  gmx_molprop_t mp,*mps=NULL;
  gau_atomprop_t gp;
  char **fns=NULL;
  int i,nmp,nfn;
  FILE *fp;
  gau_atomprop_t gaps;

  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL,&oenv);
  
  /* Read standard atom properties */
  aps = gmx_atomprop_init();
  
  /* Read polarization stuff */
  if ((pd = gmx_poldata_read(NULL,aps)) == NULL)
    gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");
    
  gaps = read_gauss_data(opt2fn_null("-f",NFILE,fnm));

  nfn = ftp2fns(&fns,efLOG,NFILE,fnm);
  nmp = 0;
  for(i=0; (i<nfn); i++) 
  {
      mp = gmx_molprop_read_log(aps,pd,fns[i],molnm,iupac,gaps,
                                th_toler,ph_toler,maxpot,bVerbose);
      if (NULL != mp) 
      {
          srenew(mps,++nmp);
          mps[nmp-1] = mp;
      }
  }
  printf("Succesfully read %d molprops from %d Gaussian files.\n",nmp,nfn);
  if (nmp > 0)
  {
      gmx_molprops_write(opt2fn("-o",NFILE,fnm),nmp,mps,(int)compress);
  }
      
  thanx(stderr);
  
  return 0;
}
