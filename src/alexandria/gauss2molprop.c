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

gmx_molprop_t gmx_molprop_read_log(gmx_atomprop_t aps,gmx_poldata_t pd,
                                   const char *fn,char *molnm,char *iupac)
{
  /* Read a gaussian log file */
  char **strings=NULL;
  char sbuf[STRLEN];
  int  nstrings,atomiccenter=-1,atomicnumber=-1;
  double ee,x,y,z,V,mass=0;
  int i,j,k,kk,zz,anumber,natom,nesp,nelprop,charge,nfitpoints=-1;
  int calcref=NOTSET,atomref=NOTSET,len,iMP;
  gmx_bool bWarnESP=FALSE;
  gmx_bool bAtomicCenter;
  gmx_molprop_t mpt;
  real mm;
  char *atomname,*ginc,*hfener,*mp2ener,*g3ener,*g4ener,*cbsener;
  char *conformation = "unknown";
  char *reference = "This Work";
  char *program=NULL,*method=NULL,*basis=NULL;
  char **ptr,**qtr;
  rvec xatom;
  rvec *esp=NULL;
  real *pot=NULL;
  
  nstrings = get_file(fn,&strings);
    
  /* Create new calculation */ 
  mpt = gmx_molprop_init();

  natom   = 0;
  nesp    = 0;
  nelprop = NOTSET;
  charge  = NOTSET;
  hfener  = NULL;
  mp2ener = NULL;
  g3ener  = NULL;
  g4ener  = NULL;
  cbsener = NULL;
  
  /* First loop over strings to deduct basic stuff */
  for(i=0; (i<nstrings); i++) 
    {
      if (NULL != strstr(strings[i],"Revision")) {
        program = strdup(strings[i]);
        trim(program);
        len = strlen(program);
        if ((len > 0) && (program[len-1] == ',')) {
          program[len-1] = '\0';
        }
      }
      else if ((NULL != strstr(strings[i],"Standard basis:")) || 
               (NULL != strstr(strings[i],"Basis read from chk"))) {
        ptr = split(' ',strings[i]);
        if (NULL != ptr[2]) 
          {
            basis = strdup(ptr[2]);
          }
      }
      else if (NULL != strstr(strings[i],"GINC")) {
          j = 0;
          while ((j<nstrings) && (strlen(strings[i+j]) > 0))
              j++;
          snew(ginc,80*(j+1));
          for(k=i; (k<i+j); k++) {
              trim(strings[k]);
              strncat(ginc,strings[k],80);
          }
          i = k-1;
          ptr = split('\\',ginc);
          for(j=0; (NULL != ptr[j]); j++) {
              if (NULL != strstr(ptr[j],"HF=")) {
                  hfener = strdup(ptr[j]+3);
              }
              if (NULL != strstr(ptr[j],"MP2=")) {
                  mp2ener = strdup(ptr[j]+4);
              }
              if (NULL != strstr(ptr[j],"G3=")) {
                  g3ener = strdup(ptr[j]+3);
              }
              if (NULL != strstr(ptr[j],"G4=")) {
                  g4ener = strdup(ptr[j]+3);
              }
              if (NULL != strstr(ptr[j],"CBSQB3=")) {
                  cbsener = strdup(ptr[j]+7);
              }
          }
      }
      else if (NULL != strstr(strings[i],"#P")) {
        ptr = split(' ',strings[i]);
        if (NULL != ptr[1]) {
          qtr = split('/',ptr[1]);
          if (NULL != qtr[0]) {
            method = strdup(qtr[0]);
          }
        }
      }
    }
    
  if ((calcref == NOTSET) && (NULL != program) && (NULL != method) && (NULL != basis))
    {
      gmx_molprop_add_calculation(mpt,program,method,basis,reference,conformation,&calcref);
      if (NULL != hfener)
      {
          ee = convert2gmx(atof(hfener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"HF","kJ/mol",ee,0);
      }
      if (NULL != mp2ener)
      {
          ee = convert2gmx(atof(mp2ener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"MP2","kJ/mol",ee,0);
      }		
      if (NULL != g3ener)
      {
          ee = convert2gmx(atof(g3ener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"G3","kJ/mol",ee,0);
      }		
      if (NULL != g4ener)
      {
          ee = convert2gmx(atof(g4ener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"G4","kJ/mol",ee,0);
      }		
      if (NULL != cbsener)
      {
          ee = convert2gmx(atof(cbsener),eg2cHartree);
          gmx_molprop_add_energy(mpt,calcref,"CBS-BQ3","kJ/mol",ee,0);
      }		
    }
  if (calcref != NOTSET) 
    {
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
                      srenew(esp,nesp);
                    }
                  esp[k-1][XX] = 100*x;
                  esp[k-1][YY] = 100*y;
                  esp[k-1][ZZ] = 100*z;
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
              snew(pot,nesp);
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
              pot[k-1] = V;
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

      for(i=0; (i<nesp); i++) {
        gmx_molprop_add_potential(mpt,calcref,"pm","Hartree/e",i,
                                  esp[i][XX],esp[i][YY],esp[i][ZZ],pot[i]);
      }
      sfree(pot);
      sfree(esp);
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
    { efDAT, "-o",    "molprop", ffWRITE }
  };
#define NFILE asize(fnm)
  static gmx_bool bVerbose = TRUE;
  static char *molnm=NULL,*iupac=NULL;
  static real th_toler=170,ph_toler=5;
  static gmx_bool compress=FALSE;
  t_pargs pa[] = {
    { "-v",      FALSE, etBOOL, {&bVerbose},
      "Generate verbose output in the top file and on terminal." },
    { "-th_toler", FALSE, etREAL, {&th_toler},
      "Minimum angle to be considered a linear A-B-C bond" },
    { "-ph_toler", FALSE, etREAL, {&ph_toler},
      "Maximum angle to be considered a planar A-B-C/B-C-D torsion" },
    { "-compress", FALSE, etBOOL, {&compress},
      "Compress output XML files" },
    { "-molnm", FALSE, etSTR, {&molnm},
      "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
    { "-iupac", FALSE, etSTR, {&iupac},
      "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." }
  };
  output_env_t oenv;
  gmx_atomprop_t aps;
  gmx_poldata_t  pd;
  gmx_molprop_t mp,*mps=NULL;
  char **fns=NULL;
  int i,nmp,nfn;
        
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL,&oenv);
  
  /* Read standard atom properties */
  aps = gmx_atomprop_init();
  
  /* Read polarization stuff */
  if ((pd = gmx_poldata_read(NULL,aps)) == NULL)
    gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");
    
  nfn = ftp2fns(&fns,efLOG,NFILE,fnm);
  nmp = 0;
  for(i=0; (i<nfn); i++) 
  {
      mp = gmx_molprop_read_log(aps,pd,fns[i],molnm,iupac);
      if (NULL != mp) 
      {
          srenew(mps,++nmp);
          mps[nmp-1] = mp;
      }
  }
  generate_composition(nmp,mps,pd,aps,TRUE,th_toler,ph_toler);
  printf("Succesfully read %d molprops from %d Gaussian files.\n",nmp,nfn);
  if (nmp > 0)
  {
      gmx_molprops_write(ftp2fn(efDAT,NFILE,fnm),nmp,mps,(int)compress);
  }
      
  thanx(stderr);
  
  return 0;
}
