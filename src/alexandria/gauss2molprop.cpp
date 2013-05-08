/*
 * $Id: gauss2molprop.cpp,v 1.26 2009/05/20 10:48:03 spoel Exp $
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
#include "poldata.h"
#include "poldata_xml.h"
#include "molprop.hpp"
#include "molprop_util.hpp"
#include "molprop_xml.hpp"
#include "gauss_io.hpp"

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
  static gmx_bool bVerbose = FALSE;
  static char *molnm=NULL,*iupac=NULL,*conf=(char *)"minimum",*basis=NULL;
  static const char *forcefield="GAFF";
  static int  maxpot=0;
  static gmx_bool compress=FALSE;
  static gmx_bool bBabel=TRUE;
  t_pargs pa[] = {
    { "-v",      FALSE, etBOOL, {&bVerbose},
      "Generate verbose terminal output." },
    { "-compress", FALSE, etBOOL, {&compress},
      "Compress output XML files" },
#ifdef HAVE_LIBOPENBABEL2
    { "-babel", FALSE, etBOOL, {&bBabel},
      "Use the OpenBabel engine to process gaussian input files" },
#endif
    { "-molnm", FALSE, etSTR, {&molnm},
      "Name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
    { "-iupac", FALSE, etSTR, {&iupac},
      "IUPAC name of the molecule in *all* input files. Do not use if you have different molecules in the input files." },
    { "-conf",  FALSE, etSTR, {&conf},
      "Conformation of the molecule" },
    { "-basis",  FALSE, etSTR, {&basis},
      "Basis-set used in this calculation for those case where it is difficult to extract from a Gaussian file" },
    { "-ff", FALSE, etSTR, {&forcefield},
      "Force field for basic atom typing available in OpenBabel" },
    { "-maxpot", FALSE, etINT, {&maxpot},
      "Max number of potential points to add to the molprop file. If 0 all points are registered, else a selection of points evenly spread over the range of values is taken" }
  };
  output_env_t oenv;
  gmx_atomprop_t aps;
  gmx_poldata_t  pd;
  std::vector<alexandria::MolProp> mp;
  alexandria::GaussAtomProp gap;
  char **fns=NULL;
  int i,nfn;
  
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL,&oenv);
  
  /* Read standard atom properties */
  aps = gmx_atomprop_init();
  
  /* Read polarization stuff */
  if ((pd = gmx_poldata_read(NULL,aps)) == NULL)
    gmx_fatal(FARGS,"Can not read the force field information. File missing or incorrect.");

  nfn = ftp2fns(&fns,efLOG,NFILE,fnm);
  for(i=0; (i<nfn); i++) 
  {
      alexandria::MolProp mmm;
      
      ReadGauss(fns[i],mmm,gap,bBabel,aps,pd,molnm,iupac,conf,basis,
                maxpot,bVerbose,gmx_poldata_get_force_field(pd));
      mp.push_back(mmm);
  }
  
  printf("Succesfully read %d molprops from %d Gaussian files.\n",
         (int)mp.size(),nfn);
  MolPropSort(mp,MPSA_MOLNAME,NULL,NULL);
  merge_doubles(mp,NULL,TRUE);
  if (mp.size() > 0)
  {
      MolPropWrite(opt2fn("-o",NFILE,fnm),mp,(int)compress);
  }
      
  thanx(stderr);
  
  return 0;
}
