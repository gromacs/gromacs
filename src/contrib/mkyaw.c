/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gromacs Runs On Most of All Computer Systems
 */
static char *SRCID_mkyaw_c = "$Id$";
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "assert.h"
#include "pdbio.h"
#include "confio.h"
#include "symtab.h"
#include "smalloc.h"
#include "symtab.h"
#include "macros.h"
#include "copyrite.h"
#include "statutil.h"
#include "string2.h"
#include "strdb.h"
#include "rdgroup.h"
#include "vec.h"
#include "typedefs.h"
#include "gbutil.h"
#include "strdb.h"
#include "rdgroup.h"
#include "physics.h"
#include "atomprop.h"

void copy_atom(t_symtab *tab,t_atoms *a1,int i1,t_atoms *a2,int i2,
	       rvec xin[],rvec xout[],rvec vin[],rvec vout[])
{
  a2->atom[i2]     = a1->atom[i1];
  a2->atomname[i2] = put_symtab(tab,*a1->atomname[i1]);
  a2->resname[a2->atom[i2].resnr] =
    put_symtab(tab,*a1->resname[a1->atom[i1].resnr]);
  copy_rvec(xin[i1],xout[i2]);
  copy_rvec(vin[i1],vout[i2]);
}

int main(int argc, char *argv[])
{
  t_symtab tab;
  static char *desc[] = {
    "mkyaw adds to an existing conf file for every OW atom an DW and SW",
    "after the hydrogens (or the inverse with the -back option)."
  };
  static bool bBack=FALSE,bDW=TRUE;
  t_pargs pa[] = {
    { "-back",   FALSE, etBOOL, &bBack, 
      "Remove SW and DW" },
    { "-dw",     FALSE, etBOOL, &bDW,
      "Use both dummy and shell" }
  };
#define NPA asize(pa)
  t_filenm fnm[] = {
    { efSTX, "-f", NULL, ffREAD },
    { efSTO, "-o", NULL, ffWRITE }
  };
#define NFILE asize(fnm)
  int  i,iout,now,natom;
  rvec *xin,*vin,*xout,*vout;
  matrix box;
  t_atoms atoms,aout;
  char *infile,*outfile,title[256];
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,NPA,pa,
		    asize(desc),desc,0,NULL);

  infile  = ftp2fn(efSTX,NFILE,fnm);
  outfile = ftp2fn(efSTO,NFILE,fnm);
  
  get_stx_coordnum(infile,&natom);
  init_t_atoms(&atoms,natom,TRUE);
  snew(xin,natom);
  snew(vin,natom);
  read_stx_conf(infile,title,&atoms,xin,vin,box);
  printf("Read %d atoms\n",atoms.nr); 
  open_symtab(&tab);
  if (!bBack) {
    now = 0;
    for(i=0; (i<natom-2); ) {
      if ((strstr(*atoms.atomname[i],"OW")   != NULL) &&
	  (strstr(*atoms.atomname[i+1],"HW") != NULL) &&
	  (strstr(*atoms.atomname[i+2],"HW") != NULL)) {
	now++;
	i+=3;
      }
      else
	i++;
    }
    fprintf(stderr,"There are %d water molecules\n",now);
    init_t_atoms(&aout,natom+2*now,TRUE);
    snew(xout,natom+2*now);
    snew(vout,natom+2*now);
    for(i=iout=0; (i<natom); i++) {
      copy_atom(&tab,&atoms,i,&aout,iout,xin,xout,vin,vout);
      iout++;
      if (i >= 2) {
	if (strstr(*atoms.atomname[i-2],"OW") != NULL) {
	  if (bDW) {
	    copy_atom(&tab,&atoms,i-2,&aout,iout,xin,xout,vin,vout);
	    aout.atomname[iout] = put_symtab(&tab,"DW");
	    iout++;
	  }
	  copy_atom(&tab,&atoms,i-2,&aout,iout,xin,xout,vin,vout);
	  aout.atomname[iout] = put_symtab(&tab,"SW");
	  iout++;
	}
      }
    }
    aout.nr = iout;
    close_symtab(&tab);
    fprintf(stderr,"iout = %d\n",iout);
    write_sto_conf(outfile,"Gravity Sucks",&aout,xout,vout,box); 
  }

  thanx(stderr);
  
  return 0;
}
