/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "smalloc.h"
#include "statutil.h"
#include "confio.h"
#include "genhydro.h"
#include "tpxio.h"
#include "index.h"
#include "vec.h"
#include "hackblock.h"

int main (int argc,char *argv[])
{
  const char *desc[] = {
    "[TT]protonate[tt] reads (a) conformation(s) and adds all missing",
    "hydrogens as defined in [TT]ffgmx2.hdb[tt]. If only [TT]-s[tt] is",
    "specified, this conformation will be protonated, if also [TT]-f[tt]",
    "is specified, the conformation(s) will be read from this file",
    "which can be either a single conformation or a trajectory.",
    "[PAR]",
    "If a pdb file is supplied, residue names might not correspond to",
    "to the GROMACS naming conventions, in which case these residues will",
    "probably not be properly protonated.",
    "[PAR]",
    "If an index file is specified, please note that the atom numbers",
    "should correspond to the [BB]protonated[bb] state."
  };
  
  char        title[STRLEN+1];  
  const char  *infile;
  char        *grpnm;
  t_topology  top;
  int         ePBC;
  t_atoms     *atoms,*iatoms;
  t_protonate protdata;
  atom_id     *index;
  t_trxstatus *status;
  t_trxstatus *out;
  t_trxframe  fr,frout;
  rvec        *x,*ix;
  int         nidx,natoms,natoms_out;
  matrix      box;
  int         i,frame,resind;
  gmx_bool        bReadMultiple;
  output_env_t oenv;
  
  t_filenm fnm[] = {
    { efTPS, NULL, NULL,         ffREAD  },
    { efTRX, "-f", NULL,         ffOPTRD },
    { efNDX, NULL, NULL,         ffOPTRD },
    { efTRO, "-o", "protonated", ffWRITE }
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME,
		    NFILE,fnm,0,NULL,asize(desc),desc,0,NULL,&oenv);
  
  infile=opt2fn("-s",NFILE,fnm);
  read_tps_conf(infile,title,&top,&ePBC,&x,NULL,box,FALSE);
  atoms=&(top.atoms);
  printf("Select group to process:\n");
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nidx,&index,&grpnm);
  bReadMultiple = opt2bSet("-f",NFILE,fnm);
  if (bReadMultiple) {
    infile = opt2fn("-f",NFILE,fnm);
    if ( !read_first_frame(oenv,&status, infile, &fr, TRX_NEED_X ) )
      gmx_fatal(FARGS,"cannot read coordinate file %s",infile);
    natoms = fr.natoms;
  } else {
    clear_trxframe(&fr,TRUE);
    fr.natoms = atoms->nr;
    fr.bTitle = TRUE;
    fr.title  = title;
    fr.bX     = TRUE;
    fr.x      = x;
    fr.bBox   = TRUE;
    copy_mat(box, fr.box);
    natoms = fr.natoms;
  }
  
  /* check input */
  if ( natoms == 0 )
    gmx_fatal(FARGS,"no atoms in coordinate file %s",infile);
  if ( natoms > atoms->nr )
    gmx_fatal(FARGS,"topology with %d atoms does not match "
		"coordinates with %d atoms",atoms->nr,natoms);
  for(i=0; i<nidx; i++)
    if (index[i] > natoms)
      gmx_fatal(FARGS,"An atom number in group %s is larger than the number of "
		  "atoms (%d) in the coordinate file %s",grpnm,natoms,infile);
  
  /* get indexed copy of atoms */
  snew(iatoms,1);
  init_t_atoms(iatoms,nidx,FALSE);
  snew(iatoms->atom, iatoms->nr);
  resind = 0;
  for(i=0; i<nidx; i++) {
    iatoms->atom[i] = atoms->atom[index[i]];
    iatoms->atomname[i] = atoms->atomname[index[i]];
    if ( i>0 && (atoms->atom[index[i]].resind!=atoms->atom[index[i-1]].resind) )
      resind++;
    iatoms->atom[i].resind = resind;
    iatoms->resinfo[resind] = atoms->resinfo[atoms->atom[index[i]].resind];
    iatoms->nres = max(iatoms->nres, iatoms->atom[i].resind+1);
  }
  
  init_t_protonate(&protdata);
  
  out = open_trx(opt2fn("-o",NFILE,fnm),"w");
  snew(ix, nidx);
  frame=0;
  do {
    if (debug) fprintf(debug,"FRAME %d (%d %g)\n",frame,fr.step,fr.time);
    /* get indexed copy of x */
    for(i=0; i<nidx; i++)
      copy_rvec(fr.x[index[i]], ix[i]);
    /* protonate */
    natoms_out = protonate(&iatoms, &ix, &protdata);
    
    /* setup output frame */
    frout = fr;
    frout.natoms = natoms_out;
    frout.bAtoms = TRUE;
    frout.atoms  = iatoms;
    frout.bV     = FALSE;
    frout.bF     = FALSE;
    frout.x      = ix;

    /* write output */
    write_trxframe(out,&frout,NULL);
    frame++;
  } while ( bReadMultiple && read_next_frame(oenv,status, &fr) );
  
  thanx(stderr);

  return 0;
}

