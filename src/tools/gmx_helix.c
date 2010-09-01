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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include "confio.h"
#include "copyrite.h"
#include "gmx_fatal.h"
#include "fitahx.h"
#include "futil.h"
#include "gstat.h"
#include "wgms.h"
#include "hxprops.h"
#include "macros.h"
#include "maths.h"
#include "pbc.h"
#include "tpxio.h"
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "xvgr.h"
#include "gmx_ana.h"


void dump_ahx(int nres,
	      t_bb bb[],rvec x[],matrix box,int teller)
{
  FILE *fp;
  char buf[256];
  int  i;
  
  sprintf(buf,"dump%d.gro",teller);
  fp=ffopen(buf,"w");
  fprintf(fp,"Dumping fitted helix frame %d\n",teller);
  fprintf(fp,"%5d\n",nres*5);
  for(i=0; (i<nres); i++) {
#define PR(AA) fprintf(fp,"%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",i+1,"GLY",#AA,bb[i].AA,x[bb[i].AA][XX],x[bb[i].AA][YY],x[bb[i].AA][ZZ]); fflush(fp)
    if (bb[i].bHelix) {
      PR(N);
      PR(H);
      PR(CA);
      PR(C);
      PR(O);
    }
  }
  for(i=0; (i<DIM); i++)
    fprintf(fp,"%10.5f",box[i][i]);
  fprintf(fp,"\n");
  ffclose(fp);
}

void dump_otrj(FILE *otrj,int natoms,atom_id all_index[],rvec x[],
	       real fac,rvec xav[])
{
  static FILE *fp=NULL;
  static real fac0=1.0;
  
  int  i,ai,m;
  real xa,xm;
  
  if (fp==NULL) {
    fp=ffopen("WEDGAMMA10.DAT","w");
    fac0=fac;
  }
  fac/=fac0;
  
  fprintf(fp,"%10g\n",fac);
  
  for(i=0; (i<natoms); i++) {
    ai=all_index[i];
    for(m=0; (m<DIM); m++) {
      xa = xav[ai][m];
      xm = x[ai][m]*fac;
      xav[ai][m] = xa+xm;
      x[ai][m]   = xm;
    }
  }
  write_gms_ndx(otrj,natoms,all_index,x,NULL);
}

int gmx_helix(int argc,char *argv[])
{
  const char *desc[] = {
    "g_helix computes all kind of helix properties. First, the peptide",
    "is checked to find the longest helical part. This is determined by",
    "Hydrogen bonds and Phi/Psi angles.",
    "That bit is fitted",
    "to an ideal helix around the Z-axis and centered around the origin.",
    "Then the following properties are computed:[PAR]",
    "[BB]1.[bb] Helix radius (file radius.xvg). This is merely the",
    "RMS deviation in two dimensions for all Calpha atoms.",
    "it is calced as sqrt((SUM i(x^2(i)+y^2(i)))/N), where N is the number",
    "of backbone atoms. For an ideal helix the radius is 0.23 nm[BR]",
    "[BB]2.[bb] Twist (file twist.xvg). The average helical angle per",
    "residue is calculated. For alpha helix it is 100 degrees,",
    "for 3-10 helices it will be smaller,", 
    "for 5-helices it will be larger.[BR]",
    "[BB]3.[bb] Rise per residue (file rise.xvg). The helical rise per", 
    "residue is plotted as the difference in Z-coordinate between Ca", 
    "atoms. For an ideal helix this is 0.15 nm[BR]",
    "[BB]4.[bb] Total helix length (file len-ahx.xvg). The total length", 
    "of the", 
    "helix in nm. This is simply the average rise (see above) times the",  
    "number of helical residues (see below).[BR]",
    "[BB]5.[bb] Number of helical residues (file n-ahx.xvg). The title says",
    "it all.[BR]",
    "[BB]6.[bb] Helix Dipole, backbone only (file dip-ahx.xvg).[BR]",
    "[BB]7.[bb] RMS deviation from ideal helix, calculated for the Calpha",
    "atoms only (file rms-ahx.xvg).[BR]",
    "[BB]8.[bb] Average Calpha-Calpha dihedral angle (file phi-ahx.xvg).[BR]",
    "[BB]9.[bb] Average Phi and Psi angles (file phipsi.xvg).[BR]",
    "[BB]10.[bb] Ellipticity at 222 nm according to [IT]Hirst and Brooks[it]",
    "[PAR]"
  };
  static const char *ppp[efhNR+2] = { 
    NULL, "RAD", "TWIST", "RISE", "LEN", "NHX", "DIP", "RMS", "CPHI", 
    "RMSA", "PHI", "PSI", "HB3", "HB4", "HB5", "CD222", NULL
  };
  static gmx_bool bCheck=FALSE,bFit=TRUE,bDBG=FALSE,bEV=FALSE;
  static int  rStart=0,rEnd=0,r0=1;
  t_pargs pa [] = {
    { "-r0", FALSE, etINT, {&r0},
      "The first residue number in the sequence" },
    { "-q",  FALSE, etBOOL,{&bCheck},
      "Check at every step which part of the sequence is helical" },
    { "-F",  FALSE, etBOOL,{&bFit},
      "Toggle fit to a perfect helix" },
    { "-db", FALSE, etBOOL,{&bDBG},
      "Print debug info" },
    { "-prop", FALSE, etENUM, {ppp},
      "Select property to weight eigenvectors with. WARNING experimental stuff" },
    { "-ev", FALSE, etBOOL,{&bEV},
      "Write a new 'trajectory' file for ED" },
    { "-ahxstart", FALSE, etINT, {&rStart},
      "First residue in helix" },
    { "-ahxend", FALSE, etINT, {&rEnd},
      "Last residue in helix" }
  };

  typedef struct {
    FILE *fp,*fp2;
    gmx_bool bfp2;
    const char *filenm;
    const char *title;
    const char *xaxis;
    const char *yaxis;
    real val;
  } t_xvgrfile;
  
  t_xvgrfile xf[efhNR] = {
    { NULL, NULL, TRUE,  "radius",  "Helix radius",               NULL, "r (nm)" , 0.0 },
    { NULL, NULL, TRUE,  "twist",   "Twist per residue",          NULL, "Angle (deg)", 0.0 },
    { NULL, NULL, TRUE,  "rise",    "Rise per residue",           NULL, "Rise (nm)", 0.0 },
    { NULL, NULL, FALSE, "len-ahx", "Length of the Helix",        NULL, "Length (nm)", 0.0 },
    { NULL, NULL, FALSE, "dip-ahx", "Helix Backbone Dipole",      NULL, "rq (nm e)", 0.0 },
    { NULL, NULL, TRUE,  "rms-ahx", "RMS Deviation from Ideal Helix", NULL, "RMS (nm)", 0.0 },
    { NULL, NULL, FALSE, "rmsa-ahx","Average RMSD per Residue",   "Residue", "RMS (nm)", 0.0 },
    { NULL, NULL,FALSE,  "cd222",   "Ellipticity at 222 nm", NULL, "nm", 0.0 },
    { NULL, NULL, TRUE,  "pprms",   "RMS Distance from \\8a\\4-helix", NULL, "deg" , 0.0 },
    { NULL, NULL, TRUE,  "caphi",   "Average Ca-Ca Dihedral",     NULL, "\\8F\\4(deg)", 0.0 },
    { NULL, NULL, TRUE,  "phi",     "Average \\8F\\4 angles", NULL, "deg" , 0.0 },
    { NULL, NULL, TRUE,  "psi",     "Average \\8Y\\4 angles", NULL, "deg" , 0.0 },
    { NULL, NULL, TRUE,  "hb3",     "Average n-n+3 hbond length", NULL, "nm" , 0.0 },
    { NULL, NULL, TRUE,  "hb4",     "Average n-n+4 hbond length", NULL, "nm" , 0.0 },
    { NULL, NULL, TRUE,  "hb5",     "Average n-n+5 hbond length", NULL, "nm" , 0.0 },
    { NULL, NULL,FALSE,  "JCaHa",   "J-Coupling Values",        "Residue", "Hz" , 0.0 },
    { NULL, NULL,FALSE,  "helicity","Helicity per Residue",     "Residue", "% of time" , 0.0 }
  };
 
  output_env_t oenv;
  FILE       *otrj;
  char       buf[54],prop[256];
  t_trxstatus *status;
  int        natoms,nre,nres;
  t_bb       *bb;
  int        i,j,ai,m,nall,nbb,nca,teller,nSel=0;
  atom_id    *bbindex,*caindex,*allindex;
  t_topology *top;
  int        ePBC;
  rvec       *x,*xref,*xav;
  real       t;
  real       rms,fac;
  matrix     box;
  gmx_rmpbc_t  gpbc=NULL;
  gmx_bool       bRange;
  t_filenm  fnm[] = {
    { efTPX, NULL,  NULL,   ffREAD  },
    { efNDX, NULL,  NULL,   ffREAD  },
    { efTRX, "-f",  NULL,   ffREAD  },
    { efG87, "-to", NULL,   ffOPTWR },
    { efSTO, "-cz", "zconf",ffWRITE },
    { efSTO, "-co", "waver",ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);
  
  bRange=(opt2parg_bSet("-ahxstart",asize(pa),pa) &&
	  opt2parg_bSet("-ahxend",asize(pa),pa));
		        
  top=read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);
  
  natoms=read_first_x(oenv,&status,opt2fn("-f",NFILE,fnm),&t,&x,box);

  if (opt2bSet("-to",NFILE,fnm)) {
    otrj=opt2FILE("-to",NFILE,fnm,"w");
    strcpy(prop,ppp[0]);
    fprintf(otrj,"%s Weighted Trajectory: %d atoms, NO box\n",prop,natoms);
  }
  else
    otrj=NULL;
    
  if (natoms != top->atoms.nr)
    gmx_fatal(FARGS,"Sorry can only run when the number of atoms in the run input file (%d) is equal to the number in the trajectory (%d)",
	    top->atoms.nr,natoms);
	    
  bb=mkbbind(ftp2fn(efNDX,NFILE,fnm),&nres,&nbb,r0,&nall,&allindex,
	     top->atoms.atomname,top->atoms.atom,top->atoms.resinfo);
  snew(bbindex,natoms);
  snew(caindex,nres);
  
  fprintf(stderr,"nall=%d\n",nall);
    
  /* Open output files, default x-axis is time */
  for(i=0; (i<efhNR); i++) {
    sprintf(buf,"%s.xvg",xf[i].filenm);
    remove(buf);
    xf[i].fp=xvgropen(buf,xf[i].title,
                      xf[i].xaxis ? xf[i].xaxis : "Time (ps)",
		      xf[i].yaxis,oenv);
    if (xf[i].bfp2) {
      sprintf(buf,"%s.out",xf[i].filenm);
      remove(buf);
      xf[i].fp2=ffopen(buf,"w");
    }
  }

  /* Read reference frame from tpx file to compute helix length */
  snew(xref,top->atoms.nr);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),
	   NULL,NULL,&natoms,xref,NULL,NULL,NULL);
  calc_hxprops(nres,bb,xref,box);
  do_start_end(nres,bb,xref,&nbb,bbindex,&nca,caindex,bRange,rStart,rEnd);
  sfree(xref);
  if (bDBG) {
    fprintf(stderr,"nca=%d, nbb=%d\n",nca,nbb);
    pr_bb(stdout,nres,bb);
  }
  
  gpbc = gmx_rmpbc_init(&top->idef,ePBC,natoms,box);

  snew(xav,natoms);
  teller=0;
  do {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rt=%.2f",t);
    gmx_rmpbc(gpbc,natoms,box,x);

    
    calc_hxprops(nres,bb,x,box);
    if (bCheck)
      do_start_end(nres,bb,x,&nbb,bbindex,&nca,caindex,FALSE,0,0);
    
    if (nca >= 5) {
      rms=fit_ahx(nres,bb,natoms,nall,allindex,x,nca,caindex,box,bFit);
      
      if (teller == 1) {
	write_sto_conf(opt2fn("-cz",NFILE,fnm),"Helix fitted to Z-Axis",
		       &(top->atoms),x,NULL,ePBC,box);
      }
            
      xf[efhRAD].val   = radius(xf[efhRAD].fp2,nca,caindex,x);
      xf[efhTWIST].val = twist(xf[efhTWIST].fp2,nca,caindex,x);
      xf[efhRISE].val  = rise(nca,caindex,x);
      xf[efhLEN].val   = ahx_len(nca,caindex,x,box);
      xf[efhCD222].val = ellipticity(nres,bb);
      xf[efhDIP].val   = dip(nbb,bbindex,x,top->atoms.atom);
      xf[efhRMS].val   = rms;
      xf[efhCPHI].val  = ca_phi(nca,caindex,x,box);
      xf[efhPPRMS].val = pprms(xf[efhPPRMS].fp2,nres,bb);
      
      for(j=0; (j<=efhCPHI); j++)
	fprintf(xf[j].fp,   "%10g  %10g\n",t,xf[j].val);
      
      av_phipsi(xf[efhPHI].fp,xf[efhPSI].fp,xf[efhPHI].fp2,xf[efhPSI].fp2,
		t,nres,bb);
      av_hblen(xf[efhHB3].fp,xf[efhHB3].fp2,
	       xf[efhHB4].fp,xf[efhHB4].fp2,
	       xf[efhHB5].fp,xf[efhHB5].fp2,
	       t,nres,bb);
      
      if (otrj) 
	dump_otrj(otrj,nall,allindex,x,xf[nSel].val,xav);
    }
  } while (read_next_x(oenv,status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  gmx_rmpbc_done(gpbc);

  close_trj(status);

  if (otrj) {
    ffclose(otrj);
    fac=1.0/teller;
    for(i=0; (i<nall); i++) {
      ai=allindex[i];
      for(m=0; (m<DIM); m++)
	xav[ai][m]*=fac;
    }
    write_sto_conf_indexed(opt2fn("-co",NFILE,fnm),
			   "Weighted and Averaged conformation",
			   &(top->atoms),xav,NULL,ePBC,box,nall,allindex);
  }
  
  for(i=0; (i<nres); i++) {
    if (bb[i].nrms > 0) {
      fprintf(xf[efhRMSA].fp,"%10d  %10g\n",r0+i,bb[i].rmsa/bb[i].nrms);
    }
    fprintf(xf[efhAHX].fp,"%10d  %10g\n",r0+i,(bb[i].nhx*100.0)/(real )teller);
    fprintf(xf[efhJCA].fp,"%10d  %10g\n",
	    r0+i,140.3+(bb[i].jcaha/(double)teller));
  }
  
  for(i=0; (i<efhNR); i++) {
    ffclose(xf[i].fp);
    if (xf[i].bfp2)
      ffclose(xf[i].fp2);
    do_view(oenv,xf[i].filenm,"-nxy");
  }
  
  thanx(stderr);
  
  return 0;
}

