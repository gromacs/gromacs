/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_g_helix_c = "$Id$";

#include <math.h>
#include "bondf.h"
#include "confio.h"
#include "copyrite.h"
#include "fatal.h"
#include "fitahx.h"
#include "futil.h"
#include "gstat.h"
#include "hxprops.h"
#include "macros.h"
#include "maths.h"
#include "pbc.h"
#include "tpxio.h"
#include "physics.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "statutil.h"
#include "string.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "typedefs.h"
#include "vec.h"
#include "xvgr.h"

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
  fclose(fp);
}

int get_prop(char prop[])
{
  static char *ppp[efhNR] = { 
    "RAD", "TWIST", "RISE", "LEN", "NHX", "DIP", "RMS", "CPHI", 
    "RMSA", "PHI", "PSI", "HB3", "HB4", "HB5", "CD222"
  };
  int i,n;
  
  fprintf(stderr,"Select something\n");
  for(i=0; (i<efhNR); ) {
    fprintf(stderr,"%2d: %10s",i,ppp[i]); i++;
    if (i<efhNR) {
      fprintf(stderr,"      %2d: %10s\n",i,ppp[i]);
      i++;
    }
    else
      fprintf(stderr,"\n");
  }
  do {
    scanf("%d",&n);
  } while ((n < 0) || (n >= efhNR));
  
  strcpy(prop,ppp[n]);
  
  return n;
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
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
  static bool bCheck=FALSE,bFit=TRUE,bDBG=FALSE,bEV=FALSE;
  static int  rStart=0,rEnd=0,r0=1;
  t_pargs pa [] = {
    { "-r0", FALSE, etINT, &r0,
      "The first residue number in the sequence" },
    { "-q",  FALSE, etBOOL,&bCheck,
      "Check at every step which part of the sequence is helical" },
    { "-F",  FALSE, etBOOL,&bFit,
      "Toggle fit to a perfect helix" },
    { "-db", FALSE, etBOOL,&bDBG,
      "Print debug info" },
    { "-ev", FALSE, etBOOL,&bEV,
      "Write a new 'trajectory' file for ED" },
    { "-ahxstart", FALSE, etINT, &rStart,
      "First residue in helix" },
    { "-ahxend", FALSE, etINT, &rEnd,
      "Last residue in helix" }
  };

  typedef struct {
    FILE *fp,*fp2;
    bool bfp2;
    char *filenm;
    char *title;
    char *xaxis;
    char *yaxis;
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
  
  FILE       *otrj;
  char       buf[54],prop[256];
  int        status;
  int        natoms,nre,nres,step;
  t_bb       *bb;
  int        i,j,ai,m,nall,nbb,nca,teller,nSel=0;
  atom_id    *bbindex,*caindex,*allindex;
  t_topology *top;
  rvec       *x,*xref,*xav;
  real       t,lambda;
  real       rms,fac;
  matrix     box;
  bool       bRange;
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
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  bRange=(opt2parg_bSet("-ahxstart",asize(pa),pa) &&
	  opt2parg_bSet("-ahxend",asize(pa),pa));
		        
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  
  natoms=read_first_x(&status,opt2fn("-f",NFILE,fnm),&t,&x,box);

  if (opt2bSet("-to",NFILE,fnm)) {
    otrj=opt2FILE("-to",NFILE,fnm,"w");
    nSel=get_prop(prop);
    fprintf(otrj,"%s Weighted Trajectory: %d atoms, NO box\n",prop,natoms);
  }
  else
    otrj=NULL;
    
  if (natoms != top->atoms.nr)
    fprintf(stderr,"Warning! in tpx are %d atoms, in trj %d\n",
	    top->atoms.nr,natoms);
	    
  bb=mkbbind(ftp2fn(efNDX,NFILE,fnm),&nres,&nbb,r0,&nall,&allindex,
	     top->atoms.atomname,top->atoms.atom,top->atoms.resname);
  snew(bbindex,natoms);
  snew(caindex,nres);
  
  fprintf(stderr,"nall=%d\n",nall);
    
  /* Open output files, default x-axis is time */
  for(i=0; (i<efhNR); i++) {
    sprintf(buf,"%s.xvg",xf[i].filenm);
    remove(buf);
    xf[i].fp=xvgropen(buf,xf[i].title,xf[i].xaxis ? xf[i].xaxis : "Time (ps)",
		      xf[i].yaxis);
    if (xf[i].bfp2) {
      sprintf(buf,"%s.out",xf[i].filenm);
      remove(buf);
      xf[i].fp2=ffopen(buf,"w");
    }
  }

  /* Read reference frame from tpx file to compute helix length */
  snew(xref,natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),
	   &step,&t,&lambda,NULL,NULL,
	   &natoms,xref,NULL,NULL,NULL);
  calc_hxprops(nres,bb,xref,box);
  do_start_end(nres,bb,xref,&nbb,bbindex,&nca,caindex,bRange,rStart,rEnd);
  sfree(xref);
  if (bDBG) {
    fprintf(stderr,"nca=%d, nbb=%d\n",nca,nbb);
    pr_bb(stdout,nres,bb);
  }
  
  snew(xav,natoms);
  teller=0;
  do {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rt=%.2f",t);
    rm_pbc(&(top->idef),top->atoms.nr,box,x,x);
    
    calc_hxprops(nres,bb,x,box);
    if (bCheck)
      do_start_end(nres,bb,x,&nbb,bbindex,&nca,caindex,FALSE,0,0);
    
    if (nca >= 5) {
/*#define DEBUG*/
#ifdef DEBUG
      fprintf(stderr,"nca=%d, nbb=%d\n",nca,nbb);
      dump_ahx(nres,bb,x,box,3);
#endif
      rms=fit_ahx(nres,bb,natoms,nall,allindex,x,nca,caindex,box,bFit);
      
      if (teller == 1) {
	write_sto_conf(opt2fn("-cz",NFILE,fnm),"Helix fitted to Z-Axis",
		       &(top->atoms),x,NULL,box);
      }
            
#ifdef DEBUG
      dump_ahx(nres,bb,x,box,4);
      exit(1);
#endif

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
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);

  if (otrj) {
    fclose(otrj);
    fac=1.0/teller;
    for(i=0; (i<nall); i++) {
      ai=allindex[i];
      for(m=0; (m<DIM); m++)
	xav[ai][m]*=fac;
    }
    write_sto_conf_indexed(opt2fn("-co",NFILE,fnm),
			   "Weighted and Averaged conformation",
			   &(top->atoms),xav,NULL,box,nall,allindex);
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
    fclose(xf[i].fp);
    if (xf[i].bfp2)
      fclose(xf[i].fp2);
    xvgr_file(xf[i].filenm,NULL);
  }
  
  thanx(stdout);
  
  return 0;
}

