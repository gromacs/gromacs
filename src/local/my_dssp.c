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
 * GROtesk MACabre and Sinister
 */
static char *SRCID_my_dssp_c = "$Id$";

#include "sysstuff.h"
#include "typedefs.h"
#include "string2.h"
#include "strdb.h"
#include "macros.h"
#include "smalloc.h"
#include "mshift.h"
#include "statutil.h"
#include "copyrite.h"
#include "confio.h"
#include "fatal.h"
#include "xvgr.h"
#include "matio.h"
#include "rdgroup.h"
#include "gstat.h"

extern void dssp_main(bool bDoAcc, bool bVerbose);

extern FILE *tapein, *tapeout;

static void strip_dssp(int nres,int r0,
		       bool bPhobres[],real t,real dt,
		       FILE *ss,FILE *acc,FILE *acct,
		       t_matrix *mat,int average_area[])
{
  static bool bFirst=TRUE;
  static char *ssbuf;
  static int xsize,frame;
  char buf[STRLEN+1];
  char SSTP;
  int  i,nr,iacc;
  real iaccf,iaccb;
  t_xpmelmt c;
  
  /* Skip header */
  do {
    fgets2(buf,STRLEN,tapeout);
  } while (strstr(buf,"KAPPA") == NULL);
  if (bFirst) {
    snew(ssbuf,nres+10);
  }
  
  iaccb=iaccf=0;
  for(nr=0; (fgets2(buf,STRLEN,tapeout) != NULL); nr++) {
    SSTP=buf[16]==' ' ? '~' : buf[16];
    ssbuf[nr]=SSTP;
    
    buf[39]='\0';
    sscanf(&(buf[34]),"%d",&iacc);
    if (acc)
      fprintf(acc,"%4d ",iacc);
    average_area[nr]+=iacc;
    if (bPhobres[nr])
      iaccb+=iacc;
    else
      iaccf+=iacc;
  }
  ssbuf[nr]='\0';
  
  if (ss)
    fprintf(ss,"%10g  & {\\tt %s } \\\\\n",t,ssbuf);
  
  if (bFirst) {
    sprintf(mat->title,"Secondary structure");
    mat->legend[0]=0;
    sprintf(mat->label_x,"Time (ps)");
    sprintf(mat->label_y,"Residue");
    mat->bDiscrete=TRUE;
    mat->ny=nr;
    snew(mat->axis_y,nr);
    for(i=0; i<nr; i++)
      mat->axis_y[i]=r0+i;
    mat->axis_x=NULL;
    mat->matrix=NULL;
    xsize=0;
    frame=0;
    bFirst=FALSE;
  }
  if (frame>=xsize) {
    xsize+=10;
    srenew(mat->axis_x,xsize);
    srenew(mat->matrix,xsize);
  }
  mat->axis_x[frame]=t;
  snew(mat->matrix[frame],nr);
  c.c2=0;
  for(i=0; i<nr; i++) {
    c.c1=ssbuf[i];
    mat->matrix[frame][i]=searchcmap(mat->nmap,mat->map,c);
  }
  frame++;
  mat->nx=frame;
  
  if (acc)
    fprintf(acc,"\n");
  if (acct)
    fprintf(acct,"%10g  %10g  %10g\n",t,0.01*iaccb,0.01*iaccf);
}

bool *bPhobics(t_atoms *atoms)
{
  int  i,nb;
  char **cb;
  bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,atoms->nres+10);
  
  for(i=0; (i<atoms->nres); i++) {
    if (search_str(nb,cb,*atoms->resname[i]) != -1)
      bb[i]=TRUE;
  }
  return bb;
}
 
static void check_oo(t_atoms *atoms)
{
  static char *OOO="O";
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    if (strcmp(*(atoms->atomname[i]),"OXT")==0)
      atoms->atomname[i]=&OOO;
    else if (strcmp(*(atoms->atomname[i]),"O1")==0)
      atoms->atomname[i]=&OOO;
  }
}

static void norm_acc(char *surffn, t_atoms *atoms, 
		     real av_area[], real norm_av_area[])
{
  int i,n,n_surf;
  
  char **surf_res, **surf_lines;
  double *surf;
  
  n_surf = get_lines(surffn, &surf_lines);
  snew(surf, n_surf);
  snew(surf_res, n_surf);
  for(i=0; (i<n_surf); i++) {
    snew(surf_res[i], 5);
    sscanf(surf_lines[i],"%s %lf",surf_res[i],&surf[i]);
  }
  
  for(i=0; (i<atoms->nres); i++) {
    n = search_str(n_surf,surf_res,*atoms->resname[i]);
    if ( n != -1)
      norm_av_area[i] = av_area[i] / surf[n];
    else 
      fprintf(stderr,"Residue %s not found in surface database (%s)\n",
	      *atoms->resname[i],surffn);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "my_dssp reads a trajectory file and computes the secondary",
    "structure for each time frame (or every [BB]dt[bb] ps) by",
    "using the dssp program. Also solvent accessible surface per",
    "residue is calculated, both in absolute values (A^2) and in",
    "fractions of the maximal accessible surface of a residue.",
    "The maximal accessible surface is defined as the accessible",
    "surface of a residue in a chain of glycines.[PAR]",
    "The output of dssp may be visualized, as a table in [BB]LaTeX[bb],",
    "as an X PixMap or as color postscript using a postprocessing program",
    "called xpm2ps."
  };
  static real dt=0.0;
  static int  r0=1;
  static bool bVerbose;
  t_pargs pa[] = {
    { "-v",  FALSE, etBOOL, &bVerbose,
      "Generate miles of useless information" },
    { "-dt", FALSE, etREAL, &dt,
      "Time interval between frames." },
    { "-r0", FALSE, etINT,  &r0,
      "Starting residue for output files" }
  };
  int        status;
  FILE       *ss,*acc,*acct;
  char       *leg[] = { "Phobic", "Phylic" };
  t_topology *top;
  t_atoms    *atoms;
  t_matrix   mat;
  int        nres,nr0;
  bool       *bPhbres,bDoAcc;
  real       t,nt;
  int        i,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm;
  atom_id    *index;
  rvec       *x;
  int        *average_area;
  real       *av_area, *norm_av_area;
  char       pdbfile[L_tmpnam],tmpfile[L_tmpnam];
  char       dssp[256],*dptr;
#define MAXBUF 1000000
  char       inbuf[MAXBUF],outbuf[MAXBUF];
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,      ffREAD },
    { efTPX, NULL,   NULL,      ffREAD },
    { efNDX, NULL,   NULL,      ffOPTRD },
    { efMAP, "-map", "ss",      ffLIBRD },
    { efDAT, "-sf",  "surface", ffLIBRD },
    { efXPM, "-ss",  "ss",      ffWRITE },
    { efTEX, "-os",  "ss",      ffOPTWR },
    { efOUT, "-oa",  "area",    ffOPTWR },
    { efXVG, "-a",   "totarea", ffOPTWR },
    { efXVG, "-aa",  "averarea",ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
	 	    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  bDoAcc=(ftp2bSet(efOUT,NFILE,fnm) || 
	  opt2bSet("-a" ,NFILE,fnm) || 
	  opt2bSet("-aa",NFILE,fnm));
  if (bDoAcc) 
    printf("Will do accessible surface calc\n");
  else
    printf("Will *NOT* do accessible surface calc\n");
  top=read_top(ftp2fn(efTPX,NFILE,fnm));
  atoms=&(top->atoms);
  check_oo(atoms);
  nres=atoms->nres;
  bPhbres=bPhobics(atoms);
  
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);
  nres=0;
  nr0=-1;
  for(i=0; (i<gnx); i++) {
    if (atoms->atom[index[i]].resnr != nr0) {
      nr0=atoms->atom[index[i]].resnr;
      nres++;
    }
  }
  fprintf(stderr,"There are %d residues in your selected group\n",nres);
  
  /* Open all files read-write */
  (void) tmpnam(pdbfile);
  tapein=ffopen(pdbfile,"w+");
  setvbuf(tapein,inbuf,MAXBUF,_IOFBF);
  
  (void) tmpnam(tmpfile);
  tapeout=ffopen(tmpfile,"w+");
  setvbuf(tapeout,outbuf,MAXBUF,_IOFBF);
  
  if ((dptr=getenv("DSSP")) == NULL)
    dptr="/home/mdgroup/dssp/dssp";
  sprintf(dssp,"%s %s %s > /dev/null",dptr,pdbfile,tmpfile);
  fprintf(stderr,"dssp cmd='%s'\n",dssp);
  
  if (ftp2bSet(efOUT,NFILE,fnm))
    acc=ftp2FILE(efOUT,NFILE,fnm,"w");
  else
    acc=NULL;
    
  if (opt2bSet("-a",NFILE,fnm)) {
    acct=xvgropen(opt2fn("-a",NFILE,fnm),"Solvent Accessible Surface Area",
		  "Time (ps)","Area (nm\\S2\\N)");
    xvgr_legend(acct,2,leg);
  } else
    acct=NULL;
  
  if (ftp2bSet(efTEX,NFILE,fnm))
    ss=ffopen(ftp2fn(efTEX,NFILE,fnm),"w");
  else
    ss=NULL;
  
  mat.map=NULL;
  mat.nmap=getcmap(libopen(opt2fn("-map",NFILE,fnm)),
		   opt2fn("-map",NFILE,fnm),&(mat.map));
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms != atoms->nr) 
    fatal_error(0,"Trajectory does not match topology!");
  
  snew(average_area,atoms->nres+10);
  snew(av_area,atoms->nres+10);
  snew(norm_av_area,atoms->nres+10);
  nt=t;
  do {
    if (t >= nt) {
      rm_pbc(&(top->idef),atoms->nr,box,x,x);
      hwrite_pdb_conf_indexed(tapein,atoms,x,box,gnx,index);
      rewind(tapein);
      dssp_main(bDoAcc,bVerbose);
      rewind(tapein);
      rewind(tapeout);
      strip_dssp(nres,r0,bPhbres,t,dt,ss,acc,acct,&mat,average_area);
      rewind(tapeout);
      /*remove(pdbfile);*/
      nt+=dt;
      nframe++;
    }
  } while(read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  if (ss)
    ffclose(ss);
  if (acc)
    ffclose(acc);
  if (acct)
    ffclose(acct);
  
  ss=ffopen(ftp2fn(efXPM,NFILE,fnm),"w");
  write_xpm_m(ss,mat);
  ffclose(ss);
  remove(tmpfile);
  remove(pdbfile);

  for(i=0; (i<atoms->nres); i++)
    av_area[i] = average_area[i]/(real) nframe;
    
  norm_acc(opt2fn("-sf",NFILE,fnm), atoms, av_area, norm_av_area);
  
  if (opt2bSet("-aa",NFILE,fnm)) {
    acc=xvgropen(opt2fn("-aa",NFILE,fnm),"Average Accessible Area",
		 "Residue","A\\S2");
    for(i=0; (i<atoms->nres); i++) {
      fprintf(acc,"%5d  %10g %10g\n",i+r0,av_area[i], norm_av_area[i]);
    }
    ffclose(acc);
  }
  
  if (opt2bSet("-a",NFILE,fnm))
    xvgr_file(opt2fn("-a",NFILE,fnm),"-nxy");
  if (opt2bSet("-aa",NFILE,fnm))
    xvgr_file(opt2fn("-aa",NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}
