/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Gnomes, ROck Monsters And Chili Sauce
 */
static char *SRCID_do_dssp_c = "$Id$";

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
#include "readcmap.h"

static void strip_dssp(char *dsspfile,int nres,int r0,
		       bool bPhobres[],real t,real dt,
		       FILE *ss,FILE *acc,FILE *acct,t_matrix *mat)
{
  static bool bFirst=TRUE;
  static char *ssbuf;
  FILE *in;
  static int xsize,frame;
  char buf[STRLEN+1];
  char SSTP;
  int  i,j,nr,iacc;
  real iaccf,iaccb;
  
  in=ffopen(dsspfile,"r");
  
  /* Skip header */
  do {
    fgets2(buf,STRLEN,in);
  } while (strstr(buf,"KAPPA") == NULL);
  if (bFirst) {
    snew(ssbuf,nres+10);
  }
  
  iaccb=iaccf=0;
  for(nr=0; (fgets2(buf,STRLEN,in) != NULL); nr++) {
    SSTP=buf[16]==' ' ? '~' : buf[16];
    ssbuf[nr]=SSTP;
    
    buf[39]='\0';
    sscanf(&(buf[34]),"%d",&iacc);
    fprintf(acc,"%4d ",iacc);
    if (bPhobres[nr])
      iaccb+=iacc;
    else
      iaccf+=iacc;
  }
  ssbuf[nr]='\0';
  
  
  fprintf(ss,"%10g  & {\\tt %s } \\\\\n",t,ssbuf);
  
  if (bFirst) {
    /* We can only write the header here, because DSSP may not recognize
     * all the residues that we have (such as ACE, NAC, etc)
     * And therefore the number of residues from the topology is not
     * the same as the one we have here.
     */
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
  for(i=0; i<nr; i++)
    mat->matrix[frame][i]=ssbuf[i];
  frame++;
  mat->nx=frame;
  
  fprintf(acc,"\n");
  fprintf(acct,"%10g  %10g  %10g\n",t,0.01*iaccb,0.01*iaccf);
  fclose(in);
}

bool *bPhobics(t_atoms *atoms)
{
  int  i,nb;
  char **cb;
  bool *bb;
  
  nb=get_strings("phbres.dat",&cb);
  snew(bb,atoms->nres);
  
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

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "do_dssp reads a trajectory file and computes the secondary",
    "structure for each time frame (or every [BB]dt[bb] ps) by",
    "calling the dssp program. If you do not have the dssp program,",
    "get it. do_dssp assumes that the dssp executable is in",
    "/home/mdgroup/dssp/dssp. If that is not the case, then you should",
    "set an environment variable [BB]DSSP[bb] as in: [PAR]",
    "[TT]setenv DSSP /usr/local/bin/dssp[tt][PAR]",
    "where the right hand side should point to the dssp executable.[PAR]",
    "The output of dssp may be visualized, as a table in [BB]LaTeX[bb],",
    "as an X PixMap or as color postscript using a postprocessing program",
    "called xpm2ps."
  };
  static int  r0=1;
  static real dt=0.0;
  t_pargs pa[] = {
    { "-dt", FALSE, etREAL, &dt,
      "Time interval between frames." },
    { "-r0", FALSE, etINT,  &r0,
      "For proteins or peptides which do not start at residue 1 , this option sets the starting residue to be the correct one in the output files" }
  };
  static char *bugs[] = {
    "The program is very slow"
  };
  int        status;
  FILE       *ss,*acc,*acct;
  char       *leg[] = { "Phobic", "Phylic" };
  t_topology *top;
  t_atoms    *atoms;
  t_matrix   mat;
  int        nres,nr0;
  bool       *bPhbres;
  real       t,nt;
  int        i,natoms,nframe=0;
  matrix     box;
  int        gnx;
  char       *grpnm;
  atom_id    *index;
  rvec       *x;
  char       pdbfile[L_tmpnam],tmpfile[L_tmpnam];
  char       dssp[256],*dptr;
  t_filenm   fnm[] = {
    { efTRX, "-f",   NULL,      ffREAD },
    { efTPB, NULL,   NULL,      ffREAD },
    { efNDX, NULL,   NULL,      ffOPTRD },
    { efMAP, "-map", "ss",      ffLIBRD },
    { efTEX, "-os",  "ss",      ffWRITE },
    { efOUT, "-oa",  "area",    ffWRITE },
    { efXVG, "-a",   "totarea", ffWRITE },
    { efXPM, "-ss",  "ss",      ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
	 	    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
		    
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
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
  
  (void) tmpnam(pdbfile);
  (void) tmpnam(tmpfile);
  if ((dptr=getenv("DSSP")) == NULL)
    dptr="/home/mdgroup/dssp/dssp";
  sprintf(dssp,"%s %s %s > /dev/null",dptr,pdbfile,tmpfile);
  fprintf(stderr,"dssp cmd='%s'\n",dssp);
  
  ss=opt2FILE("-os",NFILE,fnm,"w");
  acc=opt2FILE("-oa",NFILE,fnm,"w");
  acct=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Solvent Accessible Surface Area",
		"Time (ps)","Area (nm\\S2\\N)");
  mat.map=NULL;
  mat.nmap=getcmap(libopen(opt2fn("-map",NFILE,fnm)),
		   opt2fn("-map",NFILE,fnm),&(mat.map));
  
  xvgr_legend(acct,2,leg);
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if (natoms != atoms->nr) 
    fatal_error(0,"Trajectory does not match topology!");
  
  nt=t;
  do {
    if (t >= nt) {
      fprintf(stderr,"\rt=%.2f",t);
      rm_pbc(&(top->idef),atoms->nr,box,x,x);
      write_pdb_conf_indexed(pdbfile,atoms,x,box,gnx,index);
      system(dssp);
      strip_dssp(tmpfile,nres,r0,bPhbres,t,dt,ss,acc,acct,&mat);
      remove(tmpfile);
      remove(pdbfile);
      nt+=dt;
      nframe++;
    }
  } while(read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  close_trj(status);
  fclose(ss);
  fclose(acc);
  fclose(acct);
  ss=ffopen(opt2fn("-ss",NFILE,fnm),"w");
  write_xpm_m(ss,mat);
  ffclose(ss);
  remove(tmpfile);
  remove(pdbfile);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
  xvgr_file(opt2fn("-ss",NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}
