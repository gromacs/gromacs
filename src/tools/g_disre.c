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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_disre_c = "$Id$";

#include <math.h>
#include <string.h>
#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "mshift.h"
#include "xvgr.h"
#include "smalloc.h"
#include "nrnb.h"
#include "bondf.h"
#include "disre.h"
#include "statutil.h"
#include "force.h"
#include "gstat.h"
#include "main.h"
#include "genhydro.h"
#include "pdbio.h"
#include "rdgroup.h"
#include "mdatoms.h"

typedef struct {
  int n;
  real v;
} t_toppop;

t_toppop *top=NULL;
int      ntop=0;

static void init5(int n)
{
  ntop=n;
  snew(top,ntop);
}

static void reset5(void)
{
  int i;

  for(i=0; (i<ntop); i++) {
    top[i].n=-1;
    top[i].v= 0;
  }
}

int tpcomp(const void *a,const void *b)
{
  t_toppop *tpa;
  t_toppop *tpb;

  tpa=(t_toppop *)a;
  tpb=(t_toppop *)b;

  return  (1e7*(tpb->v-tpa->v));
}

static void add5(int ndr,real viol)
{
  int i,mini;
  
  mini=0;
  for(i=1; (i<ntop); i++)
    if (top[i].v < top[mini].v) 
      mini=i;
  if (viol > top[mini].v) {
    top[mini].v=viol;
    top[mini].n=ndr;
  }
}

static void print5(FILE *fp)
{
  int i;

  qsort(top,ntop,sizeof(top[0]),tpcomp);
  fprintf(fp,"Index:");
  for(i=0; (i<ntop); i++)
    fprintf(fp," %6d",top[i].n);
  fprintf(fp,"\nViol: ");
  for(i=0; (i<ntop); i++)
    fprintf(fp," %6.3f",top[i].v);
  fprintf(fp,"\n");
}

void check_viol(FILE *log,
		t_ilist *bonds,t_iparams forceparams[],
		t_functype functype[],
		int natoms,rvec x[],rvec f[],
		t_forcerec *fr,matrix box,t_graph *g,
		real *sumv,real *averv,
		real *maxv,int *nv,
		int isize,atom_id index[],real vvindex[])
{
  t_iatom *forceatoms;
  int     i,j,type,ftype,nat,nviol,ndr;
  real    mviol,tviol,viol,lam,dvdl;
  
  lam  =0;
  dvdl =0;
  tviol=0;
  nviol=0;
  mviol=0;
  ndr=0;
  reset5();
  forceatoms=bonds->iatoms;
  for(j=0; (j<isize); j++) {
    vvindex[j]=0;
  }
  for(i=0; (i<bonds->nr); ) {
    type=forceatoms[i++];
    ftype=functype[type];
    viol=interaction_function[ftype].ifunc(log,bonds->nr,&forceatoms[i],
					   &forceparams[type],
					   x,f,fr,g,box,lam,&dvdl,
					   NULL,0,NULL,NULL);
    if (viol > 0) {
      nviol++;
      add5(forceparams[type].disres.index,viol);
      if (viol > mviol) 
	mviol=viol;
      tviol+=viol;
      for(j=0; (j<isize); j++) {
	if (index[j] == forceparams[type].disres.index)
	  vvindex[j]=viol;
	}
    }
    ndr ++;
    nat  = interaction_function[ftype].nratoms;
    i   += nat;
  }
  *nv   = nviol;
  *maxv = mviol;
  *sumv = tviol;
  *averv= tviol/ndr;
  
  print5(log);
}

void patch_viol(t_ilist *bonds,t_iparams forceparams[],
		t_functype functype[])
{
  t_iatom *forceatoms;
  int     i,j,type,ftype,nat;
  real    viol2;
  
  viol2=0;
  forceatoms=bonds->iatoms;
  for(i=j=0; (i<bonds->nr); ) {
    type=forceatoms[i++];
    ftype=functype[type];
    if (ftype == F_DISRES)
      forceparams[ftype].disres.index=j++;
    nat=interaction_function[ftype].nratoms;
    i+=nat;
  }
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_disre computes violations of distance restraints. If necessary",
    "all protons can be added to a protein molecule. The program allways",
    "computes the instantaneous violations rather than time-averaged,",
    "because this analysis is done from a trajectory file afterwards",
    "it does not make sense to use time averaging.[PAR]",
    "An index file may be used to select out specific restraints for",
    "printing."
  };
  static bool bProt=FALSE;
  static int  ntop = 6;
  t_pargs pa[] = {
    { "-prot", FALSE, etBOOL, &bProt,
      "Protonate protein every step. This currently does not add terminal hydrogens, and therefore works only when the termini are capped." },
    { "-ntop", FALSE, etINT,  &ntop,
      "Number of large violations that are stored in the log file every step" }
  };
  
  FILE        *out,*aver,*numv,*maxxv,*xvg;
  t_inputrec  ir;
  t_topology  top;
  t_atoms     *atoms=NULL;
  t_forcerec  *fr;
  t_nrnb      nrnb;
  t_commrec   *cr;
  t_graph     *g;
  int         status,idum,natoms,i,j,nv;
  real        rdum,t,sumv,averv,maxv;
  rvec        *x,*f;
  real        epot[F_NRE];
  matrix      box;
  int         isize;
  atom_id     *index=NULL;
  char        *grpname;
  char        **leg;
  real        *vvindex;
  t_mdatoms   *mdatoms;
  
  t_filenm fnm[] = {
    { efTPB, NULL, NULL, ffREAD },
    { efTRX, "-f", NULL, ffREAD },
    { efXVG, "-ds", "drsum",  ffWRITE },
    { efXVG, "-da", "draver", ffWRITE },
    { efXVG, "-dn", "drnum",  ffWRITE },
    { efXVG, "-dm", "drmax",  ffWRITE },
    { efXVG, "-dr", "restr",  ffWRITE },
    { efLOG, "-l", "disres", ffWRITE },
    { efNDX, NULL, "viol",   ffOPTRD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  init5(ntop);
		    
  fprintf(stderr,"Reading status: %s\n",
	  read_status(ftp2fn(efTPB,NFILE,fnm),&idum,&rdum,&rdum,
		      &ir,NULL,NULL,NULL,
		      &idum,NULL,NULL,NULL,&idum,NULL,&top));
  g=mk_graph(&top.idef,top.atoms.nr,0);
  cr=init_par(NULL);
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (ftp2bSet(efNDX,NFILE,fnm)) {
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
    xvg=xvgropen(opt2fn("-dr",NFILE,fnm),"Inidividual Restraints","Time (ps)",
		 "nm");
    snew(vvindex,isize);
    snew(leg,isize);
    for(i=0; (i<isize); i++) {
      snew(leg[i],12);
      sprintf(leg[i],"index %d",index[i]);
    }
    xvgr_legend(xvg,isize,leg);
  }
  else 
    isize=0;
  
  ir.dr_tau=0.0;
  init_disres(stdlog,top.idef.il[F_DISRES].nr,&ir);

  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(f,5*natoms);
		
  out =xvgropen(opt2fn("-ds",NFILE,fnm),"Sum of Violations","Time (ps)","nm");
  aver=xvgropen(opt2fn("-da",NFILE,fnm),"Average Violation","Time (ps)","nm");
  numv=xvgropen(opt2fn("-dn",NFILE,fnm),"# Violations","Time (ps)","#");
  maxxv=xvgropen(opt2fn("-dm",NFILE,fnm),"Largest Violation","Time (ps)","nm");

  snew(atoms,1);
  atoms->nr=top.atoms.nr;
  atoms->nres=top.atoms.nres;
  snew(atoms->atomname,atoms->nr);
  snew(atoms->resname,atoms->nres);
  snew(atoms->atom,atoms->nr);
  memcpy(atoms->atom,top.atoms.atom,atoms->nr*sizeof(atoms->atom[0]));
  memcpy(atoms->atomname,top.atoms.atomname,atoms->nr*sizeof(atoms->atomname[0]));
  memcpy(atoms->resname,top.atoms.resname,atoms->nres*sizeof(atoms->resname[0]));

  mdatoms=atoms2md(&top.atoms,FALSE);  
  fr=mk_forcerec();
  fprintf(stdlog,"Made forcerec...\n");
  init_forcerec(stdlog,fr,&(ir),&(top.blocks[ebMOLS]),cr,
		&(top.blocks[ebCGS]),&(top.idef),mdatoms,box,FALSE);
  init_nrnb(&nrnb);
  j=0;
  do {
    if ((j % 10) == 0)
      fprintf(stderr,"\rFrame: %d",j);
    rm_pbc(&top.idef,top.atoms.nr,box,x,x);

    if (bProt) {
      protonate(&atoms,&x);
    }
    
    check_viol(stdlog,
	       &(top.idef.il[F_DISRES]),
	       top.idef.iparams,top.idef.functype,
	       atoms->nr,x,f,fr,box,g,&sumv,&averv,&maxv,&nv,
	       isize,index,vvindex);
    if (isize > 0) {
      fprintf(xvg,"%10g",t);
      for(i=0; (i<isize); i++)
	fprintf(xvg,"  %10g",vvindex[i]);
      fprintf(xvg,"\n");
    }    
    fprintf(out,  "%10g  %10g\n",t,sumv);
    fprintf(aver, "%10g  %10g\n",t,averv);
    fprintf(maxxv,"%10g  %10g\n",t,maxv);
    fprintf(numv, "%10g  %10d\n",t,nv);

    j++;
  } while (read_next_x(status,&t,natoms,x,box));
  
  close_trj(status);
  fclose(out);
  fclose(aver);
  fclose(numv);
  fclose(maxxv);
  if (isize > 0) {
    fclose(xvg);
    xvgr_file(opt2fn("-dr",NFILE,fnm),"-nxy");
  }
  xvgr_file(opt2fn("-dn",NFILE,fnm),NULL);
  xvgr_file(opt2fn("-da",NFILE,fnm),NULL);
  xvgr_file(opt2fn("-ds",NFILE,fnm),NULL);
  xvgr_file(opt2fn("-dm",NFILE,fnm),NULL);
  
  thanx(stdout);

#ifdef USE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}


