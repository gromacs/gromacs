/*
 * $Id$
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
#include <string.h>
#include <assert.h>

#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "mshift.h"
#include "xvgr.h"
#include "vec.h"
#include "do_fit.h"
#include "confio.h"
#include "smalloc.h"
#include "nrnb.h"
#include "disre.h"
#include "statutil.h"
#include "force.h"
#include "gstat.h"
#include "main.h"
#include "pdbio.h"
#include "index.h"
#include "mdatoms.h"
#include "nsb.h"
#include "tpxio.h"
#include "init.h"
#include "names.h"

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

static void check_viol(FILE *log,t_commrec *mcr,
		       t_ilist *disres,t_iparams forceparams[],
		       t_functype functype[],
		       rvec x[],rvec f[],
		       t_forcerec *fr,matrix box,t_graph *g,
		       real *sumv,real *averv,
		       real *maxv,int *nv,
		       int isize,atom_id index[],real vvindex[],
		       t_fcdata *fcd,real aver1[],real aver2[],real aver_3[])
{
  t_iatom *forceatoms;
  int     i,j,nat,n,type,nviol,ndr,label;
  real    ener,rt,mviol,tviol,viol,lam,dvdl;
  static  bool bFirst=TRUE;
  
  lam  =0;
  dvdl =0;
  tviol=0;
  nviol=0;
  mviol=0;
  ndr=0;
  if (ntop)
    reset5();
  forceatoms=disres->iatoms;
  for(j=0; (j<isize); j++) {
    vvindex[j]=0;
  }
  nat = interaction_function[F_DISRES].nratoms+1; 
  for(i=0; (i<disres->nr); ) {
    type  = forceatoms[i];
    n     = 0;
    label = forceparams[type].disres.label;
    if (debug) 
      fprintf(debug,"DISRE: ndr = %d, label = %d  i=%d, n =%d\n",
	      ndr,label,i,n);
    if (ndr != label) 
      fatal_error(0,"tpr inconsistency. ndr = %d, label = %d\n",ndr,label);
    do {
      n += nat;
    } while (((i+n) < disres->nr) && 
	     (forceparams[forceatoms[i+n]].disres.label == label));
    
    calc_disres_R_6(mcr,n,&forceatoms[i],forceparams,x,fr->ePBC==epbcFULL,fcd);

    rt = pow(fcd->disres.Rt_6[0],-1.0/6.0);
    aver1[ndr]  += rt;
    aver2[ndr]  += sqr(rt);
    aver_3[ndr] += pow(rt,-3.0);
    
    ener=interaction_function[F_DISRES].ifunc(n,&forceatoms[i],
					   forceparams,
					   x,f,fr,g,box,lam,&dvdl,
					   NULL,0,NULL,NULL,fcd);
    viol = fcd->disres.sumviol;
    
    
    if (viol > 0) {
      nviol++;
      if (ntop)
	add5(forceparams[type].disres.label,viol);
      if (viol > mviol) 
	mviol = viol;
      tviol += viol;
      for(j=0; (j<isize); j++) {
	if (index[j] == forceparams[type].disres.label)
	  vvindex[j]=pow(fcd->disres.Rt_6[0],-1.0/6.0);
	}
    }
    ndr ++;
    i   += n;
  }
  *nv   = nviol;
  *maxv = mviol;
  *sumv = tviol;
  *averv= tviol/ndr;
  
  if (bFirst) {
    fprintf(stderr,"\nThere are %d restraints and %d pairs\n",ndr,disres->nr/nat);
    bFirst = FALSE;
  }
  if (ntop)
    print5(log);
}

typedef struct {
  int  index;
  bool bCore;
  real up1,r,rT,viol,violT;
} t_dr_stats;

static int drs_comp(const void *a,const void *b)
{
  t_dr_stats *da,*db;
  
  da = (t_dr_stats *)a;
  db = (t_dr_stats *)b;
  
  if (da->viol > db->viol)
    return -1;
  else if (da->viol < db->viol)
    return 1;
  else
    return 0;
}

static void dump_dump(FILE *log,int ndr,t_dr_stats drs[])
{
  static char *core[] = { "All restraints", "Core restraints" };
  static char *tp[]   = { "linear", "third power" };
  real viol_tot,viol_max,viol;
  bool bCore,bTP;
  int  nviol,nrestr;
  int  i;

  for(bCore = FALSE; (bCore <= TRUE); bCore++) {
    for(bTP = FALSE; (bTP <= TRUE); bTP++) {
      viol_tot  = 0;
      viol_max  = 0;
      nviol     = 0;
      nrestr    = 0;
      for(i=0; (i<ndr); i++) {
	if (!bCore || (bCore && drs[i].bCore)) {
	  viol = bTP ? drs[i].violT : drs[i].viol;
	  viol_max     = max(viol_max,viol);
	  if (viol > 0)
	    nviol++;
	  viol_tot  += viol;
	  nrestr++;
	}
      }
      if ((nrestr > 0) || (bCore && (nrestr < ndr))) {
	fprintf(log,"\n");
	fprintf(log,"+++++++ %s ++++++++\n",core[bCore]);
	fprintf(log,"+++++++ Using %s averaging: ++++++++\n",tp[bTP]);
	fprintf(log,"Sum of violations: %8.3f nm\n",viol_tot);
	fprintf(log,"Average violation: %8.3f nm\n",viol_tot/nrestr);
	fprintf(log,"Largest violation: %8.3f nm\n",viol_max);
	fprintf(log,"Number of violated restraints: %d/%d\n",nviol,nrestr);
      }
    }
  }
}

static void dump_viol(FILE *log,int ndr,t_dr_stats *drs,bool bViol)
{
  int i;
  
  fprintf(log," Restraint Core     Up1     <r>    <rT>  <viol> <violT>\n");
  for(i=0; (i<ndr); i++) {
    if (bViol  && (drs[i].viol == 0) && (drs[i].violT > 0))
      break;
    fprintf(log,"%10d%5s%8.3f%8.3f%8.3f%8.3f%8.3f\n",
	    drs[i].index,yesno_names[drs[i].bCore],
	    drs[i].up1,drs[i].r,drs[i].rT,
	    drs[i].viol,drs[i].violT);
  }
}

static bool is_core(int i,int isize,atom_id index[])
{
  int kk;
  bool bIC = FALSE;
  
  for(kk=0; !bIC && (kk<isize); kk++)
    bIC = (index[kk] == i);
    
  return bIC;
}	      

static void dump_stats(FILE *log,int nsteps,int ndr,t_ilist *disres,
		       t_iparams ip[],real aver1[],real aver2[],
		       real aver_3[],int isize,atom_id index[],
		       t_atoms *atoms)
{
  int  i,j,nra;
  t_dr_stats *drs;

  fprintf(log,"\n");
  fprintf(log,"++++++++++++++ STATISTICS ++++++++++++++++++++++++\n");  
  snew(drs,ndr);
  j         = 0;
  nra       = interaction_function[F_DISRES].nratoms+1;
  for(i=0; (i<ndr); i++) {
    /* Search for the appropriate restraint */
    for( ; (j<disres->nr); j+=nra) {
      if (ip[disres->iatoms[j]].disres.label == i)
	break;
    }
    drs[i].index = i;
    drs[i].bCore = is_core(i,isize,index);
    drs[i].up1   = ip[disres->iatoms[j]].disres.up1;
    drs[i].r     = aver1[i]/nsteps;
    drs[i].rT    = pow(aver_3[i]/nsteps,-1.0/3.0);
    drs[i].viol  = max(0,drs[i].r-drs[i].up1);
    drs[i].violT = max(0,drs[i].rT-drs[i].up1);
    if (atoms) {
      int j1 = disres->iatoms[j+1];
      int j2 = disres->iatoms[j+2];
      atoms->pdbinfo[j1].bfac += drs[i].violT*5;
      atoms->pdbinfo[j2].bfac += drs[i].violT*5;
    }
  }
  dump_viol(log,ndr,drs,FALSE);
  
  fprintf(log,"+++ Sorted by linear averaged violations: +++\n");
  qsort(drs,ndr,sizeof(drs[0]),drs_comp);
  dump_viol(log,ndr,drs,TRUE);
	      
  dump_dump(log,ndr,drs);
  
  sfree(drs);
}

int gmx_disre(int argc,char *argv[])
{
  static char *desc[] = {
    "g_disre computes violations of distance restraints.",
    "If necessary all protons can be added to a protein molecule ",
    "using the protonate program.[PAR]",
    "The program allways",
    "computes the instantaneous violations rather than time-averaged,",
    "because this analysis is done from a trajectory file afterwards",
    "it does not make sense to use time averaging. However,",
    "the time averaged values per restraint are given in the log file.[PAR]",
    "An index file may be used to select specific restraints for",
    "printing.[PAR]",
    "When the optional[TT]-q[tt] flag is given a pdb file coloured by the",
    "amount of average violations."
  };
  static int  ntop      = 0;
  t_pargs pa[] = {
    { "-ntop", FALSE, etINT,  {&ntop},
      "Number of large violations that are stored in the log file every step" }
  };
  
  FILE        *out,*aver,*numv,*maxxv,*xvg=NULL;
  t_tpxheader header;
  t_inputrec  ir;
  t_topology  top;
  rvec        *xtop;
  t_atoms     *atoms=NULL;
  t_forcerec  *fr;
  t_fcdata    *fcd;
  t_nrnb      nrnb;
  t_nsborder  *nsb;
  t_commrec   *cr;
  t_graph     *g;
  int         status,ntopatoms,natoms,i,j,nv,step,kkk,m;
  real        t,sumv,averv,maxv,lambda;
  real        *aver1,*aver2,*aver_3;
  rvec        *x,*f,*xav=NULL;
  matrix      box;
  bool        bPDB;
  int         isize;
  atom_id     *index=NULL,*ind_fit=NULL;
  char        *grpname;
  char        **leg;
  real        *vvindex=NULL,*w_rls=NULL;
  t_mdatoms   *mdatoms;
  
  t_filenm fnm[] = {
    { efTPX, NULL, NULL, ffREAD },
    { efTRX, "-f", NULL, ffREAD },
    { efXVG, "-ds", "drsum",  ffWRITE },
    { efXVG, "-da", "draver", ffWRITE },
    { efXVG, "-dn", "drnum",  ffWRITE },
    { efXVG, "-dm", "drmax",  ffWRITE },
    { efXVG, "-dr", "restr",  ffWRITE },
    { efLOG, "-l",  "disres", ffWRITE },
    { efNDX, NULL,  "viol",   ffOPTRD },
    { efPDB, "-q",  "viol",   ffOPTWR }
  };
#define NFILE asize(fnm)

  cr  = init_par(&argc,&argv);
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
		    
  if (ntop)
    init5(ntop);
  
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&header,FALSE,NULL,NULL);
  snew(xtop,header.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,
	   box,&ntopatoms,xtop,NULL,NULL,&top);
  bPDB = opt2bSet("-q",NFILE,fnm);
  if (bPDB) {
    snew(xav,ntopatoms);
    snew(ind_fit,ntopatoms);
    snew(w_rls,ntopatoms);
    for(kkk=0; (kkk<ntopatoms); kkk++) {
      w_rls[kkk] = 1;
      ind_fit[kkk] = kkk;
    }   
    if (top.atoms.pdbinfo == NULL)
      snew(top.atoms.pdbinfo,ntopatoms);
  } 
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);

  check_nnodes_top(ftp2fn(efTPX,NFILE,fnm),&top,1);

  g   = mk_graph(&top.idef,top.atoms.nr,FALSE,FALSE);  
  
  if (ftp2bSet(efNDX,NFILE,fnm)) {
    rd_index(ftp2fn(efNDX,NFILE,fnm),1,&isize,&index,&grpname);
    xvg=xvgropen(opt2fn("-dr",NFILE,fnm),"Inidividual Restraints","Time (ps)",
		 "nm");
    snew(vvindex,isize);
    snew(leg,isize);
    for(i=0; (i<isize); i++) {
      index[i]++;
      snew(leg[i],12);
      sprintf(leg[i],"index %d",index[i]);
    }
    xvgr_legend(xvg,isize,leg);
  }
  else 
    isize=0;
  
  snew(fcd,1);
  ir.dr_tau=0.0;
  init_disres(stdlog,top.idef.il[F_DISRES].nr,top.idef.il[F_DISRES].iatoms,
	      top.idef.iparams,&ir,NULL,fcd);

  snew(aver1,fcd->disres.nr);
  snew(aver2,fcd->disres.nr);
  snew(aver_3,fcd->disres.nr);
	      
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  snew(f,5*natoms);
		
  out =xvgropen(opt2fn("-ds",NFILE,fnm),"Sum of Violations","Time (ps)","nm");
  aver=xvgropen(opt2fn("-da",NFILE,fnm),"Average Violation","Time (ps)","nm");
  numv=xvgropen(opt2fn("-dn",NFILE,fnm),"# Violations","Time (ps)","#");
  maxxv=xvgropen(opt2fn("-dm",NFILE,fnm),"Largest Violation","Time (ps)","nm");

  snew(nsb,1);
  snew(atoms,1);
  atoms->nr=top.atoms.nr;
  atoms->nres=top.atoms.nres;
  snew(atoms->atomname,atoms->nr);
  snew(atoms->resname,atoms->nres);
  snew(atoms->atom,atoms->nr);
  memcpy(atoms->atom,top.atoms.atom,atoms->nr*sizeof(atoms->atom[0]));
  memcpy(atoms->atomname,top.atoms.atomname,
	 atoms->nr*sizeof(atoms->atomname[0]));
  memcpy(atoms->resname,top.atoms.resname,
	 atoms->nres*sizeof(atoms->resname[0]));
  mdatoms = atoms2md(stdlog,&top.atoms,ir.opts.nFreeze,
		     FALSE,0,0,NULL,FALSE,FALSE);  
  fr      = mk_forcerec();
  fprintf(stdlog,"Made forcerec...\n");
  calc_nsb(stdlog,&(top.blocks[ebCGS]),1,nsb,0);
  init_forcerec(stdlog,fr,&ir,&top,cr,mdatoms,nsb,box,FALSE,NULL,FALSE);
  init_nrnb(&nrnb);
  j=0;
  do {
    rm_pbc(&top.idef,natoms,box,x,x);

    check_viol(stdlog,cr,
	       &(top.idef.il[F_DISRES]),
	       top.idef.iparams,top.idef.functype,
	       x,f,fr,box,g,&sumv,&averv,&maxv,&nv,
	       isize,index,vvindex,fcd,aver1,aver2,aver_3);
    if (bPDB) {
      reset_x(top.atoms.nr,ind_fit,top.atoms.nr,NULL,x,w_rls);
      do_fit(top.atoms.nr,w_rls,x,x);
      if (j == 0) {
	/* Store the first frame of the trajectory as 'characteristic'
	 * for colouring with violations.
	 */
	for(kkk=0; (kkk<top.atoms.nr); kkk++)
	  copy_rvec(x[kkk],xav[kkk]);
      }
    }
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
  dump_stats(stdlog,j,fcd->disres.nr,&(top.idef.il[F_DISRES]),
	     top.idef.iparams,aver1,aver2,aver_3,isize,index,
	     bPDB ? (&top.atoms) : NULL);
  if (bPDB) {
    write_sto_conf(opt2fn("-q",NFILE,fnm),
		   "Coloured by average violation in Angstrom",
		   &(top.atoms),xav,NULL,box);
  }
  fclose(out);
  fclose(aver);
  fclose(numv);
  fclose(maxxv);
  if (isize > 0) {
    fclose(xvg);
    do_view(opt2fn("-dr",NFILE,fnm),"-nxy");
  }
  do_view(opt2fn("-dn",NFILE,fnm),"-nxy");
  do_view(opt2fn("-da",NFILE,fnm),"-nxy");
  do_view(opt2fn("-ds",NFILE,fnm),"-nxy");
  do_view(opt2fn("-dm",NFILE,fnm),"-nxy");
  
  thanx(stderr);

  if (gmx_parallel_env)
    gmx_finalize();
  
  return 0;
}
