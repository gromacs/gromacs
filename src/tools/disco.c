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


#include "macros.h"
#include "statutil.h"
#include "pdbio.h"
#include "smalloc.h"
#include "random.h"
#include "vec.h"
#include "princ.h"
#include "confio.h"
#include "index.h"
#include "filenm.h"
#include "do_fit.h"
#include "tpxio.h"
#include "copyrite.h"
#include "disco.h"
#include "xvgr.h"
#include "main.h"
#include "network.h"
#ifdef GMX_MPI
#include "mpi.h"
#endif

void rand_box(bool bUserBox,
	      matrix box,rvec boxsize,int nres,bool bCubic,int *seed)
{
  int  m;
  real fac;
  
  clear_mat(box);
  
  if (bUserBox) {
    for(m=0; (m<DIM); m++)
      box[m][m] = boxsize[m];
  }
  else {
    /* Generate a random box with size between 5*nres and 10*nres in nm */
    fac = 0.5*nres; /* Ca-Ca distance is 0.35 nm */
    box[XX][XX] = fac*(1+rando(seed));
    if (bCubic)
      box[YY][YY] = box[ZZ][ZZ] = box[XX][XX];
    else {
      box[YY][YY] = fac*(1+rando(seed));
      box[ZZ][ZZ] = fac*(1+rando(seed));
    }
    for(m=0; (m<DIM); m++) 
      boxsize[m] = box[m][m];
  }
}

void rand_coord(rvec x,int *seed,rvec box)
{
  int m;
  
  for(m=0; (m<DIM); m++)
    x[m] = box[m]*rando(seed);
}

void rand_coords(int natom,rvec x[],rvec xref[],real weight[],
		 bool bCenter,rvec xcenter[],rvec box,int *seed)
{
  int i;
  
  for(i=0; (i<natom); i++) {
    if (weight[i] == 0) 
      copy_rvec(xref[i],x[i]);
    else {
      rand_coord(x[i],seed,box);
      if (bCenter)
	rvec_inc(x[i],xcenter[i]);
    }
  }
}

void pr_conv_stat(FILE *fp,int ntry,int nconv,double tnit)
{
  fprintf(fp,"\n-------------------------\n");
  fprintf(fp," Convergence statistics:\n");
  fprintf(fp," # tries:     %d\n",ntry);
  fprintf(fp," # converged: %d\n",nconv);
  fprintf(fp," # nit/ntry:  %d\n",(int)(tnit/ntry));
  if (nconv > 0)
    fprintf(fp," # nit/nconv: %d\n",(int)(tnit/nconv));
  fprintf(fp,"-------------------------\n");
}

static void do_disco(FILE *log,char *outfn,char *keepfn,t_correct *c,
		     bool bVerbose,t_atoms *atoms,
		     rvec xref[],rvec xcenter[],
		     int nstruct,int *seed,
		     bool bFit,int nfit,atom_id fit_ind[],
		     bool bPrintViol,char *violfn,rvec boxsize)
{
  FILE    *fp,*gp;
  int     *nconvdist;
  int     i,k,kk,nconv,ntry,status,kstatus,natom,nres,nit,nvtest,nviol;
  double  tnit;
  rvec    *x,xcm;
  matrix  box,wrbox;
  atom_id *wr_ind;
  real    *w_rls;
  bool    bConverged;

  natom  = atoms->nr;
  nres   = atoms->nres;
  
  /* Initiate some local arrays */
  init_corr2(c,natom);
  
  clear_mat(wrbox);
  wrbox[XX][XX] = wrbox[YY][YY] = wrbox[ZZ][ZZ] = nres;  
  status = open_trx(outfn,"w");
  if (keepfn)
    kstatus = open_trx(keepfn,"w");
  else
    kstatus = -1;
  snew(x,natom);
  snew(wr_ind,natom);
  for(k=0; (k<natom); k++)
    wr_ind[k]=k;
  snew(w_rls,natom);
  for(k=0; (k<nfit); k++)
    w_rls[fit_ind[k]] = 1;

  snew(nconvdist,c->maxnit+1);
  /* Now loop over structures */
  tnit = 0;
  for(k=nconv=ntry=0; (k<nstruct); ntry++) {
    if (bVerbose)
      fprintf(stderr,"\rTry: %d, Success: %d",ntry,nconv);
      
    /* Generate random box*/
    rand_box(c->bBox,box,boxsize,nres,c->bCubic,seed);
    
    /* Generate random coords */
    rand_coords(natom,x,xref,c->weight,c->bCenter,xcenter,boxsize,seed);
    
    /* Now correct the random coords */
    nviol = shake_coords(log,bVerbose,k,natom,xref,x,seed,box,c,&nit);
    bConverged = (nviol == 0);
    tnit += nit;

    if (bConverged)
      nconvdist[nit]++;
    
    nvtest = quick_check(bVerbose ? log : NULL,natom,x,box,c);
    fprintf(stderr,"Double checking: %d violations\n",nvtest);

    if (bConverged || keepfn) {
      center_in_box(natom,x,wrbox,x);
      if (bFit)
	do_fit(natom,w_rls,xref,x);
      write_trx(bConverged ? status : kstatus,
		natom,wr_ind,atoms,k,(real) k,wrbox,x,NULL);
	
      if (bConverged) 
	nconv++;
      
      k++;
    }
    if (bPrintViol) {
      /* Print structure coloured by the violations */
      if (!atoms->pdbinfo)
	snew(atoms->pdbinfo,natom);
      for(kk=0; (kk<natom); kk++)
	atoms->pdbinfo[kk].bfac = (real) c->bViol[kk];
      gp=ffopen(violfn,"w");
      write_pdbfile(gp,"Structure coloured by violation",atoms,x,box,0,-1);
      ffclose(gp);
    }
  }
  close_trx(status);
  if (keepfn)
    close_trx(kstatus);
  gp = xvgropen("conv_stat.xvg","Iterations per converged structure",
		"nit","N");
  for(i=0; (i<c->maxnit); i++)
    fprintf(gp,"%10d  %10d\n",i,nconvdist[i]);
  ffclose(gp);
  sfree(x);
  sfree(w_rls);
  sfree(wr_ind);
  sfree(nconvdist);
  
  pr_conv_stat(log,ntry,nconv,tnit);
  pr_conv_stat(stderr,ntry,nconv,tnit);
}
	      
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "disco reads a topology (tpr) file and runs distance geometry",
    "calculations based on the distances defined in the",
    "distance-restraints section of the topology. An appropriate tpr",
    "file may be generated by the cdist program.[PAR]",
    "The algorithm is the CONCOORD algorithm of De Groot et al.,",
    "which in turn is derived from the SHAKE alogrithm.[PAR]",
    "A parallel version of disco is under development whihc uses a",
    "master-slave approach. Slaves work asynchronously, and it is no",
    "problem when nodes are not equally fast, or when a node dies,",
    "unless it is the master node."
  };
  FILE        *fp,*dp;
  char        title[256];
  bool        bCenter=FALSE,bBox;
  t_commrec   *cr;
  t_atoms     atoms,newatoms;
  t_correct   *corr=NULL;
  rvec        xcm,*xref=NULL,*xcenter=NULL;
  matrix      box;
  real        t,lambda,tot_weight;
  int         i,nfit,step,natom;
  atom_id     *fit_ind;
  char        *grpname;
  
  static int  nstruct=10,maxnit=1000,seed=1997,nbcheck=1;
  static int  nstprint=1,nstranlist=0,ngrow=0;
  static bool bVerbose=TRUE,bCubic=FALSE,bWeight=FALSE,bLower=FALSE;
  static real lowdev=0.05;
  static bool bExplicit=FALSE,bChiral=TRUE,bFit=FALSE,bDump=FALSE,bPep=TRUE;
  static bool bRanlistFirst=TRUE;
  static rvec boxsize={ 2, 2, 2 };
  t_pargs pa[] = {
    { "-nf",    FALSE, etINT,     {&nstruct},
      "Number of structures to generate" },
    { "-nit",   FALSE, etINT,     {&maxnit},
      "Max number of iterations for a structure to converge" },
    { "-v",     FALSE, etBOOL,    {&bVerbose},
      "Be verbosive" },
    { "-chiral",   FALSE, etBOOL, {&bChiral},
      "Check chirality during disco-ing" },
    { "-pep",   FALSE,  etBOOL,   {&bPep},
      "Flip all cis-peptide bonds automatically to trans" },
    { "-lower", FALSE,  etBOOL,   {&bLower},
      "Use lower bounds only for nonbondeds." },
    { "-weighted", FALSE, etBOOL, {&bWeight},
      "Use weighted disco. The STX file must be a pdb file in this case and weights are read from the occupancy field" },
    { "-dump",     FALSE, etBOOL, {&bDump},
      "Dump the trajectory of the shaking to testX.xtc file where X is the structure number." },
    { "-cubic",    FALSE, etBOOL, {&bCubic},
      "Generate coordinates in a cubic box, rather than rectangular" },
    { "-explicit", FALSE, etBOOL, {&bExplicit},
      "Use explicit updating of positions if the sum of deviations is smaller than lowdev" },
    { "-fit",      FALSE, etBOOL, {&bFit},
      "Fit output structures to reference structure in tpx file" },
    { "-nbcheck",  FALSE, etINT,  {&nbcheck},
      "Check non-bonded interactions every N steps" },
    { "-nstprint", FALSE, etINT,  {&nstprint},
      "Print number of violations every N steps" },
    { "-ranlist",  FALSE, etINT,  {&nstranlist},
      "Update list order to avoid bias every n steps" },
    { "-ranlistfirst",  FALSE, etBOOL,  {&bRanlistFirst},
      "Randomize list once before shaking" },
    { "-lowdev",   FALSE, etREAL, {&lowdev},
      "Low deviation [Sum of distance deviation per atom in nm] beyond which nonbondeds are done every step" },
    { "-seed",     FALSE, etINT,  {&seed},
      "Seed for the random number generator" },
    { "-box",      FALSE, etRVEC, {boxsize},
      "Boxsize (nm) for generating random coordinates" },
    { "-grow",     FALSE, etINT,  {&ngrow},
      "Number of steps after which Van der Waals lower bounds grow from 0 to the real lower bounds. If this is 0 (default), the Van der Waals lower bounds are in effect from the beginning" }
    
  };
#define NPA asize(pa)
    
  t_filenm fnm[] = {
    { efLOG, "-g",     "disco", ffWRITE },
    { efSTX, "-f",      NULL,   ffREAD  },
    { efDAT, "-d",     "cdist", ffREAD },
    { efDAT, "-do",    "distout",  ffOPTWR },
    { efSTO, "-c",      NULL,   ffREAD  },
    { efSTO, "-center", NULL,   ffOPTRD },
    { efNDX, "-n",      NULL,   ffOPTRD },
    { efTRX, "-o",   "structs", ffWRITE },
    { efTRX, "-keep", "unconverged", ffOPTWR },
    { efPDB, "-viol",  "vvv",   ffOPTWR }
  };
#define NFILE asize(fnm)

  cr       = init_par(&argc,&argv);
  bVerbose = bVerbose && MASTER(cr);
  
  if (MASTER(cr)) {
    CopyRight(stderr,argv[0]);

    parse_common_args(&argc,argv,PCA_BE_NICE,
		      NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);
		      
    /* Copy arguments to correct structure */
    bCenter = opt2bSet("-center",NFILE,fnm);
    bBox    = opt2parg_bSet("-box",asize(pa),pa),    
    corr    = init_corr(maxnit,nstprint,nbcheck,nstranlist,ngrow,bExplicit,
			bChiral,bPep,bDump,lowdev,bLower,bRanlistFirst,
			bCubic,bBox,bCenter);
  }
  /* Open the log file */
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  
  if (MASTER(cr)) {
    CopyRight(stdlog,argv[0]);
    please_cite(stdlog,"Ryckaert77a");
    please_cite(stdlog,"DeGroot97a");
  }

  if (MASTER(cr)) {  
    /* Get number of atoms etc. */
    get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&natom);
  
    init_t_atoms(&atoms,natom,bWeight);
    snew(xref,natom);
    read_stx_conf(ftp2fn(efSTX,NFILE,fnm),title,&atoms,xref,NULL,box);
    
    /* Translate reference and xcenter coords to C.O.M. */
    sub_xcm(xref,natom,NULL,NULL,xcm,FALSE);
  
    snew(corr->weight,natom);
    tot_weight = 0;
    for(i=0; (i<natom); i++) {
      corr->weight[i] = bWeight ? atoms.pdbinfo[i].occup : 1;
      tot_weight += corr->weight[i];
    }
  
    fprintf(stderr,"Reading distances from %s\n",opt2fn("-d",NFILE,fnm));
    read_dist(stdlog,opt2fn("-d",NFILE,fnm),natom,corr);

    /* Dump a distance file if necessary */
    if (opt2bSet("-do",NFILE,fnm)) {
      dp = fopen(opt2fn("-do",NFILE,fnm),"w");
      pr_distances(dp,corr);
      fclose(dp);
    }
  
    /* Check distances */
    check_dist(stdlog,corr);

    /* Read index if necessary */
    if (bFit) {
      fprintf(stderr,"Select group for fitting output structures:\n");
      get_index(&atoms,ftp2fn_null(efNDX,NFILE,fnm),1,
		&nfit,&fit_ind,&grpname);
    }
    else {
      nfit = 0;
      fit_ind = NULL;
    }
  
    /* Read centers for generating coordinates (optional) */
    if (bCenter) {
      snew(xcenter,natom);
      init_t_atoms(&newatoms,natom,TRUE);
      read_stx_conf(opt2fn("-center",NFILE,fnm),title,
		    &newatoms,xcenter,NULL,box);
      free_t_atoms(&newatoms);
    
      for(i=0; (i<natom); i++)
	rvec_dec(xcenter[i],xcm);
    }

    /* 
     * define improper dihedrals that are not automatically correct
     * when all distances are correct
     */
    define_impropers(stdlog,&atoms,corr);

    /* define peptide-bonds, so we can correct cis to trans
     * Adam Kirrander 990121
     */
    define_peptide_bonds(stdlog,&atoms,corr);

    /* Print parameters */
    pr_corr(stdlog,corr);
  }

  /* Now do my thing */
#ifdef GMX_MPI
  if (PAR(cr)) {
    if (MASTER(cr))
      disco_master(cr,stdlog,opt2fn("-o",NFILE,fnm),
		   opt2bSet("-keep",NFILE,fnm) ? opt2fn("-keep",NFILE,fnm) : NULL,
		   corr,bVerbose,&atoms,
		   xref,xcenter,nstruct,&seed,bFit,nfit,fit_ind,
		   opt2bSet("-viol",NFILE,fnm),opt2fn("-viol",NFILE,fnm),boxsize);
    else
      disco_slave(cr,stdlog);
  }  
  else {
#endif
    do_disco(stdlog,opt2fn("-o",NFILE,fnm),
	     opt2bSet("-keep",NFILE,fnm) ? opt2fn("-keep",NFILE,fnm) : NULL,
	     corr,bVerbose,&atoms,
	     xref,xcenter,nstruct,&seed,bFit,nfit,fit_ind,
	     opt2bSet("-viol",NFILE,fnm),opt2fn("-viol",NFILE,fnm),
	     boxsize);
#ifdef GMX_MPI
  }
#endif
  ffclose(stdlog);
#ifdef GMX_MPI
  if (PAR(cr))
    MPI_Finalize();
#endif
  if (MASTER(cr))
    thanx(stderr);
  
  return 0;  
}
