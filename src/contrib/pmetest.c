/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "copyrite.h"
#include "main.h"
#include "nrnb.h"
#include "txtdump.h"
#include "tpxio.h"
#include "statutil.h"
#include "futil.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "mdatoms.h"
#include "coulomb.h"
#include "nsb.h"
#include "rmpbc.h"
#include "pme.h"
#include "force.h"
#include "xvgr.h"
#include "pbc.h"

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREAD_MPI
#include "tmpi.h"
#endif

#include "block_tx.h"

rvec *xptr=NULL;

static int comp_xptr(const void *a,const void *b)
{
  int  va,vb;
  real dx;
  
  va = *(int *)a;
  vb = *(int *)b;
  
  if ((dx = (xptr[va][XX] - xptr[vb][XX])) < 0)
    return -1;
  else if (dx > 0)
    return 1;
  else
    return 0;
}

static void do_my_pme(FILE *fp,real tm,gmx_bool bVerbose,t_inputrec *ir,
		      rvec x[],rvec xbuf[],rvec f[],
		      real charge[],real qbuf[],real qqbuf[],
		      matrix box,gmx_bool bSort,
		      t_commrec *cr,t_nsborder *nsb,t_nrnb *nrnb,
		      t_block *excl,real qtot,
		      t_forcerec *fr,int index[],FILE *fp_xvg,
		      int ngroups,unsigned short cENER[])
{
  real   ener,vcorr,q,xx,dvdl=0,vdip,vcharge;
  tensor vir,vir_corr,vir_tot;
  rvec   mu_tot[2];
  int    i,m,ii,ig,jg;
  real   **epme,*qptr;
  
  /* Initiate local variables */  
  fr->f_el_recip = f;
  clear_mat(vir);
  clear_mat(vir_corr);
  
  if (ngroups > 1) {
    fprintf(fp,"There are %d energy groups\n",ngroups);
    snew(epme,ngroups);
    for(i=0; (i<ngroups); i++)
      snew(epme[i],ngroups);
  }
    
  /* Put x is in the box, this part needs to be parallellized properly */
  /*put_atoms_in_box(box,nsb->natoms,x);*/
  /* Here sorting of X (and q) is done.
   * Alternatively, one could just put the atoms in one of the
   * cr->nnodes slabs. That is much cheaper than sorting.
   */
  for(i=0; (i<nsb->natoms); i++)
    index[i] = i;
  if (bSort) {
    xptr = x;
    qsort(index,nsb->natoms,sizeof(index[0]),comp_xptr);
    xptr = NULL; /* To trap unintentional use of the ptr */
  }
  /* After sorting we only need the part that is to be computed on 
   * this processor. We also compute the mu_tot here (system dipole)
   */
  clear_rvec(mu_tot[0]);
  for(i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) {
    ii      = index[i];
    q       = charge[ii];
    qbuf[i] = q;
    for(m=0; (m<DIM); m++) {
      xx         = x[ii][m];
      xbuf[i][m] = xx;
      mu_tot[0][m] += q*xx;
    }
    clear_rvec(f[ii]);
  }
  copy_rvec(mu_tot[0],mu_tot[1]);
  if (debug) {
    pr_rvec(debug,0,"qbuf",qbuf,nsb->natoms,TRUE);
    pr_rvecs(debug,0,"xbuf",xbuf,nsb->natoms);
    pr_rvecs(debug,0,"box",box,DIM);
  }
  for(ig=0; (ig<ngroups); ig++) {
    for(jg=ig; (jg<ngroups); jg++) {
      if (ngroups > 1) {
	for(i=START(nsb); (i<START(nsb)+HOMENR(nsb)); i++) {
	  if ((cENER[i] == ig) || (cENER[i] == jg))
	    qqbuf[i] = qbuf[i];
	  else
	    qqbuf[i] = 0;
	}
	qptr = qqbuf;
      }
      else
	qptr = qbuf;
      ener  = do_pme(fp,bVerbose,ir,xbuf,f,qptr,qptr,box,cr,
		     nsb,nrnb,vir,fr->ewaldcoeff,FALSE,0,&dvdl,FALSE);
      vcorr = ewald_LRcorrection(fp,nsb,cr,fr,qptr,qptr,excl,xbuf,box,mu_tot,
				 ir->ewald_geometry,ir->epsilon_surface,
				 0,&dvdl,&vdip,&vcharge);
      gmx_sum(1,&ener,cr);
      gmx_sum(1,&vcorr,cr);
      if (ngroups > 1)
	epme[ig][jg] = ener+vcorr;
    }
  }
  if (ngroups > 1) {
    if (fp_xvg) 
      fprintf(fp_xvg,"%10.3f",tm);
    for(ig=0; (ig<ngroups); ig++) {
      for(jg=ig; (jg<ngroups); jg++) {
	if (ig != jg)
	  epme[ig][jg] -= epme[ig][ig]+epme[jg][jg];
	if (fp_xvg) 
	  fprintf(fp_xvg,"  %12.5e",epme[ig][jg]);
      }
    }
    if (fp_xvg) 
      fprintf(fp_xvg,"\n");
  }
  else {
    fprintf(fp,"Time: %10.3f Energy: %12.5e  Correction: %12.5e  Total: %12.5e\n",
	    tm,ener,vcorr,ener+vcorr);
    if (fp_xvg) 
      fprintf(fp_xvg,"%10.3f %12.5e %12.5e %12.5e\n",tm,ener+vcorr,vdip,vcharge);
    if (bVerbose) {
      m_add(vir,vir_corr,vir_tot);
      gmx_sum(9,vir_tot[0],cr);
      pr_rvecs(fp,0,"virial",vir_tot,DIM); 
    }
    fflush(fp);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "The [TT]pmetest[tt] program tests the scaling of the PME code. When only given",
    "a [TT].tpr[tt] file it will compute PME for one frame. When given a trajectory",
    "it will do so for all the frames in the trajectory. Before the PME",
    "routine is called the coordinates are sorted along the X-axis.[PAR]",
    "As an extra service to the public the program can also compute",
    "long-range Coulomb energies for components of the system. When the",
    "[TT]-groups[tt] flag is given to the program the energy groups",
    "from the [TT].tpr[tt] file will be read, and half an energy matrix computed."
  };
  t_commrec    *cr,*mcr;
  static t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD  },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efLOG, "-g",      "pme",      ffWRITE },
    { efTRX, "-f",      NULL,       ffOPTRD },
    { efXVG, "-x",      "ener-pme", ffWRITE }
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  static gmx_bool bVerbose=FALSE;
  static gmx_bool bOptFFT=FALSE;
  static gmx_bool bSort=FALSE;
  static int  ewald_geometry=eewg3D;
  static int  nnodes=1;
  static int  pme_order=0;
  static rvec grid = { -1, -1, -1 };
  static real rc   = 0.0;
  static real dtol = 0.0;
  static gmx_bool bGroups = FALSE;
  static t_pargs pa[] = {
    { "-np",      FALSE, etINT, {&nnodes},
      "Number of nodes, must be the same as used for [TT]grompp[tt]" },
    { "-v",       FALSE, etBOOL,{&bVerbose},  
      "Be loud and noisy" },
    { "-sort",    FALSE, etBOOL,{&bSort},  
      "Sort coordinates. Crucial for domain decomposition." },
    { "-grid",    FALSE, etRVEC,{&grid},
      "Number of grid cells in X, Y, Z dimension (if -1 use from [TT].tpr[tt])" },
    { "-order",   FALSE, etINT, {&pme_order},
      "Order of the PME spreading algorithm" },
    { "-groups",  FALSE, etBOOL, {&bGroups},
      "Compute half an energy matrix based on the energy groups in your [TT].tpr[tt] file" },
    { "-rc",      FALSE, etREAL, {&rc},
      "Rcoulomb for Ewald summation" },
    { "-tol",     FALSE, etREAL, {&dtol},
      "Tolerance for Ewald summation" }
  };
  FILE        *fp;
  t_inputrec  *ir;
  t_topology  top;
  t_tpxheader tpx;
  t_nrnb      nrnb;
  t_nsborder  *nsb;
  t_forcerec  *fr;
  t_mdatoms   *mdatoms;
  char        title[STRLEN];
  int         natoms,step,status,i,ncg,root;
  real        t,lambda,ewaldcoeff,qtot;
  rvec        *x,*f,*xbuf;
  int         *index;
  gmx_bool        bCont;
  real        *charge,*qbuf,*qqbuf;
  matrix      box;
  
  /* Start the actual parallel code if necessary */
  cr   = init_par(&argc,&argv);
  root = 0;
  
  if (MASTER(cr)) 
    CopyRight(stderr,argv[0]);
  
  /* Parse command line on all processors, arguments are passed on in 
   * init_par (see above)
   */
  parse_common_args(&argc,argv,
		    PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_BE_NICE |
		    PCA_CAN_SET_DEFFNM | (MASTER(cr) ? 0 : PCA_QUIET),
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
#ifndef GMX_MPI
  if (nnodes > 1) 
    gmx_fatal(FARGS,"GROMACS compiled without MPI support - can't do parallel runs");
#endif

  /* Open log files on all processors */
  open_log(ftp2fn(efLOG,NFILE,fnm),cr);
  snew(ir,1);
  
  if (MASTER(cr)) {
    /* Read tpr file etc. */
    read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&tpx,FALSE,NULL,NULL);
    snew(x,tpx.natoms);
    read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,ir,
	     box,&natoms,x,NULL,NULL,&top);
    /* Charges */
    qtot = 0;
    snew(charge,natoms);
    for(i=0; (i<natoms); i++) {
      charge[i] = top.atoms.atom[i].q;
      qtot += charge[i];
    }
  
    /* Grid stuff */
    if (opt2parg_bSet("-grid",asize(pa),pa)) {
      ir->nkx = grid[XX];
      ir->nky = grid[YY];
      ir->nkz = grid[ZZ];
    }
    /* Check command line parameters for consistency */
    if ((ir->nkx <= 0) || (ir->nky <= 0) || (ir->nkz <= 0))
      gmx_fatal(FARGS,"PME grid = %d %d %d",ir->nkx,ir->nky,ir->nkz);
    if (opt2parg_bSet("-rc",asize(pa),pa)) 
      ir->rcoulomb = rc;
    if (ir->rcoulomb <= 0)
      gmx_fatal(FARGS,"rcoulomb should be > 0 (not %f)",ir->rcoulomb);
    if (opt2parg_bSet("-order",asize(pa),pa)) 
      ir->pme_order = pme_order;
    if (ir->pme_order <= 0)
      gmx_fatal(FARGS,"pme_order should be > 0 (not %d)",ir->pme_order);
    if (opt2parg_bSet("-tol",asize(pa),pa))
      ir->ewald_rtol = dtol;
    if (ir->ewald_rtol <= 0)
      gmx_fatal(FARGS,"ewald_tol should be > 0 (not %f)",ir->ewald_rtol);
  }
  else {
    init_top(&top);
  }

  /* Add parallellization code here */
  snew(nsb,1);
  if (MASTER(cr)) {
    ncg = top.blocks[ebCGS].multinr[0];
    for(i=0; (i<cr->nnodes-1); i++)
      top.blocks[ebCGS].multinr[i] = min(ncg,(ncg*(i+1))/cr->nnodes);
    for( ; (i<MAXNODES); i++)
      top.blocks[ebCGS].multinr[i] = ncg;
  }
  if (PAR(cr)) {
    /* Set some variables to zero to avoid core dumps */
    ir->opts.ngtc = ir->opts.ngacc = ir->opts.ngfrz = ir->opts.ngener = 0;
#ifdef GMX_MPI
    /* Distribute the data over processors */
    MPI_Bcast(&natoms,1,MPI_INT,root,MPI_COMM_WORLD);
    MPI_Bcast(ir,sizeof(*ir),MPI_BYTE,root,MPI_COMM_WORLD);
    MPI_Bcast(&qtot,1,GMX_MPI_REAL,root,MPI_COMM_WORLD);
#endif

    /* Call some dedicated communication routines, master sends n-1 times */
    if (MASTER(cr)) {
      for(i=1; (i<cr->nnodes); i++) {
	mv_block(i,&(top.blocks[ebCGS]));
	mv_block(i,&(top.atoms.excl));
      }
    }
    else {
      ld_block(root,&(top.blocks[ebCGS]));
      ld_block(root,&(top.atoms.excl));
    }
    if (!MASTER(cr)) {
      snew(charge,natoms);
      snew(x,natoms);
    }
#ifdef GMX_MPI
    MPI_Bcast(charge,natoms,GMX_MPI_REAL,root,MPI_COMM_WORLD);
#endif
  }
  ewaldcoeff = calc_ewaldcoeff(ir->rcoulomb,ir->ewald_rtol);
  
  
  if (bVerbose)
    pr_inputrec(stdlog,0,"Inputrec",ir);

  /* Allocate memory for temp arrays etc. */
  snew(xbuf,natoms);
  snew(f,natoms);
  snew(qbuf,natoms);
  snew(qqbuf,natoms);
  snew(index,natoms);

  /* Initialize the PME code */  
  init_pme(stdlog,cr,ir->nkx,ir->nky,ir->nkz,ir->pme_order,
	   natoms,FALSE,bOptFFT,ewald_geometry);
	   
  /* MFlops accounting */
  init_nrnb(&nrnb);
  
  /* Initialize the work division */
  calc_nsb(stdlog,&(top.blocks[ebCGS]),cr->nnodes,nsb,0);
  nsb->nodeid = cr->nodeid;
  print_nsb(stdlog,"pmetest",nsb);  

  /* Initiate forcerec */
  mdatoms = atoms2md(stdlog,&top.atoms,ir->opts.nFreeze,ir->eI,
		     ir->delta_t,0,ir->opts.tau_t,FALSE,FALSE);
  snew(fr,1);
  init_forcerec(stdlog,fr,ir,&top,cr,mdatoms,nsb,box,FALSE,NULL,NULL,FALSE);
  
  /* First do PME based on coordinates in tpr file, send them to
   * other processors if needed.
   */
  if (MASTER(cr))
    fprintf(stdlog,"-----\n"
	    "Results based on tpr file %s\n",ftp2fn(efTPX,NFILE,fnm));
#ifdef GMX_MPI
  if (PAR(cr)) {
    MPI_Bcast(x[0],natoms*DIM,GMX_MPI_REAL,root,MPI_COMM_WORLD);
    MPI_Bcast(box[0],DIM*DIM,GMX_MPI_REAL,root,MPI_COMM_WORLD);
    MPI_Bcast(&t,1,GMX_MPI_REAL,root,MPI_COMM_WORLD);
  }
#endif
  do_my_pme(stdlog,0,bVerbose,ir,x,xbuf,f,charge,qbuf,qqbuf,box,bSort,
	    cr,nsb,&nrnb,&(top.atoms.excl),qtot,fr,index,NULL,
	    bGroups ? ir->opts.ngener : 1,mdatoms->cENER);

  /* If we have a trajectry file, we will read the frames in it and compute
   * the PME energy.
   */
  if (ftp2bSet(efTRX,NFILE,fnm)) {
    fprintf(stdlog,"-----\n"
	    "Results based on trx file %s\n",ftp2fn(efTRX,NFILE,fnm));
    if (MASTER(cr)) {
      sfree(x);
      natoms = read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box); 
      if (natoms != top.atoms.nr)
	gmx_fatal(FARGS,"natoms in trx = %d, in tpr = %d",natoms,top.atoms.nr);
      fp = xvgropen(ftp2fn(efXVG,NFILE,fnm),"PME Energy","Time (ps)","E (kJ/mol)");
    }
    else
      fp = NULL;
    do {
      /* Send coordinates, box and time to the other nodes */
#ifdef GMX_MPI
      if (PAR(cr)) {
	MPI_Bcast(x[0],natoms*DIM,GMX_MPI_REAL,root,MPI_COMM_WORLD);
	MPI_Bcast(box[0],DIM*DIM,GMX_MPI_REAL,root,MPI_COMM_WORLD);
	MPI_Bcast(&t,1,GMX_MPI_REAL,root,MPI_COMM_WORLD);
      }
#endif
      rm_pbc(&top.idef,nsb->natoms,box,x,x);
      /* Call the PME wrapper function */
      do_my_pme(stdlog,t,bVerbose,ir,x,xbuf,f,charge,qbuf,qqbuf,box,bSort,cr,
		nsb,&nrnb,&(top.atoms.excl),qtot,fr,index,fp,
		bGroups ? ir->opts.ngener : 1,mdatoms->cENER);
      /* Only the master processor reads more data */
      if (MASTER(cr))
          bCont = read_next_x(status,&t,natoms,x,box);
      /* Check whether we need to continue */
#ifdef GMX_MPI
      if (PAR(cr))
          MPI_Bcast(&bCont,1,MPI_INT,root,MPI_COMM_WORLD);
#endif
      
    } while (bCont);
    
    /* Finish I/O, close files */
    if (MASTER(cr)) {
      close_trx(status);
      ffclose(fp);
    }
  }
  
  if (bVerbose) {
    /* Do some final I/O about performance, might be useful in debugging */
    fprintf(stdlog,"-----\n");
    print_nrnb(stdlog,&nrnb);
  }
  
  /* Finish the parallel stuff */  
  if (gmx_parallel_env_initialized())
    gmx_finalize(cr);

  /* Thank the audience, as usual */
  if (MASTER(cr)) 
    thanx(stderr);

  return 0;
}

