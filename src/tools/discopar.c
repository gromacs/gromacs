#include "typedefs.h"
#include "network.h"
#include "smalloc.h"
#include "vec.h"
#include "statutil.h"
#include "do_fit.h"
#include "random.h"
#include "futil.h"
#include "xvgr.h"
#include "pdbio.h"
#include "disco.h"
	
static t_correct *recv_init(t_commrec *cr,int *seed,int *natom,int *nres,
			    rvec **xref,rvec **xcenter,bool *bKeep)
{
  t_correct *c;
  
  gmx_rxs(0,seed,sizeof(*seed));
  gmx_rxs(0,natom,sizeof(*natom));
  gmx_rxs(0,nres,sizeof(*nres));
  gmx_rxs(0,bKeep,sizeof(*bKeep));
  
  snew(*xref,*natom);
  snew(*xcenter,*natom);
  gmx_rxs(0,*xref,*natom*sizeof(**xref));
  gmx_rxs(0,*xcenter,*natom*sizeof(**xcenter));

  snew(c,1);
#define cget(nn) gmx_rxs(0,record(c->nn))
  cget(maxnit);
  cget(nbcheck);
  cget(nstprint);
  cget(nstranlist);
  cget(ngrow);
  cget(bExplicit);
  cget(bChiral);
  cget(bPep);
  cget(bLowerOnly);
  cget(bRanlistFirst);
  cget(bCubic);
  cget(bBox);
  cget(bCenter);
  cget(lodev);
  cget(maxdist);
  cget(ndist);
  cget(npep);
  cget(nimp);
#undef cget
  /* Get the important input data */
  snew(c->d,c->ndist);
  gmx_rxs(0,array(c->d,c->ndist));
  
  snew(c->pepbond,c->npep);
  gmx_rxs(0,array(c->pepbond,c->npep));
  
  snew(c->imp,c->nimp+1);
  gmx_rxs(0,array(c->imp,c->nimp+1));
  
  /* Other arrays can be deduced from these */
  
  return c;
}

static void send_init(t_commrec *cr,t_correct *c,int *seed,int natom,int nres,
		      rvec *xref,rvec *xcenter,bool bKeep)
{
  int pid;
  
  for(pid=1; (pid < cr->nprocs); pid++) {
    gmx_txs(pid,seed,sizeof(*seed));
    /* Update random seed */
    (void) rando(seed);
    gmx_txs(pid,record(natom));
    gmx_txs(pid,record(nres));
    gmx_txs(pid,record(bKeep));
    
    gmx_txs(pid,array(xref,natom));
    gmx_txs(pid,array(xcenter,natom));
    
#define cput(nn) gmx_txs(pid,record(c->nn))
    cput(maxnit);
    cput(nbcheck);
    cput(nstprint);
    cput(nstranlist);
    cput(ngrow);
    cput(bExplicit);
    cput(bChiral);
    cput(bPep);
    cput(bLowerOnly);
    cput(bRanlistFirst);
    cput(bCubic);
    cput(bBox);
    cput(bCenter);
    cput(lodev);
    cput(maxdist);
    cput(ndist);
    cput(npep);
    cput(nimp);
#undef cput
    /* Send the important input data */
    gmx_txs(pid,array(c->d,c->ndist));
  
    gmx_txs(pid,array(c->pepbond,c->npep));
  
    gmx_txs(pid,array(c->imp,c->nimp+1));
  
    /* Other arrays can be deduced from these */
  }
}

static void send_coords(t_commrec *cr,bool bConverged,int nit,int k,
			int natom,rvec x[],matrix box)
{
  gmx_txs(0,record(cr->pid));
  gmx_txs(0,record(bConverged));
  gmx_txs(0,record(nit));
  gmx_txs(0,record(k));
  gmx_txs(0,record(natom));
  gmx_txs(0,array(x,natom));
  gmx_txs(0,array(box,DIM));
}

static int recv_coords(t_commrec *cr,bool *bConverged,int *nit,int *k,
		       rvec x[],matrix box)
{
  MPI_Status status;
  int        pid,natom;
  
  /* Check whether there is something from anyone */
  MPI_Recv(&pid,           /* message buffer */
	   1,              /* one data item */
	   MPI_INT,        /* of type double real */
	   MPI_ANY_SOURCE, /* receive from any sender */
	   MPI_ANY_TAG,    /* any type of message */
	   MPI_COMM_WORLD, /* always use this */
	   &status);       /* received message info */

gmx_rxs(pid,record(cr->pid));
  gmx_rxs(pid,record(bConverged));
  gmx_rxs(pid,record(nit));
  gmx_rxs(pid,record(k));
  gmx_rxs(pid,record(natom));
  gmx_rxs(pid,array(x,natom));
  gmx_rxs(pid,array(box,DIM));
  
  return pid;
}

void disco_slave(t_commrec *cr,FILE *log)
{
  t_correct *c;
  int       seed,nit,natom,nres,k;
  bool      bConverged,bKeep;
  matrix    box;
  rvec      boxsize;
  rvec      *x,*xref,*xcenter;
  
  c = recv_init(cr,&seed,&natom,&nres,&xref,&xcenter,&bKeep);
  
  /* Make tags etc. */
  init_corr2(c,natom);
  
  snew(x,natom);
  for(k=0; (TRUE); k++) {
    /* Generate random box*/
    rand_box(c->bBox,box,boxsize,nres,c->bCubic,&seed);
    
    /* Generate random coords */
    rand_coords(natom,x,xref,c->weight,c->bCenter,xcenter,boxsize,&seed);
    
    /* Now correct the random coords */
    bConverged = shake_coords(log,FALSE,k,natom,xref,x,&seed,box,c,&nit);
    
    if (bConverged || bKeep)
      send_coords(cr,bConverged,nit,k,natom,x,box);
  }
}

void disco_master(t_commrec *cr,FILE *fp,char *outfn,char *keepfn,t_correct *c,
		  bool bVerbose,t_atoms *atoms,
		  rvec xref[],rvec xcenter[],
		  int nstruct,int *seed,
		  bool bFit,int nfit,atom_id fit_ind[],
		  bool bPrintViol,char *violfn,rvec boxsize)
{
  FILE    *gp;
  int     *nconvdist;
  int     i,k,kk,nconv,ntry,status,kstatus,natom,nres,nit,nvtest,pid,kpid;
  double  tnit;
  rvec    *x,xcm;
  matrix  box,wrbox;
  atom_id *wr_ind;
  real    *w_rls;
  bool    bConverged;

  natom = atoms->nr;
  nres  = atoms->nres;
  
  /* Send out the word to my disciples */
  send_init(cr,c,seed,natom,nres,xref,xcenter,(keepfn != NULL));
  
  /* Make tags etc. */
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
  nconv = 0;
  ntry  = 0;
  tnit  = 0;
  for(k=0; (k<nstruct); k++) {
    pid = recv_coords(cr,&bConverged,&nit,&kpid,x,box);
    
    if (bConverged)
      nconvdist[nit]++;
    
    nvtest = quick_check(bVerbose ? fp : NULL,natom,x,box,c);
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
      write_pdbfile(gp,"Structure coloured by violation",
		    atoms,x,box,'A',TRUE);
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
  
  pr_conv_stat(fp,ntry,nconv,tnit);
  pr_conv_stat(stderr,ntry,nconv,tnit);
}

