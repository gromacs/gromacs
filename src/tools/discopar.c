#include <unistd.h>
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

#ifdef debug_gmx
#undef debug_gmx
#endif

#define debug_gmx() do { FILE *fp=debug ? debug : (stdlog ? stdlog : stderr);\
fprintf(fp,"PID=%d, %s  %d\n",gmx_cpu_id(),__FILE__,__LINE__); fflush(fp); } while (0)
	
static t_correct *recv_init(FILE *fp,
			    t_commrec *cr,int *seed,int *natom,int *nres,
			    rvec **xref,rvec **xcenter,bool *bKeep)
{
  t_correct *c;
  
  gmx_rxs(0,seed,sizeof(*seed));
  gmx_rxs(0,natom,sizeof(*natom));
  gmx_rxs(0,nres,sizeof(*nres));
  gmx_rxs(0,bKeep,sizeof(*bKeep));
  
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
  snew(*xref,*natom);
  gmx_rxs(0,*xref,*natom*sizeof(**xref));
  if (c->bCenter) {
    snew(*xcenter,*natom);
    gmx_rxs(0,*xcenter,*natom*sizeof(**xcenter));
  }
  
  /* Get the important input data */
  snew(c->d,c->ndist);
  gmx_rxs(0,array(c->d,c->ndist));
  
  snew(c->pepbond,c->npep);
  gmx_rxs(0,array(c->pepbond,c->npep));
  
  snew(c->imp,c->nimp+1);
  gmx_rxs(0,array(c->imp,c->nimp+1));

  snew(c->weight,*natom);
  gmx_rxs(0,array(c->weight,*natom));
  
  /* Other arrays can be deduced from these */
  fprintf(fp,"Succesfully got input data\n");
  fflush(fp);
  
  return c;
}

static void send_init(FILE *fp,
		      t_commrec *cr,t_correct *c,int *seed,int natom,int nres,
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
    gmx_txs(pid,array(xref,natom));
    
    if (c->bCenter)
      gmx_txs(pid,array(xcenter,natom));
    
    /* Send the important input data */
    gmx_txs(pid,array(c->d,c->ndist));
  
    gmx_txs(pid,array(c->pepbond,c->npep));
  
    gmx_txs(pid,array(c->imp,c->nimp+1));

    gmx_txs(pid,array(c->weight,natom));
    
    /* Other arrays can be deduced from these */
    fprintf(fp,"Succesfully sent input data to pid %d\n",pid);
    fflush(fp);
  }
}

static bool send_coords(t_commrec *cr,int nviol,int nit,int k,
			int natom,rvec x[],matrix box,bool bKeep)
{
  bool bDone;

  debug_gmx();
  gmx_tx(0,record(cr->pid));
  gmx_tx_wait(0);
  
  gmx_rxs(0,record(bDone));
  if (!bDone) {
    debug_gmx();
    gmx_txs(0,record(nviol));
    gmx_txs(0,record(nit));
    gmx_txs(0,record(k));
    if (bKeep || (nviol == 0)) {
      gmx_txs(0,record(natom));
      gmx_txs(0,array(x,natom));
      gmx_txs(0,array(box,DIM));
    }
  }
  debug_gmx();
  
  return bDone;
}

static int recv_coords(t_commrec *cr,int *nviol,int *nit,int *k,
		       rvec x[],matrix box,bool bKeep,bool bDone)
{
  int        pid,natom;
  
  /* Check whether there is something from anyone */
  debug_gmx();
  gmx_rx(MPI_ANY_SOURCE,record(pid));
  do {
    usleep(1000);
  } while (!mpiio_rx_probe(MPI_ANY_SOURCE));
  
  debug_gmx();
  if ((pid >= cr->nprocs) || (pid <= 0))
    fatal_error(0,"Reading data from pid %d",pid);
      
  gmx_txs(pid,record(bDone));
      
  if (!bDone) {
    gmx_rxs(pid,record(*nviol));
    gmx_rxs(pid,record(*nit));
    gmx_rxs(pid,record(*k));
    if (bKeep || (*nviol == 0)) {
      gmx_rxs(pid,record(natom));
      gmx_rxs(pid,array(x,natom));
      gmx_rxs(pid,array(box,DIM));
    }
  }
  return pid;
}

void disco_slave(t_commrec *cr,FILE *fp)
{
  t_correct *c;
  int       seed,nit,natom,nres,k,nviol;
  bool      bKeep,bDone;
  matrix    box;
  rvec      boxsize;
  rvec      *x,*xref,*xcenter;
  
  debug_gmx();
  c = recv_init(fp,cr,&seed,&natom,&nres,&xref,&xcenter,&bKeep);
  c->nstprint = 0;
  
  /* Make tags etc. */
  debug_gmx();
  init_corr2(c,natom);
  
  debug_gmx();
  pr_corr(fp,c);
  fprintf(fp,"Random seed   = %d\n",seed);
  
  snew(x,natom);
  bDone = FALSE;
  for(k=0; (!bDone); k++) {
    /* Generate random box*/
    debug_gmx();
    rand_box(c->bBox,box,boxsize,nres,c->bCubic,&seed);
    
    /* Generate random coords */
    debug_gmx();
    rand_coords(natom,x,xref,c->weight,c->bCenter,xcenter,boxsize,&seed);
    
    /* Now correct the random coords */
    debug_gmx();
    nviol = shake_coords(fp,FALSE,k,natom,xref,x,&seed,box,c,&nit);
    
    debug_gmx();
    bDone = send_coords(cr,nviol,nit,k,natom,x,box,bKeep);
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
  int     i,k,kk,nconv,ntry,status,kstatus,natom,nres,nit;
  int     nvtest,pid,kpid,nviol;
  double  tnit;
  rvec    *x,xcm;
  matrix  box,wrbox;
  atom_id *wr_ind;
  real    *w_rls;
  bool    bConverged;
  
  natom = atoms->nr;
  nres  = atoms->nres;
  
  /* Send out the word to my disciples */
  debug_gmx();
  send_init(fp,cr,c,seed,natom,nres,xref,xcenter,(keepfn != NULL));
  
  /* Make tags etc. */
  debug_gmx();
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
  debug_gmx();
  for(k=0; (k<nstruct); ) {
    pid        = recv_coords(cr,&nviol,&nit,&kpid,x,box,keepfn,FALSE);
    bConverged = (nviol == 0);
    tnit      += nit;
    ntry++;
    
    if (bConverged) 
      nconvdist[nit]++;
    
    /*nvtest = quick_check(bVerbose ? fp : NULL,natom,x,box,c);*/
    if ((k % 20) == 0)
      fprintf(stderr,"%8s%6s%6s%8s%6s\n",
	      "Struct","CPU","Str","nviol","niter");
    fprintf(stderr,"%8d%6d%6d%8d%6d\n",k,pid,kpid,nviol,nit);
    
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
    if (bPrintViol && 0) {
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
  for(k=1; (k<cr->nprocs); k++)
    pid = recv_coords(cr,&nviol,&nit,&kpid,x,box,keepfn,TRUE);

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

