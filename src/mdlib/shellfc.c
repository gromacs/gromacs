#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "vec.h"
#include "txtdump.h"
#include "mdrun.h"
#include "partdec.h"
#include "mdatoms.h"
#include "vsite.h"
#include "network.h"
#include "names.h"
#include "constr.h"
#include "domdec.h"
#include "partdec.h"
#include "physics.h"
#include "copyrite.h"
#include "shellfc.h"


typedef struct {
  int     nnucl;
  atom_id shell;	        /* The shell id				*/
  atom_id nucl1,nucl2,nucl3;	/* The nuclei connected to the shell	*/
  /* bool    bInterCG; */       /* Coupled to nuclei outside cg?        */
  real    k;		        /* force constant		        */
  real    k_1;		        /* 1 over force constant		*/
  rvec    xold;
  rvec    fold;
  rvec    step;
} t_shell;

typedef struct gmx_shellfc {
  int     nshell_gl;       /* The number of shells in the system       */
  t_shell *shell_gl;       /* All the shells (for DD only)             */
  int     *shell_index_gl; /* Global shell index (for DD only)         */
  int     nshell;          /* The number of local shells               */
  t_shell *shell;          /* The local shells                         */
  int     shell_nalloc;    /* The allocation size of shell             */
  bool    bPredict;        /* Predict shell positions                  */
  bool    bForceInit;      /* Force initialization of shell positions  */
  int     nflexcon;        /* The number of flexible constraints       */
  rvec    *x[2];           /* Array for iterative minimization         */
  rvec    *f[2];           /* Array for iterative minimization         */
  int     x_nalloc;        /* The allocation size of x and f           */
  rvec    *acc_dir;        /* Acceleration direction for flexcon       */
  rvec    *x_old;          /* Old coordinates for flexcon              */
  int     flex_nalloc;     /* The allocation size of acc_dir and x_old */
} t_gmx_shellfc;

	
static void pr_shell(FILE *fplog,int ns,t_shell s[])
{
  int i;
  
  fprintf(fplog,"SHELL DATA\n");
  fprintf(fplog,"%5s  %8s  %5s  %5s  %5s\n",
	  "Shell","Force k","Nucl1","Nucl2","Nucl3");
  for(i=0; (i<ns); i++) {
    fprintf(fplog,"%5d  %8.3f  %5d",s[i].shell,1.0/s[i].k_1,s[i].nucl1);
    if (s[i].nnucl == 2)
      fprintf(fplog,"  %5d\n",s[i].nucl2);
    else if (s[i].nnucl == 3)
      fprintf(fplog,"  %5d  %5d\n",s[i].nucl2,s[i].nucl3);
    else
      fprintf(fplog,"\n");
  }
}

static void predict_shells(FILE *fplog,rvec x[],rvec v[],real dt,
			   int ns,t_shell s[],
			   real mass[],t_atom *atom,bool bInit)
{
  int  i,m,s1,n1,n2,n3;
  real dt_1,dt_2,dt_3,fudge,tm,m1,m2,m3;
  rvec *ptr;
  
  /* We introduce a fudge factor for performance reasons: with this choice
   * the initial force on the shells is about a factor of two lower than 
   * without
   */
  fudge = 1.0;
    
  if (bInit) {
    if (fplog)
      fprintf(fplog,"RELAX: Using prediction for initial shell placement\n");
    ptr  = x;
    dt_1 = 1;
  }
  else {
    ptr  = v;
    dt_1 = fudge*dt;
  }
    
  for(i=0; (i<ns); i++) {
    s1 = s[i].shell;
    if (bInit)
      clear_rvec(x[s1]);
    switch (s[i].nnucl) {
    case 1:
      n1 = s[i].nucl1;
      for(m=0; (m<DIM); m++)
	x[s1][m]+=ptr[n1][m]*dt_1;
      break;
    case 2:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      if (mass) {
	m1 = mass[n1];
	m2 = mass[n2];
      } else {
	/* Not the correct masses with FE, but it is just a prediction... */
	m1 = atom[n1].m;
	m2 = atom[n2].m;
      }
      tm = dt_1/(m1+m2);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m])*tm;
      break;
    case 3:
      n1 = s[i].nucl1;
      n2 = s[i].nucl2;
      n3 = s[i].nucl3;
      if (mass) {
	m1 = mass[n1];
	m2 = mass[n2];
	m3 = mass[n3];
      } else {
	/* Not the correct masses with FE, but it is just a prediction... */
	m1 = atom[n1].m;
	m2 = atom[n2].m;
	m3 = atom[n3].m;
      }
      tm = dt_1/(m1+m2+m3);
      for(m=0; (m<DIM); m++)
	x[s1][m]+=(m1*ptr[n1][m]+m2*ptr[n2][m]+m3*ptr[n3][m])*tm;
      break;
    default:
      gmx_fatal(FARGS,"Shell %d has %d nuclei!",i,s[i].nnucl);
    }
  }
}

gmx_shellfc_t init_shell_flexcon(FILE *fplog,t_commrec *cr,
				 t_topology *top,int nflexcon,
				 bool bContinuation,rvec *x)
{
  struct gmx_shellfc *shfc;
  t_shell     *shell;
  int         *shell_index,*at2cg;
  t_atom      *atom;
  t_idef      *idef;
  int         n[eptNR],ns,nsi,nshell,nshell_tot;
  int         cg0,cg1,*cgindex,start,end,cg,i,j,type,ftype,nra;
  real        qS,alpha;
  int         aS,aN=0; /* Shell and nucleus */
  int         bondtypes[] = { F_BONDS, F_HARMONIC, F_CUBICBONDS, F_POLARIZATION, F_WATER_POL };
#define NBT asize(bondtypes)
  bool        bInterCG;
  t_iatom     *ia;

  if (PARTDECOMP(cr)) {
    pd_cg_range(cr,&cg0,&cg1);
  } else {
    cg0 = 0;
    cg1 = top->blocks[ebCGS].nr;
  }

  atom = top->atoms.atom;
  idef = &top->idef;
  cgindex = top->blocks[ebCGS].index;
  
  start = cgindex[cg0];
  end   = cgindex[cg1];

  /* Count number of shells, and find their indices */
  for(i=0; (i<eptNR); i++)
    n[i]=0;
  snew(shell_index,end-start);
  snew(at2cg,end-start);
  nsi = 0;
  for(cg=cg0; (cg<cg1); cg++) {
    for(i=cgindex[cg]; i<cgindex[cg+1]; i++) {
      n[atom[i].ptype]++;
      if (atom[i].ptype == eptShell)
	shell_index[i-start] = nsi++;
      at2cg[i-start] = cg;
    }
  }
  if (nsi != n[eptShell])
    gmx_fatal(FARGS,"Your number of shells %d is not equal to the number of shells %d",
		nsi,n[eptShell]);

  if (fplog) {
    /* Print the number of each particle type */  
    for(i=0; (i<eptNR); i++)
      if (n[i]!=0)
	fprintf(fplog,"There are: %d %ss\n",n[i],ptype_str[i]);
  }
  
  ns      = n[eptShell];
  nshell  = ns;

  nshell_tot = ns;
  if (PARTDECOMP(cr))
    gmx_sumi(1,&nshell_tot,cr);

  if (nshell_tot == 0 && nflexcon == 0)
    return NULL;

  snew(shfc,1);
  shfc->nflexcon = nflexcon;

  if (nshell_tot == 0)
    return shfc;

  snew(shell,ns);
  shfc->nflexcon = nflexcon;
  
  /* Initiate the shell structures */    
  for(i=0; (i<ns); i++) {
    shell[i].shell=NO_ATID;
    shell[i].nnucl=0;
    shell[i].nucl1=NO_ATID;
    shell[i].nucl2=NO_ATID;
    shell[i].nucl3=NO_ATID;
    /* shell[i].bInterCG=FALSE; */
    shell[i].k_1=0;
    shell[i].k=0;
  }

  /* Now fill the structures */
  bInterCG = FALSE;
  ns=0;
  for(j=0; (j<NBT); j++) {
    ia=idef->il[bondtypes[j]].iatoms;
    for(i=0; (i<idef->il[bondtypes[j]].nr); ) {
      type  = ia[0];
      ftype = idef->functype[type];
      nra   = interaction_function[ftype].nratoms;
      
      /* Check whether we have a bond with a shell */
      aS = NO_ATID;
      
      switch (bondtypes[j]) {
      case F_BONDS:
      case F_HARMONIC:
      case F_CUBICBONDS:
      case F_POLARIZATION:
	if (atom[ia[1]].ptype == eptShell) {
	  aS = ia[1];
	  aN = ia[2];
	}
	else if (atom[ia[2]].ptype == eptShell) {
	  aS = ia[2];
	  aN = ia[1];
	}
	break;
      case F_WATER_POL:
	aN    = ia[4];  /* Dummy */
	aS    = ia[5];  /* Shell */
	break;
      default:
	gmx_fatal(FARGS,"Death Horror: %s, %d",__FILE__,__LINE__);
      }
      
      if (aS != NO_ATID) {	  
	qS = atom[aS].q;
	
	/* Check whether one of the particles is a shell... */
	nsi = shell_index[aS-start];
	if ((nsi < 0) || (nsi >= nshell))
	  gmx_fatal(FARGS,"nsi is %d should be within 0 - %d. aS = %d",
		    nsi,nshell,aS);
	if (shell[nsi].shell == NO_ATID) {
	  shell[nsi].shell = aS;
	  ns ++;
	}
	else if (shell[nsi].shell != aS)
	  gmx_fatal(FARGS,"Weird stuff in %s, %d",__FILE__,__LINE__);
	
	if      (shell[nsi].nucl1 == NO_ATID) {
	  shell[nsi].nucl1 = aN;
	} else if (shell[nsi].nucl2 == NO_ATID) {
	  shell[nsi].nucl2 = aN;
	} else if (shell[nsi].nucl3 == NO_ATID) {
	  shell[nsi].nucl3 = aN;
	} else {
	  if (fplog)
	    pr_shell(fplog,ns,shell);
	  gmx_fatal(FARGS,"Can not handle more than three bonds per shell\n");
	}
	if (at2cg[aS] != at2cg[aN]) {
	  /* shell[nsi].bInterCG = TRUE; */
	  bInterCG = TRUE;
	}

	switch (bondtypes[j]) {
	case F_BONDS:
	case F_HARMONIC:
	  shell[nsi].k    += idef->iparams[type].harmonic.krA;
	  break;
	case F_CUBICBONDS:
	  shell[nsi].k    += idef->iparams[type].cubic.kb;
	  break;
	case F_POLARIZATION:
	  if (qS != atom[aS].qB)
	    gmx_fatal(FARGS,"polarize can not be used with qA != qB");
	  shell[nsi].k    += sqr(qS)*ONE_4PI_EPS0/
	    idef->iparams[type].polarize.alpha;
	  break;
	case F_WATER_POL:
	  if (qS != atom[aS].qB)
	    gmx_fatal(FARGS,"water_pol can not be used with qA != qB");
	  alpha          = (idef->iparams[type].wpol.al_x+
			    idef->iparams[type].wpol.al_y+
			    idef->iparams[type].wpol.al_z)/3.0;
	  shell[nsi].k  += sqr(qS)*ONE_4PI_EPS0/alpha;
	  break;
	default:
	  gmx_fatal(FARGS,"Death Horror: %s, %d",__FILE__,__LINE__);
	}
	shell[nsi].nnucl++;
      }
      ia += nra+1;
      i  += nra+1;
    }
  }

  sfree(at2cg);
  
  /* Verify whether it's all correct */
  if (ns != nshell)
    gmx_fatal(FARGS,"Something weird with shells. They may not be bonded to something");
  
  for(i=0; (i<ns); i++)
    shell[i].k_1 = 1.0/shell[i].k;
  
  if (debug)
    pr_shell(debug,ns,shell);

  
  shfc->nshell_gl      = ns;
  shfc->shell_gl       = shell;
  if (DOMAINDECOMP(cr)) {
    shfc->shell_index_gl = shell_index;
  } else {
    shfc->nshell         = ns;
    shfc->shell          = shell;
    sfree(shell_index);
  }

  shfc->bPredict   = (getenv("GMX_NOPREDICT") == NULL);
  shfc->bForceInit = FALSE;
  if (!shfc->bPredict) {
    if (fplog)
      fprintf(fplog,"\nWill never predict shell positions\n");
  } else {
    shfc->bForceInit = (getenv("GMX_FORCEINIT") != NULL);
    if (shfc->bForceInit && fplog)
      fprintf(fplog,"\nWill always initiate shell positions\n");
  }

  if (shfc->bPredict) {
    if (!bContinuation && (!DOMAINDECOMP(cr) || MASTER(cr))) {
      predict_shells(fplog,x,NULL,0,shfc->nshell_gl,shfc->shell_gl,
		     NULL,top->atoms.atom,TRUE);
    }

    if (bInterCG) {
      if (fplog)
	fprintf(fplog,"\nNOTE: there all shells that are connected to particles outside thier own charge group, will not predict shells positions during the run\n\n");
      shfc->bPredict = FALSE;
    }
  }

  return shfc;
}

void make_local_shells(gmx_domdec_t *dd,t_mdatoms *md,
		      struct gmx_shellfc *shfc)
{
  t_shell *shell;
  int nshell,i;

  nshell = 0;
  shell  = shfc->shell; 
  for(i=0; i<dd->nat_home; i++) {
    if (md->ptype[i] == eptShell) {
      if (nshell+1 > shfc->shell_nalloc) {
	shfc->shell_nalloc = over_alloc_dd(nshell+1);
	srenew(shell,shfc->shell_nalloc);
      }
      shell[nshell] = shfc->shell_gl[shfc->shell_index_gl[dd->gatindex[i]]];
      shell[nshell].shell = i;
      /* More is required here */
      if (shell[nshell].k < 1)
	printf("node %d shell[%d/%d].k = %f\n",
	       dd->sim_nodeid,nshell,dd->gatindex[i],shell[nshell].k);
      nshell++;
    }
  }

  shfc->nshell = nshell;
  shfc->shell  = shell;
}

static void do_1pos(rvec xnew,rvec xold,rvec f,real step)
{
  real xo,yo,zo;
  real dx,dy,dz;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];

  dx=f[XX]*step;
  dy=f[YY]*step;
  dz=f[ZZ]*step;

  xnew[XX]=xo+dx;
  xnew[YY]=yo+dy;
  xnew[ZZ]=zo+dz;
}

static void do_1pos3(rvec xnew,rvec xold,rvec f,rvec step)
{
  real xo,yo,zo;
  real dx,dy,dz;
  
  xo=xold[XX];
  yo=xold[YY];
  zo=xold[ZZ];

  dx=f[XX]*step[XX];
  dy=f[YY]*step[YY];
  dz=f[ZZ]*step[ZZ];

  xnew[XX]=xo+dx;
  xnew[YY]=yo+dy;
  xnew[ZZ]=zo+dz;
}

static void directional_sd(FILE *log,rvec xold[],rvec xnew[],rvec acc_dir[],
			   int start,int homenr,real step)
{
  int  i;

  for(i=start; i<homenr; i++)
    do_1pos(xnew[i],xold[i],acc_dir[i],step);
}

static void shell_pos_sd(FILE *log,rvec xcur[],rvec xnew[],rvec f[],
			 int ns,t_shell s[],int count)
{
  int  i,shell,d;
  real dx,df,k_est;
#ifdef PRINT_STEP  
  real step_min,step_max;

  step_min = 1e30;
  step_max = 0;
#endif
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    if (count == 1) {
      for(d=0; d<DIM; d++) {
	s[i].step[d] = s[i].k_1;
#ifdef PRINT_STEP
	step_min = min(step_min,s[i].step[d]);
	step_max = max(step_max,s[i].step[d]);
#endif
      }
    } else {
      for(d=0; d<DIM; d++) {
	dx = xcur[shell][d] - s[i].xold[d];
	df =    f[shell][d] - s[i].fold[d];
	if (dx != 0 && df != 0) {
	  k_est = -dx/df;
	  if (k_est >= 2*s[i].step[d]) {
	    s[i].step[d] *= 1.2;
	  } else if (k_est <= 0) {
	    s[i].step[d] *= 0.8;
	  } else {
	    s[i].step[d] = 0.8*s[i].step[d] + 0.2*k_est;
	  }
	} else if (dx != 0) {
	  s[i].step[d] *= 1.2;
	}
#ifdef PRINT_STEP
	step_min = min(step_min,s[i].step[d]);
	step_max = max(step_max,s[i].step[d]);
#endif
      }
    }
    copy_rvec(xcur[shell],s[i].xold);
    copy_rvec(f[shell],   s[i].fold);

    do_1pos3(xnew[shell],xcur[shell],f[shell],s[i].step);

    if (gmx_debug_at) {
      fprintf(debug,"shell[%d] = %d\n",i,shell);
      pr_rvec(debug,0,"fshell",f[shell],DIM,TRUE);
      pr_rvec(debug,0,"xold",xcur[shell],DIM,TRUE);
      pr_rvec(debug,0,"step",s[i].step,DIM,TRUE);
      pr_rvec(debug,0,"xnew",xnew[shell],DIM,TRUE);
    }
  }
#ifdef PRINT_STEP
  printf("step %.3e %.3e\n",step_min,step_max);
#endif
}

static void decrease_step_size(int nshell,t_shell s[])
{
  int i;
  
  for(i=0; i<nshell; i++)
    svmul(0.8,s[i].step,s[i].step);
}

static void print_epot(FILE *fp,int mdstep,int count,real epot,real df,
		       int ndir,real sf_dir)
{
  fprintf(fp,"MDStep=%5d/%2d EPot: %12.8e, rmsF: %6.2e",
	  mdstep,count,epot,df);
  if (ndir)
    fprintf(fp,", dir. rmsF: %6.2e\n",sqrt(sf_dir/ndir));
  else
    fprintf(fp,"\n");
}


static real rms_force(t_commrec *cr,rvec f[],int ns,t_shell s[],
		      int ndir,real *sf_dir,real *Epot)
{
  int  i,shell,ntot;
  double buf[4];

  buf[0] = *sf_dir;
  for(i=0; i<ns; i++) {
    shell = s[i].shell;
    buf[0]  += norm2(f[shell]);
  }
  ntot = ns;

  if (PAR(cr)) {
    buf[1] = ntot;
    buf[2] = *sf_dir;
    buf[3] = *Epot;
    gmx_sumd(4,buf,cr);
    ntot = (int)(buf[1] + 0.5);
    *sf_dir = buf[2];
    *Epot   = buf[3];
  }
  ntot += ndir;

  return (ntot ? sqrt(buf[0]/ntot) : 0);
}

static void check_pbc(FILE *fp,rvec x[],int shell)
{
  int m,now;
  
  now = shell-4;
  for(m=0; (m<DIM); m++)
    if (fabs(x[shell][m]-x[now][m]) > 0.3) {
      pr_rvecs(fp,0,"SHELL-X",x+now,5);
      break;
    }
}

static void dump_shells(FILE *fp,rvec x[],rvec f[],real ftol,int ns,t_shell s[])
{
  int  i,shell;
  real ft2,ff2;
  
  ft2 = sqr(ftol);
  
  for(i=0; (i<ns); i++) {
    shell = s[i].shell;
    ff2   = iprod(f[shell],f[shell]);
    if (ff2 > ft2)
      fprintf(fp,"SHELL %5d, force %10.5f  %10.5f  %10.5f, |f| %10.5f\n",
	      shell,f[shell][XX],f[shell][YY],f[shell][ZZ],sqrt(ff2));
    check_pbc(fp,x,shell);
  }
}

static void init_adir(FILE *log,
		      gmx_constr_t constr,t_topology *top,t_inputrec *ir,
		      gmx_domdec_t *dd,int dd_ac1,
		      int step,t_mdatoms *md,int start,int end,
		      rvec *x_old,rvec *x_init,rvec *x,
		      rvec *f,rvec *acc_dir,matrix box,
		      real lambda,real *dvdlambda,t_nrnb *nrnb)
{
  static rvec *xnold=NULL,*xnew=NULL;
  static int x_nalloc=0;
  double w_dt;
  int    gf,ga,gt;
  real   dt,scale;
  int    n,d; 
  unsigned short *ptype;
  rvec   p,dx;
  
  if (dd)
    n = dd_ac1;
  else
    n = end - start;
  if (n > x_nalloc) {
    x_nalloc = over_alloc_dd(n);
    srenew(xnold,x_nalloc);
    srenew(xnew,x_nalloc);
  }
    
  ptype = md->ptype;

  dt = ir->delta_t;

  /* Does NOT work with freeze or acceleration groups (yet) */
  for (n=start; n<end; n++) {  
    w_dt = md->invmass[n]*dt;
    
    for (d=0; d<DIM; d++) {
      if ((ptype[n] != eptVSite) && (ptype[n] != eptShell)) {
	xnold[n-start][d] = x[n][d] - (x_init[n][d] - x_old[n][d]);
	xnew[n-start][d] = 2*x[n][d] - x_old[n][d] + f[n][d]*w_dt*dt;
      } else {
	xnold[n-start][d] = x[n][d];
	xnew[n-start][d] = x[n][d];
      }
    }
  }
  constrain(log,FALSE,FALSE,constr,top,ir,dd,step,md,
	    x,xnold-start,NULL,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,TRUE);
  constrain(log,FALSE,FALSE,constr,top,ir,dd,step,md,
	    x,xnew-start,NULL,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,TRUE);

  /* Set xnew to minus the acceleration */
  for (n=start; n<end; n++) {
    for(d=0; d<DIM; d++)
      xnew[n-start][d] =
	-(2*x[n][d]-xnold[n-start][d]-xnew[n-start][d])/sqr(dt)
	- f[n][d]*md->invmass[n];
    clear_rvec(acc_dir[n]);
  }

  /* Project the acceleration on the old bond directions */
  constrain(log,FALSE,FALSE,constr,top,ir,dd,step,md,
	    x_old,xnew-start,acc_dir,box,
	    lambda,dvdlambda,dt,NULL,NULL,nrnb,FALSE); 
}

int relax_shell_flexcon(FILE *fplog,t_commrec *cr,bool bVerbose,
			int mdstep,t_inputrec *inputrec,
			bool bDoNS,bool bStopCM,
			t_topology *top,gmx_constr_t constr,
			real ener[],t_fcdata *fcd,
			t_state *state,rvec f[],
			rvec buf[],t_mdatoms *md,
			t_nrnb *nrnb,gmx_wallcycle_t wcycle,
			t_graph *graph,t_groups *grps,
			struct gmx_shellfc *shfc,
			t_forcerec *fr,
			real t,rvec mu_tot,
			int natoms,bool *bConverged,
			gmx_vsite_t *vsite,
			FILE *fp_field)
{
  int    nshell;
  t_shell *shell;
  rvec   *pos[2],*force[2],*acc_dir=NULL,*x_old=NULL;
  real   Epot[2],df[2],Estore[F_NRE];
  tensor vir_el_recip[2];
  rvec   dx;
  real   sf_dir,invdt;
  real   ftol,xiH,xiS,dum=0;
  char   cbuf[56];
  bool   bCont,bInit;
  int    nat,dd_ac0,dd_ac1=0,i;
  int    start=md->start,homenr=md->homenr,end=start+homenr,cg0,cg1;
  int    nflexcon,g,number_steps,d,Min=0,count=0;
#define  Try (1-Min)             /* At start Try = 1 */

  bCont        = (mdstep == inputrec->init_step) && inputrec->bContinuation;
  bInit        = (mdstep == inputrec->init_step) || shfc->bForceInit;
  ftol         = inputrec->em_tol;
  number_steps = inputrec->niter;
  nshell       = shfc->nshell;
  shell        = shfc->shell;
  nflexcon     = shfc->nflexcon;

  if (DOMAINDECOMP(cr)) {
    nat = dd_natoms_vsite(cr->dd);
    if (nflexcon > 0) {
      dd_get_constraint_range(cr->dd,&dd_ac0,&dd_ac1);
      nat = max(nat,dd_ac1);
    }
  } else {
    nat = state->natoms;
  }

  if (nat > shfc->x_nalloc) {
    /* Allocate local arrays */
    shfc->x_nalloc = over_alloc_dd(nat);
    for(i=0; (i<2); i++) {
      srenew(shfc->x[i],shfc->x_nalloc);
      srenew(shfc->f[i],shfc->x_nalloc);
    }
  }
  for(i=0; (i<2); i++) {
    pos[i]   = shfc->x[i];
    force[i] = shfc->f[i];
  }
     
  /* With particle decomposition this code only works
   * when all particles involved with each shell are in the same cg.
   */

  if (bDoNS && !DOMAINDECOMP(cr)) {
    /* This is the only time where the coordinates are used
     * before do_force is called, which normally puts all
     * charge groups in the box.
     */
    if (PARTDECOMP(cr)) {
      pd_cg_range(cr,&cg0,&cg1);
    } else {
      cg0 = 0;
      cg1 = top->blocks[ebCGS].nr;
    }
    put_charge_groups_in_box(fplog,cg0,cg1,fr->ePBC,state->box,
			     &(top->blocks[ebCGS]),state->x,fr->cg_cm);
    if (graph)
      mk_mshift(fplog,graph,fr->ePBC,state->box,state->x);
  }

  /* After this all coordinate arrays will contain whole molecules */
  if (graph)
    shift_self(graph,state->box,state->x);

  if (nflexcon) {
    if (nat > shfc->flex_nalloc) {
      shfc->flex_nalloc = over_alloc_dd(nat);
      srenew(shfc->acc_dir,shfc->flex_nalloc);
      srenew(shfc->x_old,shfc->flex_nalloc);
    }
    acc_dir = shfc->acc_dir;
    x_old   = shfc->x_old;
    for(i=0; i<homenr; i++) {
      for(d=0; d<DIM; d++)
        shfc->x_old[i][d] =
	  state->x[start+i][d] - state->v[start+i][d]*inputrec->delta_t;
    }
  }
  
  /* Do a prediction of the shell positions */
  if (shfc->bPredict && !bCont)
    predict_shells(fplog,state->x,state->v,inputrec->delta_t,nshell,shell,
		   md->massT,NULL,bInit);

  /* do_force expected the charge groups to be in the box */
  if (graph)
    unshift_self(graph,state->box,state->x);

  /* Calculate the forces first time around */
  if (gmx_debug_at) {
    pr_rvecs(debug,0,"x b4 do_force",state->x + start,homenr);
  }
  do_force(fplog,cr,inputrec,mdstep,nrnb,wcycle,top,grps,
	   state->box,state->x,force[Min],buf,md,ener,fcd,
	   state->lambda,graph,
	   TRUE,bDoNS,FALSE,TRUE,fr,mu_tot,FALSE,t,fp_field,NULL);
  sum_lrforces(force[Min],fr,start,homenr);
  copy_mat(fr->vir_el_recip,vir_el_recip[Min]);

  sf_dir = 0;
  if (nflexcon) {
    init_adir(fplog,constr,top,inputrec,cr->dd,dd_ac1,mdstep,md,start,end,
	      shfc->x_old-start,state->x,state->x,force[Min],
	      shfc->acc_dir-start,state->box,state->lambda,&dum,nrnb);

    for(i=start; i<end; i++)
      sf_dir += md->massT[i]*norm2(shfc->acc_dir[i-start]);
  }

  /* Sum the potential energy terms from group contributions */
  sum_epot(&(inputrec->opts),grps,ener);
  Epot[Min]=ener[F_EPOT];

  df[Min]=rms_force(cr,shfc->f[Min],nshell,shell,nflexcon,&sf_dir,&Epot[Min]);
  df[Try]=0;
  if (debug) {
    fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);
  }

  if (gmx_debug_at) {
    pr_rvecs(debug,0,"force0",force[Min],md->nr);
  }

  if (nshell+nflexcon > 0) {
    /* Copy x to pos[Min] & pos[Try]: during minimization only the
     * shell positions are updated, therefore the other particles must
     * be set here.
     */
    memcpy(pos[Min],state->x,nat*sizeof(state->x[0]));
    memcpy(pos[Try],state->x,nat*sizeof(state->x[0]));
  }
  
  if (bVerbose && MASTER(cr))
    print_epot(stdout,mdstep,0,Epot[Min],df[Min],nflexcon,sf_dir);

  if (debug) {
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EKIN].longname, ener[F_EKIN]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_EPOT].longname, ener[F_EPOT]);
    fprintf(debug,"%17s: %14.10e\n",
	    interaction_function[F_ETOT].longname, ener[F_ETOT]);
    fprintf(debug,"SHELLSTEP %d\n",mdstep);
  }
  
  /* First check whether we should do shells, or whether the force is 
   * low enough even without minimization.
   */
  *bConverged = (df[Min] < ftol);
  
  for(count=1; (!(*bConverged) && (count < number_steps)); count++) {
    if (vsite)
      construct_vsites(fplog,vsite,pos[Min],nrnb,inputrec->delta_t,state->v,
		       &top->idef,fr->ePBC,fr->bMolPBC,graph,cr,state->box);
     
    if (nflexcon) {
      init_adir(fplog,constr,top,inputrec,cr->dd,dd_ac1,mdstep,md,start,end,
		x_old-start,state->x,pos[Min],force[Min],acc_dir-start,
		state->box,state->lambda,&dum,nrnb);
      
      directional_sd(fplog,pos[Min],pos[Try],acc_dir-start,start,end,
		     fr->fc_stepsize);
    }
    
    /* New positions, Steepest descent */
    shell_pos_sd(fplog,pos[Min],pos[Try],force[Min],nshell,shell,count); 

    /* do_force expected the charge groups to be in the box */
    if (graph)
      unshift_self(graph,state->box,pos[Try]);

    if (gmx_debug_at) {
      pr_rvecs(debug,0,"RELAX: pos[Min]  ",pos[Min] + start,homenr);
      pr_rvecs(debug,0,"RELAX: pos[Try]  ",pos[Try] + start,homenr);
    }
    /* Try the new positions */
    do_force(fplog,cr,inputrec,1,nrnb,wcycle,
	     top,grps,state->box,pos[Try],force[Try],buf,md,ener,fcd,
	     state->lambda,graph,
	     TRUE,FALSE,FALSE,TRUE,fr,mu_tot,FALSE,t,fp_field,NULL);
    if (vsite) 
      spread_vsite_f(fplog,vsite,pos[Try],force[Try],fr->fshift,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,state->box,cr);
      
    /* Calculation of the virial must be done after vsites!    */
    /* Question: Is it correct to do the PME forces after this? */
    /*    calc_virial(fplog,START(nsb),HOMENR(nsb),pos[Try],force[Try],
		my_vir[Try],pme_vir[Try],graph,state->box,nrnb,fr,FALSE);
    */	  
    /* Spread the LR force on virtual site to the other particles... 
     * This is parallellized. MPI communication is performed
     * if the constructing atoms aren't local.
     */
    if (vsite && fr->bEwald) 
      spread_vsite_f(fplog,vsite,pos[Try],fr->f_el_recip,NULL,nrnb,
		     &top->idef,fr->ePBC,fr->bMolPBC,graph,state->box,cr);
    
    sum_lrforces(force[Try],fr,start,homenr);
    copy_mat(fr->vir_el_recip,vir_el_recip[Try]);
    
    if (gmx_debug_at) {
      pr_rvecs(debug,0,"RELAX: force[Min]",force[Min] + start,homenr);
      pr_rvecs(debug,0,"RELAX: force[Try]",force[Try] + start,homenr);
    }
    sf_dir = 0;
    if (nflexcon) {
      init_adir(fplog,constr,top,inputrec,cr->dd,dd_ac1,mdstep,md,start,end,
		x_old-start,state->x,pos[Try],force[Try],acc_dir-start,
		state->box,state->lambda,&dum,nrnb);

      for(i=start; i<end; i++)
	sf_dir += md->massT[i]*norm2(acc_dir[i-start]);
    }

    /* Sum the potential energy terms from group contributions */
    sum_epot(&(inputrec->opts),grps,ener);
    Epot[Try]=ener[F_EPOT]; 
    
    df[Try]=rms_force(cr,force[Try],nshell,shell,nflexcon,&sf_dir,&Epot[Try]);

    if (debug)
      fprintf(debug,"df = %g  %g\n",df[Min],df[Try]);

    if (debug) {
      if (gmx_debug_at)
	pr_rvecs(debug,0,"F na do_force",force[Try] + start,homenr);
      if (gmx_debug_at) {
	fprintf(debug,"SHELL ITER %d\n",count);
	dump_shells(debug,pos[Try],force[Try],ftol,nshell,shell);
      }
    }

    if (bVerbose && MASTER(cr))
      print_epot(stdout,mdstep,count,Epot[Try],df[Try],nflexcon,sf_dir);
      
    *bConverged = (df[Try] < ftol);
    
    if ((df[Try] < df[Min])) {
      if (debug)
	fprintf(debug,"Swapping Min and Try\n");
      if (nflexcon) {
	/* Correct the velocities for the flexible constraints */
	invdt = 1/inputrec->delta_t;
	for(i=start; i<end; i++)
	  for(d=0; d<DIM; d++)
	    state->v[i][d] += (pos[Try][i][d] - pos[Min][i][d])*invdt;
      }
      Min  = Try;
    } else {
      decrease_step_size(nshell,shell);
    }
  }
  if (MASTER(cr) && !(*bConverged)) {
    fprintf(fplog,
	    "step %d: EM did not converge in %d iterations, RMS force %.3f\n",
	    mdstep,number_steps,df[Min]);
    fprintf(stderr,
	    "step %d: EM did not converge in %d iterations, RMS force %.3f\n",
	    mdstep,number_steps,df[Min]);
  }

  /* Parallelise this one! */
  if (EEL_FULL(fr->eeltype)) {
    for(i=start; (i<end); i++)
      rvec_dec(force[Min][i],fr->f_el_recip[i]);
  }
  memcpy(f,force[Min],nat*sizeof(f[0]));

  /* CHECK VIRIAL */
  copy_mat(vir_el_recip[Min],fr->vir_el_recip);
  
  memcpy(state->x,pos[Min],nat*sizeof(state->x[0]));

  return count; 
}

