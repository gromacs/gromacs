#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "fatal.h"
#include "dummies.h"
#include "smalloc.h"
#include "confio.h"
#include "mdrun.h"
#include "vec.h"

#define NR_END 1
#define ITMAX 200
#define EPS 1e-12
#define TOL 2.0e-6
#define CGOLD 0.3819660
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

int ncom;
real *pcom,*xicom,(*nrfunc)(real []);

real *vector(long nl,long nh)
/* allocate a real vector with subscript range v[nl..nh] */
{
  real *v;
  
  v=(real *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(real)));
  if (!v) fatal_error(0,"allocation failure in vector()");
  return v-nl+NR_END;
}

void free_vector(real *v,long nl,long nh)
/* free a real vector allocated with vector() */
{
  free((void *) (v+nl-NR_END));
}

real **f_matrix(long nrl,long nrh,long ncl,long nch)
     /* allocate a real matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  real **m;
  
  /* allocate pointers to rows */
  m=(real **) malloc((unsigned int)((nrow+NR_END)*sizeof(real*)));
  if (!m) fatal_error(0,"allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(real *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(real)));
  if (!m[nrl]) fatal_error(0,"allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

real brent(real ax,real bx,real cx,real (*f)(real),real tol,real *xmin)
{
  int iter;
  real a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  real e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm   = 0.5*(a+b);
    tol1 = tol*fabs(x)+EPS;
    tol2 = 2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  fatal_error(0,"Too many iterations in brent");
}
#undef CGOLD
#undef SHFT

real f1dim(real x)
{
  int j;
  real f,*xt;
  
  xt=vector(1,ncom);
  for (j=1;j<=ncom;j++) 
    xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_vector(xt,1,ncom);
  return f;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) { (a)=(b);(b)=(c);(c)=(d); }
static real maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))


void mnbrak(real *ax,real *bx,real *cx,real *fa,real *fb,real *fc,
	    real (*func)(real))
{
  real ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum);
    SHFT(dum,*fb,*fa,dum);
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx));
	SHFT(*fb,*fc,fu,(*func)(u));
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
  }
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT


void linmin(real p[],real xi[],int n,real *fret,real (*func)(real []))
{
  int j;
  real xx,xmin,fx,fb,fa,bx,ax;
  
  ncom   = n;
  pcom   = vector(1,n);
  xicom  = vector(1,n);
  nrfunc = func;
  for (j=1;j<=n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
  for (j=1;j<=n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free_vector(xicom,1,n);
  free_vector(pcom,1,n);
}
#undef TOL


void frprmn(real p[],int n,real ftol,int *iter,real *fret,
	    real (*func)(real[]),void (*dfunc)(real[],real[]))
{
  int j,its;
  real gg,gam,fp,dgg;
  real *g,*h,*xi;
  
  g=vector(1,n);
  h=vector(1,n);
  xi=vector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  for (j=1;j<=n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
      free_vector(xi,1,n);
      free_vector(h,1,n);
      free_vector(g,1,n);
      return;
    }
    fp=(*func)(p);
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=1;j<=n;j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      free_vector(xi,1,n);
      free_vector(h,1,n);
      free_vector(g,1,n);
      return;
    }
    gam=dgg/gg;
    for (j=1;j<=n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  fatal_error(0,"Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS

#define ALF 1.0e-4
#define TOLX 1.0e-7

void lnsrch(int n,real xold[],real fold,real g[],real p[],real x[],real *f,
	    real stpmax,int *check,real (*func)())
{
  int i;
  real a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;
  
  *check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*func)(x);
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      *check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
	tmplam = -slope/(2.0*(*f-fold-slope));
      else {
	rhs1 = *f-fold-alam*slope;
	rhs2=f2-fold2-alam2*slope;
	a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	if (a == 0.0) tmplam = -slope/(2.0*b);
	else {
	  disc=b*b-3.0*a*slope;
	  if (disc<0.0) fatal_error(0,"Roundoff problem in lnsrch.");
	  else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
	if (tmplam>0.5*alam)
	  tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    fold2=fold;
    alam=FMAX(tmplam,0.1*alam);
  }
}
#undef ALF
#undef TOLX

#define ITMAX 200
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \

void dfpmin(real p[],int n,real gtol,int *iter,real *fret,
	    real (*func)(real []),void (*dfunc)(real[],real[]))
{
  int check,i,its,j;
  real den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
  real *dg,*g,*hdg,**hessin,*pnew,*xi;

  dg=vector(1,n);
  g=vector(1,n);
  hdg=vector(1,n);
  hessin=f_matrix(1,n,1,n);
  pnew=vector(1,n);
  xi=vector(1,n);
  fp=(*func)(p);
  (*dfunc)(p,g);
  for (i=1;i<=n;i++) {
    for (j=1;j<=n;j++) hessin[i][j]=0.0;
    hessin[i][i]=1.0;
    xi[i] = -g[i];
    sum += p[i]*p[i];
  }
  stpmax=STPMX*FMAX(sqrt(sum),(real)n);
  for (its=1;its<=ITMAX;its++) {
    *iter=its;
    lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
    fp = *fret;
    for (i=1;i<=n;i++) {
      xi[i]=pnew[i]-p[i];
      p[i]=pnew[i];
    }
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) {
      FREEALL
	return;
    }
    for (i=1;i<=n;i++) dg[i]=g[i];
    (*dfunc)(p,g);
    test=0.0;
    den=FMAX(*fret,1.0);
    for (i=1;i<=n;i++) {
      temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
      if (temp > test) test=temp;
    }
    if (test < gtol) {
      FREEALL
	return;
    }
    for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
    for (i=1;i<=n;i++) {
      hdg[i]=0.0;
      for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=sumdg=sumxi=0.0;
    for (i=1;i<=n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
      sumdg += sqr(dg[i]);
      sumxi += sqr(xi[i]);
    }
    if (fac*fac > EPS*sumdg*sumxi) {
      fac=1.0/fac;
      fad=1.0/fae;
      for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
      for (i=1;i<=n;i++) {
	for (j=1;j<=n;j++) {
	  hessin[i][j] += fac*xi[i]*xi[j]
	    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	}
      }
    }
    for (i=1;i<=n;i++) {
      xi[i]=0.0;
      for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  fatal_error(0,"too many iterations in dfpmin");
  FREEALL
    }
#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL

typedef struct {
  FILE           *fp;
  int            count,start,end;
  bool           bVerbose,bDummies;
  t_topology     *top;
  t_parm         *parm;
  t_mdatoms      *mdatoms;
  t_nrnb         *nrnb;
  t_graph        *graph;
  t_commrec      *cr,*mcr;
  t_nsborder     *nsb;
  t_groups       *grps;
  t_forcerec     *fr;
  t_fcdata       *fcd;
  t_comm_dummies *dummy_comm;
  t_state        *state;
  rvec           *f,*buf;
  real           *ener,lambda;
} t_bulk;

static t_bulk bulk;

real fletch_func(real *x)
{
  tensor force_vir,shake_vir,pme_vir;
  rvec   mu_tot;
  rvec   *xp = (rvec *)(x+1);
  
  clear_mat(force_vir); 
  clear_mat(shake_vir); 
  clear_mat(pme_vir); 
  clear_rvec(mu_tot);
  
  if (bulk.bDummies)
    construct_dummies(bulk.fp,xp,bulk.nrnb,1,NULL,&bulk.top->idef,
		      bulk.graph,bulk.cr,bulk.state->box,bulk.dummy_comm);
  
  /* Calc force & energy on new positions
   * do_force always puts the charge groups in the box and shifts again
     * We do not unshift, so molecules are always whole in steep.c
     */
  do_force(bulk.fp,bulk.cr,bulk.mcr,bulk.parm,bulk.nsb,force_vir,pme_vir,
	   bulk.count,bulk.nrnb,bulk.top,bulk.grps,xp,
	   bulk.buf,bulk.f,bulk.buf,
	   bulk.mdatoms,bulk.ener,bulk.fcd,bulk.bVerbose && !(PAR(bulk.cr)), 
	   bulk.lambda,bulk.graph,(bulk.parm->ir.nstlist>0) || (bulk.count==0),
	   FALSE,bulk.fr,mu_tot,FALSE,0.0); 
  
  /* Spread the force on dummy particle to the other particles... */
  if (bulk.bDummies) 
    spread_dummy_f(bulk.fp,xp,bulk.f,bulk.nrnb,
		   &bulk.top->idef,bulk.dummy_comm,bulk.cr);
  
  /* Sum the potential energy terms from group contributions  */
  sum_epot(&(bulk.parm->ir.opts),bulk.grps,bulk.ener); 

  fprintf(stderr,"Step: %5d  Epot:  %12.5e\n",bulk.count,bulk.ener[F_EPOT]);
  bulk.count++;
  
  return bulk.ener[F_EPOT];
}

void fletch_dfunc(real *x,real *dx)
{
  int i,m;
  rvec *fptr = (rvec *)(dx+1);
  
  fletch_func(x);
  for(i=bulk.start; (i<bulk.end); i++)
    for(m=0; (m<DIM); m++)
      fptr[i][m] = -bulk.f[i][m];
}

time_t do_steep_new(FILE *log,int nfile,t_filenm fnm[], 
 		t_parm *parm,t_topology *top, 
 		t_groups *grps,t_nsborder *nsb, 
 		t_state *state,rvec grad[],rvec buf[],t_mdatoms *mdatoms, 
 		tensor ekin,real ener[],t_fcdata *fcd,t_nrnb nrnb[], 
 		bool bVerbose,bool bDummies, t_comm_dummies *dummycomm,
		t_commrec *cr,t_commrec *mcr,t_graph *graph,
		t_forcerec *fr,rvec box_size) 
{ 
  const char *SD="Fletcher!"; 
  real   constepsize,lambda,fmax; 
  real   ustep,dvdlambda,epot;
  t_vcm      *vcm;
  int        fp_ene; 
  t_mdebin   *mdebin; 
  t_nrnb     mynrnb; 
  bool   bDone,bAbort,bLR,bLJLR,bBHAM,b14; 
  time_t start_t; 
  rvec   mu_tot;
  int    nfmax,nsteps,niter=0;
  int    count=0; 
  int    i,m,start,end,gf; 
  int    Min=0; 
  int    steps_accepted=0; 
  bool   bConstrain;
  real   *xptr,*xxx;
  /* not used */
  real   terminate=0;
#define  TRY (1-Min)
  
  init_em(log,SD,parm,&lambda,&mynrnb,mu_tot,state->box,box_size,
	  fr,mdatoms,top,nsb,cr,&vcm,&start,&end);
   
  /* Print to log file  */
  start_t=print_date_and_time(log,cr->nodeid,"Started Fletcher"); 
  
  /* Set some booleans for the epot routines  */
  set_pot_bools(&(parm->ir),top,&bLR,&bLJLR,&bBHAM,&b14);
  
  /* Open the enrgy file */   
  if (MASTER(cr)) 
    fp_ene=open_enx(ftp2fn(efENX,nfile,fnm),"w"); 
  else 
    fp_ene=-1; 
  
  /* Init bin for energy stuff  */
  mdebin=init_mdebin(fp_ene,grps,&(top->atoms),&(top->idef),bLR,bLJLR,
		     bBHAM,b14,parm->ir.efep!=efepNO,parm->ir.epc,
		     parm->ir.eDispCorr,TRICLINIC(parm->ir.compress),
		     (parm->ir.etc==etcNOSEHOOVER),cr); 
  
  if (fr->ePBC != epbcNONE)
    /* Remove periodicity */
    do_pbc_first(log,state->box,box_size,fr,graph,state->x);

  /* Copy data to global arrays */
  bulk.fp       = stdlog;
  bulk.count    = 0;
  bulk.start    = 0;
  bulk.end      = mdatoms->nr;
  bulk.bVerbose = bVerbose;
  bulk.bDummies = bDummies;
  bulk.top      = top;
  bulk.parm     = parm;
  bulk.mdatoms  = mdatoms;
  bulk.nrnb     = nrnb;
  bulk.graph    = graph;
  bulk.cr       = cr; 
  bulk.mcr      = mcr;
  bulk.nsb      = nsb;
  bulk.grps     = grps;
  bulk.fr       = fr;
  bulk.fcd      = fcd;
  bulk.dummy_comm = dummycomm;
  bulk.state    = state;
  bulk.f        = grad;
  bulk.buf      = buf;
  bulk.ener     = ener;
  bulk.lambda   = lambda;
    
  /* Do the actual minimization */
  xptr = (real *)state->x[0];
  /*frprmn(xptr-1,mdatoms->nr*DIM,parm->ir.em_tol,&count,&epot,
	 fletch_func,fletch_dfunc);*/
  dfpmin(xptr-1,mdatoms->nr*DIM,parm->ir.em_tol,&count,&epot,
	 fletch_func,fletch_dfunc);
	 
  if (MASTER(cr)) 
    fprintf(stderr,"\nwriting lowest energy (%12.5e) coordinates.\n",epot); 

  write_traj(log,cr,ftp2fn(efTRN,nfile,fnm), 
	     nsb,count,(real) count, 
	     lambda,nrnb,nsb->natoms,state->x,NULL,bulk.f,state->box); 
  if (MASTER(cr)) 
    write_sto_conf(ftp2fn(efSTO,nfile,fnm),
		   *top->name, &(top->atoms),state->x,NULL,state->box);
  
  /* To print the actual number of steps we needed somewhere */
  parm->ir.nsteps=bulk.count;
  
  return start_t;
} /* That's all folks */

