#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <time.h>
#include "filenm.h"
#include "names.h"
#include "typedefs.h"
#include "smalloc.h"
#include "gmxfio.h"
#include "xdrf.h"
#include "checkpoint.h"

#define CPT_MAGIC 171817

static const cpt_version = 1;

enum { ecpdtINT, ecpdtFLOAT, ecpdtDOUBLE, ecpdtNR };

const char *ecpdt_names[ecpdtNR] = { "int", "float", "double" };

const char *est_names[estNR]=
{
  "FE-lambda",
  "box", "box-rel", "box-v", "pcoupl-mu",
  "nosehoover-xi", "nosehoover-xi-integral",
  "x", "v", "SDx", "CGp", "LD-rng", "LD-rng-i",
};

static void do_cpt_string(XDR *xd,bool bRead,char *desc,char **s,FILE *list)
{
#define CPTSTRLEN 1024
  bool_t res=0;
  int dt=ecpdtINT;

  if (bRead)
    snew(*s,CPTSTRLEN);
  res = xdr_string(xd,s,CPTSTRLEN);
  if (list) {
    fprintf(list,"%s = %s\n",desc,*s);
    sfree(*s);
  }
}

static void do_cpt_int(XDR *xd,char *desc,int *i,FILE *list)
{
  bool_t res=0;

  res = xdr_int(xd,i);
  if (list)
    fprintf(list,"%s = %d\n",desc,*i);
}

static void do_cpt_double(XDR *xd,char *desc,double *f,FILE *list)
{
  bool_t res=0;

  res = xdr_double(xd,f);
  if (list)
    fprintf(list,"%s = %f\n",desc,*f);
}

static void do_cpte_reals(XDR *xd,int ecpt,int sflags,int n,real **v,
			  FILE *list)
{
  bool_t res=0;
#ifndef GMX_DOUBLE
  int  dtc=ecpdtFLOAT;
#else
  int  dtc=ecpdtDOUBLE;
#endif
  real *vv;
  int  nf,dt,i;

  nf = n;
  res = xdr_int(xd,&nf);
  if (nf != n && list == NULL)
    gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
	    dtc==ecpdtFLOAT ? "float" : "double",
	    dt ==ecpdtFLOAT ? "float" : "double");
  if (list) {
    snew(vv,n);
  } else {
    if (*v == NULL)
      snew(*v,n);
    vv = *v;
  }
  if (dt == ecpdtFLOAT) {
    res = xdr_vector(xd,(char *)vv,nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_float);
  } else {
    res = xdr_vector(xd,(char *)vv,nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  }
  
  if (list) {
    pr_reals(list,0,est_names[ecpt],vv,nf);
    sfree(vv);
  }
}

static do_cpte_real(XDR *xd,int ecpt,int sflags,real *r,FILE *list)
{
  do_cpte_reals(xd,ecpt,sflags,1,&r,list);
}

static void do_cpte_ints(XDR *xd,int ecpt,int sflags,int n,int **v,FILE *list)
{
  bool_t res=0;
  int  dtc=ecpdtINT;
  int *vv;
  int  nf,dt,i;

  nf = n;
  res = xdr_int(xd,&nf);
  if (nf != n && list == NULL)
    gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    gmx_fatal(FARGS,"Type mismatch for state entry %s, code type is %s, file type is %s\n",est_names[ecpt],ecpdt_names[dtc],ecpdt_names[dt]);
  if (list) {
    snew(vv,n);
  } else {
    if (*v == NULL)
      snew(*v,n);
    vv = *v;
  }
  res = xdr_vector(xd,(char *)vv,nf,
		   (unsigned int)sizeof(int),(xdrproc_t)xdr_int);
  
  if (list) {
    pr_ivec(list,0,est_names[ecpt],vv,nf);
    sfree(vv);
  }
}

static void do_cpte_int(XDR *xd,int ecpt,int sflags,int *i,FILE *list)
{
  do_cpte_ints(xd,ecpt,sflags,1,&i,list);
}

static void do_cpte_doubles(XDR *xd,int ecpt,int sflags,int n,double **v,
			    FILE *list)
{
  bool_t res=0;
  int  dtc=ecpdtDOUBLE;
  double *vv;
  int  nf,dt,i;

  nf = n;
  res = xdr_int(xd,&nf);
  if (nf != n && list == NULL)
    gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    gmx_fatal(FARGS,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
	    dtc==ecpdtFLOAT ? "float" : "double",
	    dt ==ecpdtFLOAT ? "float" : "double");
  if (list) {
    snew(vv,n);
  } else {
    if (*v == NULL)
      snew(*v,n);
    vv = *v;
  }
  res = xdr_vector(xd,(char *)vv,nf,
		   (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  
  if (list) {
    pr_doubles(list,0,est_names[ecpt],vv,nf);
    sfree(vv);
  }
}

static void do_cpte_rvecs(XDR *xd,int ecpt,int sflags,int n,rvec **v,
			  FILE *list)
{
  bool_t res=0;
#ifndef GMX_DOUBLE
  int  dtc=ecpdtFLOAT;
#else
  int  dtc=ecpdtDOUBLE;
#endif
  rvec *vv;
  int  nf,dt,i;

  nf = 3*n;
  res = xdr_int(xd,&nf);
  if (nf != 3*n && list == NULL)
    gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],3*n,nf);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
	    dtc==ecpdtFLOAT ? "float" : "double",
	    dt ==ecpdtFLOAT ? "float" : "double");
  if (list) {
    snew(vv,nf/3);
  } else {
    if (*v == NULL)
      snew(*v,nf/3);
    vv = *v;
  }
  if (dt == ecpdtFLOAT) {
    res = xdr_vector(xd,(char *)(vv[0]),nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_float);
  } else {
    res = xdr_vector(xd,(char *)(vv[0]),nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  }

  if (list) {
    pr_rvecs(list,0,est_names[ecpt],vv,nf/3);
    sfree(vv);
  }
}

static void do_cpte_matrix(XDR *xd,int ecpt,int sflags,matrix v,FILE *list)
{
  bool_t res=0;
#ifndef GMX_DOUBLE
  int  dtc=ecpdtFLOAT;
#else
  int  dtc=ecpdtDOUBLE;
#endif
  int  n,dt,i;

  n = 9;
  res = xdr_int(xd,&n);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    fprintf(stderr,"Precision mismatch for state entry %s, code precision is %s, file precision is %s\n",est_names[ecpt],
	    dtc==ecpdtFLOAT ? "float" : "double",
	    dt ==ecpdtFLOAT ? "float" : "double");
  if (dt == ecpdtFLOAT) {
    res = xdr_vector(xd,(char *)(v[0]),n,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_float);
  } else {
    res = xdr_vector(xd,(char *)(v[0]),n,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  }

  if (list) {
    pr_rvecs(list,0,est_names[ecpt],v,DIM);
  }
}

static int do_checkpoint(XDR *xd,bool bRead,
			 char **version,char **btime,char **buser,char **bmach,
			 char **ftime,
			 int *eIntegrator,int *step,double *t,
			 int *nnodes,int *dd_nc,int *npme,
			 t_state *state,
			 FILE *list)
{
  int  magic,file_version;
  int  natoms,ngtc,fflags,sflags,idum;
  int  i;

  magic = CPT_MAGIC;
  xdr_int(xd,&magic);
  if (magic != CPT_MAGIC) {
    fprintf(stderr,
	    "Magic number mismatch, checkpoint file has %d, should be %d\n",
	    magic,CPT_MAGIC);
    return -1;
  }
  do_cpt_string(xd,bRead,"GROMACS version"           ,version,list);
  do_cpt_string(xd,bRead,"GROMACS build time"        ,btime,list);
  do_cpt_string(xd,bRead,"GROMACS build user"        ,buser,list);
  do_cpt_string(xd,bRead,"GROMACS build machine"     ,bmach,list);
  do_cpt_string(xd,bRead,"checkpoint generation time",ftime,list);
  file_version = cpt_version;
  do_cpt_int(xd,"checkpoint file version",&file_version,list);
  if (file_version > cpt_version)
    gmx_fatal(FARGS,"Attempting to read a checkpoint file of version %d with code of version %d\n",file_version,cpt_version);

  natoms = state->natoms;
  do_cpt_int(xd,"#atoms",&natoms,list);
  if (state->natoms == -1) {
    state->natoms = natoms;
  } else if (natoms != state->natoms) {
    gmx_fatal(FARGS,"Checkpoint file is for a system of %d atoms, while the current system consists of %d atoms",natoms,state->natoms);
  }
  ngtc = state->ngtc;
  do_cpt_int(xd,"#T-coupling groups",&ngtc,list);
  if (state->ngtc == -1) {
    state->ngtc = ngtc;
  } else if (ngtc != state->ngtc) {
    gmx_fatal(FARGS,"Checkpoint file is for a system of %d T-coupling groups, while the current system consists of %d T-coupling groups",ngtc,state->ngtc);
  }
  do_cpt_int(xd,"integrator",eIntegrator,list);
  do_cpt_int(xd,"step",step,list);
  do_cpt_double(xd,"t",t,list);
  do_cpt_int(xd,"#nodes",nnodes,list);
  idum = 1;
  do_cpt_int(xd,"dd_nc[x]",dd_nc ? &(dd_nc[0]) : &idum,list);
  do_cpt_int(xd,"dd_nc[y]",dd_nc ? &(dd_nc[1]) : &idum,list);
  do_cpt_int(xd,"dd_nc[z]",dd_nc ? &(dd_nc[2]) : &idum,list);
  do_cpt_int(xd,"#PME-only nodes",npme,list);
  if (!bRead)
    fflags = state->flags;
  do_cpt_int(xd,"state flags",&fflags,list);
  sflags = state->flags;
  for(i=0; i<estNR; i++) {
    if (fflags & (1<<i)) {
      switch (i) {
      case estLAMBDA: do_cpte_real  (xd,i,sflags,&state->lambda,list); break;
      case estBOX:    do_cpte_matrix(xd,i,sflags,state->box,list); break;
      case estBOXV:   do_cpte_matrix(xd,i,sflags,state->boxv,list); break;
      case estPC_MU:  do_cpte_matrix(xd,i,sflags,state->pcoupl_mu,list); break;
      case estNH_XI:  do_cpte_reals (xd,i,sflags,state->ngtc,&state->nosehoover_xi,list); break;
      case estNH_IXI: do_cpte_doubles(xd,i,sflags,state->ngtc,&state->nosehoover_ixi,list); break;
      case estX:      do_cpte_rvecs (xd,i,sflags,natoms,&state->x,list); break;
      case estV:      do_cpte_rvecs (xd,i,sflags,natoms,&state->v,list); break;
      case estSDX:    do_cpte_rvecs (xd,i,sflags,natoms,&state->sd_X,list); break;
      case estLD_RNG: do_cpte_ints  (xd,i,sflags,state->nrng,(int **)&state->ld_rng,list); break;
      case estLD_RNGI: do_cpte_int  (xd,i,sflags,&state->ld_rngi,list); break;
      default:
	gmx_fatal(FARGS,"Unknown state entry %d",i);
      }
    }
  }

  return 0;
}

void write_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
		      int eIntegrator,int step,double t,t_state *state)
{
  int fp;
  char *version=VERSION;
  char *btime=BUILD_TIME;
  char *buser=BUILD_USER;
  char *bmach=BUILD_MACHINE;
  char *ftime;
  time_t now;
  int nnodes;
  char buf[1024];

  if (PAR(cr)) {
    if (DOMAINDECOMP(cr)) {
      nnodes = cr->dd->nnodes;
    } else {
      nnodes = cr->nnodes;
      gmx_fatal(FARGS,"write_checkpoint not supported with particle decomposition");
    }
  } else {
    nnodes = 1;
  }

  if (fexist(fn)) {
    /* Rename the previous checkpoint file */
    strcpy(buf,fn);
    buf[strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1] = '\0';
    strcat(buf,"_prev");
    strcat(buf,fn+strlen(fn) - strlen(ftp2ext(fn2ftp(fn))) - 1);
    rename(fn,buf);
  }

  now = time(NULL);
  ftime = strdup(ctime(&now));
  ftime[strlen(ftime)-1] = '\0';

  fprintf(stderr,"\nWriting checkpoint, step %d at %s\n",step,ftime);
  if (fplog)
    fprintf(fplog,"Writing checkpoint, step %d at %s\n\n",step,ftime);

  fp = gmx_fio_open(fn,"w");
  do_checkpoint(gmx_fio_getxdr(fp),FALSE,&version,&btime,&buser,&bmach,&ftime,
		&eIntegrator,&step,&t,&nnodes,
		DOMAINDECOMP(cr) ? cr->dd->nc : NULL,&cr->npmenodes,
		state,NULL);
  gmx_fio_close(fp);

  sfree(ftime);
}

void read_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
		     int *step,double *t,
		     int *nnodes,ivec dd_nc,int *npme,
		     t_state *state)
{
  int fp;
  char *version,*btime,*buser,*bmach,*ftime;
  int eIntegrator;

  if (PARTDECOMP(cr))
    gmx_fatal(FARGS,
	      "read_checkpoint not supported with particle decomposition");

  fp = gmx_fio_open(fn,"r");
  do_checkpoint(gmx_fio_getxdr(fp),TRUE,&version,&btime,&buser,&bmach,&ftime,
		&eIntegrator,step,t,nnodes,dd_nc,npme,state,NULL);
  gmx_fio_close(fp);

  if (cr == NULL || MASTER(cr))
    fprintf(stderr,"\nRead checkpoint file %s generated: %s\n\n",
	    fn,ftime);
  if (fplog) {
    fprintf(fplog,"\n");
    fprintf(fplog,"Read checkpoint file %s\n",fn);
    fprintf(fplog,"  file generated:        %s\n",ftime);  
    fprintf(fplog,"  GROMACS build time:    %s\n",btime);  
    fprintf(fplog,"  GROMACS build user:    %s\n",buser);  
    fprintf(fplog,"  GROMACS build machine: %s\n",bmach);  
    fprintf(fplog,"  step %d\n",*step);  
    fprintf(fplog,"  time %f\n",*t);  
    fprintf(fplog,"\n");
  }
  sfree(ftime);
  sfree(btime);
  sfree(buser);
  sfree(bmach);
}

void list_checkpoint(char *fn,FILE *out)
{
  int fp;
  char *version,*btime,*buser,*bmach,*ftime;
  int eIntegrator,step,nnodes,npme;
  double t;
  ivec dd_nc;
  t_state state;
  int indent;
  int i;

  init_state(&state,-1,-1);

  fp = gmx_fio_open(fn,"r");
  do_checkpoint(gmx_fio_getxdr(fp),TRUE,&version,&btime,&buser,&bmach,&ftime,
		&eIntegrator,&step,&t,&nnodes,dd_nc,&npme,&state,out);
  gmx_fio_close(fp);

  done_state(&state);
}
