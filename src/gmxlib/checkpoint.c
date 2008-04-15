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
#include "statutil.h"
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
  real *vp,*va=NULL;
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
  if (list || !(sflags & (1<<ecpt))) {
    snew(va,n);
    vp = va;
  } else {
    if (*v == NULL)
      snew(*v,n);
    vp = *v;
  }
  if (dt == ecpdtFLOAT) {
    res = xdr_vector(xd,(char *)vp,nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_float);
  } else {
    res = xdr_vector(xd,(char *)vp,nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  }
  
  if (list)
    pr_reals(list,0,est_names[ecpt],vp,nf);
  if (va)
    sfree(va);
}

static do_cpte_real(XDR *xd,int ecpt,int sflags,real *r,FILE *list)
{
  do_cpte_reals(xd,ecpt,sflags,1,&r,list);
}

static void do_cpte_ints(XDR *xd,int ecpt,int sflags,int n,int **v,FILE *list)
{
  bool_t res=0;
  int  dtc=ecpdtINT;
  int *vp,*va=NULL;
  int  nf,dt,i;

  nf = n;
  res = xdr_int(xd,&nf);
  if (nf != n && list == NULL)
    gmx_fatal(FARGS,"Count mismatch for state entry %s, code count is %d, file count is %d\n",est_names[ecpt],n,nf);
  dt = dtc;
  res = xdr_int(xd,&dt);
  if (dt != dtc)
    gmx_fatal(FARGS,"Type mismatch for state entry %s, code type is %s, file type is %s\n",est_names[ecpt],ecpdt_names[dtc],ecpdt_names[dt]);
  if (list || !(sflags & (1<<ecpt))) {
    snew(va,nf);
    vp = va;
  } else {
    if (*v == NULL)
      snew(*v,nf);
    vp = *v;
  }
  res = xdr_vector(xd,(char *)vp,nf,
		   (unsigned int)sizeof(int),(xdrproc_t)xdr_int);
  
  if (list)
    pr_ivec(list,0,est_names[ecpt],vp,nf);
  if (va)
    sfree(va);
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
  double *vp,*va=NULL;
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
  if (list || !(sflags & (1<<ecpt))) {
    snew(va,n);
    vp = va;
  } else {
    if (*v == NULL)
      snew(*v,n);
    vp = *v;
  }
  res = xdr_vector(xd,(char *)vp,nf,
		   (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  
  if (list)
    pr_doubles(list,0,est_names[ecpt],vp,nf);
  if (va)
    sfree(va);
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
  rvec *vp,*va=NULL;
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
  if (list || !(sflags & (1<<ecpt))) {
    snew(va,nf/3);
    vp = va;
  } else {
    if (*v == NULL)
      snew(*v,nf/3);
    vp = *v;
  }
  if (dt == ecpdtFLOAT) {
    res = xdr_vector(xd,(char *)(vp[0]),nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_float);
  } else {
    res = xdr_vector(xd,(char *)(vp[0]),nf,
		     (unsigned int)sizeof(real),(xdrproc_t)xdr_double);
  }

  if (list)
    pr_rvecs(list,0,est_names[ecpt],vp,nf/3);
  if (va)
    sfree(va);
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

static int do_cpt_header(XDR *xd,bool bRead,
			 char **version,char **btime,char **buser,char **bmach,
			 char **fprog,char **ftime,
			 int *eIntegrator,int *step,double *t,
			 int *nnodes,int *dd_nc,int *npme,
			 int *natoms,int *ngtc,int *flags,
			 FILE *list)
{
  int  magic,file_version;
  int  idum;
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
  do_cpt_string(xd,bRead,"generating program"        ,fprog,list);
  do_cpt_string(xd,bRead,"generation time"           ,ftime,list);
  file_version = cpt_version;
  do_cpt_int(xd,"checkpoint file version",&file_version,list);
  if (file_version > cpt_version)
    gmx_fatal(FARGS,"Attempting to read a checkpoint file of version %d with code of version %d\n",file_version,cpt_version);

  do_cpt_int(xd,"#atoms"            ,natoms     ,list);
  do_cpt_int(xd,"#T-coupling groups",ngtc       ,list);
  do_cpt_int(xd,"integrator"        ,eIntegrator,list);
  do_cpt_int(xd,"step"              ,step       ,list);
  do_cpt_double(xd,"t"              ,t          ,list);
  do_cpt_int(xd,"#nodes"            ,nnodes     ,list);
  idum = 1;
  do_cpt_int(xd,"dd_nc[x]",dd_nc ? &(dd_nc[0]) : &idum,list);
  do_cpt_int(xd,"dd_nc[y]",dd_nc ? &(dd_nc[1]) : &idum,list);
  do_cpt_int(xd,"dd_nc[z]",dd_nc ? &(dd_nc[2]) : &idum,list);
  do_cpt_int(xd,"#PME-only nodes",npme,list);
  do_cpt_int(xd,"state flags",flags,list);

  return 0;
}

static int do_cpt_state(XDR *xd,bool bRead,
			int fflags,int nppnodes,t_state *state,
			FILE *list)
{
  int  sflags;
  int  i;

  sflags = state->flags;
  for(i=0; i<estNR; i++) {
    if (fflags & (1<<i)) {
      switch (i) {
      case estLAMBDA: do_cpte_real  (xd,i,sflags,&state->lambda,list); break;
      case estBOX:    do_cpte_matrix(xd,i,sflags,state->box,list); break;
      case estBOX_REL:do_cpte_matrix(xd,i,sflags,state->box_rel,list); break;
      case estBOXV:   do_cpte_matrix(xd,i,sflags,state->boxv,list); break;
      case estPC_MU:  do_cpte_matrix(xd,i,sflags,state->pcoupl_mu,list); break;
      case estNH_XI:  do_cpte_reals (xd,i,sflags,state->ngtc,&state->nosehoover_xi,list); break;
      case estNH_IXI: do_cpte_doubles(xd,i,sflags,state->ngtc,&state->nosehoover_ixi,list); break;
      case estX:      do_cpte_rvecs (xd,i,sflags,state->natoms,&state->x,list); break;
      case estV:      do_cpte_rvecs (xd,i,sflags,state->natoms,&state->v,list); break;
      case estSDX:    do_cpte_rvecs (xd,i,sflags,state->natoms,&state->sd_X,list); break;
      case estLD_RNG: do_cpte_ints  (xd,i,sflags,state->nrng,(int **)&state->ld_rng,list); break;
      case estLD_RNGI: do_cpte_ints (xd,i,sflags,nppnodes,&state->ld_rngi,list); break;
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
  char *fprog;
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

  fprog = Program();
  now = time(NULL);
  ftime = strdup(ctime(&now));
  ftime[strlen(ftime)-1] = '\0';

  fprintf(stderr,"\nWriting checkpoint, step %d at %s\n",step,ftime);
  if (fplog)
    fprintf(fplog,"Writing checkpoint, step %d at %s\n\n",step,ftime);

  fp = gmx_fio_open(fn,"w");
  do_cpt_header(gmx_fio_getxdr(fp),FALSE,
		&version,&btime,&buser,&bmach,&fprog,&ftime,
		&eIntegrator,&step,&t,&nnodes,
		DOMAINDECOMP(cr) ? cr->dd->nc : NULL,&cr->npmenodes,
		&state->natoms,&state->ngtc,&state->flags,NULL);
  do_cpt_state(gmx_fio_getxdr(fp),FALSE,state->flags,
	       cr->nnodes-cr->npmenodes,state,NULL);
  gmx_fio_close(fp);

  sfree(ftime);
}

static void print_flag_mismatch(FILE *fplog,int sflags,int fflags)
{
  int i;

  fprintf(fplog,"\nState entry mismatch between the simulation and the checkpoint file\n");
  fprintf(fplog,"Entries which are not present in the checkpoint file will not be updated\n");
  fprintf(fplog,"  %24s    %11s    %11s\n","","simulation","checkpoint");
  for(i=0; i<estNR; i++) {
    if ((sflags & (1<<i)) || (fflags & (1<<i))) {
      fprintf(fplog,"  %24s    %11s    %11s\n",
	      est_names[i],
	      (sflags & (1<<i)) ? "  present  " : "not present",
	      (fflags & (1<<i)) ? "  present  " : "not present");
    }
  }
}

static void check_string(FILE *fplog,char *type,char *p,char *f,bool *mm)
{
  if (strcmp(p,f) != 0) {
    if (fplog) {
      fprintf(fplog,"  %s mismatch,\n",type);
      fprintf(fplog,"    current program: %s\n",p);
      fprintf(fplog,"    checkpoint file: %s\n",f);
      fprintf(fplog,"\n");
    }
    *mm = TRUE;
  }
}

static void check_build_match(FILE *fplog,
			      char *version,
			      char *btime,char *buser,char *bmach,
			      char *fprog)
{
  bool mm;

  mm = FALSE;
  check_string(fplog,"Version"     ,VERSION      ,version,&mm);
  check_string(fplog,"Build time"  ,BUILD_TIME   ,btime  ,&mm);
  check_string(fplog,"Build user"  ,BUILD_USER   ,buser  ,&mm);
  check_string(fplog,"Build time"  ,BUILD_MACHINE,bmach  ,&mm);
  check_string(fplog,"Program name",Program()    ,fprog  ,&mm);
  if (mm) {
    if (fplog)
      fprintf(fplog,"Continuation is exact, but is not guaranteed to be binary identical\n\n");
    fprintf(stderr,"Continuation is exact, but is not guaranteed to be binary identical,\n"
	    "see the log file for details\n\n");
  }
}

void read_checkpoint(char *fn,FILE *fplog,t_commrec *cr,
		     int *step,double *t,
		     int *nnodes,ivec dd_nc,int *npme,
		     t_state *state)
{
  int  fp;
  char *version,*btime,*buser,*bmach,*fprog,*ftime;
  int  eIntegrator;
  int  natoms,ngtc,fflags;

  if (PARTDECOMP(cr))
    gmx_fatal(FARGS,
	      "read_checkpoint not supported with particle decomposition");

  fp = gmx_fio_open(fn,"r");
  do_cpt_header(gmx_fio_getxdr(fp),TRUE,
		&version,&btime,&buser,&bmach,&fprog,&ftime,
		&eIntegrator,step,t,nnodes,dd_nc,npme,
		&natoms,&ngtc,&fflags,NULL);

  if (cr == NULL || MASTER(cr)) {
    fprintf(stderr,"\nReading checkpoint file %s generated: %s\n\n",
	    fn,ftime);
  }
  if (fplog) {
    fprintf(fplog,"\n");
    fprintf(fplog,"Reading checkpoint file %s\n",fn);
    fprintf(fplog,"  file generated by:     %s\n",fprog);  
    fprintf(fplog,"  file generated at:     %s\n",ftime);  
    fprintf(fplog,"  GROMACS build time:    %s\n",btime);  
    fprintf(fplog,"  GROMACS build user:    %s\n",buser);  
    fprintf(fplog,"  GROMACS build machine: %s\n",bmach);  
    fprintf(fplog,"  step %d\n",*step);  
    fprintf(fplog,"  time %f\n",*t);  
    fprintf(fplog,"\n");
  }

  if (natoms != state->natoms) {
    gmx_fatal(FARGS,"Checkpoint file is for a system of %d atoms, while the current system consists of %d atoms",natoms,state->natoms);
  }
  if (ngtc != state->ngtc) {
    gmx_fatal(FARGS,"Checkpoint file is for a system of %d T-coupling groups, while the current system consists of %d T-coupling groups",ngtc,state->ngtc);
  }
  if (fflags != state->flags) {
    if (MASTER(cr))
      fprintf(stderr,
	      "WARNING: The checkpoint state entries do not match the simulation,\n"
	      "         see the log file for details\n");
    print_flag_mismatch(fplog,state->flags,fflags);
  } else {
    if (MASTER(cr))
      check_build_match(fplog,version,btime,buser,bmach,fprog);
  }
  do_cpt_state(gmx_fio_getxdr(fp),TRUE,fflags,(*nnodes)-(*npme),state,NULL);
  gmx_fio_close(fp);

  sfree(fprog);
  sfree(ftime);
  sfree(btime);
  sfree(buser);
  sfree(bmach);
}

void read_checkpoint_state(char *fn,int *step,double *t,t_state *state)
{
  int  fp;
  char *version,*btime,*buser,*bmach,*fprog,*ftime;
  int  eIntegrator;
  int  nnodes,npme;
  ivec dd_nc;

  fp = gmx_fio_open(fn,"r");
  do_cpt_header(gmx_fio_getxdr(fp),TRUE,
		&version,&btime,&buser,&bmach,&fprog,&ftime,
		&eIntegrator,step,t,&nnodes,dd_nc,&npme,
		&state->natoms,&state->ngtc,&state->flags,NULL);
  do_cpt_state(gmx_fio_getxdr(fp),TRUE,state->flags,nnodes-npme,state,NULL);
  gmx_fio_close(fp);

  sfree(fprog);
  sfree(ftime);
  sfree(btime);
  sfree(buser);
  sfree(bmach);
}

void list_checkpoint(char *fn,FILE *out)
{
  int  fp;
  char *version,*btime,*buser,*bmach,*fprog,*ftime;
  int  eIntegrator,step,nnodes,npme;
  double t;
  ivec dd_nc;
  t_state state;
  int  indent;
  int  i;

  init_state(&state,-1,-1);

  fp = gmx_fio_open(fn,"r");
  do_cpt_header(gmx_fio_getxdr(fp),TRUE,
		&version,&btime,&buser,&bmach,&fprog,&ftime,
		&eIntegrator,&step,&t,&nnodes,dd_nc,&npme,
		&state.natoms,&state.ngtc,&state.flags,out);
  do_cpt_state(gmx_fio_getxdr(fp),TRUE,state.flags,nnodes-npme,&state,out);
  gmx_fio_close(fp);

  done_state(&state);
}
