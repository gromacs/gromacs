/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_edsam_c = "$Id$";

#include <stdio.h>
#include <time.h>
#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "names.h"
#include "stat.h"
#include "confio.h"
#include "mvdata.h"
#include "pdebug.h"
#include "txtdump.h"
#include "vec.h"
#include "time.h"
#include "nrnb.h"
#include "mshift.h"
#include "mdrun.h"
#include "update.h"
#include "physics.h"
#include "rmpbc.h"
#include "gstat.h"
#include "edsam.h"

#define EPS  1.0e-09

void ed_open(int nfile,t_filenm fnm[],t_edsamyn *edyn)
{
  fprintf(stderr,"ED sampling will be performed!\n");
  edyn->bEdsam=TRUE;
  edyn->edinam=ftp2fn(efEDI,nfile,fnm);
  edyn->edonam=ftp2fn(efEDO,nfile,fnm); 
}

void init_edsam(FILE *log,t_topology *top,
	   t_mdatoms *md,int start,int homenr,
	   int *nbl,int **sbl,rvec x[],matrix box, 
	   t_edsamyn *edyn,t_edpar *edi)
{
  int i,j,ned,*refmasnrs;
  rvec *xdum,*transvec;
  real rmsd;
  matrix rotmat;

  fprintf(log,"Initialising ED sampling\n\n");
  
  /* first read the input. All input is stored in edi */
  read_edi(edyn,edi,top->atoms.nr);

  /* here for initialisation */

  /* evaluate masses */
  ned=edi->ned;
  edi->tmass=0.0;
  if (edi->selmas) {
    edi->nmass=edi->npro;
    snew(edi->mass,edi->nmass);
    snew(edi->masnrs,edi->nmass);
    snew(refmasnrs,edi->nmass);
    for(i=0; (i < edi->nmass); i++) {
      edi->mass[i]=top->atoms.atom[edi->sref.anrs[i]].m;
      edi->masnrs[i]=i;
      refmasnrs[i]=i;
      edi->tmass+=edi->mass[i];
    }
  }
  else {
    edi->nmass=edi->sav.nr;
    snew(edi->mass,edi->nmass);
    snew(edi->masnrs,edi->nmass);
    snew(refmasnrs,edi->nmass);
    for(i=0; (i < edi->nmass); i++) {
      edi->mass[i]=1.0;
      edi->masnrs[i]=edi->sav.anrs[i];
      for(j=0; (j < edi->sref.nr); j++) {
	if (edi->sref.anrs[j] == edi->sav.anrs[i])
	  refmasnrs[i]=j;
      }
      edi->tmass+=edi->mass[i];
    }
  }

  /* mark atoms that are to be used for rotational fit */
  edi->nfit=0;
  for(i=0; (i < ned); i++) {
    if (strcmp(*top->atoms.atomname[i],"CA") == 0) 
      edi->nfit++;
  }
  fprintf(log,"%d CA atoms found\n", edi->nfit);
  if (edi->nfit == 0) {
    edi->nfit = edi->sav.nr;
    snew(edi->fitnrs,edi->nfit);
    for(i=0; (i < edi->nfit); i++) {
      edi->fitnrs[i] = edi->sav.anrs[i];
    }
  }
  else {
    snew(edi->fitnrs,edi->nfit);
    j = 0;
    for(i=0; (i < ned); i++) {
      if (strcmp(*top->atoms.atomname[i],"CA") == 0) {
	edi->fitnrs[j] = i;
	/*fprintf(stderr,"Found CA at position %d\n",i);*/
	j++;
      }
    }
  }

  /* put reference structure in origin */
  put_in_origin(edi->sref.nr,edi->sref.x,edi->nmass,refmasnrs,edi->mass,
		edi->tmass);

  /* reduce the reference structure */
  snew(xdum,edi->nfit);
  j=0;
  for(i=0; (i < edi->sref.nr); i++) {
    if (j < edi->nfit) {
      if (edi->sref.anrs[i] == edi->fitnrs[j]) {
	copy_rvec(edi->sref.x[i],xdum[j]);
	j++;
      }
    }
  }
  if (j != edi->nfit) fatal_error(0,"Counted %d, should be %d in init_edsam",
				  j,edi->nfit);
  edi->sref.nr=edi->nfit;
  sfree(edi->sref.x);
  snew(edi->sref.x,edi->nfit);
  for(i=0; (i < edi->nfit); i++)
    copy_rvec(xdum[i],edi->sref.x[i]);
  sfree(xdum);

  /* remove pbc */
  snew(xdum,top->atoms.nr);
  rm_pbc(&(top->idef),top->atoms.nr,box,x,xdum);

  /* fit starting positions to reference structure */
  snew(transvec,ned);
  rmsd=fitit(ned,xdum,edi,transvec,rotmat);
  fprintf(log,"Initial RMSD from reference structure = %10.5f nm\n\n",rmsd);
  sfree(transvec);

  /* calculate initial projections */
  project(xdum,edi,"x");
  fprintf(log,"Initial projections:\n");
  write_edidx(log,edi);

  /* process target structure, if required */
  if (edi->star.nr > 0) {
    rmsd=fitit(ned,edi->star.x,edi,transvec,rotmat);
    projectx(edi,edi->star.x,&edi->vecs.radcon);
  }

  /* process structure that will serve as origin of expansion circle */
  if (edi->sori.nr > 0) {
    rmsd=fitit(ned,edi->sori.x,edi,transvec,rotmat);
    projectx(edi,edi->sori.x,&edi->vecs.radacc);
    projectx(edi,edi->sori.x,&edi->vecs.radfix);
  }
  else {
    projectx(edi,xdum,&edi->vecs.radacc);
    projectx(edi,xdum,&edi->vecs.radfix);
  }
  
  /* set starting projections for linsam */
  projectx(edi,xdum,&edi->vecs.linacc);
  projectx(edi,xdum,&edi->vecs.linfix);
  sfree(xdum);

  /* calculate initial radii */
  fprintf(log,"Initial fixed increment radius=%f\n",edi->vecs.radfix.radius);
  fprintf(log,"Initial   acceptance    radius=%f\n",edi->vecs.radacc.radius);
  fprintf(log,"Initial   contracting   radius=%f\n",edi->vecs.radcon.radius);
  fflush(log);
  
  /* open output file */
  edi->edo=ffopen(edyn->edonam,"w");
  write_edidx(edi->edo,edi);
}

void read_edi(t_edsamyn *edyn,t_edpar *edi,int nr_mdatoms)
{
  FILE *in;
  int i,j,idum,magic=666,readmagic;
  rvec *xdum;

  /* the edi file is not free format, so expect problems if the input is 
     corrupt :-) */
  in=ffopen(edyn->edinam,"r");

  /* check the magic number */
  readmagic=read_edint(in);
  if (readmagic != magic)
    fatal_error(0,"Wrong magic number in %s",edyn->edinam);
  
  /* check the number of atoms */
   edi->nini=read_edint(in);
  if (edi->nini != nr_mdatoms)
    fatal_error(0,"Nr of atoms in %s (%d) does not match nr of md atoms (%d)",
		edyn->edinam,edi->nini,nr_mdatoms); 

  /* Done checking. For the rest we blindly trust the input */
  edi->npro=read_edint(in);
  idum=read_edint(in);
  if (idum == 0) edi->selmas=FALSE;
  else edi->selmas=TRUE;
  edi->outfrq=read_edint(in);
  edi->logfrq=read_edint(in);
  edi->maxedsteps=read_edint(in);
  edi->slope=read_edreal(in);
  edi->sref.nr=read_edint(in);

  /* allocate space for reference positions and read them */
  snew(edi->sref.anrs,edi->sref.nr);
  snew(edi->sref.x,edi->sref.nr);
  read_edx(in,edi->sref.nr,edi->sref.anrs,edi->sref.x);

  /* average positions. they define which atoms will be used for ED sampling */
  edi->sav.nr=read_edint(in);
  snew(edi->sav.anrs,edi->sav.nr);
  snew(edi->sav.x,edi->sav.nr);
  read_edx(in,edi->sav.nr,edi->sav.anrs,edi->sav.x);

  /* Nr of essdyn atoms */
  edi->ned=read_edint(in);

  /* eigenvectors */
  read_edvecs(in,edi->sav.nr,&edi->vecs);

  /* target positions */
  edi->star.nr=read_edint(in);
  if (edi->star.nr > 0) {
    snew(edi->star.anrs,edi->star.nr);
    snew(xdum,edi->star.nr);
    read_edx(in,edi->star.nr,edi->star.anrs,xdum);
    snew(edi->star.x,edi->ned);
    j=0;
    for(i=0; (i < edi->ned); i++) {
      if (edi->star.anrs[j] == i) {
	copy_rvec(xdum[j],edi->star.x[i]);
	j++;
      }
    }
    sfree(xdum);
  }

  /* positions defining origin of expansion circle */
  edi->sori.nr=read_edint(in);
  if (edi->sori.nr > 0) {
    snew(edi->sori.anrs,edi->sori.nr);
    snew(xdum,edi->sori.nr);
    read_edx(in,edi->sori.nr,edi->sori.anrs,xdum);
    snew(edi->sori.x,edi->ned);
    j=0;
    for(i=0; (i < edi->ned); i++) {
      if (edi->sori.anrs[j] == i) {
	copy_rvec(xdum[j],edi->sori.x[i]);
	j++;
      }
    }
    sfree(xdum);
  }
  
  /* all done */
  ffclose (in);
}

int read_edint(FILE *file)
{
  char line[STRLEN+1];
  int idum;

  fgets2 (line,STRLEN,file);
  fgets2 (line,STRLEN,file);
  sscanf (line,"%d",&idum);
  return idum;
}

real read_edreal(FILE *file)
{
  char line[STRLEN+1];
  double rdum;

  fgets2 (line,STRLEN,file);
  fgets2 (line,STRLEN,file);
  sscanf (line,"%lf",&rdum);
  return (real) rdum; /* always read as double and convert to single */
}

int read_edint2(FILE *file)
{
  char line[STRLEN+1];
  int idum;

  fgets2 (line,STRLEN,file);
  sscanf (line,"%d",&idum);
  return idum;
}

void read_edx(FILE *file,int number,int *anrs,rvec *x)
{
  int i,j;
  char line[STRLEN+1];
  double d[3];

  for(i=0; (i < number); i++) {
    fgets2 (line,STRLEN,file);
    sscanf (line,"%d%lf%lf%lf",&anrs[i],&d[0],&d[1],&d[2]);
    anrs[i]--; /* we are reading FORTRAN indices */
    for(j=0; (j < 3); j++)
      x[i][j]=d[j]; /* always read as double and convert to single */
  }
}

void read_edvecs(FILE *in,int nr,t_edvecs *vecs)
{
  read_edvec(in,nr,&vecs->mon);
  read_edvec(in,nr,&vecs->linfix);
  read_edvec(in,nr,&vecs->linacc);
  read_edvec(in,nr,&vecs->radfix);
  read_edvec(in,nr,&vecs->radacc);
  read_edvec(in,nr,&vecs->radcon);
}

void read_edvec(FILE *in,int nr,t_eigvec *tvec)
{
  int i,idum;
  double rdum;
  char line[STRLEN+1];

  tvec->neig=read_edint(in);
  snew(tvec->ieig,tvec->neig);
  snew(tvec->stpsz,tvec->neig);
  snew(tvec->vec,tvec->neig);
  snew(tvec->xproj,tvec->neig);
  snew(tvec->vproj,tvec->neig);
  snew(tvec->fproj,tvec->neig);
  snew(tvec->refproj,tvec->neig);
  for(i=0; (i < tvec->neig); i++) {
    fgets2 (line,STRLEN,in);
    sscanf (line,"%d%lf",&idum,&rdum);
    tvec->ieig[i]=idum;
    tvec->stpsz[i]=rdum;
  }
  for(i=0; (i < tvec->neig); i++) {
    snew(tvec->vec[i],nr);
    scan_edvec(in,nr,tvec->vec[i]);
  }
}
  
void scan_edvec(FILE *in,int nr,rvec *vec)
{
  char line[STRLEN+1];
  int i;
  double x,y,z;

  for(i=0; (i < nr); i++) {
    fgets2 (line,STRLEN,in);
    sscanf (line,"%le%le%le",&x,&y,&z);
    vec[i][0]=x;
    vec[i][1]=y;
    vec[i][2]=z;
  }
}
	  
real fitit(int nr, rvec *x,t_edpar *edi,rvec *transvec,matrix rmat)
{
  rvec *xdum,x_old;
  int i,j,k;
  real rmsd;

  snew(xdum,nr);
  for(i=0; (i<nr); i++)
    copy_rvec(x[i],xdum[i]); 

  /* first do translational fit */
  put_in_origin(nr,x,edi->nmass,edi->masnrs,edi->mass,edi->tmass);
  for(i=0; (i<nr); i++) 
    rvec_sub(x[i],xdum[i],transvec[i]);
  sfree(xdum);

  /* now rotational fit */
  snew(xdum,edi->nfit);
  for(i=0; (i<edi->nfit); i++)
    copy_rvec(x[edi->fitnrs[i]],xdum[i]);
  do_edfit(edi->nfit,edi->sref.x,xdum,rmat);
  sfree(xdum);

  /* apply the rotation matrix */
  for(i=0;(i<nr);i++) {
    for(j=0;(j<3);j++)
      x_old[j]=x[i][j];
    for(j=0;(j<3);j++) {
      x[i][j]=0;
      for(k=0;(k<3);k++)
        x[i][j]+=rmat[k][j]*x_old[k];
    }
  }

  /* calculate RMSD */
  rmsd=0.0;
  for(i=0; (i<edi->nfit); i++) {
    for(j=0; (j<DIM); j++)
      rmsd+=pow((edi->sref.x[i][j]-x[edi->fitnrs[i]][j]),2);
  }
  rmsd/=(real) edi->nfit;
  rmsd=sqrt(rmsd);
  /*fprintf(stderr,"RMSD to reference structure=%15.10f\n",rmsd);*/
  return rmsd;
}

void do_edfit(int natoms,rvec *xp,rvec *x,matrix R)
{
  /* this is a copy of do_fit with some modifications */
  int    c,r,n,j,i,irot;
  double omega[7][7],om[7][7],d[7],xnr,xpc;
  matrix vh,vk,vh_d,vk_d,u;
  real   du,sigd;
  int    index;
  real   max_d;

  for(i=0;(i<7);i++) {
    d[i]=0;
    for(j=0;(j<7);j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }

  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++) {
    for(c=0; (c<DIM); c++) {
      xpc=xp[n][c];
      for(r=0; (r<DIM); r++) {
	xnr=x[n][r];
	u[c][r]+=xnr*xpc;
      }
    }
  }
  
  /*calculate its determinant*/
  du=det(u);
    
  if (fabs(du)<EPS) {
    fprintf(stderr,"\n");
    pr_rvecs(stderr,0,"U",u,DIM);
    fatal_error(0,"Determinant of U = 0\n");
  }  
  sigd=du/fabs(du);

  /*construct omega*/
  /*omega is symmetric -> omega==omega' */
  for(r=0;(r<6);r++)
    for(c=0;(c<=r);c++)
      if ((r>=3) && (c<3)) {
        omega[r+1][c+1]=u[r-3][c];
        omega[c+1][r+1]=u[r-3][c];
      }
      else {
        omega[r+1][c+1]=0;
        omega[c+1][r+1]=0;
      }

  /*determine h and k*/
  jacobi(omega,6,d,om,&irot);

  if (irot==0) {
    fprintf(stderr,"IROT=0\n");
  }

  for(c=0;(c<3);c++)
    for(r=0;(r<3);r++) {
      vh_d[r][c]=om[c+4][r+1];
      vk_d[r][c]=om[c+4][r+4];
    }

  index=0; /* For the compiler only */

  for(j=0;(j<3);j++) {
    max_d=-1000;
    for(i=0;(i<6);i++)
      if (d[i+1]>max_d) {
        max_d=d[i+1];
        index=i;
      }
    d[index+1]=-10000;
    for(i=0;(i<3);i++) {
      vh[j][i]=M_SQRT2*om[i+1][index+1];
      vk[j][i]=M_SQRT2*om[i+4][index+1];
      vh_d[j][i]=om[i+1][index+1];
      vk_d[j][i]=om[i+4][index+1];
    }
  }

  /*determine R*/
  for(c=0;(c<3);c++)
    for(r=0;(r<3);r++)
      R[c][r]=vk[0][r]*vh[0][c]+
              vk[1][r]*vh[1][c]+sigd*
              vk[2][r]*vh[2][c];
}

void put_in_origin(int nr,rvec *x,int nmass,int *masnrs,real *mass,real tmass)
{
  int i,j;
  rvec cm;

  /* calculate CM */
  for (i=0;(i<DIM);i++) 
    cm[i]=0.0;	 
  for (i=0;(i<nmass);i++) {
    for (j=0;(j<DIM);j++)
      cm[j]+=x[masnrs[i]][j]*mass[i];
  }

  for (i=0;(i<DIM);i++) 
    cm[i]/=tmass;

  /* and subtract it */
  for (i=0;(i<nr);i++) 
    rvec_dec(x[i],cm);
}

void project(rvec *x,t_edpar *edi,char *mode)
{
  int i;

  /* subtract average positions */
  if (strcmp(mode,"x") == 0) {
    for (i=0;(i<edi->sav.nr);i++) 
      rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);
  }

  /* make the projections */
  do_project(x,&edi->vecs.mon,edi,mode);
  do_project(x,&edi->vecs.linfix,edi,mode);
  do_project(x,&edi->vecs.linacc,edi,mode);
  do_project(x,&edi->vecs.radfix,edi,mode);
  do_project(x,&edi->vecs.radacc,edi,mode);
  do_project(x,&edi->vecs.radcon,edi,mode);

  /* add average positions */
  if (strcmp(mode,"x") == 0) {
  for (i=0;(i<edi->sav.nr);i++) 
    rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
  }
}

void do_project(rvec *x,t_eigvec *vec,t_edpar *edi,char *mode)
{
  int i,j,k;
  real proj;

  for (i=0;(i<vec->neig);i++) { 
    proj=0.0;
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	proj+=vec->vec[i][j][k]*x[edi->sav.anrs[j]][k];
    }
    /*fprintf(stderr,"Projection on ev[%d]=%15.10f\n",vec->ieig[i],proj);*/
    if (strcmp(mode,"x") == 0)  vec->xproj[i]=proj;
    else if (strcmp(mode,"v") == 0) vec->vproj[i]=proj;
    else if (strcmp(mode,"f") == 0) vec->fproj[i]=proj;
  }
}

void projectx(t_edpar *edi,rvec *x,t_eigvec *vec)
{
  int i,j,k;
  real rad=0.0;

  /* subtract average positions */
  for (i=0;(i<edi->sav.nr);i++)
    rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);

  for (i=0;(i<vec->neig);i++) {
    vec->refproj[i]=0.0;
    for (j=0;(j<edi->sav.nr);j++) { 
      for (k=0;(k<DIM);k++)
	vec->refproj[i]+=vec->vec[i][j][k]*x[edi->sav.anrs[j]][k];
    }
    rad+=pow((vec->refproj[i]-vec->xproj[i]),2);
    /*fprintf(stderr,"Proj[%d]=%f\n",vec->ieig[i],vec->refproj[i]);*/
  }
  vec->radius=sqrt(rad);
    
  /* add average positions */
  for (i=0;(i<edi->sav.nr);i++) 
    rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
}

real do_projectx(t_edpar *edi,rvec *x,rvec *vec)
{
  int i,j;
  real proj=0.0;

  for (i=0;(i<edi->sav.nr);i++)
    for (j=0;(j<DIM);j++)
      proj+=vec[i][j]*x[edi->sav.anrs[i]][j];

  return proj;
}
	
real calc_radius(t_eigvec *vec)
{
  int i;
  
  vec->radius = 0.0;
  for (i=0;(i<vec->neig);i++) 
    vec->radius+=pow((vec->refproj[i]-vec->xproj[i]),2);
  return vec->radius=sqrt(vec->radius);

}

void do_edsam(FILE *log,t_topology *top,t_inputrec *ir,int step,
	      t_mdatoms *md,int start,int homenr,int *nbl,int **sbl,
              rvec x[],rvec xold[],rvec x_unc[],rvec force[],matrix box,
	      t_edsamyn *edyn,t_edpar *edi,bool bHave_force)
{
  int i,j,ned=edi->ned,edstep=step,iupdate=500;
  rvec *transvec,*vdum,*fdum;
  matrix rotmat;
  real rmsd,mas,rad;
  static real oldrad;
  static bool bFirst=TRUE;
  real dt   = ir->delta_t;
  real dt_1 = 1.0/dt;
  real dt_2 = 1.0/(dt*dt);

  /* initialise radacc radius for slope criterion */
  if (bFirst) {
    oldrad=calc_radius(&edi->vecs.radacc);
    bFirst=FALSE;
  }
    
  /* calculate full conservative forces & velocities when required */
  if (do_per_step(edstep,edi->outfrq) && bHave_force) {
    snew(vdum,ned);
    for(i=0; (i<ned); i++) {
      for(j=0; (j<DIM); j++) 
	vdum[i][j]=(x[i][j]-xold[i][j])*dt_1;
    }
    snew(fdum,ned);
    for(i=0; (i<ned); i++) {
      mas=top->atoms.atom[i].m;
      for(j=0; (j<DIM); j++) 
	fdum[i][j]=force[i][j]+(x[i][j]-x_unc[i][j])*dt_2*mas;
    }
  }

  /* fit the structure */
  snew(transvec,ned);
  rmsd=fitit(ned,x,edi,transvec,rotmat);

  /* update radsam references, when required */
  if (do_per_step(edstep,edi->maxedsteps) && edstep > 0) {
    project(x,edi,"x");
    projectx(edi,x,&edi->vecs.radacc);
    projectx(edi,x,&edi->vecs.radfix);
    oldrad=-1.e5;
  }

  /* update radacc references, when required */
  if (do_per_step(edstep,iupdate) && edstep > 0) {
    rad=calc_radius(&edi->vecs.radacc);
    if ((rad-oldrad) < edi->slope) {
      project(x,edi,"x");
      projectx(edi,x,&edi->vecs.radacc);
      oldrad=0.0;
    }
    else
      oldrad=rad;
  }

  /* apply the constraints */
  ed_cons(x,edi,edstep);

  /* produce output, when required */
  if (do_per_step(edstep,edi->outfrq) && bHave_force) {
    /* rotate forces and velocities */
    rotate_vec(ned,vdum,rotmat);
    rotate_vec(ned,fdum,rotmat);
    project(vdum,edi,"v");
    project(fdum,edi,"f");
    project(x,edi,"x");
    sfree(vdum);
    sfree(fdum);
    write_edo(edi,edstep,rmsd);
    fflush(edi->edo);
  }

  /* write to log, when required */
  if ((edstep > 0) && do_per_step(edstep,edi->logfrq)) {
    fprintf(log,"ED sampling information, step %d\n",edstep);
    project(x,edi,"x");
    write_edidx(log,edi);
    fprintf(log,"acceptance radius = %f\n",
	    calc_radius(&edi->vecs.radacc));
    fprintf(log,"fixed increment radius = %f\n",
	    calc_radius(&edi->vecs.radfix));
    fprintf(log,"contracting radius = %f\n",
	    calc_radius(&edi->vecs.radcon));
    fflush(log);
  }
  
  /* undo fit */
  rmfit(ned,x,transvec,rotmat);

  if (edstep == ir->nsteps) ffclose(edi->edo);
  sfree(transvec);
}

void rmfit(int ned,rvec *x,rvec *transvec,matrix rotmat)
{
  int i,j,k;
  matrix r_inv;
  rvec xdum;

  /* invert the rotation matrix and apply */
  m_inv(rotmat,r_inv);
  for(i=0;(i<ned);i++) {
    for(j=0;(j<3);j++)
      xdum[j]=x[i][j];
    for(j=0;(j<3);j++) {
      x[i][j]=0;
      for(k=0;(k<3);k++)
        x[i][j]+=r_inv[k][j]*xdum[k];
    }
  }
 
  /* subtract the translation vector */
  for(i=0;(i<ned);i++)
    rvec_dec(x[i],transvec[i]);
}

void rotate_vec(int nr,rvec *x,matrix rotmat)
{
  int i,j,k;
  rvec xdum;

  /* apply the rotation matrix */
  for(i=0;(i<nr);i++) {
    for(j=0;(j<3);j++)
      xdum[j]=x[i][j];
    for(j=0;(j<3);j++) {
      x[i][j]=0;
      for(k=0;(k<3);k++)
        x[i][j]+=rotmat[k][j]*xdum[k];
    }
  }
}

void ed_cons(rvec *x,t_edpar *edi,int step)
{
  int i;
  
  /* subtract the average positions */
  for (i=0;(i<edi->sav.nr);i++) 
      rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);
  
  /* apply the constraints */
  if (step >= 0) do_linfix(x,edi,step);
  do_linacc(x,edi);
  if (step >= 0) do_radfix(x,edi,step);
  do_radacc(x,edi);
  do_radcon(x,edi);
  
  /* add average positions */
  for (i=0;(i<edi->sav.nr);i++) 
    rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
}

void do_linfix(rvec *x,t_edpar *edi,int step)
{
  int i,j,k;
  real proj,add;

  /* loop over linfix vectors */
  for (i=0;(i < edi->vecs.linfix.neig);i++) {
    /* calculate the projection */
    proj=do_projectx(edi,x,edi->vecs.linfix.vec[i]);
    /*fprintf(stderr,"Proj[%d]=%f\n",edi->vecs.linfix.ieig[i],proj);*/
    /* calculate the correction */
    add=edi->vecs.linfix.refproj[i]+step*edi->vecs.linfix.stpsz[i]-proj;
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=add*edi->vecs.linfix.vec[i][j][k];
    }
  }
}

void do_linacc(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real proj,add;

  /* loop over linacc vectors */
  for (i=0;(i<edi->vecs.linacc.neig);i++) {
    /* calculate the projection */
    proj=do_projectx(edi,x,edi->vecs.linacc.vec[i]);
    /* calculate the correction */
    add=0.0;
    if (edi->vecs.linacc.stpsz[i] > 0.0) {
      if ((proj-edi->vecs.linacc.refproj[i]) < 0.0)
	add=edi->vecs.linacc.refproj[i]-proj;
    }
    if (edi->vecs.linacc.stpsz[i] < 0.0) {
      if ((proj-edi->vecs.linacc.refproj[i]) > 0.0)
	add=edi->vecs.linacc.refproj[i]-proj;
    }
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=add*edi->vecs.linacc.vec[i][j][k];
    }
    /* new positions will act as reference */
    edi->vecs.linacc.refproj[i]=proj+add;
  }
}

void do_radfix(rvec *x,t_edpar *edi,int step)
{
  int i,j,k;
  real *proj,rad=0.0,ratio;

  if (edi->vecs.radfix.neig == 0) return;
  snew(proj,edi->vecs.radfix.neig);
  /* loop over radfix vectors */
  for (i=0;(i<edi->vecs.radfix.neig);i++) {
    /* calculate the projections, radius */
    proj[i]=do_projectx(edi,x,edi->vecs.radfix.vec[i]);
    rad+=pow((proj[i]-edi->vecs.radfix.refproj[i]),2);
  }
  rad=sqrt(rad);
  ratio=(edi->vecs.radfix.stpsz[0]+edi->vecs.radfix.radius)/rad-1.0;
  edi->vecs.radfix.radius+=edi->vecs.radfix.stpsz[0];

  /* loop over radfix vectors */
  for (i=0;(i<edi->vecs.radfix.neig);i++) {
    proj[i]-=edi->vecs.radfix.refproj[i];
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=proj[i]*ratio*edi->vecs.radfix.vec[i][j][k];
    }
  }  
  sfree(proj);
}

void do_radacc(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real *proj,rad=0.0,ratio=0.0;

  if (edi->vecs.radacc.neig == 0) return;
  snew(proj,edi->vecs.radacc.neig);
  /* loop over radacc vectors */
  for (i=0;(i<edi->vecs.radacc.neig);i++) {
    /* calculate the projections, radius */
    proj[i]=do_projectx(edi,x,edi->vecs.radacc.vec[i]);
    /*fprintf(stderr,"Proj[%d]=%f\n",edi->vecs.radacc.ieig[i],proj[i]);*/
    rad+=pow((proj[i]-edi->vecs.radacc.refproj[i]),2);
  }
  rad=sqrt(rad);
  /*fprintf(stderr,"Radius=%f\n",rad);*/

  /* only correct when radius decreased */
  if (rad < edi->vecs.radacc.radius) {
    ratio=edi->vecs.radacc.radius/rad-1.0;
    rad=edi->vecs.radacc.radius;
  }
  else
    edi->vecs.radacc.radius=rad;

  /* loop over radacc vectors */
  for (i=0;(i<edi->vecs.radacc.neig);i++) {
    proj[i]-=edi->vecs.radacc.refproj[i];
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=proj[i]*ratio*edi->vecs.radacc.vec[i][j][k];
    }
  }  
  sfree(proj);
}

void do_radcon(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real *proj,rad=0.0,ratio=0.0;

  if (edi->vecs.radcon.neig == 0) return;
  snew(proj,edi->vecs.radcon.neig);
  /* loop over radcon vectors */
  for (i=0;(i<edi->vecs.radcon.neig);i++) {
    /* calculate the projections, radius */
    proj[i]=do_projectx(edi,x,edi->vecs.radcon.vec[i]);
    rad+=pow((proj[i]-edi->vecs.radcon.refproj[i]),2);
  }
  rad=sqrt(rad);

  /* only correct when radius increased */
  if (rad > edi->vecs.radcon.radius) {
    ratio=edi->vecs.radcon.radius/rad-1.0;
    rad=edi->vecs.radcon.radius;
  }
  else
    edi->vecs.radcon.radius=rad;

  /* loop over radcon vectors */
  for (i=0;(i<edi->vecs.radcon.neig);i++) {
    proj[i]-=edi->vecs.radcon.refproj[i];
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=proj[i]*ratio*edi->vecs.radcon.vec[i][j][k];
    }
  }  
  sfree(proj);
}

void write_edo(t_edpar *edi,int step,real rmsd)
{
  fprintf(edi->edo,"%d\n",step);
  fprintf(edi->edo,"%f\n",rmsd);
  write_proj(edi->edo,edi,"x");
  write_proj(edi->edo,edi,"v");
  write_proj(edi->edo,edi,"f");
}

void write_proj(FILE *out,t_edpar *edi,char *mode)
{
  do_write_proj(out,&edi->vecs.mon,mode);
  do_write_proj(out,&edi->vecs.linfix,mode);
  do_write_proj(out,&edi->vecs.linacc,mode);
  do_write_proj(out,&edi->vecs.radfix,mode);
  do_write_proj(out,&edi->vecs.radacc,mode);
  do_write_proj(out,&edi->vecs.radcon,mode);
}

void do_write_proj(FILE *out,t_eigvec *vec,char *mode)
{
  int i;
  
  for (i=0;(i<vec->neig);i++) {
    if (strcmp(mode,"x") == 0)
      fprintf(out,"%f ",vec->xproj[i]);
    else if (strcmp(mode,"v") == 0)
      fprintf(out,"%f ",vec->vproj[i]);
    else if (strcmp(mode,"f") == 0)
      fprintf(out,"%f ",vec->fproj[i]);
  }
  if (vec->neig > 0) fprintf(out,"\n");
}

void write_edidx(FILE *out,t_edpar *edi)
{
  int i;

  fprintf(out,"monitor eigenvectors");
  for (i=0;(i<edi->vecs.mon.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.mon.ieig[i],edi->vecs.mon.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"linfix  eigenvectors");
  for (i=0;(i<edi->vecs.linfix.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.linfix.ieig[i],edi->vecs.linfix.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"linacc  eigenvectors");
  for (i=0;(i<edi->vecs.linacc.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.linacc.ieig[i],edi->vecs.linacc.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radfix  eigenvectors");
  for (i=0;(i<edi->vecs.radfix.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radfix.ieig[i],edi->vecs.radfix.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radacc  eigenvectors");
  for (i=0;(i<edi->vecs.radacc.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radacc.ieig[i],edi->vecs.radacc.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radcon  eigenvectors");
  for (i=0;(i<edi->vecs.radcon.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radcon.ieig[i],edi->vecs.radcon.xproj[i]);
  fprintf(out,"\n");  
}

