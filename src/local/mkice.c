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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_mkice_c = "$Id$";

#include <stdio.h>
#include <math.h>
#include "typedefs.h"
#include "statutil.h"
#include "copyrite.h"
#include "fatal.h"
#include "pdbio.h"
#include "macros.h"
#include "smalloc.h"
#include "vec.h"
#include "pbc.h"
#include "physics.h"
#include "names.h"
#include "txtdump.h"
#include "trnio.h"
#include "symtab.h"
#include "confio.h"

#define ODIST 0.274
#define HDIST 0.1
#define TET   109.47
#define DCONS 0.117265878

static char *watname[] = { "OW1", "HW2", "HW3", "DW", "SW" };
static real qspc[]     = { -0.82, 0.41, 0.41 };
static real qyaw[]     = { 1.24588, 0.62134, 0.62134, 0.0, -2.48856 };
static real spc_lj[3][6] = {
  { 2.6171e-3, 2.6331e-6, 0, 0, 0, 0 },
  { 0,      0,      0, 0, 0, 0 },
  { 0,      0,      0, 0, 0, 0 }
};
#define CHH 2e-9
#define CHS 2e-9
#define COS 2e-6
static real yaw_lj[5][10] = {
  { 0, 0,   0, 0,   0,   0, 0, 0, 0,      COS },
  { 0, 0,   0, CHH, 0, CHH, 0, 0, 0,      CHS },
  { 0, 0,   0, CHH, 0, CHH, 0, 0, 0,      CHS },
  { 0, 0,   0, 0,   0,   0, 0, 0, 0,      0   },
  { 0, COS, 0, CHS, 0, CHS, 0, 0, 2.6e-3, 0   }
};

void unitcell(rvec x[],rvec box,bool bYaw)
{
#define cx  0.81649658
#define cy  0.47140452
#define cy2 0.94280904
#define cz  0.33333333
  
  rvec xx[24] = {
    { 0,   0,         0 }, /* O1 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx, cy,       -cz },
    { cx, cy,       -cz }, /* O2 */
    { 0, 0,       -1    }, /* H relative to Oxygen */
    { cx,-cy,       +cz },
    { cx, cy+cy2,     0 }, /* O3 */
    { -cx, cy,    -cz   }, /* H relative to Oxygen */
    { 0,   -cy2,    -cz },
    { 0,  2*cy+cy2, -cz }, /* O4 */
    {-cx,-cy,       +cz }, /* H relative to Oxygen */
    { 0 , cy2,      +cz },
    { 0,   0,         1 }, /* O5 */
    {-cx, cy,       +cz }, /* H relative to Oxygen */
    { 0 , -cy2,     +cz },
    { cx, cy,      1+cz }, /* O6 */
    { -cx, -cy,     -cz }, /* H relative to Oxygen */
    { 0,   cy2,     -cz },
    { cx, cy+cy2,     1 }, /* O7 */
    { 0,  0,       -1   }, /* H relative to Oxygen */
    { cx, cy,       +cz },
    { 0,  2*cy+cy2,1+cz }, /* O8 */
    { 0,   0,         1 }, /* H relative to Oxygen */
    { cx,   -cy,    -cz }
  };
  int  i,iin,iout,j,m;
  rvec tmp,t2,dip;
  
  clear_rvec(dip);
  for(i=0; (i<8); i++) {
    iin = 3*i;
    if (bYaw)
      iout = 5*i;
    else
      iout = iin;
    svmul(ODIST,xx[iin],x[iout]);
    svmul(-0.82,x[iout],t2);
    rvec_inc(dip,t2);
    for(j=1; (j<=2); j++) {
      svmul(HDIST,xx[iin+j],tmp);
      rvec_add(x[iout],tmp,x[iout+j]);
      svmul(0.41,x[iout+j],t2);
      rvec_inc(dip,t2);
    }
    if (bYaw)
      for(m=0; (m<DIM); m++) 
	x[iout+3][m] = x[iout+4][m] = 
	  (1-2*DCONS)*x[iout][m]+DCONS*(x[iout+1][m]+x[iout+2][m]);
  }
  
  box[XX] = 2*cx;
  box[YY] = 2*(cy2+cy);
  box[ZZ] = 2*(1+cz);
  for(i=0; (i<DIM); i++)
    box[i] *= ODIST;
    
  printf("Unitcell:  %10.5f  %10.5f  %10.5f\n",box[XX],box[YY],box[ZZ]);
  printf("Dipole:    %10.5f  %10.5f  %10.5f (e nm)\n",dip[XX],dip[YY],dip[ZZ]);
}

static real calc_ener(real c6,real c12,rvec dx,tensor vir)
{
  real r,e,f;
  int  m,n;
  
  r  = norm(dx);
  e  = c12*pow(r,-12.0) - c6*pow(r,-6.0);
  f  = 12*c12*pow(r,-14.0) - 6*c6*pow(r,-8.0);
  for(m=0; (m<DIM); m++) 
    for(n=0; (n<DIM); n++)
      vir[m][n] -= 0.5*f*dx[m]*dx[n];
      
  return e;
}

void interactions(FILE *fp,rvec x[],matrix box,real rcut,bool bYaw)
{
  /* List of h bonds in the crystal unitcell */
  typedef struct {
    int d,h,a;
  } t_hb;
  t_hb hb[16] = {
    {  0,  1, 12 },
    {  0,  2,  3 },
    {  3,  4, 15 },
    {  3,  5,  0 },
    {  6,  7,  9 },
    {  6,  8,  3 },
    {  9, 10,  6 },
    {  9, 11,  0 },
    { 12, 13, 15 },
    { 12, 14, 21 },
    { 15, 16, 12 },
    { 15, 17, 18 },
    { 18, 19,  6 },
    { 18, 20, 21 },
    { 21, 22,  9 },
    { 21, 23, 18 }
  };
  int    i,ad,ah,aa;
  rvec   dx[16];
  real   elj;
  tensor vir;
  
  init_pbc(box,FALSE);
  /* Initiate the periodic boundary conditions. Set bTruncOct to
   * TRUE when using a truncated octahedron box.
   */

  clear_mat(vir);
  elj = 0;
  for(i=0; (i<asize(hb)); i++) {
    ad = hb[i].d;
    ah = hb[i].h;
    aa = hb[i].a;
    pbc_dx(x[ah],x[aa],dx[i]);
    elj += calc_ener(0.0,1.0e-8,dx[i],vir);
  }
  pr_rvecs(fp,0,"H-A dist",dx,16);
  pr_rvecs(fp,0,"Virial",vir,3);
  fprintf(fp,"LJ energy is %g\n",elj);
}

void virial(FILE *fp,bool bFull,int nmol,rvec x[],matrix box,real rcut,
	    bool bYaw,real q[],bool bLJ)
{
  int  i,j,im,jm,natmol,ik,jk,m,ninter;
  rvec dx,f,ftot,dvir,vir,pres,xcmi,xcmj,*force;
  real dx6,dx2,dx1,fscal,c6,c12,vcoul,v12,v6,vctot,v12tot,v6tot;
  
  init_pbc(box,FALSE);
  /* Initiate the periodic boundary conditions. Set bTruncOct to
   * TRUE when using a truncated octahedron box.
   */
  fprintf(fp,"%3s   -  %3s: %6s %6s %6s  %6s %8s %8s %8s\n",
	  "ai","aj","dx","dy","dz","|d|","virx","viry","virz");
  clear_rvec(ftot);
  clear_rvec(vir);
  ninter = 0;
  vctot  = 0;
  v12tot = 0;
  v6tot  = 0;
  natmol = bYaw ? 5 : 3;
  snew(force,nmol*natmol);
  
  for(i=0; (i<nmol); i++) {
    im = natmol*i;
    /* Center of geometry */
    clear_rvec(xcmi);
    for(ik=0; (ik<natmol); ik++)
      rvec_inc(xcmi,x[im+ik]);
    for(m=0; (m<DIM); m++)
      xcmi[m] /= natmol;

    for(j=i+1; (j<nmol); j++) {
      jm = natmol*j;
      /* Center of geometry */
      clear_rvec(xcmj);
      for(jk=0; (jk<natmol); jk++)
	rvec_inc(xcmj,x[jm+jk]);
      for(m=0; (m<DIM); m++)
	xcmj[m] /= natmol;

      /* First check COM-COM distance */
      pbc_dx(xcmi,xcmj,dx);
      if (norm(dx) < rcut) {
	ninter++;
	/* Neirest neighbour molecules! */
	for(ik=0; (ik<natmol); ik++) 
	  for(jk=0; (jk<natmol); jk++) {
	    pbc_dx(x[im+ik],x[jm+jk],dx);
	    dx2    = iprod(dx,dx);
	    dx1    = sqrt(dx2);
	    vcoul  = q[ik]*q[jk]*ONE_4PI_EPS0/dx1;
	    vctot += vcoul;
	    
	    if (bLJ) {
	      if (bYaw) {
		c6  = yaw_lj[ik][2*jk];
		c12 = yaw_lj[ik][2*jk+1];
	      }
	      else {
		c6  = spc_lj[ik][2*jk];
		c12 = spc_lj[ik][2*jk+1];
	      }
	      dx6    = dx2*dx2*dx2;
	      v6     = c6/dx6;
	      v12    = c12/(dx6*dx6);
	      v6tot -= v6;
	      v12tot+= v12;
	      fscal  = (vcoul+12*v12-6*v6)/dx2;
	    }
	    else
	      fscal  = vcoul/dx2;
	    for(m=0; (m<DIM); m++) {
	      f[m]     = dx[m]*fscal;
	      dvir[m]  = -0.5*dx[m]*f[m];
	      vir[m]  += dvir[m];
	    }
	    rvec_inc(force[ik+im],f);
	    rvec_dec(force[jk+jm],f);
	    if (bFull)
	      fprintf(fp,"%3s%4d-%3s%4d: %6.3f %6.3f %6.3f %6.3f"
		      " %8.3f %8.3f %8.3f\n",
		      watname[ik],im+ik,watname[jk],jm+jk,
		      dx[XX],dx[YY],dx[ZZ],norm(dx),
		      dvir[XX],dvir[YY],dvir[ZZ]);
	  }
      }
    }
  }
  fprintf(fp,"There were %d interactions between the %d molecules\n",
	  ninter,nmol);
  fprintf(fp,"Vcoul: %10.4e  V12: %10.4e  V6: %10.4e  Vtot: %10.4e\n",
	  vctot,v12tot,v6tot,vctot+v12tot+v6tot);
  pr_rvec(fp,0,"vir ",vir,DIM);
  
  for(m=0; (m<DIM); m++) 
    pres[m] = -2*PRESFAC/(det(box))*vir[m];
  pr_rvec(fp,0,"pres",pres,DIM);
  pr_rvecs(fp,0,"force",force,natmol*nmol);
  sfree(force);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "mkice generates an ice crystal in the Ih crystal form which is the",
    "most stable form. The rectangular unitcell contains eight molecules",
    "and all oxygens are tetrahedrally coordinated."
  };
  static int nx=1,ny=1,nz=1;
  static bool bYaw=FALSE,bLJ=TRUE,bFull=TRUE,bSeries=FALSE;
  static real rcut=0.3;
  t_pargs pa[] = {
    { "-nx", FALSE, etINT, {&nx}, "nx" },
    { "-ny", FALSE, etINT, {&ny}, "ny" },
    { "-nz", FALSE, etINT, {&nz}, "nz" },
    { "-yaw",FALSE, etBOOL,{&bYaw},"Generate dummies and shell positions" },
    { "-lj", FALSE, etBOOL,{&bLJ},"Use LJ as well as coulomb for virial calculation" },
    { "-rcut",FALSE,etREAL,{&rcut},"Cut-off for virial calculations" },
    { "-full",FALSE,etBOOL,{&bFull},"Full virial output" },
    { "-series",FALSE, etBOOL, {&bSeries}, "Do a series of virial calculations with different cut-off (from 0.3 up till the specified one)" }
  };
  t_filenm fnm[] = {
    { efPDB, "-p", "ice", ffWRITE },
    { efSTO, "-c", NULL,  ffOPTRD },
    { efTRN, "-o", "ice", ffOPTWR }
  };
#define NFILE asize(fnm)

  FILE      *fp;
  int       i,j,k,n,nmax,m,natom,natmol;
  t_atoms   *pdba;
  t_atoms   atoms;
  t_symtab  symtab;
  rvec      box,tmp,*xx;
  matrix    boxje;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,asize(desc),
		    desc,0,NULL);
  if (debug) {
    fprintf(debug,"nx  = %3d, ny  = %3d,  nz   = %3d\n",nx,ny,nz);
    fprintf(debug,"YAW = %3s, LJ  = %3s,  rcut = %g\n",yesno_names[bYaw],
	    yesno_names[bLJ],rcut);
  }

  natmol = bYaw ? 5 : 3;
  natom = natmol*8;
  nmax = natom*nx*ny*nz;
  snew(xx,nmax);
  snew(pdba,1);
  init_t_atoms(pdba,nmax,TRUE);
  pdba->nr = nmax;
  open_symtab(&symtab);
  for(n=0; (n<nmax); n++) {
    pdba->pdbinfo[n].type   = epdbATOM;
    pdba->pdbinfo[n].atomnr = n;
    pdba->atom[n].resnr     = n/natmol;
    pdba->atomname[n] = put_symtab(&symtab,watname[(n % natmol)]);
    pdba->resname[n] = put_symtab(&symtab,"SOL");

    sprintf(pdba->pdbinfo[n].pdbresnr,"%d",n);
    pdba->atom[n].chain = ' ';
  }
  
  /* Generate the unit cell */
  unitcell(xx,box,bYaw);
  if (debug) {
    clear_mat(boxje);
    boxje[XX][XX] = box[XX];
    boxje[YY][YY] = box[YY];
    boxje[ZZ][ZZ] = box[ZZ];
    /*interactions(debug,xx,boxje,0.9,bYaw);*/
  }
  n=0;
  for(i=0; (i<nx); i++) {
    tmp[XX] = i*box[XX];
    for(j=0; (j<ny); j++) {
      tmp[YY] = j*box[YY];
      for(k=0; (k<nz); k++) {
	tmp[ZZ] = k*box[ZZ];
	for(m=0; (m<natom); m++,n++) {
	  rvec_add(xx[n % natom],tmp,xx[n]);
	}
      }
    }
  }
    
  clear_mat(boxje);
  boxje[XX][XX] = box[XX]*nx;
  boxje[YY][YY] = box[YY]*ny;
  boxje[ZZ][ZZ] = box[ZZ]*nz;
  
  printf("Crystal:   %10.5f  %10.5f  %10.5f\n",
	 nx*box[XX],ny*box[YY],nz*box[ZZ]);
  
  if (debug) {
    if (bSeries)
      for(i=3; (i<=10*rcut); i++) {
	fprintf(debug,"This is with rcut = %g\n",i*0.1);
	virial(debug,bFull,nmax/natmol,xx,boxje,
	       0.1*i,bYaw,bYaw ? qyaw : qspc,bLJ);
      }    
    else
      virial(debug,bFull,nmax/natmol,xx,boxje,
	     rcut,bYaw,bYaw ? qyaw : qspc,bLJ);
  }
  write_sto_conf(ftp2fn(efPDB,NFILE,fnm),
		 "Ice Ih crystal generated by mkice",pdba,xx,NULL,boxje);
		 
  if (ftp2bSet(efTRN,NFILE,fnm))
    write_trn(ftp2fn(efTRN,NFILE,fnm),0,0,0,boxje,nmax,xx,NULL,NULL);
	       
  return 0;
}

