/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Gromacs Runs On Most of All Computer Systems
 */
#include "copyrite.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"

void get_params(char *gcin,char *gcout,
		int *nx,int *ny,int *nz,real *dx,real *dy,real *dz,
		int *seed,bool *bRandom,rvec max_rot,rvec max_tr)
{
  t_inpfile *inp;
  int       ninp;

  inp=read_inpfile(gcin,&ninp);

  ITYPE("nx",                 *nx,         1);
  ITYPE("ny",                 *ny,         1);
  ITYPE("nz",                 *nz,         1);
  RTYPE("dx",                 *dx,      0.15);
  RTYPE("dy",                 *dy,      0.15);
  RTYPE("dz",                 *dz,      0.15);
  ITYPE("seed",               *seed,    1993);
  ETYPE("randomize",          *bRandom,   yesno_names);
  RTYPE("max_x_rotation",     max_rot[0], 0.0);
  RTYPE("max_y_rotation",     max_rot[1], 0.0);
  RTYPE("max_z_rotation",     max_rot[2], 0.0);
  RTYPE("max_x_translation",  max_tr[0],  0.0);
  RTYPE("max_y_translation",  max_tr[1],  0.0);
  RTYPE("max_z_translation",  max_tr[2],  0.0);

  write_inpfile(gcout,ninp,inp);
}

static void rand_rot(int natoms,rvec x[],rvec v[],vec4 xrot[],vec4 vrot[],
                     int *seed,rvec max_rot)
{
  mat4 mt1,mt2,mr[DIM],mtemp1,mtemp2,mtemp3,mxtot,mvtot;
  rvec xcm;
  real phi;
  int  i,m;
  
  clear_rvec(xcm);
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++) {
      xcm[m]+=x[i][m]/natoms;   /* get center of mass of one molecule  */
    }
  fprintf(stderr,"center of mass: %f, %f, %f\n",xcm[0],xcm[1],xcm[2]);
  
  translate(-xcm[XX],-xcm[YY],-xcm[ZZ],mt1);  /* move c.o.ma to origin */
  for(m=0; (m<DIM); m++) {
    phi=M_PI*max_rot[m]*rando(seed)/180;
    rotate(m,phi,mr[m]);
  }
  translate(xcm[XX],xcm[YY],xcm[ZZ],mt2);

  mult_matrix(mtemp1,mt2,mr[XX]);
  mult_matrix(mtemp2,mr[YY],mr[ZZ]);
  mult_matrix(mtemp3,mtemp1,mtemp2);
  mult_matrix(mxtot,mtemp3,mt1);
  mult_matrix(mvtot,mr[XX],mtemp2);
  
  for(i=0; (i<natoms); i++) {
    m4_op(mxtot,x[i],xrot[i]);
    m4_op(mvtot,v[i],vrot[i]);
  }
}

static void move_x(int natoms,rvec x[],matrix box)
{
  int  i,m;
  rvec xcm;

  clear_rvec(xcm);
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++)
      xcm[m]+=x[i][m];
  for(m=0; (m<DIM); m++)
    xcm[m] = 0.5*box[m][m]-xcm[m]/natoms;
  for(i=0; (i<natoms); i++)
    for(m=0; (m<DIM); m++)
      x[i][m] += xcm[m];
}

void main(int argc, char *argv[])
{
  static char *desc[] = {
    "genconf multiplies a given coordinate file by simply stacking them",
    "on top of each other, like a small child playing with wooden blocks.",
    "The program makes a grid of [IT]user defined[it]",
    "proportions (nx, ny, nz in the input file), ",
    "and interspaces the grid point with an extra space dx,dy and dz."
  };
  static char *bugs[] = {
    "When option -r is used (randomize) the program does not check for overlap between molecules on grid points. It is recommended to make the box in the input file at least as big as the coordinates + Vander Waals radius. ",
    "The program should be more flexible, to allow for random displacement off lattice points (for each cartesian coordinate), and specify the (maximum) random rotation per coordinate to be useful for building membranes."
  };

  int     nx,ny,nz,vol; /* nx,ny,nz number of gridpoints in x,y,z direction */
  t_atoms *atoms;       /* list with all atoms */
  char    title[STRLEN];
  real    dx,dy,dz;     /* space added between molecules ? */
  rvec    *x,*v;        /* coordinates? */
  vec4    *xrot,*vrot;  
  matrix  box;          /* box length matrix */
  bool    bRandom;      /* False: no random rotations */
  rvec    shift;         
  rvec    max_rot;      /* maximum rotation */
  rvec    max_tr,tr;    /* maximum translation, actual translation */
  int     seed;         /* seed for random number generator */
  int     natoms;       /* number of atoms in one molecule  */
  int     nres;         /* number of molecules? */
  int     i,j,k,l,m,ndx,nrdx;
  t_filenm fnm[] = {
    { efGCP, "-f",  NULL,     ffREAD  },
    { efGCP, "-po", "gcout",  ffWRITE },
    { efGRO, "-ci", "confin", ffREAD  },
    { efGRO, "-co", "confout",ffWRITE }
  };
#define NFILE asize(fnm)
    
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,asize(bugs),bugs);

  get_params(opt2fn("-f",NFILE,fnm),opt2fn("-po",NFILE,fnm),
	     &nx,&ny,&nz,&dx,&dy,&dz,&seed,&bRandom,max_rot,max_tr);

  if ((nx <= 0) || (ny <= 0) || (nz <= 0)) {
    fprintf(stderr,"(nx <= 0) || (ny <= 0) || (nz <= 0)\n");
    exit(1);
  }
  vol=nx*ny*nz;     /* calculate volume in grid points (= nr. molecules) */

  get_coordnum(opt2fn("-ci",NFILE,fnm),&natoms); 
  snew(atoms,1); 
  snew(atoms->resname,natoms*vol); /* get space for _all_ atoms */
  snew(atoms->atomname,natoms*vol);
  snew(atoms->atom,natoms*vol);
  snew(x,natoms*vol);              /* get space for coordinates of all atoms */
  snew(xrot,natoms);               /* get space for rotation matrix? */
  snew(v,natoms*vol);              /* velocities. not really needed? */ 
  snew(vrot,natoms); 
  read_whole_conf(opt2fn("-ci",NFILE,fnm),title,atoms,x,v,box);
 
  nres=atoms->nres;                /* nr of residues in one element? */

  for(i=0; (i<nx); i++) {          /* loop over all gridpositions    */
    shift[XX]=i*(dx+box[XX][XX]);
    
    for(j=0; (j<ny); j++) {
      shift[YY]=j*(dy+box[YY][YY]);
      
      for(k=0; (k<nz); k++)  {
	shift[ZZ]=k*(dz+box[ZZ][ZZ]);
	
	for(l=0; (l<DIM); l++)  /* to get negative values too */
	    tr[l] = (rando(&seed)-0.5)*2*max_tr[l];
 
	ndx=(i*ny*nz+j*nz+k)*natoms;
	nrdx=(i*ny*nz+j*nz+k)*nres;
	
	if (ndx > 0) {

	  /* Random rotation on input coords */
	  if (bRandom)
	   rand_rot(natoms,x,v,xrot,vrot,&seed,max_rot);
	  
	  for(l=0; (l<natoms); l++) {
	    for(m=0; (m<DIM); m++) {
	      if (bRandom) {
		x[ndx+l][m]=xrot[l][m]+shift[m];
		v[ndx+l][m]=vrot[l][m];
	      }
	      else {
		x[ndx+l][m]=x[l][m]+shift[m];
		v[ndx+l][m]=v[l][m];
	      }
	    }
	    atoms->atom[ndx+l].resnr=nrdx+atoms->atom[l].resnr;
	    atoms->atomname[ndx+l]=atoms->atomname[l];
	  }

	  for(l=0; (l<nres); l++)
	    atoms->resname[nrdx+l]=atoms->resname[l];
	}
      }
    }
  }

  box[XX][XX] = nx*(box[XX][XX]+dx); /* make box bigger */
  box[YY][YY] = ny*(box[YY][YY]+dy);
  box[ZZ][ZZ] = nz*(box[ZZ][ZZ]+dz);

  move_x(natoms*vol,x,box);          /* put atoms in box? */
    
  atoms->nr*=vol;
  atoms->nres*=vol;
  write_conf(opt2fn("-co",NFILE,fnm),title,atoms,x,v,box);
  
  thanx(stdout);
}






