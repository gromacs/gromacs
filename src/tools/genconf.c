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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_genconf_c = "$Id$";

#include "maths.h"
#include "macros.h"
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

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "genconf multiplies a given coordinate file by simply stacking them",
    "on top of each other, like a small child playing with wooden blocks.",
    "The program makes a grid of [IT]user defined[it]",
    "proportions ([TT]-nbox[tt]), ",
    "and interspaces the grid point with an extra space [TT]-dist[tt].[PAR]",
    "When option [TT]-rot[tt] is used the program does not check for overlap",
    "between molecules on grid points. It is recommended to make the box in",
    "the input file at least as big as the coordinates + Vander Waals radius."
  };
  /* The program should be more flexible, to allow for random displacement 
     off lattice points (for each cartesian coordinate), and specify the 
     (maximum) random rotation per coordinate to be useful for building 
     membranes.
     */

  int     vol;          
  t_atoms *atoms;       /* list with all atoms */
  char    title[STRLEN];
  rvec    *x,*v;        /* coordinates? */
  vec4    *xrot,*vrot;  
  matrix  box;          /* box length matrix */
  rvec    shift;         
  int     natoms;       /* number of atoms in one molecule  */
  int     nres;         /* number of molecules? */
  int     i,j,k,l,m,ndx,nrdx,nx,ny,nz;
  
  t_filenm fnm[] = {
    { efSTX, "-f", "conf", ffREAD  },
    { efSTO, "-o", "out",  ffWRITE }
  };
#define NFILE asize(fnm)
  static rvec nrbox={1,1,1};
  static int  seed;               /* seed for random number generator */
  static bool bRandom = FALSE;    /* False: no random rotations */
  static rvec dist={0,0,0};       /* space added between molecules ? */
  static rvec max_rot={90,90,90}; /* maximum rotation */
/*  static rvec max_tr={0,0,0}; */    /* maximum translation */
  t_pargs pa[] = {
    { "-nbox",  FALSE, etRVEC, &nrbox,  "Number of boxes" },
    { "-dist",  FALSE, etRVEC, &dist,   "Distance between boxes" },
    { "-seed",  FALSE, etINT,  &seed,   "Random generator seed" },
    { "-rot",   FALSE, etBOOL, &bRandom,"Randomly rotate conformations" },
    { "-maxrot",FALSE, etRVEC, &max_rot,"Maximum random rotation" },
    /*
    { "-maxtrans",FALSE, etRVEC, &max_tr),  "Max translation" },
    */
  };
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  nx=(int)(nrbox[XX]+0.5);
  ny=(int)(nrbox[YY]+0.5);
  nz=(int)(nrbox[ZZ]+0.5);
  
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    fatal_error(0,"Number of boxes (-nbox) should be positive");
  
  vol=nx*ny*nz;     /* calculate volume in grid points (= nr. molecules) */

  get_stx_coordnum(opt2fn("-f",NFILE,fnm),&natoms); 
  snew(atoms,1);
  init_t_atoms(atoms,natoms*vol,FALSE);
  snew(x,natoms*vol);              /* get space for coordinates of all atoms */
  snew(xrot,natoms);               /* get space for rotation matrix? */
  snew(v,natoms*vol);              /* velocities. not really needed? */ 
  snew(vrot,natoms); 
  read_stx_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,v,box);

  atoms->nr/=vol;
  nres=atoms->nres;                /* nr of residues in one element? */

  for(i=0; (i<nx); i++) {          /* loop over all gridpositions    */
    shift[XX]=i*(dist[XX]+box[XX][XX]);
    
    for(j=0; (j<ny); j++) {
      shift[YY]=j*(dist[YY]+box[YY][YY]);
      
      for(k=0; (k<nz); k++)  {
	shift[ZZ]=k*(dist[ZZ]+box[ZZ][ZZ]);
	
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

  box[XX][XX] = nx*(box[XX][XX]+dist[XX]); /* make box bigger */
  box[YY][YY] = ny*(box[YY][YY]+dist[YY]);
  box[ZZ][ZZ] = nz*(box[ZZ][ZZ]+dist[ZZ]);

  move_x(natoms*vol,x,box);          /* put atoms in box? */
    
  atoms->nr*=vol;
  atoms->nres*=vol;
  write_sto_conf(opt2fn("-o",NFILE,fnm),title,atoms,x,v,box);
  
  thanx(stdout);
  
  return 0;
}
