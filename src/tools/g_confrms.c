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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_g_confrms_c = "$Id$";

/* 
 * fit coordinates of file 1 to coordinates of file 2
 * fit is done on the atoms as indicated by the NDX- file 
 */

#include "filenm.h"
#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "confio.h"
#include "txtdump.h"
#include "gstat.h"

#define EPS  1.0e-09

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,d;
  
  /* check periodic boundary */
  for(d=0;(d<DIM);d++) {
    for(n=1;(n<atoms->nr);n++) {
      dist = x[n][d]-x[n-1][d];
      if ( fabs(dist) > 0.9 * box[d][d]  ) {
	if ( dist >  0 )
	  x[n][d]-=box[d][d];
	else
	  x[n][d]+=box[d][d];
      } 	
    }
  }
}

real my_fit(int natoms1,atom_id *index1,atom_id *index2,
	    rvec *xp,
	    int natoms2,rvec *x)
{
  int    c,r,n,j,m,i,irot;
  double omega[7][7],om[7][7],d[7],xnr,xpc;
  matrix vh,vk,R,vh_d,vk_d,u;
  real   du,sigd,mn;
  int    index;
  real   max_d;
  rvec   x_old;

  for(i=0;(i<7);i++) {
    d[i]=0;
    for(j=0;(j<7);j++) {
      omega[i][j]=0;
      om[i][j]=0;
    }
  }


#ifdef OLDCODE
  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms1);n++) {
    if ((mn = w_rls[n]) != 0.0) {
      for(c=0; (c<DIM); c++) {
	xpc=xp[n][c];
	for(r=0; (r<DIM); r++) {
	  xnr=x[n][r];
	  u[c][r]+=mn*xnr*xpc;
	}
      }
    }
  }
#else
  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms1);n++) {
    mn = 1.0;
    
    for(c=0; (c<DIM); c++) {
      xpc=xp[index1[n]][c];
      for(r=0; (r<DIM); r++) {
	xnr=x[index2[n]][r];
	u[c][r]+=mn*xnr*xpc;
      }
    }
  }

#endif
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
  /*real   **omega = input matrix a[1..n][1..n] must be symmetric
   *int     natoms = number of rows and columns
   *real      NULL = d[1]..d[n] are the eigenvalues of a[][]
   *real       **v = v[1..n][1..n] contains the vectors in columns
   *int      *irot = number of jacobi rotations
   */

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
  
  /*rotate X*/
  for(j=0;(j<natoms2);j++) {
    for(m=0;(m<3);m++)
      x_old[m]=x[j][m];
    for(r=0;(r<3);r++) {
      x[j][r]=0;
      for(c=0;(c<3);c++)
        x[j][r]+=R[c][r]*x_old[c];
    }
  }

  /* calculate deviation of fit */
  {
    int  aid1,aid2,i,m;
    real e=0,tmas=0;
    
    /*calculate energy */
    for(i=0; (i<natoms1); i++) {
      aid1=index1[i];
      aid2=index2[i];
      for(m=0;(m<DIM);m++)
      e+=(x[aid2][m]-xp[aid1][m])*(x[aid2][m]-xp[aid1][m]);
    }
    
    /*return energy*/
    return (sqrt(e/natoms1));
  }
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_confrms computes the root mean square deviation (RMSD) of two",
    "structure files after LSQ fitting the structures on top of each other.",
    "The two structures do NOT need to have the same number of atoms,",
    "only the two index groups used for the fit need to be identical",
    "[PAR]",
    "Option -o2 produces a Gromos file with the second structure fitted on",
    "the first structure,",
    "[PAR]",
    "Option -op produces a PDB file with the superimposed structures"
  };
  static bool bRmpbc=FALSE;
  
  t_pargs pa[] = {
    { "-pbc", FALSE, etBOOL, &bRmpbc, "Remove periodic boundary conditions" }
  };
  t_filenm fnm[] = {
    { efSTX, "-f1",  "conf1", ffREAD  },
    { efSTX, "-f2",  "conf2", ffREAD  },
    { efGRO, "-o2", "confout2.gro", ffOPTWR },
    { efPDB, "-op",  "conf.pdb",  ffOPTWR },
    { efNDX, "-n1" ,  "fit1.ndx",   ffOPTRD  },
    { efNDX, "-n2" ,  "fit2.ndx",   ffOPTRD  }
  };

#define NFILE asize(fnm)
  
  /* the two gromos files */
  char    title_1[STRLEN],title_2[STRLEN],*name1,*name2;
  t_atoms atoms_1,atoms_2;
  int     natoms_1,natoms_2; 
  rvec    *x_1,*v_1,*x_2,*v_2,*xf_1,*xf_2;
  matrix  box_1,box_2;
  
  /* counters */
  int   i,j,m;
  
  /* center of mass calculation */
  real tmas_1,tmas_2;
  rvec xcm_1,xcm_2;

  /* variables for fit */
  char *groupnames_1,*groupnames_2;
  int isize_1,isize_2;
  atom_id *index_1,*index_2;
  real rms;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  /* reading first gromos file nb this is going 
   * to be the reference structure*/
  get_stx_coordnum(opt2fn("-f1",NFILE,fnm),&(atoms_1.nr));
  snew(x_1,atoms_1.nr);
  snew(v_1,atoms_1.nr);
  snew(atoms_1.resname,atoms_1.nr);
  snew(atoms_1.atom,atoms_1.nr);
  snew(atoms_1.atomname,atoms_1.nr);
  fprintf(stderr,"\nReading first structure file\n");
  read_stx_conf(opt2fn("-f1",NFILE,fnm),title_1,&atoms_1,x_1,v_1,box_1);
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title_1,atoms_1.nr,atoms_1.nres);
  srenew(atoms_1.resname,atoms_1.nres);
    
  if ( bRmpbc ) 
    rm_gropbc(&atoms_1,x_1,box_1);


  /* check if we have an NDX file */
  /*
  if ( opt2bSet("-n1",NFILE,fnm) ) {
    fprintf(stderr,"\nSelect group for root least square fit\n");
    rd_index(opt2fn("-n1",NFILE,fnm),1,&isize_1,&index_1,&groupnames_1);
  } else {
    isize_1 = atoms_1.nr;
    snew(index_1,isize_1);
    for(i=0;(i<isize_1);i++)
      index_1[i]=i;
    groupnames_1 = title_1;
  }
  */
  get_index(&atoms_1,opt2fn_null("-n1",NFILE,fnm),
	    1,&isize_1,&index_1,&groupnames_1);
 
  if (isize_1 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");


  /* reading second gromos file */
  get_stx_coordnum(opt2fn("-f2",NFILE,fnm),&(atoms_2.nr));
  snew(x_2,atoms_2.nr);
  snew(v_2,atoms_2.nr);
  snew(atoms_2.resname,atoms_2.nr);
  snew(atoms_2.atom,atoms_2.nr);
  snew(atoms_2.atomname,atoms_2.nr);
  fprintf(stderr,"\nReading second structure file\n");
  read_stx_conf(opt2fn("-f2",NFILE,fnm),title_2,&atoms_2,x_2,v_2,box_2);
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title_2,atoms_2.nr,atoms_2.nres);
  srenew(atoms_2.resname,atoms_2.nres);

  if ( bRmpbc ) 
    rm_gropbc(&atoms_2,x_2,box_2);


  
  /* check if we have an NDX file */
  /*
  if ( opt2bSet("-n2",NFILE,fnm) ) {
    fprintf(stderr,"\nSelect group for root least square fit\n");
    rd_index(opt2fn("-n2",NFILE,fnm),1,&isize_2,&index_2,&groupnames_2);
  } else {
    isize_2 = atoms_2.nr;
    snew(index_2,isize_2);
    for(i=0;(i<isize_2);i++)
      index_2[i]=i;
    groupnames_2 = title_2;
  }
  */
  get_index(&atoms_2,opt2fn_null("-n2",NFILE,fnm),
	    1,&isize_2,&index_2,&groupnames_2);

  /* check isizes, must be equal */
  if ( isize_2 != isize_1 )
    fatal_error(0,"isize_2 != isize_2\n");

  if (isize_2 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");

  for(i=0; (i<isize_1); i++) {
    name1=*atoms_1.atomname[index_1[i]];
    name2=*atoms_2.atomname[index_2[i]];
    if (strcmp(name1,name2))
      fprintf(stderr,"Warning: atomnames at index %d don't match: %d %s, %d %s\n",i+1,index_1[i]+1,name1,index_2[i]+1,name2);
  }

  /* calculate center of mass of reference structure */
  for(m=0;(m<3);m++)
    xcm_1[m]=0;
  for(i=0;(i<isize_1);i++)
    for(m=0;(m<DIM);m++) 
      xcm_1[m]+=x_1[index_1[i]][m];
  for(m=0;(m<3);m++)
    xcm_1[m]/=isize_1;
  for(i=0;(i<atoms_1.nr);i++)
    for(m=0;(m<DIM);m++)
      x_1[i][m]-=xcm_1[m];
  
  /* calculate center of mass of structure to be fitted */
  for(m=0;(m<3);m++)
    xcm_2[m]=0;
  for(i=0;(i<isize_2);i++)
    for(m=0;(m<DIM);m++) 
      xcm_2[m]+=x_2[index_2[i]][m];
  for(m=0;(m<3);m++)
    xcm_2[m]/=isize_2;
  for(i=0;(i<atoms_2.nr);i++)
    for(m=0;(m<DIM);m++)
      x_2[i][m]-=xcm_2[m];
  
  /*do the least squares fit to the reference structure*/
  rms = my_fit(isize_1,index_1,index_2,x_1,atoms_2.nr,x_2);
  
  fprintf(stderr,"RMS = %8.3f\n",rms);

  /* reset coordinates of reference and fitted structure */
  for(i=0;(i<atoms_1.nr);i++) {
    for(m=0;(m<3);m++)
      x_1[i][m]+=xcm_1[m];
  }
  for(i=0;(i<atoms_2.nr);i++) {
    for(m=0;(m<3);m++)
      x_2[i][m]+=xcm_1[m];
  }

  /* calculate the rms deviation */
  
  

  /* write gromos file of fitted structure */
  if ( opt2bSet("-o2",NFILE,fnm) )
    write_conf(opt2fn("-o2",NFILE,fnm),title_2,&atoms_2,x_2,v_2,box_2);

  /* write pdb conf to output file */
  if ( opt2bSet("-op",NFILE,fnm) ) {
    t_atoms **temp_at;
    rvec    **temp_x;
    snew(temp_at,2);
    temp_at[0]=&atoms_1;
    temp_at[1]=&atoms_2;
    snew(temp_x,2);
    temp_x[0]=x_1;
    temp_x[1]=x_2;
    
    write_pdb_confs(opt2fn("-op",NFILE,fnm),temp_at,temp_x,2);
    sfree(temp_at);
    sfree(temp_x);
  }
  thanx(stdout);
  
  return 0;
}









