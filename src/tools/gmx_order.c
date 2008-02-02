/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <ctype.h>

#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "tpxio.h"
#include "cmat.h"

/****************************************************************************/
/* This program calculates the order parameter per atom for an interface or */
/* bilayer, averaged over time.                                             */
/* S = 1/2 * (3 * cos(i)cos(j) - delta(ij))                                 */
/* It is assumed that the order parameter with respect to a box-axis        */
/* is calculated. In that case i = j = axis, and delta(ij) = 1.             */
/*                                                                          */
/* Peter Tieleman,  April 1995                                              */
/* P.J. van Maaren, November 2005     Added tetrahedral stuff               */
/****************************************************************************/

static void find_nearest_neighbours(t_topology top,int natoms, matrix box,
				    rvec x[],int maxidx,atom_id index[], 
				    real time,
				    real *sgmean, real *skmean)
{
  FILE    *fpoutdist;
  char    fnsgdist[32];
  int     ix,jx,nsgbin, *sgbin;
  int     i1,i2,i,ibin,j,k,l,n,*nn[4];
  rvec    dx,dx1,dx2,rj,rk,urk,urj;
  real    cost,cost2,*sgmol,*skmol,rmean,rmean2,r2,box2,*r_nn[4];
  t_pbc   pbc;
  t_mat   *dmat;
  t_dist  *d;
  int     m1,mm;
  int    **nnb;
  real   onethird=1.0/3.0;
  
  /*  dmat = init_mat(maxidx, FALSE); */

  box2 = box[XX][XX] * box[XX][XX];
  
  for (i=0; (i<4); i++) {
    snew(r_nn[i],natoms);
    snew(nn[i],natoms);
    
    for (j=0; (j<natoms); j++) {
      r_nn[i][j] = box2;
    }
  }
  
  snew(sgmol,maxidx);
  snew(skmol,maxidx);

  /* Must init pbc every step because of pressure coupling */
  set_pbc(&pbc,box);
  rm_pbc(&(top.idef),natoms,box,x,x);

  nsgbin = 1 + 1/0.0005;
  snew(sgbin,nsgbin);

  *sgmean = 0.0;
  *skmean = 0.0;
  l=0;
  for (i=0; (i<maxidx); i++) { /* loop over index file */
    ix = index[i];
    for (j=0; (j<maxidx); j++) {
      if (i == j) continue;

      jx = index[j];
      
      pbc_dx(&pbc,x[ix],x[jx],dx);
      r2=iprod(dx,dx);

      /* set_mat_entry(dmat,i,j,r2); */

      /* determine the nearest neighbours */
      if (r2 < r_nn[0][i]) {
	r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
	r_nn[2][i] = r_nn[1][i]; nn[2][i] = nn[1][i];
	r_nn[1][i] = r_nn[0][i]; nn[1][i] = nn[0][i];
	r_nn[0][i] = r2;         nn[0][i] = j; 
      } else if (r2 < r_nn[1][i]) {
	r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
	r_nn[2][i] = r_nn[1][i]; nn[2][i] = nn[1][i];
	r_nn[1][i] = r2;         nn[1][i] = j;
      } else if (r2 < r_nn[2][i]) {
	r_nn[3][i] = r_nn[2][i]; nn[3][i] = nn[2][i];
	r_nn[2][i] = r2;         nn[2][i] = j;
      } else if (r2 < r_nn[3][i]) {
	r_nn[3][i] = r2;         nn[3][i] = j;
      }
    }


    /* calculate mean distance between nearest neighbours */
    rmean = 0;
    for (j=0; (j<4); j++) {
      r_nn[j][i] = sqrt(r_nn[j][i]);
      rmean += r_nn[j][i];
    }
    rmean /= 4;
    
    n = 0;
    sgmol[i] = 0.0;
    skmol[i] = 0.0;

    /* Chau1998a eqn 3 */
    /* angular part tetrahedrality order parameter per atom */
    for (j=0; (j<3); j++) {
      for (k=j+1; (k<4); k++) {
	pbc_dx(&pbc,x[ix],x[index[nn[k][i]]],rk);
	pbc_dx(&pbc,x[ix],x[index[nn[j][i]]],rj);

	unitv(rk,urk);
	unitv(rj,urj);
	
	cost = iprod(urk,urj) + onethird;
	cost2 = cost * cost;

	/* sgmol[i] += 3*cost2/32; 
	   sgmol[i] += cost2; */

	/* determine distribution */
	ibin = nsgbin * cost2;
	if (ibin < nsgbin)
	  sgbin[ibin]++;
	/* printf("%d %d %f %d %d\n", j, k, cost * cost, ibin, sgbin[ibin]);*/
	l++;
	n++;
      }
    }

    /* normalize sgmol between 0.0 and 1.0 */
    sgmol[i] = 3*sgmol[i]/32;
    *sgmean += sgmol[i];

    /* distance part tetrahedrality order parameter per atom */
    rmean2 = 4 * 3 * rmean * rmean;
    for (j=0; (j<4); j++) {
      skmol[i] += (rmean - r_nn[j][i]) * (rmean - r_nn[j][i]) / rmean2;
      /*      printf("%d %f (%f %f %f %f) \n",
	      i, skmol[i], rmean, rmean2, r_nn[j][i], (rmean - r_nn[j][i]) );
      */
    }
    
    *skmean += skmol[i];
  } /* loop over entries in index file */
  
  *sgmean /= maxidx;
  *skmean /= maxidx;

  sfree(sgbin);
  sfree(sgmol);
  sfree(skmol);
  for (i=0; (i<4); i++) {
    sfree(r_nn[i]);
    sfree(nn[i]);
  }
}


static void calc_tetra_order_parm(char *fnNDX,char *fnTPS,char *fnTRX,
			    char *sgfn,char *skfn)
{
  FILE       *fpsg=NULL,*fpsk=NULL;
  t_topology top;
  char       title[STRLEN],fn[STRLEN],subtitle[STRLEN];
  int        status;
  int        natoms;
  real       t;
  rvec       *xtop,*x;
  matrix     box;
  real       sg,sk;
  atom_id    **index;
  char       **grpname;
  int        i,*isize,ng,nframes;
  
  read_tps_conf(fnTPS,title,&top,&xtop,NULL,box,FALSE);

  ng = 1;
  /* get index groups */
  printf("Select the group that contains the atoms you want to use for the tetrahedrality order parameter calculation:\n");
  snew(grpname,ng);
  snew(index,ng);
  snew(isize,ng);
  get_index(&top.atoms,fnNDX,ng,isize,index,grpname);

  /* Analyze trajectory */
  natoms=read_first_x(&status,fnTRX,&t,&x,box);
  if ( natoms > top.atoms.nr )
    gmx_fatal(FARGS,"Topology (%d atoms) does not match trajectory (%d atoms)",
	      top.atoms.nr,natoms);
  check_index(NULL,ng,index[0],NULL,natoms);

  fpsg=xvgropen(sgfn,"S\\sg\\N Angle Order Parameter","Time (ps)","S\\sg\\N");
  fpsk=xvgropen(skfn,"S\\sk\\N Distance Order Parameter","Time (ps)","S\\sk\\N");

  /* loop over frames */
  nframes = 0;
  do {
    find_nearest_neighbours(top,natoms,box,x,isize[0],index[0],t,
			    &sg,&sk);
    fprintf(fpsg,"%f %f\n", t, sg);
    fprintf(fpsk,"%f %f\n", t, sk);
    nframes++;
  } while (read_next_x(status,&t,natoms,x,box));
  close_trj(status);
 
  sfree(grpname);
  sfree(index);
  sfree(isize);

  fclose(fpsg);
  fclose(fpsk);
}


/* Print name of first atom in all groups in index file */
static void print_types(atom_id index[], atom_id a[], int ngrps, 
			char *groups[], t_topology *top)
{
  int i;

  fprintf(stderr,"Using following groups: \n");
  for(i = 0; i < ngrps; i++)
    fprintf(stderr,"Groupname: %s First atomname: %s First atomnr %u\n", 
	    groups[i], *(top->atoms.atomname[a[index[i]]]), a[index[i]]);
  fprintf(stderr,"\n");
}

static void check_length(real length, int a, int b)
{
  if (length > 0.3)
    fprintf(stderr,"WARNING: distance between atoms %d and "
	    "%d > 0.3 nm (%f). Index file might be corrupt.\n", 
	    a, b, length);
}

void calc_order(char *fn, atom_id *index, atom_id *a, rvec **order,
		real ***slOrder, real *slWidth, int nslices, bool bSliced, 
		bool bUnsat, t_topology *top, int ngrps, int axis)
{ 
  rvec *x0,          /* coordinates with pbc                           */
    *x1,             /* coordinates without pbc                        */
    dist;            /* vector between two atoms                       */
  matrix box;        /* box (3x3)                                      */
  int   status;  
  rvec  cossum,      /* sum of vector angles for three axes            */
    Sx, Sy, Sz,      /* the three molecular axes                       */
    tmp1, tmp2,      /* temp. rvecs for calculating dot products       */
    frameorder;      /* order parameters for one frame                 */
  real *slFrameorder; /* order parameter for one frame, per slice      */
  real length,       /* total distance between two atoms               */
    t,               /* time from trajectory                           */
    z_ave,z1,z2;     /* average z, used to det. which slice atom is in */
  int natoms,        /* nr. atoms in trj                               */
    nr_tails,        /* nr tails, to check if index file is correct    */
    size=0,          /* nr. of atoms in group. same as nr_tails        */  
    i,j,m,k,l,teller = 0,
    slice,           /* current slice number                           */
    nr_frames = 0,
    *slCount;        /* nr. of atoms in one slice                      */
   real dbangle = 0, /* angle between double bond and  axis            */ 
        sdbangle = 0;/* sum of these angles                            */

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0) 
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");

  snew(slCount, nslices);
  snew(*slOrder, nslices);
  for(i = 0; i < nslices; i++)
    snew((*slOrder)[i],ngrps);
  snew(*order,ngrps);
  snew(slFrameorder, nslices);
  snew(x1, natoms);
  
  if (bSliced) {
    *slWidth = box[axis][axis]/nslices;
    fprintf(stderr,"Box divided in %d slices. Initial width of slice: %f\n",
	    nslices, *slWidth);
  } 

  nr_tails = index[1] - index[0];
  fprintf(stderr,"Number of elements in first group: %d\n",nr_tails);
  /* take first group as standard. Not rocksolid, but might catch error in index*/

  teller = 0; 

  /*********** Start processing trajectory ***********/
  do {
    if (bSliced)
      *slWidth = box[axis][axis]/nslices;
    teller++;
    
    rm_pbc(&(top->idef),top->atoms.nr,box,x0,x1);

    /* Now loop over all groups. There are ngrps groups, the order parameter can
       be calculated for grp 1 to grp ngrps - 1. For each group, loop over all 
       atoms in group, which is index[i] to (index[i+1] - 1) See block.h. Of 
       course, in this case index[i+1] -index[i] has to be the same for all 
       groups, namely the number of tails. i just runs over all atoms in a tail,
       so for DPPC ngrps = 16 and i runs from 1 to 14, including 14
     */
    
    for (i = 1; i < ngrps - 1; i++) {
      clear_rvec(frameorder);
      
      size = index[i+1] - index[i];
      if (size != nr_tails)
	gmx_fatal(FARGS,"grp %d does not have same number of"
		" elements as grp 1\n",i); 
      
      for (j = 0; j < size; j++) {
	if (bUnsat) {
	  /* Using convention for unsaturated carbons */
	  /* first get Sz, the vector from Cn to Cn+1 */
	  rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i]+j]], dist); 
	  length = norm(dist);
	  check_length(length, a[index[i]+j], a[index[i+1]+j]);
	  svmul(1/length, dist, Sz);

	  /* this is actually the cosine of the angle between the double bond
	     and axis, because Sz is normalized and the two other components of
	     the axis on the bilayer are zero */
	  sdbangle += acos(Sz[axis]);  
	} else {
	  /* get vector dist(Cn-1,Cn+1) for tail atoms */
	  rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i-1]+j]], dist);
	  length = norm(dist);      /* determine distance between two atoms */
	  check_length(length, a[index[i-1]+j], a[index[i+1]+j]);
	  
	  svmul(1/length, dist, Sz);
	  /* Sz is now the molecular axis Sz, normalized and all that */
	}

	/* now get Sx. Sx is normal to the plane of Cn-1, Cn and Cn+1 so
	   we can use the outer product of Cn-1->Cn and Cn+1->Cn, I hope */
	rvec_sub(x1[a[index[i+1]+j]], x1[a[index[i]+j]], tmp1);
	rvec_sub(x1[a[index[i-1]+j]], x1[a[index[i]+j]], tmp2);
	cprod(tmp1, tmp2, Sx);
	svmul(1/norm(Sx), Sx, Sx);
	
	/* now we can get Sy from the outer product of Sx and Sz   */
	cprod(Sz, Sx, Sy);
	svmul(1/norm(Sy), Sy, Sy);

	/* the square of cosine of the angle between dist and the axis.
	   Using the innerproduct, but two of the three elements are zero
	   Determine the sum of the orderparameter of all atoms in group 
	   */
	cossum[XX] = sqr(Sx[axis]); /* this is allowed, since Sa is normalized */
	cossum[YY] = sqr(Sy[axis]);
	cossum[ZZ] = sqr(Sz[axis]);

	for (m = 0; m < DIM; m++)
          frameorder[m] += 0.5 * (3 * cossum[m] - 1);
	
	if (bSliced) {
	  /* get average coordinate in box length for slicing,
	     determine which slice atom is in, increase count for that
	     slice. slFrameorder and slOrder are reals, not
	     rvecs. Only the component [axis] of the order tensor is
	     kept, until I find it necessary to know the others too 
	   */
	  
	  z1 = x1[a[index[i-1]+j]][axis]; 
	  z2 = x1[a[index[i+1]+j]][axis];
	  z_ave = 0.5 * (z1 + z2);
	  if (z_ave < 0)
	    z_ave += box[axis][axis];
	  if (z_ave > box[axis][axis])
	    z_ave -= box[axis][axis];

	  slice  = (int)(0.5 + (z_ave / (*slWidth))) - 1;
          slCount[slice]++;               /* determine slice, increase count */

	  slFrameorder[slice] += 0.5 * (3 * cossum[axis] - 1);
	}
      }   /* end loop j, over all atoms in group */
      
      for (m = 0; m < DIM; m++)
	(*order)[i][m] += (frameorder[m]/size);
      
      for (k = 0; k < nslices; k++) {
	if (slCount[k]) {     /* if no elements, nothing has to be added */
	  (*slOrder)[k][i] += slFrameorder[k]/slCount[k];
	  slFrameorder[k] = 0; slCount[k] = 0;
	}
      }   /* end loop i, over all groups in indexfile */
    }
    nr_frames++;
    
  } while (read_next_x(status,&t,natoms,x0,box));
  /*********** done with status file **********/
  
  fprintf(stderr,"\nRead trajectory. Printing parameters to file\n");
  
  /* average over frames */
  for (i = 1; i < ngrps - 1; i++) {
    svmul(1.0/nr_frames, (*order)[i], (*order)[i]);
    fprintf(stderr,"Atom %d Tensor: x=%g , y=%g, z=%g\n",i,(*order)[i][XX],
	    (*order)[i][YY], (*order)[i][ZZ]);
    if (bSliced) {
      for (k = 0; k < nslices; k++)
	(*slOrder)[k][i] /= nr_frames;
    }
  }

  if (bUnsat)
    fprintf(stderr,"Average angle between double bond and normal: %f\n", 
	    180*sdbangle/(nr_frames * size*M_PI));

  sfree(x0);  /* free memory used by coordinate arrays */
  sfree(x1);
}


void order_plot(rvec order[], real *slOrder[], char *afile, char *bfile, 
		char *cfile, int ngrps, int nslices, real slWidth, bool bSzonly)
{
  FILE       *ord, *slOrd;           /* xvgr files with order parameters  */
  int        atom, slice;            /* atom corresponding to order para.*/
  char       buf[256];               /* for xvgr title */
  real      S;                      /* order parameter averaged over all atoms */

  if (bSzonly) {
    sprintf(buf,"Orderparameters Sz per atom");
    ord = xvgropen(afile,buf,"Atom","S");
    fprintf(stderr,"ngrps = %d, nslices = %d",ngrps, nslices);

    sprintf(buf, "Orderparameters per atom per slice");
    slOrd = xvgropen(bfile, buf, "Slice", "S");
    
    for (atom = 1; atom < ngrps - 1; atom++)
      fprintf(ord,"%12d       %12g\n", atom, order[atom][ZZ]);

    for (slice = 0; slice < nslices; slice++) {
      S = 0;
      for (atom = 1; atom < ngrps - 1; atom++)
	S += slOrder[slice][atom];
      fprintf(slOrd,"%12g     %12g\n", slice*slWidth, S/atom);
    }

  } else {
    sprintf(buf,"Order tensor diagonal elements");
    ord = xvgropen(afile,buf,"Atom","S");
    sprintf(buf,"Deuterium order parameters");
    slOrd = xvgropen(cfile,buf, "Atom", "Scd");

    for (atom = 1; atom < ngrps - 1; atom++) {
      fprintf(ord,"%12d   %12g   %12g   %12g\n", atom, order[atom][XX],
	      order[atom][YY], order[atom][ZZ]);
      fprintf(slOrd,"%12d   %12g\n", atom, -1 * (0.6667 * order[atom][XX] + 
						 0.333 * order[atom][YY]));
    }
    
    fclose(ord);
    fclose(slOrd);
  }
}

int gmx_order(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute the order parameter per atom for carbon tails. For atom i the",
    "vector i-1, i+1 is used together with an axis. The index file has to contain",
    "a group with all equivalent atoms in all tails for each atom the",
    "order parameter has to be calculated for. The program can also give all",
    "diagonal elements of the order tensor and even calculate the deuterium",
    "order parameter Scd (default). If the option -szonly is given, only one",
    "order tensor component (specified by the -d option) is given and the",
    "order parameter per slice is calculated as well. If -szonly is not",
    "selected, all diagonal elements and the deuterium order parameter is",
    "given.[PAR]"
    "The tetrahedrality order parameters can be determined",
    "around an atom. Both angle an distance order parameters are calculated. See",
    "P.-L. Chau and A.J. Hardwick, Mol. Phys., 93, (1998), 511-518.",
    "for more details.[BR]",
    ""
  };

  static int  nslices = 1;                    /* nr of slices defined       */
  static bool bSzonly = FALSE;                /* True if only Sz is wanted  */
  static bool bUnsat = FALSE;                 /* True if carbons are unsat. */
  static char *normal_axis[] = { NULL, "z", "x", "y", NULL };
  t_pargs pa[] = {
    { "-d",      FALSE, etENUM, {normal_axis}, 
      "Direction of the normal on the membrane" },
    { "-sl",     FALSE, etINT, {&nslices},
      "Calculate order parameter as function of boxlength, dividing the box"
      " in #nr slices." },
    { "-szonly", FALSE, etBOOL,{&bSzonly},
      "Only give Sz element of order tensor. (axis can be specified with -d)" },
    { "-unsat",  FALSE, etBOOL,{&bUnsat},
      "Calculate order parameters for unsaturated carbons. Note that this can"
      "not be mixed with normal order parameters." }
  };

  rvec      *order;                         /* order par. for each atom   */
  real      **slOrder;                      /* same, per slice            */
  real      slWidth = 0.0;                  /* width of a slice           */
  char      **grpname;            	    /* groupnames                 */
  int       ngrps,                          /* nr. of groups              */
            i,
            axis=0;                         /* normal axis                */
  t_topology *top;                	    /* topology 		  */ 
  atom_id   *index,             	    /* indices for a              */
            *a;                             /* atom numbers in each group */
  t_blocka  *block;                         /* data from index file       */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efTPX, NULL, NULL,  ffREAD },    	    /* topology file           	  */
    { efXVG,"-o","order", ffWRITE }, 	    /* xvgr output file 	  */
    { efXVG,"-od","deuter", ffWRITE },      /* xvgr output file           */
    { efXVG,"-os","sliced", ffWRITE },      /* xvgr output file           */
    { efXVG,"-Sg","sg-ang", ffWRITE },      /* xvgr output file           */
    { efXVG,"-Sk","sk-dist", ffWRITE },     /* xvgr output file           */
  };
  bool      bSliced = FALSE;                /* True if box is sliced      */
#define NFILE asize(fnm)
  char *sgfnm,*skfnm,*ndxfnm,*tpsfnm,*trxfnm;

  CopyRight(stderr,argv[0]);
  
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0, NULL);

  sgfnm = opt2fn_null("-Sg",NFILE,fnm);
  skfnm = opt2fn_null("-Sk",NFILE,fnm);
  ndxfnm = ftp2fn(efNDX,NFILE,fnm);
  tpsfnm = ftp2fn(efTPX,NFILE,fnm);
  trxfnm = ftp2fn(efTRX,NFILE,fnm);
  
  /* tetraheder order parameter */
  if (skfnm || sgfnm) {
    calc_tetra_order_parm(ndxfnm,tpsfnm,trxfnm,sgfnm,skfnm);
    /* view xvgr files */
    do_view(opt2fn("-Sg",NFILE,fnm), NULL);
    do_view(opt2fn("-Sk",NFILE,fnm), NULL);
  } 
  else {  
    /* tail order parameter */
    /* Calculate axis */
    if (strcmp(normal_axis[0],"x") == 0) axis = XX;
    else if (strcmp(normal_axis[0],"y") == 0) axis = YY;
    else if (strcmp(normal_axis[0],"z") == 0) axis = ZZ;
    else gmx_fatal(FARGS,"Invalid axis, use x, y or z");
    
    switch (axis) {
    case 0:
      fprintf(stderr,"Taking x axis as normal to the membrane\n");
      break;
    case 1:
      fprintf(stderr,"Taking y axis as normal to the membrane\n");
      break;
    case 2:
      fprintf(stderr,"Taking z axis as normal to the membrane\n");
      break;
    }
    
    if (nslices > 1) {
      bSliced = TRUE;
      fprintf(stderr,"Dividing box in %d slices.\n\n", nslices);
    }
    
    if (bSzonly)
      fprintf(stderr,"Only calculating Sz\n");
    if (bUnsat)
      fprintf(stderr,"Taking carbons as unsaturated!\n");
    
    top = read_top(ftp2fn(efTPX,NFILE,fnm));     /* read topology file */
    
    block = init_index(ftp2fn(efNDX,NFILE,fnm),&grpname);
    index = block->index;                       /* get indices from t_block block */
    a = block->a;                               /* see block.h                    */
    ngrps = block->nr;           
    
    /* show atomtypes, to check if index file is correct */
    print_types(index, a, ngrps, grpname, top);
    
    calc_order(ftp2fn(efTRX,NFILE,fnm), index, a, &order, 
	       &slOrder, &slWidth, nslices, bSliced, bUnsat,
	       top, ngrps, axis); 
    
    order_plot(order, slOrder, opt2fn("-o",NFILE,fnm), opt2fn("-os",NFILE,fnm), 
	       opt2fn("-od",NFILE,fnm), ngrps, nslices, slWidth, bSzonly);
    
    do_view(opt2fn("-o",NFILE,fnm), NULL);      /* view xvgr file */
    do_view(opt2fn("-os",NFILE,fnm), NULL);     /* view xvgr file */
    do_view(opt2fn("-od",NFILE,fnm), NULL);     /* view xvgr file */
  }
  
  thanx(stderr);
  
  return 0;
}
