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
 * Good ROcking Metal Altar for Chronical Sinners
 */
#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "rmpbc.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"

/****************************************************************************/
/* This program calculates the diffusion coefficients of water per slice, as*/
/* function of the z-coordinate, from the MSD displacements using the       */
/* Einstein relation                                                        */
/*                                                                          */
/* Peter Tieleman, January 1996                                             */
/****************************************************************************/

#define TRESH 0.2    /* percentage of msd graph to neglect for fitting */
#define NRESTART 5  /* restart every 5 times per 5 ps */
#define LOWEST 20    /* lowert nr of waters in slice that stilll counts */

static int axis = 2; /* use Z-axis as default. */
static matrix box;   /* box (3x3)                               */


/* lsq things stolen from g_msd.cc. */
typedef struct {
  double yx,xx,sx,sy;
  double sigw_y;     /* sum of the weight factors */
  int    np;
} t_lsq;

void init_lsq(t_lsq *lsq)
{
  lsq->yx=lsq->xx=lsq->sx=lsq->sy=0.0;
  lsq->sigw_y=0.0;
  lsq->np=0;
}

void add_lsq(t_lsq *lsq,real x,real y,real w_y)
{
#ifdef MYDEBUG
  fprintf(bug,"Adding x,y,w_y: %f %f %f\n",x,y,w_y);
#endif
  lsq->yx+=w_y*y*x;
  lsq->xx+=x*x*w_y;
  lsq->sx+=x*w_y;
  lsq->sy+=w_y*y;
  lsq->np++;
  lsq->sigw_y+=w_y;
}

void get_lsq_ab(t_lsq *lsq,real *a,real *b)
{
  real yx,xx,sx,sy;
 
  yx=lsq->yx/lsq->sigw_y;
  xx=lsq->xx/lsq->sigw_y;
  sx=lsq->sx/lsq->sigw_y;
  sy=lsq->sy/lsq->sigw_y;

#ifdef MYDEBUG
  fprintf(stderr,"sigw_y: %f\n",lsq->sigw_y);
  fprintf(bug,"average yx,xx,sx,sy:%f %f %f %f\n",yx,xx,sx,sy);
#endif
  (*a)=(yx-sx*sy)/(xx-sx*sx);
  (*b)=(sy)-(*a)*(sx);
}

int find_slice(rvec cm, int nslices)
{
  int i, slice;
  rvec tmp;
  
  for (i = 0; i < DIM; i++)
    tmp[i] = cm[i];

  while (tmp[axis] < 0)
    tmp[axis] += box[axis][axis];
  while (tmp[axis] > box[axis][axis])
    tmp[axis] -= box[axis][axis];
  
  slice = (int)(tmp[axis] * nslices / box[axis][axis] ); 
  if (slice == nslices)
  {
#ifdef MYDEBUG
    fprintf(stderr,"Putting slice with max. number in 0\n");
#endif
    return 0;
  }
  if (slice < 0 || slice >= nslices)
  {
    fprintf(stderr,"Coordinate: %f\n", cm[axis]);
    fprintf(stderr,"HELP PANIC! slice = %d, OUT OF RANGE!\n",slice);
    slice = -1;
  }
  return slice;
}

void make_slicetab(rvec x0[], atom_id *index, int ngx, int nslices, 
		   int *slCount, int *slTab)
{
  int i; 
  int slice;

  for (i = 0; i < nslices; i++)
    slCount[i] = 0;

  for (i = 0; i < ngx; i++) /* put all atoms in index file in slice */
  {
    slice = find_slice(x0[index[i]], nslices);
    if (slice != -1)
    {
      slTab[i] = slice;
      slCount[slice]++;
    }
    else
      fprintf(stderr,"Warning: find_slice returned -1\n");
  }

#ifdef MYDEBUG
  for (i = 0; i < nslices; i++)
    fprintf(stderr,"%d molecules in slice %d\n", slCount[i], i); 
#endif
}

void calc_diffusion(char *fn, atom_id *index, int ngx, real interval, 
		    int nslices, rvec **diff, rvec **stdev)
{
  rvec *x0,          /* coordinate frame                        */
    *xprev, hbox,    /* previous frame                          */
    **ref;           /* reference frame                         */
  int   status;
  real w[NRESTART] ; /* weight factors for lsq averages         */
  rvec dr,r2,        /* displacement, squared displacement      */
      **slGframe,    /* per slic, diffusion in three directions */
      **slSigGframe, /* sum of slGframe squared                 */
      *slD,*slD2;    /* average dif. co.per slice, square of that*/
  real t, tin[NRESTART]; /* time from trajectory, interval time    */
  int natoms,        /* nr. atoms in trj                        */
    slice, i,ii,m,teller = 0,
     **slTab,        /* slTab[restart][atom] is the slice atom is in*/
     **slCount,      /* slCount[restart][slice] is the nr. of atoms in slice */
    *slDcount;  
  t_lsq ***slTot;   /* sucks. slTot[restart][slice][dimensin] */
  real a, b;
  real weight;

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0) 
  {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }

  fprintf(stderr,"Box divided in %d slices along axis %c. "
	  "Init. width of slice: %f\n", nslices, axis + 'X', 
	  box[axis][axis]/nslices);
  
  for (i = 0; i < NRESTART; i++)
  {
    tin[i] = t+i*NRESTART/interval; /* set initial restart times */
  }

  /* initialize tons of memory */
  snew(xprev, natoms);
  snew(ref, NRESTART);
  snew(slGframe, NRESTART);
  snew(slCount, NRESTART);
  snew(slTot, NRESTART);
  snew(slTab, NRESTART);
  snew(slSigGframe, NRESTART);
  for (i = 0; i < NRESTART; i++)
  {
    snew(ref[i], natoms);
    snew(slCount[i], nslices);
    snew(slTot[i], nslices);
    snew(slTab[i], natoms);
    snew(slGframe[i], nslices);
    snew(slSigGframe[i], nslices);
    for (m = 0; m < nslices; m++)
      snew(slTot[i][m], DIM);
  }
  snew(slD, nslices);
  snew(slD2, nslices);
  snew(slDcount, nslices);
  snew(*diff, nslices); 
  snew(*stdev, nslices);

  for (ii = 0; ii < NRESTART; ii++)
  {
    memcpy(ref[ii],x0,natoms*sizeof(rvec)); /* this is the reference frame  */
    for (i = 0; i < nslices; i++)
    {
      for (m = 0; m < DIM; m++)
	init_lsq(&(slTot[ii][i][m]));
    }
  }

  memcpy(xprev,x0,natoms*sizeof(rvec)); /* save this is the previous frame */

  teller = 0; 

  /* make initial slicetabels. This is a problem, because the table
     for restart 5 will be updated only after 10 ps. After the first
     time they'll be updated every 5 ps 
   */
  for (ii = 0; ii < NRESTART; ii++)
    make_slicetab(x0, index, ngx, nslices, slCount[ii], slTab[ii]);

  /*********** Start processing trajectory ***********/
  do 
  {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rFrame: %d",teller-1); 
    
    for (i = 0; i < nslices; i++)
      for (ii = 0; ii < NRESTART; ii++)
	for (m = 0; m < DIM; m++)
	{
	  slGframe[ii][i][m] = 0;
	  slSigGframe[ii][i][m] = 0;
	}
    
    for(m=0; (m<DIM); m++) 
      hbox[m]=0.5*box[m][m];
    
    /* prepare data so pbc is gone */
    for (i = 0; i < natoms; i++)  
    {
      for(m=0; (m<DIM); m++) 
      {
	while(x0[i][m]-xprev[i][m] <= hbox[m])
	  x0[i][m] += box[m][m];
	while(x0[i][m]-xprev[i][m] >  hbox[m])
	  x0[i][m] -= box[m][m];
      }
    }

    /* Loop over all atoms in the group (ngx), so typically all oxygen
       atoms. Calculate de displacement of each oxygen with respect to the
       previous frame, square it, and add it to the total msd for the 
       slice the atom belongs in. The index for the slice atom i is in 
       is just slTab[i], so slGframe[slTab[i]] is the msd for the atoms 
       in slice slTab[i], for one frame 
    */
    for (ii = 0; ii < NRESTART; ii++)
    {
      for (i = 0; i < ngx; i++) 
      {
	for (m = 0; m < DIM; m++)
	{
	  dr[m] = x0[index[i]][m] - ref[ii][index[i]][m];
	  r2[m] = dr[m]*dr[m];
	  slGframe[ii][slTab[ii][i]][m] += r2[m];
	  slSigGframe[ii][slTab[ii][i]][m] += r2[m]*r2[m];
	}
	/* ok, this sucks too. slGframe[ii] is the array with the r2
	 *for all _slices_, in all three dimensions, for restart point
	 *ii. slTab[ii][i] is the slice atom i is in in restart point
	 *ii, and [m] is the dimension we're looking at.  
	 */
	
      }
#ifdef MYDEBUG
      fprintf(stderr,"Calculated displacements for frame %d, sp. %d\n",
	      teller, ii);
#endif
      /* now we have calculated the msd for all atoms in this frame.
       Check if for certain starting points we have reach 'interval',
       and add the diffusion coeffecient to a bin (slD[i])
       */
      if ((t-tin[ii]) > TRESH*interval)  
      {
	for (i = 0; i < nslices; i++)
	  for (m = 0; m < DIM; m++)
	  {
	    if (slCount[ii][i])
	    {
	      slSigGframe[ii][i][m] /= slCount[ii][i];
	      slGframe[ii][i][m] /= slCount[ii][i];
	      weight = slSigGframe[ii][i][m] - 
		slGframe[ii][i][m]*slGframe[ii][i][m];
	      add_lsq(&(slTot[ii][i][m]), t-tin[ii], 
		      slGframe[ii][i][m],1/weight); 
#ifdef MYDEBUG	   
	    fprintf(stderr,"%f\t",1/weight);
#endif
	    }
	  }
      }

      if ( (t - tin[ii]) > interval )
      {
	tin[ii] = t; 
	for (i = 0; i < nslices; i++)
	  if (slCount[ii][i])
	  {
	    /* get diffusion coeffecient, and zero the lsq for next time */
	    for (m = 0; m < DIM; m++)
	    {
	      get_lsq_ab(&(slTot[ii][i][m]),&a,&b);
	      slD[i][m] += a*500; 
	      slD2[i][m] += a*a*500*500;  /* keep square */
	      init_lsq(&(slTot[ii][i][m]));
	    }
	    slDcount[i]++;
	  }
	/* now make new slTab and new ref. frame for this ii */
	memcpy(ref[ii],x0,natoms*sizeof(rvec));
	make_slicetab(x0, index, ngx, nslices, slCount[ii], slTab[ii]);
#ifdef MYDEBUG
	fprintf(stderr,"Making sltab for restart %d, time %f, frame%d\n\n",
		ii,t,teller);
	fprintf(stderr,"tin[%d] = %f\t",ii,tin[ii]);
#endif
      }
    }   /* end loop over NRESTARTS */
    memcpy(xprev,x0,natoms*sizeof(rvec)); /* current frame becomes previous */
    
  } while (read_next_x(status,&t,natoms,x0,box));
  /*********** done with status file **********/
  
  for (i = 0; i < nslices; i++)
    for (m = 0; m < DIM; m++)
      {
	slD[i][m] /= slDcount[i];
	slD2[i][m] /= slDcount[i];
	slD2[i][m] = (slD2[i][m] - slD[i][m]*slD[i][m]); /* get rms error */ 
     }
  
  *stdev = slD2;
  *diff = slD;
  close_trj(status);
}

void diffusion_plot(rvec diffusion[], rvec variance[], char *afile, 
		    int nslices)
{
  FILE       *dif;                   /* xvgr files with order parameters  */
  int        slice;                  /* atom corresponding to order para.*/
  char       buf[256];               /* for xvgr title */

  sprintf(buf,"Diffusion per slice");
  dif = xvgropen(afile,buf,"z(nm)","D");
  
  for (slice = 0; slice < nslices; slice++)
    fprintf(dif,"%12g     %12g %12g %12g %12g %12g %12g\n", 
	    slice*box[axis][axis]/nslices, 
	    diffusion[slice][XX], variance[slice][XX],
	    diffusion[slice][YY], variance[slice][YY],
	    diffusion[slice][ZZ], variance[slice][ZZ]);
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Simple program to calculate the diffusion coefficients of water"
    "as function of the z-coordinate"
  };

  static char *opts[] = {
    "-t",
    "-sl",
    "-d",
  };
  static char *odesc[] = {
    "Length of a time-interval",
    "Nr. of slices",
    "-d X|Y|Z, axis",
  };
 
  t_manual man = {asize(desc),desc,asize(opts),opts,odesc,NULL,NULL};
  char      **grpname;            	    /* groupnames                 */
  real      interval = 5.0;                 /* interval to use            */
  rvec      *diff, *variance;
  int       i, nslices = 1, ngx;
  atom_id   **index;             	    /* indices for a              */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL, ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL, ffREAD },    	    /* index file 		  */
    { efXVG,"-o",  NULL, ffWRITE },  	    /* xvgr output file 	  */
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE,
		    fnm, TRUE, &man);

  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-t") == 0) 
    {
      if (i < argc - 1)
      {
       sscanf(argv[i+1],"%f",&interval);
       i++;
      }
    }
    if (strcmp(argv[i],"-sl") == 0)
    {
      if (i < argc -1)
      {
	sscanf(argv[i+1],"%d",&nslices);
	i++;
      }
    }
    if (strcmp(argv[i],"-d") == 0)
    {
      if (i < argc -1)
      {
	axis = (int)argv[i+1][0] - 'X';
	i++;
      }
    }
  }

  fprintf(stderr,"Dividing box in %d slices.\n\n", nslices);
  fprintf(stderr,"Taking %c as axis\n",'X'+axis);

  snew(grpname,1);
  snew(index,1);
  rd_index(ftp2fn(efNDX,NFILE,fnm),1,&ngx,index,grpname);    /* read index */

  calc_diffusion(ftp2fn(efTRX,NFILE,fnm), index[0], ngx, interval, 
		 nslices, &diff, &variance);

  diffusion_plot(diff, variance, opt2fn("-o",NFILE,fnm), nslices);

  thanx(stdout);
}








