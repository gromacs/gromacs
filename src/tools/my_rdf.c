/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"

#define RESOLUTION 0.01   /* 0.01 nm resolution */
#define NR_HIST 200      /* compile histogram of 200 bins (from 0 to 2.0 nm)*/

static real hist[NR_HIST];
static real dens_tot; 

/****************************************************************************/
/* This program calculates density corrected and uncorrected radial         */
/* distribution functions in inhomogeneous systems like membranes.          */
/* Peter Tieleman, January 1996                                             */
/****************************************************************************/

real sphere_vol(real r)
{
  return (4.0*M_PI/3.0)*(r*r*r);
}

void calc_rdf(char *fn, atom_id **index, int gnx[], real l)
{
  rvec *x0;                /* coordinates without pbc */
  matrix box;              /* box (3x3) */
  int natoms,              /* nr. atoms in trj */
    status,   
    i,j,m,                 /* loop indices */
    teller = 0,       
    nr_frames = 0,         /* number of frames */
    nwater, nlipid;        /* number of waters, number of lipids */
  real t, 
    rho_local; 
  int bin;                 /* bin to put atom in */
  atom_id *water, *lipid;  /* the index numbers for lipid and water */
  rvec tmp;
  real dist;
  real count = 0; real dens = 0; 

  for (i = 0; i < NR_HIST; i++)
    hist[i] = 0;

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");
  
  fprintf(stderr,"Cut-off for counting is %f\n",l);
  nwater = gnx[1]; nlipid = gnx[0];

  rho_local = nwater/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
  fprintf(stderr,"Average water density: %8.3g\n"
	  "Number of lipids: %d. Number of waters: %d\n",
	  rho_local, nlipid, nwater);

  lipid = index[0];
  water = index[1];

  dens_tot = 0;

  /*********** Start processing trajectory ***********/
  do 
  {
    if ((teller++ % 10) == 0)
      fprintf(stderr,"\rFrame: %d",teller-1); 

    for (i = 0; i < nlipid; i++)
    {
      for (j = 0; j < nwater; j++)  /* nice double loop */
      {
	rvec_sub(x0[lipid[i]], x0[water[j]], tmp); 
	for (m = 0; m < DIM; m++)
	{
	  if (tmp[m] > 0.5*box[m][m])
	    tmp[m] -= box[m][m];
	  else if (tmp[m] < -0.5*box[m][m])
	    tmp[m] += box[m][m];
	}
	dist = norm(tmp);
	/*	fprintf(stderr,"dist: %f",dist); */
	if (dist < RESOLUTION*NR_HIST) 
	{
	  bin = trunc(dist/RESOLUTION);
	  if (bin < 0 || bin > NR_HIST)
	  {
	    fprintf(stderr,"Warning: bin out of range: %d\n",bin);
	    bin = NR_HIST;
	  }
	  hist[bin] += 1.0;
	}
	if (dist < 2.0)
	  dens += 1.0;

	if (dist < l)
	{
	  count += 1.0;
	}
      }
      dens_tot += dens;
      dens = 0;
    }
  } while (read_next_x(status,&t,natoms,x0,box));
  
  /*********** done with status file **********/
  close_trj(status);
  dens_tot /= nlipid * teller * sphere_vol(2);

  fprintf(stderr,"Local density of water around lipid: %f\n",dens_tot);

  for(i=1;i<NR_HIST;i++)
    hist[i] /= nlipid*teller*dens_tot;

  count /= nlipid*teller;
  fprintf(stderr,"Counted %g OW's\n", count); 
  sfree(x0);  /* free memory used by coordinate array */
}

real *integrate(real data[])
{
  int i,j;  
  real *result,
    sum;

  snew(result, NR_HIST);

  for (i = 1; i < NR_HIST ; i++)
  {
    sum = 0;
    for (j = 1; j < i; j++)
      sum += data[j] + 0.5 * (data[j+1] - data[j]);
    result[i] = sum * dens_tot;
  }
  return result;
}

void plot_rdf(char *afile, char *grpname[])
{
  FILE       *uncor;     /* xvgr file corrected rdf   */
  char       buf[256]; /* for xvgr title */
  int        i, n;
  real       *integ;

  sprintf(buf,"Uncorrected rdf");
  uncor = xvgropen(afile, buf, "r", "g(r)");

  integ = integrate(hist);

  for (i = 0; i < NR_HIST; i++)
  {
    hist[i]/= 4*M_PI*RESOLUTION*RESOLUTION*RESOLUTION*i*i; 
  /* r2 = RESOLUTION*RESOLUTION*i*i. The third comes from the change of the x-ax
     from i to i*RESOLUTION. sucks */
    fprintf(uncor,"%12g   %12g   %12g\n", i*RESOLUTION, hist[i], integ[i]); 
  }

  fclose(uncor);
}
 
void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute rdf's"
  };
  static char *opts[] = {
    "-l distance"
  };
  static char *odesc[] = {
    "distance taken as minimum for counting"
  };

  t_manual man = {asize(desc),desc,asize(opts),opts,NULL,NULL};
  real      l = 0;
  char      **grpname;            	    /* groupnames                 */
  int       i,ngrps = 0,                    /* nr. of groups              */
    *ngx;                                   /* sizes of groups            */
  atom_id   **index;             	    /* indices for all groups     */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efXVG, "-o", "graph",ffWRITE }, 	    /* xvgr output file 	  */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc, argv, PCA_CAN_TIME, NFILE, fnm, TRUE, &man);
  
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-l") == 0)  {
      if (i < argc - 1) {
	sscanf(argv[i+1],"%f",&l);
	i++;
      }
    }
  }

  snew(grpname,2);
  snew(index,2);
  snew(ngx,2);
 
  rd_index(ftp2fn(efNDX,NFILE,fnm),2,ngx,index,grpname); 

  fprintf(stderr,"Assuming group %s is the solvent.\n", grpname[1]); 

  calc_rdf(ftp2fn(efTRX,NFILE,fnm),index, ngx, l);
  
  plot_rdf(opt2fn("-o",NFILE,fnm), grpname);
}
