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
 * GROningen MAchine for Chemical Simulation
 */
#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

void calc_dist(char *xvg, char *fn, atom_id **index, int ngx[], int axis, int ngrps,
	       char *grpname[])
{
  FILE  *out;            /* output file */
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  int natoms,            /* nr. atoms in trj */
      status,  
      nframes=0,         /* frame counter */
      i,j,n;             /* loop indices */
  real t,                /* four z's per frame */ 
       z1_tot,           /* total of all coords of first group of pair */
       z2_tot,           /* total of all coords of second group of pair*/
       z1_sqr,           /* sum of squared values for z1 */
       z2_sqr;           /* sum of squared values for z2 */
  char buf[256],**legend;         /* legend for xvg-file */

  sprintf(buf,"Distances between groups in lipid bilayer");
  out = xvgropen(xvg, buf, "Time (t)","Distance (nm)");
  snew(legend,ngrps/2);
  for(i=0;i<ngrps/2;i++)
  {
    sprintf(buf,"%s-%s",grpname[2*i],grpname[2*i+1]);
    legend[i] = strdup(buf);
  }
  xvgr_legend(out,ngrps/2,legend);

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  
  /*********** Start processing trajectory ***********/
  do 
  {
    if ((nframes++ % 10) == 0)
      fprintf(stderr,"\rFrame: %d",nframes-1); 
  
    fprintf(out,"%g   ",t);
    for (n = 0; n < ngrps; n=n+2)
    {      
      z1_tot = 0; z2_tot = 0; z1_sqr = 0; z2_sqr = 0;
      for (i = 0; i < ngx[n]; i++)   /* loop over all atoms in index file */
      {
	z1_tot += x0[index[n][i]][axis];
	z2_tot += x0[index[n+1][i]][axis];
	z1_sqr += x0[index[n][i]][axis] * x0[index[n][i]][axis];
	z2_sqr += x0[index[n+1][i]][axis] * x0[index[n+1][i]][axis]; 
      }
      fprintf(out," %e   ",fabs((z1_tot - z2_tot)/ngx[n]));
    }
    fprintf(out,"\n");
  } while (read_next_x(status,&t,natoms,x0,box));
  
  /*********** done with status file **********/
  close_trj(status);
  fclose(out);
  sfree(x0);  /* free memory used by coordinate array */
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute the distances between the average coordinate of two groups of",
    "atoms. This is useful for computing the P-P or N-N distances in lipid",
    "bilayers.",
  };
  static char *opts[] = {
    "-d X | Y | Z",
  };
  static char *odesc[] = {
    "Take the normal on the membrane in direction X, Y or Z.",
  };
  static char *bugs[] = {
    ""
  };

  t_manual man = {asize(desc),desc,asize(opts),opts,odesc,NULL,NULL};
  char      **grpname;            	    /* groupnames                 */
  int       ngrps = 0,                      /* nr. of groups in index file*/
            i,
            *ngx,                           /* sizes of groups            */
            axis = 2;                       /* normal to memb. default z  */
  atom_id   **index;             	    /* indices for all groups     */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efXVG, "-o", NULL,  ffWRITE }, 	    /* xvgr output file 	  */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE,
		    fnm, TRUE, &man);
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i],"-d") == 0)  {
      if (i < argc - 1) {
	axis = (int)(argv[i+1][0] - 'X');
	i++;
      }
    }
  }
  
  printf("You can pick groups from the index. Distances will be calculated\n"
	 "between groups 1 and 2, 3 and 4, 5 and 6 etc.\n");
  printf("How many groups? ");
  do { scanf("%d",&ngrps); 
       if (ngrps % 2 == 1)
       {	 
	 fprintf(stderr,"Please given an even number of groups so I can make pairs.\n");
	 ngrps=-1;
       }
     } while (ngrps <= 0);
 
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);
 
  rd_index(ftp2fn(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 

  calc_dist(ftp2fn(efXVG,NFILE,fnm), ftp2fn(efTRX,NFILE,fnm), index, ngx, axis, 
	    ngrps, grpname); 

  thanx(stdout);
}









