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
 * Grunge ROck MAChoS
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
/* This program creates a graph of the projection of the vectors between    */
/* the middle of the C2-C3 bond and C11-C12 bond of lipids on a plane z=c   */
/* where z = average z-coordinate of the middle of the bonds involved in    */
/* defining the vector.                                                     */
/* The program first reads the trajectory to determine the average z of all */
/* tails involved. Then it reads the trajectory again to determine the      */
/* intersection of the vector for each bond with the previously determined  */
/* average z. These values are averaged and printed out.                    */
/*                                                                          */
/* Peter Tieleman, April 1995                                               */
/****************************************************************************/


/* Print name of first atom in all groups in index file */
static void print_types(atom_id index[], atom_id a[], int ngrps, 
			char *groups[], t_topology *top)
{
  int i;

  fprintf(stderr,"Using following groups: \n");
  for(i = 0; i < ngrps; i++)
    fprintf(stderr,"Groupname: %s First atomname: %s First atomnr %d\n", 
	    groups[i], *(top->atoms.atomname[a[index[i]]]), a[index[i]]);
  fprintf(stderr,"\n");
}

void find_orientations(char *fn, atom_id *index, atom_id *a, real *x, 
		       real *y, t_topology *top, int axis, int nr_tails, 
		       matrix box)
{
  rvec *x0,              /* coordinates with pbc */
       *x1,              /* coordinates without pbc */
       dist,dist2;       /* vector between two bonds */
  int  status,teller;
  real z=0,z_ave=0,      /* average z (used for plane) */
       z2_ave=0,t;                /* time (from coordinate frame) */
  int natoms,            /* nr. atoms in trj */
      i,j   = 0,
     nr_frames = 0;

  /* First find average z of all tails, over all frames  */
  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  fprintf(stderr,"Number of elements in first group: %d\n",nr_tails);
  /* first group as standard. Not rocksolid, but might catch error in index*/

  snew(x1, natoms);
  teller = 0; 

  /*********** Start processing trajectory ***********/
  do 
  {
    teller++;
    fprintf(stderr,"\rFrame: %d",teller-1); 
    /* remove pbc */
    rm_pbc(&(top->idef),top->atoms.nr,box,x0,x1);

    /* Now loop over all groups. There are 4 groups, take middle from vector
       0[i]-1[i] and middle of 2[i]-3[i].
       over all 
       atoms in group, which is index[i] to (index[i+1] - 1) See block.h. Of 
       course, in this case index[i+1] -index[i] has to be the same for all 
       groups, namely the number of tails. 
     */
    
    for (i = 1; i < 4; i++)
      if (index[i+1] - index[i] != nr_tails)
      {
	fprintf(stderr,"ERROR: grp %d does not have same number of"
		" elements as grp 1\n",i); 
	exit(1);
      }
    
    for (j = 0; j < nr_tails; j++)
    {
      /* get vector dist(Cn-1,Cn+1) for tail atoms */
      rvec_add(x1[a[index[0]+j]], x1[a[index[1]+j]], dist);
      rvec_add(x1[a[index[2]+j]], x1[a[index[3]+j]], dist2);
      svmul(0.5, dist, dist);
      svmul(0.5, dist2, dist2);
      z = 0.5*(dist[2]+dist2[2]); 
      z_ave+=z;
      z2_ave+=(z*z);
    }   /* end loop j, over all atoms in group */
  } while (read_next_x(status,&t,natoms,x0,box));
  /*********** done with status file **********/

  z_ave = z_ave/(nr_tails * (teller+1));
  z2_ave = z2_ave/(nr_tails * (teller+1));

  fprintf(stderr,"\nThe average coordinate is: %f\n"
	  "<(Delta x)2> = %f\n"
	  "error: %f\n",z_ave,z2_ave-z_ave*z_ave,
	  sqrt((z2_ave-z_ave*z_ave)/sqrt(nr_tails*(teller+1))));

  sfree(x0);  /* free memory used by coordinate arrays */
  sfree(x1);
}

void make_plot(real x[], real y[], char *afile, matrix box, int nr_tails)
{
  FILE       *out;          /* xvgr files with order parameters  */
  char       buf[256];               /* for xvgr title */
  int        i;

  sprintf(buf,"Tail projection");
  out = xvgropen(afile,buf,"x (nm)","y (nm)");

  xvgr_line_props(out,0,elNone,ecFrank);
  xvgr_view(out,0.2,0.2,0.8,0.8);
  xvgr_world(out,0,0,box[XX][XX],box[YY][YY]);
  fprintf(out,"@ s0 symbol 2\n@ s0 symbol size 0.4\n@ s0 symbol fill 1\n");

  for (i = 0; i < nr_tails; i++)
    fprintf(out,"%8.3f       %8.3f\n", x[i], y[i]);

  fclose(out);
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute the orientation of a lipid tail with respect to an interface.",
    "Output is the projection of all tails along a vector defined by the",
    "bonds C2-C3 and C11-C12, where C1 is the carbonyl atom. This gives a",
    "plot in which hexagonal packings are easily recognized."
  };
  static char *opts[] = {
    "-d X | Y | Z",
  };
  static char *odesc[] = {
    "Take the normal on the membrane in direction X, Y or Z.",
  };

  t_manual man = {asize(desc),desc,asize(opts),opts,odesc,NULL,NULL};
  FILE      *status;              	    /* trajectory file  	  */
  char      **grpname;            	    /* groupnames                 */
  int       nr_tails,ngrps,                /* nr. of groups              */
            i,
            axis = 2;                       /* normal to memb. default z  */
  real      *x, *y;                         /* values for graph           */
  matrix box;                               /* box sizes                  */
  t_topology *top;                	    /* topology 		  */ 
  atom_id   *index,             	    /* indices for a              */
            *a;                             /* atom numbers in each group */
  t_block   *block;                         /* data from index file       */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL, ffREAD },    	    /* trajectory file 	          */
    { efNDX, "-n", NULL, ffREAD },    	    /* index file 		  */
    { efTPB, "-s", NULL, ffREAD },    	    /* topology file           	  */
    { efXVG, "-o", NULL, ffWRITE }, 	    /* xvgr output file 	  */
  };

#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, NFILE,
		    fnm, TRUE, &man);
  for (i = 1; i < argc; i++)
  {
    if (strcmp(argv[i],"-d") == 0) 
    {
      if (i < argc - 1)
      {
	axis = (int)(argv[i+1][0] - 'X');
        i++;
      }
    }
  }
  
  fprintf(stderr,"Reading topology\n");
  top = read_top(ftp2fn(efTPB,NFILE,fnm));     /* read topology file */

  fprintf(stderr,"Reading indexfile\n");
  block = init_index(ftp2fn(efNDX,NFILE,fnm),&grpname);
  index = block->index;                  /* get indices from t_block block */
  a = block->a;                          /* see block.h                    */
  ngrps = block->nr;           
  nr_tails = index[1] - index[0];

  if (ngrps != 4)
  {
    fprintf(stderr,"I want 4 groups in the indexfile! Bye!\n");
    exit(1); 
  }

  /* show atomtypes, to check if index file is correct */ 
  print_types(index, a, ngrps, grpname, top);

  /* real program */
  find_orientations(ftp2fn(efTRX,NFILE,fnm), index, a, 
		    x, y, top, axis, nr_tails, box); 
  make_plot(x, y, opt2fn("-o",NFILE,fnm), box, nr_tails);
  xvgr_file(opt2fn("-o",NFILE,fnm), NULL);      /* view xvgr file */
  thanx(stdout);
}









