/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_density_c = "$Id$";

#include <math.h>
#include <ctype.h>
#include "sysstuff.h"
#include "string.h"
#include "string2.h"
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
#include "rdgroup.h"
#include "tpxio.h"

typedef struct {
  char *atomname;
  int nr_el;
} t_electron;

/****************************************************************************/
/* This program calculates the partial density across the box.              */
/* Peter Tieleman, Mei 1995                                                 */
/****************************************************************************/

/* used for sorting the list */
int compare(void *a, void *b)
{
  t_electron *tmp1,*tmp2;
  tmp1 = (t_electron *)a; tmp2 = (t_electron *)b;

  return strcmp(tmp1->atomname,tmp2->atomname);
}

int get_electrons(t_electron **eltab, char *fn)
{
  char buffer[256];  /* to read in a line   */
  char tempname[80]; /* buffer to hold name */
  int tempnr; 

  FILE *in;
  int nr;            /* number of atomstypes to read */
  int i;

  if ( !(in = fopen(fn,"r")))
    fatal_error(0,"Couldn't open %s. Exiting.\n",fn);

  fgets(buffer, 255, in);
  if (sscanf(buffer, "%d", &nr) != 1)
    fatal_error(0,"Invalid number of atomtypes in datafile\n");

  snew(*eltab,nr);

  for (i=0;i<nr;i++) {
    if (fgets(buffer, 255, in) == NULL)
      fatal_error(0,"reading datafile. Check your datafile.\n");
    if (sscanf(buffer, "%s = %d", tempname, &tempnr) != 2)
      fatal_error(0,"Invalid line in datafile at line %d\n",i+1);
    (*eltab)[i].nr_el = tempnr;
    (*eltab)[i].atomname = strdup(tempname);
  }
  
  /* sort the list */
  fprintf(stderr,"Sorting list..\n");
  qsort ((void*)*eltab, nr, sizeof(t_electron), 
	 (int(*)(const void*, const void*))compare);

  return nr;
}

void calc_electron_density(char *fn, atom_id **index, int gnx[], 
			   real ***slDensity, int *nslices, t_topology *top, 
			   int axis, int nr_grps, real *slWidth, 
			   t_electron eltab[], int nr)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  int natoms,            /* nr. atoms in trj */
      status,  
      i,n,               /* loop indices */
      ax1=0, ax2=0,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  t_electron *found;     /* found by bsearch */
  t_electron sought;     /* thingie thought by bsearch */
 
  real t, 
        z;

  switch(axis) {
  case 0:
    ax1 = 1; ax2 = 2;
    break;
  case 1:
    ax1 = 0; ax2 = 2;
    break;
  case 2:
    ax1 = 0; ax2 = 1;
    break;
  default:
    fatal_error(0,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0)
    fatal_error(0,"Could not read coordinates from statusfile\n");
  
  if (! *nslices)
    *nslices = (int)(box[axis][axis] * 10); /* default value */
  fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);

  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);
  
  /*********** Start processing trajectory ***********/
  do {
    rm_pbc(&(top->idef),top->atoms.nr,box,x0,x0);

    *slWidth = box[axis][axis]/(*nslices);
    for (n = 0; n < nr_grps; n++) {      
      for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
	  z = x0[index[n][i]][axis];
	  if (z < 0) 
	    z += box[axis][axis];
	  if (z > box[axis][axis])
	    z -= box[axis][axis];
      
	  /* determine which slice atom is in */
	  slice = (z / (*slWidth)); 
	  sought.nr_el = 0;
	  sought.atomname = strdup(*(top->atoms.atomname[index[n][i]]));

	  /* now find the number of electrons. This is not efficient. */
	  found = (t_electron *)
	    bsearch((const void *)&sought,
		    (const void *)eltab, nr, sizeof(t_electron), 
		    (int(*)(const void*, const void*))compare);

	  if (found == NULL)
	    fprintf(stderr,"Couldn't find %s. Add it to the .dat file\n",
		    *(top->atoms.atomname[index[n][i]]));
	  else  
	    (*slDensity)[n][slice] += found->nr_el - 
	                              top->atoms.atom[index[n][i]].q;
	  /* take into account partial charge */
	  free(sought.atomname);
	}
    }
      nr_frames++;
  } while (read_next_x(status,&t,natoms,x0,box));
  
  /*********** done with status file **********/
  close_trj(status);
  
/* slDensity now contains the total number of electrons per slice, summed 
   over all frames. Now divide by nr_frames and volume of slice 
*/

  fprintf(stderr,"\nRead %d frames from trajectory. Counting electrons\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++)
      (*slDensity)[n][i] = (*slDensity)[n][i] * (*nslices) /
	( nr_frames * box[axis][axis] * box[ax1][ax1] * box[ax2][ax2]);
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void calc_density(char *fn, atom_id **index, int gnx[], 
		  real ***slDensity, int *nslices, t_topology *top, 
		  int axis, int nr_grps, real *slWidth, bool bNumber,
		  bool bCount)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  int natoms,            /* nr. atoms in trj */
      status,  
      **slCount,         /* nr. of atoms in one slice for a group */
      i,j,n,               /* loop indices */
      teller = 0,      
      ax1=0, ax2=0,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  real t, 
        z;
  char *buf;             /* for tmp. keeping atomname */

  switch(axis) {
  case 0:
    ax1 = 1; ax2 = 2;
    break;
  case 1:
    ax1 = 0; ax2 = 2;
    break;
  case 2:
    ax1 = 0; ax2 = 1;
    break;
  default:
    fatal_error(0,"Invalid axes. Terminating\n");
  }

  if ((natoms = read_first_x(&status,fn,&t,&x0,box)) == 0)
    fatal_error(0,"Could not read coordinates from statusfile\n");
  
  if (! *nslices) {
    *nslices = (int)(box[axis][axis] * 10); /* default value */
    fprintf(stderr,"\nDividing the box in %d slices\n",*nslices);
  }
  
  snew(*slDensity, nr_grps);
  for (i = 0; i < nr_grps; i++)
    snew((*slDensity)[i], *nslices);
  
  /*********** Start processing trajectory ***********/
  do {
    rm_pbc(&(top->idef),top->atoms.nr,box,x0,x0);

    *slWidth = box[axis][axis]/(*nslices);
    teller++;
    
    for (n = 0; n < nr_grps; n++) {      
      for (i = 0; i < gnx[n]; i++) {   /* loop over all atoms in index file */
	z = x0[index[n][i]][axis];
	if (z < 0) 
	  z += box[axis][axis];
	if (z > box[axis][axis])
	  z -= box[axis][axis];
      
	/* determine which slice atom is in */
	slice = (int)(z / (*slWidth)); 
	if (bNumber || bCount) {
	  buf = strdup(*(top->atoms.atomname[index[n][i]]));
	  trim(buf);
	  if (buf[0] != 'H')
	    (*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m;
	  free(buf);
	} else
	  (*slDensity)[n][slice] += top->atoms.atom[index[n][i]].m;
      }
    }

    nr_frames++;
  } while (read_next_x(status,&t,natoms,x0,box));
  
  /*********** done with status file **********/
  close_trj(status);
  
  /* slDensity now contains the total mass per slice, summed over all
     frames. Now divide by nr_frames and volume of slice 
     */
  
  fprintf(stderr,"\nRead %d frames from trajectory. Calculating density\n",
	  nr_frames);

  for (n =0; n < nr_grps; n++) {
    for (i = 0; i < *nslices; i++) {
      if (bCount) 
	(*slDensity)[n][i] = (*slDensity)[n][i]/nr_frames;
      else
	(*slDensity)[n][i] = (*slDensity)[n][i] * (*nslices) * 1.66057 /
	(nr_frames * box[axis][axis] * box[ax1][ax1] * box[ax2][ax2]);
    }
  }

  sfree(x0);  /* free memory used by coordinate array */
}

void plot_density(real *slDensity[], char *afile, int nslices,
		  int nr_grps, char *grpname[], real slWidth, 
		  bool bElectron, bool bNumber, bool bCount)
{
  FILE       *den;     /* xvgr file with density   */
  char       buf[256]; /* for xvgr title */
  int        slice, n;

  sprintf(buf,"Partial densities");
  if (bElectron)
    den = xvgropen(afile, buf, "Box (nm)", "Electron density (e/nm\\S3\\N)");
  else if (bNumber)
    den = xvgropen(afile, buf, "Box (nm)","Density (atoms/nm\\S3\\N)");
  else if (bCount)
    den = xvgropen(afile, buf, "Slice","Absolute numbers");
  else
    den = xvgropen(afile, buf, "Box (nm)","Density (kg/m\\S3\\N)");

  xvgr_legend(den,nr_grps,grpname);

  for (slice = 0; slice < nslices; slice++) { 
    fprintf(den,"%12g  ", slice * slWidth);
    for (n = 0; n < nr_grps; n++)
      if (bNumber)
	fprintf(den,"   %12g", slDensity[n][slice]/1.66057);
      else
	fprintf(den,"   %12g", slDensity[n][slice]);
    fprintf(den,"\n");
  }

  fclose(den);
}
 
int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute partial densities across the box, using an index file. Densities",
    "in gram/cubic centimeter, number densities or electron densities can be",
    "calculated. For electron densities, each atom is weighed by its atomic",
    "partial charge."
  };
  
  static bool bNumber=FALSE;                  /* calculate number density   */
  static bool bElectron=FALSE;                /* calculate electron density */
  static bool bCount=FALSE;                   /* just count                 */
  static int  axis = 2;                       /* normal to memb. default z  */
  static char *axtitle="Z"; 
  static int  nslices = 10;                    /* nr of slices defined       */
  t_pargs pa[] = {
    { "-d", FALSE, etSTR, {&axtitle}, 
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-sl",  FALSE, etINT, {&nslices},
      "Divide the box in #nr slices." },
    { "-number",  FALSE, etBOOL, {&bNumber},
      "Calculate number density instead of mass density. Hydrogens are not counted!" },
    { "-ed",      FALSE, etBOOL, {&bElectron},
      "Calculate electron density instead of mass density" },
    { "-count",   FALSE, etBOOL, {&bCount},
      "Only count atoms in slices, no densities. Hydrogens are not counted"}
  };

  static char *bugs[] = {
    "When calculating electron densities, atomnames are used instead of types. This is bad.",
    "When calculating number densities, atoms with names that start with H are not counted. This may be surprising if you use hydrogens with names like OP3."
  };
  
  real     **density,                       /* density per slice          */
    slWidth;                        /* width of one slice         */
  char      **grpname;            	    /* groupnames                 */
  int       ngrps = 0,                      /* nr. of groups              */
    nr_electrons,                   /* nr. electrons              */
    *ngx;                           /* sizes of groups            */
  t_topology *top;                	    /* topology 		  */ 
  atom_id   **index;             	    /* indices for all groups     */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f", NULL,  ffREAD },    	    /* trajectory file 	          */
    { efNDX, NULL, NULL,  ffOPTRD },    	    /* index file 		  */
    { efTPX, NULL, NULL,  ffREAD },    	    /* topology file           	  */
    { efDAT, "-ei", "electrons", ffWRITE },   /* file with nr. of electrons */
    { efXVG,"-o","density",ffWRITE }, 	    /* xvgr output file 	  */
  };
  t_electron *el_tab;                       /* tabel with nr. of electrons*/
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  /* Calculate axis */
  axis = toupper(axtitle[0]) - 'X';
  
  top = read_top(ftp2fn(efTPX,NFILE,fnm));     /* read topology file */
  if ( bNumber  || bCount) {
    int n;
    for(n=0;(n<top->atoms.nr);n++)
      top->atoms.atom[n].m=1;  
    /* this sucks! mass is in amu, so a factor 1.66 is missing */
  }

  printf("How many groups? ");
  do { scanf("%d",&ngrps); } while (ngrps <= 0);
  
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);
 
  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 

  if (bElectron) {
    if (bCount)
      fatal_error(0,"I don't feel like counting electrons. Bye.\n");

    nr_electrons =  get_electrons(&el_tab,ftp2fn(efDAT,NFILE,fnm));
    fprintf(stderr,"Read %d atomtypes from datafile\n", nr_electrons);

    calc_electron_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, 
			  &nslices, top, axis, ngrps, &slWidth, el_tab, 
			  nr_electrons);
  } else
    calc_density(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, &nslices, top, 
		 axis, ngrps, &slWidth, bNumber,bCount); 
  
  plot_density(density, opt2fn("-o",NFILE,fnm),
	       nslices, ngrps, grpname, slWidth, bElectron, bNumber, bCount);
  
  xvgr_file(opt2fn("-o",NFILE,fnm), NULL);       /* view xvgr file */
  thanx(stdout);
  return 0;
}



