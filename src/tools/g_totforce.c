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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include "sysstuff.h"
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

#define EPS0 8.85419E-12
#define ELC 1.60219E-19


/* ************************************************* */
/* */
/* */ 

static int ce=0, cb=0;

/* this routine integrates the array data and returns the resulting array */
/* routine uses simple trapezoid rule                                     */
void integrate(real *result, real data[], int ndata, real slWidth)
{
  int i,j,  
      slice;
  real sum;
  
  if (ndata <= 2) 
    fprintf(stderr,"Warning: nr of slices very small. This will result"
	    "in nonsense.\n");

  fprintf(stderr,"Integrating from slice %d to slice %d\n",cb, ndata-ce);

  for (slice = cb; slice < (ndata-ce); slice ++)
  {
    sum = 0;
    for (i = cb; i < slice; i++)
      sum += slWidth * (data[i] + 0.5 * (data[i+1] - data[i]));
    result[slice] = sum;
  }
  return;
}

void calc_netforce(char *fn1,char *fn2, atom_id **index, int gnx[], 
		    real ***slPotential, real ***slCharge, 
		    real ***slField, int *nslices, 
		    t_topology *top, int axis, int nr_grps, real *slWidth,
		    real fudge_z, bool bSpherical)
{
  /*  rvec *x0; */             /* coordinates without pbc */
  rvec *f1,*f2;          /* forces in traj 1 and 2 */
  rvec df,*dftot;
  rvec *dfftot;
  matrix box;            /* box (3x3) */
  int natoms1,natoms2,   /* nr. atoms in trj */
      status1,status2,
        **slCount,       /* nr. of atoms in one slice for a group */
      i,j,n,             /* loop indices */
      teller = 0,      
      ax1, ax2,
      nr_frames = 0,     /* number of frames */
      slice;             /* current slice */
  real slVolume;         /* volume of slice for spherical averaging */
  real t, z, tm;
  rvec xcm;
  t_filenm  fnm[] = {
    { efTRJ, "-f1", NULL,  FALSE },
    { efTRJ, "-f2", NULL, FALSE },
  };
 
  FILE *fp1,*fp2;
  t_statheader sh1,sh2;

  
  fp1=opt2FILE("-f1",FILE,fnm,"r");
  fp2=opt2FILE("-f2",FILE,fnm,"r");
  rd_header(fp1,&sh1);
  rd_header(fp2,&sh2);
  
  /*
  if ((natoms1 = read_first_f(&status1,fn1,&t,&f1,box)) == 0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  if ((natoms2 = read_first_f(&status2,fn2,&t,&f2,box)) == 0) {
    fprintf(stderr,"Could not read coordinates from statusfile\n");
    exit(1);
  }
  */
  natoms1 = sh1.natoms;
  natoms2 = sh2.natoms;
  
  
  /* allocate memory */
  snew(df,1);
  snew(dftot,nr_grps);
  snew(dfftot,nr_grps);
  
  snew(f1,natoms1);
  snew(f2,natoms2);
  
  rewind(fp1);
  rewind(fp2);
  
  /*********** Start processing trajectory ***********/
  /* Berk code */
  do {
    rd_header(fp1,&sh1);
    rd_header(fp2,&sh2);
    
    /* Check whether there are forces in this frame... */
    if (sh.f_size != 0) {
      rd_hstatus(fp,&sh,&step,&t,&lambda,NULL,NULL,NULL,NULL,&natoms,
		 NULL,NULL,force,&nre,NULL,NULL);
      /* Analyse them */
      
      
      
      
      for (n = 0; n < nr_grps; n++)
        {      
          dftot[n]=0.0;
          for (i = 0; i < gnx[n]; i++)  /* loop over all atoms in index file */
            {
              rvec_sub(f1[index[n][i]],f2[index[n][i]],df);
              rvec_inc(dftot,df);
            }
        }
      nr_frames++;
      
    }
    else {
      /* Else skip the frame */
      rd_hstatus(fp,&sh,&step,&t,&lambda,NULL,NULL,NULL,NULL,&natoms,
		 NULL,NULL,NULL,&nre,NULL,NULL);
    }
  } while ((!eof(fp1)) && (!eof(fp2)));
  
  /*********** done with status file **********/
  
  fclose(fp1);
  fclose(fp2);
  
  /* slCharge now contains the total charge per slice, summed over all
     frames. Now divide by nr_frames and integrate twice 
     */
  
  fprintf(stderr,"\n\nRead %d frames from trajectory. Calculating force\n",
          nr_frames-1);
  
  for (n =0; n < nr_grps; n++)
    {
      dfftot[n]/=(nr_frames-1);
      
    }
  
  sfree(f1);  /* free memory used by coordinate array */
  sfree(f2);
  
  fprintf(stderr,"\n\nTotal forces are:\n");
  for (n =0; n < nr_grps; n++)
    {
     fprintf(stderr,"Group %d: %12.5f",n,dfftot[n]);
    } 
}


void plot_potential(real *potential[], real *charge[], real *field[], 
		    char *afile, char *bfile, char *cfile, int nslices,
		    int nr_grps, char *grpname[], real slWidth)
{
  FILE       *pot,     /* xvgr file with potential */
             *cha,     /* xvgr file with charges   */
             *fie;     /* xvgr files with fields   */
  char       buf[256]; /* for xvgr title */
  int        slice, n;

  sprintf(buf,"Electrostatic Potential");
  pot = xvgropen(afile, buf, "Box (nm)","Potential (V)");
  xvgr_legend(pot,nr_grps,grpname);

  sprintf(buf,"Charge Distribution");
  cha = xvgropen(bfile, buf, "Box (nm)", "Charge density (q/nm\\S3\\N)");
  xvgr_legend(cha,nr_grps,grpname);

  sprintf(buf, "Electric Field");
  fie = xvgropen(cfile, buf, "Box (nm)", "Field (V/nm)");
  xvgr_legend(fie,nr_grps,grpname);

  for (slice = cb; slice < (nslices - ce); slice++)
  { 
    fprintf(pot,"%12g  ", slice * slWidth);
    fprintf(cha,"%12g  ", slice * slWidth);
    fprintf(fie,"%12g  ", slice * slWidth);
    for (n = 0; n < nr_grps; n++)
    {
      fprintf(pot,"   %12g", potential[n][slice]);
      fprintf(fie,"   %12g", field[n][slice]);
      fprintf(cha,"   %12g", charge[n][slice]);
    }
    fprintf(pot,"\n");
    fprintf(cha,"\n");
    fprintf(fie,"\n");
  }

  fclose(pot);
  fclose(cha);
  fclose(fie);
}

void main(int argc,char *argv[])
{
  static char *desc[] = {
    "Compute the electrostatical potential across the box. The potential is"
    "calculated by first summing the charges per slice and then integrating"
    "twice of this charge distribution. Periodic boundaries are not taken  "
    "into account. Reference of potential is taken to be the left side of"
    "the box. It's also possible to calculate the potential in spherical"
    "coordinates as function of r by calculating a charge distribution in"
    "spherical slices and twice integrating them. epsilon_r is taken as 1,"
    "2 is more appropriate in many cases"
  };
  static int  axis = 2;                      /* normal to memb. default z  */
  static char *axtitle="Z"; 
  static int  nslices = 10;                  /* nr of slices defined       */
  static bool bSpherical = FALSE;            /* default is bilayer types   */
  static real fudge_z = 0;                    /* translate coordinates      */
  t_pargs pa [] = {
    { "-d",   FALSE, etSTR, &axtitle, 
      "Take the normal on the membrane in direction X, Y or Z." },
    { "-sl",  FALSE, etINT, &nslices,
      "Calculate potential as function of boxlength, dividing the box"
      " in #nr slices." } ,
    { "-cb",  FALSE, etINT, &cb,
      "Discard first #nr slices of box for integration" },
    { "-ce",  FALSE, etINT, &ce,
      "Discard last #nr slices of box for integration" },
    { "-tz",  FALSE, etREAL, &fudge_z,
      "Translate all coordinates <distance> in the direction of the box" },
    { "-spherical", FALSE, etBOOL,
      "Calculate spherical thingie" },
  };
  static char *bugs[] = {
    "Discarding slices for integration should not be necessary."
  };

  real     **potential,                    /* potential per slice        */
            **charge,                       /* total charge per slice     */
            **field,                        /* field per slice            */
            slWidth;                        /* width of one slice         */
  char      **grpname;            	    /* groupnames                 */
  int       ngrps = 0,                      /* nr. of groups              */
            i,
            *ngx;                           /* sizes of groups            */
  t_topology *top;                	    /* topology 		  */ 
  atom_id   **index;             	    /* indices for all groups     */
  t_filenm  fnm[] = {             	    /* files for g_order 	  */
    { efTRX, "-f1", NULL,  ffREAD },        /* trajectory file 1          */
    { efTRX, "-f2", NULL,  ffREAD },        /* trajectory file 2          */
    { efNDX, NULL, NULL,  ffREAD },    	    /* index file 		  */
    { efTPB, NULL, NULL,  ffREAD },    	    /* topology file           	  */
    { efXVG,"-o","force", ffWRITE },        /* xvgr output file 	  */
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME, TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);

  /* Calculate axis */
  /*  axis = toupper(axtitle[0]) - 'X';
   */
  top = read_top(ftp2fn(efTPB,NFILE,fnm));     /* read topology file */

  printf("How many groups? ");
  do { scanf("%d",&ngrps); } while (ngrps <= 0);
  
  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);
 
  rd_index(ftp2fn(efNDX,NFILE,fnm),ngrps,ngx,index,grpname); 

  
  calc_netforce(opt2fn("-f1",NFILE,fnm),opt2fn("-f2",NFILE,fnm), index, ngx, 
		 &potential, &charge, &field,
		 &nslices, top, axis, ngrps, &slWidth, fudge_z,
		 bSpherical); 

  plot_potential(potential, charge, field, opt2fn("-o",NFILE,fnm),
		 opt2fn("-oc",NFILE,fnm), opt2fn("-of",NFILE,fnm),
		 nslices, ngrps, grpname, slWidth);

  xvgr_file(opt2fn("-o",NFILE,fnm), NULL);       /* view xvgr file */
  xvgr_file(opt2fn("-oc",NFILE,fnm), NULL);      /* view xvgr file */  
  xvgr_file(opt2fn("-of",NFILE,fnm), NULL);      /* view xvgr file */

  thanx(stdout);
}









