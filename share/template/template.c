/*
 *       $Id$
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_template_c = "$Id$";
#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"


int main(int argc,char *argv[])
{
  static char *desc[] = {
    "this is a small test program meant to serve as a template ",
    "when writing your own analysis tools. The advantage of ",
    "using gromacs for this is that you have access to all ",
    "information in the topology, and your program will be ",
    "able to handle all types of coordinates and trajectory ",
    "files supported by gromacs. Go ahead and try it! ",
    "This test version just writes the coordinates of an ",
    "arbitrary atom to standard out for each frame. You can ",
    "select which atom you want to examine with the -n argument."
  };
  
  static int n=1;

  /* Extra arguments - but note how you always get the begin/end
   * options when running the program, without mentioning them here!
   */
  
  t_pargs pa[] = {
    { "-n", FALSE, etINT, {&n},
      "Plot data for atom number n (starting on 1)"
    }
  };
  
  t_topology top;
  char       title[STRLEN];
  t_trxframe fr;
  rvec       *xtop;
  matrix     box;
  int        status,flags;
 
  t_filenm fnm[] = {
    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    { efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  };
  
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);

  /* This is the routine responsible for adding default options,
   * calling the X/motif interface, etc. */
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  /* We don't need any topology information to write the coordinates,
   * but to show how it works we start by writing the name and
   * charge of the selected atom. It returns a boolean telling us
   * whether the topology was found and could be read
   */
  
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);
  sfree(xtop);

  n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  /* check that this atom exists */
  if(n<0 || n>(top.atoms.nr)) 
  {
    printf("Error: Atom number %d is out of range.\n",n);
    exit(1);
  }
  
  printf("Atom name: %s\n",*(top.atoms.atomname[n]));
  printf("Atom charge: %f\n",top.atoms.atom[n].q);
  
  /* The first time we read data is a little special */
  read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);
  
  /* This is the main loop over frames */
  do {
    /* coordinates are available in the vector fr.x
     * you can find this and all other structures in
     * the types directory under the gromacs include dir.
     */
    printf("Coordinates at t=%8.3f : %8.5f %8.5f %8.5f\n",fr.time,fr.x[n][XX],fr.x[n][YY],fr.x[n][ZZ]);
  } while(read_next_frame(status,&fr));

  thanx(stderr);
  
  return 0;
}

