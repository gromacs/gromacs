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
#include <copyrite.h>
#include <statutil.h>
#include <statusio.h>
#include <smalloc.h>


void main(int argc,char *argv[])
{
  bool bCont;
  rvec *x,*v;
  matrix box;
  t_statheader header;

  static char *desc[0];
  static char *opts[0];
  static char *odesc[0];

  t_manual man = { asize(desc),desc,asize(opts),opts,odesc,0,NULL};

  t_filenm fnm[] = {
    { efTPB, NULL, NULL, ffREAD },
    { efTRJ, NULL, NULL, ffREAD },
    { efNDX, NULL, NULL, ffREAD }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		    NFILE,fnm,TRUE,&man);

  if (argc > 1) {
    if (strcmp(argv[1],"-k") == 0)
      bPBC=FALSE;
    else
      usage(argv[0],argv[1]);
  }
  read_status_header(ftp2fn(efTPB,NFILE,fnm),&header);

  snew(x,header.natoms);
  snew(v,header.natoms);

  read_status(ftp2fn(efTPB,NFILE,fnm),
              &step,&t,&lambda,&ir,
              box,NULL,NULL,
              &natom,
              xp,NULL,NULL,&nre,NULL,
              &top);

  /*set box type*/
  if (ir.eBox==ebtTRUNCOCT)
    bTruncOct=TRUE;
  else
    bTruncOct=FALSE;
  init_pbc(box,bTruncOct);


  status=ftp2FILE(efTRJ,NFILE,fnm,"r");

  /* */
  read_first_x_v(status,&t,&x,&v,box);
  do {
    /* */
    fprintf(stderr,"hallo");
    bCont = read_next_x_v(status,&t,natom,x,v,box);
  } while (bCont);
      
  fclose(status);
}
