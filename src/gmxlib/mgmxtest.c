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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_mgmxtest_c = "$Id$";

#include "copyrite.h"
#include "statutil.h"
#include "macros.h"

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "g_chi computes phi, psi and chi dihedrals for all your sidechains.",
    "Output is in form of xvgr files, as well as a LaTeX table of the",
    "number of transitions per nanosecond.[PAR]",
    "Order parameters S2 for each of the dihedrals are calculated and",
    "output as xvgr file and optionally as a pdb file with the S2",
    "values as B-factor.[PAR]",
    "If option [TT]-c[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity.[PAR]",
    "g_chi computes phi, psi and chi dihedrals for all your sidechains.",
    "Output is in form of xvgr files, as well as a LaTeX table of the",
    "number of transitions per nanosecond.[PAR]",
    "Order parameters S2 for each of the dihedrals are calculated and",
    "output as xvgr file and optionally as a pdb file with the S2",
    "values as B-factor.[PAR]",
    "If option [TT]-c[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity.[PAR]",
    "g_chi computes phi, psi and chi dihedrals for all your sidechains.",
    "Output is in form of xvgr files, as well as a LaTeX table of the",
    "number of transitions per nanosecond.[PAR]",
    "Order parameters S2 for each of the dihedrals are calculated and",
    "output as xvgr file and optionally as a pdb file with the S2",
    "values as B-factor.[PAR]",
    "If option [TT]-c[tt] is given, the program will",
    "calculate dihedral autocorrelation functions. The function used",
    "is C(t) = < cos(chi(tau)) cos(chi(tau+t)) >. The use of cosines",
    "rather than angles themselves, resolves the problem of periodicity."
  };
  static char *bugs[] = {
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output.",
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output.",
    "Produces MANY output files (up to about 4 times the number of residues in the protein, twice that if autocorrelation functions are calculated). Typically several hundred files are output."
  };
  int    dummy;
  static int  r0=1,ndeg=1,maxchi=2,nf=10;
  static bool bAll=FALSE;
  static bool bPhi=FALSE,bPsi=FALSE,bChi=TRUE;
  static real bfac_init=-1.0;
  static bool bRama=FALSE,bShift=TRUE;
  static char *maxchistr[] = { "0", "1", "2", "3",  "4", "5", "6", NULL };

  t_pargs pa[] = {
    { "-r0",  FALSE, etINT, &r0,
      "Starting residue number" },
    { "-phi",  FALSE, etBOOL, &bPhi,
      "Output for Phi dihedral angles" },
    { "-psi",  FALSE, etBOOL, &bPsi,
      "Output for Psi dihedral angles" },
    { "-chi",  FALSE, etBOOL, &bChi,
      "Output for Chi dihedral angles" },
    { "-rama", FALSE, etBOOL, &bRama,
      "Generate Phi/Psi and Chi1/Chi2 ramachandran plots" },
    { "-all",  FALSE, etBOOL, &bAll,
      "Output separate files for every dihedral." },
    { "-nframes", FALSE, etINT, &nf,
      "Number of frames in your trajectory" },
    { "-shift", FALSE, etBOOL, &bShift,
      "Compute chemical shifts from Phi/Psi angles" },
    { "-run", FALSE, etINT, &ndeg,
      "Perform running average over ndeg degrees for histograms" },
    { "-maxchi", FALSE, etENUM, maxchistr,
      "Calculate first ndih Chi dihedrals" },
    { "-bfact", FALSE, etREAL, &bfac_init,
      "B-Factor value for pdb file for atoms with no calculated dihedral order parameter"}
  };

  t_filenm  fnm[] = {
    { efTPX, NULL,  NULL,     ffREAD  },
    { efTRX, "-f",  NULL,     ffREAD  },
    { efXVG, "-o",  "order",  ffWRITE },
    { efPDB, "-p",  "order",  ffOPTWR },
    { efXVG, "-jc", "Jcoupling", ffWRITE },
    { efXVG, "-c",  "dihcorr",ffOPTWR },
    { efTEX, "-t",  "trans",  ffWRITE },
    { efLOG, "-g",  "chi",    ffWRITE }
  };
#define NFILE asize(fnm)

  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs);
  
  fprintf(stderr,"You made it back in main of mgmxtest!\n");
  fprintf(stderr,"Please enter an integer and press return:");
  scanf("%d",&dummy);
  
  thanx(stdout);
    
  return 0;
}
