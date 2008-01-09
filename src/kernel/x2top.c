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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "vec.h"
#include "x2top_nm2type.h"
#include "x2top_core.h"
#include "x2top_qgen.h"
#include "atomprop.h"
#include "grompp.h"
#include "add_par.h"
#include "gmx_random.h"

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "x2top generates a primitive topology from a coordinate file.",
    "The program assumes all hydrogens are present when defining",
    "the hybridization from the atom name and the number of bonds.",
    "The program can also make an rtp entry, which you can then add",
    "to the rtp database.[PAR]",
    "When [TT]-param[tt] is set, equilibrium distances and angles",
    "and force constants will be printed in the topology for all",
    "interactions. The equilibrium distances and angles are taken",
    "from the input coordinates, the force constant are set with",
    "command line options."
    "The force fields supported currently are:[PAR]",
    "G43a1  GROMOS96 43a1 Forcefield (official distribution)[PAR]",
    "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
    "G43b1  GROMOS96 43b1 Vacuum Forcefield (official distribution)[PAR]",
    "gmx    Gromacs Forcefield (a modified GROMOS87, see manual)[PAR]",
    "G43a2  GROMOS96 43a2 Forcefield (development) (improved alkane dihedrals)[PAR]",
    "The corresponding data files can be found in the library directory",
    "with names like ffXXXX.YYY. Check chapter 5 of the manual for more",
    "information about file formats. By default the forcefield selection",
    "is interactive, but you can use the [TT]-ff[tt] option to specify",
    "one of the short names above on the command line instead. In that",
    "case pdb2gmx just looks for the corresponding file.[PAR]",
    "An optional file containing atomname charge polarizability can be",
    "given with the [TT]-d[tt] flag."
  };
  static char *bugs[] = {
    "The atom type selection is primitive. Virtually no chemical knowledge is used",
    "No improper dihedrals are generated",
    "The atoms to atomtype translation table is incomplete (ffG43a1.n2t",
    "file in the $GMXLIB directory). Please extend it and send the results",
    "back to the GROMACS crew.[PAR]",
    "Charges can be read in (when based on a QM calculation) or generated",
    "using a variety of empirical algorithms (option [TT]-qgen[tt])."
  };
  FILE       *fp;
  t_params   plist[F_NRE];
  t_excls    *excls;
  t_atoms    *atoms;       /* list with all atoms */
  t_atomtype *atype;
  t_nextnb   nnb;
  t_nm2type  *nm2t;
  t_mols     mymol;
  void       *atomprop;
  int        nnm;
  char       title[STRLEN],forcefield[32];
  rvec       *x;        /* coordinates? */
  int        *nbonds,*cgnr;
  int        bts[] = { 1,1,1,2 };
  matrix     box;          /* box length matrix */
  int        natoms;       /* number of atoms in one molecule  */
  int        nres;         /* number of molecules? */
  int        i,j,k,l,m,ndih,alg;
  real       mu;
  bool       bRTP,bTOP,bOPLS,bCharmm;
  t_symtab   symtab;
  t_q_alpha  *qa=NULL;
  int        nqa=0;
  real       cutoff,qtot,mtot;
  char       rtp[STRLEN];

  
  t_filenm fnm[] = {
    { efSTX, "-f", "conf", ffREAD  },
    { efTOP, "-o", "out",  ffOPTWR },
    { efRTP, "-r", "out",  ffOPTWR },
    { efDAT, "-d", "qpol", ffOPTRD }
  };
#define NFILE asize(fnm)
  static real scale = 1.1, kb = 4e5,kt = 400,kp = 5;
  static real btol=0.1,qtol=1e-3,fac=5.0;
  static real qtotref=0,muref=0;
  static int  nexcl = 3;
  static int  maxiter=100;
  static bool bRemoveDih = FALSE;
  static bool bParam = TRUE, bH14 = TRUE,bAllDih = FALSE,bRound = TRUE;
  static bool bPairs = TRUE, bPBC = TRUE;
  static bool bUsePDBcharge = FALSE,bVerbose=FALSE;
  static char *molnm = "ICE";
  static char *ff = "select";
  static char *qgen[] = { NULL, "None", "Linear", "Yang", "Bultinck", "SM", NULL };
  t_pargs pa[] = {
    { "-ff",     FALSE, etSTR, {&ff},
      "Select the force field for your simulation." },
    { "-v",      FALSE, etBOOL, {&bVerbose},
      "Generate verbose output in the top file." },
    { "-nexcl", FALSE, etINT,  {&nexcl},
      "Number of exclusions" },
    { "-H14",    FALSE, etBOOL, {&bH14}, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-alldih", FALSE, etBOOL, {&bAllDih}, 
      "Generate all proper dihedrals" },
    { "-remdih", FALSE, etBOOL, {&bRemoveDih}, 
      "Remove dihedrals on the same bond as an improper" },
    { "-pairs",  FALSE, etBOOL, {&bPairs},
      "Output 1-4 interactions (pairs) in topology file" },
    { "-name",   FALSE, etSTR,  {&molnm},
      "Name of your molecule" },
    { "-pbc",    FALSE, etBOOL, {&bPBC},
      "Use periodic boundary conditions." },
    { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
      "Use the B-factor supplied in a pdb file for the atomic charges" },
    { "-btol",  FALSE, etREAL, {&btol},
      "Relative tolerance for determining whether two atoms are bonded." },
    { "-param", FALSE, etBOOL, {&bParam},
      "Print parameters in the output" },
    { "-round",  FALSE, etBOOL, {&bRound},
      "Round off measured values" },
    { "-qgen",   FALSE, etENUM, {qgen},
      "HIDDENAlgorithm used for charge generation" },
    { "-qtol",   FALSE, etREAL, {&qtol},
      "HIDDENTolerance for assigning charge generation algorithm" },
    { "-qtot",   FALSE, etREAL, {&qtotref},
      "HIDDENNet charge on molecule when generating a charge" },
    { "-muref",  FALSE, etREAL, {&muref},
      "HIDDENReference dipole of molecule when optimizing charges in the SM (Van der Spoel and Van Maaren) model" },
    { "-maxiter",FALSE, etINT, {&maxiter},
      "HIDDENMax number of iterations for charge generation algorithm" },
    { "-fac",    FALSE, etREAL, {&fac},
      "HIDDENNot yet understood factor for generating charges" },
    { "-kb",    FALSE, etREAL, {&kb},
      "Bonded force constant (kJ/mol/nm^2)" },
    { "-kt",    FALSE, etREAL, {&kt},
      "Angle force constant (kJ/mol/rad^2)" },
    { "-kp",    FALSE, etREAL, {&kp},
      "Dihedral angle force constant (kJ/mol/rad^2)" }
  };
  
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,asize(bugs),bugs);
  /* Check the options */
  bRTP = opt2bSet("-r",NFILE,fnm);
  bTOP = opt2bSet("-o",NFILE,fnm);
  
  if (!bRTP && !bTOP)
    gmx_fatal(FARGS,"Specify at least one output file");
  
  if ((btol < 0) || (btol > 1)) 
    gmx_fatal(FARGS,"Bond tolerance should be between 0 and 1 (not %g)",
	      btol);
  if ((qtol < 0) || (qtol > 1)) 
    gmx_fatal(FARGS,"Charge tolerance should be between 0 and 1 (not %g)",
	      qtol);
    
  /* Read standard atom properties */
  atomprop = get_atomprop();
    
  /* Force field selection */
  if(!strncmp(ff,"select",6)) {
    /* Interactive forcefield selection */
    choose_ff(forcefield,sizeof(forcefield));
  } else {
    sprintf(forcefield,"ff%s",ff);
  }
  sprintf(rtp,"%s.rtp",forcefield);
  printf("Looking whether force field file %s exists\n",rtp);
  fclose(libopen(rtp));

  bOPLS   = (strcmp(forcefield,"ffoplsaa") == 0);
  bCharmm = (strcmp(forcefield,"ffcharmm") == 0);
    
  mymol.name = strdup(molnm);
  mymol.nr   = 1;
	
  /* Init parameter lists */
  init_plist(plist);
  
  /* Read coordinates */
  get_stx_coordnum(opt2fn("-f",NFILE,fnm),&natoms); 
  snew(atoms,1);
  
  /* make space for all the atoms */
  init_t_atoms(atoms,natoms,TRUE);
  snew(x,natoms);              

  read_stx_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,NULL,box);

  nm2t = rd_nm2type(forcefield,&nnm);
  printf("There are %d name to type translations\n",nnm);
  if (debug)
    dump_nm2type(debug,nnm,nm2t);
  
  printf("Generating bonds from distances...\n");
  snew(nbonds,atoms->nr);
  mk_bonds(nnm,nm2t,atoms,x,&(plist[F_BONDS]),nbonds,forcefield,
	   bPBC,box,atomprop,btol);

  open_symtab(&symtab);
  atype = set_atom_type(&symtab,atoms,&(plist[F_BONDS]),nbonds,nnm,nm2t);
  
  /* Read charges */
  if (opt2bSet("-d",NFILE,fnm))
    qa = rd_q_alpha(opt2fn("-d",NFILE,fnm),&nqa);

  alg = eqgNone;
  if (qgen[0]) {
    for(alg=1; (alg<=eqgNR); alg++) {
      if (strcmp(qgen[0],qgen[alg]) == 0)
	break;
    }
    if (alg > eqgNR)
      alg = eqgNone;
    else
      alg--;
  }
  assign_charge_alpha(alg,atoms,x,nqa,qa,&(plist[F_BONDS]),qtol,fac,
		      maxiter,atomprop,qtotref,muref);

  /* Make Angles and Dihedrals */
  snew(excls,atoms->nr);
  printf("Generating angles and dihedrals from bonds...\n");
  init_nnb(&nnb,atoms->nr,5);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,nexcl,plist,excls,NULL,bAllDih,bRemoveDih,TRUE);
  delete_shell_interactions(plist,atoms,atype,&nnb,excls);
  done_nnb(&nnb);
  mu = calc_dip(atoms,x);
  
  if (!bPairs)
    plist[F_LJ14].nr = 0;
  fprintf(stderr,
	  "There are %4d %s dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d atoms\n",
	  plist[F_PDIHS].nr, 
	  bOPLS ? "Ryckaert-Bellemans" : "proper",
	  plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	  plist[F_LJ14].nr, plist[F_BONDS].nr,atoms->nr);

  calc_angles_dihs(&plist[F_ANGLES],&plist[F_PDIHS],x,bPBC,box);
  
  set_force_const(plist,kb,kt,kp,bRound,bParam);

  cgnr = set_cgnr(atoms,bUsePDBcharge,&qtot,&mtot);
  printf("Total charge is %g, total mass is %g, dipole is %g D\n",
	 qtot,mtot,mu);
  if (bOPLS) {
    bts[2] = 3;
    bts[3] = 1;
  }
  if (bCharmm) {
    bts[1] = 5;
    bts[2] = 3;
    bts[3] = 1;
  }
  reset_q(atoms);
  
  if (bTOP) {    
    fp = ftp2FILE(efTOP,NFILE,fnm,"w");
    print_top_header(fp,ftp2fn(efTOP,NFILE,fnm),
		     "Generated by x2top",TRUE, forcefield,1.0);
    /*switch (alg) {
    case eqgYang:
      please_cite(fp,"Yang2006b");
      break;
    case eqgBultinck:
      please_cite(fp,"Bultinck2002a");
      break;
    default:
      break;
      } */   
    write_top(fp,NULL,mymol.name,atoms,bts,plist,excls,atype,
	      cgnr,nexcl);
    print_top_mols(fp,mymol.name,NULL,0,NULL,1,&mymol);
    
    fclose(fp);
  }
  if (bRTP)
    print_rtp(ftp2fn(efRTP,NFILE,fnm),"Generated by x2top",
	      atoms,plist,cgnr,asize(bts),bts);
  
  if (debug) {
    dump_hybridization(debug,atoms,nbonds);
  }
  close_symtab(&symtab);
    
  thanx(stderr);
  
  return 0;
}
