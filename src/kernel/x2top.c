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
#include "pdbio.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "vec.h"
#include "gmx_random.h"
#include "gmx_elements.h"
#include "x2top_eemprops.h"
#include "x2top_nm2type.h"
#include "x2top_qalpha.h"
#include "x2top_core.h"
#include "x2top_qgen.h"
#include "atomprop.h"
#include "grompp.h"
#include "add_par.h"

enum { edihNo, edihOne, edihAll, edihNR };

static int get_option(char **opts)
{
  int val = 0;
  
  if (!opts)
    return NOTSET;
  if (opts[val] != NULL)
    for(val=1; (opts[val] != NULL); val++)
      if (strcasecmp(opts[0],opts[val]) == 0)
	break;
  if (opts[val] == NULL)
    val = 0;
  else
    val--;
    
  return val;
}

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
    "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
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
  };
  FILE       *fp;
  t_params   plist[F_NRE];
  t_excls    *excls;
  t_atoms    *atoms;       /* list with all atoms */
  t_atomtype atype;
  t_nextnb   nnb;
  x2top_nm2t nm2t;
  x2top_qat  qa;
  t_mols     mymol;
  void       *atomprop;
  char       title[STRLEN],forcefield[32];
  rvec       *x;        /* coordinates? */
  int        *nbonds,*cgnr;
  int        bts[] = { 1,1,1,2 };
  int        ePBC;
  matrix     box;          /* box length matrix */
  int        natoms;       /* number of atoms in one molecule  */
  int        nres;         /* number of molecules? */
  int        i,j,k,l,m,ndih,alg,dih,cgtp;
  real       mu;
  bool       bRTP,bTOP,bOPLS,bCharmm;
  t_symtab   symtab;
  int        nqa=0;
  real       cutoff,qtot,mtot,hardness=1;
  char       *fn,rtp[STRLEN];
  gmx_conect gc;

  t_filenm fnm[] = {
    { efPDB, "-f", "conf", ffREAD  },
    { efTOP, "-o", "out",  ffWRITE },
    { efSTO, "-c", "out",  ffWRITE },
    { efRTP, "-r", "out",  ffOPTWR },
    { efDAT, "-d", "qpol", ffOPTRD }
  };
#define NFILE asize(fnm)
  static real scale = 1.1, kb = 4e5,kt = 400,kp = 5;
  static real btol=0.1,qtol=1e-3,fac=5.0;
  static real qtotref=0;
  static int  nexcl = 2;
  static int  maxiter=100;
  static bool bRemoveDih = FALSE,bQsym = TRUE;
  static bool bParam = TRUE, bH14 = TRUE,bRound = TRUE;
  static bool bPairs = TRUE, bPBC = TRUE;
  static bool bUsePDBcharge = FALSE,bVerbose=FALSE,bCONECT=TRUE;
  static char *molnm = "ICE";
  static char *ff = "select";
  static char *qgen[] = { NULL, "None", "Yang", "Bultinck", "Rappe", 
			  "SMp", "SMpp", "SMs", "SMps", "SMg", "SMpg", NULL };
  static char *dihopt[] = { NULL, "No", "Single", "All", NULL };
  static char *cgopt[] = { NULL, "Group", "Atom", "Neutral", NULL };
  t_pargs pa[] = {
    { "-ff",     FALSE, etSTR, {&ff},
      "Select the force field for your simulation." },
    { "-v",      FALSE, etBOOL, {&bVerbose},
      "Generate verbose output in the top file." },
    { "-nexcl", FALSE, etINT,  {&nexcl},
      "Number of exclusions" },
    { "-H14",    FALSE, etBOOL, {&bH14}, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-dih",    FALSE, etSTR,  {dihopt}, 
      "Which proper dihedrals to generate: none, one per rotatable bond, or all possible." },
    { "-remdih", FALSE, etBOOL, {&bRemoveDih}, 
      "Remove dihedrals on the same bond as an improper" },
    { "-pairs",  FALSE, etBOOL, {&bPairs},
      "Output 1-4 interactions (pairs) in topology file" },
    { "-name",   FALSE, etSTR,  {&molnm},
      "Name of your molecule" },
    { "-pbc",    FALSE, etBOOL, {&bPBC},
      "Use periodic boundary conditions." },
    { "-conect", FALSE, etBOOL, {&bCONECT},
      "Use CONECT records in the pdb file to signify bonds" },
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
    { "-maxiter",FALSE, etINT, {&maxiter},
      "HIDDENMax number of iterations for charge generation algorithm" },
    { "-fac",    FALSE, etREAL, {&fac},
      "HIDDENNot yet understood factor for generating charges" },
    { "-qsymm",  FALSE, etBOOL, {&bQsym},
      "HIDDENSymmetrize the charges on methyl and NH3 groups." },
    { "-cgsort", FALSE, etSTR, {cgopt},
      "HIDDENOption for assembling charge groups: based on Group (default, e.g. CH3 groups are kept together), Atom, or Neutral sections" },
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
  bTOP = TRUE;
  
  if (!bRTP && !bTOP)
    gmx_fatal(FARGS,"Specify at least one output file");
  
  if ((btol < 0) || (btol > 1)) 
    gmx_fatal(FARGS,"Bond tolerance should be between 0 and 1 (not %g)",
	      btol);
  if ((qtol < 0) || (qtol > 1)) 
    gmx_fatal(FARGS,"Charge tolerance should be between 0 and 1 (not %g)",
	      qtol);
  /* Check command line options of type enum */
  dih  = get_option(dihopt);
  cgtp = get_option(cgopt);
  alg  = get_option(qgen);
    
  /* Read standard atom properties */
  atomprop = get_atomprop();
    
  /* Force field selection */
  if (!strncmp(ff,"select",6)) {
    /* Interactive forcefield selection */
    choose_ff(forcefield,sizeof(forcefield));
  } 
  else {
    sprintf(forcefield,"ff%s",ff);
  }
  if (0) {
    sprintf(rtp,"%s.rtp",forcefield);
    if ((fp = libopen(rtp)) == NULL) 
      gmx_fatal(FARGS,"Force field file %s does not exist\n",rtp);
    fclose(fp);
  }
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

  if (bCONECT)
    gc = init_gmx_conect();
  else
    gc = NULL;
  read_pdb_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,&ePBC,box,FALSE,gc);
  get_pdb_atomnumber(atoms,atomprop);
  if (bCONECT && debug)
    dump_conection(debug,gc);
  
  nm2t = rd_nm2type(forcefield,atomprop);
  if (debug) 
    dump_nm2type(debug,nm2t);
  
  if (bVerbose)  
    printf("Generating bonds from distances...\n");
  snew(nbonds,atoms->nr);
  mk_bonds(nm2t,atoms,x,gc,&(plist[F_BONDS]),nbonds,forcefield,
	   bPBC,box,atomprop,btol);

  /* Setting the atom types: this depends on the bonding */
  open_symtab(&symtab);
  atype = set_atom_type(&symtab,atoms,&(plist[F_BONDS]),
			nbonds,nm2t,atomprop);
  if (debug) 
    dump_hybridization(debug,atoms,nbonds);
  sfree(nbonds);

  /* Read charges and polarizabilities, if provided */
  if (opt2bSet("-d",NFILE,fnm))
    qa = rd_q_alpha(opt2fn("-d",NFILE,fnm));

  /* Check which algorithm to use for charge generation */
  if (alg == eqgNone) {
    if (nqa == atoms->nr) {
      /* Use values from file */
      for(i=0; (i<nqa); i++) {
	atoms->atom[i].q  = get_qa_q(*atoms->atomname[i],qa);
	atoms->atom[i].qB = get_qa_alpha(*atoms->atomname[i],qa);
      }
    }
    else
      printf("Using zero charges\n");
  }
  else {
    /* Hardness */
    if (alg  == eqgBultinck)
      hardness = 2;
    
    generate_charges(stdout,molnm,alg,atoms,x,&(plist[F_BONDS]),qtol,fac,
		     maxiter,atomprop,qtotref,hardness);
    if (bQsym)
      symmetrize_charges(atoms,atype,&(plist[F_BONDS]),atomprop);
  }
    
  /* Make Angles and Dihedrals */
  snew(excls,atoms->nr);
  if (bVerbose)
    printf("Generating angles and dihedrals from bonds...\n");
  init_nnb(&nnb,atoms->nr,nexcl);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,nexcl,plist,excls,NULL,
	  (dih == edihAll),bRemoveDih,TRUE);
  generate_excls(&nnb,nexcl,excls);
  done_nnb(&nnb);
  
  if (!bPairs)
    plist[F_LJ14].nr = 0;
  if (dih == edihNo)
    plist[F_PDIHS].nr = 0;
  
  if ((alg == eqgSMpp) || (alg == eqgSMps) ||(alg == eqgSMpg))
    add_shells(nm2t,&atoms,atype,plist,&x,&symtab,&excls);
  
  mu = calc_dip(atoms,x);
  
  calc_angles_dihs(&plist[F_ANGLES],&plist[F_PDIHS],x,bPBC,box);
  
  set_force_const(plist,kb,kt,kp,bRound,bParam);

  if ((cgnr = generate_charge_groups(cgtp,atoms,atype,
				     &plist[F_BONDS],&plist[F_POLARIZATION],
				     bUsePDBcharge,&qtot,&mtot)) == NULL)
    gmx_fatal(FARGS,"Error generating charge groups");
  sort_on_charge_groups(cgnr,atoms,plist,atype,x);
         
  if (bVerbose) {
    printf("There are %4d %s dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d atoms\n"
	   "          %4d polarizations\n",
	   plist[F_PDIHS].nr, 
	   bOPLS ? "Ryckaert-Bellemans" : "proper",
	   plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	   plist[F_LJ14].nr, plist[F_BONDS].nr,atoms->nr,
	   plist[F_POLARIZATION].nr);
    
    printf("Total charge is %g, total mass is %g, dipole is %g D\n",
	   qtot,mtot,mu);
  }
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
    /* Write topology file */
    fn = ftp2fn(efTOP,NFILE,fnm);
    fp = ffopen(fn,"w");
    print_top_header(fp,fn,"Generated by x2top",TRUE, forcefield,1.0);
    write_top(fp,NULL,mymol.name,atoms,bts,plist,excls,atype,cgnr,nexcl);
    print_top_mols(fp,mymol.name,NULL,0,NULL,1,&mymol);
    fclose(fp);
  }
  if (bRTP) {
    /* Write force field component */
    snew(atoms->atomtype,atoms->nr);
    for(i=0; (i<atoms->nr); i++)
      atoms->atomtype[i] = put_symtab(&symtab,
				      get_atomtype_name(atoms->atom[i].type,atype));
    print_rtp(ftp2fn(efRTP,NFILE,fnm),"Generated by x2top",
	      atoms,plist,cgnr,asize(bts),bts);
  }
  /* Write coordinates */ 
  sprintf(title,"%s processed by %s",molnm,ShortProgram());
  write_sto_conf(opt2fn("-c",NFILE,fnm),title,atoms,x,NULL,ePBC,box);

  close_symtab(&symtab);
    
  thanx(stderr);
  
  return 0;
}
