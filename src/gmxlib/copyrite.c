/*
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_THREADS
#include <thread_mpi.h>
#endif

/* This file is completely threadsafe - keep it that way! */

#include <string.h>
#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "macros.h"
#include "time.h"
#include "random.h"
#include "statutil.h"
#include "copyrite.h"
#include "strdb.h"
#include "futil.h"

static void pr_two(FILE *out,int c,int i)
{
  if (i < 10)
    fprintf(out,"%c0%1d",c,i);
  else
    fprintf(out,"%c%2d",c,i);
}

void pr_difftime(FILE *out,double dt)
{
  int    ndays,nhours,nmins,nsecs;
  gmx_bool   bPrint,bPrinted;

  ndays = dt/(24*3600);
  dt    = dt-24*3600*ndays;
  nhours= dt/3600;
  dt    = dt-3600*nhours;
  nmins = dt/60;
  dt    = dt-nmins*60;
  nsecs = dt;
  bPrint= (ndays > 0);
  bPrinted=bPrint;
  if (bPrint) 
    fprintf(out,"%d",ndays);
  bPrint=bPrint || (nhours > 0);
  if (bPrint) {
    if (bPrinted)
      pr_two(out,'d',nhours);
    else 
      fprintf(out,"%d",nhours);
  }
  bPrinted=bPrinted || bPrint;
  bPrint=bPrint || (nmins > 0);
  if (bPrint) {
    if (bPrinted)
      pr_two(out,'h',nmins);
    else 
      fprintf(out,"%d",nmins);
  }
  bPrinted=bPrinted || bPrint;
  if (bPrinted)
    pr_two(out,':',nsecs);
  else
    fprintf(out,"%ds",nsecs);
  fprintf(out,"\n");
}


gmx_bool be_cool(void)
{
  /* Yes, it is bad to check the environment variable every call,
   * but we dont call this routine often, and it avoids using 
   * a mutex for locking the variable...
   */
#ifdef GMX_FAHCORE
  /*be uncool*/
  return FALSE;
#else
  return (getenv("GMX_NO_QUOTES") == NULL);
#endif
}

void space(FILE *out, int n)
{
  fprintf(out,"%*s",n,"");
}

static void sp_print(FILE *out,const char *s)
{
  int slen;
  
  slen=strlen(s);
  space(out,(80-slen)/2);
  fprintf(out,"%s\n",s);
}

static void ster_print(FILE *out,const char *s)
{
  int  slen;
  char buf[128];
  
  snprintf(buf,128,":-)  %s  (-:",s);
  slen=strlen(buf);
  space(out,(80-slen)/2);
  fprintf(out,"%s\n",buf);
}


static void pukeit(const char *db,const char *defstring, char *retstring, 
		   int retsize, int *cqnum)
{
  FILE *fp;
  char **help;
  int  i,nhlp;
  int  seed;
 
  if (be_cool() && ((fp = low_libopen(db,FALSE)) != NULL)) {
    nhlp=fget_lines(fp,&help);
    /* for libraries we can use the low-level close routines */
    ffclose(fp);
    seed=time(NULL);
    *cqnum=nhlp*rando(&seed);
    if (strlen(help[*cqnum]) >= STRLEN)
      help[*cqnum][STRLEN-1] = '\0';
    strncpy(retstring,help[*cqnum],retsize);
    for(i=0; (i<nhlp); i++)
      sfree(help[i]);
    sfree(help);
  }
  else 
    strncpy(retstring,defstring,retsize);
}

void bromacs(char *retstring, int retsize)
{
  int dum;

  pukeit("bromacs.dat",
	 "Groningen Machine for Chemical Simulation",
	 retstring,retsize,&dum);
}

void cool_quote(char *retstring, int retsize, int *cqnum)
{
  char *tmpstr;
  char *s,*ptr;
  int tmpcq,*p;
  
  if (cqnum!=NULL)
    p = cqnum;
  else
    p = &tmpcq;
  
  /* protect audience from explicit lyrics */
  snew(tmpstr,retsize+1);
  pukeit("gurgle.dat","Thanx for Using GROMACS - Have a Nice Day",
	 tmpstr,retsize-2,p);

  if ((ptr = strchr(tmpstr,'_')) != NULL) {
    *ptr='\0';
    ptr++;
    sprintf(retstring,"\"%s\" %s",tmpstr,ptr);
  }
  else {
    strcpy(retstring,tmpstr);
  }
  sfree(tmpstr);
}

void CopyRight(FILE *out,const char *szProgram)
{
  /* Dont change szProgram arbitrarily - it must be argv[0], i.e. the 
   * name of a file. Otherwise, we won't be able to find the library dir.
   */
#define NCR (int)asize(CopyrightText)
#ifdef GMX_FAHCORE
#define NGPL 0 /*FAH has an exception permission from GPL to allow digital signatures in Gromacs*/
#else
#define NGPL (int)asize(GPLText)
#endif

  char buf[256],tmpstr[1024];
  int i;

#ifdef GMX_FAHCORE
  set_program_name("Gromacs");
#else
  set_program_name(szProgram);
#endif

  ster_print(out,"G  R  O  M  A  C  S");
  fprintf(out,"\n");
  
  bromacs(tmpstr,1023);
  sp_print(out,tmpstr); 
  fprintf(out,"\n");

  ster_print(out,GromacsVersion());
  fprintf(out,"\n");

  /* fprintf(out,"\n");*/

  /* sp_print(out,"PLEASE NOTE: THIS IS A BETA VERSION\n");
  
  fprintf(out,"\n"); */

  for(i=0; (i<NCR); i++) 
    sp_print(out,CopyrightText[i]);
  for(i=0; (i<NGPL); i++)
    sp_print(out,GPLText[i]);

  fprintf(out,"\n");

  snprintf(buf,256,"%s",Program());
#ifdef GMX_DOUBLE
  strcat(buf," (double precision)");
#endif
  ster_print(out,buf);
  fprintf(out,"\n");
}


void thanx(FILE *fp)
{
  char cq[1024];
  int  cqnum;

  /* protect the audience from suggestive discussions */
  cool_quote(cq,1023,&cqnum);
  
  if (be_cool()) 
    fprintf(fp,"\ngcq#%d: %s\n\n",cqnum,cq);
  else
    fprintf(fp,"\n%s\n\n",cq);
}

typedef struct {
  const char *key;
  const char *author;
  const char *title;
  const char *journal;
  int volume,year;
  const char *pages;
} t_citerec;

void please_cite(FILE *fp,const char *key)
{
  static const t_citerec citedb[] = {
    { "Allen1987a",
      "M. P. Allen and D. J. Tildesley",
      "Computer simulation of liquids",
      "Oxford Science Publications",
      1, 1987, "1" },
    { "Berendsen95a",
      "H. J. C. Berendsen, D. van der Spoel and R. van Drunen",
      "GROMACS: A message-passing parallel molecular dynamics implementation",
      "Comp. Phys. Comm.",
      91, 1995, "43-56" },
    { "Berendsen84a",
      "H. J. C. Berendsen, J. P. M. Postma, A. DiNola and J. R. Haak",
      "Molecular dynamics with coupling to an external bath",
      "J. Chem. Phys.",
      81, 1984, "3684-3690" },
    { "Ryckaert77a",
      "J. P. Ryckaert and G. Ciccotti and H. J. C. Berendsen",
      "Numerical Integration of the Cartesian Equations of Motion of a System with Constraints; Molecular Dynamics of n-Alkanes",
      "J. Comp. Phys.",
      23, 1977, "327-341" },
    { "Miyamoto92a",
      "S. Miyamoto and P. A. Kollman",
      "SETTLE: An Analytical Version of the SHAKE and RATTLE Algorithms for Rigid Water Models",
      "J. Comp. Chem.",
      13, 1992, "952-962" },
    { "Cromer1968a",
      "D. T. Cromer & J. B. Mann",
      "X-ray scattering factors computed from numerical Hartree-Fock wave functions",
      "Acta Cryst. A",
      24, 1968, "321" },
    { "Barth95a",
      "E. Barth and K. Kuczera and B. Leimkuhler and R. D. Skeel",
      "Algorithms for Constrained Molecular Dynamics",
      "J. Comp. Chem.",
      16, 1995, "1192-1209" },
    { "Essmann95a",
      "U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee and L. G. Pedersen ",
      "A smooth particle mesh Ewald method",
      "J. Chem. Phys.",
      103, 1995, "8577-8592" },
    { "Torda89a",
      "A. E. Torda and R. M. Scheek and W. F. van Gunsteren",
      "Time-dependent distance restraints in molecular dynamics simulations",
      "Chem. Phys. Lett.",
      157, 1989, "289-294" },
    { "Tironi95a",
      "I. G. Tironi and R. Sperb and P. E. Smith and W. F. van Gunsteren",
      "Generalized reaction field method for molecular dynamics simulations",
      "J. Chem. Phys",
      102, 1995, "5451-5459" },
    { "Hess97a",
      "B. Hess and H. Bekker and H. J. C. Berendsen and J. G. E. M. Fraaije",
      "LINCS: A Linear Constraint Solver for molecular simulations",
      "J. Comp. Chem.",
      18, 1997, "1463-1472" },
    { "Hess2008a",
      "B. Hess",
      "P-LINCS: A Parallel Linear Constraint Solver for molecular simulation",
      "J. Chem. Theory Comput.",
      4, 2008, "116-122" },
    { "Hess2008b",
      "B. Hess and C. Kutzner and D. van der Spoel and E. Lindahl",
      "GROMACS 4: Algorithms for highly efficient, load-balanced, and scalable molecular simulation",
      "J. Chem. Theory Comput.",
      4, 2008, "435-447" },
    { "Hub2010",
      "J. S. Hub, B. L. de Groot and D. van der Spoel",
      "g_wham - A free weighted histogram analysis implementation including robust error and autocorrelation estimates",
      "J. Chem. Theory Comput.",
      6, 2010, "3713-3720"}, 
    { "In-Chul99a",
      "Y. In-Chul and M. L. Berkowitz",
      "Ewald summation for systems with slab geometry",
      "J. Chem. Phys.",
      111, 1999, "3155-3162" },
    { "DeGroot97a",
      "B. L. de Groot and D. M. F. van Aalten and R. M. Scheek and A. Amadei and G. Vriend and H. J. C. Berendsen",
      "Prediction of Protein Conformational Freedom From Distance Constrains",
      "Proteins",
      29, 1997, "240-251" },
    { "Spoel98a",
      "D. van der Spoel and P. J. van Maaren and H. J. C. Berendsen",
      "A systematic study of water models for molecular simulation. Derivation of models optimized for use with a reaction-field.",
      "J. Chem. Phys.",
      108, 1998, "10220-10230" },
    { "Wishart98a",
      "D. S. Wishart and A. M. Nip",
      "Protein Chemical Shift Analysis: A Practical Guide",
      "Biochem. Cell Biol.",
      76, 1998, "153-163" },
    { "Maiorov95",
      "V. N. Maiorov and G. M. Crippen",
      "Size-Independent Comparison of Protein Three-Dimensional Structures",
      "PROTEINS: Struct. Funct. Gen.",
      22, 1995, "273-283" },
    { "Feenstra99",
      "K. A. Feenstra and B. Hess and H. J. C. Berendsen",
      "Improving Efficiency of Large Time-scale Molecular Dynamics Simulations of Hydrogen-rich Systems",
      "J. Comput. Chem.",
      20, 1999, "786-798" },
    { "Timneanu2004a",
      "N. Timneanu and C. Caleman and J. Hajdu and D. van der Spoel",
      "Auger Electron Cascades in Water and Ice",
      "Chem. Phys.",
      299, 2004, "277-283" },
    { "Pascal2011a",
      "T. A. Pascal and S. T. Lin and W. A. Goddard III",
      "Thermodynamics of liquids: standard molar entropies and heat capacities of common solvents from 2PT molecular dynamics",
      "Phys. Chem. Chem. Phys.",
      13, 2011, "169-181" },
    { "Caleman2011b",
      "C. Caleman and M. Hong and J. S. Hub and L. T. da Costa and P. J. van Maaren and D. van der Spoel",
      "Force Field Benchmark 1: Density, Heat of Vaporization, Heat Capacity, Surface Tension and Dielectric Constant of 152 Organic Liquids",
      "Submitted",
      0, 2011, "" },
    { "Lindahl2001a",
      "E. Lindahl and B. Hess and D. van der Spoel",
      "GROMACS 3.0: A package for molecular simulation and trajectory analysis",
      "J. Mol. Mod.",
      7, 2001, "306-317" },
    { "Wang2001a",
      "J. Wang and W. Wang and S. Huo and M. Lee and P. A. Kollman",
      "Solvation model based on weighted solvent accessible surface area",
      "J. Phys. Chem. B",
      105, 2001, "5055-5067" },
    { "Eisenberg86a",
      "D. Eisenberg and A. D. McLachlan",
      "Solvation energy in protein folding and binding",
      "Nature",
      319, 1986, "199-203" },
    { "Eisenhaber95",
      "Frank Eisenhaber and Philip Lijnzaad and Patrick Argos and Chris Sander and Michael Scharf",
      "The Double Cube Lattice Method: Efficient Approaches to Numerical Integration of Surface Area and Volume and to Dot Surface Contouring of Molecular Assemblies",
      "J. Comp. Chem.",
      16, 1995, "273-284" },
    { "Hess2002",
      "B. Hess, H. Saint-Martin and H.J.C. Berendsen",
      "Flexible constraints: an adiabatic treatment of quantum degrees of freedom, with application to the flexible and polarizable MCDHO model for water",
      "J. Chem. Phys.",
      116, 2002, "9602-9610" },
    { "Hetenyi2002b",
      "Csaba Hetenyi and David van der Spoel",
      "Efficient docking of peptides to proteins without prior knowledge of the binding site.",
      "Prot. Sci.",
      11, 2002, "1729-1737" },
    { "Hess2003",
      "B. Hess and R.M. Scheek",
      "Orientation restraints in molecular dynamics simulations using time and ensemble averaging",
      "J. Magn. Res.",
      164, 2003, "19-27" },
    { "Rappe1991a",
      "A. K. Rappe and W. A. Goddard III",
      "Charge Equillibration for Molecular Dynamics Simulations",
      "J. Phys. Chem.",
      95, 1991, "3358-3363" },
    { "Mu2005a",
      "Y. Mu, P. H. Nguyen and G. Stock",
      "Energy landscape of a small peptide revelaed by dihedral angle principal component analysis",
      "Prot. Struct. Funct. Bioinf.",
      58, 2005, "45-52" },
    { "Okabe2001a",
      "T. Okabe and M. Kawata and Y. Okamoto and M. Mikami",
      "Replica-exchange {M}onte {C}arlo method for the isobaric-isothermal ensemble",
      "Chem. Phys. Lett.",
      335, 2001, "435-439" },
    { "Hukushima96a",
      "K. Hukushima and K. Nemoto",
      "Exchange Monte Carlo Method and Application to Spin Glass Simulations",
      "J. Phys. Soc. Jpn.",
      65, 1996, "1604-1608" },
    { "Tropp80a",
      "J. Tropp",
      "Dipolar Relaxation and Nuclear Overhauser effects in nonrigid molecules: The effect of fluctuating internuclear distances",
      "J. Chem. Phys.",
      72, 1980, "6035-6043" },
    { "Bultinck2002a",
       "P. Bultinck and W. Langenaeker and P. Lahorte and F. De Proft and P. Geerlings and M. Waroquier and J. P. Tollenaere",
      "The electronegativity equalization method I: Parametrization and validation for atomic charge calculations",
      "J. Phys. Chem. A",
      106, 2002, "7887-7894" },
    { "Yang2006b",
      "Q. Y. Yang and K. A. Sharp",
      "Atomic charge parameters for the finite difference Poisson-Boltzmann method using electronegativity neutralization",
      "J. Chem. Theory Comput.",
      2, 2006, "1152-1167" },
    { "Spoel2005a",
      "D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A. E. Mark and H. J. C. Berendsen",
      "GROMACS: Fast, Flexible and Free",
      "J. Comp. Chem.",
      26, 2005, "1701-1719" },
    { "Spoel2006b",
      "D. van der Spoel, P. J. van Maaren, P. Larsson and N. Timneanu",
      "Thermodynamics of hydrogen bonding in hydrophilic and hydrophobic media",
      "J. Phys. Chem. B",
      110, 2006, "4393-4398" },
    { "Spoel2006d",
      "D. van der Spoel and M. M. Seibert",
      "Protein folding kinetics and thermodynamics from atomistic simulations",
      "Phys. Rev. Letters",
      96, 2006, "238102" },
    { "Palmer94a",
      "B. J. Palmer",
      "Transverse-current autocorrelation-function calculations of the shear viscosity for molecular liquids",
      "Phys. Rev. E",
      49, 1994, "359-366" },
    { "Bussi2007a",
      "G. Bussi, D. Donadio and M. Parrinello",
      "Canonical sampling through velocity rescaling",
      "J. Chem. Phys.",
      126, 2007, "014101" },
    { "Hub2006",
      "J. S. Hub and B. L. de Groot",
      "Does CO2 permeate through Aquaporin-1?",
      "Biophys. J.",
      91, 2006, "842-848" },
    { "Hub2008",
      "J. S. Hub and B. L. de Groot",
      "Mechanism of selectivity in aquaporins and aquaglyceroporins",
      "PNAS",
      105, 2008, "1198-1203" },
    { "Friedrich2009",
      "M. S. Friedrichs, P. Eastman, V. Vaidyanathan, M. Houston, S. LeGrand, A. L. Beberg, D. L. Ensign, C. M. Bruns, and V. S. Pande",
      "Accelerating Molecular Dynamic Simulation on Graphics Processing Units",
      "J. Comp. Chem.",
      30, 2009, "864-872" },
    { "Engin2010",
      "O. Engin, A. Villa, M. Sayar and B. Hess",
      "Driving Forces for Adsorption of Amphiphilic Peptides to Air-Water Interface",
      "J. Phys. Chem. B",
      0, 2010, "???" },
    { "Wang2010",
      "H. Wang, F. Dommert, C.Holm",
      "Optimizing working parameters of the smooth particle mesh Ewald algorithm in terms of accuracy and efficiency",
      "J. Chem. Phys. B",
      133, 2010, "034117"
    },
    { "Sugita1999a",
      "Y. Sugita, Y. Okamoto",
      "Replica-exchange molecular dynamics method for protein folding",
      "Chem. Phys. Lett.",
      314, 1999, "141-151" },
  };
#define NSTR (int)asize(citedb)
  
  int  j,index;
  char *author;
  char *title;
#define LINE_WIDTH 79
  
  if (fp == NULL)
    return;

  for(index=0; (index<NSTR) && (strcmp(citedb[index].key,key) != 0); index++)
    ;
  
  fprintf(fp,"\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
  if (index < NSTR) {
    /* Insert newlines */
    author = wrap_lines(citedb[index].author,LINE_WIDTH,0,FALSE);
    title  = wrap_lines(citedb[index].title,LINE_WIDTH,0,FALSE);
    fprintf(fp,"%s\n%s\n%s %d (%d) pp. %s\n",
	    author,title,citedb[index].journal,
	    citedb[index].volume,citedb[index].year,
	    citedb[index].pages);
    sfree(author);
    sfree(title);
  }
  else {
    fprintf(fp,"Entry %s not found in citation database\n",key);
  }
  fprintf(fp,"-------- -------- --- Thank You --- -------- --------\n\n");
  fflush(fp);
}

#ifdef USE_VERSION_H
/* Version information generated at compile time. */
#include "version.h"
#else
/* Fall back to statically defined version. */
static const char _gmx_ver_string[]="VERSION " VERSION;
#endif

/* This routine only returns a static (constant) string, so we use a 
 * mutex to initialize it. Since the string is only written to the
 * first time, there is no risk with multiple calls overwriting the
 * output for each other.
 */
const char *GromacsVersion()
{
  return _gmx_ver_string;
}


void gmx_print_version_info(FILE *fp)
{
    fprintf(fp, "Version:          %s\n", _gmx_ver_string);
#ifdef USE_VERSION_H
    fprintf(fp, "GIT SHA1 hash:    %s\n", _gmx_full_git_hash);
    /* Only print out the branch information if present.
     * The generating script checks whether the branch point actually
     * coincides with the hash reported above, and produces an empty string
     * in such cases. */
    if (_gmx_central_base_hash[0] != 0)
    {
        fprintf(fp, "Branched from:    %s\n", _gmx_central_base_hash);
    }
#endif

#ifdef GMX_DOUBLE
    fprintf(fp, "Precision:        double\n");
#else
    fprintf(fp, "Precision:        single\n");
#endif

#ifdef GMX_THREADS
    fprintf(fp, "Parallellization: thread_mpi\n");
#elif defined(GMX_MPI)
    fprintf(fp, "Parallellization: MPI\n");
#else
    fprintf(fp, "Parallellization: none\n");
#endif

#ifdef GMX_FFT_FFTPACK
    fprintf(fp, "FFT Library:      fftpack\n");
#elif defined(GMX_FFT_FFTW2)
    fprintf(fp, "FFT Library:      fftw2\n");
#elif defined(GMX_FFT_FFTW3)
    fprintf(fp, "FFT Library:      fftw3\n");
#elif defined(GMX_FFT_MKL)
    fprintf(fp, "FFT Library:      MKL\n");
#else
    fprintf(fp, "FFT Library:      unknown\n");
#endif
}
