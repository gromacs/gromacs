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
static char *SRCID_g_analyze_c = "$Id$";

#include <math.h>
#include <string.h>
#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "txtdump.h"
#include "enxio.h"
#include "gstat.h"
#include "xvgr.h"

typedef struct {
  int nlj,nqq;
  int *lj;
  int *qq;
} t_liedata;

static t_liedata *analyze_names(int nre,char *names[],char *ligand)
{
  int       i;
  t_liedata *ld;
  
  /* Skip until we come to pressure */
  for(i=0; (i<F_NRE); i++)
    if (strcmp(names[i],interaction_function[F_PRES].longname) == 0)
      break;
      
  /* Now real analysis: find components of energies */
  snew(ld,1);
  for( ; (i<nre); i++) {
    if (strstr(names[i],ligand) != NULL) {
      if (strstr(names[i],"LJ") != NULL) {
	ld->nlj++;
	srenew(ld->lj,ld->nlj);
	ld->lj[ld->nlj-1] = i;
      }
      else if (strstr(names[i],"Coul") != NULL) {
	ld->nqq++;
	srenew(ld->qq,ld->nqq);
	ld->qq[ld->nqq-1] = i;
      }
    }
  }
  printf("Using the following energy terms:\n");
  printf("LJ:  ");
  for(i=0; (i<ld->nlj); i++)
    printf("  %12s",names[ld->lj[i]]);
  printf("\nCoul:");
  for(i=0; (i<ld->nqq); i++)
    printf("  %12s",names[ld->qq[i]]);
  printf("\n");
  
  return ld;
}

real calc_lie(t_liedata *ld,t_energy ee[],real lie_lj,real lie_qq,
	      real fac_lj,real fac_qq)
{
  int i;
  real lj_tot,qq_tot;
  
  lj_tot = 0;
  for(i=0; (i<ld->nlj); i++)
    lj_tot += ee[ld->lj[i]].e;
  qq_tot = 0;
  for(i=0; (i<ld->nqq); i++)
    qq_tot += ee[ld->qq[i]].e;
    
  /* And now the great LIE formula: */
  return fac_lj*(lj_tot-lie_lj)+fac_qq*(qq_tot-lie_qq);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_lie computes a free energy estimate based on an energy analysis",
    "from. One needs an energy file with the following components:",
    "Coul (A-B) LJ-SR (A-B) etc."
  };
  static real lie_lj=0,lie_qq=0,fac_lj=0.181,fac_qq=0.4;
  static char *ligand=NULL;
  t_pargs pa[] = {
    { "-Elj",  FALSE, etREAL, {&lie_lj},
      "Lennard-Jones interaction between ligand and solvent" },
    { "-Eqq",  FALSE, etREAL, {&lie_qq},
      "Coulomb interaction between ligand and solvent" },
    { "-Clj",  FALSE, etREAL, {&fac_lj},
      "Factor in the LIE equation for Lennard-Jones component of energy" },
    { "-Cqq",  FALSE, etREAL, {&lie_qq},
      "Factor in the LIE equation for Coulomb component of energy" },
    { "-ligand",  FALSE, etSTR, {&ligand},
      "Name of the ligand in the energy file" }
  };
#define NPA asize(pa)

  FILE      *out;
  int       fp,nre,ndr,step;
  char      **enm=NULL;
  bool      bCont;
  t_liedata *ld;
  t_energy  *ee;
  real      t,lie;
  
  t_filenm fnm[] = { 
    { efENX, "-f",    "ener",     ffREAD   },
    { efXVG, "-o",    "lie",      ffWRITE  }
  }; 
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL); 
  
  fp = open_enx(ftp2fn(efENX,NFILE,fnm),"r");
  do_enxnms(fp,&nre,&enm);
  
  ld = analyze_names(nre,enm,ligand);
  snew(ee,nre);
  out = xvgropen(ftp2fn(efXVG,NFILE,fnm),"LIE free energy estimate","Time (ps)","DGbind (kJ/mol)");
  do {
    bCont = do_enx(fp,&t,&step,&nre,ee,&ndr,NULL);
    lie   = calc_lie(ld,ee,lie_lj,lie_qq,fac_lj,fac_qq);
    fprintf(out,"%10g  %10g\n",t,lie);
  } while (bCont);
  close_enx(fp);
  fclose(out);

  do_view(ftp2fn(efXVG,NFILE,fnm),NULL);
    
  thanx(stderr);

  return 0;
}
  
