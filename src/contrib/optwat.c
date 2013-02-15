/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "macros.h"
#include "txtdump.h"
#include "tpxio.h"
#include "enxio.h"
#include "names.h"
#include "statutil.h"
#include "copyrite.h"
#include "random.h"

real ener(matrix P,real e,real e0,int nmol,real kp,real ke,gmx_bool bPScal)
{
  if (bPScal)
    return (kp*(sqr(P[XX][XX]+P[YY][YY]+P[ZZ][ZZ]-3))+
	    ke*sqr(e/nmol-e0));
  else
    return (kp*(sqr(P[XX][XX]-1)+sqr(P[YY][YY]-1)+sqr(P[ZZ][ZZ]-1)+
		sqr(P[XX][YY])+sqr(P[XX][ZZ])+sqr(P[YY][ZZ])) +
	    ke*sqr(e/nmol-e0));
}

void do_sim(char *enx,
	    t_topology *top,rvec *x,rvec *v,t_inputrec *ir,matrix box)
{
  char *tpx = "optwat.tpr";
  char buf[128];
  
  write_tpx(tpx,0,0.0,0.0,ir,box,top->atoms.nr,x,v,NULL,top);
  sprintf(buf,"xmdrun -s %s -e %s >& /dev/null",tpx,enx);
  system(buf);
}

void get_results(char *enx,real P[],real *epot,int pindex,int eindex)
{
  ener_file_t fp_ene;
  int      nre,step,ndr,i;
  real     t;
  gmx_enxnm_t *nms=NULL;
  t_enxframe fr;
  
  fp_ene = open_enx(enx,"r");
  
  do_enxnms(fp_ene,&nre,&nms);
  
  /* Read until the last frame */
  while (do_enx(fp_ene,&fr));
  
  close_enx(fp_ene);
  
  *epot = fr.ener[eindex].e;
  for(i=pindex; (i<pindex+9); i++)
    P[i-pindex] = fr.ener[i].e;
    
  sfree(ener);
}

void copy_iparams(int nr,t_iparams dest[],t_iparams src[])
{
  memcpy(dest,src,nr*sizeof(dest[0]));
}

void rand_step(FILE *fp,int nr,t_iparams ip[],int *seed,real frac)
{
  int  i;
  real ff;
  
  do {
    i = (int) (rando(seed)*nr);
  } while (ip[i].lj.c12 == 0.0);

  do {
    ff = frac*rando(seed);
  } while (ff == 0.0);
  
  if (rando(seed) > 0.5) {
    ip[i].lj.c12 *= 1.0+ff;
    fprintf(fp,"Increasing c12[%d] by %g%% to %g\n",i,100*ff,ip[i].lj.c12);
  }
  else {
    ip[i].lj.c12 *= 1.0-ff;
    fprintf(fp,"Decreasing c12[%d] by %g%% to %g\n",i,100*ff,ip[i].lj.c12);
  }
}

void pr_progress(FILE *fp,int nit,tensor P,real epot,real eFF,
		 double mc_crit,gmx_bool bConv,gmx_bool bAccept)
{
  fprintf(fp,"Iter %3d, eFF = %g, Converged = %s,  Accepted = %s\n",
	  nit,eFF,yesno_names[bConv],yesno_names[bAccept]);
  fprintf(fp,"Epot = %g  Pscal = %g,  mc_crit = %g\n",epot,
	  trace(P)/3,mc_crit);
  pr_rvecs(fp,0,"Pres",P,DIM);
  fprintf(fp,"-----------------------------------------------------\n");
  fflush(fp);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "[TT]optwat[tt] optimizes the force field parameter set of a molecular crystal",
    "to reproduce the pressure tensor and experimental energy.[PAR]",
    "Note that for good results the [TT].tpx[tt] file must contain input for a",
    "simulated annealing run, or a single point energy calculation at 0 K"
  };
  t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efEDR, "-e",      NULL,       ffRW },
    { efLOG, "-g",      NULL,       ffWRITE }
  };
#define NFILE asize(fnm)

  static real epot0    = -57, tol    = 1,   kT     = 0.0;
  static real kp       = 1,   ke     = 100, frac   = 0.1;
  static int  maxnit   = 100, eindex = 5,   pindex = 19;
  static int  seed     = 1993;
  static gmx_bool bPScal   = FALSE;
  static t_pargs pa[] = {
    { "-epot0",  FALSE, etREAL, {&epot0},
      "Potential energy in kJ/mol" },
    { "-tol",    FALSE, etREAL, {&tol},
      "Tolerance for converging" },
    { "-nit",    FALSE, etINT,  {&maxnit},
      "Max number of iterations" },
    { "-seed",   FALSE, etINT,  {&seed},
      "Random seed for MC steps" },
    { "-frac",   FALSE, etREAL, {&frac},
      "Maximum fraction by which to change parameters. Actual fraction is random between 0 and this parameter" },
    { "-pindex", FALSE, etINT,  {&pindex},
      "Index of P[X][X] in the energy file (check with [TT]g_energy[tt] and subtract 1)" },
    { "-eindex", FALSE, etINT,  {&pindex},
      "Index of Epot in the energy file (check with [TT]g_energy[tt] and subtract 1)" },
    { "-kp",     FALSE, etREAL, {&kp},
      "Force constant for pressure components"},
    { "-ke",     FALSE, etREAL, {&ke},
      "Force constant for energy component"},
    { "-kT",     FALSE, etREAL, {&kT},
      "Boltzmann Energy for Monte Carlo" },
    { "-pscal",  FALSE, etBOOL, {&bPScal},
      "Optimize params for scalar pressure, instead of tensor" }
  };

  FILE        *fp;  
  t_topology  top;
  t_tpxheader sh;
  t_inputrec  ir;
  t_iparams   *ip[2];
  int         cur=0;
#define next  (1-cur)
  rvec        *xx,*vv;
  matrix      box;
  int         i,step,natoms,nmol,nit,atnr2;
  real        t,lambda,epot,eFF[2];
  double      mc_crit=0;
  gmx_bool        bConverged,bAccept;
  tensor      P;
  
  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  /* Read initial topology and coordaintes etc. */
  read_tpxheader(ftp2fn(efTPX,NFILE,fnm),&sh,TRUE,NULL,NULL);
  snew(xx,sh.natoms);
  snew(vv,sh.natoms);
  read_tpx(ftp2fn(efTPX,NFILE,fnm),&step,&t,&lambda,&ir,box,&natoms,
	   xx,vv,NULL,&top);

  /* Open log file and print options */   
  fp    = ftp2FILE(efLOG,NFILE,fnm,"w");
  fprintf(fp,"%s started with the following parameters\n",argv[0]);
  fprintf(fp,"epot   = %8g  ke     = %8g  kp     = %8g\n",epot0,ke,kp);
  fprintf(fp,"maxnit = %8d  tol    = %8g  seed   = %8d\n",maxnit,tol,seed);
  fprintf(fp,"frac   = %8g  pindex = %8d  eindex = %8d\n",frac,pindex,eindex);
  fprintf(fp,"kT     = %8g  pscal  = %8s\n",kT,gmx_bool_names[bPScal]);
  
  /* Unpack some topology numbers */
  nmol  = top.blocks[ebMOLS].nr;
  atnr2 = top.idef.atnr*top.idef.atnr;
  
  /* Copy input params */
  snew(ip[cur],atnr2);
  snew(ip[next],atnr2);
  copy_iparams(atnr2,ip[cur],top.idef.iparams);
  copy_iparams(atnr2,ip[next],top.idef.iparams);
  
  /* Loop over iterations */
  nit = 0;
  do {
    if (nit > 0) {
      /* Do random step */
      rand_step(fp,atnr2,ip[next],&seed,frac);
      copy_iparams(atnr2,top.idef.iparams,ip[next]);
    }
    do_sim(ftp2fn(efEDR,NFILE,fnm),&top,xx,vv,&ir,box);

    get_results(ftp2fn(efEDR,NFILE,fnm),P[0],&epot,pindex,eindex);

    /* Calculate penalty */
    eFF[(nit > 0) ? next : cur]  = ener(P,epot,epot0,nmol,kp,ke,bPScal);
    
    bConverged = (eFF[(nit > 0) ? next : cur] < tol);
    
    if (nit > 0) {
      /* Do Metropolis criterium */
      if (kT > 0)
	mc_crit = exp(-(eFF[next]-eFF[cur])/kT);
      bAccept = ((eFF[next] < eFF[cur]) || 
		 ((kT > 0) && (mc_crit > rando(&seed))));
      pr_progress(fp,nit,P,epot/nmol,eFF[next],mc_crit,
		  bConverged,bAccept);
      if (bAccept) {
	/* Better params! */
	cur = next;
      }
      else {
	/* Restore old parameters */
	copy_iparams(atnr2,ip[next],ip[cur]);
      }
    }
    else
      pr_progress(fp,nit,P,epot/nmol,eFF[cur],mc_crit,bConverged,FALSE);
		
    nit++;
  } while (!bConverged && (nit < maxnit));

  for(i=0; (i<atnr2); i++)
    pr_iparams(fp,F_LJ,&ip[cur][i]);
  
  ffclose(fp);
    
  thanx(stderr);
  
  return 0;
}
