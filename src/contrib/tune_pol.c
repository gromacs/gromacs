/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maths.h"
#include "futil.h"
#include "smalloc.h"
#include "tune_pol.h"
#include "vec.h"
#include "gentop_matrix.h"

static int tabnum = 0;
static int toler  = 5;  /* Tolerance in % */
static int tolerb = 50; /* Tolerance in % */
static int my_order[eelemNR] = { 
  eelemC, 
  eelemH,  eelemHe, 
  eelemLi, eelemBe, eelemB,  eelemN,  eelemO, eelemF,  eelemNe,
  eelemNa, eelemMg, eelemAl, eelemSi, eelemP, eelemS, eelemCl, eelemAr,
  eelemK,  eelemCa, eelemGe, eelemSe, eelemBr, eelemKr, eelemRb,
  eelemI,  eelemXe, eelemCs
};

/*static char *sbasis[eqmNR] = {
  "6-31G","pVDZ","pVTZ","pVTZ","d-pVTZ", "Ahc", "Ahp", "BS", "FP" 
  };*/

static char *sbasis[eqmNR] = {
  "6-31G","aug-DZ","aug-TZ","aug-TZ","d-aug-TZ", "Ahc", "Ahp", "BS", "FP" 
};

static void check_formula(int np,t_molprop pd[])
{
  int  i,j,jj;
  char formula[1280],number[32];
  
  for(i=0; (i<np); i++) {
    formula[0] = '\0';
    for(j=0; (j<eelemNR); j++) {
      jj = my_order[j];
      if (pd[i].elem_comp[jj] > 0) {
	strcat(formula,bosque[jj].name);
	if (pd[i].elem_comp[jj] > 1) {
	  sprintf(number,"%d",pd[i].elem_comp[jj]);
	  strcat(formula,number);
	}
      }
    }
    if (0 && (strcasecmp(pd[i].formula,formula) != 0) )
      fprintf(stderr,"Formula '%s' does match '%s' based on composition for %s.\n",
	      pd[i].formula,formula,pd[i].molname);
  }
}
  
static void table_number(FILE *fp)
{
  fprintf(fp,"\\setcounter{table}{%d}\n",tabnum);
}

static void table1_header(FILE *fp,int caption)
{
  table_number(fp);
  fprintf(fp,"\\begin{table}[H]\n\\centering\n");
  if (caption)
    fprintf(fp,"\\caption{Definition of fragment-based atom types. The naming is such that the first digit represents the hybridization state and the second digit the number of hydrogen atoms attached. The number of occurences  of the fragments N is given along with the fragment polarizabilities (FP). The columns Ahc and Ahp contain group polarizabilites computed using Miller's equation~\\protect\\cite{Miller1979a} and parameters~\\protect\\cite{Miller90a} and Kang and Jhon's method~\\cite{Kang1982a} with the parametrization of Miller~\\protect\\cite{Miller90a} respectively. The column BS contains the equivalent using the polarizabilities of Bosque and Sales~\\protect\\cite{Bosque2002a}.}\n\\label{fragments}\n");
  else
    fprintf(fp,"\\caption{Continued}");
  fprintf(fp,"\\begin{tabular}{cccccccc}\n\\hline\n");
  fprintf(fp,"Group & Name  & Hybridiz. & \\multicolumn{1}{c}{N} & \\multicolumn{4}{c}{Polarizability} \\\\\n");
  fprintf(fp,"& & & & Ahc & Ahp & BS & FP \\\\\n");
  fprintf(fp,"\\hline\n");
}

static void table1_end(FILE *fp)
{
  fprintf(fp,"\\hline\n");
  fprintf(fp,"\\end{tabular}\n\\end{table}\n\n");
}

static int miller_index(char *mil)
{
  int imil;
  
  for(imil=0; (imil<emlNR); imil++)
    if (strcmp(mil,miller[imil].name) == 0)
      break;
  if (imil == emlNR) {
    fprintf(stderr,"No such Miller atom %s\n",mil);
    exit(1);
  }
  return imil;
}

static int bosque_index(char *bos)
{
  int ielem;
  
  for(ielem=0; (ielem<eelemNR); ielem++)
    if (strcmp(bos,bosque[ielem].name) == 0)
      break;
  if (ielem == eelemNR) {
    fprintf(stderr,"No such element %s",bos);
    exit(1);
  }
  return ielem;
}

static void dump_table1(FILE *fp,int np,t_molprop pd[])
{
  int    i,j,nfit,nh,ielem,imil;
  char   hybrid[8];
  double ahc,bos,ahp;
  
  table1_header(fp,1);
  strcpy(hybrid,"SP ");
  
  for(i=0; (i<eatNR); i++) {
    int len = strlen(spoel[i].name);
    if (len < 3) {
      fprintf(stderr,"No such fragment %s",spoel[i].name);
    }
    hybrid[2] = spoel[i].name[len-2];
    nh        = spoel[i].nh;
    ielem = bosque_index(spoel[i].bosque);
    imil  = miller_index(spoel[i].miller);
    bos   = bosque[ielem].value + nh*bosque[eelemH].value;
    nfit  = 0;
    for(j=0; (j<np); j++) {
      if ((mp_num_polar(&(pd[j])) > 0) && (mp_get_polar(&(pd[j])) > 0))
	nfit += pd[j].frag_comp[i];
    }
    ahc = (4.0/(nh+miller[imil].nelec))*sqr(miller[imil].tau_ahc + 
					    nh*miller[emlH].tau_ahc);
    ahp = (nh*miller[emlH].alpha_ahp + miller[imil].alpha_ahp);
    fprintf(fp,"%s & %s & %s & %d & %8.2f & %8.2f & %8.2f & %8.2f \\\\\n",
	    spoel[i].atom,spoel[i].name,hybrid,nfit,
	    ahc,ahp,bos,spoel[i].all);
    if ((((i+1) % 16) == 0) && (i != eatNR-1)) {
      table1_end(fp);
      table1_header(fp,0);
    }
  }
  table1_end(fp);
}

static void table2_header(FILE *fp,int caption,int bTrain,int bQM)
{
  int i;
  char *ptr;
  
  ptr = bTrain ? "training" : "test";
  table_number(fp);
  fprintf(fp,"\\begin{%stable}[H]\n",bQM ? "sideways" : "");
  if (caption)
    fprintf(fp,"\\caption{Molecules with experimental polarizabilities as well as %sfragment or atom based polarizabilities of Miller {\\em et al.}~\\protect\\cite{Miller1979a,Miller90a} (Ahc, Ahp), Bosque and Sales~\\protect\\cite{Bosque2002a} (BS) or fragment polarizabilities (FP, this work). %sCalculated numbers that are more than %d\\%% off are printed in bold.%s A (*) behind the experimental value indicates that these were derived in this work. %s}\n\\label{%s}",
	    bQM ? "quantum chemical and " : "",
	    bQM ? "Basissets used are aug-cc-pVXZ and d-aug-cc-pVXZ. " : "",
	    toler,
	    bQM ? " Numbering is identical to that in Supplementary material, Table 1." : "",
	    bQM ? "$\\dagger$: for this calculation no convergence was achieved and hence no polarizability could be calculated." : "",
	    bQM ? "qm" : "empirical");
  else
    fprintf(fp,"\\caption{Continued}\n");
  if (bQM) {
    fprintf(fp,"\\begin{tabular}{lcccccccccc}\n\\hline\n");
    fprintf(fp,"Molecule & Experiment & HF & \\multicolumn{2}{c}{MP2} & \\multicolumn{2}{c}{B3LYP} & \\multicolumn{4}{c}{Empirical}\\\\\n");
    fprintf(fp,"   &  ");
    for(i=0; (i<eqmNR); i++)
      fprintf(fp,"& %s ",sbasis[i]);
    fprintf(fp,"\\\\\n\\hline\n");
  }
  else {
    fprintf(fp,"\\begin{tabular}{lccccc}\n\\hline\n");
    fprintf(fp,"Molecule & Experiment ");
    for(i=NQUANT; (i<eqmNR); i++)
      fprintf(fp,"& %s ",sbasis[i]);
    fprintf(fp,"\\\\\n\\hline\n");
  }
}

static void table2_end(FILE *fp,int bQM)
{
  fprintf(fp,"\\hline\n\\end{tabular}\n\\end{%stable}\n",
	  bQM ? "sideways" : "");
}

static void dump_table2(FILE *fp,int bVerbose,int bTrain,int bQM,
			int np,t_molprop pd[],int bDS)
{
  int i,j,j0,k,header,iline,iqm,caption=1,maxline;
  char buf[32],buf2[32],*ref;
  double val;
    
  table2_header(fp,caption,bTrain,bQM);  
  header  = 1;
  iline   = 0;
  j0 = bQM ? 0 : NQUANT;
  for(i=0; (i<np); i++) {
    iqm = 0;
    for(k=0; (k<NQUANT); k++)
      if (pd[i].qm[k] > 0)
	iqm++;
    if (((bTrain && (pd[i].weight > 0)) ||
	 (!bTrain && (pd[i].weight == 0))) &&
	(mp_num_polar(&(pd[i])) > 0) &&
	((bQM && (iqm > 0)) || !bQM)) {
      fprintf(fp,"%3d. %-15s",pd[i].nr,pd[i].molname);
      header = 0;
      fprintf(fp," &");
      for(k=0; (k<mp_num_polar(&(pd[i]))); k++,iline++) {
	if (k > 0)
	  fprintf(fp," &");
	val = mp_get_prop(&(pd[i]),"Polarizability",k);
	ref = mp_get_ref_prop(&(pd[i]),"Polarizability",k);
	if (val > 0) {
	  fprintf(fp," %8.3f",val);
	  if (strcmp(ref,"Spoel2007a") == 0)
	    fprintf(fp," (*)");
	  else
	    fprintf(fp,"~\\cite{%s} ",ref ? ref : "XXX");
	}
	if (k == 0) {
	  for(j=j0; (j<eqmNR); j++) 
	    if (pd[i].qm[j] > 0) {
	      if ((fabs(pd[i].qm[j] - val) > (val*toler*0.01))) {
		if (fabs(pd[i].qm[j] - val) > (val*tolerb*0.01))
		  fprintf(fp,"& {\\LARGE\\bf %8.2f} ",pd[i].qm[j]);
		else
		  fprintf(fp,"& {\\bf %8.2f} ",pd[i].qm[j]);
	      }
	      else
		fprintf(fp,"& %8.2f ",pd[i].qm[j]);
	    }
	    else
	      fprintf(fp,"& $\\dagger$");
	  fprintf(fp,"\\\\\n");
	}
	else {
	  if (bQM) 
	    fprintf(fp,"&&&&&&&&&\\\\\n"); 
	  else
	    fprintf(fp,"&&&&\\\\\n"); 
	}
      }
      maxline = 28;
      if (!bQM)
	maxline += 8;
      if (caption == 1)
	maxline -= 5;
      if (bDS)
	maxline = maxline*0.6;
	
      if ((iline >= maxline) && (i+1<np)) {
	caption = 0;
	table2_end(fp,bQM);
	table2_header(fp,caption,bTrain,bQM);
	header  = 1;
	iline   = 0;
      }
    }
  }
    table2_end(fp,bQM);
}

static void dump_table3(FILE *fp,int np,t_molprop pd[],int ntot)
{
  int    i,j,k,ns,n,nn[eqmNR];
  double aa[eqmNR],bb[eqmNR],R[eqmNR],rms,diff,ssx,ssy,ssxy,ssxx,ssyy;
  double ra[eqmNR],pp[eqmNR],ratio[eqmNR];
  double val;
  
  table_number(fp);
  fprintf(fp,"\\begin{table}[H]\n\\centering\n");
  fprintf(fp,"\\caption{Performance of the different methods for predicting molecular polarizability for different subsets of the test molecules. RMSD from experimental values (\\AA$^3$), in brackets the number of molecules in this particular subset. For the empirical methods (Ahc,Ahp~\\cite{Miller90a}, BS~\\cite{Bosque2002a} and FP (this work) all molecules were taken into account. At the bottom the \\%% difference, the correlation coefficient R, the regression coefficient a and the intercept b are given.}\n\\label{frag_rmsd}\n");
  fprintf(fp,"\\begin{tabular}{lcccccccccc}\n\\hline\n");
  fprintf(fp,"Method & HF & \\multicolumn{2}{c}{MP2} & \\multicolumn{2}{c}{B3LYP} & \\multicolumn{4}{c}{Empirical}\\\\\n");
  
  for(i=0; (i<eqmNR); i++) 
    fprintf(fp,"& %s ",sbasis[i]);
  fprintf(fp,"\\\\\n\\hline\n");
  
  for(i=0; (i<ecNR); i++) {
    n = 0;
    for(j=0; (j<np); j++) 
      if (pd[j].category[i] > 0)
	n++;
    fprintf(fp," %s (%d) ",ec_name[i],n);
    for(k=0; (k<eqmNR); k++) {
      rms = 0;
      n   = 0;
      for(j=0; (j<np); j++) {
	if (mp_num_polar(&(pd[j])) > 0) {
	  val = mp_get_polar(&(pd[j]));
	  if ((val > 0) && (pd[j].qm[k] > 0) && (pd[j].category[i] > 0)) {
	    rms += sqr(val - pd[j].qm[k]);
	    n++;
	  }
	}
      }
      if (n > 0) {
	if (k >= NQUANT)
	  fprintf(fp,"& %8.2f",sqrt(rms/n));
	else
	  fprintf(fp,"& %8.2f(%d)",sqrt(rms/n),n);
      }
      else
	fprintf(fp,"& -");
    }
    fprintf(fp,"\\\\\n");
  }
  for(k=0; (k<eqmNR); k++) {
    ra[k] = 0;
    pp[k] = 0;
    nn[k] = 0;
    ratio[k] = 0;
    for(j=0; (j<np); j++) {
      if (mp_num_polar(&(pd[j])) > 0) {
	val = mp_get_polar(&(pd[j]));
	if ((val > 0) && (pd[j].qm[k] > 0)) {
	  diff   = val - pd[j].qm[k];
	  ra[k] += sqr(diff);
	  pp[k] += 100*diff/val;
	  ratio[k] += pd[j].qm[k]/val;
	  nn[k]++;
	}
      }
    }
  }
  fprintf(fp,"All(%d)",nn[eqmNR-1]);
  for(k=0; (k<eqmNR); k++) {
    if (k < NQUANT) 
      fprintf(fp,"& %8.2f(%d)",sqrt(ra[k]/nn[k]),nn[k]);
    else
      fprintf(fp,"& %8.2f",sqrt(ra[k]/nn[k]));
  }
  fprintf(fp,"\\\\\n\\hline\n");
  fprintf(fp,"\\%% Diff ");
  for(k=0; (k<eqmNR); k++) {
    if (nn[k] > 0) {
      fprintf(fp,"& %8.2f",pp[k]/nn[k]);
    }
    else
      fprintf(fp,"& -");
  }
  fprintf(fp,"\\\\\n");
  for(k=0; (k<eqmNR); k++) {
    ssx = ssy = ssxx = ssxy = ssyy = 0;
    ns = 0;
    for(j=0; (j<np); j++) {
      if (mp_num_polar(&(pd[j])) > 0) {
	val = mp_get_polar(&(pd[j]));
	if ((val > 0) && (pd[j].qm[k] > 0)) {
	  ssx  += val;
	  ssy  += pd[j].qm[k];
	  ssxx += sqr(val);
	  ssyy += sqr(pd[j].qm[k]);
	  ssxy += val*pd[j].qm[k];
	  ns++;
	}
      }
    }
    ssx /= ns;
    ssy /= ns;
    ssxy -= ns*ssx*ssy;
    ssxx -= ns*sqr(ssx);
    ssyy -= ns*sqr(ssy);
    R[k] = 100*(ssxy*ssxy)/(ssxx*ssyy);
    bb[k] = ssxy/ssxx;
    aa[k] = ssy-bb[k]*ssx;
  }
  fprintf(fp,"R (\\%%) ");
  for(k=0; (k<eqmNR); k++) {
    fprintf(fp,"& %8.1f",R[k]);
  }
  fprintf(fp,"\\\\\n");
  fprintf(fp,"a ");
  for(k=0; (k<eqmNR); k++) {
    fprintf(fp,"& %8.2f",bb[k]);
  }
  fprintf(fp,"\\\\\n");
  fprintf(fp,"b ");
  for(k=0; (k<eqmNR); k++) {
    fprintf(fp,"& %8.2f",aa[k]);
  }
  fprintf(fp,"\\\\\n");
  fprintf(fp,"\\hline\n\\end{tabular}\n\\end{table}\n");
}

static void table4_header(FILE *fp,int caption)
{
  table_number(fp);
  fprintf(fp,"\\begin{sidewaystable}[H]\n\\centering\n");
  if (caption)
    fprintf(fp,"\\caption{Decomposition of molecules into fragments or Miller atoms~\\protect\\cite{Miller1979a,Miller90a}}\n\\label{frag_defs}\n");
  else
    fprintf(fp,"\\caption{Continued}\n");
  fprintf(fp,"\\begin{tabular}{llll}\n\\hline\n");
  fprintf(fp,"Molecule & Formula  & Fragments & Miller \\\\\n");
  fprintf(fp,"\\hline\n");
}

static void dump_table4(FILE *fp,int np,t_molprop pd[])
{
  int i,j,k,iline;
  char buf[32],buf2[32];
  
  table4_header(fp,1);  
  iline = 0;
  for(i=0; (i<np); i++) {
    if ((i == 0) || 
	((i > 0) && strcmp(pd[i].molname,pd[i-1].molname) != 0)) {
      fprintf(fp,"%3d. %s & %s &",++iline,pd[i].molname,pd[i].formula);
      for(j=0; (j<eatNR); j++)
	if (pd[i].frag_comp[j] > 0)
	  fprintf(fp,"%s$_{%d}$ ",spoel[j].name,pd[i].frag_comp[j]);
      fprintf(fp," &");
      for(j=0; (j<emlNR); j++)
	if (pd[i].emil_comp[j] > 0)
	  fprintf(fp,"%s$_{%d}$ ",miller[j].name,pd[i].emil_comp[j]);
      fprintf(fp," \\\\\n");
      if (((iline % 28) == 0) && (i+1<np)) {
	fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
	table4_header(fp,0);
      }
    }
  }
  fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
}

static void calc_frag_miller(int bTrain,int np,t_molprop pd[])
{
  int    j,k,Nelec;
  double ahc,ahp,bos,spol;
  
  for(j=0; (j<np); j++) {
    spol = 0;
    for(k=0; (k<eatNR); k++)
      spol += pd[j].frag_comp[k] * (bTrain ? spoel[k].train : spoel[k].all);
    pd[j].qm[eqmSpoel] = spol;
    bos = bosque[eelemNR].value;
    for(k=0; (k<eelemNR); k++)
      bos += bosque[k].value*pd[j].elem_comp[k];
    if ((pd[j].qm[eqmBosque] != 0) &&
	(fabs(pd[j].qm[eqmBosque] - bos) >= 0.01*toler*bos)) {
      fprintf(stderr,"Possible error in Bosque composition for %s. Computed %g iso %g\n",pd[j].molname,bos,pd[j].qm[eqmBosque]);
    }
    pd[j].qm[eqmBosque] = bos;
    ahc   = 0;
    Nelec = 0;
    ahp   = 0;
    for(k=0; (k<emlNR); k++) {
      ahc   += miller[k].tau_ahc*pd[j].emil_comp[k];
      ahp   += miller[k].alpha_ahp*pd[j].emil_comp[k];
      Nelec += miller[k].nelec*pd[j].emil_comp[k];
    }
    pd[j].qm[eqmMiller] = 4*sqr(ahc)/Nelec;
    pd[j].qm[eqmKang]   = ahp;
  }
}
  
static void calc_rmsd(int bTrain,double fit[],double test[],int np,
		      t_molprop pd[])
{
  int i,j,k,nfit[eqmNR],ntest[eqmNR],Nelec;
  double eval,qval;
  double fpol,apol;
  
  calc_frag_miller(bTrain,np,pd);
  
  for(i=0; (i<eqmNR); i++) {
    nfit[i]  = 0;
    ntest[i] = 0;
    fit[i]   = 0;
    test[i]  = 0;
    for(j=0; (j<np); j++) {
      if (mp_num_polar(&(pd[j])) > 0) {
	eval = mp_get_polar(&(pd[j]));
	if (eval > 0) {
	  qval = pd[j].qm[i];
	  
	  if (qval > 0) {
	    if (pd[j].weight > 0) {
	      nfit[i]++;
	      fit[i] += sqr(qval-eval);
	    }
	    else {
	      ntest[i]++;
	      test[i] += sqr(qval-eval);
	    }
	  }
	}
      }
    }
    if (nfit[i] > 0) 
      fit[i] = sqrt(fit[i]/nfit[i]);
    if (ntest[i] > 0)
      test[i] = sqrt(test[i]/ntest[i]);
  }
}

static void decompose_frag(FILE *fp,int bTrain,int np,t_molprop pd[])
{
  double *x,*atx;
  double **a,**at,**ata,fpp;
  int i,j,n,m=eatNR,nn;
  int test[eatNR];

  for(i=0; (i<eatNR); i++) {
    test[i] = 0;
    for(j=0; (j<np); j++) 
      if (pd[j].weight > 0) 
	test[i] += pd[j].frag_comp[i];
    if (test[i] == 0) {
      fprintf(stderr,"You have no molecules with group %s\n",
	      spoel[i].name);
      exit(1);
    }
  }
  snew(x,np);
  for(i=n=0; (i<np); i++) {
    if ((pd[i].weight > 0) && 
	(mp_num_polar(&(pd[i])) > 0)) {
      x[n] = mp_get_polar(&(pd[i]));
      n++;
    }
  }
  a      = alloc_matrix(n,m);
  at     = alloc_matrix(m,n);
  ata    = alloc_matrix(m,m);
  for(i=nn=0; (i<np); i++) {
    if ((pd[i].weight > 0) && 
	(mp_num_polar(&(pd[i])) > 0)){
      for(j=0; (j<m); j++) {
	a[nn][j]  = pd[i].frag_comp[j];
	at[j][nn] = pd[i].frag_comp[j];
      }
      nn++;
    }
  }
  if (n != nn) {
    fprintf(stderr,"Consistency error %s, line %d\n",__FILE__,__LINE__);
    exit(1);
  }
  matrix_multiply(fp,n,m,a,at,ata);
  matrix_invert(fp,m,ata);
  snew(atx,m);
  for(i=0; (i<m); i++) 
    for(j=0; (j<n); j++)
      atx[i] += at[i][j]*x[j];
  
  for(i=0; (i<m); i++) {
    fpp = 0;
    for(j=0; (j<m); j++)
      fpp += ata[i][j]*atx[j];
    if (bTrain)
      spoel[i].train = fpp;
    else
      spoel[i].all = fpp;
  }
}

static void write_xvg(int np,t_molprop pd[])
{
  char *fn[2] = { "qm.xvg", "empir.xvg" };
  FILE *fp;
  int i,j,n,i0,i1;
  double val; 
  
  for(n=0; (n<2); n++) {
    fp = fopen(fn[n],"w");
    if (n == 0) {
      i0  = 0;
      i1  = NQUANT;
      val = 19.5;
    }
    else {
      i0  = NQUANT;
      i1  = eqmNR;
      val = 51.5;
    }
    for(i=i0; (i<i1); i++) {
      fprintf(fp,"@type xy\n");
      for(j=0; (j<np); j++) {
	if ((mp_num_polar(&(pd[j])) > 0) &&
	    (mp_get_polar(&(pd[j])) > 0) && 
	    (pd[j].qm[i] > 0))
	  fprintf(fp,"%8.3f  %8.3f\n",mp_get_polar(&(pd[j])),pd[j].qm[i]);
      }
      fprintf(fp,"&\n");
    }
    fprintf(fp,"@type xy\n");
    fprintf(fp,"1.5 1.5\n%6.1f %6.1f\n",val,val);
    fprintf(fp,"&\n");
    fclose(fp);
  }
}

static void write_csv(int np,t_molprop pd[])
{
  FILE *fp,*gp;
  int i,j;
  
  fp = fopen("matrix.csv","w");
  gp = fopen("vector.csv","w");
  for(i=0; (i<np); i++) {
    if ((pd[i].weight > 0) && 
	(mp_num_polar(&(pd[i])) > 0) &&
	(mp_get_polar(&(pd[i])) > 0)) {
      fprintf(gp,"%10.3f\n",mp_get_polar(&(pd[i])));
      for(j=0; (j<eatNR-1); j++)
	fprintf(fp,"  %5d,",pd[i].frag_comp[j]);
      fprintf(fp,"  %5d\n",pd[i].frag_comp[eatNR-1]);
    }
  }
  fclose(fp);
  fclose(gp);
}

static int comp_pd_name(const void *a,const void *b)
{
  t_molprop *pa,*pb;
  pa = (t_molprop *)a;
  pb = (t_molprop *)b;
  
  return strcasecmp(pa->molname,pb->molname);
}

static int comp_pd_elem(const void *a,const void *b)
{
  t_molprop *pa,*pb;
  pa = (t_molprop *)a;
  pb = (t_molprop *)b;
  
  return strcmp(pa->formula,pb->formula);
}

static int comp_pd_elem2(const void *a,const void *b)
{
  int i,oo=0,nn=0;
  t_molprop *pa,*pb;
  pa = (t_molprop *)a;
  pb = (t_molprop *)b;
  
  for(i=2; ((i<eelemNR) && (nn == 0)); i++)
    if ((pa->elem_comp[i] == 0) && (pb->elem_comp[i] > 0)) {
      nn = -1;
      break;
    }
    else if ((pb->elem_comp[i] == 0) && (pa->elem_comp[i] > 0)) {
      nn = 1;
      break;
    }
      
  if (nn != 0)
    return nn;
    
  for(i=0; (i<eelemNR); i++) 
    if ((oo = (pa->elem_comp[my_order[i]] - pb->elem_comp[my_order[i]])) != 0)
      return oo;
 
  return strcasecmp(pa->molname,pb->molname);
}

static void dump_n2t(char *fn)
{
  FILE   *fp;
  int    i,j,ielem;
  double q = 0, m = 0;
  
  fp = fopen(fn,"w");
  for(i=0; (i<eatNR+eatExtra); i++) {
    ielem  = bosque_index(spoel[i].bosque);
    fprintf(fp,"%3s  %4s  %9g  %9g  %1d",
	    spoel[i].bosque,spoel[i].name,
	    spoel[i].all,bosque[ielem].mass,
	    spoel[i].maxbond);
    for(j=0; (j<spoel[i].nh); j++)
      fprintf(fp," H %.3f",spoel[eatNR].blen);
    for( ; (j<spoel[i].maxbond); j++)
      fprintf(fp," * %.3f",spoel[i].blen);
    fprintf(fp,"\n");
  }
  
  fclose(fp);
}

static int read_dipoles(char *fn) 
{
  FILE *fp;
  
  fp = fopen(fn,"r");
  
  fclose(fp);
  
  return 0;
}

int main(int argc,char *argv[])
{
  FILE  *fp,*gp;
  int    i,j,np,nspoel,nbosque,ntot,nqm,nqqm,nqmtot=0;
  double fit1[eqmNR],test1[eqmNR];
  double fit2[eqmNR],test2[eqmNR];
  double *w;
  t_molprop *mp=NULL;
  
  if (argc < 2) 
    fatal_error("Please give names of database files!","");
  mp = merge_xml(argc,argv,NULL,NULL,"double_dip.dat",&np);
    
  check_formula(np,mp);
  snew(w,np);
  qsort(mp,np,sizeof(mp[0]),comp_pd_elem2);
  
  for(i=0; (i<np); i++) {
    mp[i].nr = 1+i;
    w[i] = mp[i].weight;
    if ((mp[i].weight > 0) && 
	(mp_num_polar(&(mp[i])) > 0) &&
	(mp_get_polar(&(mp[i])) == 0)) {
      fprintf(stderr,"Inconsistency for %s: weight = %g, experiment = 0\n",
	      mp[i].molname,mp[i].weight);
      exit(1);
    }
  }
  
  nbosque = nspoel = nqm = ntot = 0;
  for(i=0; (i<np); i++) {
    if ((mp_num_polar(&(mp[i])) > 0) && (mp_get_polar(&(mp[i])) > 0)) {
      char *ref = mp_get_ref_polar(&(mp[i]));
      mp[i].weight = 1;
      if (strcasecmp(ref,"Spoel2007a") == 0)
	nspoel++;
      else if (strcasecmp(ref,"Bosque2002a") == 0)
	nbosque++;
      nqqm = 0;
      for(j=0; (j<NQUANT); j++) 
	if (mp[i].qm[j] > 0)
	  nqqm++;
      if (nqqm)
	nqm++;
      nqmtot += nqqm;
      ntot++;
    }
    else {
      /*fprintf(stderr,"No experimental data for %s\n",mp[i].molname);*/
    }
    if ((i > 0) && (strcasecmp(mp[i].molname,mp[i-1].molname) == 0))
      fprintf(stderr,"Double entry %s\n",mp[i].molname);
  }
  printf("--------------------------------------------------\n");
  printf("              Some Statistics\n");
  printf("There are %d entries with experimental data\n",ntot);
  printf("There are %d experiments with Spoel2007a  as reference\n",nspoel);
  printf("There are %d experiments with Bosque2002a as reference\n",nbosque);
  printf("There are %d entries with quantum chemistry results, %d missing\n",
	 nqm,nqm*5-nqmtot);
  printf("There are %d Spoel atom types\n",eatNR);
  printf("--------------------------------------------------\n");
  
  gp = ffopen("tune_pol.log","w");
  decompose_frag(gp,0,np,mp);
  fclose(gp);
  calc_rmsd(0,fit2,test2,np,mp);

  if ((fp=fopen("molecules.tex","w")) == NULL)
    exit(1);
  
  dump_table1(fp,np,mp);
  tabnum++;
  dump_table2(fp,0,1,1,np,mp,1);
  tabnum++;
  dump_table3(fp,np,mp,ntot);
  tabnum = 0;
  fclose(fp);
  
  fp = ffopen("molecules-sup.tex","w");

  dump_table2(fp,0,1,0,np,mp,0);
  tabnum++;
  dump_table4(fp,np,mp);
  fclose(fp);
  dump_n2t("ffsm.n2t");    
  write_xvg(np,mp);
    
  return 0;
}
