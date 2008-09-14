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
#include "vec.h"
#include "gmx_fatal.h"
#include "molprop.h"
#include "molprop_util.h"
#include "gentop_matrix.h"
#include "poldata_xml.h"

static int tabnum = 0;
static int toler  = 5;  /* Tolerance in % */
static int tolerb = 50; /* Tolerance in % */

/*static char *sbasis[eqmNR] = {
  "6-31G","pVDZ","pVTZ","pVTZ","d-pVTZ", "Ahc", "Ahp", "BS", "FP" 
  };*/

static char *sbasis[] = {
  "6-31G","aug-DZ","aug-TZ","aug-TZ","d-aug-TZ", "Ahc", "Ahp", "BS", "FP" 
};
#define eqmNR asize(sbasis)
#define NQUANT 5

double mp_get_polar(gmx_molprop_t mp)
{
  double val,err;
  char   *ref;
  
  if (gmx_molprop_search_property(mp,eMOLPROP_Exp,"Polarizability",
				  &val,&err,NULL,&ref) == 1) {
    sfree(ref);
    return val;
  }
  return 0;
}

int mp_num_polar(gmx_molprop_t mp) 
{
  int npol=0,eMP;
  double val,err;
  char *prop_name,*prop_method,*prop_ref;
  
  while (gmx_molprop_get_property(mp,&eMP,&prop_name,&val,&err,
				  &prop_method,&prop_ref) == 1) {
    if (strcasecmp(prop_name,"Polarizability") == 0)
      npol++;
    sfree(prop_name);
    sfree(prop_method);
    sfree(prop_ref);
  }
  return npol;
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

static void dump_table1(FILE *fp,gmx_poldata_t pd,
			int np,gmx_molprop_t mp[])
{
  int    i,j,nfit,nh,ielem,imil,atomnumber;
  int    nhydrogen,charge,hybridization;
  double ahc,ahp,ahc_H,ahp_H,mil_pol,bos_pol,bos_H,spoel_pol,blength,pol_tmp,pol_error;
  char   *name,*elem,*miller_equiv,*prop_reference,group[32];
  
  table1_header(fp,1);
  
  if (gmx_poldata_get_bosque(pd,"H",&bos_H) == NULL) 
    gmx_fatal(FARGS,"Can not find Bosque polarizability for H");
  
  if (gmx_poldata_get_miller(pd,"H",&atomnumber,&ahc_H,&ahp_H,NULL) != NULL) 
    gmx_fatal(FARGS,"Can not find Miller polarizability for H");
  
  while(name = gmx_poldata_get_spoel(pd,NULL,&elem,
				     &miller_equiv,&nhydrogen,
				     &charge,&hybridization,
				     &spoel_pol,&blength)) {
    
    /* Determine Miller and Bosque polarizabilities for this Spoel element */
    ahc = ahp = 0;
    if (gmx_poldata_get_miller(pd,miller_equiv,&atomnumber,&ahc,&ahp,NULL) != NULL) {
      ahc = (4.0/(nhydrogen+atomnumber))*sqr(ahc + nhydrogen*ahc_H);
      ahp = (nhydrogen*ahp_H + ahp);
    }
    bos_pol = 0;
    if (gmx_poldata_get_bosque(pd,elem,&bos_pol) != NULL)
      bos_pol += nhydrogen*bos_H;
    
    /* Compute how many molecules contributed to the optimization 
     * of this polarizability.
     */
    nfit  = 0;
    for(j=0; (j<np); j++) {
      if (gmx_molprop_search_property(mp[j],eMOLPROP_Exp,"Polarizability",
				      &pol_tmp,&pol_error,
				      NULL,&prop_reference) == 1)
	nfit += gmx_molprop_count_composition_atoms(mp[j],"spoel",name);
    }
    /* Construct group name from element composition */
    if (nhydrogen > 1)
      sprintf(group,"%sH%d",elem,nhydrogen);
    else if (nhydrogen == 1)
      sprintf(group,"%sH",elem);
    else
      strcpy(group,elem);
      
    fprintf(fp,"%s & %s & SP%d & %d & %8.2f & %8.2f & %8.2f & %8.2f \\\\\n",
	    group,name,hybridization,nfit,ahc,ahp,bos_pol,spoel_pol);
    if ((((j+1) % 16) == 0)) {
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
    for(i=0; (i<eqmNR); i++)
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
			int np,gmx_molprop_t mp[],int bDS)
{
  int i,j,j0,k,header,iline,iqm,caption=1,maxline,result;
  char buf[32],buf2[32],*exp_ref,*calc_ref,*prop_name;
  double val,exp_val,exp_error,weight;
  double val_calc[eqmNR],err_calc[eqmNR];
  int    eMP,found_calc[eqmNR];
  
  table2_header(fp,caption,bTrain,bQM);  
  header  = 1;
  iline   = 0;
  for(i=0; (i<np); i++) {
    iqm = 0;
    
    if (gmx_molprop_search_property(mp[i],eMOLPROP_Exp,"Polarizability",
				    &exp_val,&exp_error,NULL,&exp_ref) == 1) {
      sfree(exp_ref);
      for(j=0; (j<eqmNR); j++) {
	found_calc[j] = 
	  gmx_molprop_search_property(mp[i],eMOLPROP_Calc,"Polarizability",
				      &(val_calc[j]),&(err_calc[j]),
				      sbasis[j],&calc_ref);
	sfree(calc_ref);
	if (found_calc[j])
	  iqm++;
      }
    
      weight = gmx_molprop_get_weight(mp[i]);
      if (((bTrain && (weight > 0)) || (!bTrain && (weight == 0))) &&
	  ((bQM && (iqm > 0)) || !bQM)) {
	fprintf(fp,"%-15s",gmx_molprop_get_molname(mp[i]));
	header = 0;
	fprintf(fp," &");
	for(k=0; (k<mp_num_polar(mp[i])); k++,iline++) {
	  if (k > 0)
	    fprintf(fp," &");
	  
	  fprintf(fp," %8.3f",exp_val);
	  if (strcmp(exp_ref,"Spoel2007a") == 0)
	    fprintf(fp," (*)");
	  else
	    fprintf(fp,"~\\cite{%s} ",exp_ref ? exp_ref : "XXX");
	  sfree(exp_ref);
	}
	if (k == 0) {
	  for(j=0; (j<eqmNR); j++) 
	    if (val_calc[j] > 0) {
	      if ((fabs(val_calc[j] - exp_val) > (exp_val*toler*0.01))) {
		if (fabs(val_calc[j] - exp_val) > (exp_val*tolerb*0.01))
		  fprintf(fp,"& {\\LARGE\\bf %8.2f} ",val_calc[j]);
		else
		  fprintf(fp,"& {\\bf %8.2f} ",val_calc[j]);
	      }
	      else
		fprintf(fp,"& %8.2f ",val_calc[j]);
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

typedef struct {
  int ncat;
  char **cat;
  int  *count;
  real *rms;
} t_cats;

static void add_cat(t_cats *cats,char *catname)
{
  int i;
  
  for(i=0; (i<cats->ncat); i++) 
    if (strcasecmp(catname,cats->cat[i]) == 0) {
      cats->count[i]++;
      break;
    }
  if (i == cats->ncat) {
    cats->ncat++;
    srenew(cats->cat,cats->ncat);
    srenew(cats->count,cats->ncat);
    cats->cat[cats->ncat-1] = strdup(catname);
    cats->count[cats->ncat-1]++;
  }
}

static void dump_table3(FILE *fp,int np,gmx_molprop_t pd[],int ntot)
{
  int    i,j,k,ns,n,nn[eqmNR];
  double aa[eqmNR],bb[eqmNR],R[eqmNR],rms,diff,ssx,ssy,ssxy,ssxx,ssyy;
  double ra[eqmNR],pp[eqmNR],ratio[eqmNR];
  double exp_val,qm_val,qm_err;
  char   *cname,*qm_ref;
  t_cats *cats;
  
  table_number(fp);
  fprintf(fp,"\\begin{table}[H]\n\\centering\n");
  fprintf(fp,"\\caption{Performance of the different methods for predicting molecular polarizability for different subsets of the test molecules. RMSD from experimental values (\\AA$^3$), in brackets the number of molecules in this particular subset. For the empirical methods (Ahc,Ahp~\\cite{Miller90a}, BS~\\cite{Bosque2002a} and FP (this work) all molecules were taken into account. At the bottom the \\%% difference, the correlation coefficient R, the regression coefficient a and the intercept b are given.}\n\\label{frag_rmsd}\n");
  fprintf(fp,"\\begin{tabular}{lcccccccccc}\n\\hline\n");
  fprintf(fp,"Method & HF & \\multicolumn{2}{c}{MP2} & \\multicolumn{2}{c}{B3LYP} & \\multicolumn{4}{c}{Empirical}\\\\\n");
  
  for(i=0; (i<eqmNR); i++) 
    fprintf(fp,"& %s ",sbasis[i]);
  fprintf(fp,"\\\\\n\\hline\n");
  
  snew(cats,1);
  
  for(j=0; (j<np); j++) 
    while ((cname = gmx_molprop_get_category(pd[j])) != NULL)
      add_cat(cats,cname);

  for(i=0; (i<cats->ncat); i++) {      
    fprintf(fp," %s (%d) ",cats->cat[i],cats->count[i]);
    for(k=0; (k<eqmNR); k++) {
      rms = 0;
      n   = 0;
      for(j=0; (j<np); j++) {
	if ((gmx_molprop_search_category(pd[j],cats->cat[i]) == 1) &&
	    ((exp_val = mp_get_polar(pd[j])) > 0) &&
	    (gmx_molprop_search_property(pd[j],eMOLPROP_Calc,"Polarizability",
					 &qm_val,&qm_err,sbasis[k],&qm_ref) == 1)) {
	  rms += sqr(exp_val - qm_val);
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
    fprintf(fp,"\\\\\n");
  }
  for(k=0; (k<eqmNR); k++) {
    ra[k] = 0;
    pp[k] = 0;
    nn[k] = 0;
    ratio[k] = 0;
    for(j=0; (j<np); j++) {
      if (((exp_val = mp_get_polar(pd[j])) > 0) &&
	  (gmx_molprop_search_property(pd[j],eMOLPROP_Calc,"Polarizability",
				       &qm_val,&qm_err,sbasis[k],&qm_ref) == 1)) {
	diff   = exp_val - qm_val;
	ra[k] += sqr(diff);
	pp[k] += 100*diff/exp_val;
	ratio[k] += qm_val/exp_val;
	nn[k]++;
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
      if (((exp_val = mp_get_polar(pd[j])) > 0) &&
	  (gmx_molprop_search_property(pd[j],eMOLPROP_Calc,"Polarizability",
				       &qm_val,&qm_err,sbasis[k],&qm_ref) == 1)) {
	ssx  += exp_val;
	ssy  += qm_val;
	ssxx += sqr(exp_val);
	ssyy += sqr(qm_val);
	ssxy += exp_val*qm_val;
	ns++;
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

static void dump_table4(FILE *fp,int np,gmx_molprop_t pd[])
{
  int i,j,k,iline,natom;
  char buf[32],buf2[32];
  char *atomname,*comp;
  
  table4_header(fp,1);  
  iline = 0;
  for(i=0; (i<np); i++) {
    fprintf(fp,"%3d. %s & %s &",++iline,
	    gmx_molprop_get_molname(pd[i]),
	    gmx_molprop_get_formula(pd[i]));
    while((comp = gmx_molprop_get_composition(pd[i])) != NULL) {
      if (strcasecmp(comp,"spoel") == 0) {
	while(gmx_molprop_get_composition_atom(pd[i],comp,&atomname,&natom) == 1) {
	  fprintf(fp,"%s$_{%d}$ ",atomname,natom);
	}
	fprintf(fp," &");
      }
      else if (strcasecmp(comp,"miller") == 0) {
	while(gmx_molprop_get_composition_atom(pd[i],comp,&atomname,&natom) == 1) {
	  fprintf(fp,"%s$_{%d}$ ",atomname,natom);
	}
	fprintf(fp," \\\\\n");
      }
    }
    if (((iline % 28) == 0) && (i+1<np)) {
      fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
      table4_header(fp,0);
    }
  }
  fprintf(fp,"\\hline\n\\end{tabular}\n\\end{sidewaystable}\n");
}

static void calc_frag_miller(int bTrain,gmx_poldata_t pd,
			     int np,gmx_molprop_t mp[])
{
  int    j,k,Nelec,natom,charge,hybridization,nhydrogen,atomnumber;
  double ahc,ahp,bos,spoel,bos0,polar,blength;
  double tau_ahc,alpha_ahp;
  char   *comp,*atomname,*elem,*miller_equiv,*ref;
  int    spoel_support,miller_support,bosque_support;
  
  if (gmx_poldata_get_bosque(pd,"0",&bos0) == NULL)
    gmx_fatal(FARGS,"Can not find Bosque polarizability for 0");
  
  for(j=0; (j<np); j++) {
    ahc = ahp = spoel = 0;
    bos = bos0;
    Nelec = 0;
    spoel_support = 1;
    miller_support = 1;
    bosque_support = 1;
      
    while ((comp = gmx_molprop_get_composition(mp[j])) != NULL) {
      while (gmx_molprop_get_composition_atom(mp[j],comp,&atomname,&natom) == 1) {
	if (strcasecmp(comp,"spoel") == 0) {
	  if (gmx_poldata_get_spoel(pd,atomname,NULL,NULL,NULL,NULL,NULL,
				    &polar,NULL) != NULL) {
	    spoel += polar*natom;
	  }
	  else 
	    spoel_support = 0;
	}
	else if (strcasecmp(comp,"miller") == 0) {
	  if (gmx_poldata_get_miller(pd,atomname,&atomnumber,
				     &tau_ahc,&alpha_ahp,NULL) != NULL) {
	    ahc   += tau_ahc*natom; 
	    ahp   += alpha_ahp*natom;
	    Nelec += atomnumber*Nelec;
	  }
	  else
	    miller_support = 0;
	}
	else if (strcasecmp(comp,"bosque") == 0) {
	  if (gmx_poldata_get_bosque(pd,atomname,&polar) != NULL) {
	    bos += polar*natom;
	  }
	  else
	    bosque_support = 0;
	}
	else
	  gmx_fatal(FARGS,"Unknown composition %s\n",comp);
      }
    }

    if (spoel_support == 1) 
      gmx_molprop_add_property(mp[j],eMOLPROP_Calc,"Polarizability",
			       spoel,0,"spoel","Spoel2009a");
    if (bosque_support == 1)
      gmx_molprop_add_property(mp[j],eMOLPROP_Calc,"Polarizability",
			       bos,0,"bosque","Bosque2002a");

    if (miller_support == 1) {
      ahc = 4*sqr(ahc)/Nelec;
      gmx_molprop_add_property(mp[j],eMOLPROP_Calc,"Polarizability",
			       ahc,0,"ahc","Miller");
      gmx_molprop_add_property(mp[j],eMOLPROP_Calc,"Polarizability",
			       ahp,0,"ahp","Khan");
    }
  }
}
  
static void calc_rmsd(int bTrain,double fit[],double test[],
		      gmx_poldata_t pd,int np,gmx_molprop_t mp[])
{
  int i,j,k,nfit[eqmNR],ntest[eqmNR],Nelec;
  double eval,qval,qerr;
  double fpol,apol;
  char   *ref;
  
  calc_frag_miller(bTrain,pd,np,mp);
  
  for(i=0; (i<eqmNR); i++) {
    nfit[i]  = 0;
    ntest[i] = 0;
    fit[i]   = 0;
    test[i]  = 0;
    for(j=0; (j<np); j++) {
      if ((eval = mp_get_polar(mp[j])) > 0) {
	if (gmx_molprop_search_property(mp[j],eMOLPROP_Calc,
					"Polarizability",&qval,&qerr,
					sbasis[i],&ref) == 1) {
	
	  if (qval > 0) {
	    if (gmx_molprop_get_weight(mp[j]) > 0) {
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

static int decompose_frag(FILE *fp,int bTrain,
			  gmx_poldata_t pd,int np,gmx_molprop_t mp[])
{
  double *x,*atx;
  double **a,**at,**ata,fpp;
  double pol,blength;
  char *elem,*miller_equiv,*name,**atype=NULL;
  int i,j,n,nn,nhydrogen,charge,hybridization;
  int *test=NULL,ntest=0;

  while ((name = gmx_poldata_get_spoel(pd,NULL,&elem,&miller_equiv,
				       &nhydrogen,&charge,
				       &hybridization,&pol,&blength)) != NULL) {
    srenew(test,++ntest);
    srenew(atype,ntest);
    test[ntest-1] = 0;
    atype[ntest-1] = strdup(name);
    for(j=0; (j<np); j++) 
      if (gmx_molprop_get_weight(mp[j]) > 0) 
	test[ntest-1] += gmx_molprop_count_composition_atoms(mp[j],"spoel",name);
    
    if (test[ntest-1] == 0) 
      gmx_fatal(FARGS,"You have no molecules with group %s\n",name);
  }
  snew(x,np);
  for(i=n=0; (i<np); i++) {
    if ((gmx_molprop_get_weight(mp[i]) > 0) && 
	((pol = mp_get_polar(mp[i])) > 0)) {
      x[n] = pol;
      n++;
    }
  }
  a      = alloc_matrix(n,ntest);
  at     = alloc_matrix(ntest,n);
  ata    = alloc_matrix(ntest,ntest);
  for(i=nn=0; (i<np); i++) {
    if ((gmx_molprop_get_weight(mp[i]) > 0) && 
	(mp_num_polar(mp[i]) > 0)) {
      for(j=0; (j<ntest); j++) {
	a[nn][j] = at[j][nn] = gmx_molprop_count_composition_atoms(mp[j],"spoel",atype[j]);
      }
      nn++;
    }
  }
  if (n != nn) {
    fprintf(stderr,"Consistency error %s, line %d\n",__FILE__,__LINE__);
    exit(1);
  }
  matrix_multiply(fp,n,ntest,a,at,ata);
  matrix_invert(fp,ntest,ata);
  snew(atx,ntest);
  for(i=0; (i<ntest); i++) 
    for(j=0; (j<n); j++)
      atx[i] += at[i][j]*x[j];
  
  for(i=0; (i<ntest); i++) {
    fpp = 0;
    for(j=0; (j<ntest); j++)
      fpp += ata[i][j]*atx[j];
    /*if (bTrain)
      spoel[i].train = fpp;
      else*/
    gmx_poldata_set_spoel(pd,atype[j],fpp);
  }
  return ntest;
}

static void write_xvg(int np,gmx_molprop_t pd[])
{
  char *fn[2] = { "qm.xvg", "empir.xvg" };
  FILE *fp;
  int i,j,n,i0,i1;
  double exp_val,qm_val,qm_err,pol_max; 
  
  for(n=0; (n<2); n++) {
    fp = fopen(fn[n],"w");
    if (n == 0) {
      i0  = 0;
      i1  = NQUANT;
      pol_max = 19.5;
    }
    else {
      i0  = NQUANT;
      i1  = eqmNR;
      pol_max = 51.5;
    }
    for(i=i0; (i<i1); i++) {
      fprintf(fp,"@type xy\n");
      for(j=0; (j<np); j++) {
	if ((exp_val = (mp_get_polar(pd[j])) > 0) &&
	    (gmx_molprop_search_property(pd[j],eMOLPROP_Calc,"Polarizability",
					 &qm_val,&qm_err,sbasis[i],NULL) == 1)) {
	  fprintf(fp,"%8.3f  %8.3f\n",exp_val,qm_val);
	}
	fprintf(fp,"&\n");
      }
    }
    fprintf(fp,"@type xy\n");
    fprintf(fp,"1.5 1.5\n%6.1f %6.1f\n",pol_max,pol_max);
    fprintf(fp,"&\n");
    fclose(fp);
  }
}

int main(int argc,char *argv[])
{
  FILE  *fp,*gp;
  int    i,j,np,nspoel,nbosque,ntot,nqm,nqqm,nqmtot=0,nspoel_atypes;
  double fit1[eqmNR],test1[eqmNR];
  double fit2[eqmNR],test2[eqmNR];
  double *w;
  char   *ref,*molname[2];
  int    cur = 0;
#define prev (1-cur)
  gmx_molprop_t *mp=NULL;
  gmx_atomprop_t ap;
  gmx_poldata_t  pd;
  
  if (argc < 2) 
    gmx_fatal(FARGS,"Please give names of database files!","");
  mp = merge_xml(argc,argv,NULL,NULL,"double_dip.dat",&np);
  ap = gmx_atomprop_init();
  pd = gmx_poldata_read("gentop.xml");
    
  generate_composition(np,mp,pd,ap);
  generate_formula(np,mp,ap);

  snew(w,np);
  gmx_molprop_sort(np,mp,2,ap);
  
  for(i=0; (i<np); i++) {
    w[i] = gmx_molprop_get_weight(mp[i]);
    molname[cur] = gmx_molprop_get_molname(mp[i]);
    if ((w[i] > 0) && 
	(mp_num_polar(mp[i]) > 0) &&
	(mp_get_polar(mp[i]) == 0)) {
      fprintf(stderr,"Inconsistency for %s: weight = %g, experiment = 0\n",
	      molname[cur],w[i]);
      exit(1);
    }
  }
  
  nbosque = nspoel = nqm = ntot = 0;
  for(i=0; (i<np); i++) {
    if (gmx_molprop_search_property(mp[i],eMOLPROP_Exp,
				    "Polarizability",NULL,NULL,
				    NULL,&ref) == 1) {
      gmx_molprop_set_weight(mp[i],1);
      if (strcasecmp(ref,"Spoel2007a") == 0)
	nspoel++;
      else if (strcasecmp(ref,"Bosque2002a") == 0)
	nbosque++;
	
      nqqm = 0;
      while (gmx_molprop_search_property(mp[i],eMOLPROP_Calc,
					 "Polarizability",NULL,NULL,
					 NULL,NULL) == 1) 
	nqqm++;
	
      if (nqqm)
	nqm++;
      nqmtot += nqqm;
      ntot++;
    }
    else {
      /*fprintf(stderr,"No experimental data for %s\n",mp[i].molname);*/
    }
    if ((i > 0) && (strcasecmp(molname[cur],molname[prev]) == 0))
      fprintf(stderr,"Double entry %s\n",molname[cur]);
    cur=prev;
  }
  
  gp = ffopen("tune_pol.log","w");
  nspoel_atypes = decompose_frag(gp,0,pd,np,mp);
  fclose(gp);
  calc_rmsd(0,fit2,test2,pd,np,mp);

  printf("--------------------------------------------------\n");
  printf("              Some Statistics\n");
  printf("There are %d entries with experimental data\n",ntot);
  printf("There are %d experiments with Spoel2007a  as reference\n",nspoel);
  printf("There are %d experiments with Bosque2002a as reference\n",nbosque);
  printf("There are %d entries with quantum chemistry results, %d missing\n",
	 nqm,nqm*5-nqmtot);
  printf("There are %d Spoel atom types\n",nspoel_atypes);
  printf("--------------------------------------------------\n");
  
  fp = ffopen("molecules.tex","w");
  
  dump_table1(fp,pd,np,mp);
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
  gmx_poldata_write("tunepol.xml",pd);
  write_xvg(np,mp);
    
  return 0;
}
