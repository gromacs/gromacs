#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "string2.h"
#include "typedefs.h"
#include "pdb2top.h"
#include "atomprop.h"
#include "poldata.h"
#include "poldata_xml.h"

static int maxline = 30;

static void begin_table(FILE *tp,char *caption,char *label,char *colinfo,char *hbuf)
{
  fprintf(tp,"\\begin{table}[H]\n\\centering\n");
  fprintf(tp,"\\caption{%s}\n",caption);
  fprintf(tp,"\\label{%s}\n",label);
  fprintf(tp,"\\begin{tabular}{%s}\n\\hline\n",colinfo);
  fprintf(tp,"%s\\\\\n\\hline\n",hbuf);
}

static void end_table(FILE *tp)
{
  fprintf(tp,"\\hline\n\\end{tabular}\n\\end{table}\n\n");
}

static void do_atypes(FILE *fp,FILE *tp,gmx_poldata_t pd,gmx_atomprop_t aps)
{
  char hbuf[1024],colinfo[1024];
  double polarizability,sig_pol,valence;
  char *smname,*elem,*desc,*gt_type,*gt_old,*charge;
  char *neighbors,*vdwparams,*geometry;
  int  numbonds,atomnumber,nline,npage;
  real mass;
  
  strcpy(hbuf,"Type & Description & Elem & Val & $\\alpha$ & Van der Waals");
  strcpy(colinfo,"clcccc");
  begin_table(tp,"Atom types defined by the Alexandria force field","atypes",
              colinfo,hbuf);
              
  fprintf(fp,"[ defaults ]\n");
  fprintf(fp,"; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n");
  fprintf(fp,"1               1               yes             1       1\n\n");

  fprintf(fp,"[ atomtypes ]\n");
  fprintf(fp,"; name      at.num  mass     charge ptype  repulsion  dispersion\n");

  gt_old = NULL;
  nline = 2;
  npage = 0;
  while (1 == gmx_poldata_get_atype(pd,&elem,&desc,&gt_type,
                                    NULL,&charge,&valence,&polarizability,&sig_pol,
                                    &vdwparams))
    {
      if (gmx_atomprop_query(aps,epropMass,"",elem,&mass)) 
	{
	  atomnumber = gmx_atomprop_atomnumber(aps,elem);
	  if ((NULL == gt_old) || (strcmp(gt_old,gt_type) != 0)) {
	    fprintf(fp,"%5s   %3d  %12.6f  %10.4f  A  %-s\n",
		    gt_type,atomnumber,mass,0.0,vdwparams);
	    fprintf(fp,"%5ss  %3d  %12.6f  %10.4f  S  0  0\n",
		    gt_type,0,0.0,0.0);
            if (0 == (nline % maxline)) {
              end_table(tp);
              fprintf(tp,"\\addtocounter{table}{-1}\n");
              npage++;
              nline = 1;
              begin_table(tp,"Atom types, continued",gmx_itoa(npage),colinfo,hbuf);
            }
            fprintf(tp,"%s & %s & %s & %g & %.2f(%.2f) & %s\\\\\n",
                    gt_type,desc,elem,valence,polarizability,sig_pol,vdwparams);
            nline++;
	  }
	  gt_old = gt_type;
	}
    }
  end_table(tp);
}

static void do_brule(FILE *tp,gmx_poldata_t pd,gmx_atomprop_t aps)
{
  int nline,npage,numbonds;
  char hbuf[1024],colinfo[1024],pbuf[1024];
  double length,sigma,bondorder;
  char *gt_name,*gt_type,*geometry,*neighbors;
  
  /* Bondtypes */
  strcpy(colinfo,"ccccc");
  strcpy(hbuf,"Rule & Type & Geometry & \\# Bonds & Neighbors");
  begin_table(tp,"Bonding rules defined in the Alexandria force field to determine atom types.",
              "brules",colinfo,hbuf);
  
  nline = 2;
  npage = 0;
  while (0 != gmx_poldata_get_bonding_rule(pd,&gt_name,&gt_type,&geometry,&numbonds,
                                           &neighbors)) {
    if (0 == (nline % maxline)) {
      end_table(tp);
      sprintf(pbuf,"brule%d",++npage);
      fprintf(tp,"\\addtocounter{table}{-1}\n");
      begin_table(tp,"Bonding rules, continued",pbuf,colinfo,hbuf);
      nline = 1;
    }
    fprintf(tp,"%s & %s & %s & %d & %s\\\\\n",gt_name,gt_type,geometry,numbonds,
            neighbors);
    nline++;
  }
  end_table(tp);
}
  
static void do_bad(FILE *fp,FILE *tp,gmx_poldata_t pd,gmx_atomprop_t aps)
{
  int bts[ebtsNR];
  int nline,npage;
  char *ai,*aj,*ak,*al,*params,*lu;
  char hbuf[1024],colinfo[1024],pbuf[1024];
  double length,ang,sigma,bondorder;
  
  lu = gmx_poldata_get_length_unit(pd);
  bts[ebtsBONDS] = gmx_poldata_get_bond_ftype(pd);
  bts[ebtsANGLES] = gmx_poldata_get_angle_ftype(pd);
  bts[ebtsPDIHS] = gmx_poldata_get_dihedral_ftype(pd,egdPDIHS);
  bts[ebtsIDIHS] = gmx_poldata_get_dihedral_ftype(pd,egdIDIHS);
  
  /* Bondtypes */
  strcpy(colinfo,"ccccl");
  sprintf(hbuf,"i & j & Length (%s) & Bond order & Params",lu);
  begin_table(tp,"Bonds defined in the Alexandria force field. Bond length (standard deviation in parentheses), bond order and Morse potential~\\protect\\cite{Morse29} parameters.",
              "btypes",colinfo,hbuf);
  
  fprintf(fp,"\n[ bondtypes ]\n");
  fprintf(fp,"; ; i    j  func       parameters\n");
  nline = 2;
  npage = 0;
  while (0 < gmx_poldata_get_bond(pd,&ai,&aj,&length,&sigma,&bondorder,&params)) {
    fprintf(fp,"%-5s  %-5s   %d  %g  %s\n",ai,aj,bts[ebtsBONDS],length,params);
    if (0 == (nline % maxline)) {
      end_table(tp);
      sprintf(pbuf,"btype%d",++npage);
      fprintf(tp,"\\addtocounter{table}{-1}\n");
      begin_table(tp,"Bonds, continued",pbuf,colinfo,hbuf);
      nline = 1;
    }
    fprintf(tp,"%s & %s & %.1f(%.1f) & %g & %s\\\\\n",ai,aj,length,sigma,bondorder,params);
    nline++;
  }
  end_table(tp);
  
  /* Angletypes */
  strcpy(colinfo,"ccccl");
  strcpy(hbuf,"i & j & k & Angle & Params");
  begin_table(tp,"Angles defined in the Alexandria force field.",
              "angtypes",colinfo,hbuf);
  
  fprintf(fp,"\n[ angletypes ]\n");
  fprintf(fp,"; ; i    j   k  func       parameters\n");
  nline = 2;
  npage = 0;
  while (0 < gmx_poldata_get_angle(pd,&ai,&aj,&ak,&ang,&sigma,&params)) {
    fprintf(fp,"%-5s  %-5s  %-5s  %d  %g  %s\n",ai,aj,ak,bts[ebtsANGLES],length,params);
    if (0 == (nline % maxline)) {
      end_table(tp);
      sprintf(pbuf,"angtype%d",++npage);
      fprintf(tp,"\\addtocounter{table}{-1}\n");
      begin_table(tp,"Angles, continued",pbuf,colinfo,hbuf);
      nline = 1;
    }
    fprintf(tp,"%s & %s & %s & %.2f(%.2f) & %s\\\\\n",ai,aj,ak,ang,sigma,params);
    nline++;
  }
  end_table(tp);
  
  /* Dihedraltypes */
  strcpy(colinfo,"cccccl");
  strcpy(hbuf,"i & j & k & l & Angle & Params");
  begin_table(tp,"Dihedrals defined in the Alexandria force field.",
              "dihtypes",colinfo,hbuf);
  
  fprintf(fp,"\n[ dihedraltypes ]\n");
  fprintf(fp,"; ; i    j   k    l  func       parameters\n");
  nline = 2;
  npage = 0;
  while (0 < gmx_poldata_get_dihedral(pd,egdPDIHS,&ai,&aj,&ak,&al,&ang,&sigma,&params)) {
    fprintf(fp,"%-5s  %-5s  %-5s  %-5s  %d  %.1f  %s\n",
            ai,aj,ak,al,bts[ebtsPDIHS],ang,params);
    if (0 == (nline % maxline)) {
      end_table(tp);
      sprintf(pbuf,"dihtype%d",++npage);
      fprintf(tp,"\\addtocounter{table}{-1}\n");
      begin_table(tp,"Dihedrals, continued",pbuf,colinfo,hbuf);
      nline = 1;
    }
    fprintf(tp,"%s & %s & %s & %s & %.1f(%.1f) & %s\\\\\n",ai,aj,ak,al,ang,sigma,params);
    nline++;
  }
  end_table(tp);
  
  /* Impropertypes */
  strcpy(colinfo,"cccccl");
  strcpy(hbuf,"i & j & k & l & Angle & Params");
  begin_table(tp,"Impropers defined in the Alexandria force field.",
              "idihtypes",colinfo,hbuf);
              
  fprintf(fp,"\n[ dihedraltypes ]\n");
  fprintf(fp,"; ; i    j   k    l  func       parameters\n");
  nline = 2;
  npage = 0;
  while (0 < gmx_poldata_get_dihedral(pd,egdIDIHS,&ai,&aj,&ak,&al,&ang,&sigma,&params)) {
    fprintf(fp,"%-5s  %-5s  %-5s  %-5s  %d  %.1f  %s\n",
            ai,aj,ak,al,bts[ebtsIDIHS],ang,params);
    if (0 == (nline % maxline)) {
      end_table(tp);
      sprintf(pbuf,"idihtype%d",++npage);
      fprintf(tp,"\\addtocounter{table}{-1}\n");
      begin_table(tp,"Impropers, continued",pbuf,colinfo,hbuf);
      nline = 1;
    }
    fprintf(tp,"%s & %s & %s & %s & %.1f(%.1f) & %s\\\\\n",ai,aj,ak,al,ang,sigma,params);
    nline++;
  }
  end_table(tp);
}

int main(int argc,char *argv[]) 
{
  gmx_poldata_t pd;
  gmx_atomprop_t aps;
  int numbonds,atomnumber,nline,maxline,npage;
  FILE *fp,*tp;
  
  aps = gmx_atomprop_init();
  pd  = gmx_poldata_read("gentop.dat",aps);
  tp = fopen("forcefield.tex","w");
  fprintf(tp,"%% Generated by gen_ff\n");
  fprintf(tp,"%% Copyright 2012 David van der Spoel\n");
    
  fp = fopen("forcefield.itp","w");
  fprintf(fp,"; This file is generated from gentop.dat by %s\n",argv[0]);
  fprintf(fp,"; Do not edit this file!\n");
  fprintf(fp,"; This is the force field file for the Alexandria FF\n");
  fprintf(fp,"; Paul J. van Maaren and David van der Spoel\n\n");
  
  do_atypes(fp,tp,pd,aps);  
  do_brule(tp,pd,aps);  
  do_bad(fp,tp,pd,aps);  
  
  fclose(fp);
  fclose(tp);  
  return 0;
}
