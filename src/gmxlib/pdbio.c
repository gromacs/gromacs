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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>

#include "sysstuff.h"
#include "string2.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "symtab.h"
#include "pdbio.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "physics.h"
#include "pbc.h"
	
static const char *pdbtp[epdbNR]={
  "ATOM  ","HETATM", "ANISOU", "CRYST1",
  "COMPND", "MODEL", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK" 
};

static bool bTER=FALSE;
static bool bWideFormat=FALSE;
#define REMARK_SIM_BOX "REMARK    THIS IS A SIMULATION BOX"

void set_pdb_wide_format(bool bSet)
{
  bWideFormat = bSet;
}

void pdb_use_ter(bool bSet)
{
  bTER=bSet;
}

static void gromacs_name(char *name)
{
  int i,length;
  char temp;

  length=strlen(name);
  if (isdigit(name[0])) {
    temp=name[0]; 
    for(i=1; i<length; i++)
      name[i-1]=name[i];
    name[length-1]=temp;
  }
  if(strcmp(name,"OXT")==0)
    strcpy(name,"O2");
}

void write_pdbfile_indexed(FILE *out,char *title,
			   t_atoms *atoms,rvec x[],matrix box,char chain,
			   int model_nr, atom_id nindex, atom_id index[])
{
  char resnm[6],nm[6],ch,pdbform[128],pukestring[100];
  atom_id i,ii;
  int  resnr,type;
  real occup,bfac;
  real alpha,beta,gamma;
  bool bOccup;

  fprintf(out,"TITLE     %s\n",(title && title[0])?title:bromacs(pukestring,99));
  if (bWideFormat) {
    fprintf(out,"REMARK    This file does not adhere to the PDB standard\n");
    fprintf(out,"REMARK    As a result of, some programs may not like it\n");
  }
  if (box && ( norm2(box[XX]) || norm2(box[YY]) || norm2(box[ZZ]) ) ) {
    if (norm2(box[YY])*norm2(box[ZZ])!=0)
      alpha = RAD2DEG*acos(cos_angle_no_table(box[YY],box[ZZ]));
    else
      alpha = 90;
    if (norm2(box[XX])*norm2(box[ZZ])!=0)
      beta  = RAD2DEG*acos(cos_angle_no_table(box[XX],box[ZZ]));
    else
      beta  = 90;
    if (norm2(box[XX])*norm2(box[YY])!=0)
      gamma = RAD2DEG*acos(cos_angle_no_table(box[XX],box[YY]));
    else
      gamma = 90;
    fprintf(out,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	    10*norm(box[XX]),10*norm(box[YY]),10*norm(box[ZZ]),
	    alpha,beta,gamma);
  }
  if (atoms->pdbinfo) {
    /* Check whether any occupancies are set, in that case leave it as is,
     * otherwise set them all to one
     */
    bOccup = TRUE;
    for (ii=0; (ii<nindex) && bOccup; ii++) {
      i      = index[ii];
      bOccup = bOccup && (atoms->pdbinfo[i].occup == 0.0);
    }
  } 
  else
    bOccup = FALSE;

  if (!bTER)
    fprintf(out,"MODEL %8d\n",model_nr>=0 ? model_nr : 1);
  for (ii=0; ii<nindex; ii++) {
    i=index[ii];
    resnr=atoms->atom[i].resnr;
    strcpy(resnm,*atoms->resname[resnr]);
    strcpy(nm,*atoms->atomname[i]);
    resnr++;
    if (resnr>=10000)
      resnr = resnr % 10000;
    if (chain)
      ch=chain;
    else
      if (atoms->atom[i].chain)
	ch=atoms->atom[i].chain;
      else
	  ch=' ';
    if (atoms->pdbinfo) {
      type  = atoms->pdbinfo[i].type;
      occup = bOccup ? 1.0 : atoms->pdbinfo[i].occup;
      bfac  = atoms->pdbinfo[i].bfac;
    }
    else {
      type  = 0;
      occup = 1.0;
      bfac  = 0.0;
    }
    if (bWideFormat)
      strcpy(pdbform,"%-6s%5u %-4.4s %3.3s %c%4d    %10.5f%10.5f%10.5f%8.4f%8.4f\n");
    else {
      if (strlen(nm)<4)
	strcpy(pdbform,pdbformat);
      else {
	strcpy(pdbform,pdbformat4);
	if (strlen(nm)>4)
	  fprintf(stderr,"WARNING: Writing out atom name (%s) longer than 4 characters to .pdb file\n",nm);
      }
      strcat(pdbform,"%6.2f%6.2f\n");
    }
    fprintf(out,pdbform,pdbtp[type],(i+1)%100000,nm,resnm,ch,resnr,
	    10*x[i][XX],10*x[i][YY],10*x[i][ZZ],occup,bfac);
    if (atoms->pdbinfo && atoms->pdbinfo[i].bAnisotropic) {
      fprintf(out,"ANISOU%5u  %-4.4s%3.3s %c%4d  %7d%7d%7d%7d%7d%7d\n",
	      (i+1)%100000,nm,resnm,ch,resnr,
	      atoms->pdbinfo[i].uij[0],atoms->pdbinfo[i].uij[1],
	      atoms->pdbinfo[i].uij[2],atoms->pdbinfo[i].uij[3],
	      atoms->pdbinfo[i].uij[4],atoms->pdbinfo[i].uij[5]);
    }
  }
 
  fprintf(out,"TER\n");
  if (!bTER)
    fprintf(out,"ENDMDL\n");
}

void write_pdbfile(FILE *out,char *title, t_atoms *atoms,rvec x[],
		   matrix box,char chain,int model_nr)
{
  atom_id i,*index;

  snew(index,atoms->nr);
  for(i=0; i<atoms->nr; i++)
    index[i]=i;
  write_pdbfile_indexed(out,title,atoms,x,box,chain,model_nr,atoms->nr,index);
  sfree(index);
}

int line2type(char *line)
{
  int  k;
  char type[8];
  
  for(k=0; (k<6); k++) 
    type[k]=line[k];
  type[k]='\0';
  
  for(k=0; (k<epdbNR); k++)
    if (strncmp(type,pdbtp[k],strlen(pdbtp[k])) == 0)
      break;
  
  return k;
}

static void read_anisou(char line[],int natom,t_atoms *atoms)
{
  int  i,j,k,atomnr;
  char nc='\0';
  char anr[12],anm[12];

  /* Skip over type */  
  j=6;
  for(k=0; (k<5); k++,j++) anr[k]=line[j];
  anr[k]=nc;
  j++;
  for(k=0; (k<4); k++,j++) anm[k]=line[j];
  anm[k]=nc;
  j++;
  
  /* Strip off spaces */
  trim(anm);
  
  /* Search backwards for number and name only */
  atomnr = atoi(anr);
  for(i=natom-1; (i>=0); i--)
    if ((strcmp(anm,*(atoms->atomname[i])) == 0) && 
	(atomnr == atoms->pdbinfo[i].atomnr))
      break;
  if (i < 0)
    fprintf(stderr,"Skipping ANISOU record (atom %s %d not found)\n",
	    anm,atomnr);
  else {
    if (sscanf(line+29,"%d%d%d%d%d%d",
	       &atoms->pdbinfo[i].uij[U11],&atoms->pdbinfo[i].uij[U22],
	       &atoms->pdbinfo[i].uij[U33],&atoms->pdbinfo[i].uij[U12],
	       &atoms->pdbinfo[i].uij[U13],&atoms->pdbinfo[i].uij[U23])
	         == 6) {
      atoms->pdbinfo[i].bAnisotropic = TRUE;
    }
    else {
      fprintf(stderr,"Invalid ANISOU record for atom %d\n",i);
      atoms->pdbinfo[i].bAnisotropic = FALSE;
    }     
  }
}

static int read_atom(t_symtab *symtab,char line[],int type,int natom,
		     t_atoms *atoms,rvec x[],bool bChange)
{
  t_atom *atomn;
  int  j,k;
  char nc='\0';
  char anr[12],anm[12],altloc,resnm[12],chain[12],resnr[12];
  char xc[12],yc[12],zc[12],occup[12],bfac[12],pdbresnr[12];
  static char oldresnm[12],oldresnr[12];
  int  newres;

  if (natom>=atoms->nr)
    gmx_fatal(FARGS,"\nFound more atoms (%d) in pdb file than expected (%d)",
	      natom+1,atoms->nr);

  /* Skip over type */  
  j=6;
  for(k=0; (k<5); k++,j++) anr[k]=line[j];
  anr[k]=nc;
  trim(anr);
  j++;
  for(k=0; (k<4); k++,j++) anm[k]=line[j];
  anm[k]=nc;
  trim(anm);
  altloc=line[j];
  j++;
  for(k=0; (k<4); k++,j++) 
    resnm[k]=line[j];
  resnm[k]=nc;
  trim(resnm);

  for(k=0; (k<1); k++,j++)
    chain[k]=line[j];
  chain[k]=nc;
  
  for(k=0; (k<4); k++,j++) {
    resnr[k]=line[j];
    pdbresnr[k]=line[j];
  }
  resnr[k]=nc;
  trim(resnr);
  pdbresnr[k]=line[j];
  pdbresnr[k+1]=nc;
  trim(pdbresnr);
  j+=4;

  /* X,Y,Z Coordinate */
  for(k=0; (k<8); k++,j++) xc[k]=line[j];
  xc[k]=nc;
  for(k=0; (k<8); k++,j++) yc[k]=line[j];
  yc[k]=nc;
  for(k=0; (k<8); k++,j++) zc[k]=line[j];
  zc[k]=nc;
  
  /* Weight */
  for(k=0; (k<6); k++,j++) occup[k]=line[j];
  occup[k]=nc;
  
  /* B-Factor */
  for(k=0; (k<7); k++,j++) bfac[k]=line[j];
  bfac[k]=nc;

  if (atoms->atom) {
    atomn=&(atoms->atom[natom]);
    if ((natom==0) || (strcmp(oldresnr,pdbresnr)!=0) || 
	(strcmp(oldresnm,resnm)!=0)) {
      strcpy(oldresnr,pdbresnr);
      strcpy(oldresnm,resnm);
      if (natom==0)
	newres=0;
      else
	newres=atoms->atom[natom-1].resnr+1;
      atoms->nres=newres+1;
      atoms->resname[newres]=put_symtab(symtab,resnm);
    }
    else
      newres=atoms->atom[natom-1].resnr;
    if (bChange)
      gromacs_name(anm); 
    atoms->atomname[natom]=put_symtab(symtab,anm);
    atomn->chain=chain[0];
    atomn->resnr=newres;
    atomn->m = 0.0;
    atomn->q = 0.0;
  }
  x[natom][XX]=atof(xc)*0.1;
  x[natom][YY]=atof(yc)*0.1;
  x[natom][ZZ]=atof(zc)*0.1;
  if (atoms->pdbinfo) {
    atoms->pdbinfo[natom].type=type;
    atoms->pdbinfo[natom].atomnr=atoi(anr);
    atoms->pdbinfo[natom].altloc=altloc;
    strcpy(atoms->pdbinfo[natom].pdbresnr,pdbresnr);
    atoms->pdbinfo[natom].bfac=atof(bfac);
    atoms->pdbinfo[natom].occup=atof(occup);
  }
  natom++;
  
  return natom;
}

bool is_hydrogen(char *nm)
{
  char buf[30];
  
  strcpy(buf,nm);
  trim(buf);
  
  if (buf[0] == 'H')
    return TRUE;
  else if ((isdigit(buf[0])) && (buf[1] == 'H'))
    return TRUE;
  return FALSE;
}

bool is_dummymass(char *nm)
{
  char buf[30];
  
  strcpy(buf,nm);
  trim(buf);
  
  if ((buf[0] == 'M') && isdigit(buf[strlen(buf)-1]))
    return TRUE;
      
  return FALSE;
}

int read_pdbfile(FILE *in,char *title,int *model_nr,
		 t_atoms *atoms,rvec x[],matrix box,bool bChange)
{
  static t_symtab symtab;
  static bool bFirst=TRUE;
  bool bCOMPND;
  char line[STRLEN+1];
  char sa[12],sb[12],sc[12];
  double fa,fb,fc,alpha,beta,gamma,cosa,cosb,cosg,sing;
  int  line_type;
  char *c,*d;
  int  natom;
  bool bStop=FALSE;

  if (box != NULL) 
    clear_mat(box);

  if (bFirst) {
    open_symtab(&symtab);
    bFirst=FALSE;
  }

  bCOMPND=FALSE;
  title[0]='\0';
  natom=0;
  while (!bStop && (fgets2(line,STRLEN,in) != NULL)) {
    line_type = line2type(line);
    
    switch(line_type) {
    case epdbATOM:
    case epdbHETATM:
      natom = read_atom(&symtab,line,line_type,natom,atoms,x,bChange);
      break;
      
    case epdbANISOU:
      if (atoms->pdbinfo)
	read_anisou(line,natom,atoms);
      break;

    case epdbCRYST1:      
      if (box) {
	sscanf(line,"%*s%s%s%s%lf%lf%lf",sa,sb,sc,&alpha,&beta,&gamma);
	fa = atof(sa)*0.1;
	fb = atof(sb)*0.1;
	fc = atof(sc)*0.1;
	clear_mat(box);
	box[XX][XX] = fa;
	if ((alpha!=90.0) || (beta!=90.0) || (gamma!=90.0)) {
	  if (alpha != 90.0) {
	    cosa = cos(alpha*DEG2RAD);
	  } else {
	    cosa = 0;
	  }
	  if (beta != 90.0) {
	    cosb = cos(beta*DEG2RAD);
	  } else {
	    cosb = 0;
	  }
	  if (gamma != 90.0) {
	    cosg = cos(gamma*DEG2RAD);
	    sing = sin(gamma*DEG2RAD);
	  } else {
	    cosg = 0;
	    sing = 1;
	  }
	  box[YY][XX] = fb*cosg;
	  box[YY][YY] = fb*sing;
	  box[ZZ][XX] = fc*cosb;
	  box[ZZ][YY] = fc*(cosa - cosb*cosg)/sing;
	  box[ZZ][ZZ] = sqrt(fc*fc
			     -box[ZZ][XX]*box[ZZ][XX]-box[ZZ][YY]*box[ZZ][YY]);
	} else {
	  box[YY][YY] = fb;
	  box[ZZ][ZZ] = fc;
	}
      }
      break;

    case epdbTITLE:
    case epdbHEADER:      
      c=line+6;
      /* skip HEADER or TITLE and spaces */
      while (c && (c[0]!=' ')) c++;
      while (c && (c[0]==' ')) c++;
      /* truncate after title */
      d=strstr(c,"      ");
      if (d) {
	d[0]='\0';
      }
      if (strlen(c)>0)
	strcpy(title,c);
      break;
      
    case epdbCOMPND:
      if ((!strstr(line,": ")) || (strstr(line+6,"MOLECULE:"))) {
	if ( !(c=strstr(line+6,"MOLECULE:")) )
	  c=line;
	/* skip 'MOLECULE:' and spaces */
	while (c && (c[0]!=' ')) c++;
	while (c && (c[0]==' ')) c++;
	/* truncate after title */
	d=strstr(c,"   ");
	if (d) {
	  while ( (d[-1]==';') && d>c)  d--;
	  d[0]='\0';
	}
	if (strlen(c) > 0) {
	  if (bCOMPND) {
	    strcat(title,"; ");
	    strcat(title,c);
	  } else
	    strcpy(title,c);
	}
	bCOMPND=TRUE;
      } 
      break;
      
    case epdbTER:
      if (bTER)
	bStop=TRUE;
      break;
    case epdbMODEL:
      if(model_nr)
	sscanf(line,"%*s%d",model_nr);
      break;
    case epdbENDMDL:
      bStop=TRUE;
      break;
    default:
      break;
    }
  }
  
  return natom;
}

void get_pdb_coordnum(FILE *in,int *natoms)
{
  char line[STRLEN];
   
  *natoms=0;
  while (fgets2(line,STRLEN,in)) {
    if ( ( bTER && (strncmp(line,"TER",3) == 0)) ||
	 (!bTER && (strncmp(line,"ENDMDL",6) == 0)) ) 
      break;
    if ((strncmp(line,"ATOM  ",6) == 0) || (strncmp(line,"HETATM",6) == 0))
      (*natoms)++;
  }
}

void read_pdb_conf(char *infile,char *title, 
		   t_atoms *atoms,rvec x[],matrix box,bool bChange)
{
  FILE *in;
  
  in = ffopen(infile,"r");
  read_pdbfile(in,title,NULL,atoms,x,box,bChange);
  ffclose(in);
}
