/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_pdbio_c = "$Id$";

#include "sysstuff.h"
#include "string2.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "symtab.h"
#include "assert.h"
#include "pdbio.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
	
static char *pdbtp[epdbNR]={"ATOM  ","HETATM"};
static char *pdbformat="%6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
/* This is THE format */
static bool bTER=FALSE;

void pdb_use_ter(bool bSet)
{
  bTER=bSet;
}

static void decompose_line(char line[],char type[],char anr[],char anm[],
			   char resnm[],char chain[],char resnr[],char xc[],
			   char yc[],char zc[],char dumc[],char bfac[],
			   char pdbresnr[])
{
  int j,k;
  char nc='\0';
  
  j=0;
  for(k=0; (k<6); k++,j++) type[k]=line[j];
  type[k]=nc;
  for(k=0; (k<5); k++,j++) anr[k]=line[j];
  anr[k]=nc;
  j++;
  for(k=0; (k<4); k++,j++) anm[k]=line[j];
  anm[k]=nc;
  j++;
  for(k=0; (k<4); k++,j++) 
    resnm[k]=line[j];
  resnm[k]=nc;
  /*j++;*/
  for(k=0; (k<1); k++,j++)
    chain[k]=line[j];
  chain[k]=nc;
  
  for(k=0; (k<4); k++,j++) {
    resnr[k]=line[j];
    pdbresnr[k]=line[j];
  }
  resnr[k]=nc;
  for( ;(k<8); k++,j++)
    pdbresnr[k]=line[j];
  pdbresnr[k]=nc;
  /*j+=4;*/
  for(k=0; (k<8); k++,j++) xc[k]=line[j];
  xc[k]=nc;
  for(k=0; (k<8); k++,j++) yc[k]=line[j];
  yc[k]=nc;
  for(k=0; (k<8); k++,j++) zc[k]=line[j];
  zc[k]=nc;
  for(k=0; (k<6); k++,j++) dumc[k]=line[j];
  dumc[k]=nc;
  for(k=0; (k<6); k++,j++) bfac[k]=line[j];
  bfac[k]=nc;
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

int read_pdbatoms(FILE *in,t_pdbatom **pdbaptr,matrix box,bool bFilterH)
{
  t_pdbatom *pdba=NULL;
  int  maxpdba=0;
  char line[STRLEN+1];
  char type[12],anr[12],anm[12],resnm[12],chain[12],resnr[12],
       xc[12],yc[12],zc[12],dumc[12],bfac[12],pdbresnr[12];
  int  natom;
  int  j,k;

  if (box != NULL) 
    clear_mat(box);
  
  natom=0;
  while (fgets2(line,STRLEN,in) != NULL) {
    decompose_line(line,type,anr,anm,resnm,chain,resnr,xc,yc,zc,dumc,bfac,
		   pdbresnr);
    for(j=0; (j<epdbNR); j++)
      if (strncmp(type,pdbtp[j],6) == 0) 
	break;
    if (j < epdbNR) {
      /* Should we filter out the hydrogens ? */
      if (bFilterH && (is_hydrogen(anm)))
	continue;
	
      if (natom >= maxpdba) {
	maxpdba+=256;
	srenew(pdba,maxpdba);
	/* Clear new entries! */
	for(k=natom; (k<maxpdba); k++)
	  memset(&(pdba[k]),0,sizeof(pdba[k]));
      }
      pdba[natom].pdbtp=j;
      pdba[natom].atomnr=atoi(anr);
      strcpy(pdba[natom].atomnm,anm);
      strcpy(pdba[natom].resnm,resnm);
      strcpy(pdba[natom].pdbresnr,pdbresnr);
      pdba[natom].chain=chain[0];
      pdba[natom].resnr=atoi(resnr);
      pdba[natom].x[XX]=atof(xc)*0.1;
      pdba[natom].x[YY]=atof(yc)*0.1;
      pdba[natom].x[ZZ]=atof(zc)*0.1;
      pdba[natom].bfac=atof(bfac);
      pdba[natom].dummy=atof(dumc);
      pdba[natom].m = 0.0;
      pdba[natom].q = 0.0;
      pdba[natom].type = 0;
      natom++;
    }
    else {
      trim(type);
      if (strncmp(type,"CRYST1",6)==0) {
	if (box != NULL) {
	  sscanf(line,"%*s%s%s%s",xc,yc,zc);
	  box[XX][XX] = atof(xc)*0.1;
	  box[YY][YY] = atof(yc)*0.1;
	  box[ZZ][ZZ] = atof(zc)*0.1;
	}
      }
      else
	if (bTER) {
	  if (strcmp(type,"TER") == 0) 
	    break;
	} else {
	  if (strcmp(type,"ENDMDL") == 0) 
	    break;
	}
    }
  }
  *pdbaptr=pdba;
  return natom;
}

void get_pdb_coordnum(char *infile,int *natoms)
{
  FILE *in;
  char line[STRLEN];
  int  j;
   
  in=ffopen(infile,"r");
  *natoms=0;
  while (fgets2(line,STRLEN,in)) {
    for(j=0; (j<epdbNR); j++)
      if (strncmp(line,pdbtp[j],6) == 0) {
	(*natoms)++;
	break;
      }
  }
  ffclose (in);
}

void print_pdbatoms(FILE *out,int natom,t_pdbatom pdba[],matrix box)
{
  int i,resnr;
  char buf[12];

  fprintf(out,"HEADER    %s\n",bromacs());
  if (box != NULL) {
    fprintf(out,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	    10*box[XX][XX],10*box[YY][YY],10*box[ZZ][ZZ],90.0,90.0,90.0);
  }
  for(i=0; (i<natom); i++) {
    sprintf(buf,"%s",pdba[i].resnm);
    buf[3]='\0';
    resnr=pdba[i].resnr + 1;
    if (resnr>=10000)
      resnr = resnr % 10000;
    fprintf(out,pdbformat,
	    pdbtp[pdba[i].pdbtp],pdba[i].atomnr + 1,pdba[i].atomnm,
	    buf,pdba[i].chain,resnr,
	    10*pdba[i].x[XX],10*pdba[i].x[YY],10*pdba[i].x[ZZ],
	    pdba[i].dummy,pdba[i].bfac);
  }
  fprintf(out,"TER\n");
}

void renumber_pdb(int natom,t_pdbatom pdba[])
{
  int i,nres=-1;
  char pdbnr[12],pdbnm[12];
  
  pdbnr[0]='\0';
  pdbnm[0]='\0';
  for(i=0; (i<natom); i++) {
    if ((strcmp(pdba[i].pdbresnr,pdbnr) != 0) ||
	(strcmp(pdba[i].resnm,pdbnm) != 0)) {
      strcpy(pdbnr,pdba[i].pdbresnr);
      strcpy(pdbnm,pdba[i].resnm);
      nres++;
    }
    pdba[i].atomnr=i;
    pdba[i].resnr=nres;
  }
}

void pdba_trimnames(int natom,t_pdbatom pdba[])
{
  int i;
  
  for(i=0; (i<natom); i++) {
    trim(pdba[i].atomnm);
    trim(pdba[i].resnm);
  }
}

int pdbasearch_atom(char *name,int resnr,int natom,t_pdbatom pdba[])
{
  int  i;
  char anmbuf[24];
  
  if (name[0] == '-') {
    name++;
    resnr--;
  }
  for(i=0; (i<natom) && (pdba[i].resnr != resnr); i++)
    ;
  for( ; (i<natom) && (pdba[i].resnr == resnr); i++) {
    strcpy(anmbuf,pdba[i].atomnm);
    trim(anmbuf);
    if (strcmp(name,pdba[i].atomnm) == 0)
      return i;
  }
  return -1;
}

void pdb2atoms(int natom,t_pdbatom pdba[],t_atoms *atoms,rvec **x,
	       t_symtab *symtab)
{
  int i,m,nres,rnr=-83;

  renumber_pdb(natom,pdba);
  nres=pdba[natom-1].resnr+1;
  
  snew(atoms->atom,natom);
  snew(atoms->atomname,natom);
  snew(atoms->resname,nres);
  snew(*x,natom);
  for(i=0; (i<natom); i++) {
    if (pdba[i].resnr != rnr) {
      rnr=pdba[i].resnr;
      atoms->resname[rnr]=put_symtab(symtab,pdba[i].resnm);
    }
    atoms->atom[i].resnr=pdba[i].resnr;
    atoms->atom[i].m=pdba[i].m;
    atoms->atom[i].q=pdba[i].q;
    atoms->atom[i].type=pdba[i].type;
    atoms->atomname[i]=put_symtab(symtab,pdba[i].atomnm);
    for(m=0; (m<DIM); m++)
      (*x)[i][m]=pdba[i].x[m];
  }
  assert(nres=rnr+1);
  atoms->nres=nres;
  atoms->nr=natom;
}

t_pdbatom *atoms2pdba(t_atoms *atoms,rvec x[])
{
  int       i;
  t_pdbatom *pdba;
  
  snew(pdba,atoms->nr);
  for(i=0; (i<atoms->nr); i++) {
    copy_rvec(x[i],pdba[i].x);
    pdba[i].pdbtp=epdbATOM;
    pdba[i].atomnr=i;
    pdba[i].resnr=atoms->atom[i].resnr;
    sprintf(pdba[i].pdbresnr,"%d",atoms->atom[i].resnr);
    strcpy(pdba[i].atomnm,*(atoms->atomname[i]));
    strcpy(pdba[i].resnm,*(atoms->resname[pdba[i].resnr]));
    pdba[i].chain=' ';
    pdba[i].m=atoms->atom[i].m;
    pdba[i].q=atoms->atom[i].q;
    pdba[i].type=atoms->atom[i].type;
  }
  
  return pdba;
}
