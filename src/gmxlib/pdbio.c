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
 * Green Red Orange Magenta Azure Cyan Skyblue
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
#include "physics.h"
#include "pbc.h"
	
static char *pdbtp[epdbNR]={
  "ATOM  ","HETATM", "ANISOU", "CRYST1",
  "COMPND", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK" 
};

static char *pdbformat ="%-6s%5u  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
static char *pdbformat4="%-6s%5u %-4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
static bool bTER=FALSE;
#define REMARK_SIM_BOX "REMARK    THIS IS A SIMULATION BOX"


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
			   bool bEndmodel,
			   atom_id nindex, atom_id index[])
{
  char resnm[6],nm[6],ch,*pdbform;
  atom_id i,ii;
  int  resnr,type;
  real occup,bfac;
  real alpha,beta,gamma;

  fprintf(out,"HEADER    %s\n",(title && title[0])?title:bromacs());
  if (box) {
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
      type=atoms->pdbinfo[i].type;
      occup=atoms->pdbinfo[i].occup;
      bfac=atoms->pdbinfo[i].bfac;
      if (occup==0.0)
	occup=1.0;
    }
    else {
      type=0;
      occup=1.0;
      bfac=0.0;
    }
    if (strlen(nm)==4)
      pdbform=pdbformat4;
    else
      pdbform=pdbformat;
    fprintf(out,pdbform,pdbtp[type],(i+1)%100000,nm,resnm,ch,resnr,
    10*x[i][XX],10*x[i][YY],10*x[i][ZZ],occup,bfac);
  }
  
  fprintf(out,"TER\n");
  if (bEndmodel && !bTER)
    fprintf(out,"ENDMDL\n");
}

void write_pdbfile(FILE *out,char *title,
		   t_atoms *atoms,rvec x[],matrix box,char chain,
		   bool bEndmodel)
{
  atom_id i,*index;

  snew(index,atoms->nr);
  for(i=0; i<atoms->nr; i++)
    index[i]=i;
  write_pdbfile_indexed(out,title,atoms,x,box,chain,bEndmodel,
			atoms->nr,index);
  sfree(index);
}

void hwrite_pdb_conf_indexed(FILE *out,char *title, 
			     t_atoms *atoms,rvec x[],matrix box,
			     atom_id nindex,atom_id index[])
{
  write_pdbfile_indexed(out,title,atoms,x,box,0,TRUE,nindex,index);
}

void write_pdb_confs(char *outfile,t_atoms **atoms,rvec *x[],int number)
{
  FILE *out;
  int n;
  char chain,str[STRLEN];

  out=ffopen(outfile,"w");
  
  for(n=0;(n<number);n++) {
    chain='A'+n;
    fprintf(stderr,"writing chain %c\n",chain);
    sprintf(str,"Chain %c",chain);
    write_pdbfile(out,str,atoms[n],x[n],NULL,chain,(n==number-1));
  }

  ffclose(out);
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
    fatal_error(0,"\nFound more atoms (%d) in pdb file than expected (%d)",
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
  for(k=0; (k<6); k++,j++) bfac[k]=line[j];
  bfac[k]=nc;

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
  atomn->m = 0.0;
  atomn->q = 0.0;
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

int read_pdbfile(FILE *in,char *title, 
		 t_atoms *atoms,rvec x[],matrix box,bool bChange)
{
  static t_symtab symtab;
  static bool bFirst=TRUE;
  bool bCOMPND;
  char line[STRLEN+1];
  char sa[12],sb[12],sc[12];
  double fa,fb,fc,alpha,beta,gamma;
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
      if (atoms->pdbinfo != NULL)
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
	  alpha *= DEG2RAD;
	  beta  *= DEG2RAD;
	  gamma *= DEG2RAD;
	  box[YY][XX] = fb*cos(gamma);
	  box[YY][YY] = fb*sin(gamma);
	  box[ZZ][XX] = fc*cos(beta);
	  box[ZZ][YY] = fc*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma);
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
  FILE      *in;
  
  in = ffopen(infile,"r");
  read_pdbfile(in,title,atoms,x,box,bChange);
  ffclose(in);
}
