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
	
static char *pdbtp[epdbNR]={
  "ATOM  ","HETATM", "ANISOU", "CRYST1",
  "COMPND", "ENDMDL", "TER", "HEADER", "TITLE", "REMARK" 
};

static char *pdbformat="%6s%5d  %-4.4s%3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
static char *pdbformat4="%6s%5d %-4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
/* This is THE format */
static char *anisouformat="%6s%5d %-4.4s %3.3s %c%4d  %-6d%-6d%-6d%-6d%-6d%-6d  %3s%-2s%2s\n";
static bool bTER=FALSE;
#define REMARK_SIM_BOX "REMARK    THIS IS A SIMULATION BOX"


void pdb_use_ter(bool bSet)
{
  bTER=bSet;
}

static void change_name(char *name)
{
  int i,length;
  char temp;
  bool bH;
  length=strlen(name);
  if (isdigit(name[length-1])&&isdigit(name[length-2])) {
    bH=FALSE;
    for (i=0;(i<length);i++)
      if (name[i] == 'H') 
        bH=TRUE;
    if (bH) {
      temp=name[length-1]; 
      for(i=length-1;(i>0);i--)
	name[i]=name[i-1];
      name[0]=temp;
    }
  }
  else {
    if(strcmp(name,"O2")==0)
    strcpy(name,"OXT");
  }
}

void write_pdbfile_indexed(FILE *out,char *title,
			   t_atoms *atoms,rvec x[],matrix box,char chain,
			   bool bEndmodel,
			   atom_id nindex, atom_id index[])
{
  char resnm[6],nm[6],ch,*pdbform;
  atom_id i,ii;
  int  resnr,type;
  real occup;

  fprintf(out,"HEADER    %s\n",title[0]?title:bromacs());
  if (box) {
    fprintf(out,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n",
	    10*box[XX][XX],10*box[YY][YY],10*box[ZZ][ZZ],90.0,90.0,90.0);
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
    occup=atoms->atom[i].occup;
    if (occup==0.0)
      occup=1.0;
    if (strlen(nm)==4)
      pdbform=pdbformat4;
    else
      pdbform=pdbformat;
    if (atoms->pdbinfo != NULL)
      type=atoms->pdbinfo[i].type;
    else
      type=0;
    fprintf(out,pdbform,
	    pdbtp[type],i+1,nm,resnm,ch,resnr,
	    10*x[i][XX],10*x[i][YY],10*x[i][ZZ],occup,atoms->atom[i].bfac);
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

extern void write_pdb_conf(char *outfile,char *title,
			   t_atoms *atoms,rvec x[],matrix box)
{
  FILE *out;
  out=ffopen(outfile,"w");
  write_pdbfile(out,title,atoms,x,box,0,TRUE);
  fclose(out);
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

static int read_atom(char line[],int natom,
		     t_atoms *atoms,rvec x[],bool bChange)
{
  static t_symtab symtab;
  t_atom *atomn;
  int  j,k;
  char nc='\0';
  char anr[12],anm[12],resnm[12],chain[12],resnr[12];
  char xc[12],yc[12],zc[12],occup[12],bfac[12],pdbresnr[12];
  static int oldres;
  int  resnri,newres;

  open_symtab(&symtab);

  /* Skip over type */  
  j=6;
  for(k=0; (k<5); k++,j++) anr[k]=line[j];
  anr[k]=nc;
  j++;
  for(k=0; (k<4); k++,j++) anm[k]=line[j];
  anm[k]=nc;
  j++;
  for(k=0; (k<4); k++,j++) 
    resnm[k]=line[j];
  resnm[k]=nc;

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
  if (atoms->pdbinfo != NULL) {
    atoms->pdbinfo[natom].type=j;
    atoms->pdbinfo[natom].atomnr=atoi(anr);
    strcpy(atoms->pdbinfo[natom].pdbresnr,pdbresnr);
  }
  resnri=atoi(resnr);
  if ((natom==0) || (resnri != oldres)) {
    oldres=resnri;
    if (natom==0)
      newres=0;
    else
      newres=atoms->atom[natom-1].resnr+1;
    atoms->nres=newres+1;
    atoms->resname[newres]=put_symtab(&symtab,resnm);
  }
  else
    newres=atoms->atom[natom-1].resnr;
  if (bChange)
    change_name(anm); 
  atoms->atomname[natom]=put_symtab(&symtab,anm);
  atomn->chain=chain[0];
  atomn->resnr=newres;
  x[natom][XX]=atof(xc)*0.1;
  x[natom][YY]=atof(yc)*0.1;
  x[natom][ZZ]=atof(zc)*0.1;
  atomn->bfac=atof(bfac);
  atomn->occup=atof(occup);
  atomn->m = 0.0;
  atomn->q = 0.0;
  natom++;
  
  close_symtab(&symtab);

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

int read_pdbfile(FILE *in,char *title, 
		 t_atoms *atoms,rvec x[],matrix box,bool bChange)
{
  bool bCOMPND,bSimBox;
  char line[STRLEN+1];
  char xc[12],yc[12],zc[12];
  int  line_type;
  char *c,*d;
  int  natom;
  int  j,k;

  if (box != NULL) 
    clear_mat(box);

  bCOMPND=FALSE;
  bSimBox=FALSE;
  title[0]='\0';
  natom=0;
  while (fgets2(line,STRLEN,in) != NULL) {
    line_type = line2type(line);
    
    switch(line_type) {
    case epdbATOM:
    case epdbHETATM:
      natom = read_atom(line,natom,atoms,x,bChange);
      break;
      
    case epdbANISOU:
      if (atoms->pdbinfo != NULL)
	read_anisou(line,natom,atoms);
      break;

    case epdbCRYST1:      
      if (!box[0][0] && bSimBox) {
	sscanf(line,"%*s%s%s%s",xc,yc,zc);
	box[XX][XX] = atof(xc)*0.1;
	box[YY][YY] = atof(yc)*0.1;
	box[ZZ][ZZ] = atof(zc)*0.1;
	
	bSimBox = FALSE;
      } else 
	fprintf(stderr,"WARNING: ignoring data in CRYST1 entry "
		"(probably not a simulation box)\n");
      break;

    case epdbTITLE:
    case epdbHEADER:      
      c=line+6;
      /* skip HEADER or TITLE and spaces */
      while (c && (c[0]!=' ')) c++;
      while (c && (c[0]==' ')) c++;
      /* truncate after title */
      d=strstr(c,"   ");
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
	if (strlen(c)>0)
	  if (bCOMPND) {
	    strcat(title,"; ");
	    strcat(title,c);
	  } else
	    strcpy(title,c);
	bCOMPND=TRUE;
      } 
      break;
      
    case epdbREMARK:
      if (strcmp(line,REMARK_SIM_BOX)==0)
	bSimBox=TRUE;
      break;
      
    case epdbTER:
    case epdbENDMDL:
    default:
      break;
    }
  }
  
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
    if ( ( bTER && (strncmp(line,"TER",3) == 0)) ||
	 (!bTER && (strncmp(line,"ENDMDL",6) == 0)) ) 
      break;
    if ((strncmp(line,"ATOM  ",6) == 0) || (strncmp(line,"HETATM",6) == 0))
      (*natoms)++;
  }
  ffclose (in);
}

void read_pdb_conf(char *infile,char *title, 
		   t_atoms *atoms,rvec x[],matrix box,bool bChange)
{
  FILE      *in;
  
  in = ffopen(infile,"r");
  read_pdbfile(in,title,atoms,x,box,bChange);
  ffclose(in);
}

void print_pdbatoms(FILE *out,char *title, 
		    int natom,t_pdbatom pdba[],matrix box)
{
  int i,resnr;
  char buf[12];

  if (title && (title[0] != '\0'))
    fprintf(out,"HEADER    %s\n",title);
  else
    fprintf(out,"HEADER    %s\n",bromacs());
  if (box != NULL) {
    fprintf(out,REMARK_SIM_BOX"\n");
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
    if (pdba[i].bAnisotropic)
      fprintf(out,anisouformat,
	      pdbtp[pdba[i].pdbtp],pdba[i].atomnr + 1,pdba[i].atomnm,
	      buf,pdba[i].chain,resnr,
	      pdba[i].uij[U11],pdba[i].uij[U22],pdba[i].uij[U33],
	      pdba[i].uij[U12],pdba[i].uij[U13],pdba[i].uij[U23],
	      "A","C","0");
	      
  }
  fprintf(out,"%s\n",bTER?"TER":"ENDMDL");
}

void renumber_pdb(int natom,t_pdbatom pdba[])
{
  int i,nres;
  char pdbnr[12],pdbnm[12];
  
  pdbnr[0]='\0';
  pdbnm[0]='\0';
  nres=-1;
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
  snew(atoms->chain,natom);
  snew(*x,natom);
  for(i=0; (i<natom); i++) {
    if (pdba[i].resnr != rnr) {
      rnr=pdba[i].resnr;
      atoms->resname[rnr]=put_symtab(symtab,pdba[i].resnm);
    }
    atoms->atom[i].resnr=pdba[i].resnr;
    atoms->chain[i]=pdba[i].chain;
    atoms->atom[i].m=pdba[i].m;
    atoms->atom[i].q=pdba[i].q;
    atoms->atom[i].type=pdba[i].type;
    sfree(atoms->atomname[i]);
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
    if (atoms->chain)
      pdba[i].chain=atoms->atom[i].chain;
    else
      pdba[i].chain=' ';
    pdba[i].m=atoms->atom[i].m;
    pdba[i].q=atoms->atom[i].q;
    pdba[i].type=atoms->atom[i].type;
  }
  
  return pdba;
}
