/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Great Red Oystrich Makes All Chemists Sane
 */

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "errno.h"
#include "macros.h"
#include "string2.h"
#include "confio.h"
#include "vec.h"
#include "symtab.h"
#include "assert.h"
#include "futil.h"
#include "xdrf.h"
#include "pdbio.h"
#include "fatal.h"
#include "copyrite.h"

static void get_coordnum_fp (FILE *in,char *title, int *natoms)
{
  char line[STRLEN+1];

  fgets2 (title,STRLEN,in);
  fgets2 (line,STRLEN,in);
  sscanf (line,"%d",natoms);
}

void get_coordnum (char *infile,int *natoms)
{
  FILE *in;
  char title[STRLEN];
  
  in=ffopen(infile,"r");
  get_coordnum_fp(in,title,natoms);
  ffclose (in);
}

static void get_w_conf(FILE *in, char *infile, char *title,
		       t_atoms *atoms, rvec x[],rvec v[], matrix box)
{
  static t_symtab symtab;
  char   name[6];
  char   line[STRLEN+1];
  char   buf[256];
  char   format[30];
  double x1,y1,z1,x2,y2,z2;
  rvec   xmin,xmax;
  int    i,m,resnr,newres,oldres,prec;
  bool   bFirst;
  char   *p1,*p2;
  
  newres=0;
  oldres=0;
  
  open_symtab(&symtab);
  name[5]=0;

  fgets2(title,STRLEN,in);

  /* read the number of atoms */
  fgets2(line,STRLEN,in);
  sscanf(line,"%d",&(atoms->nr));

  atoms->nres=0;

  bFirst=TRUE;
  /* just pray the arrays are big enough */
  for (i=0; (i < atoms->nr) ; i++) {
    if ((fgets2 (line,STRLEN,in)) == NULL) {
      unexpected_eof(infile,i+2);
    }
    
    /* determine read precision from distance between periods 
       (decimal points) */
    if (bFirst) {
      bFirst=FALSE;
      p1=strchr(line,'.');
      p2=strchr(&p1[1],'.');
      if (p1 || p2)
	prec=p2-p1;
      else
	prec=8;
      if (prec<0)
	prec=8;
      if (prec>30)
	prec=30;
      sprintf(format,"%%%dlf%%%dlf%%%dlf",prec,prec,prec);
    }
    
    /* residue number*/
    memcpy(name,line,5);
    sscanf(name,"%d",&resnr);
    memcpy(name,line+5,5);
    name[5]='\0';
    if (resnr != oldres) {
      oldres=resnr;
      newres++;
      atoms->nres=newres;
      atoms->resname[newres-1]=put_symtab(&symtab,name);
    }
    resnr=newres;
    atoms->atom[i].resnr = resnr-1;

    /* atomname */
    memcpy(name,line+10,5);
    atoms->atomname[i]=put_symtab(&symtab,name);
   
    /* eventueel controle atomnumber met i+1 */

    /* coordinates (start after residue shit) */
    if (sscanf (line+20,format,&x1,&y1,&z1) != 3) {
      too_few();
    }
    else {
      x[i][XX]=x1;
      x[i][YY]=y1;

            x[i][ZZ]=z1;
    }

    /* velocities (start after residues and coordinates) */
    if (sscanf (line+20+(3*prec),format,&x1,&y1,&z1) != 3) {
      v[i][XX] = 0.0;
      v[i][YY] = 0.0;
      v[i][ZZ] = 0.0;
    }
    else {
      v[i][XX]=x1;
      v[i][YY]=y1;
      v[i][ZZ]=z1;
    }
  }

  /* box */
  fgets2 (line,STRLEN,in);
  if (sscanf (line,"%lf%lf%lf",&x1,&y1,&z1) != 3) {
    sprintf(buf,"Bad box in file %s",infile);
    warning(buf);
    
    /* Generate a cubic box */
    for(m=0; (m<DIM); m++)
      xmin[m]=xmax[m]=x[0][m];
    for(i=1; (i<atoms->nr); i++)
      for(m=0; (m<DIM); m++) {
	xmin[m]=min(xmin[m],x[i][m]);
	xmax[m]=max(xmax[m],x[i][m]);
      }
    for (i=0; i<DIM; i++) for (m=0; m<DIM; m++) box[i][m]=0.0;
    for(m=0; (m<DIM); m++)
      box[m][m]=(xmax[m]-xmin[m]);
    fprintf(stderr,"Generated a cubic box %8.3f x %8.3f x %8.3f\n",
	    box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
  }
  else {
    /* We found the first three values, the diagonal elements */
    box[XX][XX]=x1;
    box[YY][YY]=y1;
    box[ZZ][ZZ]=z1;
    if (sscanf (line,"%*f%*f%*f%lf%lf%lf%lf%lf%lf",
		&x1,&y1,&z1,&x2,&y2,&z2) != 6) 
      x1=y1=z1=x2=y2=z2=0.0;
    box[XX][YY] = x1;
    box[XX][ZZ] = y1;
    box[YY][XX] = z1;
    box[YY][ZZ] = x2;
    box[ZZ][XX] = y2;
    box[ZZ][YY] = z2;
  }
  close_symtab(&symtab);
}

void read_whole_conf(char *infile, char *title,
		     t_atoms *atoms, rvec x[],rvec v[], matrix box)
{
  FILE   *in;
  
  /* open file */
  in=ffopen(infile,"r");

  get_w_conf(in, infile, title, atoms, x, v, box);
  
  fclose(in);
}

void read_conf(char *infile,char *title,int *natoms,
	       rvec x[],rvec v[],matrix box)
{
  t_atoms  atoms;

  atoms.nr=*natoms;
  snew(atoms.atom,*natoms);
  atoms.nres=*natoms;
  snew(atoms.resname,*natoms);
  snew(atoms.atomname,*natoms);
  
  read_whole_conf(infile,title,&atoms,x,v,box);
  
  sfree(atoms.atom);
  sfree(atoms.resname);
  sfree(atoms.atomname);
}

static void get_conf(FILE *in, char *title, int *natoms, 
		     rvec x[],rvec v[],matrix box)
{
  t_atoms  atoms;

  atoms.nr=*natoms;
  snew(atoms.atom,*natoms);
  atoms.nres=*natoms;
  snew(atoms.resname,*natoms);
  snew(atoms.atomname,*natoms);
  
  get_w_conf(in,title,title,&atoms,x,v,box);
  
  sfree(atoms.atom);
  sfree(atoms.resname);
  sfree(atoms.atomname);
}

bool gro_next_x_or_v(FILE *status,real *t,int natoms,
		     rvec x[],rvec v[],matrix box)
{
  char title[STRLEN],*p;
  
  if (eof(status))
    return FALSE;
    
  get_conf(status, title, &natoms, x, v, box);

  if (p=strstr(title,"t=")) {
    p+=2;
    if (sscanf(p,"%f",t)!=1)
      *t=0.0;
  }
  return TRUE;
}

int gro_first_x_or_v(FILE *status, real *t, 
		     rvec **x, rvec **v, matrix box)
{
  int natoms;
  char title[STRLEN];
  
  *t=0.0;
  frewind(status);
  get_coordnum_fp(status, title, &natoms);
  frewind(status);
  fprintf(stderr,"Reading frames from gro file '%s'\n",title);
  if (natoms==0)
    fatal_error(1,"No coordinates in gro file\n");
  fprintf(stderr,"No of atoms: %d.\n", natoms);
  
  snew(*x,natoms);
  snew(*v,natoms);
  gro_next_x_or_v(status, t, natoms, *x, *v, box);
  
  return natoms;
}

bool gro_next_v(FILE *status,real *t,int natoms,rvec v[],matrix box)
{
  int i,d;
  bool result,vel;
  rvec *x;
  
  snew(x,natoms);
  do {
    result = gro_next_x_or_v(status,t,natoms,x,v,box);
    vel=FALSE;
    for (i=0; (i<natoms); i++)
      for (d=0; (d<DIM); d++)
	vel=vel || v[d][XX];
  } while(result && !vel);
  sfree(x);
  
  return result;
}

int gro_first_v(FILE *status, real *t, rvec **v, matrix box)
{
  bool result,vel;
  rvec *x;
  int i,d,natoms;
  
  natoms = gro_first_x_or_v(status,t,&x,v,box);
  result = natoms;
  vel=FALSE;
  for (i=0; (i<natoms); i++)
    for (d=0; (d<DIM); d++)
      vel=vel || v[d][XX];
  while(result && !vel) {
    result = gro_next_x_or_v(status,t,natoms,x,*v,box);
    vel=FALSE;
    for (i=0; (i<natoms); i++)
      for (d=0; (d<DIM); d++)
	vel=vel || v[d][XX];
  }
  return result;
}


bool gro_next_x(FILE *status,real *t,int natoms,rvec x[],matrix box)
{
  rvec *v;
  bool result;
  
  snew(v,natoms);
  result = gro_next_x_or_v(status,t,natoms,x,v,box);
  sfree(v);
  
  return result;
}

int gro_first_x(FILE *status, real *t, rvec **x, matrix box)
{
  rvec *v;
  int result;
  
  result = gro_first_x_or_v(status,t,x,&v,box);
  sfree(v);
  
  return result;
}

bool gro_next_x_v(FILE *status,real *t,int natoms,rvec x[],rvec v[],matrix box)
{
  bool result,vel;
  int i,d;
    
  do {
    result = gro_next_x_or_v(status,t,natoms,x,v,box);
    vel=FALSE;
    for (i=0; (i<natoms); i++)
      for (d=0; (d<DIM); d++)
	vel=vel || v[d][XX];
  } while(result && !vel);
  
  return result;
}

int gro_first_x_v(FILE *status, real *t, rvec **x, rvec **v, matrix box)
{
  bool result,vel;
  int  i,d,natoms;
  
  natoms = gro_first_x_or_v(status,t,x,v,box);
  result = natoms;
  vel=FALSE;
  for (i=0; (i<natoms); i++)
    for (d=0; (d<DIM); d++)
      vel=vel || v[d][XX];
  while(result && !vel) {
    result = gro_next_x_or_v(status,t,natoms,*x,*v,box);
    vel=FALSE;
    for (i=0; (i<natoms); i++)
      for (d=0; (d<DIM); d++)
	vel=vel || v[d][XX];
  }
  return result;
}

void write_hconf_indexed_p(FILE *out,char *title,t_atoms *atoms,
			   int nx,atom_id index[], int pr,
			   rvec *x,rvec *v,matrix box)
{
  char resnm[6],nm[6],format[100];
  int  ai,i,resnr,l;

  fprintf (out,"%s\n",title);
  fprintf (out,"%5d\n",nx);
  /* build format string for printing, 
     something like "%8.3f%8.3f%8.3f\n" */
  if (pr<0)
    pr=0;
  if (pr>30)
    pr=30;
  l=pr+5;
  if (v)
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
	    l,pr,l,pr,l,pr,l,pr,l,pr,l,pr);
  else
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df\n",l,pr,l,pr,l,pr);
  
  for (i=0; (i<nx); i++) {
    ai=index[i];
    
    resnr=atoms->atom[ai].resnr;
    strcpy(resnm," ??? ");
    if (resnr<atoms->nres)
      strcpy(resnm,*atoms->resname[resnr]);
    
    if (atoms->atom)
      strcpy(nm,*atoms->atomname[ai]);
    else
      strcpy(nm," ??? ");

    fprintf(out,"%5d%-5.5s%5.5s%5d",resnr+1,resnm,nm,ai+1);
    /* next fprintf uses built format string */
    if (v)
      fprintf(out,format,
	      x[ai][XX], x[ai][YY], x[ai][ZZ], v[ai][XX],v[ai][YY],v[ai][ZZ]);
    else
      fprintf(out,format,
	      x[ai][XX], x[ai][YY], x[ai][ZZ]);
  }

  if (pr<5) 
    pr=5;
  l=pr+5;
  
  if (box[XX][YY] || box[XX][ZZ] || box[YY][XX] || box[YY][ZZ] ||
      box[ZZ][XX] || box[ZZ][YY]) {
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df"
	    "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
	    l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr,l,pr);
    fprintf(out,format,
	    box[XX][XX],box[YY][YY],box[ZZ][ZZ],
	    box[XX][YY],box[XX][ZZ],box[YY][XX],
	    box[YY][ZZ],box[ZZ][XX],box[ZZ][YY]);
  } else {
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df\n",l,pr,l,pr,l,pr);
    fprintf(out,format,
	    box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
  }
  fflush(out);
}
		 
void write_hconf_indexed(FILE *out,char *title,t_atoms *atoms,
			 int nx,atom_id index[],
			 rvec *x,rvec *v,matrix box)
{
  write_hconf_indexed_p(out,title,atoms,nx,index,3,x,v,box);
}

void write_hconf_p(FILE *out,char *title,t_atoms *atoms, int pr,
		   rvec *x,rvec *v,matrix box)
{
  atom_id *aa;
  int     i;
  
  snew(aa,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    aa[i]=i;
  write_hconf_indexed_p(out,title,atoms,atoms->nr,aa,pr,x,v,box);
  sfree(aa);
}

void write_hconf(FILE *out,char *title,t_atoms *atoms,
		 rvec *x,rvec *v,matrix box)
{
  write_hconf_p(out, title, atoms, 3, x, v, box);
}

void write_conf_p(char *outfile, char *title, t_atoms *atoms, int pr,
		  rvec *x, rvec *v,matrix box)
{
  FILE *out;

  out=ffopen(outfile,"w");
  write_hconf_p(out,title,atoms,pr,x,v,box);

  ffclose (out);
}

void write_conf(char *outfile, char *title, t_atoms *atoms,
		rvec *x, rvec *v,matrix box)
{
  write_conf_p(outfile, title, atoms, 3, x, v, box);
}

void read_pdb_conf(char *infile,t_atoms *atoms,rvec x[],matrix box)
{
  FILE      *in;
  t_symtab  tab;
  t_pdbatom *pdba=NULL;
  rvec      *xx=NULL;
  int       i,natom;
  
  in    = ffopen(infile,"r");
  natom = read_pdbatoms(in,&pdba,box,FALSE);
  ffclose(in);

  open_symtab(&tab);  
  pdb2atoms(natom,pdba,atoms,&xx,&tab);
  close_symtab(&tab);
  for(i=0; (i<natom); i++)
    copy_rvec(xx[i],x[i]);
  sfree(xx);
  sfree(pdba);
}

void change_name(char *name)
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
    strcpy(name,"OXT");;
  }
}
    
void write_pdb_conf(char *outfile,t_atoms *atoms,rvec x[],matrix box,
		    bool bChange)
{
  char resnm[6],nm[6];
  int  i,resnr;
  FILE *out;

  out=ffopen(outfile,"w");

  fprintf(out,"HEADER    %s\n",bromacs());
  if (box != NULL) {
    fprintf(out,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f %6.2f%6.2f%6.2f P 1            1\n",
	    10*box[XX][XX],10*box[YY][YY],10*box[ZZ][ZZ],90.0,90.0,90.0);
  }
  for (i=0; (i<atoms->nr); i++) {
    resnr=atoms->atom[i].resnr;
    strcpy(resnm,*atoms->resname[resnr]);
    strcpy(nm,*atoms->atomname[i]);
    if (bChange)
      change_name(nm);
    if (strlen(nm)==4)
      fprintf(out,"ATOM  %5d %-4.4s %3.3s  %4d    ",i+1,nm,resnm,resnr+1);
    else
      fprintf(out,"ATOM  %5d  %-4.4s%3.3s  %4d    ",i+1,nm,resnm,resnr+1);
    fprintf(out,"%8.3f%8.3f%8.3f  1.00  0.00\n",10*x[i][XX],10*x[i][YY],10*x[i][ZZ]);
  }
  fprintf(out,"TER\n");
  ffclose(out);
}

void write_pdb_confs(char *outfile,t_atoms **atoms,rvec *x[],int number)
{
  char resnm[6],nm[6];
  int  i,resnr;
  FILE *out;
  int n;
  char chain;
  out=ffopen(outfile,"w");
  
  for(n=0;(n<number);n++) {
    chain=n+65;
    fprintf(stderr,"writing chain %c\n",chain);
    for (i=0; (i<atoms[n]->nr); i++) {
      resnr=atoms[n]->atom[i].resnr;
      strcpy(resnm,*atoms[n]->resname[resnr]);
      strcpy(nm,*atoms[n]->atomname[i]);
      change_name(nm);
      if (strlen(nm)==4)
	fprintf(out,"ATOM  %5d %-4.4s %3.3s %c%4d    ",i+1,nm,resnm,chain,resnr+1);
      else
	fprintf(out,"ATOM  %5d  %-4.4s%3.3s %c%4d    ",i+1,nm,resnm,chain,resnr+1);
      fprintf(out,"%8.3f%8.3f%8.3f  1.00  0.00\n",10*x[n][i][XX],10*x[n][i][YY],10*x[n][i][ZZ]);
    }
    fprintf(out,"TER\n");
  }
  ffclose(out);
}

void hwrite_pdb_conf_indexed(FILE *out,t_atoms *atoms,rvec x[],matrix box,
			     int gnx,atom_id index[])
{
  char resnm[6],nm[6];
  int  ii,i,resnr;

  fprintf(out,"HEADER    %s\n",bromacs());
  if (box != NULL) {
    fprintf(out,"REMARK    THIS IS A SIMULATION BOX\n");
    fprintf(out,"CRYST1%9.3f%9.3f%9.3f %6.2f%6.2f%6.2f P 1            1\n",
	    10*box[XX][XX],10*box[YY][YY],10*box[ZZ][ZZ],90.0,90.0,90.0);
  }
  for (ii=0; (ii<gnx); ii++) {
    i=index[ii];
    resnr=atoms->atom[i].resnr;
    strcpy(resnm,*atoms->resname[resnr]);
    strcpy(nm,*atoms->atomname[i]);
    change_name(nm);
    if (strlen(nm)==4)
      fprintf(out,"ATOM  %5d %-4.4s %3.3s  %4d    ",i+1,nm,resnm,resnr+1);
    else
      fprintf(out,"ATOM  %5d  %-4.4s%3.3s  %4d    ",i+1,nm,resnm,resnr+1);
    fprintf(out,"%8.3f%8.3f%8.3f  1.00  0.00\n",10*x[i][XX],10*x[i][YY],10*x[i][ZZ]);
  }
  fprintf(out,"TER\n");
}

void write_pdb_conf_indexed(char *outfile,t_atoms *atoms,rvec x[],matrix box,
			    int gnx,atom_id index[])
{
  FILE *out;

  out=ffopen(outfile,"w");
  hwrite_pdb_conf_indexed(out,atoms,x,box,gnx,index);
  fflush(out);
  ffclose(out);
}


void write_xdr_conf(char *outfile,char *title,t_atoms *atoms,
		    rvec x[],rvec v[],matrix box)
{
  XDR xd;
  int i,j,t;
  real prec;

  int num_of_coord;
  
  xdropen(&xd,outfile,"w");
  
  xdr_string(&xd,&title,strlen(title));
  xdr_int(&xd,&atoms->nr);
  
  for(i=0;(i<atoms->nr);i++) {
    xdr_string(&xd,&(*atoms->atomname[i]),6);
    xdr_int(&xd,&(atoms->atom[i].resnr));
    xdr_string(&xd,&(*atoms->resname[atoms->atom[i].resnr]),6);
  }
 
  /* write the coordinates */
  prec=1000.0;
  xdr3drcoord(&xd, x[0], &atoms->nr, &prec);

  /* write the velocities */
  prec=10000.0;
  xdr3drcoord(&xd, v[0], &atoms->nr, &prec);

  /* write the box */
  num_of_coord = DIM;
  prec=1000.0;
  xdr3drcoord(&xd, box[0], &num_of_coord, &prec);

  xdrclose(&xd);
}

void read_xdr_coordnum(char *infile,int *natoms)
{
  XDR xd;
  char *title;
  
  /* */
  snew(title,STRLEN);

  /* read the xdr file */
  xdropen(&xd, infile,"r");
  xdr_string(&xd, &title, STRLEN);
  xdr_int(&xd,natoms);
  xdrclose(&xd);
}

void read_xdr_conf(char *infile,char *title,t_atoms *atoms,rvec x[],rvec v[],matrix box)
{
  XDR xd;
  int n;
  real prec;
  int num_of_coord;
  char *line;
  
  static t_symtab symtab;
  char name[6];
  
  /* */
  snew(line,STRLEN);
  open_symtab(&symtab);
  
  /* read the xdr file */
  xdropen(&xd, infile,"r");
  *line = '\0';
  xdr_string(&xd, &title, STRLEN);
  xdr_int(&xd, &atoms->nr);
  
  /* read the titles and strings */
  atoms->nres=0;
  for(n=0;(n<atoms->nr);n++) {
    /* atomname */
    xdr_string(&xd,&line,STRLEN);
    memcpy(name,line,5);
    name[5]='\0';
    atoms->atomname[n]=put_symtab(&symtab,name);
    
    /* residue number */
    xdr_int(&xd,&(atoms->atom[n].resnr));
    if ( atoms->atom[n].resnr > atoms->nres) 
      atoms->nres=atoms->atom[n].resnr;

    /* residue name */
    xdr_string(&xd,&line,STRLEN);
    memcpy(name,line,5);
    name[5]='\0';
    atoms->resname[atoms->atom[n].resnr]=put_symtab(&symtab,name);

  }
  atoms->nres++;

  /* read coordinates */
  xdr3drcoord(&xd, x[0], &atoms->nr, &prec);


  /* read velocities */
  xdr3drcoord(&xd, v[0], &atoms->nr, &prec);

  /* read the box */
  num_of_coord = DIM;
  xdr3drcoord(&xd, box[0], &num_of_coord, &prec);
  
  xdrclose(&xd);
  close_symtab(&symtab);
}

void get_stx_coordnum (char *infile,int *natoms)
{
  int ftp;
  int i1,i2;
  real r1,r2;

  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    get_coordnum (infile, natoms);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    get_pdb_coordnum (infile, natoms);
    break;
  case efTPB:
    read_status(infile,&i1,&r1,&r2,NULL,
		NULL,NULL,NULL,
		natoms,NULL,NULL,NULL,&i2,
		NULL,NULL);
    break;
  default:
    fatal_error(0,"Not supported in read_stx_conf: %s",infile);
  }
}

void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
		   rvec x[],rvec v[], matrix box)
{
  int ftp;
  t_topology *top;
  int natoms,i1,i2;
  real r1,r2;

  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    read_whole_conf(infile, title, atoms, x, v, box);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    read_pdb_conf(infile, atoms, x, box);
    break;
  case efTPB:
    snew(top,1); 
    read_status(infile,&i1,&r1,&r2,NULL,
		box,NULL,NULL,
		&natoms,x,v,NULL,&i2,
		NULL,top);
    strcpy(title,*(top->name));
    *atoms=top->atoms;
    sfree(top);
    break;
  default:
    fatal_error(0,"Not supported in read_stx_conf: %s",infile);
  }
}
