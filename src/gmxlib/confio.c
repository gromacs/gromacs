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
static char *SRCID_confio_c = "$Id$";

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
#include "filenm.h"
#include "pdbio.h"
#include "tpxio.h"
#include "fatal.h"
#include "copyrite.h"
#include "filenm.h"
/*#include "statusio.h"*/

void clear_g96info(t_g96info *info)
{
  info->bTitle = FALSE;
  info->bTime  = FALSE;
  info->bAtoms = FALSE;
  info->bPos   = FALSE;
  info->bVel   = FALSE;
  info->bBox   = FALSE;
}

#define CHAR_SHIFT 24

static int read_g96_pos(char line[],t_symtab *symtab,FILE *fp,char *infile,
			int nwanted,t_g96info *info,
			t_atoms *atoms,rvec *x)
{
  bool   bEnd;
  int    natoms,atnr,resnr,oldres,newres,shift;
  char   anm[STRLEN],resnm[STRLEN];
  char   c1,c2;
  double db1,db2,db3;

  natoms = 0;

  if (info->bPos) {
    if (info->bAtoms)
      shift = CHAR_SHIFT;
    else
      shift = 0;
    newres  = 0;
    oldres  = -666; /* Unlikely number for the first residue! */
    bEnd    = FALSE;
    while (!bEnd && fgets2(line,STRLEN,fp)) {
      bEnd = (strncmp(line,"END",3) == 0);
      if (!bEnd  && (line[0] != '#')) {
	if (sscanf(line+shift,"%15lf%15lf%15lf",&db1,&db2,&db3) != 3)
	  fatal_error(0,"Did not find 3 coordinates for atom %d in %s\n",
		      natoms+1,infile);
	if ((nwanted != -1) && (natoms >= nwanted))
	  fatal_error(0,
		      "Found more coordinates (%d) in %s than expected %d\n",
		      natoms,infile,nwanted);
	if (atoms) {
	  if (atoms && info->bAtoms &&
	      (sscanf(line,"%5d%c%5s%c%5s%7d",&resnr,&c1,resnm,&c2,anm,&atnr) 
	       != 6)) {
	    if (oldres>=0)
	      resnr = oldres;
	    else {
	      resnr    = 1;
	      strcpy(resnm,"???"); 
	    }
	    strcpy(anm,"???"); 
	  }
	  atoms->atomname[natoms]=put_symtab(symtab,anm);
	  if (resnr != oldres) {
	    oldres = resnr;
	    if (newres >= atoms->nr)
	      fatal_error(0,"More residues than atoms in %s (natoms = %d)",
			  infile,atoms->nr);
	    atoms->resname[newres] = put_symtab(symtab,resnm);
	    newres++;
	    if (newres > atoms->nres)
	      atoms->nres = newres;
	  }
	  resnr = newres;
	  atoms->atom[natoms].resnr = resnr-1;
	}
	if (x) {
	  x[natoms][0] = db1;
	  x[natoms][1] = db2;
	  x[natoms][2] = db3;
	}
	natoms++;
      }
    }
    if ((nwanted != -1) && natoms != nwanted)
      fprintf(stderr,
	      "Warning: found less coordinates (%d) in %s than expected %d\n",
	      natoms,infile,nwanted);
  }

  return natoms;
}

static int read_g96_vel(char line[],FILE *fp,char *infile,
			int nwanted,t_g96info *info,rvec *v)
{
  bool   bEnd;
  int    natoms=-1,shift;
  double db1,db2,db3;

  if (v && info->bVel) {
    if (strcmp(line,"VELOCITYRED") == 0)
      shift = 0;
    else
      shift = CHAR_SHIFT;
    natoms = 0;
    bEnd = FALSE;
    while (!bEnd && fgets2(line,STRLEN,fp)) {
      bEnd = (strncmp(line,"END",3) == 0);
      if (!bEnd && (line[0] != '#')) {
	if (sscanf(line+shift,"%15lf%15lf%15lf",&db1,&db2,&db3) != 3)
	  fatal_error(0,"Did not find 3 velocities for atom %d in %s",
			natoms+1,infile);
	if ((nwanted != -1) && (natoms >= nwanted))
	  fatal_error(0,"Found more velocities (%d) in %s than expected %d\n",
		      natoms,infile,nwanted);
	if (v) {
	  v[natoms][0] = db1;
	  v[natoms][1] = db2;
	  v[natoms][2] = db3;
	}
	natoms++;
      }
    }
    if ((nwanted != -1) && (natoms != nwanted))
      fprintf(stderr,
	      "Warning: found less velocities (%d) in %s than expected %d\n",
	      natoms,infile,nwanted);
  }
  
  return natoms;
}

int read_g96_conf(FILE *fp,char *infile,int nwanted,t_g96info *info,
		  char *title,t_atoms *atoms,rvec *x, rvec *v,matrix box)
{
  static t_symtab *symtab=NULL;
  static char line[STRLEN+1]; /* VERY DIRTY, you can not read two       *
		               * Gromos96 trajectories at the same time */  
  bool   bAtStart,bTime,bAtoms,bPos,bVel,bBox,bEnd,bFinished;
  int    natoms;
  double db1,db2,db3;

  bAtStart = (ftell(fp) == 0);

  clear_g96info(info);
  
  if (!symtab) {
    snew(symtab,1);
    open_symtab(symtab);
  }
  
  natoms=0;

  if (bAtStart) {
    while ( !info->bTitle && fgets2(line,STRLEN,fp))
      info->bTitle = (strcmp(line,"TITLE") == 0);
    if (title)
      fgets2(title,STRLEN,fp);
    bEnd = FALSE;
    while (!bEnd && fgets2(line,STRLEN,fp))
      bEnd = (strcmp(line,"END") == 0);
    fgets2(line,STRLEN,fp);
  }
  
  /* Do not get a line if we are not at the start of the file, *
   * because without a parameter file we don't know what is in *
   * the trajectory and we have already read the line in the   *
   * previous call (VERY DIRTY).                               */
  bFinished = FALSE;
  do {
    bTime  = (strcmp(line,"TIMESTEP") == 0);
    bAtoms = (strcmp(line,"POSITION") == 0);
    bPos   = (bAtoms || (strcmp(line,"POSITIONRED") == 0));
    bVel   = (strncmp(line,"VELOCITY",8) == 0);
    bBox   = (strcmp(line,"BOX") == 0);
    if (bTime) {
      if (!info->bTime && !info->bPos) {
	info->bTime = bTime;
	do 
	  bFinished = (fgets2(line,STRLEN,fp) == NULL);
	while (!bFinished && (line[0] == '#'));
	sscanf(line,"%15d%15lf",&(info->step),&db1);
	info->time = db1;
      } else
	bFinished = TRUE;
    }
    if (bPos) {
      if (!info->bPos) {
	info->bAtoms = bAtoms;
	info->bPos   = bPos;
	natoms = read_g96_pos(line,symtab,fp,infile,nwanted,info,atoms,x);
      } else
	bFinished = TRUE;
    }
    if (v && bVel) {
      info->bVel = bVel;
      natoms = read_g96_vel(line,fp,infile,nwanted,info,v);
    }
    if (bBox) {
      info->bBox = bBox;
      clear_mat(box);
      bEnd = FALSE;
      while (!bEnd && fgets2(line,STRLEN,fp)) {
	bEnd = (strncmp(line,"END",3) == 0);
	if (!bEnd && (line[0] != '#')) {
	  if (sscanf(line,"%15lf%15lf%15lf",&db1,&db2,&db3) != 3)
	    fatal_error(0,"Found a BOX line, but no box in %s",infile);
	  box[XX][XX] = db1;
	  box[YY][YY] = db2;
	  box[ZZ][ZZ] = db3;
	}
      }
      bFinished = TRUE;
    }
  } while (!bFinished && fgets2(line,STRLEN,fp));
  
  close_symtab(symtab);
  
  return natoms;
}

void write_g96_conf(FILE *out,char *title,t_atoms *atoms,
		    rvec *x,rvec *v,matrix box,
		    int nindex,atom_id *index)
{
  int nout,i,a;
  
  if (index)
    nout = nindex;
  else
    nout = atoms->nr; 

  fprintf(out,"TITLE\n%s\nEND\n",title);
  if (x) {
    fprintf(out,"POSITION\n");
    for(i=0; i<nout; i++) {
      if (index)
	a = index[i];
      else
	a = i;
      fprintf(out,"%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
	      atoms->atom[a].resnr+1,*atoms->resname[atoms->atom[a].resnr],
	      *atoms->atomname[a],i+1,x[a][XX],x[a][YY],x[a][ZZ]);
    }
    fprintf(out,"END\n");
  }
  if (v) {
    fprintf(out,"VELOCITY\n");
    for(i=0; i<nout; i++) {
      if (index)
	a = index[i];
      else
	a = i;
      fprintf(out,"%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
	      atoms->atom[a].resnr+1,*atoms->resname[atoms->atom[a].resnr],
	      *atoms->atomname[a],i+1,v[a][XX],v[a][YY],v[a][ZZ]);
    }
    fprintf(out,"END\n");
  }
  fprintf(out,"BOX\n");
  fprintf(out,"%15.9f%15.9f%15.9f\n",box[XX][XX],box[YY][YY],box[ZZ][ZZ]);
  fprintf(out,"END\n");
}

static void get_coordnum_fp (FILE *in,char *title, int *natoms)
{
  char line[STRLEN+1];

  fgets2 (title,STRLEN,in);
  fgets2 (line,STRLEN,in);
  sscanf (line,"%d",natoms);
}

static void get_coordnum (char *infile,int *natoms)
{
  FILE *in;
  char title[STRLEN];
  
  in=ffopen(infile,"r");
  get_coordnum_fp(in,title,natoms);
  ffclose (in);
}

static void get_w_conf(FILE *in, char *infile, char *title,
		       t_atoms *atoms, rvec x[],rvec *v, matrix box)
{
  static t_symtab *symtab=NULL;
  char   name[6];
  char   line[STRLEN+1];
  char   buf[256];
  char   format[30];
  double x1,y1,z1,x2,y2,z2;
  rvec   xmin,xmax;
  int    natoms,i,m,resnr,newres,oldres,prec;
  bool   bFirst;
  char   *p1,*p2;
  
  newres  = 0;
  oldres  = -666; /* Unlikely number for the first residue! */
  prec    = 0;
  
  if (!symtab) {
    snew(symtab,1);
    open_symtab(symtab);
  }
  fgets2(title,STRLEN,in);

  /* read the number of atoms */
  fgets2(line,STRLEN,in);
  sscanf(line,"%d",&natoms);
  if (natoms > atoms->nr)
    fatal_error(0,"gro file contains more atoms (%d) than expected (%d)",
		natoms,atoms->nr);
  else if (natoms <  atoms->nr)
    fprintf(stderr,"Warning: gro file contains less atoms (%d) than expected"
	    " (%d)\n",natoms,atoms->nr);
  
  bFirst=TRUE;
  
  /* just pray the arrays are big enough */
  for (i=0; (i < natoms) ; i++) {
    if ((fgets2 (line,STRLEN,in)) == NULL) {
      unexpected_eof(infile,i+2);
    }
    if (strlen(line) < 39)
      fatal_error(0,"Invalid line in %s for atom %d:\n%s",infile,i+1,line);

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
      /* this will be something like "%8lf%8lf%8lf" */
    }
    
    /* residue number*/
    memcpy(name,line,5);
    name[5]='\0';
    sscanf(name,"%d",&resnr);
    memcpy(name,line+5,5);
    name[5]='\0';
    if (resnr != oldres) {
      oldres = resnr;
      if (newres >= natoms)
	fatal_error(0,"More residues than atoms in %s (natoms = %d)",
		    infile,natoms);
      atoms->resname[newres] = put_symtab(symtab,name);
      newres++;
    }
    resnr = newres;
    atoms->atom[i].resnr = resnr-1;

    /* atomname */
    memcpy(name,line+10,5);
    atoms->atomname[i]=put_symtab(symtab,name);
   
    /* eventueel controle atomnumber met i+1 */

    /* coordinates (start after residue shit) */
    /* 'format' was built previously */
    if (sscanf (line+20,format,&x1,&y1,&z1) != 3) {
      too_few();
    }
    else {
      x[i][XX]=x1;
      x[i][YY]=y1;
      x[i][ZZ]=z1;
    }

    /* velocities (start after residues and coordinates) */
    /* 'format' was built previously */
    if (v) {
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
  }
  atoms->nres=newres;

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
  close_symtab(symtab);
}

static void read_whole_conf(char *infile, char *title,
			    t_atoms *atoms, rvec x[],rvec *v, matrix box)
{
  FILE   *in;
  
  /* open file */
  in=ffopen(infile,"r");

  get_w_conf(in, infile, title, atoms, x, v, box);
  
  fclose(in);
}

static void get_conf(FILE *in, char *title, int *natoms, 
		     rvec x[],rvec *v,matrix box)
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
		     rvec x[],rvec *v,matrix box)
{
  t_atoms  atoms;
  char   title[STRLEN],*p;
  double tt;

  if (eof(status))
    return FALSE;

  atoms.nr=natoms;
  snew(atoms.atom,natoms);
  atoms.nres=natoms;
  snew(atoms.resname,natoms);
  snew(atoms.atomname,natoms);
  
  get_w_conf(status,title,title,&atoms,x,v,box);
  
  sfree(atoms.atom);
  sfree(atoms.resname);
  sfree(atoms.atomname);

  if ((p=strstr(title,"t=")) != NULL) {
    p+=2;
    if (sscanf(p,"%lf",&tt)==1)
      *t=tt;
    else
      *t=0.0;
  }
  
  if (atoms.nr != natoms)
    fatal_error(0,"Number of atoms in gro frame (%d) doesn't match the number in the previous frame (%d)",atoms.nr,natoms);
  
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

bool gro_next_v(FILE *status,real *t,int natoms,rvec *v,matrix box)
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

bool gro_next_x_v(FILE *status,real *t,int natoms,rvec x[],rvec *v,matrix box)
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
  int  ai,i,resnr,l,vpr;

  fprintf (out,"%s\n",title[0]?title:bromacs());
  fprintf (out,"%5d\n",nx);
  /* build format string for printing, 
     something like "%8.3f" for x and "%8.4f" for v */
  if (pr<0)
    pr=0;
  if (pr>30)
    pr=30;
  l=pr+5;
  vpr=pr+1;
  if (v)
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
	    l,pr,l,pr,l,pr,l,vpr,l,vpr,l,vpr);
  else
    sprintf(format,"%%%d.%df%%%d.%df%%%d.%df\n",l,pr,l,pr,l,pr);
  
  for (i=0; (i<nx); i++) {
    ai=index[i];
    
    resnr=atoms->atom[ai].resnr;
    strcpy(resnm," ??? ");
    if (resnr < atoms->nres)
      strcpy(resnm,*atoms->resname[resnr]);
    
    if (atoms->atom)
      strcpy(nm,*atoms->atomname[ai]);
    else
      strcpy(nm," ??? ");

    fprintf(out,"%5d%-5.5s%5.5s%5d",resnr+1,resnm,nm,(ai+1)%100000);
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

static void write_conf(char *outfile, char *title, t_atoms *atoms,
		       rvec *x, rvec *v,matrix box)
{
  write_conf_p(outfile, title, atoms, 3, x, v, box);
}

void write_xdr_conf(char *outfile,char *title,t_atoms *atoms,
		    rvec x[],rvec *v,matrix box)
{
  XDR xd;
  int i;
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

void read_xdr_conf(char *infile,char *title,t_atoms *atoms,rvec x[],rvec *v,matrix box)
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

void write_sto_conf_indexed(char *outfile,char *title,t_atoms *atoms, 
			    rvec x[],rvec *v,matrix box,
			    atom_id nindex,atom_id index[])
{
  FILE       *out;
  int        ftp;

  ftp=fn2ftp(outfile);
  switch (ftp) {
  case efGRO:
    out=ffopen(outfile,"w");
    write_hconf_indexed(out, title, atoms, nindex, index, x, v, box);
    fclose(out);
    break;
  case efG96:
    out=ffopen(outfile,"w");
    write_g96_conf(out, title, atoms, x, v, box, nindex, index);
    fclose(out);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    out=ffopen(outfile,"w");
    hwrite_pdb_conf_indexed(out, title, atoms, x, box, nindex, index);
    fclose(out);
    break;
  case efTPR:
  case efTPB:
  case efTPA:
    fatal_error(0,"Sorry, can not write a topology to %s",outfile);
    break;
  default:
    fatal_error(0,"Not supported in write_sto_conf_indexed: %s",outfile);
  }
}

void write_sto_conf(char *outfile, char *title,t_atoms *atoms, 
		   rvec x[],rvec *v, matrix box)
{
  FILE       *out;
  int        ftp;

  ftp=fn2ftp(outfile);
  switch (ftp) {
  case efGRO:
    write_conf(outfile, title, atoms, x, v, box);
    break;
  case efG96:
    out=ffopen(outfile,"w");
    write_g96_conf(out, title, atoms, x, v, box, -1, NULL);
    fclose(out);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    out=ffopen(outfile,"w");
    write_pdbfile(out, title, atoms, x, box, 0, TRUE);
    fclose(out);
    break;
  case efTPR:
  case efTPB:
  case efTPA:
    fatal_error(0,"Sorry, can not write a topology to %s",outfile);
    break;
  default:
    fatal_error(0,"Not supported in write_sto_conf: %s",outfile);
  }
}

void get_stx_coordnum(char *infile,int *natoms)
{
  FILE *in;
  int ftp;
  t_g96info g96info;
  matrix dumbox;

  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    get_coordnum(infile, natoms);
    break;
  case efG96:
    in=ffopen(infile,"r");
    *natoms=read_g96_conf(in,infile,-1,&g96info,NULL,NULL,NULL,NULL,dumbox);
    fclose(in);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    in=ffopen(infile,"r");
    get_pdb_coordnum(in, natoms);
    fclose(in);
    break;
  case efTPA:
  case efTPB:
  case efTPR: {
    t_tpxheader tpx;
    
    read_tpxheader(infile,&tpx);
    *natoms = tpx.natoms;
    break;
  }
  default:
    fatal_error(0,"Not supported in get_stx_coordnum: %s",infile);
  }
}

void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
		   rvec x[],rvec *v, matrix box)
{
  FILE *in;
  t_topology *top;
  t_g96info  g96info;
  int        ftp,natoms,i1;
  real       r1,r2;

  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    read_whole_conf(infile, title, atoms, x, v, box);
    break;
  case efG96:
    clear_mat(box);
    in = ffopen(infile,"r");
    read_g96_conf(in, infile, atoms->nr, &g96info, title, atoms, x, v, box);
    fclose(in);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    read_pdb_conf(infile, title, atoms, x, box, TRUE);
    break;
  case efTPR:
  case efTPB:
  case efTPA: 
    snew(top,1);
    read_tpx(infile,&i1,&r1,&r2,NULL,box,&natoms,x,v,NULL,top);
    strcpy(title,*(top->name));
    *atoms=top->atoms;
    break;
  default:
    fatal_error(0,"Not supported in read_stx_conf: %s",infile);
  }
}


