/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.2
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2007, The GROMACS development team,
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
#include "futil.h"
#include "xdrf.h"
#include "filenm.h"
#include "pdbio.h"
#include "tpxio.h"
#include "gmx_fatal.h"
#include "copyrite.h"
#include "filenm.h"
#include "statutil.h"
#include "pbc.h"

#define CHAR_SHIFT 24

static int read_g96_pos(char line[],t_symtab *symtab,FILE *fp,char *infile,
			t_trxframe *fr)
{
  t_atoms *atoms;
  bool   bEnd;
  int    nwanted,natoms,atnr,resnr,oldres,newres,shift;
  char   anm[STRLEN],resnm[STRLEN];
  char   c1,c2;
  double db1,db2,db3;
  
  nwanted = fr->natoms;

  atoms = fr->atoms;

  natoms = 0;

  if (fr->bX) {
    if (fr->bAtoms)
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
	  gmx_fatal(FARGS,"Did not find 3 coordinates for atom %d in %s\n",
		      natoms+1,infile);
	if ((nwanted != -1) && (natoms >= nwanted))
	  gmx_fatal(FARGS,
		      "Found more coordinates (%d) in %s than expected %d\n",
		      natoms,infile,nwanted);
	if (atoms) {
	  if (atoms && fr->bAtoms &&
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
	      gmx_fatal(FARGS,"More residues than atoms in %s (natoms = %d)",
			  infile,atoms->nr);
	    atoms->resname[newres] = put_symtab(symtab,resnm);
	    newres++;
	    if (newres > atoms->nres)
	      atoms->nres = newres;
	  }
	  resnr = newres;
	  atoms->atom[natoms].resnr = resnr-1;
	}
	if (fr->x) {
	  fr->x[natoms][0] = db1;
	  fr->x[natoms][1] = db2;
	  fr->x[natoms][2] = db3;
	}
	natoms++;
      }
    }
    if ((nwanted != -1) && natoms != nwanted)
      fprintf(stderr,
	      "Warning: found less coordinates (%d) in %s than expected %d\n",
	      natoms,infile,nwanted);
  }

  fr->natoms = natoms;

  return natoms;
}

static int read_g96_vel(char line[],FILE *fp,char *infile,
			t_trxframe *fr)
{
  bool   bEnd;
  int    nwanted,natoms=-1,shift;
  double db1,db2,db3;

  nwanted = fr->natoms;

  if (fr->v && fr->bV) {
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
	  gmx_fatal(FARGS,"Did not find 3 velocities for atom %d in %s",
		      natoms+1,infile);
	if ((nwanted != -1) && (natoms >= nwanted))
	  gmx_fatal(FARGS,"Found more velocities (%d) in %s than expected %d\n",
		      natoms,infile,nwanted);
	if (fr->v) {
	  fr->v[natoms][0] = db1;
	  fr->v[natoms][1] = db2;
	  fr->v[natoms][2] = db3;
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

int read_g96_conf(FILE *fp,char *infile,t_trxframe *fr)
{
  static t_symtab *symtab=NULL;
  static char line[STRLEN+1]; /* VERY DIRTY, you can not read two       *
		               * Gromos96 trajectories at the same time */  
  bool   bAtStart,bTime,bAtoms,bPos,bVel,bBox,bEnd,bFinished;
  int    natoms,nbp;
  double db1,db2,db3,db4,db5,db6,db7,db8,db9;

  bAtStart = (ftell(fp) == 0);

  clear_trxframe(fr,FALSE);
  
  if (!symtab) {
    snew(symtab,1);
    open_symtab(symtab);
  }
  
  natoms=0;

  if (bAtStart) {
    while ( !fr->bTitle && fgets2(line,STRLEN,fp))
      fr->bTitle = (strcmp(line,"TITLE") == 0);
    if (fr->title)
      fgets2(fr->title,STRLEN,fp);
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
      if (!fr->bTime && !fr->bX) {
	fr->bStep = bTime;
	fr->bTime = bTime;
	do 
	  bFinished = (fgets2(line,STRLEN,fp) == NULL);
	while (!bFinished && (line[0] == '#'));
	sscanf(line,"%15d%15lf",&(fr->step),&db1);
	fr->time = db1;
      } else
	bFinished = TRUE;
    }
    if (bPos) {
      if (!fr->bX) {
	fr->bAtoms = bAtoms;
	fr->bX     = bPos;
	natoms = read_g96_pos(line,symtab,fp,infile,fr);
      } else
	bFinished = TRUE;
    }
    if (fr->v && bVel) {
      fr->bV = bVel;
      natoms = read_g96_vel(line,fp,infile,fr);
    }
    if (bBox) {
      fr->bBox = bBox;
      clear_mat(fr->box);
      bEnd = FALSE;
      while (!bEnd && fgets2(line,STRLEN,fp)) {
	bEnd = (strncmp(line,"END",3) == 0);
	if (!bEnd && (line[0] != '#')) {
	  nbp = sscanf(line,"%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf%15lf",
		       &db1,&db2,&db3,&db4,&db5,&db6,&db7,&db8,&db9);
	  if (nbp < 3)
	    gmx_fatal(FARGS,"Found a BOX line, but no box in %s",infile);
	  fr->box[XX][XX] = db1;
	  fr->box[YY][YY] = db2;
	  fr->box[ZZ][ZZ] = db3;
	  if (nbp == 9) {
	    fr->box[XX][YY] = db4;
	    fr->box[XX][ZZ] = db5;
	    fr->box[YY][XX] = db6;
	    fr->box[YY][ZZ] = db7;
	    fr->box[ZZ][XX] = db8;
	    fr->box[ZZ][YY] = db9; 
	  }
	}
      }
      bFinished = TRUE;
    }
  } while (!bFinished && fgets2(line,STRLEN,fp));
  
  close_symtab(symtab);

  fr->natoms = natoms;
  
  return natoms;
}

void write_g96_conf(FILE *out,t_trxframe *fr,
		    int nindex,atom_id *index)
{
  t_atoms *atoms;
  int nout,i,a;
  
  atoms = fr->atoms;

  if (index)
    nout = nindex;
  else
    nout = fr->natoms; 

  if (fr->bTitle)
    fprintf(out,"TITLE\n%s\nEND\n",fr->title);
  if (fr->bStep || fr->bTime)
    /* Officially the time format is %15.9, which is not enough for 10 ns */
    fprintf(out,"TIMESTEP\n%15d%15.6f\nEND\n",fr->step,fr->time);
  if (fr->bX) {
    if (fr->bAtoms) {
      fprintf(out,"POSITION\n");
      for(i=0; i<nout; i++) {
	if (index) a = index[i]; else a = i;
	fprintf(out,"%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
		(atoms->atom[a].resnr+1) % 100000,
		*atoms->resname[atoms->atom[a].resnr],
		*atoms->atomname[a],(i+1) % 10000000,
		fr->x[a][XX],fr->x[a][YY],fr->x[a][ZZ]);
      }
    } else {
      fprintf(out,"POSITIONRED\n");
      for(i=0; i<nout; i++) {
	if (index) a = index[i]; else a = i;
	fprintf(out,"%15.9f%15.9f%15.9f\n",
		fr->x[a][XX],fr->x[a][YY],fr->x[a][ZZ]);
      }
    }
    fprintf(out,"END\n");
  }
  if (fr->bV) {
    if (fr->bAtoms) {
      fprintf(out,"VELOCITY\n");
      for(i=0; i<nout; i++) {
	if (index) a = index[i]; else a = i;
	fprintf(out,"%5d %-5s %-5s%7d%15.9f%15.9f%15.9f\n",
		(atoms->atom[a].resnr+1) % 100000,
		*atoms->resname[atoms->atom[a].resnr],
		*atoms->atomname[a],(i+1) % 10000000,
		fr->v[a][XX],fr->v[a][YY],fr->v[a][ZZ]);
      }
    } else {
      fprintf(out,"VELOCITYRED\n");
      for(i=0; i<nout; i++) {
	if (index) a = index[i]; else a = i;
	fprintf(out,"%15.9f%15.9f%15.9f\n",
		fr->v[a][XX],fr->v[a][YY],fr->v[a][ZZ]);
      }
    }
    fprintf(out,"END\n");
  }
  if (fr->bBox) {
    fprintf(out,"BOX\n");
    fprintf(out,"%15.9f%15.9f%15.9f",
	    fr->box[XX][XX],fr->box[YY][YY],fr->box[ZZ][ZZ]);
    if (fr->box[XX][YY] || fr->box[XX][ZZ] || fr->box[YY][XX] ||
	fr->box[YY][ZZ] || fr->box[ZZ][XX] ||fr->box[ZZ][YY])
      fprintf(out,"%15.9f%15.9f%15.9f%15.9f%15.9f%15.9f",
	      fr->box[XX][YY],fr->box[XX][ZZ],fr->box[YY][XX],
	      fr->box[YY][ZZ],fr->box[ZZ][XX],fr->box[ZZ][YY]);
    fprintf(out,"\n");
    fprintf(out,"END\n");
  }
}

static int get_espresso_word(FILE *fp,char word[])
{
  int  ret,nc,i;
  
  ret = 0;
  nc = 0;
  
  do {
    i = fgetc(fp);
    if (i != EOF) {
      if (i == ' ' || i == '\n' || i == '\t') {
	if (nc > 0)
	  ret = 1;
      } else if (i == '{') {
	if (nc == 0)
	  word[nc++] = '{';
	ret = 2;
      } else if (i == '}') {
	if (nc == 0)
	  word[nc++] = '}';
	ret = 3;
      } else {
	word[nc++] = (char)i;
      }
    }
  } while (i != EOF && ret == 0);
  
  word[nc] = '\0';

  /*  printf("'%s'\n",word); */
  
  return ret;
}

static int check_open_parenthesis(FILE *fp,int r,
				  const char *infile,const char *keyword)
{
  int level_inc;
  char word[STRLEN];

  level_inc = 0;
  if (r == 2) {
    level_inc++;
  } else {
    r = get_espresso_word(fp,word);
    if (r == 2)
      level_inc++;
    else
      gmx_fatal(FARGS,"Expected '{' after '%s' in file '%s'",
		keyword,infile);
  }

  return level_inc;
}

static int check_close_parenthesis(FILE *fp,int r,
				  const char *infile,const char *keyword)
{
  int level_inc;
  char word[STRLEN];
  
  level_inc = 0;
  if (r == 3) {
    level_inc--;
  } else {
    r = get_espresso_word(fp,word);
    if (r == 3)
      level_inc--;
    else
      gmx_fatal(FARGS,"Expected '}' after section '%s' in file '%s'",
		keyword,infile);
  }

  return level_inc;
}

enum { espID, espPOS, espTYPE, espQ, espV, espF, espNR };
const char *esp_prop[espNR] = { "id", "pos", "type", "q", "v", "f" };

static void read_espresso_conf(char *infile,
			       t_atoms *atoms,rvec x[],rvec *v,matrix box)
{
  static t_symtab *symtab=NULL;
  FILE *fp;
  char word[STRLEN],buf[STRLEN];
  int  natoms,level,npar,r,nprop,p,i,m;
  int  prop[32];
  double d;
  bool bFoundParticles,bFoundProp,bFoundVariable;

  if (!symtab) {
    snew(symtab,1);
    open_symtab(symtab);
  }

  clear_mat(box);
  
  fp = ffopen(infile,"r");
  
  bFoundParticles = FALSE;
  bFoundVariable = FALSE;
  level = 0;
  while ((r=get_espresso_word(fp,word))) {
    if (level==1 && strcmp(word,"particles")==0 && !bFoundParticles) {
      bFoundParticles = TRUE;
      level += check_open_parenthesis(fp,r,infile,"particles");
      nprop = 0;
      while (level == 2 && (r=get_espresso_word(fp,word))) {
	bFoundProp = FALSE;
	for(p=0; p<espNR; p++) {
	  if (strcmp(word,esp_prop[p]) == 0) {
	    bFoundProp = TRUE;
	    prop[nprop++] = p;
	    /* printf("  prop[%d] = %s\n",nprop-1,esp_prop[prop[nprop-1]]); */
	  }
	}
	if (!bFoundProp && word[0] != '}') {
	  gmx_fatal(FARGS,"Can not read Espresso files with particle property '%s'",word);
	}
	if (r == 3)
	  level--;
      }
      
      i = 0;
      while (level > 0 && (r=get_espresso_word(fp,word))) {
	if (r == 2) {
	  level++;
	} else if (r == 3) {
	  level--;
	}
	if (level == 2) {
	  for(p=0; p<nprop; p++) {
	    switch (prop[p]) {
	    case espID:
	      r = get_espresso_word(fp,word);
	      /* Not used */
	      break;
	    case espPOS:
	      for(m=0; m<3; m++) {
		r = get_espresso_word(fp,word);
		sscanf(word,"%lf",&d);
		x[i][m] = d;
	      }
	      break;
	    case espTYPE:
	      r = get_espresso_word(fp,word);
	      atoms->atom[i].type = atoi(word);
	      break;
	    case espQ:
	      r = get_espresso_word(fp,word);
	      sscanf(word,"%lf",&d);
	      atoms->atom[i].q = d;
	      break;
	    case espV:
	      for(m=0; m<3; m++) {
		r = get_espresso_word(fp,word);
		sscanf(word,"%lf",&d);
		v[i][m] = d;
	      }
	      break;
	    case espF:
	      for(m=0; m<3; m++) {
		r = get_espresso_word(fp,word);
		/* not used */
	      }
	      break;
	    }
	  }
	  atoms->atom[i].resnr = i;
	  sprintf(buf,"T%d",atoms->atom[i].type);
	  atoms->atomname[i] = put_symtab(symtab,buf);
	  if (atoms->atom[i].type < 26) {
	    sprintf(buf,"T%c",'A'+atoms->atom[i].type);
	  } else {
	    sprintf(buf,"T%c%c",
		    'A'+atoms->atom[i].type/26,'A'+atoms->atom[i].type%26);
	  }
	  atoms->resname[atoms->atom[i].resnr] = put_symtab(symtab,buf);
	  
	  if (r == 3)
	    level--;
	  i++;
	}
      }
      atoms->nres = atoms->nr;

      if (i != atoms->nr) {
	gmx_fatal(FARGS,"Internal inconsistency in Espresso routines, read %d atoms, expected %d atoms",i,atoms->nr);
      }
    } else if (level==1 && strcmp(word,"variable")==0 && !bFoundVariable) {
      bFoundVariable = TRUE;
      level += check_open_parenthesis(fp,r,infile,"variable");
      while (level==2 && (r=get_espresso_word(fp,word))) {
	if (level==2 && strcmp(word,"box_l") == 0) {
	  for(m=0; m<3; m++) {
	    r = get_espresso_word(fp,word);
	    sscanf(word,"%lf",&d);
	    box[m][m] = d;
	  }
	  level += check_close_parenthesis(fp,r,infile,"box_l");
	}
      }
    } else if (r == 2) {
      level++;
    } else if (r == 3) {
      level--;
    }
  }
  
  if (!bFoundParticles) {
    fprintf(stderr,"Did not find a particles section in Espresso file '%s'\n",
	    infile);
  }
  
  ffclose(fp);
}

static int get_espresso_coordnum(char *infile)
{
  FILE *fp;
  char word[STRLEN];
  int  natoms,level,r;
  bool bFoundParticles;

  natoms = 0;
  
  fp = ffopen(infile,"r");
  
  bFoundParticles = FALSE;
  level = 0;
  while ((r=get_espresso_word(fp,word)) && !bFoundParticles) {
    if (level==1 && strcmp(word,"particles")==0 && !bFoundParticles) {
      bFoundParticles = TRUE;
      level += check_open_parenthesis(fp,r,infile,"particles");
      while (level > 0 && (r=get_espresso_word(fp,word))) {
	if (r == 2) {
	  level++;
	  if (level == 2)
	    natoms++;
	} else if (r == 3) {
	  level--;
	}
      }
    } else if (r == 2) {
      level++;
    } else if (r == 3) {
      level--;
    }
  }
  if (!bFoundParticles) {
    fprintf(stderr,"Did not find a particles section in Espresso file '%s'\n",
	    infile);
  }
  
  ffclose(fp);
  
  return natoms;
}

static void write_espresso_conf_indexed(FILE *out,const char *title,
					t_atoms *atoms,int nx,atom_id *index,
					rvec *x,rvec *v,matrix box)
{
  int i,j;

  fprintf(out,"# %s\n",title);
  if (TRICLINIC(box))
    fprintf(stderr,"\nWARNING: the Espresso format does not support triclinic unit-cells\n\n");
  fprintf(out,"{variable {box_l %f %f %f}}\n",box[0][0],box[1][1],box[2][2]);
  
  fprintf(out,"{particles {id pos type q%s}\n",v ? " v" : "");
  for(i=0; i<nx; i++) {
    if (index)
      j = index[i];
    else
      j = i;
    fprintf(out,"\t{%d %f %f %f %d %g",
	    j,x[j][XX],x[j][YY],x[j][ZZ],
	    atoms->atom[j].type,atoms->atom[j].q);
    if (v)
      fprintf(out," %f %f %f",v[j][XX],v[j][YY],v[j][ZZ]);
    fprintf(out,"}\n");
  }
  fprintf(out,"}\n");
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

static bool get_w_conf(FILE *in, char *infile, char *title,
		       t_atoms *atoms, int *ndec, rvec x[],rvec *v, matrix box)
{
  static t_symtab *symtab=NULL;
  char   name[6];
  char   line[STRLEN+1];
  char   buf[256];
  char   format[30];
  double x1,y1,z1,x2,y2,z2;
  rvec   xmin,xmax;
  int    natoms,i,m,resnr,newres,oldres,ddist;
  bool   bFirst,bVel;
  char   *p1,*p2;
  
  newres  = 0;
  oldres  = NOTSET; /* Unlikely number for the first residue! */
  ddist   = 0;
  
  if (!symtab) {
    snew(symtab,1);
    open_symtab(symtab);
  }
  fgets2(title,STRLEN,in);

  /* read the number of atoms */
  fgets2(line,STRLEN,in);
  sscanf(line,"%d",&natoms);
  if (natoms > atoms->nr)
    gmx_fatal(FARGS,"gro file contains more atoms (%d) than expected (%d)",
		natoms,atoms->nr);
  else if (natoms <  atoms->nr)
    fprintf(stderr,"Warning: gro file contains less atoms (%d) than expected"
	    " (%d)\n",natoms,atoms->nr);
  
  bFirst=TRUE;
  
  bVel = FALSE;

  /* just pray the arrays are big enough */
  for (i=0; (i < natoms) ; i++) {
    if ((fgets2 (line,STRLEN,in)) == NULL) {
      unexpected_eof(infile,i+2);
    }
    if (strlen(line) < 39)
      gmx_fatal(FARGS,"Invalid line in %s for atom %d:\n%s",infile,i+1,line);

    /* determine read precision from distance between periods 
       (decimal points) */
    if (bFirst) {
      bFirst=FALSE;
      p1=strchr(line,'.');
      if (p1 == NULL)
	gmx_fatal(FARGS,"A coordinate in file %s does not contain a '.'",infile);
      p2=strchr(&p1[1],'.');
      if (p1 || p2)
	ddist=p2-p1;
      else
	ddist=8;
      if (ddist<0)
	ddist=8;
      if (ddist>30)
	ddist=30;
      sprintf(format,"%%%dlf%%%dlf%%%dlf",ddist,ddist,ddist);
      /* this will be something like "%8lf%8lf%8lf" */
      *ndec = ddist-5;
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
	gmx_fatal(FARGS,"More residues than atoms in %s (natoms = %d)",
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
      if (sscanf (line+20+(3*ddist),format,&x1,&y1,&z1) != 3) {
	v[i][XX] = 0.0;
	v[i][YY] = 0.0;
	v[i][ZZ] = 0.0;
      }
      else {
	v[i][XX]=x1;
	v[i][YY]=y1;
	v[i][ZZ]=z1;
	bVel = TRUE;
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

  return bVel;
}

static void read_whole_conf(char *infile, char *title,
			    t_atoms *atoms, rvec x[],rvec *v, matrix box)
{
  FILE   *in;
  int    ndec;
  
  /* open file */
  in=ffopen(infile,"r");

  get_w_conf(in, infile, title, atoms, &ndec, x, v, box);
  
  fclose(in);
}

static void get_conf(FILE *in, char *title, int *natoms, 
		     rvec x[],rvec *v,matrix box)
{
  t_atoms  atoms;
  int      ndec;

  atoms.nr=*natoms;
  snew(atoms.atom,*natoms);
  atoms.nres=*natoms;
  snew(atoms.resname,*natoms);
  snew(atoms.atomname,*natoms);
  
  get_w_conf(in,title,title,&atoms,&ndec,x,v,box);
  
  sfree(atoms.atom);
  sfree(atoms.resname);
  sfree(atoms.atomname);
}

bool gro_next_x_or_v(FILE *status,t_trxframe *fr)
{
  t_atoms atoms;
  char    title[STRLEN],*p;
  double  tt;
  int     ndec,i;

  if (eof(status))
    return FALSE;

  atoms.nr=fr->natoms;
  snew(atoms.atom,fr->natoms);
  atoms.nres=fr->natoms;
  snew(atoms.resname,fr->natoms);
  snew(atoms.atomname,fr->natoms);
  
  fr->bV = get_w_conf(status,title,title,&atoms,&ndec,fr->x,fr->v,fr->box);
  fr->bPrec = TRUE;
  fr->prec = 1;
  /* prec = 10^ndec: */
  for(i=0; i<ndec; i++)
    fr->prec *= 10;
  fr->title = title;
  fr->bTitle = TRUE;
  fr->bX = TRUE;
  fr->bBox = TRUE;

  sfree(atoms.atom);
  sfree(atoms.resname);
  sfree(atoms.atomname);

  if ((p=strstr(title,"t=")) != NULL) {
    p+=2;
    if (sscanf(p,"%lf",&tt)==1) {
      fr->time = tt;
      fr->bTime = TRUE;
    } else {
      fr->time = 0;
      fr->bTime = FALSE;
    }
  }
  
  if (atoms.nr != fr->natoms)
    gmx_fatal(FARGS,"Number of atoms in gro frame (%d) doesn't match the number in the previous frame (%d)",atoms.nr,fr->natoms);
  
  return TRUE;
}

int gro_first_x_or_v(FILE *status,t_trxframe *fr)
{
  char title[STRLEN];
  
  frewind(status);
  fprintf(stderr,"Reading frames from gro file");
  get_coordnum_fp(status, title, &fr->natoms);
  frewind(status);
  fprintf(stderr," '%s', %d atoms.\n",title, fr->natoms);
  fr->bTitle = TRUE;
  fr->title = title;
  if (fr->natoms==0)
    gmx_file("No coordinates in gro file");
  
  snew(fr->x,fr->natoms);
  snew(fr->v,fr->natoms);
  gro_next_x_or_v(status, fr);
  
  return fr->natoms;
}

void write_hconf_indexed_p(FILE *out,char *title,t_atoms *atoms,
			   int nx,atom_id index[], int pr,
			   rvec *x,rvec *v,matrix box)
{
  char resnm[6],nm[6],format[100];
  int  ai,i,resnr,l,vpr;

  bromacs(format,99);
  fprintf (out,"%s\n",(title && title[0])?title:format);
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

    fprintf(out,"%5d%-5.5s%5.5s%5d",(resnr+1)%100000,resnm,nm,(ai+1)%100000);
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

void write_sto_conf_indexed(char *outfile,char *title,t_atoms *atoms, 
			    rvec x[],rvec *v,matrix box,
			    atom_id nindex,atom_id index[])
{
  FILE       *out;
  int        ftp;
  t_trxframe fr;

  ftp=fn2ftp(outfile);
  switch (ftp) {
  case efGRO:
    out=ffopen(outfile,"w");
    write_hconf_indexed_p(out, title, atoms, nindex, index, 3, x, v, box);
    fclose(out);
    break;
  case efG96:
    clear_trxframe(&fr,TRUE);
    fr.bTitle = TRUE;
    fr.title = title;
    fr.natoms = atoms->nr;
    fr.bAtoms = TRUE;
    fr.atoms = atoms;
    fr.bX = TRUE;
    fr.x = x;
    if (v) {
      fr.bV = TRUE;
      fr.v = v;
    }
    fr.bBox = TRUE;
    copy_mat(box,fr.box);
    out=ffopen(outfile,"w");
    write_g96_conf(out, &fr, nindex, index);
    fclose(out);
    break;
  case efPDB:
  case efBRK:
  case efENT:
  case efPQR:
    out=ffopen(outfile,"w");
    write_pdbfile_indexed(out, title, atoms, x, box, 0, -1, nindex, index);
    fclose(out);
    break;
  case efESP:
    out=ffopen(outfile,"w"); 
    write_espresso_conf_indexed(out, title, atoms, nindex, index, x, v, box);
    fclose(out);
    break;
  case efTPR:
  case efTPB:
  case efTPA:
    gmx_fatal(FARGS,"Sorry, can not write a topology to %s",outfile);
    break;
  default:
    gmx_incons("Not supported in write_sto_conf_indexed");
  }
}

void write_sto_conf(char *outfile, char *title,t_atoms *atoms, 
		   rvec x[],rvec *v, matrix box)
{
  FILE       *out;
  int        ftp;
  t_trxframe fr;

  ftp=fn2ftp(outfile);
  switch (ftp) {
  case efGRO:
    write_conf(outfile, title, atoms, x, v, box);
    break;
  case efG96:
    clear_trxframe(&fr,TRUE);
    fr.bTitle = TRUE;
    fr.title = title;
    fr.natoms = atoms->nr;
    fr.bAtoms = TRUE;
    fr.atoms = atoms;
    fr.bX = TRUE;
    fr.x = x;
    if (v) {
      fr.bV = TRUE;
      fr.v = v;
    }
    fr.bBox = TRUE;
    copy_mat(box,fr.box);
    out=ffopen(outfile,"w");
    write_g96_conf(out, &fr, -1, NULL);
    fclose(out);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    out=ffopen(outfile,"w");
    write_pdbfile(out, title, atoms, x, box, 0, -1);
    fclose(out);
    break;
  case efESP:
    out=ffopen(outfile,"w");
    write_espresso_conf_indexed(out, title, atoms, atoms->nr, NULL, x, v, box);
    fclose(out);
    break;
  case efTPR:
  case efTPB:
  case efTPA:
    gmx_fatal(FARGS,"Sorry, can not write a topology to %s",outfile);
    break;
  default:
    gmx_incons("Not supported in write_sto_conf");
  }
}

void get_stx_coordnum(char *infile,int *natoms)
{
  FILE *in;
  int ftp,tpxver,tpxgen;
  t_trxframe fr;

  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    get_coordnum(infile, natoms);
    break;
  case efG96:
    in=ffopen(infile,"r");
    fr.title = NULL;
    fr.natoms = -1;
    fr.atoms = NULL;
    fr.x = NULL;
    fr.v = NULL;
    fr.f = NULL;
    *natoms=read_g96_conf(in,infile,&fr);
    fclose(in);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    in=ffopen(infile,"r");
    get_pdb_coordnum(in, natoms);
    fclose(in);
    break;
  case efESP:
    *natoms = get_espresso_coordnum(infile);
    break;
  case efTPA:
  case efTPB:
  case efTPR: {
    t_tpxheader tpx;
    
    read_tpxheader(infile,&tpx,TRUE,&tpxver,&tpxgen);
    *natoms = tpx.natoms;
    break;
  }
  default:
    gmx_incons("Not supported in get_stx_coordnum");
  }
}

void read_stx_conf(char *infile, char *title,t_atoms *atoms, 
		   rvec x[],rvec *v, matrix box)
{
  FILE       *in;
  char       buf[256];
  t_topology *top;
  t_trxframe fr;
  int        i,ftp,natoms,i1;
  real       d,r1,r2;

  if (atoms->nr == 0)
    fprintf(stderr,"Warning: Number of atoms in %s is 0\n",infile);
  else if (atoms->atom == NULL)
    gmx_mem("Uninitialized array atom");
  
  ftp=fn2ftp(infile);
  switch (ftp) {
  case efGRO:
    read_whole_conf(infile, title, atoms, x, v, box);
    break;
  case efG96:
    fr.title = title;
    fr.natoms = atoms->nr;
    fr.atoms = atoms;
    fr.x = x;
    fr.v = v;
    fr.f = NULL;
    in = ffopen(infile,"r");
    read_g96_conf(in, infile, &fr);
    fclose(in);
    copy_mat(fr.box,box);
    break;
  case efPDB:
  case efBRK:
  case efENT:
    read_pdb_conf(infile, title, atoms, x, box, TRUE);
    break;
  case efESP:
    read_espresso_conf(infile,atoms,x,v,box);
    break;
  case efTPR:
  case efTPB:
  case efTPA: 
    snew(top,1);
    read_tpx(infile,&i1,&r1,&r2,NULL,box,&natoms,x,v,NULL,top);
    
    strcpy(title,*(top->name));
    /* Scalars */
    atoms->nr       = top->atoms.nr;
    atoms->nres     = top->atoms.nres;
    atoms->ngrpname = top->atoms.ngrpname;
    
    /* Arrays */
    if (!atoms->atom)
      snew(atoms->atom,atoms->nr);
    if (!atoms->atomname)
      snew(atoms->atomname,atoms->nr);
    if (!atoms->atomtype)
      snew(atoms->atomtype,atoms->nr);
    if (!atoms->atomtypeB)
      snew(atoms->atomtypeB,atoms->nr);
    for(i=0; (i<atoms->nr); i++) {
      atoms->atom[i]      = top->atoms.atom[i];
      atoms->atomname[i]  = top->atoms.atomname[i];
      atoms->atomtype[i]  = top->atoms.atomtype[i];
      atoms->atomtypeB[i] = top->atoms.atomtypeB[i];

    }
    
    if (!atoms->resname)
      snew(atoms->resname,atoms->nres);
    for(i=0; (i<atoms->nres); i++) 
      atoms->resname[i] = top->atoms.resname[i];
    
    if (!atoms->grpname)
      snew(atoms->grpname,atoms->ngrpname);
    for(i=0; (i<atoms->ngrpname); i++) 
      atoms->grpname[i] = top->atoms.grpname[i];
      
    for(i=0; (i<egcNR); i++) {
      atoms->grps[i].nr = top->atoms.grps[i].nr;
      if (atoms->grps[i].nr > 0) {
	snew(atoms->grps[i].nm_ind,atoms->grps[i].nr);
	memcpy(atoms->grps[i].nm_ind,top->atoms.grps[i].nm_ind,
	       atoms->grps[i].nr*sizeof(atoms->grps[i].nm_ind[0]));
      }
    }
      
    /* Ignore exclusions */
    
    break;
  default:
    gmx_incons("Not supported in read_stx_conf");
  }
}


