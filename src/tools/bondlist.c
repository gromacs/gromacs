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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cdist.h"

typedef struct list {
  char *residue;
  struct atombond {
    char *ai,*aj;
    real dist;
  } *atombond;
  struct atomangle {
    char *ai,*aj;
    real angle;
  } *atomangle;
  int natombond;
  int natomangle;
  struct list *next;
} bondlist;

static bondlist *lst=NULL;

static void newresidue(char *residue)
{
  bondlist *p;
  
  snew(p,1);
  p->residue    = strdup(residue);
  p->atombond   = NULL;
  p->atomangle  = NULL;
  p->natombond  = 0;
  p->natomangle = 0;
  p->next       = lst;  
  lst           = p;
  if (debug)
    fprintf(debug,"initialised element for residue %s\n",residue);
}

void read_O_dist(void)
{
#define BLEN 255
  static   char *fn= "refi_aa.dat";
  FILE     *fp;
  char     buf[BLEN+1],buf_ai[9],buf_aj[9],junk1[15],junk2[15];
  int      nres = 0;
  double   dist;
  real     junk3;

  fp = libopen(fn);
  fprintf(stderr,"Going to read %s\n",fn);
  fprintf(stderr,"Take note that the angle O C N+ will not be read from the");
  fprintf(stderr," correct residue\nThis is not a problem in the current");
  fprintf(stderr,"(980927) param-file while \nthis angle is const. 123.000\n");
  while (fgets2(buf,BLEN,fp) != NULL) {    
    switch (buf[0]) {
    case '#':
    case '.':
    case 't':
      continue;
      break;
    case 'r':
      sscanf(buf,"%s%s",buf_ai,buf_aj);
      if (strcmp(buf_ai,"residue") == 0 ) {
	nres++;
	/* O only uses HIS, Gromos has HISB... */
	if (strcmp(buf_aj,"HIS") == 0) strcpy(buf_aj,"HISB");
	newresidue(buf_aj);
      }
      break;
    case 'b':
      if ( buf[5] == 'd' ) {
	sscanf(buf,"%*s%s%s%lf",buf_ai,buf_aj,&dist);
	/*fprintf(stderr,"1: %s %s %g\n",buf_ai,buf_aj,dist);*/
	if (strcmp(buf_ai,"C-") == 0) strcpy(buf_ai,"C");
	if (strcmp(buf_ai,"CA-") == 0) strcpy(buf_ai,"CA");
	if (strcmp(buf_ai,"N+") == 0) strcpy(buf_ai,"N");
	if (strcmp(buf_aj,"C-") == 0) strcpy(buf_aj,"C");
	if (strcmp(buf_aj,"CA-") == 0) strcpy(buf_aj,"CA");
	if (strcmp(buf_aj,"N+") == 0) strcpy(buf_aj,"N");
	/*fprintf(stderr,"2: %s %s %g\n",buf_ai,buf_aj,dist);*/
     	lst->natombond++;
	lst->atombond = srenew(lst->atombond,lst->natombond);
	lst->atombond[lst->natombond-1].ai = strdup(buf_ai);
	lst->atombond[lst->natombond-1].aj = strdup(buf_aj);
	lst->atombond[lst->natombond-1].dist = (dist/10.0);
      }
      if ( buf[5] == 'a' ) {
	sscanf(buf,"%*s%s%*s%s%lf",buf_ai,buf_aj,&dist);
	/*fprintf(stderr,"1: %s %s %g\n",buf_ai,buf_aj,dist);*/
	if (strcmp(buf_ai,"C-") == 0) strcpy(buf_ai,"C");
	if (strcmp(buf_ai,"CA-") == 0) strcpy(buf_ai,"CA");
	if (strcmp(buf_ai,"N+") == 0) strcpy(buf_ai,"N");
	if (strcmp(buf_aj,"C-") == 0) strcpy(buf_aj,"C");
	if (strcmp(buf_aj,"CA-") == 0) strcpy(buf_aj,"CA");
	if (strcmp(buf_aj,"N+") == 0) strcpy(buf_aj,"N");
	/*fprintf(stderr,"2: %s %s %g\n",buf_ai,buf_aj,dist);*/
	lst->natomangle++;
	lst->atomangle = srenew(lst->atomangle,lst->natomangle);
	lst->atomangle[lst->natomangle-1].ai = strdup(buf_ai);
	lst->atomangle[lst->natomangle-1].aj = strdup(buf_aj);
	lst->atomangle[lst->natomangle-1].angle = dist;
      }
      break;
    default:
      fprintf(stderr,"Ooops! Unexpected first character in line:\n%s",buf);
      continue;
      break;
    }
  }
  fclose(fp);
  fprintf(stderr,"Read distances for %d residues from file %s \n",nres,fn);
}

static real do_specialbond(char *name1,char *name2,char *res)
{
  bondlist *p; 
  int      i;
  
  if (debug)
    fprintf(debug,"Looking for dist between %s %s in res. %s\n",
	    name1,name2,res);

  for ( p = lst; p != NULL; p = p->next)
    if (strcmp(res,p->residue) == 0) 
      break;
    
  /* Here either p == NULL or found */
  for (i=0; p!= NULL && (i<p->natombond) ; i++) {
    if (((strcmp(name1,p->atombond[i].ai) == 0) && 
	 (strcmp(name2,p->atombond[i].aj) == 0)) ||
	((strcmp(name1,p->atombond[i].aj) == 0) && 
	 (strcmp(name2,p->atombond[i].ai) == 0))) {
      if (debug) {
	fprintf(debug,"Found dist. %s %s %s",name1,name2,res);
	fprintf(debug,"   %g\n",p->atombond[i].dist);
      }
      return p->atombond[i].dist;
    }
  }
  if (debug)
    fprintf(debug,"Didn't find dist. %s %s %s\n",name1,name2,res);
    
  return 0.0;
}	  

static real do_specialangle(char *name1,char *name2,char *res)
{

  bondlist *p; 
  int      i;

  if (debug)
    fprintf(debug,"Looking for angle between %s %s in res. %s\n",
	    name1,name2, res);

  for ( p = lst; p != NULL; p = p->next)
    if (strcmp(res,p->residue) == 0) 
      break;
      
  /* Here either p == NULL or found */
  for (i=0; p != NULL && (i<p->natomangle) ; i++) {
    if (((strcmp(name1,p->atomangle[i].ai) == 0) && 
	 (strcmp(name2,p->atomangle[i].aj) == 0)) ||
	((strcmp(name1,p->atomangle[i].aj) == 0) && 
	 (strcmp(name2,p->atomangle[i].ai) == 0))) {
      if (debug) {
	fprintf(debug,"Found angle. %s %s %s",name1,name2,res);
	fprintf(debug,"   %g\n",p->atomangle[i].angle);
      }
      return p->atomangle[i].angle;
    }
  }
  if (debug)
    fprintf(debug,"Didn't find angle %s %s %s\n",name1,name2,res);
    
  return 0.0;
}	  

real lookup_bondlength_(int ai,int aj,t_ilist ilist[],
			t_iparams iparams[],bool bFail,t_atoms *atoms,
			char *file,int line)
{
  int  dodist[] = { F_BONDS, F_MORSE, F_SHAKE, F_G96BONDS };
  int  i,j,type,ftype,aai,aaj,nra;
  real blen=NOBOND;
  real spclen;
  
  if (debug) 
    fprintf(debug,"call do_spec with %s (%d) %s (%d) res %s \n",
	    *atoms->atomname[ai],ai,*atoms->atomname[aj],aj,
	    *(atoms->resname[atoms->atom[ai].resnr]));
	    
  spclen = do_specialbond(*(atoms->atomname[ai]),*(atoms->atomname[aj]),
			  *(atoms->resname[atoms->atom[aj].resnr]));
			  
  /* Comment; if you're using distances from e.g. O, it is non-trivial
     which residuename (ai or aj) you call do_special with, in order
     to get the peptide bonds correct. */
  if (spclen != 0) {
    if (debug) 
      fprintf(debug,"%s %d %d %g\n",*atoms->resname[atoms->atom[ai].resnr],
	      ai,aj,spclen);
    blen = spclen;
  }
  else {
    if (debug)
      fprintf(debug,"No spclen found for %s (%d) %s (%d) res %s \n",
	      *atoms->atomname[ai],ai,*atoms->atomname[aj],aj,
	      *(atoms->resname[atoms->atom[ai].resnr]));
    for(j=0; (j<asize(dodist)) && (blen == NOBOND); j++) {
      ftype = dodist[j];
      nra   = interaction_function[ftype].nratoms;
    
      for(i=0; (i<ilist[ftype].nr) && (blen == NOBOND); i+=nra+1) {
	type  = ilist[ftype].iatoms[i];
	aai    = ilist[ftype].iatoms[i+1];
	aaj    = ilist[ftype].iatoms[i+2];
      
	if (((aai == ai) && (aaj == aj)) || ((aaj == ai) && (aai == aj))) {
	  switch (ftype) {
	  case F_G96BONDS:
	    blen = sqrt(iparams[type].harmonic.rA);
	    break;
	  case F_BONDS: 
	    blen = iparams[type].harmonic.rA;
	    break;
	  case F_MORSE:
	    blen = iparams[type].morse.b0;
	    break;
	  case F_SHAKE:
	    blen = iparams[type].shake.dA;
	    break;
	  default:
	    break;
	  }
	}
      }
    }
  }
  if (blen == NOBOND) {
    if (bFail)
      gmx_fatal(FARGS,"No bond between atoms %d and %d (called from %s line %d)\n",
		  ai,aj,file,line);
    else
      return NOBOND;
  }
  
  return blen*10;
  /* This should be the only conversion from nanometer to Angstrom */
}

real lookup_angle_(int ai,int aj,int ak,t_ilist ilist[],
		   t_iparams iparams[],t_atoms *atoms,
		   char *file,int line)
{
  int  ft[] = { F_ANGLES, F_G96ANGLES };
  int  i,j,type,ff,ftype,aai,aaj,aak,nra;
  real angle;

  /* First check the Engh & Huber database */
  angle = DEG2RAD*do_specialangle(*atoms->atomname[ai],*atoms->atomname[ak],
				  *atoms->resname[atoms->atom[ak].resnr]);
  /* Now check the topology */
  for(ff=0; (ff<asize(ft)) && (angle == 0); ff++) {
    ftype = ft[ff];
    nra   = interaction_function[ftype].nratoms;
    
    for(i=0; (i<ilist[ftype].nr) && (angle == 0); i+=nra+1) {
      type  = ilist[ftype].iatoms[i];
      aai   = ilist[ftype].iatoms[i+1];
      aaj   = ilist[ftype].iatoms[i+2];
      aak   = ilist[ftype].iatoms[i+3];
      
      if (((aai == ai) && (aaj == aj) && (aak == ak)) || 
	  ((aak == ai) && (aaj == aj) && (aai == ak))) {
	if (ftype == F_ANGLES)
	  angle = DEG2RAD*iparams[type].harmonic.rA;
	else if (ftype == F_G96ANGLES)
	  angle = acos(iparams[type].harmonic.rA);
	else
	  gmx_fatal(FARGS,"Unknown angletype %s in %s, line %d",
		      interaction_function[ftype].longname,__FILE__,__LINE__);
      }
    }
  }
  if (angle == 0) {
    fprintf(stderr,
	    "Warning: no known angle between atoms %d, %d, %d. Using 109.47\n",
	    ai,aj,ak);
    angle = 109.47*DEG2RAD;
  }    
  return angle;
}

real angle_length_(int ai,int aj,int ak,real theta,
		   t_ilist ilist[],t_iparams iparams[],t_atoms *atoms,
		   char *file,int line)
{
  real rij,rjk;
  
  rij = lookup_bondlength_(ai,aj,ilist,iparams,TRUE,atoms,file,line);
  rjk = lookup_bondlength_(aj,ak,ilist,iparams,TRUE,atoms,file,line);
  
  if (debug)
    fprintf(debug,"angle_length uses %g %g and angle %g\n",rij,rjk,theta);

  return sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(DEG2RAD*theta));
}

/*
  static void special_bonds(t_dist *d,t_atoms *atoms,
  char *res,int nbonds,t_bonddef bdef[])
  {
  int i,j,k,resnr,na,ai,aj,ndist;
  
  resnr = -1;
  for(i=0; (i<atoms->nr); i++) {
    if (strcmp(*atoms->resname[atoms->atom[i].resnr],res) == 0) {
      resnr = atoms->atom[i].resnr;
      for(k=0; (k<nbonds); k++) {
        ai = aj = -1;
        for(j=i; ((j<atoms->nr) && (atoms->atom[j].resnr == resnr)); j++) {
          if (strcmp(bdef[k].ai,*atoms->atomname[j]) == 0)
            if (ai != -1)
              gmx_fatal(FARGS,"Atom %s multiply defined in res %s %d",
                          bdef[k].ai,res,resnr);
            else
              ai = j;
          if (strcmp(bdef[k].aj,*atoms->atomname[j]) == 0)
            if (aj != -1)
              gmx_fatal(FARGS,"Atom %s multiply defined in res %s %d",
                          bdef[k].aj,res,resnr);
            else
              aj = j;
        }
        if ((ai != -1) && (aj != -1)) {
	  lb=(1.0-pep_margin)*bdef[k].len;
	  ub=(1.0+pep_margin)*bdef[k].len;
	  if (((weight[ai] != 0.0) || (weight[aj] != 0.0)) &&
	      (!dist_set(d,natoms,ai,aj))) {
	    set_dist(d,natoms,ai,aj,lb,ub,len);
	    ndist++;
	  } 
        }
      }
      } */
    /* Optimizations ... */
/*  }
   fprintf(stderr,"There were %d special distances defined (GROMOS-override)",
	   ndist);
 }
*/

  
void pdih_lengths_(int ai,int aj,int ak,int al,
		   t_ilist ilist[],t_iparams iparams[],
		   real *lb,real *ub,t_atoms *atoms,char *file,int line)
{
  /* Returns the length corresponding to a cis dihedral */
  real rij,rjk,rkl,rik;
  real th1,th2,th3;
  real half_pi = M_PI*0.5;
  
  rij = lookup_bondlength_(ai,aj,ilist,iparams,TRUE,atoms,file,line);
  rjk = lookup_bondlength_(aj,ak,ilist,iparams,TRUE,atoms,file,line);
  rkl = lookup_bondlength_(ak,al,ilist,iparams,TRUE,atoms,file,line);
  th1 = lookup_angle_(ai,aj,ak,ilist,iparams,atoms,file,line);
  th2 = lookup_angle_(aj,ak,al,ilist,iparams,atoms,file,line);

  /* Compute distance from i to k using law of cosines */
  rik = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(th1));
  
  /* Compute angle th3 using law of sines */
  th3 = asin(rij*sin(th1)/rik);
  
  /* Compute cis length using law of cosines */
  *lb = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(th2-th3));
  
  /* Compute trans length using law of cosines */
  *ub = sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(th2+th3));
}

real idih_lengths(int ai,int aj,int ak,int al,
		  t_ilist ilist[],t_iparams iparams[],t_atoms *atoms)
{
  /* Returns the length corresponding to a cis dihedral */
  real rij,rjk,rkl,rik;
  real th1,th2,th3,cis;
  real half_pi = M_PI*0.5;
  real lb;

  lb = 0.0;
    
  if ((rij = lookup_bondlength(ai,aj,ilist,iparams,FALSE,atoms)) == NOBOND)
    return lb;
  if ((rjk = lookup_bondlength(aj,ak,ilist,iparams,FALSE,atoms)) == NOBOND)
    return lb;
  if ((rkl = lookup_bondlength(ak,al,ilist,iparams,FALSE,atoms)) == NOBOND)
    return lb;
    
  if (debug)
    fprintf(debug,"Found an improper: %d %d %d %d\n",ai,aj,ak,al);
  th1 = lookup_angle(ai,aj,ak,ilist,iparams,atoms);
  th2 = lookup_angle(aj,ak,al,ilist,iparams,atoms);

  /* Compute distance from i to k using law of cosines */
  rik = sqrt(rij*rij+rjk*rjk-2.0*rij*rjk*cos(th1));
  
  /* Compute angle th3 using law of sines */
  th3 = asin(rij*sin(th1)/rik);
  
  /* Compute cis length using law of cosines */
  return sqrt(rik*rik+rkl*rkl-2.0*rik*rkl*cos(th2-th3));
}
