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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_genconf_c = "$Id$";

#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "topexcl.h"
#include "x2top.h"

t_nm2type *rd_nm2type(char *ff,int *nnm)
{
  FILE      *fp;
  bool      bCont;
  char      libfilename[128];
  char      format[128],f1[128];
  char      buf[1024],elem[16],type[16],nbbuf[16],**newbuf;
  int       i,n,nb,nnnm,line=1;
  t_nm2type *nm2t=NULL;
  
  sprintf(libfilename,"%s.n2t",ff);
  fp = libopen(libfilename);
  
  nnnm = 0;
  do {
    /* Read a line from the file */
    bCont = (fgets2(buf,1023,fp) != NULL); 
    
    if (bCont) {
      /* Remove comment */
      strip_comment(buf);
      if ((n = sscanf(buf,"%s%s%d",elem,type,&nb)) == 3) {
	/* If we can read the first three, there probably is more */
	if (nb > 0) {
	  snew(newbuf,nb);
	  strcpy(format,"%*s%*s%*d");
	  for(i=0; (i<nb); i++) {
	    /* Complicated format statement */
	    strcpy(f1,format);
	    strcat(f1,"%s");
	    if (sscanf(buf,f1,nbbuf) != 1)
	      fatal_error(0,"Error on line %d of %s",line,libfilename);
	    newbuf[i] = strdup(nbbuf);
	    strcat(format,"%*s");
	  }
	}
	else
	  newbuf = NULL;
	srenew(nm2t,nnnm+1);
	nm2t[nnnm].elem   = strdup(elem);
	nm2t[nnnm].type   = strdup(type);
	nm2t[nnnm].nbonds = nb;
	nm2t[nnnm].bond   = newbuf;
	nnnm++;
      }
      line++;
    }
  } while(bCont);
  fclose(fp);
  
  *nnm = nnnm;
  
  return nm2t;
}

void dump_nm2type(FILE *fp,int nnm,t_nm2type nm2t[])
{
  int i,j;
  
  fprintf(fp,"; nm2type database\n");
  for(i=0; (i<nnm); i++) {
    fprintf(fp,"%-8s%-8s%-4d",nm2t[i].elem,nm2t[i].type,nm2t[i].nbonds);
    for(j=0; (j<nm2t[i].nbonds); j++)
      fprintf(fp,"%-5s",nm2t[i].bond[j]);
    fprintf(fp,"\n");
  }
}

char *nm2type(int nnm,t_nm2type nm2t[],char *elem,int nbonds)
{
  int i;
  
  for(i=0; (i<nnm); i++) {
    if ((strcasecmp(nm2t[i].elem,elem) == 0) &&
	(nm2t[i].nbonds == nbonds))
      return nm2t[i].type;
  }
	      
  return NULL;
}
     
