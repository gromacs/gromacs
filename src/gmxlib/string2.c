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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_string2_c = "$Id$";

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <sys/types.h>
#ifndef NO_PWUID
#include <pwd.h>
#endif
#include <time.h>

#include "typedefs.h"
#include "smalloc.h"
#include "fatal.h"
#include "macros.h"
#include "string2.h"

int continuing(char *s)
/* strip trailing spaces and if s ends with a CONTINUE remove that too.
 * returns TRUE if s ends with a CONTINUE, FALSE otherwise.
 */
{
  int sl;

  rtrim(s);
  sl = strlen(s);
  if ((sl > 0) && (s[sl-1] == CONTINUE)) {
    s[sl-1] = 0;
    return TRUE;
  }
  else
    return FALSE;
}

char *fgets2(char *line, int n, FILE *stream)
/* This routine reads a string from stream of max length n
 * and zero terminated, without newlines
 * line should be long enough (>= n)
 */
{
  char *c;
  if (fgets(line,n,stream)==NULL) return NULL;
  if ((c=strchr(line,'\n'))!=NULL) *c=0;
  return line;
}

void strip_comment (char *line)
{
  char *c;

  if (!line)
    return;

  /* search for a comment mark and replace it by a zero */
  if ((c = strchr(line,COMMENTSIGN)) != NULL) 
    (*c) = 0;
}

void upstring (char *str)
{
  int i;

  for (i=0; (i < (int)strlen(str)); i++) 
    str[i] = toupper(str[i]);
}

void ltrim (char *str)
{
  char *tr;
  int c;

  if (!str)
    return;

  tr = strdup (str);
  c  = 0;
  while ((tr[c] == ' ') || (tr[c] == '\t'))
    c++;

  strcpy (str,tr+c);
  free (tr);
}

void rtrim (char *str)
{
  int nul;

  if (!str)
    return;

  nul = strlen(str)-1;
  while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
    str[nul] = '\0';
    nul--;
  }
}

void trim (char *str)
{
  ltrim (str);
  rtrim (str);
}

void nice_header (FILE *out,char *fn)
{
  time_t clock;
#ifndef NO_PWUID
  struct passwd *pw;
  uid_t uid;
#endif

  /* Print a nice header above the file */
  clock = time (0);
  fprintf (out,"%c\n",COMMENTSIGN);
#ifndef NO_PWUID
  uid = getuid();
  if ((pw = getpwuid(uid)) != NULL)
    fprintf (out,"%c\tUser %s (%d)\n",COMMENTSIGN,pw->pw_name,(int) uid);
  else
    fprintf (out,"%c\tUser %d not in password file\n",COMMENTSIGN,(int) uid);
#endif
  fprintf (out,"%c\t%s",COMMENTSIGN,ctime(&clock));
  fprintf (out,"%c\t%s\n",COMMENTSIGN,fn);
  fprintf (out,"%c\n",COMMENTSIGN);
}

int gmx_strcasecmp(const char *str1, const char *str2)
{
  char ch1,ch2;
  
  do
    {
      do
	ch1=toupper(*(str1++));
      while ((ch1=='-') || (ch1=='_'));
      do 
	ch2=toupper(*(str2++));
      while ((ch2=='-') || (ch2=='_'));
      if (ch1!=ch2) return (ch1-ch2);
    }
  while (ch1);
  return 0; 
}

char *gmx_strdup(const char *src)
{
  char *dest;

  if ((dest=(char *)malloc((size_t)(strlen(src)+1))) != NULL)
    strcpy(dest,src);
  return dest;
}

char *wrap_lines(char *buf,int line_width, int indent)
{
  char *b2;
  int i,i0,i2,j,b2len,lspace,l2space;
  bool bFirst;

  b2=NULL;
  b2len=strlen(buf)+1;
  snew(b2,b2len);
  i0=0;
  i2=0;
  bFirst=TRUE;
  do {
    for(i=i0; (i<i0+line_width) && (buf[i]); i++) {
      b2[i2++] = buf[i];
      if (buf[i] == ' ') {
        lspace = i;
	l2space = i2-1;
      }
      if (buf[i]=='\n') 
	i0=i+1;
    }
    if (buf[i]) {
      b2[l2space] = '\n';
      i0 = lspace+1;
      i2 = l2space+1;
      if (indent) {
	if (bFirst) {
	  line_width-=indent;
	  bFirst=FALSE;
	}
	b2len+=indent+1;
	srenew(b2, b2len);
	for(j=0; (j<indent); j++)
	  b2[i2++]=' ';
      }
    }
  } while (buf[i]);
  b2[i2] = '\0';
  
  return b2;
}



