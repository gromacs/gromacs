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

#include <stdlib.h>
#include <string.h>
#include <metacode.h>
#include <stdarg.h>
#include <ctype.h>

/* This file is NOT threadsafe, but it is only used to create
 * the innerloops during the build process, so it will never be
 * executed by multiple threads.
 */

#define MAXCODESIZE 1000000
#define MAXDECL 1000
#define F77IND 6
#define CIND 2

#ifdef DOUBLE
int prec   = 8;
#else
int prec   = 4;
#endif

#define REAL_FORMAT "%.16f"

int IND=F77IND;

char header[10000]; /* Should be enough for a few comment lines
		     * and the function name
		     */
FILE *output;
char *codebuffer=NULL;
decl_t *decl_list;

int ndecl=0;
int nargs=0;  /* the first nargs entries in the declaration list
	       * are function arguments 
	       */

bool bC       = TRUE;

void init_metacode()
{
  static bool first=TRUE;
  int i;

  IND = bC ? CIND : F77IND;
  
  /* sanity check if the buffers are initialized more than once. 
   * They are also emptied upon flushing.
   */
  if(first) {
    decl_list=(decl_t *)malloc(sizeof(decl_t)*MAXDECL);
    codebuffer=(char *)malloc(MAXCODESIZE);
    first=FALSE;
  }
  ndecl=0;
  nargs=0;
  codebuffer[0]=0;
  header[0]=0;
}


void add_to_buffer(char *buffer,char *term)
{
  if(strlen(buffer)>0)
    strcat(buffer,"+");
  strcat(buffer,term);  
}


static bool findname(char *buffer,char *name)
{
  /* This routine returns true if name is found in
   * buffer, and the surrounding characters are non
   * alphanumeric. (i.e. it is not a substring of
   * a longer word). It might catch things in comments
   * if you choose to keep them with -DKEEP_COMMENTS, though,
   * in which case you'll get warnings about unused variables.
   */
  char *ch;
  
  /* the first character will be a space,
   * so it is safe to start at buffer+1.
   */
  ch=buffer;
  
  do {
    if((ch=strstr(ch+1,name))!=NULL) {
      /* Found something. But is it a full variable
       * or a substring in a larger name?
       * Check if the prev/next chars are alphanumeric!
       */
      if(!isalnum(*(ch+strlen(name))) && !isalnum(*(ch-1)))
	return TRUE; /* found it! */
    }
  } while(ch!=NULL);

  return FALSE; /* no hit */
}


void flush_buffers(void)
{
  int i;
  int nwritten;
  char *ch;
  /* scan the code output buffer for arguments and
   * variables. Remove all that are unreferenced!
   */
  for(i=0;i<ndecl;i++) 
    decl_list[i].breferenced=findname(codebuffer,decl_list[i].name);    
  
  /* write the function name (and start argument list) */
  fprintf(output,header);
  
  nwritten=0;
  /* write out all referenced FUNCTION ARGUMENTS */
  for(i=0;i<nargs;i++) {
    if(!decl_list[i].breferenced)
      continue;
    if(nwritten) /* separate from earlier arg with comma */
      fprintf(output,",");
    nwritten++;
    if(bC) 
      fprintf(output,"\n\t%15s %s%s",decl_list[i].typename,decl_list[i].name,
	      decl_list[i].bvector ? "[]" : "");
    else
      fprintf(output,"\n" FCON "  %s",decl_list[i].name);
  }

  /* finish argument list, start function */
  fprintf(output,")\n");
  if(bC)
    fprintf(output,"{\n");
  else
    newline();

  /* declare arguments in fortran */
  if(!bC) {
    fprintf(output,"%simplicit none\n",indent());
    for(i=0;i<nargs;i++) {
      if(!decl_list[i].breferenced)
	continue;
      fprintf(output,"%s%-10s %s%s\n",indent(),
	      decl_list[i].typename,decl_list[i].name,
	      decl_list[i].bvector ? "(*)" : "");
    }
  }
    
  /* declare all non-removed VARIABLES following the arguments */
  for(i=nargs;i<ndecl;i++) {
    if(!decl_list[i].breferenced)
      continue;
    if(bC) {
      fprintf(output,"%s%-10s %s%s",indent(),
	      decl_list[i].typename,decl_list[i].name,
	      decl_list[i].bvector ? "[]" : "");
      if(decl_list[i].bconst)
	fprintf(output,"  =  %s;\n",decl_list[i].constval);
      else
	fprintf(output,";\n");
    } else {
      fprintf(output,"%s%-10s %s%s\n",indent(),
	      decl_list[i].typename,decl_list[i].name,
	      decl_list[i].bvector ? "(*)" : "");
    }
  }

  /* assign fortran parameters */
  if(!bC) {
    for(i=nargs;i<ndecl;i++) 
      if(decl_list[i].breferenced && decl_list[i].bconst)
	fprintf(output,"%sparameter (%s  =  %s)\n",
		indent(),decl_list[i].name,decl_list[i].constval);
  }
  /* write one huge piece of function code */
  fprintf(output,codebuffer);
  
  /* empty the buffers to prepare for next routine... */

  ndecl=nargs=0;
  codebuffer[0]=0;
  header[0]=0;
}


/* Return the correct indentation as a string */
char *indent(void)
{
  static char indbuf[1024];
  int i,n;
  
  n = max(0,IND);
  for(i=0; (i<n); i++)
    indbuf[i] = ' ';
  indbuf[i] = '\0';
  
  return indbuf;
}

void fortran_splitline(char *line)
{
  char tmpbuf[1024],linebuf[1000];
  int i,j,maxlen;
  
  maxlen=strlen(line);
  
  i=0;
  while(i+72-IND<maxlen) {
    j=i+71-IND;
    if (j>=maxlen)
      j=maxlen;
    
    while(j>(i+1)) {
      if(line[j]=='+' ||
	 line[j]=='-' ||
	 line[j]=='/' ||
	 (line[j]=='*' && line[j-1]!='*') ||  /* dont split pows */
	 line[j]==' ')
	break;
      j--;
    }
    if(j==i) {
      printf("Error: Couldn't break this line:\n%s\n",line);
      exit(1);
    }
    strncpy(tmpbuf,line+i,j-i+1);
    tmpbuf[j-i]=0;
    
    strcat(codebuffer,tmpbuf);
    newline();
    i=j;
    strcat(codebuffer,FCON "  ");
  }
  strcat(codebuffer,line+i);
}


/* Print a line of code to the output file */
void code(char *fmt, ...)
{
  va_list ap;
  int d;
  double f;
  char *s;
  char tmpbuf[1024];
  char outbuf[1024];

  sprintf(outbuf,"%s",indent());
  
  va_start(ap,fmt);
  
  while(*fmt) {
    if(*fmt=='%') {
      switch(*(++fmt)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmpbuf,"%d",d);
	strcat(outbuf,tmpbuf);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmpbuf,REAL_FORMAT,f);
	strcat(outbuf,tmpbuf);
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmpbuf,s);
	strcat(outbuf,tmpbuf);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else {
      sprintf(tmpbuf,"%c",*fmt);
      strcat(outbuf,tmpbuf);
    }
    fmt++;
  }
  va_end(ap);
  if(bC)
    strcat(codebuffer,outbuf);
  else
    fortran_splitline(outbuf);
  newline();
}




void newline(void)
{
  strcat(codebuffer,"\n");
}

/* Add a comment - might just come in handy for debugging,
   but we don't need it in production level code */
void comment(char *s)
{
  char buf[512];
#ifdef KEEP_COMMENTS
  if (bC)
    sprintf(buf,"\n%s/* %s */\n",indent(),s);
  else {
    IND--;
    sprintf(buf,"\nC%s%s\n",indent(),s);
    IND++;
  }
  strcat(codebuffer,buf);
#endif
}


/* Define a new floating-point variable */
void declare_real(char *name)
{

  if (bC)
    strcpy(decl_list[ndecl].typename,"real");
  else
    sprintf(decl_list[ndecl].typename,"real*%d",prec);

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_real4(char *name)
{

  if (bC)
    strcpy(decl_list[ndecl].typename,"float");
  else
    sprintf(decl_list[ndecl].typename,"real*4");

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_real_vector(char *name)
{

  if (bC)
    strcpy(decl_list[ndecl].typename,"real");
  else
    sprintf(decl_list[ndecl].typename,"real*%d",prec);

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=TRUE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_const_real_vector(char *name)
{

  if (bC)
    strcpy(decl_list[ndecl].typename,"const real");
  else
    sprintf(decl_list[ndecl].typename,"real*%d",prec);

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=TRUE;
  decl_list[ndecl].bconst=TRUE;
  ndecl++;
}

void declare_intreal(char *name)
{
#ifdef DOUBLE
  declare_int8(name);
#else
  declare_int4(name);
#endif
}


void declare_const_real(char *name,double val)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"const real");
  else
    sprintf(decl_list[ndecl].typename,"real*%d",prec);

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE; /* cant have const vectors */
  decl_list[ndecl].bconst=TRUE;
  sprintf(decl_list[ndecl].constval,REAL_FORMAT,val);
  ndecl++;
}

void declare_const_int(char *name,int val)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"const int");
  else
    sprintf(decl_list[ndecl].typename,"integer*4");

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE; /* cant have const vectors */
  decl_list[ndecl].bconst=TRUE;
  sprintf(decl_list[ndecl].constval,"%d",val);
  ndecl++;
}


void declare_int(char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"int");
  else
    sprintf(decl_list[ndecl].typename,"integer*%d",(int)sizeof(int));

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_int_vector(char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"int");
  else
    sprintf(decl_list[ndecl].typename,"integer*%d",(int)sizeof(int));

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=TRUE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_const_int_vector(char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"const int");
  else
    sprintf(decl_list[ndecl].typename,"integer*%d",(int)sizeof(int));

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=TRUE;
  decl_list[ndecl].bconst=TRUE;
  ndecl++;
}

void declare_int4(char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"int");
  else
    strcpy(decl_list[ndecl].typename,"integer*4");

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}


void declare_int8(char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,"long long");
  else
    strcpy(decl_list[ndecl].typename,"integer*8");

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}

void declare_other(char *typename,char *name)
{
  if (bC)
    strcpy(decl_list[ndecl].typename,typename);
  else
    strcpy(decl_list[ndecl].typename,typename);

  strcpy(decl_list[ndecl].name,name);
  
  decl_list[ndecl].breferenced=TRUE;
  decl_list[ndecl].bvector=FALSE;
  decl_list[ndecl].bconst=FALSE;
  ndecl++;
}


/* Cray vector pragma */
void vector_pragma(void)
{
#ifdef CRAY_PRAGMA
  if (bC)
    strcat(codebuffer,"#pragma ivdep\n");
  else
    strcat(codebuffer,"cdir$ivdep\n");
#endif
}

char *_array(char *a, char *idx, ...)
{
  char arrtmp[1000],idxtmp[1000],tmp[1000];
  va_list ap;
  int d;
  char *s,c;
  
  arrtmp[0]=idxtmp[0]=0;

  va_start(ap,idx);
  
  while(*a) {
    if(*a=='%') {
      switch(*(++a)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to _array()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*a);
    a++;
    strcat(arrtmp,tmp);
  }
  
  while(*idx) {
    if(*idx=='%') {
      switch(*(++idx)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to _array()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*idx);
    idx++;
    strcat(idxtmp,tmp);
  }

  va_end(ap);

  sprintf(tmp,"%s%c%s%c",arrtmp, bC ? '[' : '(',idxtmp, bC ? ']' : ')');
  /* ok, now we got the array reference in tmp. But there might be
   * some stupid things which need to be removed. First, if we add
   * a negative offset of e.g. -1 somewhere, we will get a "+-1" which 
   * is bad... remove the minus sign:
   */
  if((s=strstr(tmp,"+-"))!=NULL) {
    strcpy(arrtmp,s+1); /* copy to tmparray */
    strcpy(s,arrtmp);   /* copy back */
  }
     
  /* It is also stupid to add a zero offset. Kill that cat! */
  if((s=strstr(tmp,"+0"))!=NULL) {
    strcpy(arrtmp,s+2); /* copy to tmparray */
    strcpy(s,arrtmp);   /* copy back */
  }

  return strdup(tmp);
}
    

void file_error(char *fn)
{
  fprintf(stderr,"Error creating file %s\n",fn);
  exit(-1);
}


void _p_state(char *left,char *right,char *symb)
{
  char buf[512];

  if (bC) {
    if (IND+16+3+strlen(right) > 78) {
      sprintf(buf,"%s%-16s %2s \n",indent(),left,symb);
      strcat(codebuffer,buf);
      IND+=2;
      sprintf(buf,"%s%s;\n",indent(),right);
      strcat(codebuffer,buf);
      IND-=2;
    }
    else {
      sprintf(buf,"%s%-16s %2s %s;\n",indent(),left,symb,right);
      strcat(codebuffer,buf);
    }
  }
  else {
    if (IND+16+3+strlen(right) > 72) {
      sprintf(buf,"%s%-16s = \n%s",indent(),left,FCON);
      strcat(codebuffer,buf);
      IND-=6-3;
      code(right);
      IND+=6-3;
    }
    else {
      sprintf(buf,"%s%-16s = %s\n",indent(),left,right); 
      strcat(codebuffer,buf);
    }
  }
}


void assign(char *left,char *right, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;

  va_start(ap,right);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }
  
  while(*right) {
    if(*right=='%') {
      switch(*(++right)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*right);
    right++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");
}


void increment(char *left,char *right, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;
 
  va_start(ap,right);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }

  strcpy(rtmp,ltmp);
  strcat(rtmp," + ");

  while(*right) {
    if(*right=='%') {
      switch(*(++right)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*right);
    right++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");

}


void decrement(char *left,char *right, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;

  va_start(ap,right);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }

  strcpy(rtmp,ltmp);
  strcat(rtmp," - ");

  while(*right) {
    if(*right=='%') {
      switch(*(++right)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*right);
    right++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");
}



void add(char *left,char *r1,char *r2, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;

  va_start(ap,r2);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }
  
  while(*r1) {
    if(*r1=='%') {
      switch(*(++r1)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r1);
    r1++;
    strcat(rtmp,tmp);
  }
  strcat(rtmp," + ");

  while(*r2) {
    if(*r2=='%') {
      switch(*(++r2)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r2);
    r2++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");
}


void subtract(char *left,char *r1,char *r2, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;

  va_start(ap,r2);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }
  
  while(*r1) {
    if(*r1=='%') {
      switch(*(++r1)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r1);
    r1++;
    strcat(rtmp,tmp);
  }
  strcat(rtmp," - ");

  while(*r2) {
    if(*r2=='%') {
      switch(*(++r2)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r2);
    r2++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");
}



void multiply(char *left,char *r1,char *r2, ...)
{
  char ltmp[1000],rtmp[1000],tmp[1000];
  va_list ap;
  int d;
  double f;
  char *s,c;
  
  ltmp[0]=rtmp[0]=0;

  va_start(ap,r2);
  
  while(*left) {
    if(*left=='%') {
      switch(*(++left)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*left);
    left++;
    strcat(ltmp,tmp);
  }
  
  while(*r1) {
    if(*r1=='%') {
      switch(*(++r1)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r1);
    r1++;
    strcat(rtmp,tmp);
  }
  strcat(rtmp," * ");

  while(*r2) {
    if(*r2=='%') {
      switch(*(++r2)) {
      case 'd':
	d = va_arg(ap, int);
	sprintf(tmp,"%d",d);
	break;
      case 'f':
	f = va_arg(ap, double);
	sprintf(tmp,REAL_FORMAT,f);		
	break;
      case 's':
	s = va_arg(ap, char *);
	sprintf(tmp,s);
	break;
      case 'c':
	c = va_arg(ap, int);
	sprintf(tmp,"%c",c);
	break;
      default:
	fprintf(stderr,"Error, unsupported format supplied to code()\n");
	exit(-1);
	break;
      }
    } else
      sprintf(tmp,"%c",*r2);
    r2++;
    strcat(rtmp,tmp);
  }

  va_end(ap);
  _p_state(ltmp,rtmp,"=");

}



void edit_warning(char *fn)
{
  if(bC) 
    fprintf(output,
	  "  /**********************************************************\n" 
	  "   *    This code is generated automatically by %s\n"
	  "   *                DO NOT EDIT THIS FILE\n"
	  "   *     Erik Lindahl, David van der Spoel 1999-2000\n"
	  "   **********************************************************/"
	  "\n",fn);
  else
    fprintf(output,
	  "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
	  "C       This code is generated automatically by %s\n" 
	  "C                   DO NOT EDIT THIS FILE\n"
	  "C        Erik Lindahl, David van der Spoel 1999-2000\n"
	  "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
	  fn);
}

void closeit(void)
{
  if (bC)
    strcat(codebuffer,"}\n\n");
  else {
    code("return");
    code("end");
  }
}

void usage(int argc,char *argv[])
{
  fprintf(stderr,"Usage: %s language\n",argv[0]);
  fprintf(stderr,"\tAvailable languages:  c  fortran\n");
  exit(-1);
}

int count_lines(char *fn)
{
  FILE *fp;
  int nl=0;
  char buf[1024];

  if ((fp = fopen(fn,"r")) == NULL) {
    perror(fn);
    exit(1);
  }

  while (fgets(buf,255,fp) != NULL)
    nl++;
  fclose(fp);
  
  return nl;
}


void start_loop(char *lvar,char *from,char *to)
{
  if (bC) 
    code("for(%s=%s; (%s<%s); %s++) {", lvar,from,lvar,to,lvar);
  else 
    code("do %s=%s,%s",lvar,from,to);
  
  IND += 2;
}


void start_stride_loop(char *lvar,char *from,char *to, char *stride)
{
  if (bC) 
    code("for(%s=%s; (%s<%s); %s+=%s) {", lvar,from,lvar,to,lvar,stride);
  else 
    code("do %s=%s,%s,%s",lvar,from,to,stride);

  IND += 2;
}



void end_loop(void)
{
  IND -= 2;

  if (bC)
    code("}");
  else
    code("end do");
}

void start_if(char *cond)
{
  if (bC) 
    code("if(%s) {", cond);
  else 
    code("if (%s) then",cond);
  
  IND += 2;  
}

void do_else()
{
  IND -= 2;

  if (bC)
    code("} else {");
  else
    code("else");

  IND += 2;
}

    

void end_if(void)
{
  IND -= 2;

  if (bC)
    code("}");
  else
    code("endif");
}


void close_func()
{
  if (bC) 
    strcat(codebuffer,"}\n\n");
  else {
    code("return");
    code("end");
  }
}
