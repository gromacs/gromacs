/*
 *
 * Gromacs 4.0                         Copyright (c) 1991-2003
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 *
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>

#include <mknb_metacode.h>

/* This file is NOT threadsafe, but it is only used to create
 * the nonbonded loops during the build process, so it will never be
 * executed by multiple threads.
 */


/* Global variables: */
int                   mknb_fortran=0;  /* 1 if fortran is used  */
int                   mknb_double=0;   /* 1 if double precision */
int                   mknb_keep_comments=0;
int                   mknb_indent_level=0;
FILE *                mknb_output=NULL;




/* Return the current indentation as a string */
static char *
mknb_indent(void)
{
  static char buf[255];
  int i;
  int indent;
  
  if(mknb_fortran)
	  indent = 6 + MKNB_FORTRAN_INDENT_STEP * mknb_indent_level;
  else
      indent = MKNB_C_INDENT_STEP * (mknb_indent_level+1);
	    
  if(indent > 254)
      indent = 254;
  	
  for(i=0; i<indent; i++)
      buf[i] = ' ';

  buf[i] = '\0';
  
  return buf;
}


/* Write the provided string as a comment in the file */
void 
mknb_comment(char *s)
{
  if(mknb_keep_comments) {
    if(mknb_fortran) {
      mknb_indent_level--;
      fprintf(mknb_output,"\nC%s%s\n",mknb_indent(),s);
      mknb_indent_level++;
    } else {
      fprintf(mknb_output,"\n%s/* %s */\n",mknb_indent(),s);
    }
  }
}



/* We need to separate argument and variable declarations.
 * The arguments are written directly to the output file,
 * while 
 */



/* Define a new floating-point variable */
void
mknb_declare_real(char *name)
{
  char type_name[255];


	
  if(mknb_fortran) 
  {
#ifdef IBM_FORTRAN_CPP
  	  sprintf(type_name, "%-13s", "gmxreal");
#else
	  sprintf(type_name, "%-13s", mknb_double ? "real*8" : "real*4");
#endif
  }
  else
    sprintf(type_name, "%-13s", mknb_double ? "double" : "float");

  mknb_code("%s %s%s",type_name,name, mknb_fortran ? "" : ";");
}

/* Declare a single-precision floating point var. */
void
mknb_declare_real4(char *name)
{
  char type_name[255];

  sprintf(type_name, "%-13s", mknb_fortran ? "real*4" : "float");

  mknb_code("%s %s%s", type_name ,name, mknb_fortran ? "" : ";");
}


/* Declare a constant fp variable */
void
mknb_declare_const_real(char *name, double value)
{
  char type_name[255];

  if(mknb_fortran) 
  {
#ifdef IBM_FORTRAN_CPP
	  sprintf(type_name, "%-13s", "gmxreal");
#else
	  sprintf(type_name, "%-13s", mknb_double ? "real*8" : "real*4");
#endif
	  mknb_code("%s %s",type_name,name);
      mknb_code("    parameter (%s = %f)",name,value);
  } else {
    sprintf(type_name, "%-13s",
	mknb_double ? "const double" : "const float");
    mknb_code("%s %s = %.16f;",type_name,name,value);
  }
}

/* Integer */
void
mknb_declare_int(char *name)
{
  char type_name[255];

  sprintf(type_name, "%-13s", mknb_fortran ? "integer*4" : "int"); 

  mknb_code("%s %s%s", type_name ,name, mknb_fortran ? "" : ";");
}


/* integer constant */
void
mknb_declare_const_int(char *name, int value)
{
  char type_name[255];

  sprintf(type_name, "%-13s", mknb_fortran ? "integer*4" : "const int");

  if(mknb_fortran) {
    mknb_code("%s %s", type_name ,name);
    mknb_code("    parameter (%s = %d)",name,value);
  } else
    mknb_code("%s %s = %d;", type_name ,name, value);
}

/* 4-byte Integer (same size as single precision fp) */
void
mknb_declare_int4(char *name)
{
  char type_name[255];

  sprintf(type_name, "%-13s", mknb_fortran ? "integer*4" : "int"); 

  mknb_code("%s %s%s", type_name ,name, mknb_fortran ? "" : ";");
}


/* Arbitrary declaration */
void
mknb_declare_other(char *type_name,char *name)
{
  char tmp[255];
  sprintf(tmp,"%-13s", type_name);

  mknb_code("%s %s%s", tmp, name, mknb_fortran ? "" : ";");
}


/* Reference an element in a list */
char *
mknb_array(char *a, char *idx)
{
  static char buf[1024];
  
  sprintf(buf,"%s%c%s%c",a, mknb_fortran ? '(' : '[',
	  idx, mknb_fortran ? ')' : ']');

  return buf;
}


/* Split a line to adhere to the stupid 72 column limits of fortran,
 * and write it to the output file
 */
static void 
mknb_fortran_splitline(char *line)
{
  char tmpbuf[4096];
  int i,j,maxlen;
  int indent;
  maxlen=strlen(line);

  i=0;
  
  indent = MKNB_FORTRAN_INDENT_STEP*mknb_indent_level;
  /* we can write 72-indentation characters on each line */
  
  while(i+72-indent<maxlen) {
    /* set j to the last position */
    j=i+71-indent;
    if (j>=maxlen)
      j=maxlen;

    while(j>(i+1)) {
      if(line[j]=='+' ||
         line[j]=='-' ||
         line[j]=='/' ||
         (line[j]=='*' && line[j-1]!='*') ||  /* dont split "**" */
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

  fprintf(mknb_output,"%s",tmpbuf);
  fprintf(mknb_output,"\n");
  fprintf(mknb_output,"     &  ");
  i=j;
  }
  fprintf(mknb_output,"%s",line+i);
}




/* Print a line of code to the output file.
 * This routine does proper indentation, and also supports the
 * same type of variable-argument lists as printf, apart from
 * field widths.
 */
void
mknb_code(char *format, ...)
{
  va_list ap;
  int d;
  double f;
  char *s;
  char tmp[4096];
  char outbuf[4096];

  sprintf(outbuf,"%s",mknb_indent());
 
  va_start(ap,format);
  
  while(*format) {
    if(*format!='%') 
      sprintf(tmp,"%c",*format);
    else {
      switch(*(++format)) {
      case 'd': /* read an integer */
        d = va_arg(ap, int);
        sprintf(tmp,"%d",d);
        break;
      case 'f': /* read a (double precision) floating point number */
        f = va_arg(ap, double);
        sprintf(tmp,"%.16f",f);
        break;
      case 's': /* read a string */
        s = va_arg(ap, char *);
        sprintf(tmp,"%s",s);
        break;
      default:
        fprintf(stderr,
		"Error, unsupported format supplied to mknb_code():\nn");
	
        exit(-1);
        break;
      }
    }
    format++;
    strcat(outbuf,tmp);
  }
  va_end(ap);

  if(mknb_fortran)
    mknb_fortran_splitline(outbuf);
  else
    fprintf(mknb_output,"%s",outbuf);

  fprintf(mknb_output,"\n");
}
  


/* Prints an assignment.
 * This routine does proper indentation, and also supports the
 * same type of variable-argument lists as printf (both in
 * the left and right-hand side buffers), apart from field widths.
 *
 * a statement like mknb_assign("atom%d","data%d",3,5) will give
 *
 * atom3 = data5;
 *
 *
 * In contrast to mknb_code(), mknb_assign() appends a semicolon when the
 * language is not set to fortran.
 * 
 */
void
mknb_assign(char *left, char *right, ...)
{
  int i;
  char *format;
  char buf[4096],tmp[4096];
  char outbuf[4096];
  va_list ap;
  int d;
  double f;
  char *s;

  
  sprintf(outbuf,"%s",mknb_indent());
 
  va_start(ap,right); 
 
  for(i=0;i<=1;i++) {
    /* first we do the left buffer, then repeat everything for the right. */
    if(i==0)
      format=left;
    else
      format=right;
    
    buf[0]='\0';

    while(*format) {
      if(*format!='%') 
	      sprintf(tmp,"%c",*format);
      else {
        switch(*(++format)) {
          case 'd': /* read an integer */
            d = va_arg(ap, int);
            sprintf(tmp,"%d",d);
            break;
          case 'f': /* read a (double precision) floating point number */
            f = va_arg(ap, double);
            sprintf(tmp,"%.16f",f);
            break;
          case 's': /* read a string */
            s = va_arg(ap, char *);
            sprintf(tmp,"%s",s);
            break;
          default:
            fprintf(stderr,
		    "Error, unsupported format supplied to mknb_assign()\n");
            exit(-1);
            break;
        }
      }
      strcat(buf,tmp);
      format++;
    }
    if(i==1 && !mknb_fortran)
      strcat(buf,";");
    
    sprintf(tmp,"%-16s",buf);
    strcat(outbuf,tmp);
    
    if(i==0)
      strcat(outbuf," = ");
  }
  va_end(ap);

  if(mknb_fortran)
    mknb_fortran_splitline(outbuf);
  else
    fprintf(mknb_output,"%s",outbuf);
  
  fprintf(mknb_output,"\n");
}
  


/* Start a for loop and increase indentation.*/
void
mknb_start_loop(char *lvar,char *from,char *to)
{

    mknb_code("");

    if(mknb_fortran)
    {
	    mknb_code("do %s=%s,%s",lvar,from,to);
    }
	else
    {
        mknb_code("for(%s=%s; (%s<%s); %s++)", lvar,from,lvar,to,lvar);
        mknb_code("{");
    }
    mknb_indent_level++;
}


/* decrease indentation and close for loop */
void
mknb_end_loop(void)
{
    mknb_indent_level--;

    if(mknb_fortran)
        mknb_code("end do");
    else
        mknb_code("}");
		
	mknb_code("");
}


void
mknb_start_if(char *cond)
{

    mknb_code("");
	
    if(mknb_fortran)
        mknb_code("if (%s) then",cond);
    else
	{
        mknb_code("if(%s)", cond);
        mknb_code("{");
    }
    mknb_indent_level++;
}

void
mknb_do_else()
{
    mknb_indent_level--;

    if(mknb_fortran)
        mknb_code("else");
    else
	{
        mknb_code("}");
        mknb_code("else");
        mknb_code("{");
    }
    mknb_indent_level++;
}



void
mknb_end_if(void)
{
    mknb_indent_level--;

    if (mknb_fortran)
        mknb_code("endif");
    else
        mknb_code("}");

    mknb_code("");

}


