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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _metacode_h
#define _metacode_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <types/simple.h>

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

#define FCON "     &" /* continuation line in f77 */
#define max(a,b) (((a) > (b)) ? (a) : (b))


/* Array referencing shortcut */
#define ARRAY(a,idx)          _array(#a,#idx)

typedef struct {
  char typename[150];
  char name[150];
  bool breferenced;
  bool bvector;
  bool bconst;
  char constval[50];
} decl_t;         /* Argument and variable buffer element */

extern int prec;  /* precision (4=single, 8=double) */
extern int IND;   /* current indentation */
extern char *codebuffer; /* buffer to which code is written */
extern char header[10000]; /* buffer for info and loop name */
extern FILE *output;  /* output file */
extern decl_t *decl_list; /* list with args and vars */
extern int ndecl; /* length of above list */
extern int nargs; /* first nargs are arguments, rest is vars */

extern bool bC;   

/* Concatenate a string to a buffer with plus sign between terms. */
void add_to_buffer(char *buffer,char *term);
  
/* initiate output buffers */
void init_metacode(void);

/* write a function to file and empty buffers */
void flush_buffers(void);

/* Return the correct indentation as a string */
char *indent(void);

/* Print a line of code to the output file */
void code(char *fmt, ...);

void newline(void);

/* Add a comment */
void comment(char *s);

/* Define a new floating-point variable */
void declare_real(char *name);
void declare_real_vector(char *name);
void declare_const_real_vector(char *name);

void declare_const_real(char *name,double val);
void declare_const_int(char *name,int val);

void declare_int(char *name);
void declare_int_vector(char *name);
void declare_const_int_vector(char *name);
void declare_real4(char *name);
void declare_int4(char *name);
void declare_int8(char *name);
void declare_intreal(char *name);
void declare_other(char *typename,char *name);

/* Cray vector pragma */
void vector_pragma(void);


char *_array(char *a,char *idx, ...);

void _p_state(char *left,char *right,char *symb);

void file_error(char *fn);

void assign(char *left, char *right, ...);

void increment(char *left,char *right, ...);

void decrement(char *left,char *right, ...);

void add(char *left,char *r1,char *r2, ...);

void subtract(char *left,char *r1,char *r2, ...);

void multiply(char *left,char *r1,char *r2, ...);

void closeit(void);

void usage(int argc,char *argv[]);

int count_lines(char *fn);

void edit_warning(char *fn);

void start_loop(char *lvar,char *from,char *to);

void start_stride_loop(char *lvar,char *from,char *to, char *stride);

void end_loop(void);

void start_if(char *cond);
void do_else(void);
void end_if(void);

void close_func(void);

#endif

