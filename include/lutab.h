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

#ifndef _lutab_h
#define _lutab_h

#include <stdio.h>

static char *SRCID_lutab_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) lutab.h 1.7 2/2/97"
#endif /* HAVE_IDENT */
/*
   Fast sqrt(x) routine (actually 1.0/sqrt(x)) By A. Sijbers (c) 1992
   
   Based on Newton Rhapson iteration, initial value from lookup table

   Floating point representation:

         msb                                   lsb
         3322 2222 2222 1111 1111 1000 0000 0000
   bitnr 1098 7654 3210 9876 5432 1098 7654 3210
         |||| |||| |+-- ---- ---- ---- ---- ---+ fraction
         |+-- ---- + exponent
         +sign

   S = sign
   E = exponent
   F = fraction
                      S  E-127
   IEEE : value = (-1) (2     ) (1.F)     <------ must be float representation

                      S  E-128
   DEC  : value = (-1) (2     ) (0.1F)
*/

#define	EXP_LSB		0x00800000
#define	EXP_SEED_SIZE	256
#define	EXP_MASK	0x7f800000
#define	EXP_SHIFT	23
#define	MAX_FRACT	0x007fffff
#define	FRACT_MASK	0x007fffff
#define	FRACT_SIZE	11              /* significant part of fraction */
#define	FRACT_SHIFT	(EXP_SHIFT-FRACT_SIZE)
#define	FRACT_SEED_SIZE	(1<<(FRACT_SIZE+1))   /* one bit for even/odd */
#define	FRACT_FIRST	(0x3f800000>>FRACT_SHIFT)
#define	NOT_INITIALISED	~0
#define	EXP_ADDR(val)	(((val)&EXP_MASK)>>EXP_SHIFT)
#define	FRACT_ADDR(val)	(((val)&(FRACT_MASK|EXP_LSB))>>FRACT_SHIFT)

typedef unsigned int word;

typedef union 
{
  word bval;
  float fval;
} t_convert;

typedef struct
{
  word exp_seed[EXP_SEED_SIZE];
  word fract_seed[FRACT_SEED_SIZE];
} t_lutab;

/* Global variable, must be initiated by a call to init_lookup_table */
extern t_lutab lookup_table;

extern void init_lookup_table(FILE *log);

#endif	/* _lutab_h */
