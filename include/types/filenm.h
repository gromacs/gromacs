/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
 
/* this enum should correspond to the array deffile in gmxlib/filenm.c */
enum {
  efMDP, efGCT,
  efTRX, efTRN, efTRR, efTRJ, efXTC, efG87, 
  efENX, efEDR, efENE,
  efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT,
  efLOG, efXVG, efOUT,
  efNDX, 
  efTOP, efITP,
  efTPX, efTPS, efTPR, efTPA, efTPB,
  efTEX, efRTP, efATP, efHDB,
  efDAT, efDLG, 
  efMAP, efEPS, efMAT, efM2P,
  efMTX,
  efEDI, efEDO, 
  efPPA, efPDO,
  efHAT,
  efXPM,
  efNR
};

typedef struct {
  int  ftp;		/* File type (see enum above)		*/
  char *opt;		/* Command line option			*/
  char *fn;		/* File name (as set in source code)	*/
  unsigned long flag;	/* Flag for all kinds of info (see defs)*/
  int  nf;		/* number of files			*/
  char **fns;		/* File names				*/
} t_filenm;

#define ffSET	1<<0
#define ffREAD	1<<1
#define ffWRITE	1<<2
#define ffOPT	1<<3
#define ffLIB	1<<4
#define ffMULT	1<<5
#define ffRW	(ffREAD	| ffWRITE)
#define ffOPTRD	(ffREAD	| ffOPT)
#define ffOPTWR	(ffWRITE| ffOPT)
#define ffOPTRW	(ffRW	| ffOPT)
#define ffLIBRD	(ffREAD	| ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
#define ffRDMULT   (ffREAD  | ffMULT)
