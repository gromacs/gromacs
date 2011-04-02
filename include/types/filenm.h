/*
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifdef __cplusplus
extern "C" {
#endif

 
/* this enum should correspond to the array deffile in gmxlib/filenm.c */
enum {
  efMDP, efGCT,
  efTRX, efTRO, efTRN, efTRR, efTRJ, efXTC, efG87, 
  efEDR,
  efSTX, efSTO, efGRO, efG96, efPDB, efBRK, efENT, efESP, efPQR, efXYZ,
  efCPT,
  efLOG, efXVG, efOUT,
  efNDX, 
  efTOP, efITP,
  efTPX, efTPS, efTPR, efTPA, efTPB, 
  efTEX, efRTP, efATP, efHDB,
  efDAT, efDLG, 
  efMAP, efEPS, efMAT, efM2P,
  efMTX,
  efEDI, efEDO, 
  efHAT,
  efCUB,
  efXPM,
  efRND,
  efNR
};

typedef struct {
  int  ftp;		/* File type (see enum above)		*/
  const char *opt;	/* Command line option			*/
  const char *fn;       /* File name (as set in source code)	*/
  unsigned long flag;	/* Flag for all kinds of info (see defs)*/
  int  nfiles;		/* number of files			*/
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
#define ffOPTRDMULT   (ffRDMULT | ffOPT)
#define ffWRMULT   (ffWRITE  | ffMULT)
#define ffOPTWRMULT   (ffWRMULT | ffOPT)

#ifdef __cplusplus
}
#endif

