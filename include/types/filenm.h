/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Grunge ROck MAChoS
 */
 
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
  efHAT,
  efXPM,
  efNR
};

typedef struct {
  int  ftp;		/* File type (see enum above)		*/
  char *opt;		/* Command line option			*/
  char *fn;		/* File name				*/
  ulong flag;		/* Flag for all kinds of info (see defs)*/
} t_filenm;

#define ffSET    1<<0
#define ffREAD   1<<1
#define ffWRITE  1<<2
#define ffOPT    1<<3
#define ffLIB    1<<4
#define ffRW     (ffREAD  | ffWRITE)
#define ffOPTRD  (ffREAD  | ffOPT)
#define ffOPTWR  (ffWRITE | ffOPT)
#define ffOPTRW  (ffRW    | ffOPT)
#define ffLIBRD  (ffREAD  | ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
