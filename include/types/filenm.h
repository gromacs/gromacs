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
  efMDP, efGBP, efGCP, efGIP, efGDP, efWDP, efGCT, efGPP,
  efTRX, efTRN, efTRR, efTRJ, efXTC, efG87, 
  efENX, efEDR, efENE,
  efSTX, efGRO, efPDB, efBRK, efENT, efXDR,
  efLOG, efXVG, efOUT,
  efNDX, 
  efTOP, efITP,
  efTPX, efTPR,  efTPA, efTPB,
  efTEX, efRTP, efATP, efHDB,
  efDAT, efDLG, 
  efGLD,
  efMAP, efEPS, efMAT, efM2P,
  efVDW,
  efMTX, efVEC,
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

#define ffSET    1
#define ffREAD   2
#define ffWRITE  4
#define ffRW     (ffREAD  | ffWRITE)
#define ffOPT    8
#define ffOPTRD  (ffREAD  | ffOPT)
#define ffOPTWR  (ffWRITE | ffOPT)
#define ffOPTRW  (ffRW    | ffOPT)
#define ffLIB    16
#define ffLIBRD  (ffREAD  | ffLIB)
#define ffLIBOPTRD (ffOPTRD | ffLIB)
