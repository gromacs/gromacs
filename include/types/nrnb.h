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
 * Gyas ROwers Mature At Cryogenic Speed
 */
enum { eNR_LJC, eNR_QQ, eNR_LJCRF, eNR_QQRF, eNR_FREE, eNR_COULTAB,
       eNR_TAB, eNR_WEIGHTS, eNR_SPREADQ, eNR_GATHERF, eNR_FFT, eNR_CONV,
       eNR_LR, eNR_BHAM, eNR_BHAMRF, eNR_NS, eNR_RESETX, eNR_SHIFTX,
       eNR_CGCM, eNR_FSUM,
       eNR_BONDS, eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER, eNR_RB,
       eNR_DISRES, eNR_POSRES, eNR_MORSE,
       eNR_VIRIAL, eNR_UPDATE, eNR_STOPCM, eNR_PCOUPL, eNR_EKIN,
       eNR_SHAKE, eNR_SHAKE_V, eNR_SHAKE_RIJ, eNR_SHAKE_VIR, 
       eNR_SETTLE, eNR_PSHAKEINITLD, eNR_PSHAKEINITMD, eNR_PSHAKE,
       eNR_DUM1, eNR_DUM2, eNR_DUM3, eNR_DUM4,
       eNRNB };

typedef struct {
  double n[eNRNB];
} t_nrnb;


