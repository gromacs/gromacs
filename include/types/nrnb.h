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
enum { eNR_LJC,      eNR_QQ,           eNR_LJCRF,       eNR_QQRF, 
       eNR_BHAM,     eNR_BHAMRF,       eNR_COULTAB,     eNR_TAB, 
       eNR_BHAMTAB,  eNR_BHAMTAB_WAT,
       eNR_LJC_WAT,  eNR_QQ_WAT,       eNR_LJCRF_WAT,   eNR_QQRF_WAT, 
       eNR_BHAM_WAT, eNR_BHAMRF_WAT,   eNR_COULTAB_WAT, eNR_TAB_WAT, 
       eNR_LJC_FREE, eNR_BHAM_FREE,    eNR_LJC_EW,      eNR_QQ_EW,        
       eNR_BHAM_EW,  eNR_LJC_WAT_EW,   eNR_QQ_WAT_EW,
       eNR_BHAM_WAT_EW,
       eNR_INLOOP,       
       eNR_INL_IATOM=eNR_INLOOP,
       eNR_WEIGHTS,  eNR_SPREADQ,      eNR_SPREADQBSP,  eNR_GATHERF,
       eNR_GATHERFBSP,                 eNR_FFT, 
       eNR_CONV,     eNR_SOLVEPME,
       eNR_NS,       eNR_RESETX,       eNR_SHIFTX,
       eNR_CGCM,     eNR_FSUM,
       eNR_BONDS,    eNR_ANGLES,       eNR_PROPER,      eNR_IMPROPER, eNR_RB,
       eNR_DISRES,   eNR_POSRES,       eNR_ANGRES,      eNR_ANGRESZ,
       eNR_MORSE,    eNR_WPOL,
       eNR_VIRIAL,   eNR_UPDATE,       eNR_STOPCM,      eNR_PCOUPL, eNR_EKIN,
       eNR_SHAKE,    eNR_SHAKE_V,      eNR_SHAKE_RIJ,   eNR_SHAKE_VIR, 
       eNR_SETTLE,   eNR_PSHAKEINITLD, eNR_PSHAKEINITMD, eNR_PSHAKE,
       eNR_DUM2,     eNR_DUM3,         eNR_DUM3FD,      eNR_DUM3FAD, 
       eNR_DUM3OUT,  eNR_DUM4FD, 
       eNRNB };

typedef struct {
  double n[eNRNB];
} t_nrnb;


