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

enum {
  ebCGS,ebMOLS,ebSBLOCKS,ebNR
};

enum {
  ebtRECT, ebtTRICLINIC, ebtNONE, ebtNR
};

enum {
  epcNO, epcISOTROPIC, epcANISOTROPIC, epcTRICLINIC, epcNR
};

enum {
  eelTWIN, eelRF, eelGRF, eelPME, eelPPPM, eelSWITCH, eelSHIFT, eelNR
};

enum { 
  ensGRID, ensSIMPLE, enNR
};

enum {
  eiMD, eiSteep, eiCG, eiLD, eiNR
};

enum {
  estLINCS, estSHAKE, estNR
};
