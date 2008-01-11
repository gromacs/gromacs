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
 * GRoups of Organic Molecules in ACtion for Science
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* note: these enums should correspond to the names in gmxlib/names.c */

enum {
  ebCGS,ebMOLS,ebSBLOCKS,ebEXCLS,ebNR
};

enum {
  epbcXYZ, epbcNONE, epbcXY, epbcNR
};

enum {
  etcNO, etcBERENDSEN, etcNOSEHOOVER, etcYES, etcANDERSEN, etcANDERSENINTERVAL, etcNR
}; /* yes is an alias for berendsen */

enum {
  epcNO, epcBERENDSEN, epcPARRINELLORAHMAN, epcISOTROPIC, epcNR
}; /* isotropic is an alias for berendsen */

enum {
  epctISOTROPIC, epctSEMIISOTROPIC, epctANISOTROPIC,
  epctSURFACETENSION, epctNR
};

enum {
  erscNO, erscALL, erscCOM, erscNR
};

enum {
  eelCUT,     eelRF,     eelGRF,   eelPME,  eelEWALD,  eelPPPM, 
  eelPOISSON, eelSWITCH, eelSHIFT, eelUSER, eelGB, eelRF_NEC, eelENCADSHIFT, 
  eelPMEUSER, eelPMESWITCH, eelPMEUSERSWITCH, eelRF_ZERO, eelNR
};

/* Ewald geometry */
enum { 
  eewg3D, eewg3DC, eewgNR
};

#define EEL_RF(e) ((e) == eelRF || (e) == eelGRF || (e) == eelRF_NEC || (e) == eelRF_ZERO )

#define EEL_PME(e)  ((e) == eelPME || (e) == eelPMESWITCH || (e) == eelPMEUSER || (e) == eelPMEUSERSWITCH)
#define EEL_FULL(e) (EEL_PME(e) || (e) == eelPPPM || (e) == eelPOISSON || (e) == eelEWALD)

#define EEL_NOCUT(e) ((e) == eelSWITCH || (e) == eelSHIFT || (e) == eelENCADSHIFT || (e) == eelPMESWITCH || (e) == eelPMEUSERSWITCH || (e) == eelRF_ZERO)

enum {
  evdwCUT, evdwSWITCH, evdwSHIFT, evdwUSER, evdwENCADSHIFT, evdwNR
};

#define EVDW_NOCUT(e) ((e) == evdwSWITCH || (e) == evdwSHIFT || (e) == evdwENCADSHIFT)

enum { 
  ensGRID, ensSIMPLE, ensNR
};

enum {
  eiMD, eiSteep, eiCG, eiBD, eiSD2, eiNM, eiLBFGS, eiTPI, eiTPIC, eiSD, eiNR
};

#define EI_SD(e) ((e) == eiSD || (e) == eiSD2)
#define EI_DYNAMICS(e) ((e) == eiMD || EI_SD(e) || (e) == eiBD)
#define EI_ENERGY_MINIMIZATION(e) ((e) == eiSteep || (e) == eiCG || (e) == eiLBFGS)
#define EI_TPI(e) ((e) == eiTPI || (e) == eiTPIC)

#define EI_STATE_VELOCITY(e) ((e) == eiMD || EI_SD(e))

enum {
  estLINCS, estSHAKE, estNR
};

enum {
  edrNone, edrSimple, edrEnsemble, edrNR
};

enum {
  edrwConservative, edrwEqual, edrwNR
};

/* Combination rule things */
enum { 
  eCOMB_NONE, eCOMB_GEOMETRIC, eCOMB_ARITHMETIC, eCOMB_GEOM_SIG_EPS, eCOMB_NR 
};

/* NBF selection */
enum { 
  eNBF_NONE, eNBF_LJ, eNBF_BHAM, eNBF_NR 
};

/* FEP selection */
enum {
  efepNO, efepYES, efepNR
};

/* Solvent model */
enum {
  esolNO, esolSPC, esolTIP4P, esolNR
};

/* Neighborlist type */
enum {
  enlistATOM, enlistWATER, enlistWATERWATER, enlistNR
};

/* Dispersion correction */
enum {
  edispcNO, edispcEnerPres, edispcEner, edispcAllEnerPres, edispcAllEner, edispcNR
}; 

/* Shell types, for completion stuff */
enum {
  eshellCSH, eshellBASH, eshellZSH, eshellNR
}; 

/* Center of mass motion selection */
enum { 
  ecmLINEAR, ecmANGULAR, ecmNO, ecmNR 
};

/* New version of simulated annealing */
enum { 
  eannNO, eannSINGLE, eannPERIODIC, eannNR 
};

/* Algorithms for calculating GB radii */
enum { 
  egbSTILL, egbKARPLUS, egbNR 
};

/* Implicit solvent algorithms */
enum { 
  eisNO, eisLCPO, eisNR 
};

/* Wall types */
enum {
  ewt93, ewt104, ewtTABLE, ewtNR
};

/* Pull stuff */
enum {
  epullNO, epullUMBRELLA, epullCONSTRAINT, epullCONST_F, epullNR
};

enum {
  epullgDIST, epullgDIR, epullgCYL, epullgPOS, epullgNR
};

#define PULL_CYL(pull) ((pull)->eGeom == epullgCYL)

/* QMMM */
enum {
  eQMmethodAM1, eQMmethodPM3, eQMmethodRHF, 
  eQMmethodUHF, eQMmethodDFT, eQMmethodB3LYP, eQMmethodMP2, eQMmethodCASSCF, eQMmethodB3LYPLAN,
  eQMmethodDIRECT, eQMmethodNR
};

enum {
  eQMbasisSTO3G, eQMbasisSTO3G2, eQMbasis321G, 
  eQMbasis321Gp, eQMbasis321dGp, eQMbasis621G,
  eQMbasis631G, eQMbasis631Gp, eQMbasis631dGp, 
  eQMbasis6311G, eQMbasisNR
};

enum {
  eQMMMschemenormal,eQMMMschemeoniom,eQMMMschemeNR
};

