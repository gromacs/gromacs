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

/* Oh my god, it's full of loops!
 * There are quite a few innerloops, so they have been given numbers
 * instead of names. The first figure is the coulomb alternative, the
 * second vdw, the third the solvent opt and finally the fourth free
 * energy. 0 implies no, none or turned off. The other figugures mean:
 *                                     value
 * pos              1                   2           3             4
 * 1st Coul      Normal           Reaction-field  Table
 * 2nd Vdw       Lennard-Jones    Buckingham      Table        Bham-table
 * 3rd Sol       General solvent  Water           Water-Water
 * 4th FreeEner  Lambda           Softcore
 */

#define eNR_INLNONE -1

enum {
  eNR_INL0100, eNR_INL0110,
  eNR_INL0200, eNR_INL0210,
  eNR_INL0300, eNR_INL0301, eNR_INL0302, eNR_INL0310,
  eNR_INL0400, eNR_INL0401, eNR_INL0402, eNR_INL0410,
  eNR_INL1000, eNR_INL1010,
  eNR_INL1020, eNR_INL1030, eNR_INL1100, eNR_INL1110, eNR_INL1120,
  eNR_INL1130, eNR_INL1200, eNR_INL1210, eNR_INL1220, eNR_INL1230,
  eNR_INL1300, eNR_INL1310, eNR_INL1320, eNR_INL1330, eNR_INL1400,
  eNR_INL1410, eNR_INL1420, eNR_INL1430, eNR_INL2000, eNR_INL2010,
  eNR_INL2020, eNR_INL2030, eNR_INL2100, eNR_INL2110, eNR_INL2120,
  eNR_INL2130, eNR_INL2200, eNR_INL2210, eNR_INL2220, eNR_INL2230,
  eNR_INL2300, eNR_INL2310, eNR_INL2320, eNR_INL2330, eNR_INL2400,
  eNR_INL2410, eNR_INL2420, eNR_INL2430, eNR_INL3000, eNR_INL3001,
  eNR_INL3002, eNR_INL3010, eNR_INL3020, eNR_INL3030, eNR_INL3100,
  eNR_INL3110, eNR_INL3120, eNR_INL3130, eNR_INL3200, eNR_INL3210,
  eNR_INL3220, eNR_INL3230, eNR_INL3300, eNR_INL3301, eNR_INL3302,
  eNR_INL3310, eNR_INL3320, eNR_INL3330, eNR_INL3400, eNR_INL3401,
  eNR_INL3402, eNR_INL3410, eNR_INL3420, eNR_INL3430, eNR_INLOOP,       
  eNR_INL_IATOM=eNR_INLOOP,
  eNR_WEIGHTS,              eNR_SPREADQ,              eNR_SPREADQBSP,
  eNR_GATHERF,              eNR_GATHERFBSP,           eNR_FFT,
  eNR_CONV,                 eNR_SOLVEPME,eNR_NS,      eNR_RESETX,
  eNR_SHIFTX,               eNR_CGCM,                 eNR_FSUM,
  eNR_BONDS,                eNR_G96BONDS,             eNR_FENEBONDS,
  eNR_ANGLES,               eNR_G96ANGLES,            eNR_QANGLES,
  eNR_PROPER,               eNR_IMPROPER,
  eNR_RB,                   eNR_FOURDIH,              eNR_DISRES,               
  eNR_ORIRES,               eNR_DIHRES,
  eNR_POSRES,               eNR_ANGRES,               eNR_ANGRESZ,
  eNR_MORSE,                eNR_CUBICBONDS,
  eNR_WPOL,                 eNR_VIRIAL,
  eNR_UPDATE,               eNR_EXTUPDATE,            eNR_STOPCM,
  eNR_PCOUPL,               eNR_EKIN,                 eNR_LINCS,
  eNR_LINCSMAT,             eNR_SHAKE,                eNR_CONSTR_V,
  eNR_SHAKE_RIJ,            eNR_CONSTR_VIR,           eNR_SETTLE,
  eNR_VSITE2,               eNR_VSITE3,               eNR_VSITE3FD,
  eNR_VSITE3FAD,            eNR_VSITE3OUT,            eNR_VSITE4FD, 
  eNRNB
};


typedef struct {
  double n[eNRNB];
} t_nrnb;


