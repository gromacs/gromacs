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
#ifndef _types_nrnb_h
#define _types_nrnb_h


#ifdef __cplusplus
extern "C" {
#endif

/* The nonbonded kernels are documented in gmxlib/nonbonded_kernels, 
 * but here's a lazy version of the numbering. The first position
 * is the Coulomb interaction (0 for none), second is Van der Waals
 * (again, 0 means no interaction), and the third is the water optimization
 * (0 meaning no water optimization = standard atom-atom loop)
 *
 *                                     value
 * pos                 1                   2           3              4
 * 1st Coul        Normal,1/r       Reaction-field  Table            Generalized born
 * 2nd Vdw         Lennard-Jones    Buckingham      Table             n/a
 * 3rd Water. opt  SPC-other atom   SPC-SPC         TIP4p-other at.  TIP4p-TIP4p 
 */

#define eNR_NBKERNEL_NONE -1

enum 
{
    eNR_NBKERNEL010, eNR_NBKERNEL020, eNR_NBKERNEL030,
    eNR_NBKERNEL100, eNR_NBKERNEL101, eNR_NBKERNEL102, eNR_NBKERNEL103, eNR_NBKERNEL104,
    eNR_NBKERNEL110, eNR_NBKERNEL111, eNR_NBKERNEL112, eNR_NBKERNEL113, eNR_NBKERNEL114,
    eNR_NBKERNEL120, eNR_NBKERNEL121, eNR_NBKERNEL122, eNR_NBKERNEL123, eNR_NBKERNEL124,
    eNR_NBKERNEL130, eNR_NBKERNEL131, eNR_NBKERNEL132, eNR_NBKERNEL133, eNR_NBKERNEL134,
    eNR_NBKERNEL200, eNR_NBKERNEL201, eNR_NBKERNEL202, eNR_NBKERNEL203, eNR_NBKERNEL204,
    eNR_NBKERNEL210, eNR_NBKERNEL211, eNR_NBKERNEL212, eNR_NBKERNEL213, eNR_NBKERNEL214,
    eNR_NBKERNEL220, eNR_NBKERNEL221, eNR_NBKERNEL222, eNR_NBKERNEL223, eNR_NBKERNEL224,
    eNR_NBKERNEL230, eNR_NBKERNEL231, eNR_NBKERNEL232, eNR_NBKERNEL233, eNR_NBKERNEL234,
    eNR_NBKERNEL300, eNR_NBKERNEL301, eNR_NBKERNEL302, eNR_NBKERNEL303, eNR_NBKERNEL304,
    eNR_NBKERNEL310, eNR_NBKERNEL311, eNR_NBKERNEL312, eNR_NBKERNEL313, eNR_NBKERNEL314,
    eNR_NBKERNEL320, eNR_NBKERNEL321, eNR_NBKERNEL322, eNR_NBKERNEL323, eNR_NBKERNEL324,
    eNR_NBKERNEL330, eNR_NBKERNEL331, eNR_NBKERNEL332, eNR_NBKERNEL333, eNR_NBKERNEL334,
    eNR_NBKERNEL400, eNR_NBKERNEL410, eNR_NBKERNEL430,
    eNR_NBKERNEL010NF, eNR_NBKERNEL020NF, eNR_NBKERNEL030NF,
    eNR_NBKERNEL100NF, eNR_NBKERNEL101NF, eNR_NBKERNEL102NF, eNR_NBKERNEL103NF, eNR_NBKERNEL104NF,
    eNR_NBKERNEL110NF, eNR_NBKERNEL111NF, eNR_NBKERNEL112NF, eNR_NBKERNEL113NF, eNR_NBKERNEL114NF,
    eNR_NBKERNEL120NF, eNR_NBKERNEL121NF, eNR_NBKERNEL122NF, eNR_NBKERNEL123NF, eNR_NBKERNEL124NF,
    eNR_NBKERNEL130NF, eNR_NBKERNEL131NF, eNR_NBKERNEL132NF, eNR_NBKERNEL133NF, eNR_NBKERNEL134NF,
    eNR_NBKERNEL200NF, eNR_NBKERNEL201NF, eNR_NBKERNEL202NF, eNR_NBKERNEL203NF, eNR_NBKERNEL204NF,
    eNR_NBKERNEL210NF, eNR_NBKERNEL211NF, eNR_NBKERNEL212NF, eNR_NBKERNEL213NF, eNR_NBKERNEL214NF,
    eNR_NBKERNEL220NF, eNR_NBKERNEL221NF, eNR_NBKERNEL222NF, eNR_NBKERNEL223NF, eNR_NBKERNEL224NF,
    eNR_NBKERNEL230NF, eNR_NBKERNEL231NF, eNR_NBKERNEL232NF, eNR_NBKERNEL233NF, eNR_NBKERNEL234NF,
    eNR_NBKERNEL300NF, eNR_NBKERNEL301NF, eNR_NBKERNEL302NF, eNR_NBKERNEL303NF, eNR_NBKERNEL304NF,
    eNR_NBKERNEL310NF, eNR_NBKERNEL311NF, eNR_NBKERNEL312NF, eNR_NBKERNEL313NF, eNR_NBKERNEL314NF,
    eNR_NBKERNEL320NF, eNR_NBKERNEL321NF, eNR_NBKERNEL322NF, eNR_NBKERNEL323NF, eNR_NBKERNEL324NF,
    eNR_NBKERNEL330NF, eNR_NBKERNEL331NF, eNR_NBKERNEL332NF, eNR_NBKERNEL333NF, eNR_NBKERNEL334NF,
    eNR_NBKERNEL400NF, eNR_NBKERNEL410NF, eNR_NBKERNEL430NF, 
    eNR_NBKERNEL_NR,
    eNR_NBKERNEL_FREE_ENERGY = eNR_NBKERNEL_NR,
    eNR_NBKERNEL_ALLVSALL,
    eNR_NBKERNEL_ALLVSALLGB,
    eNR_NBKERNEL_OUTER,
    eNR_NB14,
    eNR_BORN_RADII_STILL,     eNR_BORN_RADII_HCT_OBC,
    eNR_BORN_CHAINRULE,
    eNR_BORN_AVA_RADII_STILL, eNR_BORN_AVA_RADII_HCT_OBC,
    eNR_BORN_AVA_CHAINRULE,
    eNR_WEIGHTS,              eNR_SPREADQ,              eNR_SPREADQBSP,
    eNR_GATHERF,              eNR_GATHERFBSP,           eNR_FFT,
    eNR_CONV,                 eNR_SOLVEPME,eNR_NS,      eNR_RESETX,
    eNR_SHIFTX,               eNR_CGCM,                 eNR_FSUM,
    eNR_BONDS,                eNR_G96BONDS,             eNR_FENEBONDS,
    eNR_TABBONDS,             eNR_RESTRBONDS,
    eNR_ANGLES,               eNR_G96ANGLES,            eNR_QANGLES,
    eNR_TABANGLES,            eNR_PROPER,               eNR_IMPROPER,
    eNR_RB,                   eNR_FOURDIH,              eNR_TABDIHS,
    eNR_DISRES,               eNR_ORIRES,               eNR_DIHRES,
    eNR_POSRES,               eNR_ANGRES,               eNR_ANGRESZ,
    eNR_MORSE,                eNR_CUBICBONDS,           eNR_WALLS,
    eNR_WPOL,                 eNR_THOLE,                eNR_VIRIAL,
    eNR_UPDATE,               eNR_EXTUPDATE,            eNR_STOPCM,
    eNR_PCOUPL,               eNR_EKIN,                 eNR_LINCS,
    eNR_LINCSMAT,             eNR_SHAKE,                eNR_CONSTR_V,
    eNR_SHAKE_RIJ,            eNR_CONSTR_VIR,           eNR_SETTLE,
    eNR_VSITE2,               eNR_VSITE3,               eNR_VSITE3FD,
    eNR_VSITE3FAD,            eNR_VSITE3OUT,            eNR_VSITE4FD,
    eNR_VSITE4FDN,            eNR_VSITEN,               eNR_GB,
	eNR_CMAP,
    eNRNB
};


typedef struct {
  double n[eNRNB];
} t_nrnb;


typedef struct gmx_wallcycle *gmx_wallcycle_t;

#ifdef __cplusplus
}
#endif

#endif
