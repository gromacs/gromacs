/*
 * $Id$
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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _names_h
#define _names_h

static char *SRCID_names_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) names.h 1.23 5/20/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

/* All string arrays are NULL terminated, and therefore have an
 * extra argument (the +1)
 * these should correspond to names.c and include/types/enums.h
 */
extern char *eblock_names[ebNR+1];
extern char *epbc_names[epbcNR+1];
extern char *etcoupl_names[etcNR+1];
extern char *epcoupl_names[epcNR+1];
extern char *epcoupltype_names[epctNR+1];
extern char *ens_names[ensNR+1];
extern char *ei_names[eiNR+1];
extern char *yesno_names[BOOL_NR+1];
extern char *bool_names[BOOL_NR+1];
extern char *eel_names[eelNR+1];
extern char *eewg_names[eewgNR+1];
extern char *evdw_names[evdwNR+1];
extern char *eshake_names[estNR+1];
extern char *ptype_str[eptNR+1];
extern char *egrp_nm[egNR+1];
extern char *edisre_names[edrNR+1];
extern char *edisreweighting_names[edrwNR+1];
extern char *enbf_names[eNBF_NR+1];
extern char *ecomb_names[eCOMB_NR+1];
extern char *gtypes[egcNR+1];
extern char *efep_names[efepNR+1];
extern char *esolv_names[esolNR+1];
extern char *edispc_names[edispcNR+1];
extern char *ecm_names[ecmNR+1];

#define	UNDEFINED		"UNDEFINED"
#define ENUM_NAME(e,max,names)	((((e)<0)||((e)>=(max)))?UNDEFINED:(names)[e])

#define BOOL(e)        ENUM_NAME(e,BOOL_NR,bool_names)
#define ENS(e)         ENUM_NAME(e,ensNR,ens_names)
#define EI(e)          ENUM_NAME(e,eiNR,ei_names)
#define EPBC(e)        ENUM_NAME(e,epbcNR,epbc_names)
#define ETCOUPLTYPE(e) ENUM_NAME(e,etcNR,etcoupl_names)
#define EPCOUPLTYPE(e) ENUM_NAME(e,epcNR,epcoupl_names)
#define EPCOUPLTYPETYPE(e) ENUM_NAME(e,epctNR,epcoupltype_names)
#define EBLOCKS(e)     ENUM_NAME(e,ebNR,eblock_names)
#define EPARAM(e)      ENUM_NAME(e,epNR,eparam_names)
#define EELTYPE(e)     ENUM_NAME(e,eelNR,eel_names)
#define EVDWTYPE(e)    ENUM_NAME(e,evdwNR,evdw_names)
#define ESHAKETYPE(e)  ENUM_NAME(e,estNR,eshake_names)
#define EDISRETYPE(e)  ENUM_NAME(e,edrNR,edisre_names)
#define EDISREWEIGHTING(e)  ENUM_NAME(e,edrwNR,edisreweighting_names)
#define ENBFNAME(e)    ENUM_NAME(e,eNBF_NR,enbf_names)
#define ECOMBNAME(e)   ENUM_NAME(e,eCOMB_NR,ecomb_names)
#define EFEPTYPE(e)    ENUM_NAME(e,efepNR,efep_names)
#define ESOLVTYPE(e)   ENUM_NAME(e,esolNR,esolv_names)
#define EDISPCORR(e)   ENUM_NAME(e,edispcNR,edispc_names)
#define ECOM(e)        ENUM_NAME(e,ecmNR,ecm_names)
#endif	/* _names_h */
