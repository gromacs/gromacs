/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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

#ifndef _names_h
#define _names_h

static char *SRCID_names_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) names.h 1.23 5/20/97"
#endif /* HAVE_IDENT */

#include "typedefs.h"

/* All string arrays are NULL terminated, and therefore have an
 * extra argument (the +1)
 * these should correspond to names.c and include/types/enums.h
 */
extern char *eblock_names[ebNR+1];
extern char *eboxtype_names[ebtNR+1];
extern char *epcoupl_names[epcNR+1];
extern char *ens_names[enNR+1];
extern char *ei_names[eiNR+1];
extern char *yesno_names[BOOL_NR+1];
extern char *bool_names[BOOL_NR+1];
extern char *eel_names[eelNR+1];
extern char *eshake_names[estNR+1];
extern char *ptype_str[eptNR+1];
extern char *egrp_nm[egNR+1];
extern char *edisre_names[edrNR+1];
extern char *edisreweighting_names[edrwNR+1];
extern char *enbf_names[eNBF_NR+1];
extern char *ecomb_names[eCOMB_NR+1];

#define	UNDEFINED		"UNDEFINED"
#define ENUM_NAME(e,max,names)	((((e)<0)||((e)>=(max)))?UNDEFINED:(names)[e])

#define BOOL(e)        ENUM_NAME(e,BOOL_NR,bool_names)
#define ENS(e)         ENUM_NAME(e,enNR,ens_names)
#define EI(e)          ENUM_NAME(e,eiNR,ei_names)
#define EBOXTYPE(e)    ENUM_NAME(e,ebtNR,eboxtype_names)
#define EPCOUPLTYPE(e) ENUM_NAME(e,epcNR,epcoupl_names)
#define EBLOCKS(e)     ENUM_NAME(e,ebNR,eblock_names)
#define EPARAM(e)      ENUM_NAME(e,epNR,eparam_names)
#define EELTYPE(e)     ENUM_NAME(e,eelNR,eel_names)
#define ESHAKETYPE(e)  ENUM_NAME(e,estNR,eshake_names)
#define EDISRETYPE(e)  ENUM_NAME(e,edrNR,edisre_names)
#define EDISREWEIGHTING(e)  ENUM_NAME(e,edrwNR,edisreweighting_names)
#define ENBFNAME(e)    ENUM_NAME(e,eNBF_NR,enbf_names)
#define ECOMBNAME(e)   ENUM_NAME(e,eCOMB_NR,ecomb_names)

#endif	/* _names_h */
