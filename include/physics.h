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
 * Gyas ROwers Mature At Cryogenic Speed
 */

#ifndef _physics_h
#define _physics_h

static char *SRCID_physics_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) physics.h 1.6 11/23/92"
#endif /* HAVE_IDENT */

/*
 * Physical constants to be used in Gromacs.
 * No constants (apart from 0, 1 or 2) should
 * be anywhere else in the code.
 */

#include <math.h>

#ifndef M_PI
#ifdef _PI
#define M_PI _PI
#else
#define M_PI        3.14159265358979323846
#endif
#endif

#define ANGSTROM	(1e-10)				/* Old...	*/
#define KILO 		(1e3)				/* Thousand	*/
#define NANO		(1e-9)				/* A Number	*/
#define A2NM		(ANGSTROM/NANO)			/* 0.1		*/
#define NM2A		(NANO/ANGSTROM)			/* 10.0		*/ 
#define RAD2DEG		(180.0/M_PI)			/* Conversion	*/
#define DEG2RAD		(M_PI/180.0)			/* id.		*/
#define CAL2JOULE	(4.184)				/* id.		*/
#define E_CHARGE        (1.60217733e-19)		/* Coulomb	*/

#define AMU             (1.6605402e-27)                 /* kg           */
#define BOLTZMANN	(1.380658e-23)			/* (J/K)	*/
#define AVOGADRO	(6.0221367e23)			/* ()		*/
#define RGAS            (BOLTZMANN*AVOGADRO)            /* (J/(mol K))  */
#define BOLTZ           (RGAS/KILO)	                /* (kJ/(mol K))	*/
#define FARADAY         (E_CHARGE*AVOGADRO)             /* (C/mol)      */
#define EPS0            
#define EPSILON0	(5.72765E-4)			/* (e^2 Na/(kJ nm)) */
/*#define ONE_4PI_EPS0	(1.0/(4.0*M_PI*EPSILON0))       	*/
#define FACEL		332.0636*CAL2JOULE	        /* (sqrt(ONE_4PI_EPS0)) */
#define ONE_4PI_EPS0	FACEL*0.1
#define PRESFAC         (16.6054)                       /* bar / pressure unity */

/* to convert from a acceleration in (e V)/(amu nm) */
/* FIELDFAC is also Faraday's constant and E_CHARGE/(1e6 AMU) */
#define FIELDFAC       (FARADAY/KILO)        
#endif	/* _physics_h */
