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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _physics_h
#define _physics_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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

#define ANGSTROM	 (1e-10)		/* Old...	*/
#define KILO 		 (1e3)			/* Thousand	*/
#define NANO		 (1e-9)			/* A Number	*/
#define PICO		 (1e-12)		/* A Number	*/
#define A2NM		 (ANGSTROM/NANO)	/* NANO	        */
#define NM2A		 (NANO/ANGSTROM)	/* 10.0		*/
#define RAD2DEG		 (180.0/M_PI)		/* Conversion	*/
#define DEG2RAD		 (M_PI/180.0)		/* id		*/
#define CAL2JOULE	 (4.184)		/* id		*/
#define E_CHARGE         (1.60217733e-19)	/* Coulomb	*/

#define AMU              (1.6605402e-27)        /* kg           */
#define BOLTZMANN	 (1.380658e-23)		/* (J/K)	*/
#define AVOGADRO	 (6.0221367e23)		/* ()		*/
#define RGAS             (BOLTZMANN*AVOGADRO)   /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)            /* (kJ/(mol K)) */
#define FARADAY          (E_CHARGE*AVOGADRO)    /* (C/mol)      */
#define ELECTRONVOLT     (E_CHARGE*AVOGADRO/KILO) /* (kJ/mol)   */     
#define PLANCK           (6.6262e-34*AVOGADRO/(PICO*KILO)) /* kJ/mol ps */

#define EPSILON0 	 (5.72765E-4)		/* (e^2 Na/(kJ nm))     
						   == (e^2/(kJ mol nm)) */
                                                
#define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
#define ATOMICMASS_keV   (940000.0)             /* Atomic mass in keV   */
#define ELECTRONMASS_keV (512.0)                /* Electron mas in keV  */

#define FACEL		  332.0636*CAL2JOULE    /* (sqrt(ONE_4PI_EPS0)) */
#define ONE_4PI_EPS0	  FACEL*0.1
#define PRESFAC           (16.6054)             /* bar / pressure unity */
#define ENM2DEBYE         48.0321               /* Convert electron nm  *
						 * to debye             */
#define DEBYE2ENM         0.02081941
/* to convert from a acceleration in (e V)/(amu nm) */
/* FIELDFAC is also Faraday's constant and E_CHARGE/(1e6 AMU) */
#define FIELDFAC          (FARADAY/KILO)

#endif	/* _physics_h */


