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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _physics_h
#define _physics_h

/*
 * Physical constants to be used in Gromacs.
 * No constants (apart from 0, 1 or 2) should
 * be anywhere else in the code.
 */

#include <math.h>

/* we do it anyway. */
#ifdef __cplusplus
extern "C" {
#endif

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
#define PLANCK1          (6.6262e-34)           /* J s */
#define PLANCK           (6.6262e-34*AVOGADRO/(PICO*KILO)) /* (kJ/mol) ps */

#define EPSILON0 	 (5.72765E-4)		/* (e^2 / Na (kJ nm))     
						   == (e^2 mol/(kJ nm)) */
                                                
#define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
#define ATOMICMASS_keV   (940000.0)             /* Atomic mass in keV   */
#define ELECTRONMASS_keV (512.0)                /* Electron mas in keV  */

/* Improved accuracy (PL & EL, 20090421) */
#define FACEL		 (332.0636930*CAL2JOULE)/* (10 * (ONE_4PI_EPS0)) */
#define ONE_4PI_EPS0	 (FACEL*0.1)            /* 1/(4*pi*e0)*/
#define PRESFAC           (16.6054)             /* bar / pressure unity */
#define ENM2DEBYE         48.0321               /* Convert electron nm  *
						 * to debye             */
#define DEBYE2ENM         0.02081941
/* to convert from a acceleration in (e V)/(amu nm) */
/* FIELDFAC is also Faraday's constant and E_CHARGE/(1e6 AMU) */
#define FIELDFAC          (FARADAY/KILO)

/* to convert AU to MD units: */
#define HARTREE2KJ        4.3597482e-21
#define BOHR2NM           0.0529177249
#define HARTREE_BOHR2MD   (HARTREE2KJ*AVOGADRO/BOHR2NM)


/* The four basic units */
#define unit_length   "nm"
#define unit_time     "ps"
#define unit_mass     "u"
#define unit_energy   "kJ/mol"

/* Temperature unit, T in this unit times BOLTZ give energy in unit_energy */
#define unit_temp_K   "K"

/* Charge unit, electron charge, involves ONE_4PI_EPS0 */
#define unit_charge_e "e"

/* Pressure unit, pressure in basic units times PRESFAC gives this unit */
#define unit_pres_bar "bar"

/* Dipole unit, debye, conversion from the unit_charge_e involves ENM2DEBYE */
#define unit_dipole_D "D"

/* Derived units from basic units only */
#define unit_vel      unit_length "/" unit_time
#define unit_volume   unit_length "^3"
#define unit_invtime  "1/" unit_time

/* Other derived units */
#define unit_surft_bar unit_pres_bar " " unit_length

/* SI units, conversion from basic units involves NANO, PICO and AMU */
#define unit_length_SI  "m"
#define unit_time_SI    "s"
#define unit_mass_SI    "kg"

#define unit_density_SI unit_mass_SI "/" unit_length_SI "^3"
#define unit_invvisc_SI unit_length_SI " " unit_time_SI "/" unit_mass_SI

#ifdef __cplusplus
}
#endif


#endif	/* _physics_h */


