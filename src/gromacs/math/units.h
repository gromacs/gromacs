/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_MATH_UNITS_H
#define GMX_MATH_UNITS_H

/*
 * Physical constants to be used in Gromacs.
 * No constants (apart from 0, 1 or 2) should
 * be anywhere else in the code.
 */

#include "gromacs/math/utilities.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ANGSTROM         (1e-10)                           /* Old...	*/
#define KILO             (1e3)                             /* Thousand	*/
#define NANO             (1e-9)                            /* A Number	*/
#define PICO             (1e-12)                           /* A Number	*/
#define A2NM             (ANGSTROM/NANO)                   /* NANO	        */
#define NM2A             (NANO/ANGSTROM)                   /* 10.0		*/
#define RAD2DEG          (180.0/M_PI)                      /* Conversion	*/
#define DEG2RAD          (M_PI/180.0)                      /* id		*/
#define CAL2JOULE        (4.184)                           /* id		*/
#define E_CHARGE         (1.602176565e-19)                 /* Coulomb, NIST 2010 CODATA */

#define AMU              (1.660538921e-27)                 /* kg, NIST 2010 CODATA  */
#define BOLTZMANN        (1.3806488e-23)                   /* (J/K, NIST 2010 CODATA */
#define AVOGADRO         (6.02214129e23)                   /* no unit, NIST 2010 CODATA */
#define RGAS             (BOLTZMANN*AVOGADRO)              /* (J/(mol K))  */
#define BOLTZ            (RGAS/KILO)                       /* (kJ/(mol K)) */
#define FARADAY          (E_CHARGE*AVOGADRO)               /* (C/mol)      */
#define ELECTRONVOLT     (E_CHARGE*AVOGADRO/KILO)          /* (kJ/mol)   */
#define PLANCK1          (6.62606957e-34)                  /* J s, NIST 2010 CODATA */
#define PLANCK           (PLANCK1*AVOGADRO/(PICO*KILO))    /* (kJ/mol) ps */

#define EPSILON0_SI      (8.854187817e-12)                 /* F/m,  NIST 2010 CODATA */
/* Epsilon in our MD units: (e^2 / Na (kJ nm)) == (e^2 mol/(kJ nm)) */
#define EPSILON0         (EPSILON0_SI*NANO*KILO)/(E_CHARGE*E_CHARGE*AVOGADRO)

#define SPEED_OF_LIGHT   (2.99792458E05)                   /* nm/ps, NIST 2010 CODATA */
#define ATOMICMASS_keV   (931494.061)                      /* Atomic mass in keV, NIST 2010 CODATA   */
#define ELECTRONMASS_keV (510.998928)                      /* Electron mas in keV, NIST 2010 CODATA  */

#define RYDBERG          (1.0973731568539e-02)             /* nm^-1, NIST 2010 CODATA */

#define ONE_4PI_EPS0     (1.0/(4.0*M_PI*EPSILON0))
#define FACEL            10.0*ONE_4PI_EPS0

/* Pressure in MD units is:
 * 1 bar = 1e5 Pa = 1e5 kg m^-1 s^-2 = 1e-28 kg nm^-1 ps^-2 = 1e-28 / AMU amu nm^1 ps ^2
 */
#define BAR_MDUNITS      (1e5*NANO*PICO*PICO/AMU)
#define PRESFAC          (1.0/BAR_MDUNITS)

/* DEBYE2ENM should be (1e-21*PICO)/(SPEED_OF_LIGHT*E_CHARGE*NANO*NANO),
 * but we need to factor out some of the exponents to avoid single-precision overflows.
 */
#define DEBYE2ENM        (1e-15/(SPEED_OF_LIGHT*E_CHARGE))
#define ENM2DEBYE        (1.0/DEBYE2ENM)

/* to convert from a acceleration in (e V)/(amu nm) */
/* FIELDFAC is also Faraday's constant and E_CHARGE/(1e6 AMU) */
#define FIELDFAC         (FARADAY/KILO)

/* to convert AU to MD units: */
#define HARTREE2KJ       ((2.0*RYDBERG*PLANCK*SPEED_OF_LIGHT)/AVOGADRO)
#define BOHR2NM          (0.052917721092)                  /* nm^-1, NIST 2010 CODATA */
#define HARTREE_BOHR2MD  (HARTREE2KJ*AVOGADRO/BOHR2NM)


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

/* The routines below can be used for converting units from or to GROMACS
   internal units. */
enum {
    eg2cAngstrom, eg2cNm, eg2cBohr, eg2cKcal_Mole,
    eg2cHartree, eg2cHartree_e, eg2cAngstrom3, eg2cCoulomb,
    eg2cDebye, eg2cElectron, eg2cBuckingham, eg2cNR
};

/* Convert value x to GROMACS units. Energy -> Energy, Length -> Length etc.
   The type of x is deduced from unit,
   which should be taken from the enum above. */
extern double convert2gmx(double x, int unit);

/* Convert value x from GROMACS units to the desired one.
   The type of return value is deduced from unit, see above */
extern double gmx2convert(double x, int unit);

/* Convert the string to one of the units supported. Returns -1 if not found. */
extern int string2unit(char *string);

/* Convert the unit to a string. Return NULL when unit is out of range. */
extern const char *unit2string(int unit);

#ifdef __cplusplus
}
#endif

#endif
