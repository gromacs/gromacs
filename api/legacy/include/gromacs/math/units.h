/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_MATH_UNITS_H
#define GMX_MATH_UNITS_H

#include <cmath>

/*
 * Physical constants to be used in Gromacs.
 * No constants (apart from 0, 1 or 2) should
 * be anywhere else in the code.
 */

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#    define M_PI_2 1.57079632679489661923
#endif

#ifndef M_2PI
#    define M_2PI 6.28318530717958647692
#endif

#ifndef M_SQRT2
#    define M_SQRT2 std::sqrt(2.0)
#endif

#ifndef M_1_PI
#    define M_1_PI 0.31830988618379067154
#endif

#ifndef M_FLOAT_1_SQRTPI /* used in GPU kernels */
/* 1.0 / sqrt(M_PI) */
#    define M_FLOAT_1_SQRTPI 0.564189583547756f
#endif

#ifndef M_1_SQRTPI
/* 1.0 / sqrt(M_PI) */
#    define M_1_SQRTPI 0.564189583547756
#endif

#ifndef M_2_SQRTPI
/* 2.0 / sqrt(M_PI) */
#    define M_2_SQRTPI 1.128379167095513
#endif

namespace gmx
{

constexpr double c_angstrom       = 1e-10;
constexpr double c_kilo           = 1e3;
constexpr double c_nano           = 1e-9;
constexpr double c_pico           = 1e-12;
constexpr double c_nm2A           = c_nano / c_angstrom;
constexpr double c_cal2Joule      = 4.184;           /* Exact definition of the calorie */
constexpr double c_electronCharge = 1.602176634e-19; /* Exact definition, Coulomb NIST 2018 CODATA */

constexpr double c_amu       = 1.66053906660e-27; /* kg, NIST 2018 CODATA  */
constexpr double c_boltzmann = 1.380649e-23;      /* (J/K, Exact definition, NIST 2018 CODATA */
constexpr double c_avogadro  = 6.02214076e23;     /* 1/mol, Exact definition, NIST 2018 CODATA */
constexpr double c_universalGasConstant = c_boltzmann * c_avogadro;        /* (J/(mol K))  */
constexpr double c_boltz                = c_universalGasConstant / c_kilo; /* (kJ/(mol K)) */
constexpr double c_faraday              = c_electronCharge * c_avogadro;   /* (C/mol)      */
constexpr double c_planck1 = 6.62607015e-34; /* J/Hz, Exact definition, NIST 2018 CODATA */
constexpr double c_planck  = (c_planck1 * c_avogadro / (c_pico * c_kilo)); /* (kJ/mol) ps */

constexpr double c_epsilon0Si = 8.8541878128e-12; /* F/m,  NIST 2018 CODATA */
/* Epsilon in our MD units: (e^2 / Na (kJ nm)) == (e^2 mol/(kJ nm)) */
constexpr double c_epsilon0 =
        ((c_epsilon0Si * c_nano * c_kilo) / (c_electronCharge * c_electronCharge * c_avogadro));

constexpr double c_speedOfLight =
        2.99792458e05; /* units of nm/ps, Exact definition, NIST 2018 CODATA */

constexpr double c_rydberg = 1.0973731568160e-02; /* nm^-1, NIST 2018 CODATA */

constexpr double c_one4PiEps0 = (1.0 / (4.0 * M_PI * c_epsilon0));

/* Pressure in MD units is:
 * 1 bar = 1e5 Pa = 1e5 kg m^-1 s^-2 = 1e-28 kg nm^-1 ps^-2 = 1e-28 / amu amu nm^1 ps ^2
 */
constexpr double c_barMdunits = (1e5 * c_nano * c_pico * c_pico / c_amu);
constexpr double c_presfac    = 1.0 / c_barMdunits;

/* c_debye2Enm should be (1e-21*c_pico)/(c_speedOfLight*c_electronCharge*c_nano*c_nano),
 * but we need to factor out some of the exponents to avoid single-precision overflows.
 */
constexpr double c_debye2Enm = (1e-15 / (c_speedOfLight * c_electronCharge));
constexpr double c_enm2Debye = 1.0 / c_debye2Enm;

/* to convert from a acceleration in (e V)/(amu nm) */
/* c_fieldfac is also Faraday's constant and c_electronCharge/(1e6 amu) */
constexpr double c_fieldfac = c_faraday / c_kilo;

/* to convert AU to MD units: */
constexpr double c_hartree2Kj     = ((2.0 * c_rydberg * c_planck * c_speedOfLight) / c_avogadro);
constexpr double c_bohr2Nm        = 0.0529177210903; /* nm^-1, NIST 2018 CODATA */
constexpr double c_hartreeBohr2Md = (c_hartree2Kj * c_avogadro / c_bohr2Nm);

constexpr double c_rad2Deg = 180.0 / M_PI;
constexpr double c_deg2Rad = M_PI / 180.0;
} // namespace gmx

/* The four basic units */
#define unit_length "nm"
#define unit_time "ps"
#define unit_mass "u"
#define unit_energy "kJ/mol"

/* Temperature unit, T in this unit times c_boltz give energy in unit_energy */
#define unit_temp_K "K"

/* Charge unit, electron charge, involves c_one4PiEps0 */
#define unit_charge_e "e"

/* Pressure unit, pressure in basic units times c_presfac gives this unit */
#define unit_pres_bar "bar"

/* Dipole unit, debye, conversion from the unit_charge_e involves c_enm2Debye */
#define unit_dipole_D "D"

/* Derived units from basic units only */
#define unit_vel unit_length "/" unit_time
#define unit_volume unit_length "^3"
#define unit_invtime "1/" unit_time

/* Other derived units */
#define unit_surft_bar unit_pres_bar " " unit_length

/* SI units, conversion from basic units involves c_nano, c_pico and amu */
#define unit_length_SI "m"
#define unit_time_SI "s"
#define unit_mass_SI "kg"

#define unit_density_SI unit_mass_SI "/" unit_length_SI "^3"
#define unit_invvisc_SI unit_length_SI " " unit_time_SI "/" unit_mass_SI

#endif
