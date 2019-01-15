/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief Data used by LINCS tests.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#include "gromacs/math/vectypes.h"

//! 0.1/sqrt(2), used to set coordinates
#define ONE_TENTH_OVER_SQRT_TWO  0.0707106781
//! 0.1/sqrt(3), used to set coordinates
#define TWO_TENTHS_OVER_SQRT_THREE 0.1154700538

namespace gmx
{

namespace test
{

//! A box with zero sides to check if PBC are actually disabled.
static const matrix c_infinitesimalBox = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

/*
 * System of two connected atoms (e.g. OH).
 */
//! System of two connected atoms: one constraint (type, i, j).
static const std::vector<int>  c_oneBondConstraints( {0, 0, 1} );
//! System of two connected atoms: equilibrium distance.
static const std::vector<real> c_oneBondConstraintsR0( {0.1} );
//! System of two connected atoms: starting type id for constraints.
static const int               c_oneBondConstraintsFirstType = 0;
//! System of two connected atoms: coordinates before 'timestep'.
static const std::vector<RVec> c_oneBondX(
        {{
             { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
             { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 }
         }} );
//! System of two connected atoms: coordinates after 'timestep', before constraints.
static const std::vector<RVec> c_oneBondXPrime(
        {{
             {  0.01,  0.08,  0.01 },
             {  0.06,  0.01, -0.01 }
         }} );
//! System of two connected atoms: masses.
static const std::vector<real> c_oneBondMasses ( {1.0, 12.0} );

/*
 * Two disjoint bonds (taken from the SHAKE test).
 */
//! Two disjoint bonds: two constraint (type1, i1, j1, type2, i2, j2).
static const std::vector<int>  c_twoDJBondsConstraints( {0, 0, 1, 1, 2, 3} );
//! Two disjoint bonds: equilibrium distances.
static const std::vector<real> c_twoDJBondsConstraintsR0( {2.0, 1.0} );
//! Two disjoint bonds: starting type id for constraints.
static const int               c_twoDJBondsConstraintsFirstType = 0;
//! Two disjoint bonds: coordinates.
static const std::vector<RVec> c_twoDJBondsX(
        {{
             {  2.50, -3.10, 15.70 },
             {  0.51, -3.02, 15.55 },
             { -0.50, -3.00, 15.20 },
             { -1.51, -2.95, 15.05 }
         }} );
//! Two disjoint bonds: masses.
static const std::vector<real> c_twoDJBondsMasses ( {0.5, 1.0/3.0, 0.25, 1.0} );


/*
 * Three atoms, connected longitudal (e.g. CH2).
 */
//! Three atoms, connected longitudaly: two constraints (type1, i1, j1, type2, i2, j2).
static const std::vector<int>  c_twoBondsConstraints( {0, 0, 1, 1, 1, 2} );
//! Three atoms, connected longitudaly: two distances.
static const std::vector<real> c_twoBondsConstraintsR0( {0.1, 0.2} );
//! Three atoms, connected longitudaly: starting type id for constraints.
static const int               c_twoBondsConstraintsFirstType = 0;
//! Three atoms, connected longitudaly: coordinates before timestep.
static const std::vector<RVec> c_twoBondsX(
        {{
             { ONE_TENTH_OVER_SQRT_TWO, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
             { 0.0, 0.0, 0.0 },
             { TWO_TENTHS_OVER_SQRT_THREE, TWO_TENTHS_OVER_SQRT_THREE, TWO_TENTHS_OVER_SQRT_THREE }
         }} );
//! Three atoms, connected longitudaly: coordinates after 'timestep', before constraints.
static const std::vector<RVec> c_twoBondsXPrime(
        {{
             {  0.08, 0.07,  0.01 },
             { -0.02, 0.01, -0.02 },
             {  0.10, 0.12,  0.11 }
         }} );
//! Three atoms, connected longitudaly: masses.
static const std::vector<real> c_twoBondsMasses ( {1.0, 12.0, 16.0 } );

/*
 * Three atoms, connected to the central atom (e.g. CH3).
 */
//! Three atoms, connected to the central atom: three constraints (type1, i1, j1, type2, i2, j2, type3, i3, j3).
static const std::vector<int>  c_threeBondsCentralConstraints( {0, 0, 1, 0, 0, 2, 0, 0, 3} );
//! Three atoms, connected to the central atom: two distances.
static const std::vector<real> c_threeBondsCentralConstraintsR0( {0.1} );
//! Three atoms, connected to the central atom: starting type id for constraints.
static const int               c_threeBondsCentralConstraintsFirstType = 0;
//! Three atoms, connected to the central atom: coordinates before timestep.
static const std::vector<RVec> c_threeBondsCentralX(
        {{
             { 0.00,  0.00,  0.00 },
             { 0.10,  0.00,  0.00 },
             { 0.00, -0.10,  0.00 },
             { 0.00,  0.00,  0.10 }
         }} );
//! Three atoms, connected to the central atom: coordinates after 'timestep', before constraints.
static const std::vector<RVec> c_threeBondsCentralXPrime(
        {{
             { 0.004,  0.009, -0.010 },
             { 0.110, -0.006,  0.003 },
             {-0.007, -0.102, -0.007 },
             {-0.005,  0.011,  0.102 }
         }} );
//! Three atoms, connected to the central atom: masses.
static const std::vector<real> c_threeBondsCentralMasses ( {12.0, 1.0, 1.0, 1.0 } );

/*
 * Four atoms, connected longitudaly (taken from SHAKE test).
 */
//! Four atoms, connected longitudaly: two constraints (type1, i1, j1, type2, i2, j2).
static const std::vector<int>  c_threeBondsConstraints( {0, 0, 1, 1, 1, 2, 2, 2, 3} );
//! Four atoms, connected longitudaly: two distances.
static const std::vector<real> c_threeBondsConstraintsR0( {2.0, 1.0, 1.0} );
//! Four atoms, connected longitudaly: starting type id for constraints.
static const int               c_threeBondsConstraintsFirstType = 0;
//! Four atoms, connected longitudaly: coordinates.
static const std::vector<RVec> c_threeBondsX(
        {{
             {  2.50, -3.10, 15.70 },
             {  0.51, -3.02, 15.55 },
             { -0.50, -3.00, 15.20 },
             { -1.51, -2.95, 15.05 }
         }} );
//! Four atoms, connected longitudaly: masses.
static const std::vector<real> c_threeBondsMasses ( {0.5, 1.0/3.0, 0.25, 1.0} );

/*
 * Basic triangle (tree atoms, connected with each other).
 */
//! Basic triangle: three constraints (type1, i1, j1, type2, i2, j2, type3, i3, j3).
static const std::vector<int>  c_triangleConstraints( {0, 0, 1, 2, 0, 2, 1, 1, 2} );
//! Basic triangle: euilibrium distances.
static const std::vector<real> c_triangleConstraintsR0( {0.1, 0.1, 0.1} );
//! Basic triangle: starting type id for constraints.
static const int               c_triangleConstraintsFirstType = 0;
//! Basic triangle: coordinates before 'timestep'.
static const std::vector<RVec> c_triangleX(
        {{
             { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 },
             { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
             { 0.0, 0.0, ONE_TENTH_OVER_SQRT_TWO }
         }} );
//! Basic triangle: coordinates after 'timestep', before constraining.
static const std::vector<RVec> c_triangleXPrime(
        {{
             {  0.09, -0.02,  0.01 },
             { -0.02,  0.10, -0.02 },
             {  0.03, -0.01,  0.07 }
         }} );
//! Basic triangle: masses.
static const std::vector<real> c_triangleMasses ( {1.0, 1.0, 1.0} );



/*
 * Real-life system: Cys-Val-Trp peptide.
 */
//! CVW peptide: periodic box.
static const matrix c_cvwBox = {{real(2.570950), 0, 0}, {0, real(2.570950), 0}, {0, 0, real(2.570950)}};


/*
 * Constraints only on covalent bonds with hydrogens.
 */
//! CVW peptide: constraints on bonds with hydrogens (type1, i1, j1, type2, i2, j2,...).
static const std::vector<int> c_cvwHBondsConstraints(
        {465,  0,  1, 465,  0,  2, 465,  0,  3, 466,  4,  5, 467,  6,  7, 467,  6,  8, 468,  9, 10, 469, 13, 14,
         466, 15, 16, 467, 17, 18, 467, 19, 20, 467, 19, 21, 467, 19, 22, 467, 23, 24, 467, 23, 25, 467, 23, 26,
         469, 29, 30, 466, 31, 32, 467, 33, 34, 467, 33, 35, 466, 37, 38, 470, 39, 40, 466, 43, 44, 466, 45, 46,
         466, 47, 48, 466, 49, 50} );
//! CVW peptide: equilibrium distances (one for each type of HBond).
static const std::vector<real> c_cvwHBondsConstraintsR0(
        {0.104000, 0.108000, 0.111100, 0.132500, 0.099700, 0.097600} );
//! CVW peptide: first type id used for constraints when HBonds are constrained.
static const int c_cvwHBondsConstraintsFirstType = 465;


/*
 *  Constraints on all covalent bonds.
 */
//! CVW peptide: constraints on all bonds (type1, i1, j1, type2, i2, j2,...).
static const std::vector<int> c_cvwAllBondsConstraints(
        {447,  0,  1, 447,  0,  2, 447,  0,  3, 448,  0,  4, 449,  4,  5, 450,  4,  6, 451,  4, 11, 452,  6,  7,
         452,  6,  8, 453,  6,  9, 454,  9, 10, 455, 11, 12, 456, 11, 13, 457, 13, 14, 458, 13, 15, 449, 15, 16,
         459, 15, 17, 451, 15, 27, 452, 17, 18, 450, 17, 19, 450, 17, 23, 452, 19, 20, 452, 19, 21, 452, 19, 22,
         452, 23, 24, 452, 23, 25, 452, 23, 26, 455, 27, 28, 456, 27, 29, 457, 29, 30, 458, 29, 31, 449, 31, 32,
         450, 31, 33, 460, 31, 51, 452, 33, 34, 452, 33, 35, 461, 33, 36, 462, 36, 37, 463, 36, 42, 449, 37, 38,
         464, 37, 39, 465, 39, 40, 466, 39, 41, 467, 41, 42, 468, 41, 47, 468, 42, 43, 449, 43, 44, 466, 43, 45,
         449, 45, 46, 466, 45, 49, 449, 47, 48, 466, 47, 49, 449, 49, 50, 469, 51, 52, 469, 51, 53} );
//! CVW peptide: equilibrium distances (one for each type of bond).
static const std::vector<real> c_cvwAllBondsConstraintsR0(
        {0.104000, 0.148000, 0.108000, 0.153800, 0.149000, 0.111100, 0.181800, 0.132500, 0.123000, 0.134500,
         0.099700, 0.143000, 0.150000, 0.152200, 0.151000, 0.136500, 0.144000, 0.137000, 0.097600, 0.137500,
         0.140000, 0.136800, 0.126000} );
//! CVW peptide: first type id used for constraints when all bonds are constrained.
static const int c_cvwAllBondsConstraintsFirstType = 447;


/*
 *  Constraints on all covalent bonds and all angles with hydrogens.
 */
//! CVW peptide: constraints on angles with hydrogens and all bonds (type1, i1, j1, type2, i2, j2,...).
static const std::vector<int> c_cvwHAnglesConstraints(
        {444,  1,  2, 444,  1,  3, 444,  2,  3, 445,  7,  8, 446, 20, 21, 446, 20, 22, 446, 21, 22, 446, 24, 25,
         446, 24, 26, 446, 25, 26, 445, 34, 35, 447,  0,  1, 447,  0,  2, 447,  0,  3, 448,  0,  4, 449,  4,  5,
         450,  4,  6, 451,  4, 11, 452,  6,  7, 452,  6,  8, 453,  6,  9, 454,  9, 10, 455, 11, 12, 456, 11, 13,
         457, 13, 14, 458, 13, 15, 449, 15, 16, 459, 15, 17, 451, 15, 27, 452, 17, 18, 450, 17, 19, 450, 17, 23,
         452, 19, 20, 452, 19, 21, 452, 19, 22, 452, 23, 24, 452, 23, 25, 452, 23, 26, 455, 27, 28, 456, 27, 29,
         457, 29, 30, 458, 29, 31, 449, 31, 32, 450, 31, 33, 460, 31, 51, 452, 33, 34, 452, 33, 35, 461, 33, 36,
         462, 36, 37, 463, 36, 42, 449, 37, 38, 464, 37, 39, 465, 39, 40, 466, 39, 41, 467, 41, 42, 468, 41, 47,
         468, 42, 43, 449, 43, 44, 466, 43, 45, 449, 45, 46, 466, 45, 49, 449, 47, 48, 466, 47, 49, 449, 49, 50,
         469, 51, 52, 469, 51, 53} );
//! CVW peptide: equilibrium distances (one for each type of constraint).
static const std::vector<real> c_cvwHAnglesConstraintsR0(
        {0.169861, 0.180896, 0.180218, 0.104000, 0.148000, 0.108000, 0.153800, 0.149000, 0.111100, 0.181800,
         0.132500, 0.123000, 0.134500, 0.099700, 0.143000, 0.150000, 0.152200, 0.151000, 0.136500, 0.144000,
         0.137000, 0.097600, 0.137500, 0.140000, 0.136800, 0.126000} );
//! CVW peptide: first type id used.
static const int c_cvwHAnglesConstraintsFirstType = 444;


//! Coordinates of all 54 atoms of the Cys-Val-Trp system before integration step
static const std::vector<RVec> c_cvwX(
        {{
             { 1.767684, 1.434249, 1.065166 },
             { 1.815431, 1.474586, 0.982044 },
             { 1.723099, 1.510077, 1.120647 },
             { 1.843457, 1.398358, 1.126699 },
             { 1.680907, 1.314366, 1.024421 },
             { 1.599158, 1.351943, 0.964679 },
             { 1.763957, 1.225508, 0.927973 },
             { 1.697651, 1.153768, 0.875057 },
             { 1.791425, 1.291521, 0.842938 },
             { 1.912151, 1.142162, 0.994032 },
             { 2.004733, 1.230230, 0.958979 },
             { 1.614016, 1.228077, 1.132146 },
             { 1.646264, 1.110038, 1.143347 },
             { 1.524281, 1.284690, 1.211713 },
             { 1.492912, 1.379322, 1.210727 },
             { 1.454358, 1.207665, 1.318520 },
             { 1.458889, 1.100976, 1.302365 },
             { 1.511439, 1.240173, 1.457196 },
             { 1.502893, 1.346852, 1.487025 },
             { 1.453182, 1.146234, 1.561338 },
             { 1.501521, 1.159239, 1.660522 },
             { 1.348249, 1.180627, 1.573558 },
             { 1.472060, 1.042550, 1.526174 },
             { 1.658998, 1.219855, 1.459685 },
             { 1.708764, 1.291962, 1.391369 },
             { 1.713251, 1.235827, 1.555313 },
             { 1.684297, 1.117161, 1.425668 },
             { 1.308254, 1.239545, 1.313299 },
             { 1.271975, 1.350856, 1.275704 },
             { 1.222009, 1.141459, 1.353709 },
             { 1.249863, 1.048388, 1.376113 },
             { 1.097332, 1.173901, 1.411180 },
             { 1.102543, 1.269817, 1.460547 },
             { 0.975575, 1.167806, 1.301958 },
             { 1.006804, 1.100189, 1.219521 },
             { 0.890948, 1.104405, 1.336041 },
             { 0.922556, 1.298571, 1.257180 },
             { 0.947959, 1.376065, 1.145647 },
             { 1.006637, 1.350066, 1.058786 },
             { 0.873536, 1.491546, 1.146863 },
             { 0.897865, 1.574278, 1.101155 },
             { 0.787589, 1.492447, 1.255357 },
             { 0.822571, 1.376258, 1.333269 },
             { 0.763600, 1.350945, 1.451014 },
             { 0.794500, 1.269275, 1.514568 },
             { 0.668060, 1.445360, 1.503948 },
             { 0.634763, 1.434579, 1.606120 },
             { 0.697418, 1.578386, 1.298386 },
             { 0.671498, 1.660423, 1.233101 },
             { 0.636332, 1.554670, 1.426202 },
             { 0.561673, 1.626039, 1.457772 },
             { 1.065444, 1.065935, 1.521060 },
             { 1.129701, 0.961464, 1.524790 },
             { 0.967337, 1.080518, 1.598781 }
         }});

//! Coordinates of all 54 atoms of the Cys-Val-Trp system after integration step but before constraining
static const std::vector<RVec> c_cvwXPrime(
        {{
             { 1.768339, 1.434433, 1.065362 },
             { 1.814669, 1.474604, 0.980970 },
             { 1.722634, 1.507937, 1.122880 },
             { 1.844173, 1.395584, 1.125705 },
             { 1.681784, 1.313778, 1.023859 },
             { 1.599284, 1.350122, 0.963407 },
             { 1.763891, 1.225048, 0.926888 },
             { 1.699154, 1.152251, 0.874503 },
             { 1.788218, 1.292084, 0.840873 },
             { 1.912297, 1.142312, 0.994077 },
             { 2.005435, 1.230518, 0.956869 },
             { 1.614129, 1.228172, 1.131230 },
             { 1.646479, 1.110057, 1.143719 },
             { 1.524452, 1.285291, 1.211535 },
             { 1.490126, 1.379115, 1.210074 },
             { 1.453782, 1.207371, 1.318983 },
             { 1.456185, 1.101159, 1.303672 },
             { 1.511876, 1.240363, 1.456429 },
             { 1.502018, 1.345875, 1.487895 },
             { 1.453369, 1.146519, 1.561700 },
             { 1.501137, 1.158523, 1.660962 },
             { 1.348219, 1.180476, 1.573327 },
             { 1.471519, 1.041876, 1.527539 },
             { 1.659148, 1.219680, 1.459319 },
             { 1.707492, 1.290872, 1.388623 },
             { 1.712837, 1.235675, 1.554854 },
             { 1.683602, 1.117280, 1.424265 },
             { 1.308488, 1.239062, 1.312624 },
             { 1.272906, 1.350718, 1.275423 },
             { 1.222072, 1.140970, 1.353401 },
             { 1.250822, 1.047869, 1.376124 },
             { 1.096433, 1.173753, 1.411309 },
             { 1.103983, 1.270662, 1.461157 },
             { 0.975261, 1.167491, 1.301993 },
             { 1.008031, 1.099915, 1.219796 },
             { 0.889935, 1.104291, 1.334737 },
             { 0.923996, 1.298329, 1.257339 },
             { 0.947920, 1.376093, 1.145296 },
             { 1.006564, 1.350725, 1.057899 },
             { 0.873005, 1.492125, 1.147714 },
             { 0.897790, 1.574741, 1.100132 },
             { 0.787382, 1.492639, 1.255753 },
             { 0.822208, 1.376489, 1.333382 },
             { 0.764090, 1.350932, 1.451206 },
             { 0.793943, 1.269517, 1.514956 },
             { 0.668235, 1.445627, 1.504846 },
             { 0.637204, 1.433128, 1.608086 },
             { 0.696730, 1.578846, 1.298825 },
             { 0.673728, 1.662102, 1.232573 },
             { 0.636632, 1.555425, 1.426305 },
             { 0.561369, 1.625075, 1.457715 },
             { 1.065776, 1.065940, 1.521084 },
             { 1.129488, 0.961840, 1.524568 },
             { 0.967672, 1.081063, 1.599542 }
         }});

//! Coordinates of 54 atoms of the Cys-Val-Trp system after H-Bonds were constrained with LINCS
static const std::vector<RVec> c_cvwXConstrainedHBondsLincs(
        {{
             { 1.768371, 1.434429, 1.065359 },
             { 1.814530, 1.474487, 0.981212 },
             { 1.722609, 1.507980, 1.122911 },
             { 1.843891, 1.395718, 1.125476 },
             { 1.681752, 1.313792, 1.023836 },
             { 1.599664, 1.349947, 0.963684 },
             { 1.763925, 1.225101, 0.926869 },
             { 1.698888, 1.151964, 0.874291 },
             { 1.788075, 1.291739, 0.841316 },
             { 1.912320, 1.142334, 0.994068 },
             { 2.004714, 1.229833, 0.957142 },
             { 1.614129, 1.228172, 1.131230 },
             { 1.646479, 1.110057, 1.143719 },
             { 1.524448, 1.285305, 1.211535 },
             { 1.490190, 1.378923, 1.210076 },
             { 1.453780, 1.207421, 1.318991 },
             { 1.456211, 1.100554, 1.303580 },
             { 1.511879, 1.240322, 1.456417 },
             { 1.501979, 1.346368, 1.488032 },
             { 1.453364, 1.146484, 1.561669 },
             { 1.501248, 1.158553, 1.661191 },
             { 1.348235, 1.180471, 1.573325 },
             { 1.471447, 1.042270, 1.527673 },
             { 1.659141, 1.219699, 1.459286 },
             { 1.707384, 1.290716, 1.388771 },
             { 1.712994, 1.235722, 1.555132 },
             { 1.683632, 1.117160, 1.424225 },
             { 1.308488, 1.239062, 1.312624 },
             { 1.272906, 1.350718, 1.275423 },
             { 1.222079, 1.140948, 1.353406 },
             { 1.250730, 1.048177, 1.376050 },
             { 1.096437, 1.173838, 1.411353 },
             { 1.103928, 1.269647, 1.460634 },
             { 0.975265, 1.167479, 1.301979 },
             { 1.007969, 1.100050, 1.219961 },
             { 0.889948, 1.104301, 1.334732 },
             { 0.923996, 1.298329, 1.257339 },
             { 0.947931, 1.376089, 1.145280 },
             { 1.006432, 1.350784, 1.058095 },
             { 0.873020, 1.492176, 1.147686 },
             { 0.897579, 1.574023, 1.100529 },
             { 0.787382, 1.492639, 1.255753 },
             { 0.822208, 1.376489, 1.333382 },
             { 0.764082, 1.350953, 1.451189 },
             { 0.794041, 1.269257, 1.515159 },
             { 0.668222, 1.445623, 1.504884 },
             { 0.637354, 1.433176, 1.607628 },
             { 0.696714, 1.578896, 1.298785 },
             { 0.673918, 1.661501, 1.233051 },
             { 0.636672, 1.555387, 1.426288 },
             { 0.560890, 1.625534, 1.457917 },
             { 1.065776, 1.065940, 1.521084 },
             { 1.129488, 0.961840, 1.524568 },
             { 0.967672, 1.081063, 1.599542 }
         }});

//! Virial tensor of the Cys-Val-Trp system after H-Bonds were constrained with LINCS
static const tensor c_cvwVirialHBondsLincs = {
    {  9.674206e-05,  7.297217e-05, -7.171490e-05 },
    {  7.297219e-05,  2.097286e-04, -1.490687e-04 },
    { -7.171491e-05, -1.490686e-04,  1.801554e-04 }
};

//! Coordinates of 54 atoms of the Cys-Val-Trp system after H-Bonds were constrained with SHAKE
static const std::vector<RVec> c_cvwXConstrainedHBondsShake(
        {{
             { 1.768383, 1.434433, 1.065362 },
             { 1.814521, 1.474481, 0.981190 },
             { 1.722609, 1.507947, 1.122962 },
             { 1.843899, 1.395669, 1.125452 },
             { 1.681769, 1.313781, 1.023825 },
             { 1.599672, 1.349907, 0.963664 },
             { 1.763925, 1.225093, 0.926850 },
             { 1.698908, 1.151929, 0.874265 },
             { 1.788006, 1.291744, 0.841285 },
             { 1.912323, 1.142337, 0.994069 },
             { 2.004718, 1.229828, 0.957108 },
             { 1.614129, 1.228173, 1.131214 },
             { 1.646483, 1.110058, 1.143726 },
             { 1.524453, 1.285318, 1.211529 },
             { 1.490131, 1.378913, 1.210071 },
             { 1.453768, 1.207416, 1.319001 },
             { 1.456153, 1.100543, 1.303604 },
             { 1.511887, 1.240324, 1.456402 },
             { 1.501976, 1.346359, 1.488055 },
             { 1.453367, 1.146488, 1.561675 },
             { 1.501241, 1.158539, 1.661205 },
             { 1.348239, 1.180483, 1.573317 },
             { 1.471447, 1.042266, 1.527696 },
             { 1.659142, 1.219697, 1.459279 },
             { 1.707355, 1.290687, 1.388731 },
             { 1.712999, 1.235719, 1.555120 },
             { 1.683619, 1.117162, 1.424200 },
             { 1.308492, 1.239053, 1.312611 },
             { 1.272924, 1.350716, 1.275417 },
             { 1.222079, 1.140937, 1.353401 },
             { 1.250744, 1.048171, 1.376051 },
             { 1.096425, 1.173837, 1.411357 },
             { 1.103956, 1.269650, 1.460639 },
             { 0.975257, 1.167475, 1.301978 },
             { 1.007982, 1.100038, 1.219975 },
             { 0.889930, 1.104291, 1.334694 },
             { 0.924025, 1.298321, 1.257343 },
             { 0.947931, 1.376089, 1.145273 },
             { 1.006423, 1.350805, 1.058076 },
             { 0.873012, 1.492187, 1.147702 },
             { 0.897574, 1.574021, 1.100517 },
             { 0.787374, 1.492648, 1.255760 },
             { 0.822200, 1.376494, 1.333387 },
             { 0.764092, 1.350953, 1.451189 },
             { 0.794033, 1.269258, 1.515172 },
             { 0.668223, 1.445629, 1.504903 },
             { 0.637408, 1.433150, 1.607660 },
             { 0.696703, 1.578904, 1.298791 },
             { 0.673967, 1.661517, 1.233042 },
             { 0.636679, 1.555399, 1.426292 },
             { 0.560873, 1.625526, 1.457918 },
             { 1.065784, 1.065939, 1.521085 },
             { 1.129482, 0.961848, 1.524564 },
             { 0.967678, 1.081072, 1.599556 }
         }});

//! Virial tensor of the Cys-Val-Trp system after H-Bonds were constrained with SHAKE
static const tensor c_cvwVirialHBondsShake = {
    {  9.575119e-05,  7.292318e-05, -7.163925e-05 },
    {  7.292317e-05,  2.072188e-04, -1.496644e-04 },
    { -7.163925e-05, -1.496644e-04,  1.790191e-04 }
};

//! Coordinates of 54 atoms of the Cys-Val-Trp system after all bonds were constrained with LINCS
static const std::vector<RVec> c_cvwXConstrainedAllBondsLincs(
        {{
             { 1.765786, 1.431034, 1.064172 },
             { 1.813698, 1.473691, 0.982321 },
             { 1.723303, 1.506400, 1.121876 },
             { 1.842830, 1.396068, 1.124647 },
             { 1.683956, 1.313432, 1.027175 },
             { 1.602427, 1.348627, 0.965742 },
             { 1.763739, 1.225913, 0.929087 },
             { 1.699695, 1.152954, 0.875088 },
             { 1.787759, 1.291080, 0.842383 },
             { 1.911792, 1.142630, 0.993833 },
             { 2.004630, 1.229783, 0.957208 },
             { 1.614965, 1.229995, 1.129455 },
             { 1.645988, 1.111852, 1.143557 },
             { 1.523809, 1.281873, 1.213567 },
             { 1.490984, 1.375923, 1.210248 },
             { 1.454970, 1.211261, 1.317055 },
             { 1.456086, 1.104072, 1.304092 },
             { 1.508965, 1.240387, 1.453899 },
             { 1.501968, 1.345938, 1.487809 },
             { 1.453157, 1.146206, 1.561909 },
             { 1.501343, 1.158630, 1.661237 },
             { 1.348093, 1.180439, 1.573415 },
             { 1.471460, 1.042133, 1.527614 },
             { 1.661224, 1.219419, 1.459254 },
             { 1.707876, 1.291195, 1.388439 },
             { 1.713441, 1.235850, 1.555928 },
             { 1.683914, 1.116512, 1.424065 },
             { 1.308460, 1.237802, 1.313074 },
             { 1.273262, 1.349578, 1.275825 },
             { 1.221916, 1.143209, 1.353544 },
             { 1.250210, 1.050204, 1.375553 },
             { 1.092993, 1.170530, 1.408951 },
             { 1.103608, 1.265809, 1.458622 },
             { 0.979231, 1.169195, 1.305543 },
             { 1.007326, 1.101671, 1.221951 },
             { 0.892305, 1.106067, 1.333756 },
             { 0.921962, 1.300642, 1.258356 },
             { 0.946803, 1.375744, 1.147181 },
             { 1.005333, 1.351209, 1.059828 },
             { 0.873173, 1.491251, 1.148273 },
             { 0.897314, 1.573153, 1.101035 },
             { 0.788880, 1.489787, 1.256843 },
             { 0.823490, 1.376323, 1.331072 },
             { 0.761792, 1.354559, 1.451155 },
             { 0.792950, 1.271888, 1.513212 },
             { 0.670891, 1.444455, 1.501661 },
             { 0.638591, 1.433470, 1.604113 },
             { 0.695044, 1.578519, 1.301853 },
             { 0.674295, 1.660130, 1.234255 },
             { 0.638738, 1.553832, 1.424809 },
             { 0.562849, 1.623524, 1.457135 },
             { 1.065249, 1.069404, 1.519204 },
             { 1.129809, 0.961362, 1.524578 },
             { 0.968409, 1.080956, 1.598943 }
         }});

//! Coordinates of 54 atoms of the Cys-Val-Trp system after all bonds were constrained with SHAKE
static const std::vector<RVec> c_cvwXConstrainedAllBondsShake(
        {{
             { 1.765800, 1.431069, 1.064180 },
             { 1.813715, 1.473702, 0.982321 },
             { 1.723292, 1.506446, 1.121845 },
             { 1.842831, 1.396108, 1.124675 },
             { 1.683946, 1.313470, 1.027131 },
             { 1.602375, 1.348679, 0.965727 },
             { 1.763749, 1.225906, 0.929080 },
             { 1.699659, 1.152964, 0.875096 },
             { 1.787824, 1.291090, 0.842396 },
             { 1.911793, 1.142625, 0.993834 },
             { 2.004623, 1.229785, 0.957239 },
             { 1.614972, 1.229976, 1.129452 },
             { 1.645997, 1.111797, 1.143560 },
             { 1.523787, 1.281925, 1.213569 },
             { 1.491005, 1.376022, 1.210245 },
             { 1.454995, 1.211206, 1.317074 },
             { 1.456135, 1.104001, 1.304064 },
             { 1.508970, 1.240407, 1.453943 },
             { 1.501969, 1.345989, 1.487805 },
             { 1.453145, 1.146188, 1.561920 },
             { 1.501355, 1.158645, 1.661234 },
             { 1.348089, 1.180429, 1.573422 },
             { 1.471462, 1.042128, 1.527590 },
             { 1.661238, 1.219418, 1.459259 },
             { 1.707900, 1.291218, 1.388478 },
             { 1.713438, 1.235854, 1.555946 },
             { 1.683926, 1.116510, 1.424087 },
             { 1.308458, 1.237830, 1.313073 },
             { 1.273232, 1.349627, 1.275814 },
             { 1.221928, 1.143154, 1.353546 },
             { 1.250223, 1.050125, 1.375571 },
             { 1.093001, 1.170574, 1.408984 },
             { 1.103598, 1.265879, 1.458658 },
             { 0.979225, 1.169112, 1.305528 },
             { 1.007325, 1.101636, 1.221872 },
             { 0.892256, 1.106035, 1.333799 },
             { 0.921979, 1.300614, 1.258308 },
             { 0.946837, 1.375763, 1.147116 },
             { 1.005394, 1.351168, 1.059766 },
             { 0.873197, 1.491268, 1.148215 },
             { 0.897335, 1.573213, 1.101011 },
             { 0.788870, 1.489822, 1.256806 },
             { 0.823468, 1.376328, 1.331102 },
             { 0.761771, 1.354508, 1.451228 },
             { 0.792990, 1.271804, 1.513255 },
             { 0.670853, 1.444468, 1.501686 },
             { 0.638527, 1.433485, 1.604139 },
             { 0.695019, 1.578561, 1.301844 },
             { 0.674249, 1.660155, 1.234221 },
             { 0.638711, 1.553840, 1.424819 },
             { 0.562836, 1.623556, 1.457146 },
             { 1.065250, 1.069384, 1.519225 },
             { 1.129825, 0.961334, 1.524583 },
             { 0.968377, 1.080949, 1.598953 }
         }});

//! Coordinates of 54 atoms of the Cys-Val-Trp system after H-Angles were constrained with SHAKE
static const std::vector<RVec> c_cvwXConstrainedHAnglesShake(
        {{
             { 1.766280, 1.431175, 1.064547 },
             { 1.809003, 1.472220, 0.979083 },
             { 1.717227, 1.505634, 1.118062 },
             { 1.840674, 1.389080, 1.123789 },
             { 1.684689, 1.313471, 1.027242 },
             { 1.602869, 1.348475, 0.966071 },
             { 1.763611, 1.225883, 0.928499 },
             { 1.697811, 1.148030, 0.884341 },
             { 1.794194, 1.296915, 0.848744 },
             { 1.911540, 1.142768, 0.993725 },
             { 2.004591, 1.229761, 0.957261 },
             { 1.615064, 1.230131, 1.129241 },
             { 1.645954, 1.111936, 1.143551 },
             { 1.523873, 1.281837, 1.213502 },
             { 1.491052, 1.375919, 1.210256 },
             { 1.454975, 1.211220, 1.317003 },
             { 1.456158, 1.104010, 1.304064 },
             { 1.509331, 1.240710, 1.453659 },
             { 1.501971, 1.346179, 1.487803 },
             { 1.453856, 1.146844, 1.562136 },
             { 1.501190, 1.156948, 1.662149 },
             { 1.345963, 1.171572, 1.571776 },
             { 1.462769, 1.040964, 1.529682 },
             { 1.661570, 1.219597, 1.459277 },
             { 1.705641, 1.291096, 1.386569 },
             { 1.711297, 1.234313, 1.557520 },
             { 1.682879, 1.116333, 1.424266 },
             { 1.308430, 1.237807, 1.313089 },
             { 1.273231, 1.349611, 1.275821 },
             { 1.221804, 1.143240, 1.353615 },
             { 1.250180, 1.050223, 1.375547 },
             { 1.092608, 1.170398, 1.408561 },
             { 1.103566, 1.265544, 1.458468 },
             { 0.978802, 1.168248, 1.305148 },
             { 1.017841, 1.111860, 1.217754 },
             { 0.888800, 1.116321, 1.344450 },
             { 0.922141, 1.300200, 1.258493 },
             { 0.946725, 1.375608, 1.147414 },
             { 1.005252, 1.351224, 1.059984 },
             { 0.873212, 1.491204, 1.148247 },
             { 0.897322, 1.573151, 1.101044 },
             { 0.788903, 1.489750, 1.256848 },
             { 0.823622, 1.376185, 1.330981 },
             { 0.761818, 1.354573, 1.451088 },
             { 0.792960, 1.271885, 1.513193 },
             { 0.670880, 1.444444, 1.501657 },
             { 0.638517, 1.433493, 1.604106 },
             { 0.695050, 1.578522, 1.301840 },
             { 0.674239, 1.660133, 1.234233 },
             { 0.638718, 1.553828, 1.424810 },
             { 0.562845, 1.623555, 1.457144 },
             { 1.065282, 1.069466, 1.519154 },
             { 1.129797, 0.961381, 1.524584 },
             { 0.968422, 1.080939, 1.598913 }
         }});

//! Masses of atoms in the Cys-Val-Trp system
static const std::vector<real> c_cvwMasses ( {
                                                 14.007000,
                                                 1.008000,
                                                 1.008000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 1.008000,
                                                 32.060001,
                                                 1.008000,
                                                 12.011000,
                                                 15.999000,
                                                 14.007000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 1.008000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 1.008000,
                                                 1.008000,
                                                 12.011000,
                                                 15.999000,
                                                 14.007000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 1.008000,
                                                 12.011000,
                                                 12.011000,
                                                 1.008000,
                                                 14.007000,
                                                 1.008000,
                                                 12.011000,
                                                 12.011000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 1.008000,
                                                 12.011000,
                                                 15.999400,
                                                 15.999400
                                             } );

}  // namespace test

}  // namespace gmx
