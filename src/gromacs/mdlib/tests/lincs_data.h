/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#define ONE_TENTH_OVER_SQRT_TWO  0.0707106781 
#define TWO_TENTH_OVER_SQRT_TREE 0.1154700538 

namespace gmx
{

namespace test
{

static const matrix c_infinitesimalBox = { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} };

/*
 * Two connected atoms.
 */

static const std::vector<int>  c_oneBondConstraints( {0, 0, 1} );
static const std::vector<real> c_oneBondConstraintsR0( {0.1} );
static const int c_oneBondConstraintsFirstType = 0;
static const std::vector<RVec> c_oneBondCoordinates( 
        {{
            { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
            { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 }
        }} );
static const std::vector<real> c_oneBondMasses ( {1.0, 12.0} );

/*
 * Three atoms, connected longitudaly.
 */

static const std::vector<int>  c_twoBondsConstraints( {0, 0, 1, 1, 1, 2} );
static const std::vector<real> c_twoBondsConstraintsR0( {0.1, 0.2} );
static const int c_twoBondsConstraintsFirstType = 0;
static const std::vector<RVec> c_twoBondsCoordinates( 
        {{
            { ONE_TENTH_OVER_SQRT_TWO, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
            { 0.0, 0.0, 0.0 },
            { TWO_TENTH_OVER_SQRT_TREE, TWO_TENTH_OVER_SQRT_TREE, TWO_TENTH_OVER_SQRT_TREE }
        }} ); 
static const std::vector<real> c_twoBondsMasses ( {1.0, 12.0, 16.0 } );


/*
 * Basic triangle
 */

static const std::vector<int>  c_triangleConstraints( {0, 0, 1, 2, 0, 2, 1, 1, 2} );
static const std::vector<real> c_triangleConstraintsR0( {0.1, 0.1, 0.1} );
static const int c_triangleConstraintsFirstType = 0;
static const std::vector<RVec> c_triangleCoordinates( 
        {{
            { ONE_TENTH_OVER_SQRT_TWO, 0.0, 0.0 },
            { 0.0, ONE_TENTH_OVER_SQRT_TWO, 0.0 },
            { 0.0, 0.0, ONE_TENTH_OVER_SQRT_TWO }
        }} );
static const std::vector<real> c_triangleMasses ( {1.0, 1.0, 1.0} );



/*
 * Real-life system: Cys-Val-Trp peptide.
 */

static const matrix c_cvwBox = {{real(2.570950), 0, 0}, {0, real(2.570950), 0}, {0, 0, real(2.570950)}};


// Constraints only on covalent bonds with hydrogens

static const std::vector<int> c_cvwHBondsConstraints(
        {465,  0,  1, 465,  0,  2, 465,  0,  3, 466,  4,  5, 467,  6,  7, 467,  6,  8, 468,  9, 10, 469, 13, 14, 
         466, 15, 16, 467, 17, 18, 467, 19, 20, 467, 19, 21, 467, 19, 22, 467, 23, 24, 467, 23, 25, 467, 23, 26, 
         469, 29, 30, 466, 31, 32, 467, 33, 34, 467, 33, 35, 466, 37, 38, 470, 39, 40, 466, 43, 44, 466, 45, 46, 
         466, 47, 48, 466, 49, 50} );

static const std::vector<real> c_cvwHBondsConstraintsR0(
        {0.104000, 0.108000, 0.111100, 0.132500, 0.099700, 0.097600} );
        
static const int c_cvwHBondsConstraintsFirstType = 465;


// Constraints on all covalent bonds

static const std::vector<int> c_cvwAllBondsConstraints(
        {447,  0,  1, 447,  0,  2, 447,  0,  3, 448,  0,  4, 449,  4,  5, 450,  4,  6, 451,  4, 11, 452,  6,  7, 
         452,  6,  8, 453,  6,  9, 454,  9, 10, 455, 11, 12, 456, 11, 13, 457, 13, 14, 458, 13, 15, 449, 15, 16, 
         459, 15, 17, 451, 15, 27, 452, 17, 18, 450, 17, 19, 450, 17, 23, 452, 19, 20, 452, 19, 21, 452, 19, 22, 
         452, 23, 24, 452, 23, 25, 452, 23, 26, 455, 27, 28, 456, 27, 29, 457, 29, 30, 458, 29, 31, 449, 31, 32, 
         450, 31, 33, 460, 31, 51, 452, 33, 34, 452, 33, 35, 461, 33, 36, 462, 36, 37, 463, 36, 42, 449, 37, 38, 
         464, 37, 39, 465, 39, 40, 466, 39, 41, 467, 41, 42, 468, 41, 47, 468, 42, 43, 449, 43, 44, 466, 43, 45, 
         449, 45, 46, 466, 45, 49, 449, 47, 48, 466, 47, 49, 449, 49, 50, 469, 51, 52, 469, 51, 53} );

static const std::vector<real> c_cvwAllBondsConstraintsR0(   
        {0.104000, 0.148000, 0.108000, 0.153800, 0.149000, 0.111100, 0.181800, 0.132500, 0.123000, 0.134500, 
         0.099700, 0.143000, 0.150000, 0.152200, 0.151000, 0.136500, 0.144000, 0.137000, 0.097600, 0.137500, 
         0.140000, 0.136800, 0.126000} );
         
static const int c_cvwAllBondsConstraintsFirstType = 447;



//! Coordinates of all 54 atoms of the Cys-Val-Trp system before constraints were applied
static const std::vector<RVec> c_cvwInitialCoordinates( 
          {{
            { 1.746000, 1.404000, 1.052000 },
            { 1.806000, 1.450000, 0.979000 },
            { 1.672000, 1.473000, 1.086000 },
            { 1.814000, 1.368000, 1.123000 },
            { 1.681000, 1.283000, 0.992000 },
            { 1.612000, 1.309000, 0.913000 },
            { 1.801000, 1.204000, 0.919000 },
            { 1.772000, 1.100000, 0.896000 },
            { 1.814000, 1.266000, 0.827000 },
            { 1.974000, 1.193000, 0.998000 },
            { 2.058000, 1.204000, 0.896000 },
            { 1.594000, 1.204000, 1.095000 },
            { 1.568000, 1.083000, 1.091000 },
            { 1.533000, 1.270000, 1.201000 },
            { 1.540000, 1.372000, 1.213000 },
            { 1.458000, 1.206000, 1.310000 },
            { 1.457000, 1.094000, 1.307000 },
            { 1.532000, 1.230000, 1.433000 },
            { 1.527000, 1.340000, 1.449000 },
            { 1.478000, 1.148000, 1.544000 },
            { 1.544000, 1.157000, 1.633000 },
            { 1.378000, 1.176000, 1.581000 },
            { 1.461000, 1.043000, 1.518000 },
            { 1.685000, 1.181000, 1.411000 },
            { 1.747000, 1.243000, 1.345000 },
            { 1.730000, 1.190000, 1.515000 },
            { 1.688000, 1.075000, 1.377000 },
            { 1.312000, 1.253000, 1.324000 },
            { 1.268000, 1.367000, 1.318000 },
            { 1.224000, 1.149000, 1.344000 },
            { 1.264000, 1.058000, 1.349000 },
            { 1.095000, 1.170000, 1.405000 },
            { 1.103000, 1.263000, 1.454000 },
            { 0.975000, 1.166000, 1.297000 },
            { 1.025000, 1.113000, 1.214000 },
            { 0.887000, 1.110000, 1.333000 },
            { 0.924000, 1.306000, 1.265000 },
            { 0.952000, 1.386000, 1.165000 },
            { 1.017000, 1.361000, 1.079000 },
            { 0.882000, 1.507000, 1.174000 },
            { 0.911000, 1.589000, 1.129000 },
            { 0.819000, 1.515000, 1.296000 },
            { 0.840000, 1.389000, 1.355000 },
            { 0.777000, 1.364000, 1.480000 },
            { 0.785000, 1.264000, 1.522000 },
            { 0.698000, 1.458000, 1.538000 },
            { 0.640000, 1.433000, 1.626000 },
            { 0.746000, 1.617000, 1.350000 },
            { 0.719000, 1.711000, 1.307000 },
            { 0.682000, 1.582000, 1.472000 },
            { 0.631000, 1.668000, 1.515000 },
            { 1.070000, 1.060000, 1.508000 },
            { 1.150000, 0.960000, 1.515000 },
            { 0.966000, 1.066000, 1.569000 }
	    }});        
        
//! Coordinates of 54 atoms of the Cys-Val-Trp system after H-Bonds constaints were constrained 
static const std::vector<RVec> c_cvwFinalCoordinatesHBonds( 
          {{        
            { 1.745950, 1.404137, 1.052042 },
            { 1.805376, 1.449522, 0.979759 },
            { 1.673804, 1.471318, 1.085171 },
            { 1.813515, 1.368257, 1.122494 },
            { 1.680997, 1.283001, 0.991996 },
            { 1.612038, 1.308986, 0.913044 },
            { 1.801019, 1.204075, 0.918974 },
            { 1.771832, 1.099398, 0.895867 },
            { 1.813938, 1.265703, 0.827441 },
            { 1.974002, 1.193000, 0.997998 },
            { 2.057943, 1.203992, 0.896070 },
            { 1.594000, 1.204000, 1.095000 },
            { 1.568000, 1.083000, 1.091000 },
            { 1.533015, 1.270216, 1.201025 },
            { 1.539794, 1.369004, 1.212648 },
            { 1.457997, 1.205687, 1.309992 },
            { 1.457033, 1.097730, 1.307100 },
            { 1.531999, 1.230013, 1.433002 },
            { 1.527007, 1.339845, 1.448977 },
            { 1.478082, 1.148102, 1.544007 },
            { 1.543998, 1.157000, 1.632997 },
            { 1.377262, 1.176207, 1.581273 },
            { 1.460769, 1.041573, 1.517647 },
            { 1.685022, 1.180937, 1.411234 },
            { 1.747673, 1.243673, 1.344284 },
            { 1.729067, 1.189813, 1.512844 },
            { 1.687993, 1.075258, 1.377083 },
            { 1.312000, 1.253000, 1.324000 },
            { 1.268000, 1.367000, 1.318000 },
            { 1.223995, 1.149011, 1.343999 },
            { 1.264064, 1.057854, 1.349008 },
            { 1.094985, 1.169824, 1.404907 },
            { 1.103180, 1.265097, 1.455105 },
            { 0.975024, 1.166056, 1.297020 },
            { 1.025283, 1.112700, 1.213531 },
            { 0.886430, 1.109637, 1.333233 },
            { 0.924000, 1.306000, 1.265000 },
            { 0.952121, 1.385953, 1.164840 },
            { 1.015558, 1.361555, 1.080908 },
            { 0.882007, 1.507018, 1.173990 },
            { 0.910909, 1.588743, 1.129141 },
            { 0.819000, 1.515000, 1.296000 },
            { 0.840000, 1.389000, 1.355000 },
            { 0.777004, 1.363946, 1.480023 },
            { 0.784949, 1.264642, 1.521730 },
            { 0.697987, 1.457994, 1.538020 },
            { 0.640158, 1.433068, 1.625761 },
            { 0.746023, 1.616921, 1.350036 },
            { 0.718729, 1.711945, 1.306568 },
            { 0.681970, 1.582051, 1.472026 },
            { 0.631363, 1.667388, 1.514694 },
            { 1.070000, 1.060000, 1.508000 },
            { 1.150000, 0.960000, 1.515000 },
            { 0.966000, 1.066000, 1.569000 }
        }});

//! Coordinates of 54 atoms of the Cys-Val-Trp system after All bonds constraints were applied 
static const std::vector<RVec> c_cvwFinalCoordinatesAllBonds( 
          {{
            { 1.744192, 1.400876, 1.050408 },
            { 1.804620, 1.448942, 0.980679 },
            { 1.674799, 1.470390, 1.084714 },
            { 1.812763, 1.368655, 1.121708 },
            { 1.683948, 1.276865, 0.996242 },
            { 1.616170, 1.307429, 0.917774 },
            { 1.803010, 1.207954, 0.927155 },
            { 1.773417, 1.105083, 0.897124 },
            { 1.813404, 1.263159, 0.831216 },
            { 1.970550, 1.193223, 0.996389 },
            { 2.057327, 1.203912, 0.896817 },
            { 1.595071, 1.207624, 1.094104 },
            { 1.568932, 1.087338, 1.091143 },
            { 1.533602, 1.264980, 1.199359 },
            { 1.539437, 1.363790, 1.212034 },
            { 1.455744, 1.209606, 1.306025 },
            { 1.457067, 1.101535, 1.307202 },
            { 1.536585, 1.230178, 1.430802 },
            { 1.527027, 1.339398, 1.448912 },
            { 1.477277, 1.146856, 1.545698 },
            { 1.544449, 1.157061, 1.633605 },
            { 1.376383, 1.176453, 1.581598 },
            { 1.460628, 1.040705, 1.517432 },
            { 1.681448, 1.182094, 1.411715 },
            { 1.746716, 1.242717, 1.345302 },
            { 1.728696, 1.189739, 1.511987 },
            { 1.687955, 1.076581, 1.377507 },
            { 1.313454, 1.250569, 1.323999 },
            { 1.268736, 1.365093, 1.318100 },
            { 1.224082, 1.152044, 1.344591 },
            { 1.262989, 1.060300, 1.348874 },
            { 1.093660, 1.168520, 1.401086 },
            { 1.102950, 1.262423, 1.453696 },
            { 0.977797, 1.168815, 1.299821 },
            { 1.024283, 1.113760, 1.215191 },
            { 0.888605, 1.111021, 1.332343 },
            { 0.921772, 1.305848, 1.269492 },
            { 0.951364, 1.387941, 1.164422 },
            { 1.015236, 1.361678, 1.081333 },
            { 0.882251, 1.505804, 1.174960 },
            { 0.910507, 1.587606, 1.129765 },
            { 0.819059, 1.513487, 1.296914 },
            { 0.840064, 1.387880, 1.355340 },
            { 0.778667, 1.365803, 1.475722 },
            { 0.784710, 1.267623, 1.520478 },
            { 0.698397, 1.459946, 1.535866 },
            { 0.641444, 1.433622, 1.623810 },
            { 0.745363, 1.614637, 1.352384 },
            { 0.719425, 1.709521, 1.307677 },
            { 0.683959, 1.580100, 1.470510 },
            { 0.632793, 1.664977, 1.513489 },
            { 1.072860, 1.060117, 1.505945 },
            { 1.149574, 0.960533, 1.514963 },
            { 0.964441, 1.066090, 1.569914 }
        }});

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
              15.999400 } );


}  // namespace test

}  // namespace gmx
