/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#ifndef GMX_TRAJECTORYANALYSIS_MODULES_DOSUTILS_H
#define GMX_TRAJECTORYANALYSIS_MODULES_DOSUTILS_H

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_output_env_t;

//! Type of velocities: trans, angular, vibrational, rotational, total
enum Velocity {
    V_T, V_A, V_V, V_R, V_ATOM, V_SUM, V_NR
};

//! Component phase of the density of states calculationk
enum Component {
    C_SOLID, C_GAS, C_TOTAL, C_NR
};

//! Thermodynamics algorithm
enum ThermoMethod {
    TM_1PT, TM_2PT, TM_MD, TM_NR
};

//! Level of theory, pure classical, quantum or difference
enum LevelOfTheory {
    LOT_CLASSICAL, LOT_QUANTUM, LOT_QUANTUM_CORRECTION, LOT_NR
};

//! Thermodynamics property to evaluate
enum Density {
    DOS_CV, DOS_S, DOS_A, DOS_E, DOS_NR
};

const char *lotName(LevelOfTheory lot);

void dsub_xcm(rvec x[], int gnx, int *index, t_atom atom[], rvec xcm);

/* eq 34 JCP 119, 11792 (2003) */
double FD(double Delta, double f);

/* fluidicity eq. for hard spheres packing, eq 31 JCP 119, 11792 (2003) */
double YYY(double f, double y);

/* compressibility of hard spheres z=(1+y+y^2-y^3)/(1-y)^3 */
double calc_compress(double y);

double bisector(double Delta, double tol,
                double ff0, double ff1,
                double ff(double, double));

/* calculate fluidicity f */
double calc_fluidicity(double Delta, double tol);

/* hard sphere packing fraction, eq 32 JCP 119, 11792 (2003) */
double calc_y(double f, double Delta, double toler);

/* entropy for hard spheres */
double calc_Shs(double fy);

real weight(int Density, int LevelOfTheory,
            int ThermoMethod, int Component,
            real nu, real beta);

#define NDIM 4
void principal(int n, int index[], t_atom atom[], rvec x[],
               matrix trans, rvec d, double **inten,
               double **ev);

struct t_energy_term
{
    int    ftype, fff, nesum;
    bool   bCheckDrift;
    double energy, stddev;
    double esum, e2sum, esum0, e2sum0;
};

void get_energy_terms2(const char *enx_fn, int net, t_energy_term et[]);

void dump_w(gmx_output_env_t *oenv, real beta);

#endif
