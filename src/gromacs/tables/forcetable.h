/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
#ifndef GMX_TABLES_FORCETABLE_H
#define GMX_TABLES_FORCETABLE_H

/*! \libinternal \file
 * \brief
 * Old routines for table generation (will eventually be replaced)
 *
 * \inlibraryapi
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"

struct EwaldCorrectionTables;
struct bondedtable_t;
struct interaction_const_t;

/*! \brief Flag to select user tables for make_tables */
#define GMX_MAKETABLES_FORCEUSER (1 << 0)
/*! \brief Flag to only make 1,4 pair tables for make_tables */
#define GMX_MAKETABLES_14ONLY (1 << 1)

//! \brief The types of interactions contained in the table
enum class TableInteraction : int
{
    VdwRepulsionVdwDispersion,
    ElectrostaticVdwRepulsionVdwDispersion,
    Count
};

/* Different formats for table data. Cubic spline tables are typically stored
 * with the four Y,F,G,H intermediate values (check tables.c for format), which
 * makes it easy to load with a single 4-way SIMD instruction too.
 * Linear tables only need one value per table point, or two if both V and F
 * are calculated. However, with SIMD instructions this makes the loads unaligned,
 * and in that case we store the data as F, D=F(i+1)-F(i), V, and then a blank value,
 * which again makes it possible to load as a single instruction.
 */
enum class TableFormat : int
{
    CubicsplineYfgh,
    Count
};

//! \internal \brief Structure describing the data in a single table
struct t_forcetable
{
    t_forcetable(TableInteraction interaction, TableFormat format);

    ~t_forcetable();

    //! Types of interactions stored in this table
    TableInteraction interaction_;
    //! Interpolation type and data format
    TableFormat format_;
    //! range of the table
    real interactionRange;
    //! n+1 is the number of table points
    int numTablePoints;
    //! distance (nm) between two table points
    real scale;
    //! The actual table data
    std::vector<real, gmx::AlignedAllocator<real>> data;

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it
     * makes it much easier to access the tables in the nonbonded kernels when we can set the data
     * from variables. It is always true that stride = formatsize*ninteractions
     */

    //! Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.), only YFGH is implemented
    static constexpr int formatsize = 4;
    //! Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp.
    int numInteractions;
    //! Distance to next table point (number of fp variables per table point in total)
    int stride;
};

/*! \brief Enumerated type to describe the interaction types in a table */
enum
{
    etiCOUL, //!< Coulomb
    etiLJ6,  //!< Dispersion
    etiLJ12, //!< Repulsion
    etiNR    //!< Total number of interaction types
};

/*! \brief Function pointer to calculate the grid contribution for coulomb/LJ
 *
 * Used to tell table_spline3_fill_ewald_lr whether it
 * should calculate the grid contribution for electrostatics or LJ.
 */
typedef double (*real_space_grid_contribution_computer)(double, double);


/*! \brief Construct tables with the Ewald long-range force interaction
 *
 * Creates and fills tables of numPoints points with the spacing
 * set to 1/tableScaling with the Ewald long-range (mesh) force.
 * There are three separate tables with format F, V, FDV0.
 * This function interpolates the Ewald mesh potential contribution
 * with coefficient beta using a quadratic spline.
 * The force is then be interpolated linearly.
 *
 * \param numPoints     Number of points in the tables
 * \param tableScaling  1/spacing, units 1/nm
 * \param beta          Ewald splitting parameter, units 1/nm
 * \param v_lr          Pointer to function calculating real-space grid contribution
 * \returns a set of Ewald correction tables
 */
EwaldCorrectionTables generateEwaldCorrectionTables(int    numPoints,
                                                    double tableScaling,
                                                    real   beta,
                                                    real_space_grid_contribution_computer v_lr);

/*! \brief Compute scaling for the Ewald quadratic spline tables.
 *
 * \param ic                     Pointer to interaction constant structure
 * \param generateCoulombTables  Take the spacing for Coulomb Ewald corrections into account
 * \param generateVdwTables      Take the spacing for Van der Waals Ewald corrections into account
 * \return The scaling factor in units 1/nm
 */
real ewald_spline3_table_scale(const interaction_const_t& ic, bool generateCoulombTables, bool generateVdwTables);

/*! \brief Return the real space grid contribution for Ewald
 *
 *  \param beta  Ewald splitting parameter
 *  \param r     Distance for which to calculate the real-space contrib
 *  \return      Real space grid contribution for Ewald electrostatics
 */
double v_q_ewald_lr(double beta, double r);

/*! \brief Return the real space grid contribution for LJ-Ewald
 *
 *  \param beta  Ewald splitting parameter
 *  \param r     Distance for which to calculate the real-space contrib
 *  \return      Real space grid contribution for Ewald Lennard-Jones interaction
 */
double v_lj_ewald_lr(double beta, double r);

/*! \brief Return tables for inner loops.
 *
 * \param fp     Log file pointer
 * \param ic     Non-bonded interaction constants
 * \param fn     File name from which to read user tables
 * \param rtab   Largest interaction distance to tabulate
 * \param flags  Flags to select table settings
 *
 * \return Pointer to inner loop table structure
 */
std::unique_ptr<t_forcetable>
make_tables(FILE* fp, const interaction_const_t* ic, const char* fn, real rtab, int flags);

/*! \brief Return a table for bonded interactions,
 *
 * \param  fplog   Pointer to log file
 * \param  fn      File name
 * \param  angle   Type of angle: bonds 0, angles 1, dihedrals 2
 * \return New bonded table datatype
 */
bondedtable_t make_bonded_table(FILE* fplog, const char* fn, int angle);

/*! \brief Construct and return tabulated dispersion and repulsion interactions
 *
 * This table can be used to compute long-range dispersion corrections.
 * Returns pointer owning nothing when tabfn=nullptr.
 */
std::unique_ptr<t_forcetable>
makeDispersionCorrectionTable(FILE* fp, const interaction_const_t* ic, real rtab, const char* tabfn);

#endif /* GMX_TABLES_FORCETABLE_H */
