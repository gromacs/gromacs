/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014,2015,2019,2020,2021, by the GROMACS development team, led by
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
#ifndef GMX_MDTYPES_NBLIST_H
#define GMX_MDTYPES_NBLIST_H

#include <vector>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/alignedallocator.h"
#include "gromacs/utility/real.h"

/* The interactions contained in a (possibly merged) table
 * for computing electrostatic, VDW repulsion and/or VDW dispersion
 * contributions.
 */
enum gmx_table_interaction
{
    GMX_TABLE_INTERACTION_ELEC,
    GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWEXPREP_VDWDISP,
    GMX_TABLE_INTERACTION_ELEC_VDWDISP,
    GMX_TABLE_INTERACTION_NR
};

/* Different formats for table data. Cubic spline tables are typically stored
 * with the four Y,F,G,H intermediate values (check tables.c for format), which
 * makes it easy to load with a single 4-way SIMD instruction too.
 * Linear tables only need one value per table point, or two if both V and F
 * are calculated. However, with SIMD instructions this makes the loads unaligned,
 * and in that case we store the data as F, D=F(i+1)-F(i), V, and then a blank value,
 * which again makes it possible to load as a single instruction.
 */
enum gmx_table_format
{
    GMX_TABLE_FORMAT_CUBICSPLINE_YFGH,
    GMX_TABLE_FORMAT_LINEAR_VF,
    GMX_TABLE_FORMAT_LINEAR_V,
    GMX_TABLE_FORMAT_LINEAR_F,
    GMX_TABLE_FORMAT_LINEAR_FDV0,
    GMX_TABLE_FORMAT_NR
};

struct t_nblist
{
    int              nri    = 0; /* Current number of i particles	   */
    int              maxnri = 0; /* Max number of i particles	   */
    int              nrj    = 0; /* Current number of j particles	   */
    int              maxnrj = 0; /* ;Max number of j particles	   */
    std::vector<int> iinr;       /* The i-elements                        */
    std::vector<int> gid;        /* Index in energy arrays                */
    std::vector<int> shift;      /* Shift vector index                    */
    std::vector<int> jindex;     /* Index in jjnr                         */
    std::vector<int> jjnr;       /* The j-atom list                       */
    std::vector<int> excl_fep;   /* Exclusions for FEP with Verlet scheme */
};

/* Structure describing the data in a single table */
struct t_forcetable
{
    t_forcetable(enum gmx_table_interaction interaction, enum gmx_table_format format);

    ~t_forcetable();

    enum gmx_table_interaction interaction; /* Types of interactions stored in this table */
    enum gmx_table_format      format;      /* Interpolation type and data format */

    real r;     /* range of the table */
    int  n;     /* n+1 is the number of table points */
    real scale; /* distance (nm) between two table points */
    std::vector<real, gmx::AlignedAllocator<real>> data; /* the actual table data */

    /* Some information about the table layout. This can also be derived from the interpolation
     * type and the table interactions, but it is convenient to have here for sanity checks, and it
     * makes it much easier to access the tables in the nonbonded kernels when we can set the data
     * from variables. It is always true that stride = formatsize*ninteractions
     */
    int formatsize; /* Number of fp variables for each table point (1 for F, 2 for VF, 4 for YFGH, etc.) */
    int ninteractions; /* Number of interactions in table, 1 for coul-only, 3 for coul+rep+disp. */
    int stride; /* Distance to next table point (number of fp variables per table point in total) */
};

#endif /* GMX_MDTYPES_NBLIST_H */
