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
#ifndef GMX_TOPOLOGY_IFUNC_H
#define GMX_TOPOLOGY_IFUNC_H

#include "gromacs/libgromacs_export.h"
#include "gromacs/math/vectypes.h"

struct t_fcdata;
struct t_graph;
union t_iparams;
struct t_mdatoms;
struct t_pbc;

/* TODO: Remove this typedef when t_ilist is removed */
typedef int t_iatom;

/* Real vector type with an additional, unused 4th element.
 * This type is used to allow aligned 4-wide SIMD loads and stores.
 */
typedef real rvec4[4];

/*
 * The function type t_ifunc() calculates one interaction, using iatoms[]
 * and iparams. Within the function the number of atoms to be used is
 * known. Within the function only the atomid part of the iatoms[] array
 * is supplied, not the type field (see also t_ilist). The function
 * returns the potential energy. If pbc==NULL the coordinates in x are
 * assumed to be such that no calculation of PBC is necessary,
 * If pbc!=NULL a full PBC calculation is performed.
 * If g!=NULL it is used for determining the shift forces.
 * With domain decomposition ddgatindex can be used for getting global
 * atom numbers for warnings and error messages.
 * ddgatindex is NULL when domain decomposition is not used.
 */

constexpr unsigned int IF_NULL       = 0;
constexpr unsigned int IF_BOND       = 1 << 0;
constexpr unsigned int IF_VSITE      = 1 << 1;
constexpr unsigned int IF_CONSTRAINT = 1 << 2;
constexpr unsigned int IF_CHEMBOND   = 1 << 3;
constexpr unsigned int IF_BTYPE      = 1 << 4;
constexpr unsigned int IF_ATYPE      = 1 << 5;
constexpr unsigned int IF_DIHEDRAL   = 1 << 6;
constexpr unsigned int IF_PAIR       = 1 << 7;
constexpr unsigned int IF_TABULATED  = 1 << 8;
constexpr unsigned int IF_LIMZERO    = 1 << 9;
/* These flags tell to some of the routines what can be done with this
 * item in the list.
 * With IF_BOND a bonded interaction will be calculated.
 * With IF_BTYPE grompp can convert the bond to a Morse potential.
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 * IF_LIMZERO indicates that for a bonded interaction the potential
 * does goes to zero for large distances, thus if such an interaction
 * it not assigned to any node by the domain decomposition, the simulation
 * still continue, if mdrun has been told so.
 */

struct t_interaction_function // NOLINT (clang-analyzer-optin.performance.Padding)
{
    const char* name;         /* the name of this function			*/
    const char* longname;     /* The name for printing etc.                   */
    int         nratoms;      /* nr of atoms needed for this function		*/
    int         nrfpA, nrfpB; /* number of parameters for this function.      */
                              /* this corresponds to the number of params in  */
                              /* iparams struct! (see idef.h)                 */
    /* A and B are for normal and free energy components respectively.    */
    unsigned int flags; /* Flags (see above)                            */
};

#define NRFPA(ftype) (interaction_function[(ftype)].nrfpA)
#define NRFPB(ftype) (interaction_function[(ftype)].nrfpB)
#define NRFP(ftype) (NRFPA(ftype) + NRFPB(ftype))
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)

#define IS_CHEMBOND(ftype) \
    (interaction_function[(ftype)].nratoms == 2 && (interaction_function[(ftype)].flags & IF_CHEMBOND))
/* IS_CHEMBOND tells if function type ftype represents a chemical bond */

/* IS_ANGLE tells if a function type ftype represents an angle
 * Per Larsson, 2007-11-06
 */
#define IS_ANGLE(ftype) \
    (interaction_function[(ftype)].nratoms == 3 && (interaction_function[(ftype)].flags & IF_ATYPE))
#define IS_VSITE(ftype) (interaction_function[(ftype)].flags & IF_VSITE)

#define IS_TABULATED(ftype) (interaction_function[(ftype)].flags & IF_TABULATED)

/* this MUST correspond to the
   t_interaction_function[F_NRE] in src/gromacs/topology/ifunc.cpp */
enum
{
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12_NOLONGERUSED,
    F_GB13_NOLONGERUSED,
    F_GB14_NOLONGERUSED,
    F_GBPOL_NOLONGERUSED,
    F_NPSOLVATION_NOLONGERUSED,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR_NOLONGERUSED,
    F_BHAM_LR_NOLONGERUSED,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR_NOLONGERUSED,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE1,
    F_VSITE2,
    F_VSITE2FD,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_DENSITYFITTING,
    F_EQM,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE /* This number is for the total number of energies      */
};

static inline bool IS_RESTRAINT_TYPE(int ifunc)
{
    return ifunc == F_POSRES || ifunc == F_FBPOSRES || ifunc == F_DISRES || ifunc == F_RESTRBONDS
           || ifunc == F_DISRESVIOL || ifunc == F_ORIRES || ifunc == F_ORIRESDEV
           || ifunc == F_ANGRES || ifunc == F_ANGRESZ || ifunc == F_DIHRES;
}

/* Maximum allowed number of atoms, parameters and terms in interaction_function.
 * Check kernel/toppush.c when you change these numbers.
 */
constexpr int MAXATOMLIST   = 6;
constexpr int MAXFORCEPARAM = 12;
constexpr int NR_RBDIHS     = 6;
constexpr int NR_CBTDIHS    = 6;
constexpr int NR_FOURDIHS   = 4;

LIBGROMACS_EXPORT extern const t_interaction_function interaction_function[F_NRE];
/* initialised interaction functions descriptor				*/

#endif
