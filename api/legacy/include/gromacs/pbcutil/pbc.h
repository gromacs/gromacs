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
#ifndef GMX_PBCUTIL_PBC_H
#define GMX_PBCUTIL_PBC_H

#include <cstdio>

#include <string>

#include "gromacs/libgromacs_export.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

struct gmx_domdec_t;
struct gmx_mtop_t;

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

//! Names for all values in PBC types enumeration
extern LIBGROMACS_EXPORT const gmx::EnumerationArray<PbcType, std::string> c_pbcTypeNames;

/* Maximum number of combinations of single triclinic box vectors
 * required to shift atoms that are within a brick of the size of
 * the diagonal of the box to within the maximum cut-off distance.
 */
#define MAX_NTRICVEC 12

/*! \brief Structure containing info on periodic boundary conditions */
typedef struct t_pbc
{
    //! The PBC type
    PbcType pbcType;
    //! Number of dimensions in which PBC is exerted
    int ndim_ePBC;
    /*! \brief Determines how to compute distance vectors.
     *
     *  Indicator of how to compute distance vectors, depending
     *  on PBC type (depends on pbcType and dimensions with(out) DD)
     *  and the box angles.
     */
    int pbcTypeDX;
    /*! \brief Used for selecting which dimensions to use in PBC.
     *
     *  In case of 1-D PBC this indicates which dimension is used,
     *  in case of 2-D PBC this indicates the opposite
     */
    int dim;
    //! The simulation box
    matrix box;
    //! The lengths of the diagonal of the full box
    rvec fbox_diag;
    //! Halve of the above
    rvec hbox_diag;
    //! Negative of the above
    rvec mhbox_diag;
    //! Maximum allowed cutoff squared for the box and PBC used
    real max_cutoff2;
    /*! \brief Number of triclinic shift vectors.
     *
     *  Number of triclinic shift vectors depends on the skewedness
     *  of the box, that is mostly on the angles. For triclinic boxes
     *  we first take the closest image along each Cartesian dimension
     *  independently. When the resulting distance^2 is larger than
     *  max_cutoff2, up to ntric_vec triclinic shift vectors need to
     *  be tried. Because of the restrictions imposed on the unit-cell
     *  by GROMACS, ntric_vec <= MAX_NTRICVEC = 12.
     */
    int ntric_vec;
    //! The triclinic shift vectors in grid cells. Internal use only.
    ivec tric_shift[MAX_NTRICVEC];
    //!  The triclinic shift vectors in length units
    rvec tric_vec[MAX_NTRICVEC];
} t_pbc;

#define TRICLINIC(box) ((box)[YY][XX] != 0 || (box)[ZZ][XX] != 0 || (box)[ZZ][YY] != 0)

#define NTRICIMG 14
#define NCUCVERT 24
#define NCUCEDGE 36

enum
{
    ecenterTRIC, /* 0.5*(a+b+c)                  */
    ecenterRECT, /* (0.5*a[x],0.5*b[y],0.5*c[z]) */
    ecenterZERO, /* (0,0,0)                      */
    ecenterDEF = ecenterTRIC
};

/*! \brief Returns the number of dimensions that use pbc
 *
 * \param[in] pbcType The periodic boundary condition type
 * \return the number of dimensions that use pbc, starting at X
 */
int numPbcDimensions(PbcType pbcType);

/*! \brief Dump the contents of the pbc structure to the file
 *
 * \param[in] fp  The file pointer to write to
 * \param[in] pbc The periodic boundary condition information structure
 */
void dump_pbc(FILE* fp, t_pbc* pbc);

/*! \brief Check the box for consistency
 *
 * When \p pbcType=PbcTypes::Unset, the type of pbc is guessed from the box matrix.
 *
 * \param[in] pbcType The pbc identifier
 * \param[in] box     The box matrix
 * \return NULL if the box is supported by Gromacs.
 *         Otherwise returns a string with the problem.
 */
const char* check_box(PbcType pbcType, const matrix box);

/*! \brief Creates box matrix from edge lengths and angles.
 *
 * \param[in,out] box        The box matrix
 * \param[in] vec            The edge lengths
 * \param[in] angleInDegrees The angles
 */
void matrix_convert(matrix box, const rvec vec, const rvec angleInDegrees);

/*! \brief Compute the maximum cutoff for the box

 * Returns the square of the maximum cut-off allowed for the box,
 * taking into account that the grid neighborsearch code and pbc_dx
 * only check combinations of single box-vector shifts.
 *
 * \param[in] pbcType The pbc identifier
 * \param[in] box  The box matrix
 * \return the maximum cut-off.
 */
real max_cutoff2(PbcType pbcType, const matrix box);

/*! \brief Guess PBC type
 *
 * Guesses the type of periodic boundary conditions using the box
 *
 * \param[in] box  The box matrix
 * \return The pbc type identifier
 */
PbcType guessPbcType(const matrix box);

/*! \brief Corrects the box if necessary
 *
 * Checks for un-allowed box angles and corrects the box.
 *
 * \param[in] fplog File for debug output
 * \param[in] step  The MD step number
 * \param[in] box   The simulation cell
 * \return TRUE when the box was corrected.
 */
bool correct_box(FILE* fplog, int64_t step, tensor box);

/*! \brief Initiate the periodic boundary condition algorithms.
 *
 * pbc_dx will not use pbc and return the normal difference vector
 * when one or more of the diagonal elements of box are zero.
 * When \p pbcType=PbcType::Unset, the type of pbc is guessed from the box matrix.
 *
 * \param[in,out] pbc The pbc information structure
 * \param[in] pbcType The PBC identifier
 * \param[in] box     The box tensor
 */
void set_pbc(t_pbc* pbc, PbcType pbcType, const matrix box);

/*! \brief Initiate the periodic boundary condition algorithms.
 *
 * As set_pbc, but additionally sets that correct distances can
 * be obtained using (combinations of) single box-vector shifts.
 * Should be used with pbc_dx_aiuc.
 * If domdecCells!=NULL pbc is not used for directions
 * with dd->nc[i]==1 with bSingleDir==TRUE or
 * with dd->nc[i]<=2 with bSingleDir==FALSE.
 * Note that when no PBC is required only pbc->pbcType is set,
 * the rest of the struct will be invalid.
 *
 * \param[in,out] pbc     The pbc information structure
 * \param[in] pbcType     The PBC identifier
 * \param[in] domdecCells 3D integer vector describing the number of DD cells
 *                        or nullptr if not using DD.
 * \param[in] bSingleDir  TRUE if DD communicates only in one direction along dimensions
 * \param[in] box         The box tensor
 * \return the pbc structure when pbc operations are required, NULL otherwise.
 */
t_pbc* set_pbc_dd(t_pbc* pbc, PbcType pbcType, const gmx::IVec* domdecCells, bool bSingleDir, const matrix box);

/*! \brief Compute distance with PBC
 *
 * Calculate the correct distance vector from x2 to x1 and put it in dx.
 * set_pbc must be called before ever calling this routine.
 *
 * Note that for triclinic boxes that do not obey the GROMACS unit-cell
 * restrictions, pbc_dx and pbc_dx_aiuc will not correct for PBC.
 * \param[in,out] pbc The pbc information structure
 * \param[in]    x1  Coordinates for particle 1
 * \param[in]    x2  Coordinates for particle 2
 * \param[out]   dx  Distance vector
 */
void pbc_dx(const t_pbc* pbc, const rvec x1, const rvec x2, rvec dx);

/*! \brief Compute distance vector for simple PBC types
 *
 * Calculate the correct distance vector from x2 to x1 and put it in dx,
 * This function can only be used when all atoms are in the rectangular
 * or triclinic unit-cell.
 * set_pbc_dd or set_pbc must be called before ever calling this routine.
 * \param[in,out] pbc The pbc information structure
 * \param[in]    x1  Coordinates for particle 1
 * \param[in]    x2  Coordinates for particle 2
 * \param[out]   dx  Distance vector
 * \return the ishift required to shift x1 at closest distance to x2;
 * i.e. if 0<=ishift<c_numShiftVectors then x1 - x2 + shift_vec[ishift] = dx
 * (see calc_shifts below on how to obtain shift_vec)
 */
int pbc_dx_aiuc(const t_pbc* pbc, const rvec x1, const rvec x2, rvec dx);

/*! \brief Compute distance with PBC
 *
 * As pbc_dx, but for double precision vectors.
 * set_pbc must be called before ever calling this routine.
 * \param[in,out] pbc The pbc information structure
 * \param[in]    x1  Coordinates for particle 1
 * \param[in]    x2  Coordinates for particle 2
 * \param[out]   dx  Distance vector
 */
void pbc_dx_d(const t_pbc* pbc, const dvec x1, const dvec x2, dvec dx);

/*! \brief Computes shift vectors
 *
 * This routine calculates ths shift vectors necessary to use the
 * neighbor searching routine.
 * \param[in]  box       The simulation box
 * \param[out] shift_vec The shifting vectors
 */
void calc_shifts(const matrix box, gmx::ArrayRef<gmx::RVec> shift_vec);

/*! \brief Calculates the center of the box.
 *
 * See the description for the enum ecenter above.
 * \param[in]  ecenter    Description of center type
 * \param[in]  box        The simulation box
 * \param[out] box_center The center of the box
 */
void calc_box_center(int ecenter, const matrix box, rvec box_center);

/*! \brief Calculates the NTRICIMG box images
 *
 * \param[in]  box The simulation box
 * \param[out] img The triclinic box images
 */
void calc_triclinic_images(const matrix box, rvec img[]);

/*! \brief Calculates the NCUCVERT vertices of a compact unitcell
 *
 * \param[in]  ecenter The center type
 * \param[in]  box     The simulation box
 * \param[out] vert    The vertices
 */
void calc_compact_unitcell_vertices(int ecenter, const matrix box, rvec vert[]);

/*! \brief Compute unitcell edges
 *
 * \return an array of unitcell edges of length NCUCEDGE*2,
 * this is an index in vert[], which is calculated by calc_unitcell_vertices.
 * The index consists of NCUCEDGE pairs of vertex indices.
 * The index does not change, so it needs to be retrieved only once.
 */
int* compact_unitcell_edges();

/*! \brief Put atoms inside the simulations box
 *
 * These routines puts ONE or ALL atoms in the box, not caring
 * about charge groups!
 * Also works for triclinic cells.
 *
 * \param[in]     pbcType The pbc type
 * \param[in]     box     The simulation box
 * \param[in,out] x       The coordinates of the atoms
 */
void put_atoms_in_box(PbcType pbcType, const matrix box, gmx::ArrayRef<gmx::RVec> x);

/*! \brief Parallellizes put_atoms_in_box()
 *
 * This wrapper function around put_atoms_in_box() with the ugly manual
 * workload splitting is needed to avoid silently introducing multithreading
 * in tools.
 *
 * \param[in]     pbcType    The pbc type
 * \param[in]     box        The simulation box
 * \param[in]     haveBoxDeformation  Whether the box is being continously deformed
 * \param[in]     boxDeformation      The deformation speed of the box components in units of nm/ps
 * \param[in,out] x          The coordinates of the atoms
 * \param[in,out] v          The velocities of the atoms
 * \param[in]     nth        number of threads to be used in the given module
 */
void put_atoms_in_box_omp(PbcType                  pbcType,
                          const matrix             box,
                          bool                     haveBoxDeformation,
                          const matrix             boxDeformation,
                          gmx::ArrayRef<gmx::RVec> x,
                          gmx::ArrayRef<gmx::RVec> v,
                          gmx_unused int           nth);

/*! \brief Put atoms inside triclinic box
 *
 * This puts ALL atoms in the triclinic unit cell, centered around the
 * box center as calculated by calc_box_center.
 * \param[in]    ecenter The pbc center type
 * \param[in]    box     The simulation box
 * \param[in,out] x       The coordinates of the atoms
 */
void put_atoms_in_triclinic_unitcell(int ecenter, const matrix box, gmx::ArrayRef<gmx::RVec> x);

/*! \brief Put atoms inside the unitcell
 *
 * This puts ALL atoms at the closest distance for the center of the box
 * as calculated by calc_box_center.
 * When \p pbcType=PbcType::Unset, the type of pbc is guessed from the box matrix.
 *
 * \param[in]    pbcType The pbc type
 * \param[in]    ecenter The pbc center type
 * \param[in]    box     The simulation box
 * \param[in,out] x      The coordinates of the atoms
 */
void put_atoms_in_compact_unitcell(PbcType pbcType, int ecenter, const matrix box, gmx::ArrayRef<gmx::RVec> x);

/*! \brief Make all molecules whole by shifting positions
 *
 * \param[in]     fplog     Log file
 * \param[in]     pbcType   The PBC type
 * \param[in]     correctVelocitiesForBoxDeformation  Whether to correct the velocities for continous box deformation
 * \param[in]     boxDeformation  The box deformation velocity
 * \param[in]     box       The simulation box
 * \param[in]     mtop      System topology definition
 * \param[in,out] x         The coordinates of the atoms
 * \param[in,out] v         The velocities of the atoms, needed with correctVelocitiesForBoxDeformation
 */
void do_pbc_first_mtop(FILE*                    fplog,
                       PbcType                  pbcType,
                       bool                     correctVelocitiesForBoxDeformation,
                       const matrix             boxDeformation,
                       const matrix             box,
                       const gmx_mtop_t*        mtop,
                       gmx::ArrayRef<gmx::RVec> x,
                       gmx::ArrayRef<gmx::RVec> v);

/*! \brief Make molecules consisting of multiple charge groups whole by shifting positions
 *
 * \param[in]     pbcType   The PBC type
 * \param[in]     box       The simulation box
 * \param[in]     mtop      System topology definition
 * \param[in,out] x         The coordinates of the atoms
 */
void do_pbc_mtop(PbcType pbcType, const matrix box, const gmx_mtop_t* mtop, rvec x[]);


/*! \brief Sets the box deformation rate
 *
 * \param[in]      boxDeformation      The deformation speed of the box components in units of nm/ps
 * \param[in]      box                 The box
 * \param[out]     boxDeformationRate  The deformation rate of the box in units of 1/ps
 */
void setBoxDeformationRate(const matrix boxDeformation, const matrix box, matrix boxDeformationRate);

/*! \brief Correct the velocity of a particle for displacement along a flow field
 *
 * \tparam         invertDisplacement  When true, correct for -displacement
 * \param[in]      boxDeformationRate  The deformation rate of the box in units of 1/ps
 * \param[in,out]  v                   The velocity field to be corrected
 * \param[in]      displacement        The coordinate displacement to correct for
 */
template<bool invertDisplacement>
static inline void correctVelocityForDisplacement(const matrix boxDeformationRate, rvec v, const rvec displacement)
{
    for (int d1 = 0; d1 < DIM; d1++)
    {
        for (int d2 = 0; d2 <= d1; d2++)
        {
            if constexpr (invertDisplacement)
            {
                v[d2] -= boxDeformationRate[d1][d2] * displacement[d1];
            }
            else
            {
                v[d2] += boxDeformationRate[d1][d2] * displacement[d1];
            }
        }
    }
}

#endif
