/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_TRAJECTORYANALYSIS_SURFACEAREA_H
#define GMX_TRAJECTORYANALYSIS_SURFACEAREA_H

#include "gromacs/legacyheaders/types/simple.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"

struct t_pbc;

#define FLAG_DOTS       01
#define FLAG_VOLUME     02
#define FLAG_ATOM_AREA  04

namespace gmx
{

/*! \internal
 * \brief
 * Computes surface areas for a group of atoms/spheres.
 *
 * This class provides a surface area/volume calculator.
 *
 * The algorithm is based on representing each atom/sphere surface as a set of
 * dots, and determining which dots are on the surface (not covered by any
 * other atom/sphere).  The dots are distributed evenly using an icosahedron- or
 * a dodecahedron-based method (see the original reference cited in the code).
 * The area is then estimated from the area represented by each dot.
 * The volume is calculated by selecting a fixed point and integrating over the
 * surface dots, summing up the cones whose tip is at the fixed point and base
 * at the surface points.
 *
 * The default dot density per sphere is 32, which gives quite inaccurate
 * areas and volumes, but a reasonable number of surface points.  According to
 * original documentation of the method, a density of 600-700 dots gives an
 * accuracy of 1.5 A^2 per atom.
 *
 * \ingroup module_trajectoryanalysis
 */
class SurfaceAreaCalculator
{
    public:
        /*! \brief
         * Initializes a surface area calculator.
         *
         * \throws std::bad_alloc if out of memory.
         */
        SurfaceAreaCalculator();
        ~SurfaceAreaCalculator();

        /*! \brief
         * Sets the number of surface dots per sphere to use.
         *
         * This function must be called before calculate() to set the desired
         * accuracy/computational cost.
         */
        void setDotCount(int dotCount);
        /*! \brief
         * Sets the radii of spheres to use in the calculation.
         *
         * \param[in]  radius  Radius for each atom/sphere.
         *
         * This function must be called before calculate() to set the radii for
         * the spheres.  All calculations must use the same set of radii to
         * share the same grid search.
         * These radii are used as-is, without adding any probe radius.
         * The passed array must remain valid for the lifetime of this object.
         *
         * Does not throw.
         */
        void setRadii(const ConstArrayRef<real> &radius);

        /*! \brief
         * Requests calculation of volume.
         *
         * If not called, and FLAG_VOLUME is not passed to calculate(), the
         * volume output is not produced.
         *
         * Does not throw.
         */
        void setCalculateVolume(bool bVolume);
        /*! \brief
         * Requests output of per-atom areas.
         *
         * If not called, and FLAG_ATOM_AREA is not passed to calculate(), the
         * atom area output is not produced.
         *
         * Does not throw.
         */
        void setCalculateAtomArea(bool bAtomArea);
        /*! \brief
         * Requests output of all surface dots.
         *
         * If not called, and FLAG_DOTS is not passed to calculate(), the
         * surface dot output is not produced.
         *
         * Does not throw.
         */
        void setCalculateSurfaceDots(bool bDots);

        /*! \brief
         * Calculates the surface area for a set of positions.
         *
         * \param[in]  x       Atom positions (sphere centers).
         * \param[in]  pbc     PBC information (if `NULL`, calculation is done
         *     without PBC).
         * \param[in]  nat     Number of atoms to calculate.
         * \param[in]  index   Atom indices to include in the calculation.
         * \param[in]  flags   Additional flags for the calculation.
         * \param[out] area    Total surface area (must be non-`NULL`).
         * \param[out] volume  Total volume (can be `NULL`).
         * \param[out] at_area Surface area for each atom in \p index
         *     (\p nat values) (can be `NULL`).
         * \param[out] lidots  Surface dots as x,y,z triplets (`3*lidots` values)
         *     (can be `NULL`).
         * \param[out] n_dots Number of surface dots in \p lidots
         *     (can be `NULL`).
         *
         * Calculates the surface area of spheres centered at `x[index[0]]`,
         * ..., `x[index[nat-1]]`, with radii `radii[index[0]]`, ..., where
         * `radii` is the array passed to setRadii().
         *
         * If \p flags is 0, the calculation is done for the items specified
         * with setCalculateVolume(), setCalculateAtomArea(), and
         * setCalculateSurfaceDots().  Flags can specify FLAG_VOLUME,
         * FLAG_ATOM_AREA, and/or FLAG_DOTS to request additional output for
         * this particular calculation.  If any output is `NULL`, that output
         * is not calculated, irrespective of the calculation mode set.
         *
         * \todo
         * Make the output options more C++-like, in particular for the array
         * outputs.
         */
        void calculate(const rvec *x, const t_pbc *pbc,
                       int nat, atom_id index[], int flags, real *area,
                       real *volume, real **at_area,
                       real **lidots, int *n_dots) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
