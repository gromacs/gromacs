/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

 #ifndef GMX_APPLIEDFORCES_DENSITYFITTING_LOCALDENSFITDATA_H_
 #define GMX_APPLIEDFORCES_DENSITYFITTING_LOCALDENSFITDATA_H_

 #include "gromacs/domdec/localatomset.h"
 #include "gromacs/domdec/localatomsetmanager.h"

 #include "densfitdata.h"
 #include "gromacs/utility/alignedallocator.h"
 #include "gromacs/math/griddata/operations/densityspreader.h"
 #include "gromacs/math/paddedvector.h"

#include <vector>
namespace gmx
{

//! \brief local density fitting data
class LocalDensfitData
{
    public:
        LocalDensfitData(LocalAtomSet localAtomSet, const DensfitData &parameters, const std::vector<RVec> &x_densfit_whole);
        ~LocalDensfitData() = default;
        void updateTimeDependentData(real time, const DensfitData &parameters);

        void do_forces();

        const GridDataReal3D &simulatedMap() const;

        std::string infoString() const;

        void add_forces(rvec * f) const;

        void triggerShiftUpdate();

        void communicate(PaddedArrayRef<RVec> x, const t_commrec * cr, matrix box);

        void spreadAtoms(const t_commrec * cr, bool normalize);

        void sumSimulatedMap(const t_commrec * cr);
        template <class T> void densityDensityDerivative(const GridDataReal3D &referenceMap);
        void densityDensityDerivative(DensityPotential potential, const GridDataReal3D &referenceMap);

        template <class T> void goodnessOfFit(const GridDataReal3D &referenceMap);
        void goodnessOfFit(DensityPotential potential, const GridDataReal3D &referenceMap);
    private:
        real forceConstant_;                        /**< Current value k of strength constant of map
                                                       potential V_fit = k*(1 - cc) with cc being the
                                                       correlation coefficient                        */
        std::vector<real> weights;                  /**< weight of the local atoms for spreadingn and forces */
        real              sigma;                    /**< Current value of Gaussian width used for
                                                         spreading atomic positions on the map          */
        real              goodnessOfFit_;           /**< Correlation coefficient of the two maps        */
        std::vector<RVec> x_assembled;              /**< Just the positions that are spread             */
        std::vector<RVec> f_loc;                    /**< Forces on the local atoms of the
                                                       density fitting group                          */
        std::vector<RVec> xWholeMoleculeReference_; /**< x_assembled from last step                     */
        std::vector<IVec> x_shifts;                 /**< To make the molecule whole                     */
        std::vector<IVec> extra_shifts;             /**< Shifts added since last NS                     */
        gmx_bool          bUpdateShifts;            /**< Do we have to calculate new shift vectors
                                                         from x_assembled and x_old?                    */
        int               voxrange;                 /**< Max. number of voxels to be computed for a
                                                         single atom in a single dimension x, y, or z   */
        int               vox_per_atom;             /**< Total number of voxels for one atom (x,y, + z) */
        LocalAtomSet      localAtoms;

        GridDataReal3D    map_sim_;       /* The map simulated from atomic position data     */
        GridDataReal3D    derivativeMap_; /* The derivative of goodness of fit of simulated Map to reference map */

        DensitySpreader   densitySpreader_;

        /* The following two temporary vectors (one for each OpenMP thread) store
         * erf values around a single atoms, thus we can compute them all in one go,
         * and with SIMD acceleration */
        using AlignedRealVector = std::vector < real, gmx::AlignedAllocator < real>>;
        std::vector<AlignedRealVector> erfVector;  /**< vector of vectors of erf values */
        std::vector<AlignedRealVector> expVector;  /**< same for exp values             */
};
}


 #endif /* end of include guard: GMX_APPLIEDFORCES_DENSITYFITTING_LOCALDENSFITDATA_H_ */
