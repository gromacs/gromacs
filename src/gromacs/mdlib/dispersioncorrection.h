/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
#ifndef GMX_MDLIB_DISPERSIONCORRECTION_H
#define GMX_MDLIB_DISPERSIONCORRECTION_H

#include <cstdio>

#include <array>
#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct interaction_const_t;
struct t_forcetable;
struct t_inputrec;
enum class DispersionCorrectionType : int;
enum class VanDerWaalsType : int;
enum class FreeEnergyPerturbationType : int;
namespace gmx
{
class MDLogger;
} // namespace gmx

class DispersionCorrection
{
public:
    /*! \brief Constructor
     *
     * \param[in] mtop           The global topology
     * \param[in] inputrec       The input record
     * \param[in] useBuckingham  True when Buckingham is used instead of LJ
     * \param[in] ic             The nonbonded interaction parameters
     * \param[in] tableFileName  Table file name, should != nullptr (checked)
     */
    DispersionCorrection(const gmx_mtop_t&          mtop,
                         const t_inputrec&          inputrec,
                         bool                       useBuckingham,
                         const interaction_const_t& ic,
                         const char*                tableFileName);

    /*! \brief Print dispersion correction information to the log file
     *
     * \param[in] mdlog  The MD logger
     */
    void print(const gmx::MDLogger& mdlog) const;

    /*! \brief Computes and sets energy and virial correction parameters
     *
     * Sets all parameters that are affected by the cut-off and/or the
     * LJ-Ewald coefficient. Should be called before calling calculate()
     * and whenever interaction settings change, e.g. PME tuning.
     *
     * \param[in] ic  The nonbonded interaction parameters
     */
    void setParameters(const interaction_const_t& ic);

    /*! \internal
     * \brief Struct for returning all dispersion correction quantities
     */
    struct Correction
    {
        /*! \brief Correct the virial tensor for the missing dispersion
         *
         * \param[in,out] virialTensor  The virial tensor to correct
         */
        void correctVirial(tensor virialTensor) const
        {
            for (int m = 0; m < DIM; m++)
            {
                virialTensor[m][m] += virial;
            }
        }

        real virial   = 0; //!< Scalar correction to the virial
        real pressure = 0; //!< Scalar correction to the pressure
        real energy   = 0; //!< Correction to the energy
        real dvdl     = 0; //!< Correction to dH/dlambda
    };

    /*! \brief Computes and returns the dispersion correction for the pressure and energy
     *
     * \param[in]  box       The simulation unit cell
     * \param[in]  lambda    The free-energy coupling parameter
     */
    Correction calculate(const matrix box, real lambda) const;

private:
    /*! \internal \brief Parameters that depend on the topology only
     */
    class TopologyParams
    {
    public:
        TopologyParams(const gmx_mtop_t& mtop, const t_inputrec& inputrec, bool useBuckingham);

        //! The number of atoms for computing the atom density
        int numAtomsForDensity_;
        //! The number of interactions to correct for, usually num. atoms/2
        real numCorrections_;
        //! Average C6 coefficient for for topology A/B ([0]/[1])
        std::array<real, 2> avcsix_;
        //! Average C12 coefficient for for topology A/B ([0]/[1])
        std::array<real, 2> avctwelve_;
    };

    /*! \internal \brief Parameters that depend on the interaction functions and topology
     */
    struct InteractionParams
    {
    public:
        ~InteractionParams();

        //! Table used for correcting modified LJ interactions
        std::unique_ptr<t_forcetable> dispersionCorrectionTable_;

        //! Dispersion energy shift constant
        real enershiftsix_ = 0;
        //! Repulsion energy shift constant
        real enershifttwelve_ = 0;
        //! Dispersion energy difference per atom per unit of volume
        real enerdiffsix_ = 0;
        //! Repulsion energy difference per atom per unit of volume
        real enerdifftwelve_ = 0;
        //! Dispersion virial difference per atom per unit of volume
        real virdiffsix_ = 0;
        //! Repulsion virial difference per atom per unit of volume
        real virdifftwelve_ = 0;
    };

    //! Sets the interaction parameters
    static void setInteractionParameters(DispersionCorrection::InteractionParams* iParams,
                                         const interaction_const_t&               ic,
                                         const char*                              tableFileName);

    //! Returns whether we correct both dispersion and repulsion
    bool correctFullInteraction() const;

    //! Type of dispersion correction
    DispersionCorrectionType eDispCorr_;
    //! Type of Van der Waals interaction
    VanDerWaalsType vdwType_;
    //! Free-energy perturbation
    FreeEnergyPerturbationType eFep_;
    //! Topology parameters
    TopologyParams topParams_;
    //! Interaction parameters
    InteractionParams iParams_;
};

#endif
