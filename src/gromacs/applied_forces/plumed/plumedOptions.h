/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Declares options for PLUMED. This class handles parameters set during
 * pre-processing time.
 *
 * \author Daniele Rapetti <drapetti@sissa.it>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_PLUMEDOPTIONPROVIDER_H
#define GMX_APPLIED_FORCES_PLUMEDOPTIONPROVIDER_H

#include <optional>
#include <string>

#include "gromacs/utility/real.h"

struct gmx_mtop_t;
struct t_commrec;

namespace gmx
{
enum class StartingBehavior;
struct EnsembleTemperature;

struct PlumedOptions
{
    std::string         plumedFile_;
    int                 natoms_;
    const t_commrec*    cr_;
    real                simulationTimeStep_;
    std::optional<real> ensembleTemperature_{};
    StartingBehavior    startingBehavior_{};
    bool                active_{ false };
};

class PlumedOptionProvider
{
public:
    /*! @brief Sets the needed informations from the topology object
     *
     * As now oly hte number of atoms is fetched
     * @param mtop topology object
     */
    void setTopology(const gmx_mtop_t& mtop);
    /*! @brief Sets the (eventual) ensemble temperature
     *  @param temp the object with the optional temperature
     */
    void setEnsembleTemperature(const EnsembleTemperature& temp);
    /*! @brief Sets the name of the PLUMED file to read
     *
     * When called, with a non empty string, it activates the PLUMED module
     * this simulation.
     *  @param  fname the (optional) name of the file
     */
    void setPlumedFile(const std::optional<std::string>& fname);
    /*! @brief Sets the timestep
     * @param timeStep the timestep value
     */
    void setSimulationTimeStep(double timeStep);
    /*! @brief Sets the starting beahviour of the simulation
     * @param startingBehavior the starting behaviopur object
     */
    void setStartingBehavior(const StartingBehavior& startingBehavior);
    /*! @brief Sets the address to the communication record object
     * @param cr  the Communication Record object
     */
    void setComm(const t_commrec& cr);
    //! @brief returns the active status of the module
    bool active() const;

    //! @brief returns a reference to the internal PlumedOptions element
    const PlumedOptions& options() const;

private:
    PlumedOptions opts_;
};
} // namespace gmx
#endif // GMX_APPLIED_FORCES_PLUMEDOPTIONPROVIDER_H
