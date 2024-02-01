/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Implements the Colvars GROMACS proxy class
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#include "colvarproxygromacs.h"

#include <sstream>

#include "colvarproxy_gromacs_version.h"


namespace gmx
{

ColvarProxyGromacs::ColvarProxyGromacs(const std::string& colvarsConfigString,
                                       t_atoms            atoms,
                                       PbcType            pbcType,
                                       const MDLogger*    logger,
                                       bool               doParsing,
                                       const std::map<std::string, std::string>& inputStrings,
                                       real ensembleTemperature,
                                       int  seed) :
    gmxAtoms_(atoms), pbcType_(pbcType), logger_(logger), doParsing_(doParsing)
{
    engine_name_ = "GROMACS";
    version_int  = get_version_from_string(COLVARPROXY_VERSION);

    //! From colvarproxy
    //! The 5 variables below are defined in the `colvarproxy` base class

    // Retrieve masses and charges from input file
    updated_masses_ = updated_charges_ = true;

    // User-scripted forces are not available in GROMACS
    have_scripts = false;

    angstrom_value_ = 0.1;

    // From Gnu units
    // $ units -ts 'k' 'kJ/mol/K/avogadro'
    boltzmann_ = 0.0083144621;

    // Get the thermostat temperature.
    set_target_temperature(ensembleTemperature);

    // GROMACS random number generation.
    // Used the number defined in the mdp options for the seed.
    // -1 (default value) stands for random
    if (seed == -1)
    {
        rng_.seed(makeRandomSeed());
    }
    else
    {
        rng_.seed(seed);
    }


    // Read configuration file and set up the proxy during pre-processing
    // and during simulation phase but only on the master node.
    if (doParsing)
    {

        // Retrieve input files stored as string in the key-value-tree (KVT) of the TPR.
        // Add them to the map of colvars input data.
        for (const auto& [inputName, content] : inputStrings)
        {
            // input_streams_ defined in colvarproxy
            input_streams_[inputName] = new std::istringstream(content);
        }

        colvars = new colvarmodule(this);
        cvm::log(cvm::line_marker);
        cvm::log("Start colvars Initialization.");


        colvars->cite_feature("GROMACS engine");
        colvars->cite_feature("Colvars-GROMACS interface");

        if (cvm::debug())
        {
            cvm::log("Initializing the colvars proxy object.\n");
        }

        int errorCode = colvarproxy::setup();
        errorCode |= colvars->read_config_string(colvarsConfigString);
        errorCode |= colvars->update_engine_parameters();
        errorCode |= colvars->setup_input();

        if (errorCode != COLVARS_OK)
        {
            error("Error when initializing Colvars module.");
        }

        // Citation Reporter
        cvm::log(std::string("\n") + colvars->feature_report(0) + std::string("\n"));

        // TODO get initial step number from MDModules
        // colvars->set_initial_step(static_cast<cvm::step_number>(0L));
    }
}


cvm::real ColvarProxyGromacs::rand_gaussian()
{
    return normalDistribution_(rng_);
}

void ColvarProxyGromacs::log(std::string const& message)
{
    if (logger_)
    {
        std::istringstream is(message);
        std::string        line;
        while (std::getline(is, line))
        {
            GMX_LOG(logger_->info).appendText("colvars: " + line + "\n");
        }
    }
}

void ColvarProxyGromacs::error(std::string const& message)
{
    log(message);
    GMX_THROW(InternalError("Error in collective variables module.\n"));
}


int ColvarProxyGromacs::set_unit_system(std::string const& unitsIn, bool /*colvarsDefined*/)
{
    if (unitsIn != "gromacs")
    {
        cvm::error(
                "Specified unit system \"" + unitsIn
                + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
        return COLVARS_ERROR;
    }
    return COLVARS_OK;
}


// **************** ATOMS ****************

int ColvarProxyGromacs::check_atom_id(int atomNumber)
{
    // GROMACS uses zero-based arrays.
    int const aid = (atomNumber - 1);

    if (cvm::debug())
    {
        log("Adding atom " + cvm::to_str(atomNumber) + " for collective variables calculation.\n");
    }
    if ((aid < 0) || (aid >= gmxAtoms_.nr))
    {
        cvm::error("Error: invalid atom number specified, " + cvm::to_str(atomNumber) + "\n",
                   COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
    }

    return aid;
}


int ColvarProxyGromacs::init_atom(int atomNumber)
{
    // GROMACS uses zero-based arrays.
    int aid = atomNumber - 1;

    // atoms_ids & atoms_refcount declared in `colvarproxy_atoms` class
    for (size_t i = 0; i < atoms_ids.size(); i++)
    {
        if (atoms_ids[i] == aid)
        {
            // this atom id was already recorded
            atoms_refcount[i] += 1;
            return i;
        }
    }

    aid = check_atom_id(atomNumber);

    if (aid < 0)
    {
        return COLVARS_INPUT_ERROR;
    }

    int const index = add_atom_slot(aid);
    updateAtomProperties(index);
    return index;
}

void ColvarProxyGromacs::updateAtomProperties(int index)
{

    // update mass
    double const mass = gmxAtoms_.atom[atoms_ids[index]].m;
    if (mass <= 0.001)
    {
        this->log("Warning: near-zero mass for atom " + cvm::to_str(atoms_ids[index] + 1)
                  + "; expect unstable dynamics if you apply forces to it.\n");
    }

    // atoms_masses & atoms_charges declared in `colvarproxy_atoms` class
    atoms_masses[index]  = mass;
    atoms_charges[index] = gmxAtoms_.atom[atoms_ids[index]].q;
}

ColvarProxyGromacs::~ColvarProxyGromacs()
{
    if (colvars != nullptr)
    {
        delete colvars;
        colvars = nullptr;
    }
}


cvm::rvector ColvarProxyGromacs::position_distance(cvm::atom_pos const& pos1, cvm::atom_pos const& pos2) const
{
    rvec r1, r2, dr;
    r1[0] = pos1.x;
    r1[1] = pos1.y;
    r1[2] = pos1.z;
    r2[0] = pos2.x;
    r2[1] = pos2.y;
    r2[2] = pos2.z;

    pbc_dx(&gmxPbc_, r2, r1, dr);
    return cvm::atom_pos(dr[0], dr[1], dr[2]);
}


} // namespace gmx
