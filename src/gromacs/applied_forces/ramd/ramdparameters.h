/*! \internal \file
 * \brief
 * Declares RAMD parameters
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDPARAMETERS_H
#define GMX_APPLIED_FORCES_RAMDPARAMETERS_H

#include <cstdint>

#include <string>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \internal
 * \brief Parameters for RAMD group
 */
struct RAMDGroup
{
    //! Force to be applied in kcal/mol/Angstrom
    real force_ = 600.0;

    std::string        receptor_         = "Protein";
    std::vector<Index> receptor_indices_ = {};
    int                receptor_pbcatom_ = 0;

    std::string        ligand_         = "Ligand";
    std::vector<Index> ligand_indices_ = {};
    int                ligand_pbcatom_ = 0;

    //! Specifies the distance in Angstrom between the COMs of the ligand
    //! and the receptor when the simulation is stopped
    real max_dist_ = 4.0;

    //! Specifies the minimum distance in Angstrom
    //! to be traveled by the ligand in one RAMD step
    real r_min_dist_ = 0.0025;
};

/*! \internal
 * \brief Holding all directly user-provided parameters for Random Acceleration Molecular Dynamics (RAMD)
 *
 * Also used for setting all default parameters.
 */
struct RAMDParameters
{
    //! Indicate if RAMD is active
    bool active_ = false;

    //! Initialization number for pseudo random number generator
    int64_t seed_ = 1234;

    //! Interval for evaluating the COM distance and possibly changing the force direction
    int eval_freq_ = 50;

    //! Number of RAMD groups
    int ngroups_ = 0;

    //! List of RAMD receptor-ligand pairs
    std::vector<RAMDGroup> groups_ = {};

    //! Interval for writing out the COM distances of all RAMD groups
    int out_freq_ = 100;

    //! Use previous step COM as PBC reference
    bool pbc_ref_prev_step_com_ = false;

    //! Use old angle distribution
    bool old_angle_dist_ = false;

    //! Behavior of re-entering ligands into the dissociation radius
    bool connected_ligands_ = false;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_RAMDPARAMETERS_H
