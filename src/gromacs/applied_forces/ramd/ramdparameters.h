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

namespace gmx
{

/*! \internal
 * \brief Parameters for RAMD group
 */
struct RAMDGroup
{
    real force_;      ///< Force to be applied in kcal/mol/Angstrom
    real max_dist_;   ///< Specifies the distance in Angstrom between the COMs of the ligand
                      ///  and the receptor when the simulation is stopped
    real r_min_dist_; ///< Specifies the minimum distance in Angstrom
                      ///  to be traveled by the ligand in one RAMD step
};

/*! \internal
 * \brief Holding all directly user-provided parameters for Random Acceleration Molecular Dynamics (RAMD)
 *
 * Also used for setting all default parameters.
 */
struct RAMDParameters
{
    //! Indicate if density fitting is active
    bool active_ = false;

    int64_t                seed_;      ///< Initialization number for pseudo random number generator
    int                    ngroup_;    ///< Number of RAMD groups
    std::vector<RAMDGroup> group_;     ///< List of RAMD receptor-ligand pairs
    int                    eval_freq_; ///< Number of MD steps in one RAMD step
    int         force_out_freq_;       ///< Every 'force_out_freq' steps detailed output of forces will be written
    gmx_bool    old_angle_dist_;       ///< Use old angle distribution
    gmx_bool    connected_ligands_;    ///< Behavior of re-entering ligands into the dissociation radius
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_RAMDPARAMETERS_H
