/*! \internal \file
 * \brief
 * Declares options for RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_RAMDOPTIONS_H
#define GMX_APPLIED_FORCES_RAMDOPTIONS_H

#include <string>

#include "gromacs/mdtypes/imdpoptionprovider.h"

#include "ramdparameters.h"

namespace gmx
{

class EnergyCalculationFrequencyErrors;
class IndexGroupsAndNames;
class KeyValueTreeObject;
class KeyValueTreeBuilder;

/*! \internal
 * \brief Input data storage for RAMD
 */
class RAMDOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for RAMD to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     * \note This should be symmetrical to option initialization without
     *       employing manual prefixing with the section name string once
     *       the legacy code blocking this design is removed.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Report if this set of options is active
    bool active() const;

    //! Get parameters_ instance
    const RAMDParameters& parameters();

private:

    //! Parameter values for force & energy evaluation
    RAMDParameters parameters_;

};

} // namespace gmx

#endif
