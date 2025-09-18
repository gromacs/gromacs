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

    //! Process input options to parameters, including input file reading.
    const RAMDParameters& buildParameters();

    /*! \brief Evaluate and store atom indices.
     *
     * During pre-processing, use the group string from the options to
     * evaluate the indices of the atoms to be subject to forces from this
     * module.
     */
    void setFitGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames);

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    //! Check if input parameters are consistent with other simulation parameters
    void checkEnergyCaluclationFrequency(EnergyCalculationFrequencyErrors* energyCalculationFrequencyErrors) const;

private:

    //! Parameter values for force & energy evaluation
    RAMDParameters parameters_;

};

} // namespace gmx

#endif
