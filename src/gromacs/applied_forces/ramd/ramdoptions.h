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

class IndexGroupsAndNames;
class KeyValueTreeObject;
class KeyValueTreeBuilder;
class MDLogger;

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

    //! Set the MDLogger instance
    void setLogger(const MDLogger& logger);

    //! Get the logger instance
    const MDLogger& logger() const;

    //! Store the parameters that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set atom groups
    void setInputGroupIndices(const IndexGroupsAndNames&);

private:

    //! Parameter values for force & energy evaluation
    RAMDParameters parameters_;

    /*! \brief MDLogger during preprocessing
     *
     * This is a pointer only because we need an "optional reference"
     * to a const MDLogger before the notification always provides the
     * actual reference. */
    const MDLogger* logger_ = nullptr;

};

} // namespace gmx

#endif
