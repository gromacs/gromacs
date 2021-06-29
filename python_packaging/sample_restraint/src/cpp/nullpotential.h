/*! \file
 * \brief Provide a minimal pluggable restraint potential for illustration and testing purposes.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#ifndef GMXAPI_EXTENSION_NULLPOTENTIAL_H
#define GMXAPI_EXTENSION_NULLPOTENTIAL_H


#include <array>
#include <memory>
#include <mutex>
#include <vector>

#include "gmxapi/gromacsfwd.h"
#include "gmxapi/session.h"
#include "gmxapi/md/mdmodule.h"

#include "gromacs/restraint/restraintpotential.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "sessionresources.h"


// Ultimately, the shared object library for the Python module should not export any symbols,
// but we use a generic namespace here for tidiness.
namespace plugin
{

struct null_input_param_type
{
    std::vector<int> sites_;
    int              count_;
};

//! Creation function for NullRestraint input and state data.
null_input_param_type makeNullParams(std::vector<int>&& sites);

//! Support the IRestraintPotential protocol.
std::vector<int> sites(const null_input_param_type& input);

//! Implement the NullRestraint force-provider.
gmx::PotentialPointData evaluate(gmx::Vector r1, gmx::Vector r2, double t, null_input_param_type* input);

//! Get the number of times evaluate has been called.
int count(const null_input_param_type& input);

class NullRestraint : public ::gmx::IRestraintPotential
{
public:
    using input_param_type = null_input_param_type;

    NullRestraint(std::vector<int>                            sites,
                  const input_param_type&                     params,
                  std::shared_ptr<Resources> resources);

    ~NullRestraint() override = default;
    [[nodiscard]] std::vector<int> sites() const override;
    gmx::PotentialPointData        evaluate(gmx::Vector r1, gmx::Vector r2, double t) override;

    input_param_type data_;
};

// Important: Just declare the template instantiation here for client code.
// We will explicitly instantiate a definition in the .cpp file where the input_param_type is defined.
extern template class RestraintModule<NullRestraint>;

} // end namespace plugin

#endif // GMXAPI_EXTENSION_NULLPOTENTIAL_H
