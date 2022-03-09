/*! \file
 * \brief Code to implement the potential declared in nullpotential.h
 *
 * This restraint is implemented in a transitional style. We are moving in the direction of
 * callback based data flow. There is also a preference amongst the GROMACS developers for
 * stateless objects or free functions. State can be provided by library-managed facilities
 * rather than stored in long-lived objects.
 *
 * The IRestraintPotential framework is due for revision in conjunction with ongoing evolution of
 * the gmx::MDModules interactions. Until then, we try to use a forward looking design.
 *
 * \author M. Eric Irrgang <ericirrgang@gmail.com>
 */

#include "nullpotential.h"

#include <utility>

namespace plugin
{

null_input_param_type makeNullParams(std::vector<int>&& sites)
{
    return null_input_param_type{ std::move(sites), 0 };
}

std::vector<int> sites(const null_input_param_type& input)
{
    return input.sites_;
}

gmx::PotentialPointData evaluate(gmx::Vector /*r1*/, gmx::Vector /*r2*/, double /*t*/, null_input_param_type* input)
{
    ++input->count_;
    return gmx::PotentialPointData();
}

int count(const null_input_param_type& input)
{
    return input.count_;
}

gmx::PotentialPointData NullRestraint::evaluate(gmx::Vector r1, gmx::Vector r2, double t)
{
    return ::plugin::evaluate(r1, r2, t, &data_);
}

std::vector<int> NullRestraint::sites() const
{
    return ::plugin::sites(data_);
}

NullRestraint::NullRestraint(std::vector<int>                       sites,
                             const NullRestraint::input_param_type& params,
                             std::shared_ptr<Resources> /*resources*/) :
    data_{ std::move(sites), params.count_ }
{
}

// Important: Explicitly instantiate a definition for the templated class declared in
// ensemblepotential.h. Failing to do this will cause a linker error.
template class ::plugin::RestraintModule<NullRestraint>;

} // end namespace plugin
