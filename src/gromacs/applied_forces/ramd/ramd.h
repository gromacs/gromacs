#pragma once

#include <memory>
#include <string>

namespace gmx
{

class IMDModule;

/*! \internal
 *  \brief Information about the RAMD module.
 *
 * Provides name and method to create a RAMD module.
 */
struct RAMDModuleInfo
{
    //! Creates a module for applying forces according to a RAMD.
    static std::unique_ptr<IMDModule> create();
    //! The name of the module
    static const std::string name_;
};

} // namespace gmx
