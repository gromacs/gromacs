#ifndef GMX_COMMANDLINE_CMDLINEMODULEMANAGER_IMPL_H
#define GMX_COMMANDLINE_CMDLINEMODULEMANAGER_IMPL_H

#include <map>
#include <string>
#include <vector>

#include "cmdlinemodule.h"
#include "cmdlinemodulemanager.h"

#include "gromacs/utility/common.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{
//! Container type for mapping module names to module objects.
typedef std::map<std::string, CommandLineModulePointer> CommandLineModuleMap;

/*! \internal
 * \brief
 * Internal data for a CommandLineModuleManager module group.
 *
 * This class contains the state of a module group.  CommandLineModuleGroup
 * provides the public interface to construct/alter the state, and
 * CommandLineModuleManager and its associated classes use it for help output.
 *
 * \ingroup module_commandline
 */
class CommandLineModuleGroupData
{
    public:
        /*! \brief
         * Shorthand for a list of modules contained in the group.
         *
         * The first element in the contained pair contains the tag
         * (gmx-something) of the module, and the second element contains the
         * description.  The second element is never NULL.
         */
        typedef std::vector<std::pair<std::string, const char *> > ModuleList;

        /*! \brief
         * Constructs an empty module group.
         *
         * \param[in] modules  List of all modules
         *     (used for checking and default descriptions).
         * \param[in] title    Title of the group.
         *
         * Does not throw.
         */
        CommandLineModuleGroupData(const CommandLineModuleMap &modules,
                                   const char                 *title)
            : allModules_(modules), title_(title)
        {
        }

        //! Returns the title for the group.
        const char *title() const { return title_; }
        //! Returns the list of modules in the group.
        const ModuleList &modules() const { return modules_; }

        /*! \brief
         * Adds a module to the group.
         *
         * \param[in] name        Name of the module.
         * \param[in] description Description of the module in this group.
         * \throws    std::bad_alloc if out of memory.
         *
         * If \p description is NULL, the description returned by the module is
         * used.
         */
        void addModule(const char *name, const char *description);

    private:
        const CommandLineModuleMap &allModules_;
        const char                 *title_;
        ModuleList                  modules_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineModuleGroupData);
};

//! Smart pointer type for managing a CommandLineModuleGroup.
typedef gmx_unique_ptr<CommandLineModuleGroupData>::type
    CommandLineModuleGroupDataPointer;
//! Container type for keeping a list of module groups.
typedef std::vector<CommandLineModuleGroupDataPointer>
    CommandLineModuleGroupList;

} // namespace gmx

#endif
