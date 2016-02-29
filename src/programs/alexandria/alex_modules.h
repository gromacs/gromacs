/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GMX_ALEXMODULES_H
#define GMX_ALEXMODULES_H

namespace gmx
{
class CommandLineModuleManager;
} // namespace gmx

/*! \internal \brief
 * Registers all alex command-line modules.
 *
 * \param[in] manager  Command-line module manager to receive the modules.
 * \throws    std::bad_alloc if out of memory.
 *
 * Registers all modules corresponding to pre-5.0 binaries such that
 * they can be run through \p manager.
 */
void registerAlexandriaModules(gmx::CommandLineModuleManager *manager);

#endif
