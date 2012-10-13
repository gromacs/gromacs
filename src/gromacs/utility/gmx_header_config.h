/*! \file
 * \brief
 * Include file for configuration macros that affect installed headers.
 *
 * This include file (or rather, one that it includes) will configured by CMake
 * and installed with GROMACS header files so that they can refer to a central
 * location for #defines that will be available for builds of projects that
 * depend on GROMACS.
 *
 * The actual defines are in gmx_header_config_gen.h to allow usage of relative
 * include paths before installation.
 *
 * \todo
 * It would be better to have the defines here such that they are not generated
 * from CMake, but instead detected using #ifdefs (possible for some of the
 * macros currently used).
 * Even better would be to not have these defines at all.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#include "gmx_header_config_gen.h"
