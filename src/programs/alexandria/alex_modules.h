/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
 
#ifndef GMX_ALEXMODULES_H
#define GMX_ALEXMODULES_H

int alex_gentop(int argc, char *argv[]);
int alex_tune_fc(int argc, char *argv[]);
int alex_tune_eem(int argc, char *argv[]);
int alex_tune_pol(int argc, char *argv[]);
int alex_poldata_test(int argc, char *argv[]);
int alex_fit_qs_zeta(int argc, char *argv[]);
int alex_gauss2molprop(int argc, char *argv[]);
int alex_bastat(int argc, char *argv[]);
int alex_analyze(int argc, char *argv[]);
int alex_gen_table(int argc, char *argv[]);
int alex_merge_mp(int argc, char *argv[]);
int alex_merge_pd(int argc, char *argv[]);
int alex_mp2csv(int argc, char *argv[]);
int alex_molprop_test(int argc, char *argv[]);
int alex_molprop_check(int argc, char *argv[]);
int alex_tune_zeta(int argc, char *argv[]);
int alex_molselect(int argc, char *argv[]);

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
