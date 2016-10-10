/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Describes tester functions and common routines for PME tests.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_ewald
 */
#ifndef GMX_EWALD_PME_BASE_TESTER_H
#define GMX_EWALD_PME_BASE_TESTER_H

#include <functional>
#include <memory>

#include "gromacs/mdrunutility/mdmodules.h"

// Testers

//! A typedef for a function being tested
typedef std::function<void()> TestFunc;

//! A tester function should accept the TestFunc reference and work with it
typedef void (*TesterFuncPtr)(const TestFunc &);
//! A global pointer to the death tester
extern TesterFuncPtr g_TestForFatalErrorPtr;
//! A global pointer to the success tester
extern TesterFuncPtr g_TestForSuccessPtr;

// PME-specific stuff

struct gmx_pme_t;
struct t_inputrec;

//! Environment for getting the t_inputrec structure easily
static gmx::MDModules s_mdModules;

/*! \brief A safe pointer type for PME. */
typedef std::unique_ptr<gmx_pme_t, void(*)(gmx_pme_t *)> pmeSafePointer;

//! Simple PME initialization based on input
pmeSafePointer PmeDummyInit(const t_inputrec *inputRec);


#endif // GMX_EWALD_PME_BASE_TESTER_H
