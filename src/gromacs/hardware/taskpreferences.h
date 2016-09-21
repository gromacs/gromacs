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
#ifndef GMX_HARDWARE_TASKPREFERENCES_H
#define GMX_HARDWARE_TASKPREFERENCES_H

/*! \libinternal \file
 * \brief
 * This file declares class TaskPreferences,
 * which stores the user input and program decision
 * on assigning GPU-aware tasks to hardware resources.
 *
 * \author Aleksei Iupinov <a.yupinov@gmail.com>
 * \ingroup module_hardware
 */

#include <algorithm>
#include "gromacs/utility/exceptions.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"

//! GPU-aware tasks
enum class GpuTask
{
    NB,
    PME
};

/*! \brief
 * Options for ssignment to hardware.
 * This should have the same order and meaningful element count as nbpu_opt definition in mdrun.cpp.
 */
enum class DevicePreference : int
{
    Auto = 1,  // Up to scheduler, user trusts the program
    Cpu,
    Gpu,
    Hybrid // Using GPU locally, CPU non-locally, only implemented for NB
};

/*! \brief
 * Information about tasks' preferences concerning the hardware to run on
 */
class TaskPreferences
{
    public:
        bool definitelyUseGpu(GpuTask task)
        {
            return (get(task) == DevicePreference::Gpu) || (get(task) == DevicePreference::Hybrid);
        }

        bool maybeUseGpu(GpuTask task)
        {
            return definitelyUseGpu(task) || (get(task) == DevicePreference::Auto);
        }

        bool definitelyUseGpuForSomething()
        {
            bool r = false;
            for (const auto &ref : preferences)
            {
                r = r || definitelyUseGpu(ref.first);
            }
            return r;
        }

        bool maybeUseGpuForSomething()
        {
            bool r = false;
            for (const auto &ref : preferences)
            {
                r = r || maybeUseGpu(ref.first);
            }
            return r;
        }

        void set(GpuTask task, DevicePreference preference)
        {
            preferences[task] = preference;
        }

        DevicePreference get(GpuTask task)
        {
            return preferences.at(task); // can throw!
        }


        void finalize(t_commrec *cr)
        {
            bool undecided = false;
            for (const auto &ref : preferences)
            {
                undecided = undecided || (ref.second == DevicePreference::Auto);
            }
            GMX_RELEASE_ASSERT(!undecided, "The GPU task should have been assigned to a hardware resource");

            if (PAR(cr))
            {
                for (const auto &ref : preferences)
                {
                    GpuTask          temp1 = ref.first;
                    DevicePreference temp2 = ref.second;
                    gmx_bcast_sim(sizeof(temp1), &temp1, cr);
                    gmx_bcast_sim(sizeof(temp2), &temp2, cr);
                    set(temp1, temp2);
                }
            }
        }

    private:
        std::map<GpuTask, DevicePreference> preferences;
};

#endif
