/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#ifndef GMX_HARDWARE_HWINFO_H
#define GMX_HARDWARE_HWINFO_H

#include <string>
#include <vector>

#include "gromacs/hardware/gpu_hw_info.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
namespace gmx
{
class CpuInfo;
class HardwareTopology;
} // namespace

/* Hardware information structure with CPU and GPU information.
 * It is initialized by gmx_detect_hardware().
 * NOTE: this structure may only contain structures that are globally valid
 *       (i.e. must be able to be shared among all threads) */
struct gmx_hw_info_t
{
    /* Data for our local physical node */
    struct gmx_gpu_info_t gpu_info;                /* Information about GPUs detected in the system */

    int                   nthreads_hw_avail;       /* Number of hardware threads available; this number
                                                      is based on the number of CPUs reported as available
                                                      by the OS at the time of detection. */

    const gmx::CpuInfo *         cpuInfo;          /* Information about CPU capabilities */
    const gmx::HardwareTopology *hardwareTopology; /* Information about hardware topology */

    /* Data reduced through MPI over all physical nodes */
    int                 nphysicalnode;       /* Number of physical nodes */
    int                 ncore_tot;           /* Sum of #cores over all nodes, can be 0 */
    int                 ncore_min;           /* Min #cores over all nodes */
    int                 ncore_max;           /* Max #cores over all nodes */
    int                 nhwthread_tot;       /* Sum of #hwthreads over all nodes */
    int                 nhwthread_min;       /* Min #hwthreads over all nodes */
    int                 nhwthread_max;       /* Max #hwthreads over all nodes */
    int                 ngpu_compatible_tot; /* Sum of #GPUs over all nodes */
    int                 ngpu_compatible_min; /* Min #GPUs over all nodes */
    int                 ngpu_compatible_max; /* Max #GPUs over all nodes */

    int                 simd_suggest_min;    /* Highest SIMD instruction set supported by all ranks */
    int                 simd_suggest_max;    /* Highest SIMD instruction set supported by at least one rank */

    gmx_bool            bIdenticalGPUs;      /* TRUE if all ranks have the same type(s) and order of GPUs */
    bool                haveAmdZenCpu;       /* TRUE when at least one CPU in any of the nodes is AMD Zen */
};


/* The options for the thread affinity setting, default: auto */
enum {
    threadaffSEL, threadaffAUTO, threadaffON, threadaffOFF, threadaffNR
};

/*! \internal \brief Threading and GPU options, can be set automatically or by the user
 *
 * \todo During mdrunner(), if the user has left any of these values
 * at their defaults (which tends to mean "choose automatically"),
 * then those values are over-written with the result of such
 * automation. This creates problems for the subsequent code in
 * knowing what was done, why, and reporting correctly to the
 * user. Find a way to improve this.
 */
struct gmx_hw_opt_t
{
    //! Total number of threads requested (thread-MPI + OpenMP).
    int           nthreads_tot = 0;
    //! Number of thread-MPI threads requested.
    int           nthreads_tmpi = 0;
    //! Number of OpenMP threads requested.
    int           nthreads_omp = 0;
    //! Number of OpenMP threads to use on PME_only ranks.
    int           nthreads_omp_pme = 0;
    //! Thread affinity switch, see enum above.
    int           thread_affinity = threadaffSEL;
    //! Logical core pinning stride.
    int           core_pinning_stride = 0;
    //! Logical core pinning offset.
    int           core_pinning_offset = 0;
    //! Empty, or a string provided by the user declaring (unique) GPU IDs available for mdrun to use.
    std::string   gpuIdsAvailable = "";
    //! Empty, or a string provided by the user mapping GPU tasks to devices.
    std::string   userGpuTaskAssignment = "";
    //! Tells whether mdrun is free to choose the total number of threads (by choosing the number of OpenMP and/or thread-MPI threads).
    bool          totNumThreadsIsAuto;
};

/*! \brief Contains and manages changes to hardware options.
 *
 * Keeps track of whether the user has set hardware options at the command line, and
 * disallows modifications if so. The modification state of the options can be queried,
 * which is clearer than checking if an option is set to 0.
 */
class hardwareOptionsManager
{
    public:
        /*! \brief Container class for a single hardware option
         *
         * At construction, saves whether or not a value was explicitly
         * set by the user. Continually keeps track of whether or not the
         * option has been set (by the user or the program). Disallows any
         * attempts to modify a user-set parameter
         */
        template <typename T> class hardwareOption
        {
            public:
                //! No construction without a reference.
                hardwareOption() = delete;
                //! Compare to default value - not set by user if equal
                hardwareOption<T>(T value, T defaultValue);

                //! Access parameter value
                T operator()() const {return value_; }
                //! Check if user has set value explicitly
                bool isSetByUser() const {return isSetByUser_; }
                //! Check is value is currently set
                bool isSet() const {return isSet_; }
                //! Attempt to set value. Will update isSet_ if successful
                void set(T newValue)
                {
                    if (isSetByUser_)
                    {
                        GMX_THROW(gmx::APIError("User-defined hardware options should not be overwritten"));
                    }
                    value_ = newValue;
                    isSet_ = true;
                }
            private:
                //! Parameter is currently set, either by the user or Gromacs
                bool isSet_;
                //! Parameter is set by the user
                bool isSetByUser_;
                //! Hardware option parameter
                T    value_;
        };

        //! Default constructor
        hardwareOptionsManager();
        //! construct from initial options
        explicit hardwareOptionsManager(const gmx_hw_opt_t &user_hw_opt);
        //! copy construct from another instance
        hardwareOptionsManager(const hardwareOptionsManager &optionsManager) = default;
        //! copy assign from another instance
        hardwareOptionsManager &operator=(const hardwareOptionsManager &rhs) = default;
        //! move construct
        hardwareOptionsManager(hardwareOptionsManager &&optionsManager) = default;
        //! move assign
        hardwareOptionsManager &operator=(hardwareOptionsManager &&rhs) = default;
        //! destructor
        ~hardwareOptionsManager() = default;

        //! Total threads
        hardwareOption<int>         nthreads_tot;
        //! Total mpi or thread-mpi ranks
        hardwareOption<int>         nthreads_tmpi;
        //! openMP threads per rank
        hardwareOption<int>         nthreads_omp;
        //! openMP threads per pme rank. Defaults to nthreads_omp
        hardwareOption<int>         nthreads_omp_pme;
        //! Setting of thread affinity, options are in enum
        hardwareOption<int>         threadAffinity;
        //! Stride of core pinning. Usually 1
        hardwareOption<int>         corePinningStride;
        //! Core offset from 0. Useful when running multiple simulations on a node
        hardwareOption<int>         corePinningOffset;
        //! String of GPU IDs for Gromacs to use
        hardwareOption<std::string> gpuIdsAvailable;
        //! String mapping tasks to GPU IDs
        hardwareOption<std::string> gpuTaskAssignment;

};


inline hardwareOptionsManager::hardwareOptionsManager() :
    nthreads_tot(0, 0),
    nthreads_tmpi(0, 0),
    nthreads_omp(0, 0),
    nthreads_omp_pme(0, 0),
    // SEL = 0, AUTO = 1, sel and auto are equivalent
    threadAffinity(threadaffSEL, threadaffAUTO),
    corePinningStride(0, 0),
    corePinningOffset(0, 0),
    gpuIdsAvailable("", ""),
    gpuTaskAssignment("", "")
{}

inline hardwareOptionsManager::hardwareOptionsManager(const gmx_hw_opt_t &user_hw_opt) :
    nthreads_tot(user_hw_opt.nthreads_tot, 0),
    nthreads_tmpi(user_hw_opt.nthreads_tmpi, 0),
    nthreads_omp(user_hw_opt.nthreads_omp, 0),
    nthreads_omp_pme(user_hw_opt.nthreads_omp_pme, 0),
    threadAffinity(user_hw_opt.thread_affinity, threadaffAUTO),
    corePinningStride(user_hw_opt.core_pinning_stride, 0),
    corePinningOffset(user_hw_opt.core_pinning_offset, 0),
    gpuIdsAvailable(user_hw_opt.gpuIdsAvailable, ""),
    gpuTaskAssignment(user_hw_opt.userGpuTaskAssignment, "")
{}



template<typename T> hardwareOptionsManager::hardwareOption<T>::hardwareOption(T value, T defaultValue)
{
    value_     = value;
    // check if value is greater than default (mostly zeros), to account for
    // user defined -1's, which aren't really "set" by the user, but an indication
    // that Gromacs should set them. Also works for the enums and strings
    isSetByUser_ = (value > defaultValue);
    isSet_       = isSetByUser_;
}


//! \Brief Checks if the total thread count is fixed or can be set by Gromacs
inline bool totNumThreadsIsAuto(const hardwareOptionsManager &hardwareOptions)
{
    return (!hardwareOptions.nthreads_tot.isSet() &&
            !hardwareOptions.nthreads_tmpi.isSet() &&
            !hardwareOptions.nthreads_omp.isSet());
}

#endif
