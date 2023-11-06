/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief
 * Define basic tracing API for manual instrumentation.
 *
 * This header implements a simple set of tracing range start/stop functions
 * which can be used for manual instrumentation of application code.
 * Since current use is only through the wallcycle module, we define two
 * sets of start/stop functions corresponding to the main and sub-counters
 * in the wallcycle module.
 *
 * The current implementation supports the following tracing APIs:
 * - NVIDIA NVTX
 * - AMD ROCTX
 * - Intel ITT
 *
 * \author Szilárd Páll <pall.szilard@gmail.com
 * \author Andrey Alekseenko <al42and@gmail.com>
 *
 */

#include "gromacs/utility/basedefinitions.h"


//
// Forward declarations of the tracing functions with the inlineable definitions for each tracing API below.
//

/*! \brief Start a main tracing region.
 *
 * Note that the \p rangeId argument is currently only used with NVTX for aiding
 * in coloring of trace regions.
 *
 * \param[in] rangeName   String containing the name of the traced range.
 * \param[in] rangeId     Numeric ID of the range.
 */
static void traceRangeStart(const char* rangeName, int rangeId);


/*! \brief Start a tracing sub-region.
 *
 * Note that the \p rangeId argument is currently only used with NVTX for aiding
 * in coloring of trace regions.
 *
 * \param[in] rangeName   String containing the name of the traced range.
 * \param[in] rangeId     Numeric ID of the range.
 */
static void traceSubRangeStart(const char* rangeName, int rangeId);

/*! \brief End a main tracing region.
 *
 * Note that this should always be paired with a traceRangeStart().
 */
static void traceRangeEnd();

/*! \brief End a tracing sub-region.
 *
 * Note that this should always be paired with a traceSubRangeStart().
 */
static void traceSubRangeEnd();

#if (GMX_USE_NVTX + GMX_USE_ROCTX + GMX_USE_ITT) > 1
#    error "Cannot have multiple instrumentation flavors enabled at the same time"
#endif


#if GMX_USE_NVTX

#    include "nvToolsExt.h"

//! List of colors for main ranges
static constexpr uint32_t c_rangeColors[] = { 0xff00ff00, 0xff0000ff, 0xffffff00, 0xffff00ff,
                                              0xff00ffff, 0xffff0000, 0xffffffff };
//! Number of colors for main ranges
static constexpr int c_numRangeColors = sizeof(c_rangeColors) / sizeof(uint32_t);

//! List of colors for sub-ranges
static constexpr uint32_t c_subRangeColors[] = { 0x9900ff00, 0x990000ff, 0x99ffff00, 0x99ff00ff,
                                                 0x9900ffff, 0x99ff0000, 0x99ffffff };
//! Number of colors for sub-ranges
static constexpr int c_numSubRangeColors = sizeof(c_subRangeColors) / sizeof(uint32_t);

static void traceRangeStart(const char* rangeName, int rangeId)
{
    int                   colorId     = rangeId % c_numRangeColors;
    nvtxEventAttributes_t eventAttrib = { 0 };
    eventAttrib.version               = NVTX_VERSION;
    eventAttrib.size                  = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType             = NVTX_COLOR_ARGB;
    eventAttrib.color                 = c_rangeColors[colorId];
    eventAttrib.messageType           = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii         = rangeName;
    nvtxRangePushEx(&eventAttrib);
}

static void traceSubRangeStart(const char* rangeName, int rangeId)
{
    int                   colorId     = rangeId % c_numSubRangeColors;
    nvtxEventAttributes_t eventAttrib = { 0 };
    eventAttrib.version               = NVTX_VERSION;
    eventAttrib.size                  = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType             = NVTX_COLOR_ARGB;
    eventAttrib.color                 = c_subRangeColors[colorId];
    eventAttrib.messageType           = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii         = rangeName;
    nvtxRangePushEx(&eventAttrib);
}

static void traceRangeEnd()
{
    nvtxRangePop();
}

static void traceSubRangeEnd()
{
    nvtxRangePop();
}

#elif GMX_USE_ROCTX

#    include "roctracer/roctx.h"

static void traceRangeStart(const char* rangeName, int /*rangeId*/)
{
    roctxRangePush(rangeName);
}

static void traceSubRangeStart(const char* rangeName, int /*rangeId*/)
{
    roctxRangePush(rangeName);
}

static void traceRangeEnd()
{
    roctxRangePop();
}

static void traceSubRangeEnd()
{
    roctxRangePop();
}

#elif GMX_USE_ITT

#    ifdef __clang__
#        pragma clang diagnostic push
#        pragma clang diagnostic ignored "-Wold-style-cast"
#        pragma clang diagnostic ignored "-Wnewline-eof"
#    endif
#    include <ittnotify.h>
#    ifdef __clang__
#        pragma clang diagnostic pop
#    endif

// Defined in wallcycle.cpp, initialized in wallcycle_init
extern const __itt_domain*  g_ittDomain;
extern __itt_string_handle* g_ittCounterHandles[];
extern __itt_string_handle* g_ittSubCounterHandles[];

static void traceRangeStart(const char* /*rangeName*/, int rangeId)
{
    __itt_task_begin(g_ittDomain, __itt_null, __itt_null, g_ittCounterHandles[rangeId]);
}

static void traceSubRangeStart(const char* /*rangeName*/, int rangeId)
{
    __itt_task_begin(g_ittDomain, __itt_null, __itt_null, g_ittSubCounterHandles[rangeId]);
}

static void traceRangeEnd()
{
    __itt_task_end(g_ittDomain);
}

static void traceSubRangeEnd()
{
    __itt_task_end(g_ittDomain);
}

#else

gmx_unused static void traceRangeStart(gmx_unused const char* rangeName, gmx_unused int rangeId) {}
gmx_unused static void traceSubRangeStart(gmx_unused const char* rangeName, gmx_unused int rangeId)
{
}

gmx_unused static void traceRangeEnd() {}
gmx_unused static void traceSubRangeEnd() {}


#endif
