/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#include "gromacs/utility/fatalerror.h"
#include "wallcycle.h"

#ifndef GMX_TRACING_H
#define GMX_TRACING_H

#ifdef HAVE_EXTRAE
#define FUNC_QUALIFIER
#define FUNC_TERM ;
/* Set an identifier for Gromacs events.
 * It could be any number, but preferably between 1000-10000
 */
#define EXTRAE_GMX_EVENT  1001

#else
#define FUNC_QUALIFIER static
#define FUNC_TERM {}
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* start the tracer */
FUNC_QUALIFIER
void gmx_tracer_start() FUNC_TERM

/* stop the tracer */
FUNC_QUALIFIER
void gmx_tracer_stop() FUNC_TERM

/* turn-on instrumentation */
FUNC_QUALIFIER
void gmx_tracer_resume() FUNC_TERM

/* turn-off instrumentation */
FUNC_QUALIFIER
void gmx_tracer_pause() FUNC_TERM

/* set a marker for tracing a given event */
FUNC_QUALIFIER
void start_range(int epem) FUNC_TERM

/* unset the event marker */
FUNC_QUALIFIER
void stop_range(int epem) FUNC_TERM

/* read tracing preferences */
FUNC_QUALIFIER
void gmx_tracer_readconf(gmx_wallcycle_t wc) FUNC_TERM

/* push an element to a stack of events to be traced */
FUNC_QUALIFIER
void push(int) FUNC_TERM

/* pop an element*/
FUNC_QUALIFIER
int  pop(void) FUNC_TERM

/* whats on the top of the stack */
FUNC_QUALIFIER
int eventStackTop(void) FUNC_TERM
#ifdef __cplusplus
}
#endif

#undef FUNC_TERM
#undef FUNC_QUALIFIER

#endif  /* GMX_TRACING_H */
