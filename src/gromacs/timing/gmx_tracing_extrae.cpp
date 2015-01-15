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
#include "config.h"
#include <stdlib.h>
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gmx_tracing.h"
#include "extrae_user_events.h"
#include "gromacs/timing/wallcycle.h"

// TODO: remove when done with development
#include "gromacs/utility/fatalerror.h"

// Define max size for the stack of events
// It corresponds to the depth of nesting of ranges
#define TRACE_EVENT_MAXDEPTH 30

struct stack
{
    int stk[TRACE_EVENT_MAXDEPTH];
    int top;
};
typedef struct stack STACK;
STACK events_stack;

void push(int);
int  pop(void);

/* start the tracer */
void gmx_tracer_start()
{
    // Extrae is automatically initialized upon launch
    //Extrae_init();

    int            i;
    extrae_type_t  eventID   = EXTRAE_GMX_EVENT;
    char          *type_desc = (char *) "GMX_EVENT";
    unsigned int   nevents   = ewcNR;
    extrae_value_t event[ewcNR + 1];
    char          *event_labels[ewcNR+1];


    /* Define labels for the event markers */
    /* we have ewcNR+1 events, element 0 is a dummy */
    event_labels[0] = gmx_strdup("End");
    event[0]        = 0;

    for (i = 0; i < ewcNR; i++)
    {
        event[i+1]        = i+1;
        event_labels[i+1] = gmx_strdup(wcn_name_get(i));
    }

    Extrae_define_event_type(&eventID, type_desc, &nevents, event, event_labels);

/* We deallocate here because Extrae stores the info in the .sym file */
    for (i = 0; i < ewcNR+1; i++)
    {
        sfree(event_labels[i]);
    }

    pop();
};

/* stop the tracer */
void gmx_tracer_stop()
{
    Extrae_fini();
};


/* turn-on instrumentation */
void gmx_tracer_resume()
{

    Extrae_restart();
};

/* tunr-off instrumentation */
void gmx_tracer_pause()
{

    Extrae_shutdown();
};

/* set a marker for tracing a given event */
void start_range(int epem)
{

// Uncomment for debugging
// gmx_warning("PROFILER: Start range Event %d/n", epem);

    // printf(">>> Start range Event %d\n", epem);
    Extrae_event(EXTRAE_GMX_EVENT, epem+1);

};

/* unset the event marker */
void stop_range(int epem)
{

// Uncomment for debugging
// gmx_warning("PROFILER: STOP range Event %d/n", epem);

    //printf("xxx Stop range Event %d\n", epem);
    Extrae_event(EXTRAE_GMX_EVENT, 0);
};

void push (int eventID)
{
    if (events_stack.top == (TRACE_EVENT_MAXDEPTH - 1))
    {
        printf ("Stack is Full\n");
        return;
    }
    else
    {
        events_stack.top = events_stack.top + 1;

        events_stack.stk[events_stack.top] = eventID;
    }
    return;
}

int pop ()
{
    int eventID;
    if (events_stack.top == -1)
    {
        //  the stack is empty
        return (events_stack.top);
    }
    else
    {
        eventID          = events_stack.stk[events_stack.top];
        events_stack.top = events_stack.top - 1;
    }
    return(eventID);
}

int eventStackTop()
{

    return events_stack.top;

}
