/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

%Module gromacs.TrajectoryAnalysis

%DefaultEncoding "UTF-8"

%Import options/Options.sip
%Import topology/Topology.sip

%Include definitions.sip

%Include analysissettings.sip
%Include analysismodule.sip

%ModuleCode
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/exceptions.h"

class PyFactory : public gmx::TrajectoryAnalysisCommandLineRunner::ModuleFactoryFunctor {
public:
    PyFactory(PyObject *mod) : mod(mod) {} ;
    gmx::TrajectoryAnalysisModulePointer operator () () {
        int iserr = 0;
        int can = sipCanConvertToType(mod, sipType_TrajectoryAnalysisModule, 0);
        // TODO: throw if can't
        gmx::TrajectoryAnalysisModule *module = (gmx::TrajectoryAnalysisModule*) sipConvertToType(
            mod, sipType_TrajectoryAnalysisModule, NULL, 0, NULL, &iserr);
        return TrajectoryAnalysisModulePointer(module);
    }
private:
    PyObject *mod;
};

void runAsMain(PyObject *mod, PyObject *py_argv) {
    int argc = PyList_GET_SIZE(py_argv);

    char **argv = new char *[argc + 1];

    // Convert the list.
    for (int a = 0; a < argc; ++a)
    {
        PyObject *arg_obj = PyList_GET_ITEM(py_argv, a);
        const char *arg = sipString_AsLatin1String(&arg_obj);

        if (arg)
        {
            arg = strdup(arg);
            Py_DECREF(arg_obj);
        }
        else
        {
            arg = "unknown";
        }

        argv[a] = const_cast<char *>(arg);
    }

    argv[argc] = NULL;

    PyFactory factory(mod);
    TrajectoryAnalysisCommandLineRunner::runAsMain(argc, argv, &factory);
}

%End

void runAsMain(SIP_PYOBJECT mod, SIP_PYLIST argv);
