/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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

%PostInitialisationCode
    import_array();
%End

struct t_trxframe {

%TypeHeaderCode
#include <gromacs/fileio/trx.h>
#include "numpy_conv.h"
%End

    int flags;
    int not_ok;
    bool bDouble;
    int natoms;
    double t0;
    double tf;
    double tpf;
    double tppf;
    bool bTitle;
    const char *title;
    bool bStep;
    int step;
    bool bTime;
    double time;
    bool bLambda;
    bool bFepState;
    double lambda;
    int fep_state;
    bool bAtoms;
    t_atoms *atoms;
    bool bPrec;
    double prec;
    bool bX;
    SIP_PYOBJECT x {
    %GetCode
        sipPy = array2dToNumpy(sipCpp->natoms, 3, sipCpp->x);
    %End
    %SetCode

    %End
    };
    bool bV;
    SIP_PYOBJECT v {
    %GetCode
        sipPy = array2dToNumpy(sipCpp->natoms, 3, sipCpp->v);
    %End
    %SetCode

    %End
    };
    bool bF;
    SIP_PYOBJECT f {
    %GetCode
        sipPy = array2dToNumpy(sipCpp->natoms, 3, sipCpp->f);
    %End
    %SetCode

    %End
    };
    bool bBox;
    SIP_PYOBJECT box {
    %GetCode
        sipPy = array2dToNumpy(3, 3, sipCpp->box);
    %End
    %SetCode

    %End
    };
    bool bPBC;
    int ePBC;
    // TODO
    /*t_gmxvmdplugin* *vmdplugin; */
};

struct t_pbc {

%TypeHeaderCode
#include <gromacs/pbcutil/pbc.h>
#include "numpy_conv.h"
%End
    int ePBC;
    int ndim_ePBC;
    int ePBCDX;
    int dim;
    SIP_PYOBJECT box {
    %GetCode
        sipPy = array2dToNumpy(3, 3, sipCpp->box);
    %End
    %SetCode

    %End
    };
    SIP_PYOBJECT fbox_diag {
    %GetCode
        sipPy = array1dToNumpy(3, sipCpp->fbox_diag);
    %End
    %SetCode

    %End
    };
    SIP_PYOBJECT hbox_diag {
    %GetCode
        sipPy = array1dToNumpy(3, sipCpp->hbox_diag);
    %End
    %SetCode

    %End
    };
    SIP_PYOBJECT mhbox_diag {
    %GetCode
        sipPy = array1dToNumpy(3, sipCpp->mhbox_diag);
    %End
    %SetCode

    %End
    };
    double max_cutoff2;
    bool bLimitDistance;
    double limit_distance2;
    int  ntric_vec;
    /*ivec tric_shift[MAX_NTRICVEC];*/
    /*rvec tric_vec[MAX_NTRICVEC];*/
};

class TrajectoryAnalysisModuleData /NoDefaultCtors/ {

%TypeHeaderCode
#include <gromacs/trajectoryanalysis/analysismodule.h>
using namespace gmx;
%End

public:
    virtual void finish() = 0;
};

class TrajectoryAnalysisModule {

%TypeHeaderCode
#include <gromacs/trajectoryanalysis/analysismodule.h>
using namespace gmx;
%End
public:
    virtual void initOptions(Options *options, TrajectoryAnalysisSettings *settings) = 0;
    virtual void optionsFinished(Options *options, TrajectoryAnalysisSettings *settings);
    virtual void initAnalysis(const TrajectoryAnalysisSettings &settings, const TopologyInformation &top) = 0;
    virtual void initAfterFirstFrame(const TrajectoryAnalysisSettings &settings, const t_trxframe &fr);
//    virtual TrajectoryAnalysisModuleDataPointer 	startFrames (const AnalysisDataParallelOptions &opt, const SelectionCollection &selections);
    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc, TrajectoryAnalysisModuleData *pdata ) = 0;
    virtual void finishFrames(TrajectoryAnalysisModuleData *pdata);
    virtual void finishAnalysis(int nframes) = 0;
    virtual void writeOutput() = 0;
    const char* name() const;
    const char* description() const;
    int datasetCount() const;
protected:
    TrajectoryAnalysisModule(const char *name, const char *description);
private:
    TrajectoryAnalysisModule(const TrajectoryAnalysisModule &other);
};

%Exception gmx::InconsistentInputError(SIP_Exception) {

%TypeHeaderCode
#include <gromacs/utility/exceptions.h>
using namespace std;
%End

%RaiseCode
    const char *detail = sipExceptionRef.what();

    SIP_BLOCK_THREADS
    PyErr_SetString(sipException_gmx_InconsistentInputError, detail);
    SIP_UNBLOCK_THREADS
%End
};

class TrajectoryAnalysisCommandLineRunner {
%TypeHeaderCode
#include <gromacs/trajectoryanalysis/cmdlinerunner.h>
%End

public:
    TrajectoryAnalysisCommandLineRunner(TrajectoryAnalysisModule *module);
    int run(SIP_PYLIST) throw (gmx::InconsistentInputError);
    %MethodCode
    int argc = PyList_GET_SIZE(a0);

    char **argv = new char *[argc + 1];

    // Convert the list.
    for (int a = 0; a < argc; ++a)
    {
        PyObject *arg_obj = PyList_GET_ITEM(a0, a);
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

    try {
        sipCpp->run(argc, argv);
    } catch (gmx::InconsistentInputError &e) {
        sipIsErr = 1;
        PyErr_SetString(sipException_gmx_InconsistentInputError, e.what());
    }
    %End
private:
    TrajectoryAnalysisCommandLineRunner(const TrajectoryAnalysisCommandLineRunner &other);
};

