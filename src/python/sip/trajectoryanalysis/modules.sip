/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
template <T> class shared_ptr /NoDefaultCtors/ {
%TypeHeaderCode
#include <boost/shared_ptr.hpp>
using namespace boost;
%End
};

typedef shared_ptr<TrajectoryAnalysisModule> TrajectoryAnalysisModulePointer;

template <TYPE>
%MappedType py_shared_ptr<TYPE> /NoRelease/ {
%TypeHeaderCode
#include <boost/shared_ptr.hpp>
template <typename T> using py_shared_ptr = boost::shared_ptr<T>;
%End

%ConvertFromTypeCode
    PyObject *wrp = sipConvertFromType(sipCpp->get(), sipType_TYPE, NULL);
    PyObject *ptr = sipConvertFromNewType(sipCpp, sipType_TYPEPointer, wrp);

    return wrp;
%End

%ConvertToTypeCode
    // Make it impossible to convert into this type
    if (!sipIsErr)
        return 0;

    GMX_ASSERT(true, "Converting something into shared_ptr. This should not happen!");
    return 0;
%End
};

class AngleInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/angle.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class DistanceInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/distance.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class FreeVolumeInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/freevolume.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class PairDistanceInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/pairdist.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class RdfInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/rdf.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class SasaInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/sasa.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};

class SelectInfo {
%TypeHeaderCode
#include "gromacs/trajectoryanalysis/modules/select.h"
using namespace gmx::analysismodules;
%End
private:
public:
    static const char* name /NoSetter/;
    static const char* shortDescription /NoSetter/;
    static py_shared_ptr<TrajectoryAnalysisModule> create() /Factory/;
};
