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

%PostInitialisationCode
import_array();
%End

enum e_index_t {
    INDEX_UNKNOWN,
    INDEX_ATOM,
    INDEX_RES,
    INDEX_MOL,
    INDEX_ALL
};

class SelectionPosition /NoDefaultCtors/ {
%TypeHeaderCode
#include "gromacs/selection/selection.h"
using namespace gmx;
%End

public:
    e_index_t type();
    SIP_PYOBJECT x();
    %MethodCode
        sipRes = array1dToNumpy(3, sipCpp->x());
    %End
    SIP_PYOBJECT v();
    %MethodCode
        sipRes = array1dToNumpy(3, sipCpp->v());
    %End
    SIP_PYOBJECT f();
    %MethodCode
        sipRes = array1dToNumpy(3, sipCpp->f());
    %End
    double mass();
    double charge();
    int atomCount();
    SIP_PYOBJECT atomIndices();
    %MethodCode
        sipRes = iarray1dToNumpy(sipCpp->atomIndices().size(), sipCpp->atomIndices().data());
    %End
    bool selected();
    int refId();
    int mappedId();
};

enum e_coverfrac_t {
    CFRAC_NONE,
    CFRAC_SOLIDANGLE
};

class Selection {
%TypeHeaderCode
#include "gromacs/selection/selection.h"
#include "numpy_conv.h"
using namespace gmx;
%End

public:
    bool isValid();
    bool __eq__(const Selection&);
    bool __ne__(const Selection&);
    const char* name();
    const char* selectionText();
    bool isDynamic();
    e_index_t type();
    bool hasOnlyAtoms();
    int atomCount();
    SIP_PYOBJECT atomIndices();
    %MethodCode
        sipRes = iarray1dToNumpy(sipCpp->atomIndices().size(), sipCpp->atomIndices().data());
    %End
    int posCount();
    SelectionPosition position(int);
    SIP_PYOBJECT coordinates();
    %MethodCode
        sipRes = array2dToNumpy(sipCpp->coordinates().size(), 3, sipCpp->coordinates().data());
    %End
    bool hasVelocities();
    SIP_PYOBJECT velocities();
    %MethodCode
        sipRes = array2dToNumpy(sipCpp->velocities().size(), 3, sipCpp->velocities().data());
    %End
    bool hasForces();
    SIP_PYOBJECT forces();
    %MethodCode
        sipRes = array2dToNumpy(sipCpp->forces().size(), 3, sipCpp->forces().data());
    %End
    SIP_PYOBJECT masses();
    %MethodCode
        sipRes = array1dToNumpy(sipCpp->masses().size(), sipCpp->masses().data());
    %End
    SIP_PYOBJECT charges();
    %MethodCode
        sipRes = array1dToNumpy(sipCpp->charges().size(), sipCpp->charges().data());
    %End
    SIP_PYOBJECT refIds();
    %MethodCode
        sipRes = iarray1dToNumpy(sipCpp->refIds().size(), sipCpp->refIds().data());
    %End
    SIP_PYOBJECT mappedIds();
    %MethodCode
        sipRes = iarray1dToNumpy(sipCpp->mappedIds().size(), sipCpp->mappedIds().data());
    %End
    bool isCoveredFractionDynamic();
    double coveredFraction();
    bool initCoveredFraction(e_coverfrac_t);
    void setEvaluateVelocities(bool);
    void setEvaluateForces(bool);
    void setOriginalId(int, int);
};
