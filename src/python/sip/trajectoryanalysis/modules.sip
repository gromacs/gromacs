
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
