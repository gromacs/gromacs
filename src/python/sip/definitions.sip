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

%Exception gmx::GromacsException(SIP_Exception) /Default/ {
%TypeHeaderCode
#include <gromacs/utility/exceptions.h>
%End
%RaiseCode
    std::string detail = formatExceptionMessageToString(sipExceptionRef);

    SIP_BLOCK_THREADS
    PyErr_SetString(sipException_gmx_GromacsException, detail.data());
    SIP_UNBLOCK_THREADS
%End
};

// This MappedType makes sure that unique_ptr does not get destructed while wrapped pointer is still needed.
// A class called 'TYPEPointer' needs to be created for every use.

template <TYPE>
%MappedType std::unique_ptr<TYPE> /NoRelease/ {

%ConvertFromTypeCode
    PyObject *wrp = sipConvertFromType(sipCpp->get(), sipType_TYPE, NULL);
    PyObject *ptr = sipConvertFromNewType(sipCpp, sipType_TYPEPointer, wrp);
    // ptr is created only to bind wrp's ownership to it
    GMX_UNUSED_VALUE(ptr);

    return wrp;
%End

%ConvertToTypeCode
    // Make it impossible to convert into this type
    if (!sipIsErr)
        return 0;

    GMX_ASSERT(true, "Converting something into unique_ptr. This should not happen!");
    return 0;
%End
};

%VirtualErrorHandler vehandler
    SIP_RELEASE_GIL(sipGILState);

    PyErr_PrintEx(1);
    GMX_THROW(gmx::InternalError("Python virtual overload raised an exception, see traceback"));
%End

%ModuleHeaderCode
#include <functional>
int argv_wrapper(PyObject *, const std::function<int(int, char**)> &);
template <typename T> class PyFactory;
%End

// This is put into UnitPostIncludeCode so that sip* declarations are available.
%UnitPostIncludeCode
#include <memory>
#include <gromacs/utility/exceptions.h>
template <typename T> class PyFactory {
public:
    PyFactory(const sipTypeDef *sipType, PyObject *mod) : sipType_(sipType), mod(mod) {} ;
    std::unique_ptr<T> operator () () {
        int iserr = 0;
        int can = sipCanConvertToType(mod, sipType_, 0);

        if (!can)
            GMX_THROW(gmx::APIError("Provided python object can not be converted to type required by runAsMain-like function"));

        T *module = (T*) sipConvertToType(
            mod, sipType_, NULL, 0, NULL, &iserr);
        return std::unique_ptr<T>(module);
    }
private:
    const sipTypeDef *sipType_;
    PyObject *mod;
};
%End

%ModuleCode
int argv_wrapper(PyObject *py_argv, const std::function<int(int, char**)> &runner) {
    int argc = PyList_GET_SIZE(py_argv);

    char **argv = new char *[argc + 1];
    bool to_free[argc];
    std::fill(to_free, to_free+argc, false);

    // Convert the list.
    // TODO: Use something better than AsLatin1String to avoid unicode errors
    for (int a = 0; a < argc; ++a)
    {
        PyObject *arg_obj = PyList_GET_ITEM(py_argv, a);
        const char *arg = sipString_AsLatin1String(&arg_obj);

        if (arg)
        {
            arg = strdup(arg);
            Py_DECREF(arg_obj);
            to_free[a] = true;
        }
        else
        {
            arg = "unknown";
        }

        argv[a] = const_cast<char *>(arg);
    }

    argv[argc] = NULL;

    int result = runner(argc, argv);

    for (int i = 0; i < argc; i++)
        if (to_free[i])
            free(argv[i]);
    delete[] argv;

    return result;
}
%End
