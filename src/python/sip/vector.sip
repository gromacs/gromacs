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
// SIP support for std::vector
// by Giovanni Bajo <rasky <at> develer.com>
// Public domain

// ****************************************************
// SIP generic implementation for std::vector<>
// ****************************************************
// ALas, this template-based generic implementation is valid only
// if the element type is a SIP-wrapped type. For basic types (int, double, etc.)
// we are forced to cut & paste to provide a specialization.

template<TYPE*>
%MappedType std::vector<TYPE*>
{
%TypeHeaderCode
#include <vector>
%End

%ConvertFromTypeCode
    PyObject *l = PyList_New(sipCpp -> size());

    // Create the Python list of the correct length.
    if (!l)
        return NULL;

    // Go through each element in the C++ instance and convert it to a
    // wrapped P2d.
    for (int i = 0; i < (int)sipCpp->size(); ++i) {
        TYPE *cpp = sipCpp->at(i);
        PyObject *pobj = sipConvertFromInstance(cpp, sipClass_TYPE, sipTransferObj);

        // Get the Python wrapper for the Type instance, creating a new
        // one if necessary, and handle any ownership transfer.
        if (!pobj) {
            // There was an error so garbage collect the Python list.
            Py_DECREF(l);
            return NULL;
        }

        // Add the wrapper to the list.
        PyList_SET_ITEM(l, i, pobj);
    }

    // Return the Python list.
    return l;
%End

%ConvertToTypeCode
    // Check if type is compatible
    if (!sipIsErr) {
        // Must be any iterable
        PyObject *i = PyObject_GetIter(sipPy);
        bool iterable = (i != NULL);
        Py_XDECREF(i);
        return iterable;
    }

    // Iterate over the object
    PyObject *iterator = PyObject_GetIter(sipPy);
    PyObject *item;

    std::vector<TYPE*> *V = new std::vector<TYPE*>();

    while ((item = PyIter_Next(iterator)))
    {
        if (!sipCanConvertToInstance(item, sipClass_TYPE, SIP_NOT_NONE)) {
            PyErr_Format(PyExc_TypeError, "object in iterable cannot be converted to TYPE");
            *sipIsErr = 1;
            break;
        }

        int state;
        TYPE* p = reinterpret_cast<TYPE*>(
             sipConvertToInstance(item, sipClass_TYPE, 0, SIP_NOT_NONE, &state, sipIsErr));

        if (!*sipIsErr)
            V->push_back(p);

        sipReleaseInstance(p, sipClass_TYPE, state);
        Py_DECREF(item);
    }

    Py_DECREF(iterator);

    if (*sipIsErr) {
        delete V;
        return 0;
    }

    *sipCppPtr = V;
    return sipGetState(sipTransferObj);
%End
};

// ****************************************************
// Specialization for std::vector<char*>
// ****************************************************

%MappedType std::vector<char*>
{
%TypeHeaderCode
#include <vector>
%End

%ConvertFromTypeCode
    PyObject *l;

    // Create the Python list of the correct length.
    if ((l = PyList_New(sipCpp -> size())) == NULL)
        return NULL;

    // Go through each element in the C++ instance and convert it to a
    // wrapped object.
    for (int i = 0; i < (int)sipCpp -> size(); ++i)
    {
        // Add the wrapper to the list.
        PyList_SET_ITEM(l, i, sipBuildResult(NULL, "A", sipCpp->at(i)));
    }

    // Return the Python list.
    return l;
%End

%ConvertToTypeCode
    // Check if type is compatible
    if (sipIsErr == NULL)
    {
        // Must be any iterable
        PyObject *i = PyObject_GetIter(sipPy);
        bool iterable = (i != NULL);
        Py_XDECREF(i);
        return iterable;
    }

    // Iterate over the object
    PyObject *iterator = PyObject_GetIter(sipPy);
    PyObject *item;

    // Maximum number of elements
    int len = PyObject_Size(sipPy);
    std::vector<char*> *V = new std::vector<char*>();
    V->reserve(len);

    if (len)
    {
        while ((item = PyIter_Next(iterator)))
        {
            if (!PyUnicode_Check(item))
            {
                PyErr_Format(PyExc_TypeError, "object in iterable is not a string");
                *sipIsErr = 1;
                break;
            }
            //FIXME: memory leak here
            const char* str = strdup(sipString_AsLatin1String(&item));

            V->push_back(const_cast<char*>(str));

            Py_DECREF(item);
        }

        Py_DECREF(iterator);

        if (*sipIsErr)
        {
            delete V;
            return 0;
        }
    }

    *sipCppPtr = V;
    return sipGetState(sipTransferObj);
%End
};

/*
// ****************************************************
// Specialization for std::vector<double>
// ****************************************************

%MappedType std::vector<double>
{
%TypeHeaderCode
#include <vector>
%End

%ConvertFromTypeCode
    PyObject *l;

    // Create the Python list of the correct length.
    if ((l = PyList_New(sipCpp -> size())) == NULL)
        return NULL;

    // Go through each element in the C++ instance and convert it to a
    // wrapped object.
    for (int i = 0; i < (int)sipCpp -> size(); ++i)
    {
        // Add the wrapper to the list.
        PyList_SET_ITEM(l, i, PyFloat_FromDouble(sipCpp -> at(i)));
    }

    // Return the Python list.
    return l;
%End

%ConvertToTypeCode
    // Check if type is compatible
    if (sipIsErr == NULL)
    {
        // Must be any iterable
        PyObject *i = PyObject_GetIter(sipPy);
        bool iterable = (i != NULL);
        Py_XDECREF(i);
        return iterable;
    }

    // Iterate over the object
    PyObject *iterator = PyObject_GetIter(sipPy);
    PyObject *item;

    // Maximum number of elements
    int len = PyObject_Size(sipPy);
    std::vector<double> *V = new std::vector<double>();
    V->reserve(len);

    if (len)
    {
        while ((item = PyIter_Next(iterator)))
        {
            if (!PyNumber_Check(item))
            {
                PyErr_Format(PyExc_TypeError, "object in iterable is not a number");
                *sipIsErr = 1;
                break;
            }

            PyObject *f = PyNumber_Float(item);
            V->push_back(PyFloat_AsDouble(f));

            Py_DECREF(f);
            Py_DECREF(item);
        }

        Py_DECREF(iterator);

        if (*sipIsErr)
        {
            delete V;
            return 0;
        }
    }

    *sipCppPtr = V;
    return sipGetState(sipTransferObj);
%End
};

// ****************************************************
// Specialization for std::vector<int>
// ****************************************************

%MappedType std::vector<int>
{
%TypeHeaderCode
#include <vector>
%End

%ConvertFromTypeCode
    PyObject *l;

    // Create the Python list of the correct length.
    if ((l = PyList_New((SIP_SSIZE_T)sipCpp -> size())) == NULL)
        return NULL;

    // Go through each element in the C++ instance and convert it to a
    // wrapped object.
    for (int i = 0; i < (int)sipCpp -> size(); ++i)
    {
        // Add the wrapper to the list.
        PyList_SET_ITEM(l, i, PyInt_FromLong(sipCpp -> at(i)));
    }

    // Return the Python list.
    return l;
%End

%ConvertToTypeCode
    // Check if type is compatible
    if (sipIsErr == NULL)
    {
        // Must be any iterable
        PyObject *i = PyObject_GetIter(sipPy);
        bool iterable = (i != NULL);
        Py_XDECREF(i);
        return iterable;
    }

    // Iterate over the object
    PyObject *iterator = PyObject_GetIter(sipPy);
    PyObject *item;

    // Maximum number of elements
    int len = PyObject_Size(sipPy);
    std::vector<int> *V = new std::vector<int>();
    V->reserve(len);

    if (len)
    {
        while ((item = PyIter_Next(iterator)))
        {
            if (!PyInt_Check(item))
            {
                PyErr_Format(PyExc_TypeError, "object in iterable cannot be converted to float");
                *sipIsErr = 1;
                break;
            }

            int val = PyInt_AsLong(item);
            V->push_back(val);

            Py_DECREF(item);
        }

        Py_DECREF(iterator);

        if (*sipIsErr)
        {
            delete V;
            return 0;
        }
    }

    *sipCppPtr = V;
    return sipGetState(sipTransferObj);
%End
};


// ****************************************************
// Specialization for std::vector<unsigned int>
// ****************************************************

%MappedType std::vector<unsigned int>
{
%TypeHeaderCode
#include <vector>
%End

%ConvertFromTypeCode
    PyObject *l;

    // Create the Python list of the correct length.
    if ((l = PyList_New(sipCpp -> size())) == NULL)
        return NULL;

    // Go through each element in the C++ instance and convert it to a
    // wrapped object.
    for (int i = 0; i < (int)sipCpp -> size(); ++i)
    {
        // Add the wrapper to the list.
        PyList_SET_ITEM(l, i, PyInt_FromLong(sipCpp -> at(i)));
    }

    // Return the Python list.
    return l;
%End

%ConvertToTypeCode
    // Check if type is compatible
    if (sipIsErr == NULL)
    {
        // Must be any iterable
        PyObject *i = PyObject_GetIter(sipPy);
        bool iterable = (i != NULL);
        Py_XDECREF(i);
        return iterable;
    }

    // Iterate over the object
    PyObject *iterator = PyObject_GetIter(sipPy);
    PyObject *item;

    // Maximum number of elements
    Py_ssize_t size = PyObject_Size(sipPy);
    if (size == -1) {
        Py_DECREF(iterator);
        return 0;
    }

    unsigned int len = size;
    std::vector<unsigned int> *V = new std::vector<unsigned int>();
    V->reserve(len);

    if (len)
    {
        while ((item = PyIter_Next(iterator)))
        {
            if (!PyInt_Check(item))
            {
                PyErr_Format(PyExc_TypeError, "object in iterable cannot be converted to float");
                *sipIsErr = 1;
                break;
            }

            unsigned int val = PyInt_AsLong(item);
            V->push_back(val);

            Py_DECREF(item);
        }

        Py_DECREF(iterator);

        if (*sipIsErr)
        {
            delete V;
            return 0;
        }
    }

    *sipCppPtr = V;
    return sipGetState(sipTransferObj);
%End
};*/
