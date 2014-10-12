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

%ModuleHeaderCode
#include <exception>
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/selection/selectionoption.h"
#include <map>
#include <string>

class PyOptionsHolder {
public:
    class DuplicateOption : public std::exception {
    public:
        virtual const char* what() const throw();
    };
    class StringList {
    public:
        StringList(const std::vector<std::string>*);
        StringList(const std::string*, size_t);
        const char* operator[] (size_t);
    private:
        const std::string *list;
        size_t size;
        const std::vector<std::string> *vector;
    };
    gmx::DoubleOption doubleOption(const char*, size_t = 1, double = 0);
    gmx::IntegerOption integerOption(const char*, size_t = 1, int = 0);
    gmx::StringOption stringOption(const char*, size_t = 1, const char* = "");
    gmx::BooleanOption booleanOption(const char*, size_t = 1, bool = 0);
    gmx::SelectionOption selectionOption(const char*, size_t = 1);
    gmx::FileNameOption fileNameOption(const char*, size_t = 1, const char* = 0);
    PyObject* get_value(const char*);
    ~PyOptionsHolder();
private:
    struct option_value {
        void *value;
        const char *type;
        int *count;
        bool vector;
        const char *name;
    };
    std::map<std::string, option_value> storage;
    template<typename T, typename U> U createStorage(const char*, const char*, int);
    template<typename T> PyObject* buildValue(const char*, const option_value&, const sipTypeDef* = NULL, size_t = 0);
    template<typename T> void deleteStorage(const option_value&);
};
%End

%ModuleCode
#include "gromacs/utility/gmxassert.h"

const char* PyOptionsHolder::DuplicateOption::what() const throw() {
    return "This option is already defined";
}

PyOptionsHolder::StringList::StringList(const std::vector<std::string> *vector) :
    list(NULL), size(vector->size()), vector(vector)
{}

PyOptionsHolder::StringList::StringList(const std::string *list, size_t size) :
    list(list), size(size), vector(NULL)
{}

const char* PyOptionsHolder::StringList::operator[](size_t i) {
    if (i >= size)
        return NULL;
    return vector ? (*vector)[i].data() : list[i].data();
}

template<typename T, typename U> U PyOptionsHolder::createStorage(const char *name, const char *type, int count) {
    if (storage.count(name))
        throw DuplicateOption();

    U option(name);

    int *saved_count = new int;
    option.storeCount(saved_count);

    option_value value;
    value.type = type;
    value.count = saved_count;
    value.name = name;

    if (count == 0) { // std::vector of values
        std::vector<T> *store = new std::vector<T>();
        option.storeVector(store);
        option.multiValue();

        value.value = store;
        value.vector = true;
    } else { // exactly `count` values
        T *store = new T[count];
        option.store(store);
        if (count > 1)
            option.valueCount(count);

        value.value = store;
        value.vector = false;
    }

    storage[name] = value;
    return option;
}

gmx::DoubleOption PyOptionsHolder::doubleOption(const char *name, size_t count, double def) {
    gmx::DoubleOption option = createStorage<double, gmx::DoubleOption>(name, "d", count);

    return option;
}

gmx::IntegerOption PyOptionsHolder::integerOption(const char *name, size_t count, int def) {
    gmx::IntegerOption option = createStorage<int, gmx::IntegerOption>(name, "i", count);

    return option;
}

gmx::StringOption PyOptionsHolder::stringOption(const char *name, size_t count, const char* def) {
    gmx::StringOption option = createStorage<std::string, gmx::StringOption>(name, "A", count);

    return option;
}

// FIXME: Unify with the rest
gmx::BooleanOption PyOptionsHolder::booleanOption(const char *name, size_t count, bool def) {
    if (storage.count(name))
        throw DuplicateOption();

    bool *store = new bool;
    *store = def;
    gmx::BooleanOption option(name);
    option.store(store);
    option_value value = {store, "b"};
    storage[name] = value;

    return option;
}

gmx::SelectionOption PyOptionsHolder::selectionOption(const char *name, size_t count) {
    gmx::SelectionOption option = createStorage<gmx::Selection, gmx::SelectionOption>(name, "S", count);

    return option;
}

gmx::FileNameOption PyOptionsHolder::fileNameOption(const char *name, size_t count, const char* def) {
    gmx::FileNameOption option = createStorage<std::string, gmx::FileNameOption>(name, "f", count);

    return option;
}

template<typename T> PyObject* PyOptionsHolder::buildValue(const char *sipType, const PyOptionsHolder::option_value &value, const sipTypeDef *td, size_t size) {
    int count = *value.count;
    T *store = !value.vector ? (T*) value.value : ((std::vector<T>*) value.value)->data();

    if (count == 1 && !value.vector) {
        if (td)
            return sipConvertFromType(store, td, NULL);
        else
            return sipBuildResult(NULL, sipType, *store);
    } else {
        if (td)
            return sipConvertToTypedArray(store, td, NULL, size, count, SIP_READ_ONLY);
        else
            return sipConvertToArray(store, sipType, count, SIP_READ_ONLY);
    }
}


PyObject* PyOptionsHolder::get_value(const char *name) {
    if (!storage.count(name))
        return NULL;

    option_value v = storage[name];
    switch (*v.type) {
    case 'd': // double
        return buildValue<double>("d", v);
        break;
    case 'i': // int
        return buildValue<int>("i", v);
        break;
    case 'A': // std::string
    case 'f': // FileNameOption
        if (v.vector)
            return sipConvertFromType(new PyOptionsHolder::StringList(((std::vector<std::string>*) v.value)), sipType_PyOptionsHolder_StringList, Py_None);
        else if (*v.count == 1)
            return sipBuildResult(NULL, "A", ((std::string*) v.value)->data());
        else
            return sipConvertFromType(new PyOptionsHolder::StringList((std::string*) v.value, *v.count), sipType_PyOptionsHolder_StringList, Py_None);
        break;
    case 'b': // bool
        // TODO: support vector bool options
        return sipBuildResult(NULL, v.type, *((bool*) v.value));
        break;
    case 'S': // Selection
        return buildValue<gmx::Selection>(NULL, v, sipType_Selection, sizeof(gmx::Selection));
        break;
    }

    GMX_ASSERT(false, "Some type is not handled in PyOptionsHolder.get_value");
    return NULL;
}

template<typename T> void PyOptionsHolder::deleteStorage(const PyOptionsHolder::option_value &value) {
    if (value.vector)
        delete (std::vector<T> *) value.value;
    else
        delete[] (T*) value.value;

    delete value.count;
}

PyOptionsHolder::~PyOptionsHolder() {
    for (std::map<std::string, option_value>::iterator it = storage.begin(); it != storage.end(); it++)
        switch (*(it->second.type)) {
            case 'd': // double
                deleteStorage<double>(it->second);
                break;
            case 'i': // int
                deleteStorage<int>(it->second);
                break;
            case 'A': // std::string
            case 'f': // FileNameOption (the storage type is still std::string)
                deleteStorage<std::string>(it->second);
                break;
            case 'b': // bool
                delete (bool*) it->second.value;
                break;
            case 'S': // Selection
                deleteStorage<gmx::Selection>(it->second);
                break;
        }
}
%End

%Exception PyOptionsHolder::DuplicateOption {
%RaiseCode
        SIP_BLOCK_THREADS;
        PyErr_SetString(PyExc_ValueError, sipExceptionRef.what());
        SIP_UNBLOCK_THREADS;
%End
};

class PyOptionsHolder {
%TypeHeaderCode
#include "gromacs/options/basicoptions.h"
%End

public:
    class StringList /NoDefaultCtors/ {
    public:
        const char* operator[] (int);
        %MethodCode
            sipRes = (*sipCpp)[a0];

            if (!sipRes) {
                SIP_BLOCK_THREADS;
                PyErr_SetString(PyExc_IndexError, "Index out of range");
                SIP_UNBLOCK_THREADS;
                sipIsErr = 1;
            }
        %End
    };
    DoubleOption doubleOption(const char *name /KeepReference/, int count = 1, double def = 0) throw (PyOptionsHolder::DuplicateOption);
    IntegerOption integerOption(const char *name /KeepReference/, int count = 1, int def = 0) throw (PyOptionsHolder::DuplicateOption);
    StringOption stringOption(const char *name /KeepReference/, int count = 1, const char *def = 0) throw (PyOptionsHolder::DuplicateOption);
    BooleanOption booleanOption(const char *name /KeepReference/, int count = 1, bool def = 0) throw (PyOptionsHolder::DuplicateOption);
    SelectionOption selectionOption(const char *name /KeepReference/, int count = 1) throw (PyOptionsHolder::DuplicateOption);
    FileNameOption fileNameOption(const char *name /KeepReference/, int count = 1, const char *def = 0) throw (PyOptionsHolder::DuplicateOption);
    SIP_PYOBJECT __getitem__(const char *name);
    %MethodCode
        sipRes = sipCpp->get_value(a0);

        if (!sipRes) {
            SIP_BLOCK_THREADS;
            PyErr_SetString(PyExc_KeyError, "Invalid option name");
            SIP_UNBLOCK_THREADS;
            sipIsErr = 1;
        }
    %End
};
