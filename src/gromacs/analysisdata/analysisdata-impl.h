/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Declares private implementation class for gmx::AnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ANALYSISDATA_IMPL_H
#define GMX_ANALYSISDATA_ANALYSISDATA_IMPL_H

#include <vector>

#include "gromacs/utility/uniqueptr.h"

#include "analysisdata.h"
#include "datastorage.h"

namespace gmx
{

namespace internal
{
/*! \internal \brief
 * Private implementation class for AnalysisDataHandle.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataHandleImpl
{
    public:
        //! Creates a handle associated with the given data object.
        explicit AnalysisDataHandleImpl(AnalysisData *data);

        //! The data object that this handle belongs to.
        AnalysisData             &data_;
        //! Current storage frame object, or NULL if no current frame.
        AnalysisDataStorageFrame *currentFrame_;
};
} // namespace internal

/*! \internal \brief
 * Private implementation class for AnalysisData.
 *
 * \ingroup module_analysisdata
 */
class AnalysisData::Impl
{
    public:
        //! Smart pointer type to manage a data handle implementation.
        typedef gmx_unique_ptr<internal::AnalysisDataHandleImpl>::type
                HandlePointer;
        //! Shorthand for a list of data handles.
        typedef std::vector<HandlePointer> HandleList;

        Impl();
        ~Impl();

        //! Storage implementation.
        AnalysisDataStorage     storage_;
        /*! \brief
         * List of handles for this data object.
         *
         * Note that AnalysisDataHandle objects also contain (raw) pointers
         * to these objects.
         */
        HandleList              handles_;
};

} // namespace gmx

#endif
