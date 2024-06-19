/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 * \brief
 * Declares gmx::ProcessFrameConversion
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATIO_FRAMECONVERTERS_REGISTER_H
#define GMX_COORDINATIO_FRAMECONVERTERS_REGISTER_H

#include <memory>
#include <utility>
#include <vector>

#include "gromacs/coordinateio/frameconverterenums.h"
#include "gromacs/coordinateio/iframeconverter.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"

struct t_trxframe;

namespace gmx
{

/*!\brief
 * ProcessFrameConversion class for handling the running of several analysis steps
 *
 * This analysis module allows to register modules for coordinate frame manipulation
 * that are then run once this modules convertFrame method is invoked.
 * If modules are registered here, the method will take care of memory management
 * to ensure that input data is not modified.
 *
 * It is possible to chain different versions of this class together to have several
 * independend containers. In this case, only the outermost container will usually
 * own the memory, but it is possible to envision different implementations that take advantage
 * of the indivudal memory owning objects.
 *
 * \inlibraryapi
 * \ingroup module_coordinatedata
 *
 */
class ProcessFrameConversion : public IFrameConverter
{
public:
    /*! \brief
     * Default constructor for ProcessFrameConversion.
     */
    ProcessFrameConversion();

    /*! \brief
     * Default destructor for ProcessFrameConversion.
     */
    ~ProcessFrameConversion() override;

    /*! \brief
     * Change coordinate frame information for output.
     *
     * This method is used to perform the actual coordinate frame manipulation.
     * In this case, it acts as a wrapper that runs the method for all
     * modules that have been registered for the analysis chain.
     *
     * \param[in,out]  input  Coordinate frame to be modified n process chain.
     */
    void convertFrame(t_trxframe* input) override;

    //! Guarantees provided by this method and all included ones.
    unsigned long guarantee() const override { return listOfGuarantees_; }

    /*! \brief
     * Add new guarantee and check if it invalidates previous one.
     *
     * \param[in] guarantee New guarantee provided by frameconverter to add.
     */
    void addAndCheckGuarantee(unsigned long guarantee);

    /*! \brief
     * Add framemodule to analysis chain.
     *
     * Other modules derived from IFrameConverter can be registered to be analysed
     * by this method instead of having to be analysed separately.
     *
     * \param[in] module FrameConverter module to add to the chain.
     * \throws unspecified Any exception thrown by any of the \p module
     * objects in the chain during analysis.
     */
    void addFrameConverter(FrameConverterPointer module);

    //! Get number of converters registered.
    int getNumberOfConverters() const { return moduleChain_.size(); }

    /*! \brief
     * Wrapper method that allows input of a const \p inputFrame.
     *
     * Takes care of preparing a memory owning object that will contain the
     * final coordinates after all modifcations have been applied.
     * As long as this module exists, the memory and the containing object is valid.
     *
     * \param[in] inputFrame Coordinate data input.
     * \returns None owning pointer to new coordinates.
     */
    t_trxframe* prepareAndTransformCoordinates(const t_trxframe* inputFrame);

private:
    /*! \brief
     * Handle memory for coordinate frames being processed.
     *
     * Allocates new memory for methods that change coordiante frames.
     * This ensures that a methods that changes the coordiantes will never
     * change the input data.
     *
     * \param[in] inputFrame Coordinate frame to operate on in this container.
     */
    void prepareNewCoordinates(const t_trxframe* inputFrame);
    /*! \libinternal \brief
     * Storage for info about different modules chained together.
     *
     * This is storing the pointers to the individual methods in the analysis chain.
     * For each method, one pointer to the method is stored to be used in the custom
     * convertFrame method.
     *
     */
    struct FrameModule
    {
        //! Initializes module, stolen from datamodulemanager.
        explicit FrameModule(FrameConverterPointer module) : module_(std::move(module)) {}
        //! Pointer to module.
        FrameConverterPointer module_;
    };
    //! Shorthand for list of chained modules
    using FrameModuleList = std::vector<FrameModule>;

    //! List of chained modules.
    FrameModuleList moduleChain_;

    //! Internal storage object for new coordinate frame.
    std::unique_ptr<t_trxframe> frame_;
    //! Internal storage of new coordinates.
    std::vector<RVec> localX_;
    //! Internal storage of new velocities.
    std::vector<RVec> localV_;
    //! Internal storage of new forces.
    std::vector<RVec> localF_;
    //! Internal storage for number of guarantees provided by chained methods.
    unsigned long listOfGuarantees_ = convertFlag(FrameConverterFlags::NoGuarantee);
};

//! Smart pointer to manage the analyse object.
using ProcessFrameConversionPointer = std::unique_ptr<ProcessFrameConversion>;

} // namespace gmx

#endif
