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
/*!\inlibraryapi \file
 * <H3>Overview of coordinate modification (coordinate converter) routines</H3>
 * The methods based on the IFrameConverter interface provide ways to change
 * the coordinates found in input files according to a user or program specified
 * pattern that fulfils a specified requirement. The basic adapter implements a single
 * method that performs one operation on a coordinate file, which provides a guarantee
 * on the produced coordinates. Those guarantees should be simple enough so that they
 * can be fulfilled by the single operation, and avoid complex invocation patterns.
 *
 * <H3>Combination of different frame converters</H3>
 * Several individual converters can be combined and executed in sequence to provide
 * more complex operations. Here, each individual converter can provide part of a complex
 * guarantee requested for the final output structure, depending on the combination
 * and sequence of the adapters. Implementers need to take care to consider if an
 * adapter added to a particular sequence will invalidate a previous guarantee, to
 * make sure there is a way to avoid operations unknowingly invalidating each other.
 *
 * <H3>Data handling for frame conversion methods</H3>
 * Methods that use the IFrameConverter based registration and chaining method do not
 * need handle their own data. The registration method provides data storage for the
 * modified coordinates and returns a final, modified t_trxframe datastructure with
 * the method owned coordinates. No changes are applied to the velocity and force
 * fields (if present), with eventually present data used in the new datastructure.
 *
 * Methods that implement single converters without the registration machinery need to
 * implement their own memory handling.
 *
 * <H3>FrameAdapers</H3>
 * Each class implemented on top of the IFrameConverter interface implements each own
 * convertFrame method that performs its operation on the coordinates in the t_trxframe
 * input. As well as this, the method should report which kind of guarantee is provided
 * by the method, as reported by the guarantee method.
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   requirements [ label="Requirements to run method" ],
   container [ label="FrameConverterChain" ],
   frameconverters [ label="FrameConverters" ];

   analysistool => builder [ label="Requests to run analysis method" ];
   analysistool => builder [ label="Specifies required modifications needed" ];
   requirements => container [ label="Add methods needed for modifications" ];
   container => frameconverters [ label="Adds each converter element" ];
   frameconverters => container [ label="Adds to final guarantees provided by chained methods" ];
   container => analysistool [ label="Reports on final set of method that provide the needed changes" ];
   analysistool box analysistool [ label="Can now run the analysis with the preprocessing from the frameconverters" ];
 * \endmsc
 *
 *
 * \msc
   wordwraparcs=true,
   hscale="2";

   analysistool,
   analysisloop,
   inputcoordinates,
   frameconverterholder [ label="Chain of modification tools" ],
   frameconverters [ label="FrameConverters" ];

   --- [ label="Setup of frameconverters complete, analysis phase begins" ];
    analysistool   => analysisloop [ label="Starts iteration over frames" ];
    analysisloop   => frameconvertholder [ label="Provides initial unmodified coordinates" ];
    frameconverterholder => frameconverters [ label="Each successive converter modifies coordinates" ];
    frameconverters => frameconverterholder [ label="Return coordinates from modification" ];
    frameconverterholder => analysisloop [ label="Return final coordinates to analysis tool for analysis" ];
    analysisloop box analysisloop [ label="Iterates over frames" ];
    --- [ label="Analysis complete, object is destructed and files are closed" ];

 *  \endmsc
 *
 * \if libapi
 * <H3>Preparing new FrameConverters</H3>
 *
 * If additional methods are needed to perform modification of coordinate data,
 * new FrameConverters can be written that again implement the IFrameConverter interface.
 * The new method should follow the approach of the other modules that are present
 * in performing single modifications on the coordinates.
 * \endif
 * \brief
 * Interface class for frame handling, provides handles for all calls.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_coordinateio
 */
#ifndef GMX_COORDINATEIO_IFRAMECONVERTER_H
#define GMX_COORDINATEIO_IFRAMECONVERTER_H

#include <memory>

#include "frameconverterenums.h"

struct t_trxframe;

namespace gmx
{

/*!\inlibraryapi
 * \brief
 * IFrameConverter interface for manipulating coordinate information.
 *
 * This interface is aimed at providing the base methods to manipulate the
 * coordinate (usually position) data in a t_trxframe object according
 * to the requirements of an analysis module. It is similar to the
 * ICoordinateOutput interface, but instead of
 * simply passing through frames returns new t_trxframe objects with
 * changes applied to them.
 *
 * \ingroup module_coordinateio
 *
 */
class IFrameConverter
{
public:
    IFrameConverter() {}

    virtual ~IFrameConverter() {}

    //! Move constructor for old object.
    explicit IFrameConverter(IFrameConverter&& old) noexcept = default;

    /*! \brief
     * Change coordinate frame information for output.
     *
     * Takes the previously internally stored coordinates and saves them again.
     * Applies correct number of atoms, as well as changing things such as
     * frame time or affect output of velocities or forces.
     * Methods derived from this should not affect total number of atoms,
     * and should check for internal consistency of the input and output data.
     *
     * \param[in,out]  input  Coordinate frame to be modified.
     */
    virtual void convertFrame(t_trxframe* input) = 0;

    /*! \brief
     * Allows other methods to probe if a specific requirement is fulfilled by running a converter.
     *
     * When modifying coordinate frames with different frameconverters,
     * it can be important to know what kind of modifications are done by a
     * specific converter, to e.g. check if it makes a system whole or moves
     * the simulation box in a specific way.
     *
     * \returns What kind of modification is guaranteed by this converter.
     */
    virtual unsigned long guarantee() const = 0;
};
/*! \brief
 * Typedef to have direct access to the individual FrameConverter modules.
 */
using FrameConverterPointer = std::unique_ptr<IFrameConverter>;

} // namespace gmx

#endif
