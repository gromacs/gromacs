/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
/*! \internal \file
 * \brief Test helper utilities for generating frame data.
 *
 * \author Petter Johansson <pettjoha@kth.se>
 * \ingroup module_testutils
 */

#ifndef GMX_TESTUTILS_GENERATE_FRAME_DATA_H
#define GMX_TESTUTILS_GENERATE_FRAME_DATA_H

namespace gmx
{

template<typename T>
class ArrayRef;
template<typename T>
class BasicVector;
template<typename T>
class BasicMatrix3x3;

namespace test
{

//! \brief Mode for trajectory data generation.
enum class TrajectoryFrameMode
{
    Positive,    //!< Positive sign for all values.
    Negative,    //!< Negative sign for all values.
    Alternating, //!< Sign alternating every other value.
};

/*! \brief
 * Generator of per-frame trajectory (array/vector of RVecs) data.
 *
 * Each generation fills an output container of RVec-like values with
 * point data generated with increments between each RVec and each
 * value within the RVec. Subsequent generations apply a per-frame
 * increment to all point values, making it easy to create per-frame
 * trajectory data where data points are easily identifiable for
 * each frame, atom and dimension.
 *
 * Simple example (using Google Test):
 * \code
   namespace gmx
   {
   namespace test
   {
   TEST(MyTest, SimpleTest)
   {
       constexpr int nframes = 3;
       constexpr int natoms  = 4;
       std::array<std::array<RVec, natoms>, nframes> xs;
       std::for_each(xs.begin(), xs.end(), TrajectoryFrameDataGenerator());

       for (int frame = 0; frame < nframes; ++frame)
       {
           // ... do some testing with xs[frame] ...
       }
   }
   } // namespace test
   } // namespace gmx
 * \endcode
 */
class TrajectoryFrameDataGenerator
{
public:
    //! \brief Offsets for generated trajectory data.
    struct Offset
    {
        /*! \brief Constructor for offset data.
         * \param[in] start Starting offset for frame 0.
         * \param[in] frame Per-frame offset (increment).
         * \param[in] value Per-value offset (increment).
         * \param[in] dim   Per-dimension offset (increment).
         */
        Offset(double start = 0.0, double frame = 1.0, double value = 0.01, double dim = 0.001);

        //!< Accumulated offset up to the current frame.
        double current_;

        //!< Per-frame offset for all values.
        const double frame_;

        //!< Per-value offset for each vector in a frame.
        const double value_;

        //!< Per-dim (axis) offset for each value in a vector.
        const double dim_;
    };

    //! \brief Construct a generator with \p mode and \p offset incrementation.
    TrajectoryFrameDataGenerator(TrajectoryFrameMode mode   = TrajectoryFrameMode::Positive,
                                 const Offset&       offset = {});

    //! \brief Generate the next frame values and write them into \p values.
    void operator()(ArrayRef<BasicVector<float>> values);

    //! \brief Generate the next frame values and write them into \p values.
    void operator()(ArrayRef<BasicVector<double>> values);

    //! \brief Generate the next frame values and write them into \p values.
    // TODO: Remove once we are fully in BasicVector land
    void operator()(ArrayRef<float[3]> values);

    //! \brief Generate the next frame values and write them into \p values.
    // TODO: Remove once we are fully in BasicVector land
    void operator()(ArrayRef<double[3]> values);

private:
    //!< Data generation mode.
    const TrajectoryFrameMode mode_;

    //!< Offsets (increments) for generated data.
    Offset offset_;
};

/*! \brief
 * Generator of per-frame 3x3-matrix data.
 *
 * Each generation fills an output matrix with point data generated
 * with increments between each value. Subsequent generations apply
 * a per-frame increment to all point values, making it easy to create
 * per-frame matrix data where data points are easily identifiable
 * for each frame and matrix element.
 *
 * Simple example (using Google Test):
 * \code
   namespace gmx
   {
   namespace test
   {
   TEST(MyTest, SimpleTest)
   {
       constexpr int nframes = 3;
       std::array<matrix, nframes> boxes;
       std::for_each(boxes.begin(), boxes.end(), MatrixFrameDataGenerator());

       for (int frame = 0; frame < nframes; ++frame)
       {
           // ... do some testing with boxes[frame] ...
       }
   }
   } // namespace test
   } // namespace gmx
 * \endcode
 */
class MatrixFrameDataGenerator
{
public:
    //! \brief Offsets for generated matrix data.
    struct Offset
    {
        /*! \brief Constructor for offset data.
         * \param[in] start Starting offset for frame 0.
         * \param[in] frame Per-frame offset (increment).
         * \param[in] value Per-value offset (increment).
         */
        Offset(double start = 0.0, double frame = 1.0, double value = 0.1);

        //!< Accumulated offset up to the current frame.
        double current_;

        //!< Per-frame offset for all values.
        const double frame_;

        //!< Per-value offset in each frame.
        const double value_;
    };

    //! \brief Construct a generator with \p offset incrementation.
    MatrixFrameDataGenerator(const Offset& offset = {});

    //! \brief Generate the next frame values and write them into \p box.
    void operator()(BasicMatrix3x3<float>& box);

    //! \brief Generate the next frame values and write them into \p box.
    void operator()(BasicMatrix3x3<double>& box);

    //! \brief Generate the next frame values and write them into \p box.
    // TODO: Remove once we are fully in BasicMatrix land
    void operator()(double box[3][3]);

    //! \brief Generate the next frame values and write them into \p box.
    // TODO: Remove once we are fully in BasicMatrix land
    void operator()(float box[3][3]);

private:
    //!< Offsets (increments) for generated data.
    Offset offset_;
};

} // namespace test
} // namespace gmx

#endif // GMX_TESTUTILS_GENERATE_FRAME_DATA_H
