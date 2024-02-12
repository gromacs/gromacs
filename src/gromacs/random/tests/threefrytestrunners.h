/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Runners for the ThreeFry random engine.
 *
 * The test runners are used to evaluate different implementations of the
 * same test case against the same reference data.
 *
 * The abstract interface class is used to unify the interfaces of the CPU
 * and GPU implementations. This way, the same test data can be used in
 * different implementations that inherit the interfaces of the parent class.
 * The host and device implementations differ in the host-device data
 * transfers.
 *
 * To ensure full control over the key, to be able to compare to reference
 * data, the whole 128 bit key is specified using two 64-bit key values.
 *
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_TESTS_THREEFRYTESTRUNNERS_H
#define GMX_RANDOM_TESTS_THREEFRYTESTRUNNERS_H

#include "config.h"

#include "gromacs/random/threefry.h"

/*
 * ThreeFry is supported on CUDA, HIP and SYCL
 */
#define GPU_THREEFRY_SUPPORTED (GMX_GPU_CUDA || GMX_GPU_HIP || GMX_GPU_SYCL)

#if GPU_THREEFRY_SUPPORTED
#    include "gromacs/gpu_utils/devicebuffer_datatype.h"
#endif

#include "testutils/test_device.h"

namespace gmx
{

namespace test
{

/* \brief ThreeFry test runner interface. */
template<unsigned int rounds, unsigned int internalCounterBits>
class IThreeFryGeneralTestRunner
{
public:
    //! Virtual destructor
    virtual ~IThreeFryGeneralTestRunner() {}

    virtual std::string hardwareDescription() = 0;

    virtual void     restartRng(uint64_t ctr0, uint64_t ctr1) = 0;
    virtual uint64_t nextRng()                                = 0;
};

template<unsigned int rounds, unsigned int internalCounterBits>
class ThreeFryGeneralHostTestRunner : public IThreeFryGeneralTestRunner<rounds, internalCounterBits>
{
public:
    //! Constructor
    ThreeFryGeneralHostTestRunner(uint64_t key0, uint64_t key1) :
        rng_(gmx::ThreeFry2x64General<rounds, internalCounterBits>(key0, key1))
    {
    }

    ~ThreeFryGeneralHostTestRunner() override {}

    std::string hardwareDescription() override { return "CPU"; }

    void     restartRng(uint64_t ctr0, uint64_t ctr1) override { rng_.restart(ctr0, ctr1); }
    uint64_t nextRng() override { return (rng_)(); }

private:
    gmx::ThreeFry2x64General<rounds, internalCounterBits> rng_;
};

template<unsigned int internalCounterBits>
class ThreeFryHostTestRunner : public ThreeFryGeneralHostTestRunner<20, internalCounterBits>
{
public:
    ThreeFryHostTestRunner(uint64_t key0, uint64_t key1) :
        ThreeFryGeneralHostTestRunner<20, internalCounterBits>(key0, key1)
    {
    }
};

template<unsigned int internalCounterBits>
class ThreeFryFastHostTestRunner : public ThreeFryGeneralHostTestRunner<13, internalCounterBits>
{
public:
    ThreeFryFastHostTestRunner(uint64_t key0, uint64_t key1) :
        ThreeFryGeneralHostTestRunner<13, internalCounterBits>(key0, key1)
    {
    }
};


template<unsigned int rounds, unsigned int internalCounterBits>
class ThreeFryGeneralDeviceTestRunner : public IThreeFryGeneralTestRunner<rounds, internalCounterBits>
{
public:
    ThreeFryGeneralDeviceTestRunner(const TestDevice& testDevice, uint64_t key0, uint64_t key1);
    ~ThreeFryGeneralDeviceTestRunner() override;

    std::string hardwareDescription() override { return testDevice_.description(); }

    void     restartRng(uint64_t ctr0, uint64_t ctr1) override;
    uint64_t nextRng() override;


private:
    //! Test device to be used in the runner.
    const TestDevice& testDevice_;

#if GPU_THREEFRY_SUPPORTED
    DeviceBuffer<gmx::ThreeFry2x64General<rounds, internalCounterBits>> rng_;
    DeviceBuffer<uint64_t>                                              d_result;
#endif
};

template<unsigned int internalCounterBits>
class ThreeFryDeviceTestRunner : public ThreeFryGeneralDeviceTestRunner<20, internalCounterBits>
{
public:
    ThreeFryDeviceTestRunner(const TestDevice& testDevice, uint64_t key0, uint64_t key1) :
        ThreeFryGeneralDeviceTestRunner<20, internalCounterBits>(testDevice, key0, key1)
    {
    }
};

template<unsigned int internalCounterBits>
class ThreeFryFastDeviceTestRunner : public ThreeFryGeneralDeviceTestRunner<13, internalCounterBits>
{
public:
    ThreeFryFastDeviceTestRunner(const TestDevice& testDevice, uint64_t key0, uint64_t key1) :
        ThreeFryGeneralDeviceTestRunner<13, internalCounterBits>(testDevice, key0, key1)
    {
    }
};

/* Disable the "explicit template instantiation 'ThreeFryGeneralDeviceTestRunner<...>' will emit a vtable in every
 * translation unit [-Wweak-template-vtables]" warning.
 * This warning is not used anymore in clang 19 and will be removed in future releases.
 */
#ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wweak-template-vtables"
#endif
extern template class ThreeFryGeneralDeviceTestRunner<13, 0>;
extern template class ThreeFryGeneralDeviceTestRunner<20, 0>;
extern template class ThreeFryGeneralDeviceTestRunner<40, 0>;

} // namespace test
} // namespace gmx

#endif // GMX_RANDOM_TESTS_THREEFRYTESTRUNNERS_H
