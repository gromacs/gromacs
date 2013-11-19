/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010,2012 The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "smalloc.h"
#include "string2.h"
#include "types/hw_info.h"

#include "gpu_utils.h"
#include "../cuda_tools/cudautils.cuh"
#include "memtestG80_core.h"


#define QUICK_MEM       250 /*!< Amount of memory to be used in quick memtest. */
#define QUICK_TESTS     MOD_20_32BIT | LOGIC_4_ITER_SHMEM | RANDOM_BLOCKS /*!< Bit flag with type of tests
                                                                            to run in quick memtest. */
#define QUICK_ITER      3 /*!< Number of iterations in quick memtest. */

#define FULL_TESTS      0x3FFF /*!<  Bitflag with all test set on for full memetest. */
#define FULL_ITER       25 /*!< Number of iterations in full memtest. */

#define TIMED_TESTS     MOD_20_32BIT | LOGIC_4_ITER_SHMEM | RANDOM_BLOCKS /*!< Bit flag with type of tests to
                                                                            run in time constrained memtest. */

static int cuda_max_device_count = 32; /*! Max number of devices supported by CUDA (for consistency checking).
                                           In reality it 16 with CUDA <=v5.0, but let's stay on the safe side. */

/*! Dummy kernel used for sanity checking. */
__global__ void k_dummy_test(){}


/*! Bit-flags which refer to memtestG80 test types and are used in do_memtest to specify which tests to run. */
enum memtest_G80_test_types {
    MOVING_INVERSIONS_10 =      0x1,
    MOVING_INVERSIONS_RAND =    0x2,
    WALKING_8BIT_M86 =          0x4,
    WALKING_0_8BIT =            0x8,
    WALKING_1_8BIT =            0x10,
    WALKING_0_32BIT =           0x20,
    WALKING_1_32BIT =           0x40,
    RANDOM_BLOCKS =             0x80,
    MOD_20_32BIT =              0x100,
    LOGIC_1_ITER =              0x200,
    LOGIC_4_ITER =              0x400,
    LOGIC_1_ITER_SHMEM =        0x800,
    LOGIC_4_ITER_SHMEM =        0x1000
};


/*! 
  * \brief Runs GPU sanity checks.
  *
  * Runs a series of checks to determine that the given GPU and underlying CUDA
  * driver/runtime functions properly.
  * Returns properties of a device with given ID or the one that has
  * already been initialized earlier in the case if of \dev_id == -1.
  *
  * \param[in]  dev_id      the device ID of the GPU or -1 if the device has already been initialized
  * \param[out] dev_prop    pointer to the structure in which the device properties will be returned
  * \returns                0 if the device looks OK
  *
  * TODO: introduce errors codes and handle errors more smoothly.
  */
static int do_sanity_checks(int dev_id, cudaDeviceProp *dev_prop)
{
    cudaError_t cu_err;
    int         dev_count, id;

    cu_err = cudaGetDeviceCount(&dev_count);
    if (cu_err != cudaSuccess)
    {
       fprintf(stderr, "Error %d while querying device count: %s\n", cu_err,
               cudaGetErrorString(cu_err));
        return -1;
    }

    /* no CUDA compatible device at all */
    if (dev_count == 0)
        return -1;

    /* things might go horribly wrong if cudart is not compatible with the driver */
    if (dev_count < 0 || dev_count > cuda_max_device_count)
        return -1;

    if (dev_id == -1) /* device already selected let's not destroy the context */
    {
        cu_err = cudaGetDevice(&id);
        if (cu_err != cudaSuccess)
        {
            fprintf(stderr, "Error %d while querying device id: %s\n", cu_err,
                    cudaGetErrorString(cu_err));
            return -1;
        }
    }
    else
    {
        id = dev_id;
        if (id > dev_count - 1) /* pfff there's no such device */
        {
            fprintf(stderr, "The requested device with id %d does not seem to exist (device count=%d)\n",
                    dev_id, dev_count);
            return -1;
        }
    }

    memset(dev_prop, 0, sizeof(cudaDeviceProp));
    cu_err = cudaGetDeviceProperties(dev_prop, id);
    if (cu_err != cudaSuccess)
    {
        fprintf(stderr, "Error %d while querying device properties: %s\n", cu_err,
                cudaGetErrorString(cu_err));
        return -1;
    }

    /* both major & minor is 9999 if no CUDA capable devices are present */
    if (dev_prop->major == 9999 && dev_prop->minor == 9999)
        return -1;
    /* we don't care about emulation mode */
    if (dev_prop->major == 0)
        return -1;

    if (id != -1)
    {
        cu_err = cudaSetDevice(id);
        if (cu_err != cudaSuccess)
        {
            fprintf(stderr, "Error %d while switching to device #%d: %s\n",
                    cu_err, id, cudaGetErrorString(cu_err));
            return -1;
        }
    }

    /* try to execute a dummy kernel */
    k_dummy_test<<<1, 512>>>();
    if (cudaThreadSynchronize() != cudaSuccess)
    {
        return -1;
    }

    /* destroy context if we created one */
    if (id != -1)
    {
#if CUDA_VERSION < 4000
        cu_err = cudaThreadExit();
        CU_RET_ERR(cu_err, "cudaThreadExit failed");
#else
        cu_err = cudaDeviceReset();
        CU_RET_ERR(cu_err, "cudaDeviceReset failed");
#endif
    }

    return 0;
}


/*!
 * \brief Runs a set of memory tests specified by the given bit-flags.
 * Tries to allocate and do the test on \p megs Mb memory or
 * the greatest amount that can be allocated (>10Mb).
 * In case if an error is detected it stops without finishing the remaining
 * steps/iterations and returns greater then zero value.
 * In case of other errors (e.g. kernel launch errors, device querying errors)
 * -1 is returned.
 *
 * \param[in] which_tests   variable with bit-flags of the requested tests
 * \param[in] megs          amount of memory that will be tested in MB
 * \param[in] iter          number of iterations
 * \returns                 0 if no error was detected, otherwise >0
 */
static int do_memtest(unsigned int which_tests, int megs, int iter)
{
    memtestState    tester;
    int             i;
    uint            err_count; //, err_iter;

    // no parameter check as this fn won't be called externally

    // let's try to allocate the mem
    while (!tester.allocate(megs) && (megs - 10 > 0))
        { megs -= 10; tester.deallocate(); }

    if (megs <= 10)
    {
        fprintf(stderr, "Unable to allocate GPU memory!\n");
        return -1;
    }

    // clear the first 18 bits
    which_tests &= 0x3FFF;
    for (i = 0; i < iter; i++)
    {
        // Moving Inversions (ones and zeros)
        if ((MOVING_INVERSIONS_10 & which_tests) == MOVING_INVERSIONS_10)
        {
            tester.gpuMovingInversionsOnesZeros(err_count);
            if (err_count > 0)
                return MOVING_INVERSIONS_10;
        }
        // Moving Inversions (random)
        if ((MOVING_INVERSIONS_RAND & which_tests) == MOVING_INVERSIONS_RAND)
        {
            tester.gpuMovingInversionsRandom(err_count);
            if (err_count > 0)
                return MOVING_INVERSIONS_RAND;
        }
       // Memtest86 Walking 8-bit
        if ((WALKING_8BIT_M86 & which_tests) == WALKING_8BIT_M86)
        {
            for (uint shift = 0; shift < 8; shift++)
            {
                tester.gpuWalking8BitM86(err_count, shift);
                if (err_count > 0)
                    return WALKING_8BIT_M86;
            }
      }
        // True Walking zeros (8-bit)
        if ((WALKING_0_8BIT & which_tests) == WALKING_0_8BIT)
        {
            for (uint shift = 0; shift < 8; shift++)
            {
                tester.gpuWalking8Bit(err_count, false, shift);
                if (err_count > 0)
                    return WALKING_0_8BIT;
            }
        }
        // True Walking ones (8-bit)
        if ((WALKING_1_8BIT & which_tests) == WALKING_1_8BIT)
        {
            for (uint shift = 0; shift < 8; shift++)
            {
                tester.gpuWalking8Bit(err_count, true, shift);
                if (err_count > 0)
                    return WALKING_1_8BIT;
            }
        }
        // Memtest86 Walking zeros (32-bit)
        if ((WALKING_0_32BIT & which_tests) == WALKING_0_32BIT)
        {
            for (uint shift = 0; shift < 32; shift++)
            {
                tester.gpuWalking32Bit(err_count, false, shift);
                if (err_count > 0)
                    return WALKING_0_32BIT;
            }
        }
       // Memtest86 Walking ones (32-bit)
        if ((WALKING_1_32BIT & which_tests) == WALKING_1_32BIT)
        {
            for (uint shift = 0; shift < 32; shift++)
            {
                tester.gpuWalking32Bit(err_count, true, shift);
                if (err_count > 0)
                    return WALKING_1_32BIT;
            }
       }
        // Random blocks
        if ((RANDOM_BLOCKS & which_tests) == RANDOM_BLOCKS)
        {
            tester.gpuRandomBlocks(err_count,rand());
            if (err_count > 0)
                return RANDOM_BLOCKS;

        }

        // Memtest86 Modulo-20
        if ((MOD_20_32BIT & which_tests) == MOD_20_32BIT)
        {
            for (uint shift = 0; shift < 20; shift++)
            {
                tester.gpuModuloX(err_count, shift, rand(), 20, 2);
                if (err_count > 0)
                    return MOD_20_32BIT;
            }
        }
        // Logic (one iteration)
        if ((LOGIC_1_ITER & which_tests) == LOGIC_1_ITER)
        {
            tester.gpuShortLCG0(err_count,1);
            if (err_count > 0)
                return LOGIC_1_ITER;
        }
        // Logic (4 iterations)
        if ((LOGIC_4_ITER & which_tests) == LOGIC_4_ITER)
        {
            tester.gpuShortLCG0(err_count,4);
            if (err_count > 0)
                return LOGIC_4_ITER;

        }
        // Logic (shared memory, one iteration)
        if ((LOGIC_1_ITER_SHMEM & which_tests) == LOGIC_1_ITER_SHMEM)
        {
            tester.gpuShortLCG0Shmem(err_count,1);
            if (err_count > 0)
                return LOGIC_1_ITER_SHMEM;
        }
        // Logic (shared-memory, 4 iterations)
        if ((LOGIC_4_ITER_SHMEM & which_tests) == LOGIC_4_ITER_SHMEM)
        {
            tester.gpuShortLCG0Shmem(err_count,4);
            if (err_count > 0)
                return LOGIC_4_ITER_SHMEM;
        }
    }

    tester.deallocate();
    return err_count;
}

/*! \brief Runs a quick memory test and returns 0 in case if no error is detected.
 * If an error is detected it stops before completing the test and returns a
 * value greater then 0. In case of other errors (e.g. kernel launch errors,
 * device querying errors) -1 is returned.
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has already been selected
 * \returns             0 if no error was detected, otherwise >0
 */
int do_quick_memtest(int dev_id)
{
    cudaDeviceProp  dev_prop;
    int             devmem, res, time=0;

    if (debug) { time = getTimeMilliseconds(); }

    if (do_sanity_checks(dev_id, &dev_prop) != 0)
    {
        // something went wrong
        return -1;
    }

    if (debug)
    {
        devmem = dev_prop.totalGlobalMem/(1024*1024); // in MiB
        fprintf(debug, ">> Running QUICK memtests on %d MiB (out of total %d MiB), %d iterations\n",
            QUICK_MEM, devmem, QUICK_ITER);
    }

    res = do_memtest(QUICK_TESTS, QUICK_MEM, QUICK_ITER);

    if (debug)
    {
        fprintf(debug, "Q-RES = %d\n", res);
        fprintf(debug, "Q-runtime: %d ms\n", getTimeMilliseconds() - time);
    }

    /* destroy context only if we created it */
    if (dev_id !=-1) cudaThreadExit();
    return res;
}

/*! \brief Runs a full memory test and returns 0 in case if no error is detected.
 * If an error is detected  it stops before completing the test and returns a
 * value greater then 0. In case of other errors (e.g. kernel launch errors,
 * device querying errors) -1 is returned.
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has already been selected
 * \returns             0 if no error was detected, otherwise >0
 */

int do_full_memtest(int dev_id)
{
    cudaDeviceProp  dev_prop;
    int             devmem, res, time=0;

    if (debug) { time = getTimeMilliseconds(); }

    if (do_sanity_checks(dev_id, &dev_prop) != 0)
    {
        // something went wrong
        return -1;
    }

    devmem = dev_prop.totalGlobalMem/(1024*1024); // in MiB

    if (debug) 
    { 
        fprintf(debug, ">> Running FULL memtests on %d MiB (out of total %d MiB), %d iterations\n",
            devmem, devmem, FULL_ITER); 
    }

    /* do all test on the entire memory */
    res = do_memtest(FULL_TESTS, devmem, FULL_ITER);

    if (debug)
    {
        fprintf(debug, "F-RES = %d\n", res);
        fprintf(debug, "F-runtime: %d ms\n", getTimeMilliseconds() - time);
    }

    /* destroy context only if we created it */
    if (dev_id != -1) cudaThreadExit();
    return res;
}

/*! \brief Runs a time constrained memory test and returns 0 in case if no error is detected.
 * If an error is detected it stops before completing the test and returns a value greater
 * than zero. In case of other errors (e.g. kernel launch errors, device querying errors) -1
 * is returned. Note, that test iterations are not interrupted therefor the total runtime of
 * the test will always be multipple of one iteration's runtime.
 *
 * \param[in] dev_id        the device id of the GPU or -1 if the device has laredy been selected
 * \param[in] time_constr   the time limit of the testing
 * \returns                 0 if no error was detected, otherwise >0
 */
int do_timed_memtest(int dev_id, int time_constr)
{
    cudaDeviceProp  dev_prop;
    int             devmem, res=0, time=0, startt;

    if (debug) { time = getTimeMilliseconds(); }

    time_constr *= 1000;  /* convert to ms for convenience */
    startt = getTimeMilliseconds();

    if (do_sanity_checks(dev_id, &dev_prop) != 0)
    {
        // something went wrong
        return -1;
    }

    devmem = dev_prop.totalGlobalMem/(1024*1024); // in MiB

    if (debug) 
    { 
        fprintf(debug, ">> Running time constrained memtests on %d MiB (out of total %d MiB), time limit of %d s \n",
        devmem, devmem, time_constr); 
    }

    /* do the TIMED_TESTS set, one step at a time on the entire memory 
       that can be allocated, and stop when the given time is exceeded */
    while ( ((int)getTimeMilliseconds() - startt) < time_constr)
    {        
        res = do_memtest(TIMED_TESTS, devmem, 1);
        if (res != 0) break;
    }

    if (debug)
    {
        fprintf(debug, "T-RES = %d\n", res);
        fprintf(debug, "T-runtime: %d ms\n", getTimeMilliseconds() - time);
    }

    /* destroy context only if we created it */
    if (dev_id != -1) cudaThreadExit();
    return res;
}

/*! \brief Initializes the GPU with the given index.
 *
 * The varible \mygpu is the index of the GPU to initialize in the
 * gpu_info.cuda_dev array.
 *
 * \param[in]  mygpu        index of the GPU to initialize
 * \param[out] result_str   the message related to the error that occurred
 *                          during the initialization (if there was any).
 * \param[in] gpu_info      GPU info of all detected devices in the system.
 * \param[in] gpu_opt       options for using the GPUs in gpu_info
 * \returns                 true if no error occurs during initialization.
 */
gmx_bool init_gpu(int mygpu, char *result_str,
                  const gmx_gpu_info_t *gpu_info,
                  const gmx_gpu_opt_t *gpu_opt)
{
    cudaError_t stat;
    char sbuf[STRLEN];
    int gpuid;

    assert(gpu_info);
    assert(result_str);

    if (mygpu < 0 || mygpu >= gpu_opt->ncuda_dev_use)
    {
        sprintf(sbuf, "Trying to initialize an inexistent GPU: "
                "there are %d %s-selected GPU(s), but #%d was requested.",
                 gpu_opt->ncuda_dev_use, gpu_opt->bUserSet ? "user" : "auto", mygpu);
        gmx_incons(sbuf);
    }

    gpuid = gpu_info->cuda_dev[gpu_opt->cuda_dev_use[mygpu]].id;

    stat = cudaSetDevice(gpuid);
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);

    if (debug)
    {
        fprintf(stderr, "Initialized GPU ID #%d: %s\n", gpuid, gpu_info->cuda_dev[gpuid].prop.name);
    }

    return (stat == cudaSuccess);
}

/*! \brief Frees up the CUDA GPU used by the active context at the time of calling.
 *
 * The context is explicitly destroyed and therefore all data uploaded to the GPU
 * is lost. This should only be called when none of this data is required anymore.
 *
 * \param[out] result_str   the message related to the error that occurred
 *                          during the initialization (if there was any).
 * \returns                 true if no error occurs during the freeing.
 */
gmx_bool free_gpu(char *result_str)
{
    cudaError_t stat;

    assert(result_str);

    if (debug)
    {
        int gpuid;
        stat = cudaGetDevice(&gpuid);
        CU_RET_ERR(stat, "cudaGetDevice failed");
        fprintf(stderr, "Cleaning up context on GPU ID #%d\n", gpuid);
    }

#if CUDA_VERSION < 4000
    stat = cudaThreadExit();
#else
    stat = cudaDeviceReset();
#endif
    strncpy(result_str, cudaGetErrorString(stat), STRLEN);

    return (stat == cudaSuccess);
}

/*! \brief Returns true if the gpu characterized by the device properties is
 *  supported by the native gpu acceleration.
 *
 * \param[in] dev_prop  the CUDA device properties of the gpus to test.
 * \returns             true if the GPU properties passed indicate a compatible
 *                      GPU, otherwise false.
 */
static bool is_gmx_supported_gpu(const cudaDeviceProp *dev_prop)
{
    return (dev_prop->major >= 2);
}

/*! \brief Helper function that checks whether a given GPU status indicates compatible GPU.
 *
 * \param[in] stat  GPU status.
 * \returns         true if the provided status is egpuCompatible, otherwise false.
 */
static bool is_compatible_gpu(int stat)
{
    return (stat == egpuCompatible);
}

/*! \brief Checks if a GPU with a given ID is supported by the native GROMACS acceleration.
 *
 *  Returns a status value which indicates compatibility or one of the following
 *  errors: incompatibility, insistence, or insanity (=unexpected behavior).
 *  It also returns the respective device's properties in \dev_prop (if applicable).
 *
 *  \param[in]  dev_id   the ID of the GPU to check.
 *  \param[out] dev_prop the CUDA device properties of the device checked.
 *  \returns             the status of the requested device
 */
static int is_gmx_supported_gpu_id(int dev_id, cudaDeviceProp *dev_prop)
{
    cudaError_t stat;
    int         ndev;

    stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        return egpuInsane;
    }

    if (dev_id > ndev - 1)
    {
        return egpuNonexistent;
    }

    /* TODO: currently we do not make a distinction between the type of errors
     * that can appear during sanity checks. This needs to be improved, e.g if
     * the dummy test kernel fails to execute with a "device busy message" we
     * should appropriately report that the device is busy instead of insane.
     */
    if (do_sanity_checks(dev_id, dev_prop) == 0)
    {
        if (is_gmx_supported_gpu(dev_prop))
        {
            return egpuCompatible;
        }
        else
        {
            return egpuIncompatible;
        }
    }
    else
    {
        return egpuInsane;
    }
}


/*! \brief Detect all NVIDIA GPUs in the system.
 *
 *  Will detect every NVIDIA GPU supported by the device driver in use. Also
 *  check for the compatibility of each and fill the gpu_info->cuda_dev array
 *  with the required information on each the device: ID, device properties,
 *  status.
 *
 *  \param[in] gpu_info    pointer to structure holding GPU information.
 *  \param[out] err_str    The error message of any CUDA API error that caused
 *                         the detection to fail (if there was any). The memory
 *                         the pointer points to should be managed externally.
 *  \returns               non-zero if the detection encountered a failure, zero otherwise.
 */
int detect_cuda_gpus(gmx_gpu_info_t *gpu_info, char *err_str)
{
    int             i, ndev, checkres, retval;
    cudaError_t     stat;
    cudaDeviceProp  prop;
    cuda_dev_info_t *devs;

    assert(gpu_info);
    assert(err_str);

    gpu_info->ncuda_dev_compatible = 0;

    ndev    = 0;
    devs    = NULL;

    stat = cudaGetDeviceCount(&ndev);
    if (stat != cudaSuccess)
    {
        const char *s;

        /* cudaGetDeviceCount failed which means that there is something
         * wrong with the machine: driver-runtime mismatch, all GPUs being
         * busy in exclusive mode, or some other condition which should
         * result in us issuing a warning a falling back to CPUs. */
        retval = -1;
        s = cudaGetErrorString(stat);
        strncpy(err_str, s, STRLEN*sizeof(err_str[0]));
    }
    else
    {
        snew(devs, ndev);
        for (i = 0; i < ndev; i++)
        {
            checkres = is_gmx_supported_gpu_id(i, &prop);

            devs[i].id   = i;
            devs[i].prop = prop;
            devs[i].stat = checkres;

            if (checkres == egpuCompatible)
            {
                gpu_info->ncuda_dev_compatible++;
            }
        }
        retval = 0;
    }

    gpu_info->ncuda_dev = ndev;
    gpu_info->cuda_dev  = devs;

    return retval;
}

/*! \brief Select the GPUs compatible with the native GROMACS acceleration.
 *
 * This function selects the compatible gpus and initializes
 * gpu_info->cuda_dev_use and gpu_info->ncuda_dev_use.
 *
 * Given the list of GPUs available in the system check each device in
 * gpu_info->cuda_dev and place the indices of the compatible GPUs into
 * cuda_dev_use with this marking the respective GPUs as "available for use."
 * Note that \detect_cuda_gpus must have been called before.
 *
 * \param[in]     gpu_info    pointer to structure holding GPU information
 * \param[in,out] gpu_opt     pointer to structure holding GPU options
 */
void pick_compatible_gpus(const gmx_gpu_info_t *gpu_info,
                          gmx_gpu_opt_t *gpu_opt)
{
    int i, ncompat;
    int *compat;

    assert(gpu_info);
    /* cuda_dev/ncuda_dev have to be either NULL/0 or not (NULL/0) */
    assert((gpu_info->ncuda_dev != 0 ? 0 : 1) ^ (gpu_info->cuda_dev == NULL ? 0 : 1));

    snew(compat, gpu_info->ncuda_dev);
    ncompat = 0;
    for (i = 0; i < gpu_info->ncuda_dev; i++)
    {
        if (is_compatible_gpu(gpu_info->cuda_dev[i].stat))
        {
            ncompat++;
            compat[ncompat - 1] = i;
        }
    }

    gpu_opt->ncuda_dev_use = ncompat;
    snew(gpu_opt->cuda_dev_use, ncompat);
    memcpy(gpu_opt->cuda_dev_use, compat, ncompat*sizeof(*compat));
    sfree(compat);
}

/*! \brief Check the existence/compatibility of a set of GPUs specified by their device IDs.
 *
 * Given the a list of gpu->ncuda_dev_use GPU device IDs stored in
 * gpu_opt->cuda_dev_use check the existence and compatibility
 * of the respective GPUs. Also provide the caller with an array containing
 * the result of checks in \checkres.
 *
 * \param[out]  checkres    check result for each ID passed in \requested_devs
 * \param[in]   gpu_info    pointer to structure holding GPU information
 * \param[out]  gpu_opt     pointer to structure holding GPU options
 * \returns                 TRUE if every the requested GPUs are compatible
 */
gmx_bool check_selected_cuda_gpus(int *checkres,
                                  const gmx_gpu_info_t *gpu_info,
                                  gmx_gpu_opt_t *gpu_opt)
{
    int i, id;
    bool bAllOk;

    assert(checkres);
    assert(gpu_info);
    assert(gpu_opt->ncuda_dev_use >= 0);

    if (gpu_opt->ncuda_dev_use == 0)
    {
        return TRUE;
    }

    assert(gpu_opt->cuda_dev_use);

    /* we will assume that all GPUs requested are valid IDs,
       otherwise we'll bail anyways */

    bAllOk = true;
    for (i = 0; i < gpu_opt->ncuda_dev_use; i++)
    {
        id = gpu_opt->cuda_dev_use[i];

        /* devices are stored in increasing order of IDs in cuda_dev */
        gpu_opt->cuda_dev_use[i] = id;

        checkres[i] = (id >= gpu_info->ncuda_dev) ?
            egpuNonexistent : gpu_info->cuda_dev[id].stat;

        bAllOk = bAllOk && is_compatible_gpu(checkres[i]);
    }

    return bAllOk;
}

/*! \brief Frees the cuda_dev and cuda_dev_use array fields of \gpu_info.
 *
 * \param[in]    gpu_info    pointer to structure holding GPU information
 */
void free_gpu_info(const gmx_gpu_info_t *gpu_info)
{
    if (gpu_info == NULL)
    {
        return;
    }

    sfree(gpu_info->cuda_dev);
}

/*! \brief Formats and returns a device information string for a given GPU.
 *
 * Given an index *directly* into the array of available GPUs (cuda_dev)
 * returns a formatted info string for the respective GPU which includes
 * ID, name, compute capability, and detection status.
 *
 * \param[out]  s           pointer to output string (has to be allocated externally)
 * \param[in]   gpu_info    pointer to structure holding GPU information
 * \param[in]   index       an index *directly* into the array of available GPUs
 */
void get_gpu_device_info_string(char *s, const gmx_gpu_info_t *gpu_info, int index)
{
    assert(s);
    assert(gpu_info);

    if (index < 0 && index >= gpu_info->ncuda_dev)
    {
        return;
    }

    cuda_dev_info_t *dinfo = &gpu_info->cuda_dev[index];

    bool bGpuExists =
        dinfo->stat == egpuCompatible ||
        dinfo->stat == egpuIncompatible;

    if (!bGpuExists)
    {
        sprintf(s, "#%d: %s, stat: %s",
                dinfo->id, "N/A",
                gpu_detect_res_str[dinfo->stat]);
    }
    else
    {
        sprintf(s, "#%d: NVIDIA %s, compute cap.: %d.%d, ECC: %3s, stat: %s",
                dinfo->id, dinfo->prop.name,
                dinfo->prop.major, dinfo->prop.minor,
                dinfo->prop.ECCEnabled ? "yes" : " no",
                gpu_detect_res_str[dinfo->stat]);
    }
}

/*! \brief Returns the device ID of the GPU with a given index into the array of used GPUs.
 *
 * Getter function which, given an index into the array of GPUs in use
 * (cuda_dev_use) -- typically a tMPI/MPI rank --, returns the device ID of the
 * respective CUDA GPU.
 *
 * \param[in]    gpu_info   pointer to structure holding GPU information
 * \param[in]    gpu_opt    pointer to structure holding GPU options
 * \param[in]    idx        index into the array of used GPUs
 * \returns                 device ID of the requested GPU
 */
int get_gpu_device_id(const gmx_gpu_info_t *gpu_info,
                      const gmx_gpu_opt_t *gpu_opt,
                      int idx)
{
    assert(gpu_info);
    assert(gpu_opt);
    assert(idx >= 0 && idx < gpu_opt->ncuda_dev_use);

    return gpu_info->cuda_dev[gpu_opt->cuda_dev_use[idx]].id;
}

/*! \brief Returns the device ID of the GPU currently in use.
 *
 * The GPU used is the one that is active at the time of the call in the active context.
 *
 * \param[in]    gpu_info   pointer to structure holding GPU information
 * \returns                 device ID of the GPU in use at the time of the call
 */
int get_current_gpu_device_id(void)
{
    int gpuid;
    CU_RET_ERR(cudaGetDevice(&gpuid), "cudaGetDevice failed");

    return gpuid;
}

/*! \brief Returns the size of the cuda_dev_info struct.
 *
 * The size of cuda_dev_info can be used for allocation and communication.
 *
 * \returns                 size in bytes of cuda_dev_info
 */
size_t sizeof_cuda_dev_info(void)
{
    return sizeof(cuda_dev_info);
}
