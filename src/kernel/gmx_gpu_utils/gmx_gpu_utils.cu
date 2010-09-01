/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team,
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
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#include <stdio.h>
#include <stdlib.h>

#include "cuda.h"
#include "cuda_runtime_api.h"

#include "memtestG80_core.h"

/*! \cond  TEST */
#ifdef _DEBUG_
#undef _DEBUG_
#endif
#define _DEBUG_           0

#if _DEBUG_ >= 1
#define debug stderr 
#define DUPME(msg) printf("---> %s\n", msg);
#else
#define DUPME(msg) ;
#endif
/*! \endcond TEST*/

#if _DEBUG_ == 0/* no gromacs utils in debug mode */
#include "gmx_fatal.h"
#include "string2.h"
#endif

#define QUICK_MEM       250 /*!< Amount of memory to be used in quick memtest. */
#define QUICK_TESTS     MOD_20_32BIT | LOGIC_4_ITER_SHMEM | RANDOM_BLOCKS /*!< Bitflag with type of tests 
                                                                            to run in quick memtest. */
#define QUICK_ITER      3 /*!< Number of iterations in quick memtest. */

#define FULL_TESTS      0x3FFF /*!<  Bitflag with all test set on for full memetest. */
#define FULL_ITER       25 /*!< Number of iterations in full memtest. */

#define TIMED_TESTS     MOD_20_32BIT | LOGIC_4_ITER_SHMEM | RANDOM_BLOCKS /*!< Bitflag with type of tests to 
                                                                            run in time constrained memtest. */

/*! Number of supported GPUs */
#define NB_GPUS (sizeof(SupportedGPUs)/sizeof(SupportedGPUs[0]))

/*
TODO add proper gromacs logging?
*/

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

// TODO put this list into an external file and include it so that the list is easily accessible
/*! List of supported GPUs. */
static const char * const SupportedGPUs[] = {
    /* GT400 */
    "Geforce GTX 480",
    "Geforce GTX 470",
    "Geforce GTX 465",
    "Geforce GTX 460",

    "Tesla C2070",
    "Tesla C2050",
    "Tesla S2070",
    "Tesla S2050",
    "Tesla M2070",
    "Tesla M2050",

    "Quadro 5000",
    "Quadro 6000",

    /* GT200 */
    "Geforce GTX 295",
    "Geforce GTX 285",
    "Geforce GTX 280",
    "Geforce GTX 275",
    "Geforce GTX 260",
    "GeForce GTS 250",
    "GeForce GTS 150",

    "GeForce GTX 285M",
    "GeForce GTX 280M",

    "Tesla S1070",
    "Tesla C1060",
    "Tesla M1060",

    "Quadro FX 5800",
    "Quadro FX 4800",
    "Quadro CX",
    "Quadro Plex 2200 D2",
    "Quadro Plex 2200 S4",

    /* G90 */
    "GeForce 9800 G", /* GX2, GTX, GTX+, GT */
    "GeForce 9800M GTX",

    "Quadro FX 4700",
    "Quadro Plex 2100 D4"
};

/*! \cond  TEST */
#ifndef _string2_h
/* debug functions, see @the end */
void ltrim (char *);
void rtrim (char *);
void trim  (char *);
int gmx_strncasecmp(const char*, const char*, int);
#endif
/*! \endcond  TEST */


/*! 
  * \brief Runs GPU sanity checks.
  * Returnes properties of a device with given id or the one that has
  * already been initialized earlier in the case if of dev_id == -1.
  *
  * \param[in] dev_id       the device id of the GPU or -1 if the device has laredy been selected
  * \param[out] dev_prop    pointer to the structure in which the device properties will be returned
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
    if (dev_count < 0 || dev_count > 20)
        return -1;

    if (dev_id == -1) /* device already selected let's do not destroy the context */
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

    if ((dev_id != -1) && (cu_err = cudaSetDevice(dev_id)) != cudaSuccess)
    {
        fprintf(stderr, "Error %d while switching to device #%d: %s\n", cu_err, dev_id,
                cudaGetErrorString(cu_err));
        return -1;
    }

    return 0;
}

/*! 
 * \brief Checks whether the GPU with the given name is supported.
 * 
 * \param[in] gpu_name  the name of the CUDA device
 * \returns             1 if the device is supported, otherwise 0
 */
static int is_supported_gpu_n(char *gpuName)
{
    size_t i;
    for (i = 0; i < NB_GPUS; i++)
    {
        trim(gpuName);
        if (gmx_strncasecmp(gpuName, SupportedGPUs[i], strlen(SupportedGPUs[i])) == 0)
            return 1;
    }
    return 0;
}

/*! \brief Checks whether the GPU with the given device id is supported. 
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has laredy been selected
 * \param[out] gpu_name Set to contain the name of the CUDA device, if NULL passed, no device name is set. 
 * \returns             1 if the device is supported, otherwise 0
 */
int is_supported_cuda_gpu(int dev_id, char *gpu_name)
{
    cudaDeviceProp dev_prop;

    if (debug) fprintf(debug, "Checking compatibility with device #%d, %s\n", dev_id, gpu_name);

    if (do_sanity_checks(dev_id, &dev_prop) != 0)
        return -1;

    if (gpu_name != NULL)
    { 
        strcpy(gpu_name, dev_prop.name);
    }
    return is_supported_gpu_n(dev_prop.name);
}


/*!
 * \brief Runs a set of memory tests specified by the given bit-flags.
 * Tries to allocate and do the test on \p megs Mb memory or 
 * the greatest amount that can be allocated (>10Mb).
 * In case if an error is detected it stops without finishing the remainings 
 * steps/iterations and returns greater then zero value.  
 * In case of other errors (e.g. kernel launch errors, device querying erros) 
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
 * device querying erros) -1 is returned.
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has laredy been selected
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
 * device querying erros) -1 is returned.
 *
 * \param[in] dev_id    the device id of the GPU or -1 if the device has laredy been selected
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
 * than zero. In case of other errors (e.g. kernel launch errors, device querying erros) -1 
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

/*! \cond TEST */

/*******************************************************
 * The code below is for testing purposes. */
int do_custom_memtest(int dev_id)
{
    cudaDeviceProp  dev_prop;
    int             mem2test, /*devmem,*/ res;
//    memtestState    tester;
//    double          bandwidth;

#if _DEBUG_ >= 1
    int time = getTimeMilliseconds();
#endif

   if (do_sanity_checks(dev_id, &dev_prop) != 0)
        return -1;

//    if ((res=tester.allocate(100))==0)
//        printf("alloc failed\n");
//    printf("alloc res = %d\n", res);
//    res = tester.gpuMemoryBandwidth(bandwidth, tester.size(), 10);
//    printf("Bandwidth on %d (res %d)= %5.2f\n", tester.size(), res, bandwidth);
//    tester.deallocate();

//    devmem   = dev_prop.totalGlobalMem/(1024*1024); // in MiB
    mem2test = 80;

#if _DEBUG_ >= 1
    printf(">> Running CUSTOM memtests [%x] on %d MiB, %d iterations\n",
        QUICK_TESTS, mem2test, 1);
#endif

    res = do_memtest(QUICK_TESTS, mem2test, 1);
    cudaThreadExit();

#if _DEBUG_ >= 1
    printf("C-RES = %d\n", res);
    printf("C-runtime: %d ms\n", getTimeMilliseconds() - time);
#endif
    return res;
}

#if _DEBUG_ > 1
/*!
 * Only for debugging purposes, compile with:
 * nvcc -DLINUX -D_DEBUG_=2  -L -O  -Xcompiler -Wall memtestG80_core.o gmx_gpu_utils.cu  -o gmx_gpu_utils_test
 */
int main( int argc, char** argv)
{
    int dev_id = 0;
    char msg[100];
    sprintf(msg, "Device #%d supported: ", dev_id);
    switch (is_supported_cuda_gpu(dev_id, NULL))
    {
        case -1: strcat(msg, "error occured"); break;
        case  0: strcat(msg, "no"); break;
        case  1: strcat(msg, "yes"); break;
        default: strcat(msg, "\nhmmm, you should not see this!");
    }
    printf("%s\n", msg);

    printf("Doing memtest.\n");
    printf("quick memtest result: %d\n", do_quick_memtest(dev_id));
    printf("timed memtest result: %d\n", do_timed_memtest(dev_id, 15));
    printf("full memtest result: %d\n", do_full_memtest(dev_id));
    return 0;
}
#endif


#ifndef _string2_h
#include <string.h>
/* 
    Functions only used if this file is compiled in debug mode (_DEBUG_ > 0)
    when the gromacs version are not available.
    - string trimming function - duplicated from ~/src/gmxlib/string2.c 
    - case agnostic straing compare
 */
static void ltrim (char *str)
{
  char *tr;
  int c;

  if (!str)
    return;

  tr = strdup (str);
  c  = 0;
  while ((tr[c] == ' ') || (tr[c] == '\t'))
    c++;

  strcpy (str,tr+c);
  free (tr);
}

static void rtrim (char *str)
{
  int nul;

  if (!str)
    return;

  nul = strlen(str)-1;
  while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
    str[nul] = '\0';
    nul--;
  }
}

static void trim (char *str)
{
  ltrim (str);
  rtrim (str);
}

static int gmx_strncasecmp(const char* s1, const char* s2, int len)
{
  return strncasecmp(s1, s2, len);
}
#endif

/*! \endcond TEST */
