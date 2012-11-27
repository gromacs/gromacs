/*
 * This file is part of the GROMACS molecular simulation package,
 * version 4.6
 *
 * Author: Imran Haque, 2009
 * Copyright 2009, Stanford University
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
/*
 * memtestG80_core.h
 * Public API for core memory test functions for MemtestG80
 * Includes functional and OO interfaces to GPU test functions.
 */
#ifndef _MEMTESTG80_CORE_H_
#define _MEMTESTG80_CORE_H_

#if defined (WINDOWS) || defined (WINNV)
    #include <windows.h>
    inline unsigned int getTimeMilliseconds(void) {
        return GetTickCount();
    }
    #include <windows.h>
	#define SLEEPMS(x) Sleep(x)
#elif defined (LINUX) || defined (OSX)
    #include <sys/time.h>
    inline unsigned int getTimeMilliseconds(void) {
        struct timeval tv;
        gettimeofday(&tv,NULL);
        return tv.tv_sec*1000 + tv.tv_usec/1000;
    }
    #include <unistd.h>
    #define SLEEPMS(x) usleep(x*1000)
#else
    #error Must #define LINUX, WINDOWS, WINNV, or OSX
#endif

// By default the driver will spinwait when blocked on a kernel call
// Use the SOFTWAIT macro to replace this with a thread sleep and occasional poll
// limit expresses the max time we're willing to stay in the sleep loop - default = 15sec
inline int _pollStatus(unsigned length=1,unsigned limit=15000) {
    //while (cudaStreamQuery(0) != cudaSuccess) {SLEEPMS(length);}
    unsigned startTime = getTimeMilliseconds();
    while (cudaStreamQuery(0) == cudaErrorNotReady) {
        if ((getTimeMilliseconds() - startTime) > limit) return -1;
        SLEEPMS(length);
    }
    return 0;
}
#define SOFTWAIT() if (_pollStatus() != 0) {return 0xFFFFFFFE;} // -2
#define SOFTWAIT_LIM(lim) if (_pollStatus(1,lim) != 0) {return 0xFFFFFFFE;} // -2
//#define SOFTWAIT()
//#define SOFTWAIT(delay) if (_pollStatus(delay) != 0) return -2;
//#define SOFTWAIT(delay,limit) if (_pollStatus(delay,limit) != 0) return -2;
//#define SOFTWAIT() while (cudaStreamQuery(0) != cudaSuccess) {SLEEPMS(1);}
//#define SOFTWAIT(x) while (cudaStreamQuery(0) != cudaSuccess) {SLEEPMS(x);}
//#define SOFTWAIT()

// Use this macro to check for kernel errors
#define CHECK_LAUNCH_ERROR() if (cudaGetLastError() != cudaSuccess) {return 0xFFFFFFFF; /* -1 */}


typedef unsigned int uint;

// OO interface to MemtestG80 functions
class memtestState {
protected:
	const uint nBlocks;
	const uint nThreads;
    uint loopIters;
	uint megsToTest;
    int lcgPeriod;
	uint* devTestMem;
	uint* devTempMem;
	uint* hostTempMem;
	bool allocated;
public:
    uint initTime;
    memtestState() : nBlocks(1024), nThreads(512), loopIters(0), megsToTest(0), allocated(false), devTestMem(NULL),devTempMem(NULL),hostTempMem(NULL), initTime(0),lcgPeriod(1024) {};
    ~memtestState() {deallocate();}

	uint allocate(uint mbToTest);
	void deallocate();
	bool isAllocated() const {return allocated;}
	uint size() const {return megsToTest;}
    void setLCGPeriod(int period) {lcgPeriod = period;}
    int getLCGPeriod() const {return lcgPeriod;}

    bool gpuMemoryBandwidth(double& bandwidth,uint mbToTest,uint iters=5);
	bool gpuWriteConstant(const uint constant) const;
	bool gpuVerifyConstant(uint& errorCount,const uint constant) const;
	bool gpuShortLCG0(uint& errorCount,const uint repeats) const;
	bool gpuShortLCG0Shmem(uint& errorCount,const uint repeats) const;
	bool gpuMovingInversionsOnesZeros(uint& errorCount) const;
	bool gpuWalking8BitM86(uint& errorCount,const uint shift) const;
	bool gpuWalking8Bit(uint& errorCount,const bool ones,const uint shift) const;
	bool gpuMovingInversionsRandom(uint& errorCount) const;
	bool gpuWalking32Bit(uint& errorCount,const bool ones,const uint shift) const;
	bool gpuRandomBlocks(uint& errorCount,const uint seed) const;
	bool gpuModuloX(uint& errorCount,const uint shift,const uint pattern,const uint modulus,const uint overwriteIters) const;
};

// Utility functions
__host__ double gpuMemoryBandwidth(uint* src,uint* dst,uint mbToTest,uint iters);
__host__ void gpuWriteConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant);
__host__ uint gpuVerifyConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant,uint* blockErrorCount,uint* errorCounts);

__host__ void cpuWriteConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant);
__host__ uint cpuVerifyConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant);

// Logic tests
__host__ uint gpuShortLCG0(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint repeats,const int period,uint* blockErrorCounts,uint* errorCounts);
__host__ uint gpuShortLCG0Shmem(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint repeats,const int period,uint* blockErrorCounts,uint* errorCounts);

// Memtest86 Test 2: tseq=0,4
__host__ uint gpuMovingInversionsOnesZeros(const uint nBlocks,const uint nThreads,uint* base,uint N,uint* blockErrorCounts,uint* errorCounts);

// Memtest86 Test 3: tseq=1
__host__ uint gpuWalking8BitM86(const uint nBlocks,const uint nThreads,uint* base,uint N,uint shift,uint* blockErrorCounts,uint* errorCounts);
__host__ uint cpuWalking8BitM86(const uint nBlocks,const uint nThreads,uint* base,uint N,uint shift);
__host__ uint gpuWalking8Bit(const uint nBlocks,const uint nThreads,uint* base,uint N,bool ones,uint shift,uint* blockErrorCount,uint* errorCounts);

// Memtest86 Test 4: tseq=10
__host__ uint gpuMovingInversionsRandom(const uint nBlocks,const uint nThreads,uint* base,uint N,uint* blockErrorCounts,uint* errorCounts);

// Memtest86 Test 6: tseq=2
__host__ uint gpuWalking32Bit(const uint nBlocks,const uint nThreads,uint* base,uint N,bool ones,uint shift,uint* blockErrorCount,uint* errorCounts);
//
// Memtest86 Test 7: tseq=9
__host__ uint gpuRandomBlocks(const uint nBlocks,const uint nThreads,uint* base,uint N,uint seed,uint* blockErrorCount,uint* errorCounts);

// Memtest86 Test 8: tseq=3 (M86 uses modulus = 20)
__host__ uint gpuModuloX(const uint nBlocks,const uint nThreads,uint* base,const uint N,uint shift,uint pattern1,const uint modulus,const uint iters,uint* blockErrorCount,uint* errorCounts);

#endif
