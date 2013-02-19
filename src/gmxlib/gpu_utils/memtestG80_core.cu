/*
 * memtestG80_core.cu
 * MemtestG80 core memory test functions and OOP interface to tester.
 *
 * Author: Imran Haque, 2009
 * Copyright 2009, Stanford University
 *
 * This file is licensed under the terms of the LGPL. Please see
 * the COPYING file in the accompanying source distribution for
 * full license terms.
 *
 */

 /*
  * CUDA grid layout: Linear in blocks and threads.
  * Intended usage = 1k blocks, 512 t/blk, with N words (iterations) per thread
  *     -> 2*N MiB tested per grid
  * thread address at iteration i = base + blockIdx.x * N * blockDim.x + i*blockDim.x + threadIdx.x
  *
  */

// Naming convention: gpuXXX and cpuXXX functions are user-accessible; deviceXXX functions are internal
//                    gpuXXX functions execute a particular test on a block of GPU memory
//                    cpuXXX "          "      "   "         "    " "  "    "  CPU "

#define THREAD_ADDRESS(base,N,i) (base + blockIdx.x * N * blockDim.x + i * blockDim.x + threadIdx.x)
#define THREAD_OFFSET(N,i) (blockIdx.x * N * blockDim.x + i * blockDim.x + threadIdx.x)
#define BITSDIFF(x,y) __popc((x) ^ (y))


#include "memtestG80_core.h"

#include <stdio.h>




void memtestState::deallocate() {
		if (allocated) {
			cudaFree(devTestMem);
			cudaFree(devTempMem);
			free(hostTempMem);
			devTestMem = NULL;
			devTempMem = NULL;
			hostTempMem = NULL;
			allocated = false;
		}
        initTime = 0;
	}

uint memtestState::allocate(uint mbToTest) {
		deallocate();

        initTime = getTimeMilliseconds();
		
        // Round up to nearest 2MiB
		if (mbToTest % 2) mbToTest++;

		megsToTest = mbToTest;
		loopIters = megsToTest/2;

		if (megsToTest == 0) return 0;
		
		try {
			if (cudaMalloc((void**)&devTestMem,megsToTest*1048576UL) != cudaSuccess) throw 1;
			if (cudaMalloc((void**)&devTempMem,sizeof(uint)*nBlocks) != cudaSuccess) throw 2;
			if ( (hostTempMem = (uint*)malloc(sizeof(uint)*nBlocks)) == NULL) throw 3;
		} catch (...) {
            // Clear CUDA error flag for outside world
            cudaGetLastError();
			if (devTempMem) {
				cudaFree(devTempMem);
				devTempMem = NULL;
			}
			if (devTestMem) {
				cudaFree(devTestMem);
				devTestMem = NULL;
			}
			if (hostTempMem) {
				free(hostTempMem);
				hostTempMem = NULL;
			}
			return 0;
		}
		allocated = true;
		return megsToTest;
	}
bool memtestState::gpuMemoryBandwidth(double& bandwidth,uint mbToTest,uint iters) {
    if (!allocated || megsToTest < 2*mbToTest) return false;
    bandwidth = ::gpuMemoryBandwidth(devTestMem,devTestMem+mbToTest*1048576/4,mbToTest,iters);
    return cudaGetLastError() == cudaSuccess;
}
bool memtestState::gpuWriteConstant(const uint constant) const {
	if (!allocated) return false;
	::gpuWriteConstant(nBlocks,nThreads,devTestMem,loopIters,constant);
	return cudaGetLastError() == cudaSuccess;
}

bool memtestState::gpuVerifyConstant(uint& errorCount,const uint constant) const {
	if (!allocated) return false;
	errorCount = ::gpuVerifyConstant(nBlocks,nThreads,devTestMem,loopIters,constant,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}

bool memtestState::gpuShortLCG0(uint& errorCount,const uint repeats) const {
	if (!allocated) return false;
	errorCount = ::gpuShortLCG0(nBlocks,nThreads,devTestMem,loopIters,repeats,lcgPeriod,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuShortLCG0Shmem(uint& errorCount,const uint repeats) const {
	if (!allocated) return false;
	errorCount = ::gpuShortLCG0Shmem(nBlocks,nThreads,devTestMem,loopIters,repeats,lcgPeriod,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuMovingInversionsOnesZeros(uint& errorCount) const {
	if (!allocated) return false;
	errorCount = ::gpuMovingInversionsOnesZeros(nBlocks,nThreads,devTestMem,loopIters,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuWalking8BitM86(uint& errorCount,const uint shift) const {
	if (!allocated) return false;
	errorCount = ::gpuWalking8BitM86(nBlocks,nThreads,devTestMem,loopIters,shift,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuWalking8Bit(uint& errorCount,const bool ones,const uint shift) const {
	if (!allocated) return false;
	errorCount = ::gpuWalking8Bit(nBlocks,nThreads,devTestMem,loopIters,ones,shift,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuMovingInversionsRandom(uint& errorCount) const {
	if (!allocated) return false;
	errorCount = ::gpuMovingInversionsRandom(nBlocks,nThreads,devTestMem,loopIters,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuWalking32Bit(uint& errorCount,const bool ones,const uint shift) const {
	if (!allocated) return false;
	errorCount = ::gpuWalking32Bit(nBlocks,nThreads,devTestMem,loopIters,ones,shift,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuRandomBlocks(uint& errorCount,const uint seed) const {
	if (!allocated) return false;
	errorCount = ::gpuRandomBlocks(nBlocks,nThreads,devTestMem,loopIters,seed,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
bool memtestState::gpuModuloX(uint& errorCount,const uint shift,const uint pattern,const uint modulus,const uint overwriteIters) const {
	if (!allocated) return false;
	errorCount = ::gpuModuloX(nBlocks,nThreads,devTestMem,loopIters,shift,pattern,modulus,overwriteIters,devTempMem,hostTempMem);
	return ((cudaGetLastError() == cudaSuccess) && (errorCount != 0xFFFFFFFF) && (errorCount != 0xFFFFFFFE));
}
	
		

__global__ void deviceWriteConstant(uint* base, uint N, const uint constant);
__global__ void deviceVerifyConstant(uint* base,uint N,const uint constant,uint* blockErrorCount);
__global__ void deviceShortLCG0(uint* base,uint N,uint repeats,const int period);
__global__ void deviceShortLCG0Shmem(uint* base,uint N,uint repeats,const int period);
__global__ void deviceWriteRandomBlocks(uint* base,uint N,int seed);
__global__ void deviceVerifyRandomBlocks(uint* base,uint N,int seed,uint* blockErrorCount);
__global__ void deviceWriteWalking32Bit(uint* base,uint N,bool ones,uint shift);
__global__ void deviceVerifyWalking32Bit(uint* base,uint N,bool ones,uint shift,uint* blockErrorCount);
__global__ void deviceWritePairedConstants(uint* base,uint N,uint pattern0,uint pattern1);
__global__ void deviceVerifyPairedConstants(uint* base,uint N,uint pattern0,uint pattern1,uint* blockErrorCount);
__global__ void deviceWritePairedModulo(uint* base,const uint N,const uint shift,const uint pattern1,const uint pattern2,const uint modulus,const uint iters);
__global__ void deviceVerifyPairedModulo(uint* base,uint N,const uint shift,const uint pattern1,const uint modulus,uint* blockErrorCount);


// Utility function to measure memory bandwidth
__host__ double gpuMemoryBandwidth(uint* src,uint* dst,uint mbToTest,uint iters) {
       uint start = getTimeMilliseconds();
	   for (uint i = 0; i < iters; i++) {
           cudaMemcpy(dst,src,mbToTest*1048576,cudaMemcpyDeviceToDevice);
       }
       //D-to-D memory copies are non-blocking, so sync to get correct timing
       cudaThreadSynchronize();
       //SOFTWAIT();
       uint end = getTimeMilliseconds();
	   
       // Calculate bandwidth in MiB/s
	   // Multiply by 2 since we are reading and writing to the same memory
       double bw = 2.0*((double)mbToTest*iters)/((end-start)/1000.0);
	   return bw;
}

// Utility functions to write/verify pure constants in memory, CPU/GPU {{{
__host__ void gpuWriteConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant) { //{{{
    deviceWriteConstant<<<nBlocks,nThreads>>>(base,N,constant);
}

__global__ void deviceWriteConstant(uint* base, uint N, const uint constant) {
    for (uint i = 0 ; i < N; i++) {      
        *(THREAD_ADDRESS(base,N,i)) = constant;
    }
}
//}}}
__host__ uint gpuVerifyConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant,uint* blockErrorCount,uint* errorCounts) { //{{{
    // Given device arrays base (tested memory) and blockErrorCount (nBlocks uints in length of temp space)
    
	deviceVerifyConstant<<<nBlocks,nThreads,sizeof(uint)*nThreads>>>(base,N,constant,blockErrorCount);
	CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
	
    cudaMemcpy(errorCounts,blockErrorCount,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);

    // Sum-reduce block error counts on the host - it's only order of 1k numbers.
    uint totalErrors = 0;
    for (uint i = 0; i < nBlocks; i++) {
        totalErrors += errorCounts[i];
    }
    return totalErrors;
}

__global__ void deviceVerifyConstant(uint* base,uint N,const uint constant,uint* blockErrorCount) {
    // Verifies memory at base to make sure it has a constant pattern
    // Sums number of errors found in block and stores error count into blockErrorCount[blockIdx.x]
    // Sum-reduce this array afterwards to get total error count over tested region
    // Uses 4*blockDim.x bytes of shared memory
    
    extern __shared__ uint threadErrorCount[];
    threadErrorCount[threadIdx.x] = 0;

    for (uint i = 0; i < N; i++) {
        //if ( *(THREAD_ADDRESS(base,N,i)) != constant ) threadErrorCount[threadIdx.x]++;
        threadErrorCount[threadIdx.x] += BITSDIFF(*(THREAD_ADDRESS(base,N,i)),constant);
    }
    // Parallel-reduce error counts over threads in block
    for (uint stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
        if (threadIdx.x < stride)
            threadErrorCount[threadIdx.x] += threadErrorCount[threadIdx.x + stride];
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        blockErrorCount[blockIdx.x] = threadErrorCount[0];
    
    return;
}
//}}}

 __host__ void cpuWriteConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant) { //{{{
    dim3 blockDim(nThreads,0,0);
    dim3 threadIdx(0,0,0);
    dim3 blockIdx(0,0,0);
    for (blockIdx.x = 0; blockIdx.x < nBlocks; blockIdx.x++) {
        for (uint i = 0; i < N; i++) {
            for (threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
                *(THREAD_ADDRESS(base,N,i)) = constant;
            }
        }
    }
}
//}}}
__host__ uint cpuVerifyConstant(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint constant) { //{{{
    dim3 blockDim(nThreads,0,0);
    dim3 threadIdx(0,0,0);
    dim3 blockIdx(0,0,0);
    uint errorCount = 0;
    for (blockIdx.x = 0; blockIdx.x < nBlocks; blockIdx.x++) {
        for (uint i = 0; i < N; i++) {
            for (threadIdx.x = 0; threadIdx.x < blockDim.x; threadIdx.x++) {
                if (*(THREAD_ADDRESS(base,N,i)) != constant) errorCount++;
            }
        }
    }
    return errorCount;
} 
//}}}
//}}}

// Logic test 
// Idea: Run a varying number of iterations (k*N) of a short-period (per=N) LCG that returns to zero (or F's) quickly {{{
// Store only the result of the last iteration
// Compare output to the desired constant
// Compare results between varying k - memory error rate for a given pattern should be constant,
//                                     so variation should be due to logic errors in loop count
__host__ uint gpuShortLCG0(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint repeats,const int period,uint* blockErrorCounts,uint* errorCounts) { //{{{
    deviceShortLCG0<<<nBlocks,nThreads>>>(base,N,repeats,period);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
    CHECK_LAUNCH_ERROR();
    return gpuVerifyConstant(nBlocks,nThreads,base,N,0,blockErrorCounts,errorCounts);
} //}}}

__host__ uint gpuShortLCG0Shmem(const uint nBlocks,const uint nThreads,uint* base,uint N,const uint repeats,const int period,uint* blockErrorCounts,uint* errorCounts) { //{{{
    deviceShortLCG0Shmem<<<nBlocks,nThreads,sizeof(uint)*nThreads>>>(base,N,repeats,period);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
    CHECK_LAUNCH_ERROR();
    return gpuVerifyConstant(nBlocks,nThreads,base,N,0,blockErrorCounts,errorCounts);
} //}}}

// Put the LCG loop into a macro so we don't repeat code between versions of logic tester.
// The paired XOR adds diversity to the instruction stream, and is not reduced to a NOT
// as a single XOR is (verified with decuda).
// {{{
#if defined (LINUX) || defined(OSX)
#define LCGLOOP(var,repeats,period,a,c) for (uint rep = 0; rep < repeats; rep++) {\
    (var) = ~(var);\
    _Pragma("unroll 1")\
    for (uint iter = 0; iter < period; iter++) {\
        (var) = ~(var);\
        (var) = (a)*(var)+(c);\
        (var) ^= 0xFFFFFFF0;\
        (var) ^= 0xF;\
    }\
    (var) = ~(var);\
}
#elif defined (WINDOWS) || defined (WINNV)
#define LCGLOOP(var,repeats,period,a,c) for (uint rep = 0; rep < repeats; rep++) {\
    (var) = ~(var);\
    __pragma("unroll 1")\
    for (uint iter = 0; iter < period; iter++) {\
        (var) = ~(var);\
        (var) = (a)*(var)+(c);\
        (var) ^= 0xFFFFFFF0;\
        (var) ^= 0xF;\
    }\
    (var) = ~(var);\
}
#endif
//}}}

__global__ void deviceShortLCG0(uint* base,uint N,uint repeats,const int period) { //{{{
    // Pick a different block for different LCG lengths
    // Short periods are useful if LCG goes inside for i in 0..N loop
    int a,c;
    switch (period) {
        case 1024: a = 0x0fbfffff; c = 0x3bf75696; break;
        case 512:  a = 0x61c8647f; c = 0x2b3e0000; break;
        case 256:  a = 0x7161ac7f; c = 0x43840000; break;
        case 128:  a = 0x0432b47f; c = 0x1ce80000; break;
        case 2048: a = 0x763fffff; c = 0x4769466f; break;
        default:   a = 0; c = 0; break;
    }
    
    uint value = 0;
    LCGLOOP(value,repeats,period,a,c)

    for (uint i = 0 ; i < N; i++) {
        *(THREAD_ADDRESS(base,N,i)) = value;
    }
} //}}} 
// _shmem version uses shared memory to store inter-iteration values
// is more sensitive to shared memory errors from (eg) shader overclocking 
__global__ void deviceShortLCG0Shmem(uint* base,uint N,uint repeats,const int period) { //{{{
    // Pick a different block for different LCG lengths
    // Short periods are useful if LCG goes inside for i in 0..N loop
    int a,c;
    extern __shared__ uint shmem[];
    switch (period) {
        case 1024: a = 0x0fbfffff; c = 0x3bf75696; break;
        case 512:  a = 0x61c8647f; c = 0x2b3e0000; break;
        case 256:  a = 0x7161ac7f; c = 0x43840000; break;
        case 128:  a = 0x0432b47f; c = 0x1ce80000; break;
        case 2048: a = 0x763fffff; c = 0x4769466f; break;
        default:   a = 0; c = 0; break;
    }
    shmem[threadIdx.x] = 0;
    LCGLOOP(shmem[threadIdx.x],repeats,period,a,c)

    for (uint i = 0 ; i < N; i++) {
        *(THREAD_ADDRESS(base,N,i)) = shmem[threadIdx.x];

    }
} //}}} //}}}


// Memtest86 Test 2: tseq=0,4
__host__ uint gpuMovingInversionsOnesZeros(const uint nBlocks,const uint nThreads,uint* base,uint N,uint* blockErrorCounts,uint* errorCounts) { //{{{
    
    uint errorCount;
    gpuWriteConstant(nBlocks,nThreads,base,N,0xFFFFFFFF);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	errorCount = gpuVerifyConstant(nBlocks,nThreads,base,N,0xFFFFFFFF,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();

	gpuWriteConstant(nBlocks,nThreads,base,N,0x0);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	errorCount += gpuVerifyConstant(nBlocks,nThreads,base,N,0x0,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();
    return errorCount;
} //}}}

// Memtest86 Test 3: tseq=1
__host__ uint gpuWalking8BitM86(const uint nBlocks,const uint nThreads,uint* base,uint N,uint shift,uint* blockErrorCounts,uint* errorCounts) { //{{{
    // Performs the Memtest86 variation on the walking 8-bit pattern, where the same shifted pattern is
    // written into each 32-bit word in memory, verified, and its complement written and verified
    shift &= 0x7;
    uint pattern = 1 << shift;
    pattern = pattern | (pattern << 8) | (pattern << 16) | (pattern << 24);

    uint errorCount;
    gpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	errorCount = gpuVerifyConstant(nBlocks,nThreads,base,N,pattern,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();

	pattern = ~pattern;
    gpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	errorCount += gpuVerifyConstant(nBlocks,nThreads,base,N,pattern,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();
    return errorCount;
} //}}}
__host__ uint cpuWalking8BitM86(const uint nBlocks,const uint nThreads,uint* base,uint N,uint shift) { //{{{
    // Performs the Memtest86 variation on the walking 8-bit pattern, where the same shifted pattern is
    // written into each 32-bit word in memory, verified, and its complement written and verified
    shift &= 0x7;
    uint pattern = 1 << shift;
    pattern = pattern | (pattern << 8) | (pattern << 16) | (pattern << 24);

    uint errorCount;
    cpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    errorCount = cpuVerifyConstant(nBlocks,nThreads,base,N,pattern);

    pattern = ~pattern;
    cpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    errorCount += cpuVerifyConstant(nBlocks,nThreads,base,N,pattern);

    return errorCount;
} //}}}
__host__ uint gpuWalking8Bit(const uint nBlocks,const uint nThreads,uint* base,uint N,bool ones,uint shift,uint* blockErrorCount,uint* errorCounts) { //{{{
    // Implements one iteration of true walking 8-bit ones/zeros test
    uint patterns[2]={0x0,0x0};
    
    // Build the walking-ones paired pattern of 8-bits with the given shift
    shift &= 0x7;
    uint bits = 0x1 << shift;
    for (uint i = 0; i < 4; i++) {
        patterns[0] = (patterns[0] << 8) | bits;
        bits = (bits == 0x80) ? 0x01 : bits<<1;
    }
    for (uint i = 0; i < 4; i++) {
        patterns[1] = (patterns[1] << 8) | bits;
        bits = (bits == 0x80) ? 0x01 : bits<<1;
    }

    if (!ones) {
        patterns[0] = ~patterns[0];
        patterns[1] = ~patterns[1];
    }
	
	//printf("Host Patterns: %08x %08x\n",patterns[0],patterns[1]);
    deviceWritePairedConstants<<<nBlocks,nThreads>>>(base,N,patterns[0],patterns[1]);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
	//cudaMemcpy(errorCounts,base,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);
    //printf("First few words in tested RAM: %08x %08x %08x %08x %08x %08x\n",errorCounts[0],errorCounts[1],errorCounts[2],errorCounts[3],errorCounts[4],errorCounts[5]);
    // Given device arrays base (tested memory) and blockErrorCount (nBlocks uints in length of temp space)
    deviceVerifyPairedConstants<<<nBlocks,nThreads,sizeof(uint)*nThreads>>>(base,N,patterns[0],patterns[1],blockErrorCount);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
    //if (cudaGetLastError() != cudaSuccess) {
	//	return 0xFFFFFFFF; // -1
	//}
	//uint errorCounts[nBlocks];
    cudaMemcpy(errorCounts,blockErrorCount,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);

    // Sum-reduce block error counts on the host - it's only order of 1k numbers.
    uint totalErrors = 0;
    for (uint i = 0; i < nBlocks; i++) {
        totalErrors += errorCounts[i];
    }
    return totalErrors;
}

__global__ void deviceWritePairedConstants(uint* base,uint N,uint pattern0,uint pattern1) {
    // Writes paired constants to memory, such that each offset that is X mod 2 receives patterns[X]
    // Used for true walking-ones/zeros 8-bit test
    //if (threadIdx.x == 0)
    //    printf("Device Patterns Block %u: %08x %08x\n",blockIdx.x,patterns[0],patterns[1]);
    const uint pattern = (threadIdx.x & 0x1) ? pattern1 : pattern0;
    //const uint pattern = patterns[threadIdx.x & 0x1];
    for (uint i = 0 ; i < N; i++) {      
        *(THREAD_ADDRESS(base,N,i)) = pattern;
        //*(base+blockIdx.x*N*blockDim.x + i*blockDim.x + threadIdx.x) = 0;
    }

}

__global__ void deviceVerifyPairedConstants(uint* base,uint N,uint pattern0,uint pattern1,uint* blockErrorCount) {
    // Verifies memory at base to make sure it has a correct paired-constant pattern
    // Sums number of errors found in block and stores error count into blockErrorCount[blockIdx.x]
    // Sum-reduce this array afterwards to get total error count over tested region
    // Uses 4*blockDim.x bytes of shared memory
    
    extern __shared__ uint threadErrorCount[];
    threadErrorCount[threadIdx.x] = 0;
    //const uint pattern = patterns[threadIdx.x & 0x1];
    const uint pattern = (threadIdx.x & 0x1) ? pattern1 : pattern0;
    
    for (uint i = 0; i < N; i++) {
        //if ( *(THREAD_ADDRESS(base,N,i)) != pattern ) threadErrorCount[threadIdx.x]++;
        threadErrorCount[threadIdx.x] += BITSDIFF(*(THREAD_ADDRESS(base,N,i)),pattern);
    }
    // Parallel-reduce error counts over threads in block
    for (uint stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
        if (threadIdx.x < stride)
            threadErrorCount[threadIdx.x] += threadErrorCount[threadIdx.x + stride];
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        blockErrorCount[blockIdx.x] = threadErrorCount[0];
    
    return;
}
//}}}

// Memtest86 Test 4: tseq=10
__host__ uint gpuMovingInversionsRandom(const uint nBlocks,const uint nThreads,uint* base,uint N,uint* blockErrorCounts,uint* errorCounts) { //{{{
    
    uint errorCount;

    uint pattern = (uint)rand();
    gpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
	
	errorCount = gpuVerifyConstant(nBlocks,nThreads,base,N,pattern,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();
    
	pattern = ~pattern;
    gpuWriteConstant(nBlocks,nThreads,base,N,pattern);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
	
	errorCount += gpuVerifyConstant(nBlocks,nThreads,base,N,pattern,blockErrorCounts,errorCounts);
	CHECK_LAUNCH_ERROR();
    return errorCount;
} //}}}

// Memtest86 Test 6: tseq=2
__host__ uint gpuWalking32Bit(const uint nBlocks,const uint nThreads,uint* base,uint N,bool ones,uint shift,uint* blockErrorCount,uint* errorCounts) { //{{{
    // Given device arrays base (tested memory) and blockErrorCount (nBlocks uints in length of temp space)
    // Does one iteration of the walking-{ones/zeros} 32-bit test paralleling Memtest
    // With the starting pattern 1<<shift
    // NUMBER OF THREADS SHOULD BE A MULTIPLE OF 32

    deviceWriteWalking32Bit<<<nBlocks,nThreads>>>(base,N,ones,shift);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	deviceVerifyWalking32Bit<<<nBlocks,nThreads,sizeof(uint)*nThreads>>>(base,N,ones,shift,blockErrorCount);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

    cudaMemcpy(errorCounts,blockErrorCount,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);

    // Sum-reduce block error counts on the host - it's only order of 1k numbers.
    uint totalErrors = 0;
    for (uint i = 0; i < nBlocks; i++) {
        totalErrors += errorCounts[i];
    }
    return totalErrors;
    
}

__global__ void deviceWriteWalking32Bit(uint* base,uint N,bool ones,uint shift) {
    // Writes one iteration of the walking-{ones/zeros} 32-bit pattern to gpu memory

    // Want to write in a 1 << (offset from base + shift % 32)
    // Since thread indices are aligned with base, this reduces to
    // 1 << ((threadIdx.x+shift) & 0x1f)
    // With conditional inversion for walking zeros
    uint pattern = 1 << ((threadIdx.x + shift) & 0x1f);
    pattern = ones ? pattern : ~pattern;
    
    for (uint i = 0; i < N; i++) {
        *(THREAD_ADDRESS(base,N,i)) = pattern;
    }
}

__global__ void deviceVerifyWalking32Bit(uint* base,uint N,bool ones,uint shift,uint* blockErrorCount) {
    // Verifies memory at base to make sure it has a constant pattern
    // Sums number of errors found in block and stores error count into blockErrorCount[blockIdx.x]
    // Sum-reduce this array afterwards to get total error count over tested region
    // Uses 4*blockDim.x bytes of shared memory
    
    extern __shared__ uint threadErrorCount[];
    threadErrorCount[threadIdx.x] = 0;

    uint pattern = 1 << ((threadIdx.x + shift) & 0x1f);
    pattern = ones ? pattern : ~pattern;
    
    for (uint i = 0; i < N; i++) {
        //if ( *(THREAD_ADDRESS(base,N,i)) != pattern ) threadErrorCount[threadIdx.x]++;
        threadErrorCount[threadIdx.x] += BITSDIFF(*(THREAD_ADDRESS(base,N,i)),pattern);
    }
    // Parallel-reduce error counts over threads in block
    for (uint stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
        if (threadIdx.x < stride)
            threadErrorCount[threadIdx.x] += threadErrorCount[threadIdx.x + stride];
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        blockErrorCount[blockIdx.x] = threadErrorCount[0];
    
    return;
}
//}}}

// Memtest86 Test 7: tseq=9
__host__ uint gpuRandomBlocks(const uint nBlocks,const uint nThreads,uint* base,uint N,uint seed,uint* blockErrorCount,uint* errorCounts) { //{{{ {{{
    // Writes random numbers into memory and verifies pattern
    //uint errorCounts[nBlocks];
    
    deviceWriteRandomBlocks<<<nBlocks,nThreads,4*nThreads>>>(base,N,seed);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();

	//cudaMemcpy(errorCounts,base,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);
    //printf("First few words in tested RAM: %08x %08x %08x %08x %08x %08x\n",errorCounts[0],errorCounts[1],errorCounts[2],errorCounts[3],errorCounts[4],errorCounts[5]);
	
	deviceVerifyRandomBlocks<<<nBlocks,nThreads,12*nThreads>>>(base,N,seed,blockErrorCount);
    CHECK_LAUNCH_ERROR();
    SOFTWAIT();
	CHECK_LAUNCH_ERROR();
	
	
    cudaMemcpy(errorCounts,blockErrorCount,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);

    // Sum-reduce block error counts on the host - it's only order of 1k numbers.
    uint totalErrors = 0;
    for (uint i = 0; i < nBlocks; i++) {
        totalErrors += errorCounts[i];
    }
    return totalErrors;
}
//}}}
//
// Math functions modulo the Mersenne prime 2^31 -1 {{{
__device__ void deviceMul3131 (uint v1, uint v2,uint& LO, uint& HI)
{
    // Given v1, v2 < 2^31
    // Emulate a 31-bit integer multiply by doing instead a 32-bit multiply into LO and HI
    // And shifting bits around to make it look right.
    LO = v1*v2;
    HI = __umulhi(v1,v2);
    HI <<= 1;
    HI |= (LO & 0x80000000) >> 31;
    LO &= 0x7FFFFFFF;
    
}

__device__ uint deviceModMP31(uint LO,uint HI) {
    // Modulo a 62-bit number HI<<31 + LO, mod 2^31-1
    // Encyclopedia of Cryptography and Security By Henk C. A. van Tilborg
    // page 381, Mersenne Primes
    uint sum = LO+HI;
    if (sum >= 0x80000000) {
        // If a+b > 2^31, then high bit will be set
        return sum - 0x80000000 + 1;
    } else {
        return sum;
    }
}
__device__ uint deviceMulMP31(uint a,uint b) {
    // Multiplies a pair of 31-bit integers a and b mod the Mersenne prime 2^31-1
    // Takes result through a 62-bit intermediate
    uint LO,HI;
    deviceMul3131(a,b,LO,HI);
    return deviceModMP31(LO,HI);
}

__device__ uint deviceExpoModMP31(uint base,uint exponent) {
    uint result = 1;
    while (exponent > 0) {
        if (exponent & 1) {
            result = deviceMulMP31(result,base);
        }
        exponent >>= 1;
        base = deviceMulMP31(base,base);
    }
    return result;
}
//}}}
// deviceRan0p: Parallelized closed-form version of NR's ran0  {{{
__device__ uint deviceRan0p(int seed,int n) { // 
    uint an = deviceExpoModMP31(16807,n+1);
    return deviceMulMP31(an,seed);
}
//}}}
// deviceIrbit2: random bit generation, from NR {{{
__device__ int deviceIrbit2(uint& seed) {
    const uint IB1  = 1;
    const uint IB2  = 2;
    const uint IB5  = 16;
    const uint IB18 = 131072;
    const uint MASK = IB1+IB2+IB5;
    if (seed & IB18) {
        seed = ((seed ^ MASK) << 1) | IB1;
        return 1;
    } else {
        seed <<= 1;
        return 0;
    }
}
//}}}
__global__ void deviceWriteRandomBlocks(uint* base,uint N,int seed) { //{{{
    // Requires 4*nThreads bytes of shared memory
    extern __shared__ uint randomBlock[];

    // Make sure seed is not zero.
    if (seed == 0) seed = 123459876+blockIdx.x;
    uint bitSeed = deviceRan0p(seed + threadIdx.x,threadIdx.x);

    for (uint i=0; i < N; i++) {
        // Generate a block of random numbers in parallel using closed-form expression for ran0
        // OR in a random bit because Ran0 will never have the high bit set
        randomBlock[threadIdx.x] = deviceRan0p(seed,threadIdx.x) | (deviceIrbit2(bitSeed) << 31);
        __syncthreads();
        
        // Set the seed for the next round to the last number calculated in this round
        seed = randomBlock[blockDim.x-1];
        
        // Blit shmem block out to global memory
        *(THREAD_ADDRESS(base,N,i)) = randomBlock[threadIdx.x];
    }
}
//}}}
__global__ void deviceVerifyRandomBlocks(uint* base,uint N,int seed,uint* blockErrorCount) { //{{{
    // Verifies memory at base to make sure it has a correct random pattern given the seed
    // Sums number of errors found in block and stores error count into blockErrorCount[blockIdx.x]
    // Sum-reduce this array afterwards to get total error count over tested region
    // Uses 12*blockDim.x bytes of shared memory
    
    extern __shared__ uint shmem[];
    uint* threadErrorCount = shmem;
    uint* randomBlock = shmem + blockDim.x;
    // Put these into shmem to cut register count
    uint* bitSeeds = randomBlock + blockDim.x;
    
    threadErrorCount[threadIdx.x] = 0;

    // Make sure seed is not zero.
    if (seed == 0) seed = 123459876+blockIdx.x;
    //uint bitSeed = deviceRan0p(seed + threadIdx.x,threadIdx.x);
    bitSeeds[threadIdx.x] = deviceRan0p(seed + threadIdx.x,threadIdx.x);
    
    for (uint i = 0; i < N; i++) {
        // Generate a block of random numbers in parallel using closed-form expression for ran0
        // OR in a random bit because Ran0 will never have the high bit set
        //randomBlock[threadIdx.x] = deviceRan0p(seed,threadIdx.x) | (deviceIrbit2(bitSeed) << 31);
        randomBlock[threadIdx.x] = deviceRan0p(seed,threadIdx.x) | (deviceIrbit2(bitSeeds[threadIdx.x]) << 31);
        __syncthreads();
        
        // Set the seed for the next round to the last number calculated in this round
        seed = randomBlock[blockDim.x-1];
        
        //if ( randomBlock[threadIdx.x] != *(THREAD_ADDRESS(base,N,i))) threadErrorCount[threadIdx.x]++;
        threadErrorCount[threadIdx.x] += BITSDIFF(*(THREAD_ADDRESS(base,N,i)),randomBlock[threadIdx.x]);
        
    }

    // Parallel-reduce error counts over threads in block
    for (uint stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
        if (threadIdx.x < stride)
            threadErrorCount[threadIdx.x] += threadErrorCount[threadIdx.x + stride];
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        blockErrorCount[blockIdx.x] = threadErrorCount[0];
    
    return;
}
//}}}
//}}}

// Memtest86 Test 8: tseq=3 (M86 uses modulus = 20)
__host__ uint gpuModuloX(const uint nBlocks,const uint nThreads,uint* base,const uint N,uint shift,uint pattern1,const uint modulus,const uint iters,
						 uint* blockErrorCount,uint* errorCounts) { //{{{
    // Given device arrays base (tested memory) and blockErrorCount (nBlocks uints in length of temp space)
    // Given a shift, modulus, pattern to test and number of overwrite iterations
    // Performs Modulo-X test on memory
    
    //uint errorCounts[nBlocks];
    uint totalErrors = 0;
    shift %= modulus;

    // Test both the given pattern and its inverse
    for (uint i = 0; i < 2; i++, pattern1 = ~pattern1) {
        deviceWritePairedModulo<<<nBlocks,nThreads>>>(base,N,shift,pattern1,~pattern1,modulus,iters);
	    CHECK_LAUNCH_ERROR();
        SOFTWAIT();
	    CHECK_LAUNCH_ERROR();

		deviceVerifyPairedModulo<<<nBlocks,nThreads,sizeof(uint)*nThreads>>>(base,N,shift,pattern1,modulus,blockErrorCount);
		CHECK_LAUNCH_ERROR();
        SOFTWAIT();
	    CHECK_LAUNCH_ERROR();

        cudaMemcpy(errorCounts,blockErrorCount,sizeof(uint)*nBlocks,cudaMemcpyDeviceToHost);

        // Sum-reduce block error counts on the host - it's only order of 1k numbers.
        for (uint i = 0; i < nBlocks; i++) {
            totalErrors += errorCounts[i];
        }
    }
    return totalErrors;
}

__global__ void deviceWritePairedModulo(uint* base,const uint N,const uint shift,const uint pattern1,const uint pattern2,const uint modulus,const uint iters) {
    // First writes pattern1 into every offset that is 0 mod modulus
    // Next  (iters times) writes ~pattern1 into every other address
    uint offset;
    for (uint i = 0 ; i < N; i++) {      
        offset = THREAD_OFFSET(N,i);
        if ((offset % modulus) == shift) *(base+offset) = pattern1;
    }
    __syncthreads();
    for (uint j = 0; j < iters; j++) {
        for (uint i = 0 ; i < N; i++) {      
            offset = THREAD_OFFSET(N,i);
            if ((offset % modulus) != shift) *(base+offset) = pattern2;
        }
    }
}
__global__ void deviceVerifyPairedModulo(uint* base,uint N,const uint shift,const uint pattern1,const uint modulus,uint* blockErrorCount) {
    // Verifies that memory at each (offset mod modulus == shift) stores pattern1
    // Sums number of errors found in block and stores error count into blockErrorCount[blockIdx.x]
    // Sum-reduce this array afterwards to get total error count over tested region
    // Uses 4*blockDim.x bytes of shared memory
    
    extern __shared__ uint threadErrorCount[];
    threadErrorCount[threadIdx.x] = 0;
    uint offset;
    
    for (uint i = 0; i < N; i++) {
        offset = THREAD_OFFSET(N,i);
        //if (((offset % modulus) == shift) && (*(base+offset) != pattern1)) threadErrorCount[threadIdx.x]++;
        if ((offset % modulus) == shift) threadErrorCount[threadIdx.x] += BITSDIFF(*(base+offset),pattern1);
    }
    // Parallel-reduce error counts over threads in block
    for (uint stride = blockDim.x>>1; stride > 0; stride >>= 1) {
        __syncthreads();
        if (threadIdx.x < stride)
            threadErrorCount[threadIdx.x] += threadErrorCount[threadIdx.x + stride];
    }
    __syncthreads();
    
    if (threadIdx.x == 0)
        blockErrorCount[blockIdx.x] = threadErrorCount[0];
    
    return;
}
//}}}
