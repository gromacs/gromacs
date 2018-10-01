/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/


#pragma once
#if !defined( AMD_CLFFT_plan_H )
#define AMD_CLFFT_plan_H
#include <cstring>
#include "private.h"
#include "lock.h"
#include "generator.h"

std::string getKernelName(const clfftGenerators gen, const clfftPlanHandle plHandle, bool withPlHandle);

namespace ARBITRARY {
	// TODO:  These arbitrary parameters should be tuned for the type of GPU
	//	being used.  These values are probably OK for Radeon 58xx and 68xx.
	enum {
		MAX_DIMS  = 3,
			//  The clEnqueuNDRangeKernel accepts a multi-dimensional domain array.
			//  The # of dimensions is arbitrary, but limited by the OpenCL implementation
			//  usually to 3 dimensions (CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS).
			//  The kernel generator also assumes a limit on the # of dimensions.

		SIMD_WIDTH = 64,
			//  Workgroup size.  This is the # of work items that share
			//  local data storage (LDS).  This # is best for Evergreen gpus,
			//  but might change in the future.

		LDS_BANK_BITS = 5,
		LDS_BANK_SIZE = (1 << LDS_BANK_BITS),
		LDS_PADDING   = false,//true,
			//  On AMD hardware, the low-order bits of the local_id enumerate
			//  the work items that access LDS in parallel.  Ideally, we will
			//  pad our LDS arrays so that these work items access different banks
			//  of the LDS.
			//  2 ** LDS_BANK_BITS is the number of LDS banks.
			//  If LDS_PADDING is non-zero, the kernel generator should pad the
			//  LDS arrays to reduce or eliminate bank conflicts.

		LDS_FRACTION_IDEAL = 6,    // i.e., 1/6th
		LDS_FRACTION_MAX   = 4,    // i.e., 1/4
			//  For best performance, each workgroup should use 1/IDEAL'th the amount of LDS
			//  revealed by clGetDeviceInfo (.. CL_DEVICE_LOCAL_MEM_SIZE, ...)
			//  However, we can use up to 1/MAX'th of LDS per workgroup when necessary to
			//  perform the FFT in a single pass instead of multiple passes.
			//  This tuning parameter is a good value for Evergreen gpus,
			//  but might change in the future.

		LDS_COMPLEX = false,
			//  This is the default value for FFTKernelGenKeyParams::fft_LdsComplex.
			//  The generated kernels require so many bytes of LDS for each single precision
			//..complex number in the vector.
			//  If LDS_COMPLEX, then we declare an LDS array of complex numbers (8 bytes each)
			//  and swap data between workitems with a single barrier.
			//  If ! LDS_COMPLEX, then we declare an LDS array or scalar numbers (4 bytes each)
			//  and swap data between workitems in two phases, with extra barriers.
			//  The former approach uses fewer instructions and barriers;
			//  The latter uses half as much LDS space, so twice as many wavefronts can be run
			//  in parallel.

		TWIDDLE_DEE = 8,
			//  number of bits per row of matrix.
	};

};


enum BlockComputeType
{
	BCT_C2C,	// Column to column
	BCT_C2R,	// Column to row
	BCT_R2C,	// Row to column
};


//NonSquareKernelType
enum NonSquareTransposeKernelType
{
    NON_SQUARE_TRANS_PARENT,
    NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING,
    NON_SQUARE_TRANS_TRANSPOSE_BATCHED,
    NON_SQUARE_TRANS_SWAP
};

/*
There are three ways of conducting inplace transpose with 1:2 (or 2:1) dimension ratio.
A. first conduct line swapping kernels for the whole non square matrix
then conduct batched square transpose along column dim (a 'real' batched transpose)
B. first conduct batched square transpose along column dim (a 'real' batched transpose)
then conduct line swapping kernels for the whole non square matrix (for 2:1 case)
C. first conduct batched square transpose along leading dim (row dim)
then conduct line swapping kernels for the whole non square matrix
Note that the twiddle computation has to go at the begining of the first kernel or the end of the second kernel

if leading dimension is bigger, it makes more sense (faster) to swap line first and then conduct batched square transpose
if leading dimension is smaller, it makes more sense (faster) to conduct batched transpose and then swap lines.
*/
enum NON_SQUARE_KERNEL_ORDER
{
	NOT_A_TRANSPOSE,
	SWAP_AND_TRANSPOSE, // A.
	TRANSPOSE_AND_SWAP, // B.
	TRANSPOSE_LEADING_AND_SWAP, // C.
};

#define CLFFT_CB_SIZE 32
#define CLFFT_MAX_INTERNAL_DIM 16

/*! @brief Data structure to store the callback function string and other metadata passed by client 
*  @details Client sets the callback function and other required parameters through clfftSetPlanCallback() 
*  in order to register the callback function. The library populates these values into this data structure
*/ 
typedef struct clfftCallbackParam_
{
	int localMemSize;			/*!< optional local memory size if needed by callback */
	const char* funcname;		/*!< callback function name */
	const char* funcstring;		/*!< callback function in string form */
}clfftCallbackParam;

struct FFTKernelGenKeyParams {
	/*
	 *	This structure distills a subset of the fftPlan data,
	 *	including all information that is used to generate the OpenCL kernel.
	 *	This structure can be used as a key to reusing kernels that have already
	 *	been compiled.
	 */
	size_t                   fft_DataDim;       // Dimensionality of the data
	size_t                   fft_N[CLFFT_MAX_INTERNAL_DIM];          // [0] is FFT size, e.g. 1024
	                                            // This must be <= size of LDS!
	size_t                   fft_inStride [CLFFT_MAX_INTERNAL_DIM];  // input strides
	size_t                   fft_outStride[CLFFT_MAX_INTERNAL_DIM];  // output strides

	clfftResultLocation   fft_placeness;
	clfftLayout           fft_inputLayout;
	clfftLayout           fft_outputLayout;
	clfftPrecision        fft_precision;
	double                   fft_fwdScale;
	double                   fft_backScale;

	size_t                   fft_SIMD;          // Assume this SIMD/workgroup size
	size_t                   fft_LDSsize;       // Limit the use of LDS to this many bytes.
	size_t                   fft_R;             // # of complex values to keep in working registers
	                                            // SIMD size * R must be <= size of LDS!

	size_t					 fft_MaxWorkGroupSize; // Limit for work group size

	bool                     fft_3StepTwiddle;  // This is one pass of the "3-step" algorithm;
	                                            // so extra twiddles are applied on output.
	bool					 fft_twiddleFront;	// do twiddle scaling at the beginning pass

	bool					 fft_realSpecial;	// this is the flag to control the special case step (4th step)
	                                            // in the 5-step real 1D large breakdown
	size_t					 fft_realSpecial_Nr;

	bool                     fft_RCsimple; 

	bool					 transOutHorizontal;	// tiles traverse the output buffer in horizontal direction

	bool					 blockCompute;
	BlockComputeType		 blockComputeType;
	size_t					 blockSIMD;
	size_t					 blockLDS;
    
	NonSquareTransposeKernelType      nonSquareKernelType;
	// sometimes non square matrix are broken down into a number of
	// square matrix during inplace transpose
	// let's call this number transposeMiniBatchSize
	// no user of the library should set its value
	size_t transposeMiniBatchSize;
	// transposeBatchSize is the number of batchs times transposeMiniBatchSzie
	// no user of the library should set its value
	size_t transposeBatchSize;
	// no user of the library should set its value 
	NON_SQUARE_KERNEL_ORDER nonSquareKernelOrder;

	bool fft_hasPreCallback;
	clfftCallbackParam fft_preCallback;

	bool fft_hasPostCallback;
	clfftCallbackParam fft_postCallback;

	cl_ulong   limit_LocalMemSize;

	// Default constructor
	FFTKernelGenKeyParams()
	{
		fft_DataDim = 0;
		for(int i=0; i<CLFFT_MAX_INTERNAL_DIM; i++)
		{
			fft_N[i] = 0;
			fft_inStride[i] = 0;
			fft_outStride[i] = 0;
		}

		fft_placeness = CLFFT_OUTOFPLACE;
		fft_inputLayout = CLFFT_COMPLEX_INTERLEAVED;
		fft_outputLayout = CLFFT_COMPLEX_INTERLEAVED;
		fft_precision = CLFFT_SINGLE;
		fft_fwdScale = fft_backScale = 0.0;
		fft_SIMD = 0;
		fft_LDSsize = 0;
		fft_R = 0;
		fft_MaxWorkGroupSize = 0;
		fft_3StepTwiddle = false;
		fft_twiddleFront = false;

		transOutHorizontal = false;

		fft_realSpecial = false;
		fft_realSpecial_Nr = 0;

		fft_RCsimple = false;

		blockCompute = false;
		blockComputeType = BCT_C2C;
		blockSIMD = 0;
		blockLDS = 0;
        nonSquareKernelType = NON_SQUARE_TRANS_PARENT;
		transposeMiniBatchSize = 1;
		transposeBatchSize = 1;
		fft_hasPreCallback = false;
		fft_hasPostCallback = false;
		limit_LocalMemSize = 0;
	}
};


//	Sorting operator for struct FFTKernelGenKeyParams, such that it can be used in a map
bool operator<( const FFTKernelGenKeyParams& lhs, const FFTKernelGenKeyParams& rhs);

class	FFTPlan;
class   FFTRepo;

// Action ID
enum FFTActionImplID
{
    FFT_DEFAULT_STOCKHAM_ACTION,
    FFT_DEFAULT_TRANSPOSE_ACTION,
    FFT_DEFAULT_COPY_ACTION,
    FFT_STATIC_STOCKHAM_ACTION
};


// 
// FFTKernelSignatureHeader
// 
// This structure is a wrapper for the FFTKernelSignature.  
// It stores the signature size and the action ID. This ensure that every 
// FFTKernelSignature (even with an empty DATA) is unique
// 
// This class is used as the return type of FFTAction::getSignatureData()
// 
struct FFTKernelSignatureHeader
{
    int datasize;
    FFTActionImplID id;

    //clfftLayout           fft_inputLayout;
    //clfftLayout           fft_outputLayout;

    FFTKernelSignatureHeader(int size_, FFTActionImplID id_)
    {
        // Set to 0 the whole signature structure
        ::memset(this, 0, size_);
        datasize = size_;
        id = id_;
    }
};

// 
// FFTKernelSignature
// 
// This struct represents the signature of an action. An action signature 
// stores (by inheritage):
//  - the action ID
//  - its signature data size
//  - a set of bytes caracterizes a FFT action
// 
// This template class FFTKernelSignature provides a simple mechanism to
// build a proper signature (see also in src/library/repo.h) from any POD type.
//  
// It is used as a key in the different cache mecanisms (binary cache and
// dynamic cache managed by FFTRepo)
// 
template <typename DATA, FFTActionImplID ID>
struct FFTKernelSignature : public FFTKernelSignatureHeader, public DATA
{
    FFTKernelSignature()
        : FFTKernelSignatureHeader(sizeof(FFTKernelSignature<DATA, ID>), ID)
    {
    }
};



// 
// FFTAction is the base class for all actions used to implement FFT computes
// 
// An action basically implements some OpenCL related actions, for instance:
//  - Generating a kernel source code from kernel characteristics
//  - Computing the work-group local sizes according to a kernel
//  - Enqueuing arguments to the kernel call
// 
// Kernels generated and compiled by an action are stored in the different
// cache mechanism (see repo.h for the dynamic cache and fft_binary_lookup.h
// for disk cache for more information). Each cache entry can be obtained from
// the action signature which is set of byte characterizing the action itself.
// 
// At the moment, FFTAction only implements the enqueue method which is
// inherited by every action subclasses. But it may change in the time. There
// are no clear rules and the choices made until now are still subject to
// change.
// 
class FFTAction
{
public:
    FFTAction(FFTPlan * plan, clfftStatus & err);

    virtual clfftStatus enqueue(clfftPlanHandle plHandle,
                                clfftDirection dir,
                                cl_uint numQueuesAndEvents,
                                cl_command_queue* commQueues,
                                cl_uint numWaitEvents,
                                const cl_event* waitEvents,
                                cl_event* outEvents,
                                cl_mem* clInputBuffers,
                                cl_mem* clOutputBuffers);

protected:

    virtual clfftGenerators                getGenerator() = 0;

    clfftStatus                            compileKernels  ( const cl_command_queue commQueueFFT, const clfftPlanHandle plHandle, FFTPlan* fftPlan);
    clfftStatus                            writeKernel     ( const clfftPlanHandle plHandle, const clfftGenerators gen, const FFTKernelSignatureHeader* data, const cl_context& context, const cl_device_id &device);

    virtual clfftStatus                    generateKernel  ( FFTRepo & fftRepo, const cl_command_queue commQueueFFT) = 0;
    virtual clfftStatus                    getWorkSizes    ( std::vector<size_t> & globalws, std::vector<size_t> & localws) = 0;

    virtual const FFTKernelSignatureHeader * getSignatureData() = 0;

    FFTPlan * plan;

private:

    clfftStatus selectBufferArguments(FFTPlan * plan,
                                      cl_mem* clInputBuffers,
                                      cl_mem* clOutputBuffers,
                                      std::vector< cl_mem > &inputBuff,
                                      std::vector< cl_mem > &outputBuff);

    virtual bool buildForwardKernel() = 0;
    virtual bool buildBackwardKernel() = 0;
};


//	The "envelope" is a set of limits imposed by the hardware
//	This will depend on the GPU(s) in the OpenCL context.
//	If there are multiple devices, this should be the least
//	common denominators.
//
struct FFTEnvelope {
	cl_ulong   limit_LocalMemSize;
	           //  this is the minimum of CL_DEVICE_LOCAL_MEM_SIZE
	size_t     limit_Dimensions;
	           //  this is the minimum of CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
	size_t     limit_Size[8];
	           //  these are the minimima of CL_DEVICE_MAX_WORK_ITEM_SIZES[0..n]
	size_t     limit_WorkGroupSize;
	           //  this is the minimum of CL_DEVICE_MAX_WORK_GROUP_SIZE

	// ??  CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE

	FFTEnvelope ()
	:	limit_LocalMemSize (0)
	,	limit_Dimensions (0)
	,	limit_WorkGroupSize (0)
	{
		::memset( &limit_Size, 0, sizeof( limit_Size ) );
	}
};


//	This class contains objects that are specific to a particular FFT transform, and the data herein is useful
//	for us to know ahead of transform time such that we can optimize for these settings
class	FFTPlan
{

public:

	bool baked;

	//	Properties provided by the user.
	clfftDim             dim;
	clfftLayout          inputLayout;
	clfftLayout          outputLayout;
	clfftResultLocation  placeness;
	clfftResultTransposed transposed;
	clfftPrecision       precision;
	cl_context              context;
	double                  forwardScale, backwardScale;
	size_t                  iDist, oDist;
	size_t                  batchsize;

	// Note the device passed to BakePlan, assuming we are baking for one device
	// TODO, change this logic for handling multiple GPUs/devices
	cl_device_id bakeDevice;

	// Disabling devices member, plan has 1-on-1 mapping with single device as identified by bakeDevice
	//	Devices that the user specified in the context passed to the create function
	// std::vector< cl_device_id > devices;

	//	Length of the FFT in each dimension
	std::vector< size_t >	length;

	//	Stride of the FFT in each dimension
	std::vector< size_t >	inStride, outStride;

	//	Hardware Limits
	FFTEnvelope                 envelope;


	// Reserved copy for large 1d, 2d, and 3d plan
	size_t tmpBufSize;
	cl_mem intBuffer;
	bool libCreatedIntBuffer;

	// for RC copies
	size_t	tmpBufSizeRC;
	cl_mem	intBufferRC;

	// for C-to-R transforms that are OUTOFPLACE
	// we need this because the user supplied output buffer is not big enough
	// to hold intermediate results for any problem other than normal 1D
	size_t  tmpBufSizeC2R;
	cl_mem  intBufferC2R;


	size_t  large1D;
	bool    large2D;
	bool	twiddleFront;

	clfftPlanHandle planX;
	clfftPlanHandle planY;
	clfftPlanHandle planZ;

	bool transflag;
	bool transOutHorizontal;
	clfftPlanHandle planTX;
	clfftPlanHandle planTY;
	clfftPlanHandle planTZ; //reserve for 3D transpose

	clfftPlanHandle planRCcopy;
	clfftPlanHandle planCopy;

	// Plan resources
	//
	cl_mem const_buffer;

	// Generator type
	clfftGenerators gen;


	// Real-Complex simple flag
	// if this is set we do real to-and-from full complex using simple algorithm
	// where imaginary of input is set to zero in forward and imaginary not written in backward
	bool RCsimple;

	// Real FFT special flag
	// if this is set it means we are doing the 4th step in the 5-step real FFT breakdown algorithm
	bool realSpecial;
	
	size_t realSpecial_Nr; // this value stores the logical column height (N0) of matrix in the 4th step
	                       // length[1] should be 1 + N0/2

	// User created plan
	bool userPlan;


	// Allocate no extra memory
	bool allOpsInplace;

	// flag to indicate transpose placeness in 2D breakdown
	bool transpose_in_2d_inplace;


	// A flag to say that blocked FFTs are going to be performed
	// It can only be one of these: column to row, row to column or column to column
	// row to row is just the normal case where blocking is not needed
	bool blockCompute;
	BlockComputeType blockComputeType;

	bool hasPreCallback;
	bool hasPostCallback;

	clfftCallbackParam preCallback;
	clfftCallbackParam postCallbackParam;

	cl_mem precallUserData;
	cl_mem postcallUserData;

    clfftPlanHandle plHandle;

    // The action
    FFTAction * action;

    NonSquareTransposeKernelType nonSquareKernelType;
	// sometimes non square matrix are broken down into a number of
	// square matrix during inplace transpose
	// let's call this number transposeMiniBatchSize
	// no user of the library should set its value
	size_t transposeMiniBatchSize;
	NON_SQUARE_KERNEL_ORDER nonSquareKernelOrder;

	FFTPlan ()
	:	baked (false)
	,	dim (CLFFT_1D)
	,	inputLayout (CLFFT_COMPLEX_INTERLEAVED)
	,	outputLayout (CLFFT_COMPLEX_INTERLEAVED)
	,	placeness (CLFFT_INPLACE)
	,   transposed (CLFFT_NOTRANSPOSE)
	,	precision (CLFFT_SINGLE)
	,	context (NULL)
	,	forwardScale (1.0)
	,	backwardScale (1.0)
	,	iDist( 1 ), oDist( 1 )
	,	batchsize (1)
	,   tmpBufSize (0)
	,	intBuffer( NULL )
	,	libCreatedIntBuffer(false)
	,	tmpBufSizeRC (0)
	,	intBufferRC( NULL )
	,	tmpBufSizeC2R (0)
	,	intBufferC2R( NULL )
	,   large1D(0)
	,   large2D(false)
	,	twiddleFront(false)
	,   planX( 0 )
	,   planY( 0 )
	,   planZ( 0 )
	,   transflag(false)
	,	transOutHorizontal(false)
	,	RCsimple(false)
	,	realSpecial(false)
	,	realSpecial_Nr(0)
	,	userPlan(false)
	,	allOpsInplace(false)
	,	transpose_in_2d_inplace(false)
	,	blockCompute(false)
	,	blockComputeType(BCT_C2C)
	,   planTX( 0 )
	,   planTY( 0 )
	,   planTZ( 0 )
	,	planRCcopy(0)
	,	planCopy(0)
	,	const_buffer( NULL )
	,	gen(Stockham)
    ,   action(0)
    ,   nonSquareKernelType(NON_SQUARE_TRANS_PARENT)
	,   transposeMiniBatchSize(1)
	,   nonSquareKernelOrder(NOT_A_TRANSPOSE)
    ,   plHandle(0)
	,   hasPreCallback(false)
	,   hasPostCallback(false)
	{
	};


	size_t ElementSize() const;

	clfftStatus AllocateBuffers ();
	clfftStatus ReleaseBuffers ();

	clfftStatus GetMax1DLength (size_t *longest ) const;

	clfftStatus ConstructAndEnqueueConstantBuffers( cl_command_queue* commQueueFFT );

	clfftStatus GetEnvelope (const FFTEnvelope **) const;
	clfftStatus SetEnvelope ();

	clfftStatus GetMax1DLengthStockham (size_t *longest ) const;

	~FFTPlan ()
	{
		ReleaseBuffers ();

		if (action != NULL)
		{
			delete action;
			action = 0;
		}
	}
};

static bool Is1DPossible(size_t length, size_t large1DThreshold)
{
	if (length > large1DThreshold)
		return false;

	if ( (length%7 == 0) && (length%5 == 0) && (length%3 == 0) )
		return false;

	// radix 11 & 2 is ok, anything else we cannot do in 1 kernel
	if ( (length % 11 == 0) && ((length % 13 == 0) || (length % 7 == 0) || (length % 5 == 0) || (length % 3 == 0)) )
		return false;
	
	// radix 13 & 2 is ok, anything else we cannot do in 1 kernel
	if ( (length % 13 == 0) && ((length % 11 == 0) || (length % 7 == 0) || (length % 5 == 0) || (length % 3 == 0)) )
		return false;

	return true;
}

#endif // AMD_CLFFT_plan_H

