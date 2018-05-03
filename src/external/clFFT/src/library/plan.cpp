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

////////////////////////////////////////////

// clfft.plan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>
#include "private.h"
#include "repo.h"
#include "plan.h"
#include "generator.stockham.h"
#include "../include/convenienceFunctions.h"
#include "action.h"
#include "fft_binary_lookup.h"

using std::vector;

const std::string beginning_of_binary( "<[£_beginning_of_binary_£]>" );
const std::string end_of_binary( "<[£_I_may_be_a_sorry_case,_but_I_don't_write_jokes_in_base_13_£]>" );
const std::string end_of_file( "<[£_You're_off_the_edge_of_the_map,_mate._Here_there_be_monsters_£]>" );

static bool pow235(size_t num, size_t &pow2, size_t &pow3, size_t &pow5)
{
	//a helper function to decide if a number is only radix 2, 3 and 5
	if (num % 2 != 0 && num % 3 != 0 && num % 5 != 0)
		return false;

	while (num > 1)
	{
		if (num % 5 == 0)
		{
			num /= 5;
			pow5++;
			continue;
		}
		if (num % 3 == 0)
		{
			num /= 3;
			pow3++;
			continue;
		}
		if (num % 2 == 0)
		{
			num /= 2;
			pow2++;
			continue;
		}
		return false;
	}
	return true;
}

static bool split1D_for_inplace(size_t num, vector<vector<size_t> > &splitNums, clfftPrecision precision, size_t threshold)
{
	/* a helper function to split big 1D to friendly 2D sizes for inplace transpose kernels
	   currently only radix 2, 3 and 5 are supported
	   the algorithm looks for ways to split up the 1D into 2D such that one of the dimensions is multiples of the other dimension.
	   And this mupliple is radix2, 3 or 5.
	   each splited dimentsion should be further splited until that it is smaller than 4096
	*/
	if (num <= threshold)
		return true;
	if (num % 2 != 0 && num % 3 != 0 && num % 5 != 0)
		return false;

	//let's figure out pow2, pow3 and pow5 such that num = 2^pow2 * 3^pow3 * 5^pow5
	size_t pow2, pow3, pow5;
	pow2 = pow3 = pow5 = 0;
	bool status = pow235(num, pow2, pow3, pow5);
	if (!status)
		return status;

	size_t divide_factor;
	if (pow2 % 2 != 0)
	{
		//pow2 is odd
		if (pow3 % 2 != 0)
		{
			//pow2 and pow3 are odd
			if (pow5 % 2 != 0)
			{
				//pow2, pow3 and pow5 are odd
				//one dimension is 2*3*5 = 30 times bigger than the other dimension
				divide_factor = 2 * 3 * 5;
			}
			else
			{
				//pow2 and pow3 are odd, pow 5 is even
				//one dimension is 2*3 = 6 times bigger than the other dimension
				divide_factor = 2 * 3;
			}
		}
		else
		{
			//pow2 is odd, pow3 is even
			if (pow5 % 2 != 0)
			{
				//pow2, pow5 are odd pow3 is eve
				divide_factor = 2 * 5;
			}
			else
			{
				//pow2 is odd, pow3 and pow5 are even
				divide_factor = 2;
			}

		}
	}
	else
	{
		//pow2 is even
		if (pow3 % 2 != 0)
		{
			//pow3 is odd pow2 is even
			if (pow5 % 2 != 0)
			{
				//pow2 is even, pow3 and pow5 are odd
				divide_factor = 3 * 5;
			}
			else
			{
				//pow2 and pow5 are even, pow3 is odd
				divide_factor = 3;
			}
		}
		else
		{
			//pow2 and are even
			if (pow5 % 2 != 0)
			{
				//pow5 is odd pow2 pow3 is eve
				divide_factor = 5;
			}
			else
			{
				//all even
				divide_factor = 1;
			}

		}
	}
	//add some special cases
	if (num == 2687385600)
		divide_factor = 2 * 2 * 3 * 3;
	if (num == 2916000000)
		divide_factor = 2 * 2 * 3 * 3 * 5 * 5;
	if (num == 3057647616)
		divide_factor = 2 * 2 * 3 * 3;

	num = num / divide_factor;
	//now the remaining num should have even number of pow2, pow3 and pow5 and we can do sqrt
	size_t temp = (size_t)sqrt((double)num);
	vector<size_t> splitVec;
	splitVec.push_back(temp*divide_factor);
	splitVec.push_back(temp);
	splitNums.push_back(splitVec);

	status = status && split1D_for_inplace(temp*divide_factor, splitNums, precision, threshold);
	status = status && split1D_for_inplace(temp, splitNums, precision, threshold);
	return status;

}

// Returns CLFFT_SUCCESS if the fp64 is present, CLFFT_DEVICE_NO_DOUBLE if it is not found.  
clfftStatus checkDevExt( std::string ext, const cl_device_id &device )
{
	size_t deviceExtSize	= 0;
	OPENCL_V( ::clGetDeviceInfo( device, CL_DEVICE_EXTENSIONS, 0, NULL, &deviceExtSize ),
		"Getting CL_DEVICE_EXTENSIONS Platform Info string size ( ::clGetDeviceInfo() )" );

	std::vector< char > szDeviceExt( deviceExtSize );
	OPENCL_V( ::clGetDeviceInfo( device, CL_DEVICE_EXTENSIONS, deviceExtSize, &szDeviceExt[ 0 ], NULL ),
		"Getting CL_DEVICE_EXTENSIONS Platform Info string ( ::clGetDeviceInfo() )" );

	std::string strDeviceExt = &szDeviceExt[ 0 ];

	if( strDeviceExt.find( ext.c_str( ), 0 ) == std::string::npos )
		return CLFFT_DEVICE_NO_DOUBLE;


	return CLFFT_SUCCESS;
}

clfftStatus	clfftCreateDefaultPlanInternal( clfftPlanHandle* plHandle, cl_context context, const clfftDim dim,
						const size_t* clLengths )
{
	if( clLengths == NULL )
		return CLFFT_INVALID_HOST_PTR;

	size_t lenX = 1, lenY = 1, lenZ = 1;

	switch( dim )
	{
		case CLFFT_1D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 )
				return CLFFT_INVALID_ARG_VALUE;

			if( !IsASupportedLength( clLengths[ DimX ] ) )
			{
				return CLFFT_NOTIMPLEMENTED;
			}

			lenX = clLengths[ DimX ];
		}
			break;
		case CLFFT_2D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 || clLengths[ DimY ] == 0 )
				return CLFFT_INVALID_ARG_VALUE;

			if( !IsASupportedLength( clLengths[ DimX ] ) || !IsASupportedLength( clLengths[ DimY ] ) )
			{
				return CLFFT_NOTIMPLEMENTED;
			}

			lenX = clLengths[ DimX ];
			lenY = clLengths[ DimY ];
		}
			break;
		case CLFFT_3D:
		{
			//	Minimum length size is 1
			if( clLengths[ DimX ] == 0 || clLengths[ DimY ] == 0 || clLengths[ DimZ ] == 0 )
				return CLFFT_INVALID_ARG_VALUE;

			if( !IsASupportedLength( clLengths[ DimX ] ) || !IsASupportedLength( clLengths[ DimY ] ) ||
				!IsASupportedLength( clLengths[ DimZ ] ))
			{
				return CLFFT_NOTIMPLEMENTED;
			}

			lenX = clLengths[ DimX ];
			lenY = clLengths[ DimY ];
			lenZ = clLengths[ DimZ ];
		}
			break;
		default:
			return CLFFT_NOTIMPLEMENTED;
			break;
	}

	FFTPlan *fftPlan = NULL;
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	OPENCL_V( fftRepo.createPlan( plHandle, fftPlan ), _T( "fftRepo.insertPlan failed" ) );

	fftPlan->baked			= false;
	fftPlan->dim			= dim;
	fftPlan->placeness		= CLFFT_INPLACE;
	fftPlan->inputLayout	= CLFFT_COMPLEX_INTERLEAVED;
	fftPlan->outputLayout	= CLFFT_COMPLEX_INTERLEAVED;
	fftPlan->precision		= CLFFT_SINGLE;
	fftPlan->context		= context;
	fftPlan->forwardScale	= 1.0;
	fftPlan->backwardScale	= 1.0 / static_cast< double >( lenX * lenY * lenZ );
	fftPlan->batchsize		= 1;
	fftPlan->gen			= Stockham; //default setting

	OPENCL_V(fftPlan->SetEnvelope(), _T("SetEnvelope failed"));

	clRetainContext( fftPlan->context );

#if 0
	/////////////////////////////////////////////////////////////////
	// Detect OpenCL devices
	/////////////////////////////////////////////////////////////////
	// First, get the size of device list data
	size_t deviceListSize;
	OPENCL_V( ::clGetContextInfo( context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize ),
		"Getting device array size ( ::clGetContextInfo() )" );

	//	Allocate memory for the devices
	fftPlan->devices.resize( deviceListSize / sizeof( cl_device_id ) );

	/* Now, get the device list data */
	OPENCL_V( ::clGetContextInfo( context, CL_CONTEXT_DEVICES, deviceListSize, &fftPlan->devices[ 0 ], NULL ),
		"Getting device array ( ::clGetContextInfo() )" );
#endif

	//	Need to devise a way to generate better names
	tstringstream	tstream;
	tstream << _T( "plan_" ) << *plHandle;

	lockRAII* planLock	= NULL;
	OPENCL_V( fftRepo.getPlan( *plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	planLock->setName( tstream.str( ) );

	//	Set the lengths and default strides/pitches depending on the dim that the user passes to us
	switch( dim )
	{
		case CLFFT_1D:
		{
			fftPlan->length.push_back( lenX );
			fftPlan->inStride.push_back( 1 );
			fftPlan->outStride.push_back( 1 );
			fftPlan->iDist		= lenX;
			fftPlan->oDist		= lenX;
		}
			break;
		case CLFFT_2D:
		{
			fftPlan->length.push_back( lenX );
			fftPlan->length.push_back( lenY );
			fftPlan->inStride.push_back( 1 );
			fftPlan->inStride.push_back( lenX );
			fftPlan->outStride.push_back( 1 );
			fftPlan->outStride.push_back( lenX );
			fftPlan->iDist		= lenX*lenY;
			fftPlan->oDist		= lenX*lenY;
		}
			break;
		case CLFFT_3D:
		{
			fftPlan->length.push_back( lenX );
			fftPlan->length.push_back( lenY );
			fftPlan->length.push_back( lenZ );
			fftPlan->inStride.push_back( 1 );
			fftPlan->inStride.push_back( lenX );
			fftPlan->inStride.push_back( lenX*lenY );
			fftPlan->outStride.push_back( 1 );
			fftPlan->outStride.push_back( lenX );
			fftPlan->outStride.push_back( lenX*lenY );
			fftPlan->iDist		= lenX*lenY*lenZ;
			fftPlan->oDist		= lenX*lenY*lenZ;
		}
			break;
	}

	fftPlan->plHandle = *plHandle;

	return	CLFFT_SUCCESS;
}

// This external entry-point should not be called from within the library. Use clfftCreateDefaultPlanInternal instead.
clfftStatus	clfftCreateDefaultPlan( clfftPlanHandle* plHandle, cl_context context, const clfftDim dim,
						const size_t* clLengths )
{
	clfftStatus ret = clfftCreateDefaultPlanInternal(plHandle, context, dim, clLengths);

	if(ret == CLFFT_SUCCESS)
	{
		FFTRepo& fftRepo	= FFTRepo::getInstance( );
		FFTPlan *fftPlan = NULL;
		lockRAII* planLock	= NULL;
		OPENCL_V( fftRepo.getPlan( *plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );

		fftPlan->userPlan = true;
	}
	
	return ret;

}

std::string getKernelName(const clfftGenerators gen, const clfftPlanHandle plHandle, bool withPlHandle)
{
    //	Logic to define a sensible filename
    const std::string kernelPrefix( "clfft.kernel." );
    std::string generatorName;
    std::stringstream kernelPath;

    switch( gen )
    {

    case Stockham:			    generatorName = "Stockham"; break;
	case Transpose_GCN:		    generatorName = "Transpose"; break;
	case Transpose_SQUARE:	    generatorName = "Transpose"; break;
    case Transpose_NONSQUARE:	generatorName = "TransposeNonSquare"; break;
	case Copy:				    generatorName = "Copy"; break;

    }

    kernelPath << kernelPrefix << generatorName ;

    if (withPlHandle)
        kernelPath << plHandle;

    kernelPath << ".cl";

    return kernelPath.str();
}


clfftStatus selectAction(FFTPlan * fftPlan, FFTAction *& action, cl_command_queue* commQueueFFT)
{
    // set the action we are baking a leaf
    clfftStatus err;
    
    switch (fftPlan->gen)
    {
    case Stockham:  
		{
			// Instantiate the default stockham generator
			action = new FFTGeneratedStockhamAction (fftPlan->plHandle, fftPlan, *commQueueFFT, err);
			OPENCL_V( err, "FFTGeneratedStockhamAction() failed");
		}
		break;

    case Transpose_GCN: 
		{
			action = new FFTGeneratedTransposeGCNAction(fftPlan->plHandle, fftPlan, *commQueueFFT, err);
			OPENCL_V( err, "FFTGeneratedTransposeGCNAction() failed");
		}
		break;


    case Copy:
		{
			action = new FFTGeneratedCopyAction     (fftPlan->plHandle, fftPlan, *commQueueFFT, err);
			OPENCL_V( err, "FFTGeneratedCopyAction() failed");
		}
		break;

    default:
		{
			assert(false);
			OPENCL_V( CLFFT_NOTIMPLEMENTED, "selectAction() failed");
		}
    }

	return CLFFT_SUCCESS;
}


inline size_t PrecisionWidth(clfftPrecision pr)
{
	switch(pr)
	{
	case CLFFT_SINGLE:	return 1;
	case CLFFT_DOUBLE:	return 2;
	default:		assert(false); return 1;
	}
}



clfftStatus	clfftBakePlan( clfftPlanHandle plHandle, cl_uint numQueues, cl_command_queue* commQueueFFT,
							void (CL_CALLBACK *pfn_notify)( clfftPlanHandle plHandle, void *user_data ), void* user_data )
{
	//	We do not currently support multi-GPU transforms
	if( numQueues > 1 )
		return CLFFT_NOTIMPLEMENTED;

	//	Notification mechanism is not set up yet; BakePlan can be called recursively to decompose higher dimension FFT's into
	//	arrays of 1d transforms, and this must be implemented to make only a single callback to the user.
	if( pfn_notify != NULL )
		return CLFFT_NOTIMPLEMENTED;

	if( user_data != NULL )
		return CLFFT_NOTIMPLEMENTED;

	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftBakePlan" ) );

	// if we have already baked the plan and nothing has changed since, we're done here
	if( fftPlan->baked == true )
	{
		return CLFFT_SUCCESS;
	}

	// Store the device for which we are baking
	clGetCommandQueueInfo(*commQueueFFT, CL_QUEUE_DEVICE, sizeof(cl_device_id), &fftPlan->bakeDevice, NULL);

	//find product of lengths
	size_t maxLengthInAnyDim = 1;
	switch(fftPlan->dim)
	{
		case CLFFT_3D: maxLengthInAnyDim = maxLengthInAnyDim > fftPlan->length[DimZ] ? maxLengthInAnyDim : fftPlan->length[DimZ];
		case CLFFT_2D: maxLengthInAnyDim = maxLengthInAnyDim > fftPlan->length[DimY] ? maxLengthInAnyDim : fftPlan->length[DimY];
		case CLFFT_1D: maxLengthInAnyDim = maxLengthInAnyDim > fftPlan->length[DimX] ? maxLengthInAnyDim : fftPlan->length[DimX];
	}

	const bool rc = (fftPlan->inputLayout == CLFFT_REAL) || (fftPlan->outputLayout == CLFFT_REAL);

	// upper bounds on transfrom lengths - address this in the next release
	size_t SP_MAX_LEN = 1 << 24;
	size_t DP_MAX_LEN = 1 << 22;
	if((fftPlan->precision == CLFFT_SINGLE) && (maxLengthInAnyDim > SP_MAX_LEN) && rc) return CLFFT_NOTIMPLEMENTED;
	if((fftPlan->precision == CLFFT_DOUBLE) && (maxLengthInAnyDim > DP_MAX_LEN) && rc) return CLFFT_NOTIMPLEMENTED;


	// release buffers, as these will be created only in EnqueueTransform
	if( NULL != fftPlan->intBuffer ) { OPENCL_V( clReleaseMemObject( fftPlan->intBuffer ), _T( "Failed to release internal temporary buffer" ) ); fftPlan->intBuffer = NULL; }
	if( NULL != fftPlan->intBufferRC ) { OPENCL_V( clReleaseMemObject( fftPlan->intBufferRC ), _T( "Failed to release internal temporary buffer" ) ); fftPlan->intBufferRC = NULL; }
	if( NULL != fftPlan->intBufferC2R ) { OPENCL_V( clReleaseMemObject( fftPlan->intBufferC2R ), _T( "Failed to release internal temporary buffer" ) ); fftPlan->intBufferC2R = NULL; }


    if( fftPlan->userPlan ) // confirm it is top-level plan (user plan)
	{
		if(fftPlan->placeness == CLFFT_INPLACE)
		{
			if( (fftPlan->inputLayout == CLFFT_HERMITIAN_PLANAR) || (fftPlan->outputLayout == CLFFT_HERMITIAN_PLANAR) )
				return CLFFT_INVALID_PLAN;
		}

		// Make sure strides & distance are same for C-C transforms
		if(fftPlan->placeness == CLFFT_INPLACE)
		{
			if( (fftPlan->inputLayout != CLFFT_REAL) && (fftPlan->outputLayout != CLFFT_REAL) )
			{
				// check strides
				for(size_t i=0; i<fftPlan->dim; i++)
					if(fftPlan->inStride[i] != fftPlan->outStride[i])
						return CLFFT_INVALID_PLAN;

				// check distance
				if(fftPlan->iDist != fftPlan->oDist)
					return CLFFT_INVALID_PLAN;
			}
		}
	}

	if(fftPlan->gen == Copy)
	{
        clfftStatus err;
        fftPlan->action = new FFTGeneratedCopyAction(plHandle, fftPlan, *commQueueFFT, err);
        OPENCL_V( err, _T( "FFTGeneratedCopyAction() failed" ) );
		fftPlan->baked		= true;
		return	CLFFT_SUCCESS;
	}


	if( fftPlan->userPlan )
	{
		//	If the user specifies double precision, check that the device supports double precision first
		if( fftPlan->precision == CLFFT_DOUBLE || fftPlan->precision == CLFFT_DOUBLE_FAST )
		{
			clfftStatus retAmdFp64 = checkDevExt( "cl_amd_fp64", fftPlan->bakeDevice );
			if( retAmdFp64 != CLFFT_SUCCESS )
			{
				//	If AMD's extention is not supported, check for Khronos extention
				clfftStatus retKhrFp64 = checkDevExt( "cl_khr_fp64", fftPlan->bakeDevice );
				if( retKhrFp64 != CLFFT_SUCCESS )
					return retKhrFp64;
			}
		}
	}

	// Compress the plan by discarding length '1' dimensions
	// decision to pick generator
	if( fftPlan->userPlan && !rc ) // confirm it is top-level plan (user plan)
	{
		size_t dmnsn = fftPlan->dim;
		bool pow2flag = true;

		// switch case flows with no 'break' statements
		switch(fftPlan->dim)
		{
		case CLFFT_3D:

			if(fftPlan->length[DimZ] == 1)
			{
				dmnsn -= 1;
				fftPlan-> inStride.erase(fftPlan-> inStride.begin() + 2);
				fftPlan->outStride.erase(fftPlan->outStride.begin() + 2);
				fftPlan->   length.erase(fftPlan->   length.begin() + 2);
			}
			else
			{
				if( !IsPo2(fftPlan->length[DimZ])) pow2flag=false;
			}
		case CLFFT_2D:

			if(fftPlan->length[DimY] == 1)
			{
				dmnsn -= 1;
				fftPlan-> inStride.erase(fftPlan-> inStride.begin() + 1);
				fftPlan->outStride.erase(fftPlan->outStride.begin() + 1);
				fftPlan->   length.erase(fftPlan->   length.begin() + 1);
			}
			else
			{
				if( !IsPo2(fftPlan->length[DimY])) pow2flag=false;
			}

		case CLFFT_1D:

			if( (fftPlan->length[DimX] == 1) && (dmnsn > 1) )
			{
				dmnsn -= 1;
				fftPlan-> inStride.erase(fftPlan-> inStride.begin());
				fftPlan->outStride.erase(fftPlan->outStride.begin());
				fftPlan->   length.erase(fftPlan->   length.begin());
			}
			else
			{
				if( !IsPo2(fftPlan->length[DimX])) pow2flag=false;
			}
		}

		fftPlan->dim = (clfftDim)dmnsn;
	}

	// first time check transposed
	if (fftPlan->transposed != CLFFT_NOTRANSPOSE && fftPlan->dim != CLFFT_2D &&
		fftPlan->dim == fftPlan->length.size())
		return CLFFT_TRANSPOSED_NOTIMPLEMENTED;

	//	The largest vector we can transform in a single pass
	//	depends on the GPU caps -- especially the amount of LDS
	//	available
	//
	size_t Large1DThreshold = 0;


	OPENCL_V(fftPlan->GetMax1DLength (&Large1DThreshold), _T("GetMax1DLength failed"));
	BUG_CHECK (Large1DThreshold > 1);

	//	Verify that the data passed to us is packed
	switch( fftPlan->dim )
	{
	case CLFFT_1D:
		{
			if ( !Is1DPossible(fftPlan->length[0], Large1DThreshold) )
			{
				size_t clLengths[] = { 1, 1, 0 };
				size_t in_1d, in_x, count;

				BUG_CHECK (IsPo2 (Large1DThreshold))


				if( IsPo2(fftPlan->length[0]) )
				{
					// Enable block compute under these conditions
					if( (fftPlan->inStride[0] == 1) && (fftPlan->outStride[0] == 1) && !rc
						&& (fftPlan->length[0] <= 262144/PrecisionWidth(fftPlan->precision)) && (fftPlan->length.size() <= 1)
						&& (!clfftGetRequestLibNoMemAlloc() || (fftPlan->placeness == CLFFT_OUTOFPLACE)) )
					{
						fftPlan->blockCompute = true;

						if(1 == PrecisionWidth(fftPlan->precision))
						{
							switch(fftPlan->length[0])
							{
							case 8192:		clLengths[1] = 64;	break;
							case 16384:		clLengths[1] = 64;	break;
							case 32768:		clLengths[1] = 128;	break;
							case 65536:		clLengths[1] = 256;	break;
							case 131072:	clLengths[1] = 64;	break;
							case 262144:	clLengths[1] = 64;	break;
							case 524288:	clLengths[1] = 256; break;
							case 1048576:	clLengths[1] = 256; break;
							default:		assert(false);
							}
						}
						else
						{
							switch(fftPlan->length[0])
							{
							case 4096:		clLengths[1] = 64;	break;
							case 8192:		clLengths[1] = 64;	break;
							case 16384:		clLengths[1] = 64;	break;
							case 32768:		clLengths[1] = 128;	break;
							case 65536:		clLengths[1] = 64;	break;
							case 131072:	clLengths[1] = 64;	break;
							case 262144:	clLengths[1] = 128;	break;
							case 524288:	clLengths[1] = 256; break;
							default:		assert(false);
							}
						}
					}
					else
					{
						if( clfftGetRequestLibNoMemAlloc() && !rc && (fftPlan->placeness == CLFFT_INPLACE) )
						{
							in_x = BitScanF(fftPlan->length[0]);
							in_x /= 2;
							clLengths[1] = (size_t)1 << in_x;
						}
						else if( fftPlan->length[0] > (Large1DThreshold * Large1DThreshold) )
						{
							clLengths[1] = fftPlan->length[0] / Large1DThreshold;
						}
						else
						{
							in_1d = BitScanF (Large1DThreshold);	// this is log2(LARGE1D_THRESHOLD)
							in_x  = BitScanF (fftPlan->length[0]);	// this is log2(length)
							BUG_CHECK (in_1d > 0)
							count = in_x/in_1d;
							if (count*in_1d < in_x)
							{
								count++;
								in_1d = in_x / count;
								if (in_1d * count < in_x) in_1d++;
							}
							clLengths[1] = (size_t)1 << in_1d;
						}
					}
				}
				else
				{
					// This array must be kept sorted in the ascending order

					size_t supported[] = {	1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 20, 21, 22, 24,
											25, 26, 27, 28, 30, 32, 33, 35, 36, 39, 40, 42, 44, 45, 48, 49, 50, 52, 54,
											55, 56, 60, 63, 64, 65, 66, 70, 72, 75, 77, 78, 80, 81, 84, 88, 90, 91, 96,
											98, 99, 100, 104, 105, 108, 110, 112, 117, 120, 121, 125, 126, 128, 130, 132,
											135, 140, 143, 144, 147, 150, 154, 156, 160, 162, 165, 168, 169, 175, 176,
											180, 182, 189, 192, 195, 196, 198, 200, 208, 210, 216, 220, 224, 225, 231,
											234, 240, 242, 243, 245, 250, 252, 256, 260, 264, 270, 273, 275, 280, 286,
											288, 294, 297, 300, 308, 312, 315, 320, 324, 325, 330, 336, 338, 343, 350,
											351, 352, 360, 363, 364, 375, 378, 384, 385, 390, 392, 396, 400, 405, 416,
											420, 429, 432, 440, 441, 448, 450, 455, 462, 468, 480, 484, 486, 490, 495,
											500, 504, 507, 512, 520, 525, 528, 539, 540, 546, 550, 560, 567, 572, 576,
											585, 588, 594, 600, 605, 616, 624, 625, 630, 637, 640, 648, 650, 660, 672,
											675, 676, 686, 693, 700, 702, 704, 715, 720, 726, 728, 729, 735, 750, 756,
											768, 770, 780, 784, 792, 800, 810, 819, 825, 832, 840, 845, 847, 858, 864,
											875, 880, 882, 891, 896, 900, 910, 924, 936, 945, 960, 968, 972, 975, 980,
											990, 1000, 1001, 1008, 1014, 1024, 1029, 1040, 1050, 1053, 1056, 1078, 1080,
											1089, 1092, 1100, 1120, 1125, 1134, 1144, 1152, 1155, 1170, 1176, 1183, 1188,
											1200, 1210, 1215, 1225, 1232, 1248, 1250, 1260, 1274, 1280, 1287, 1296, 1300,
											1320, 1323, 1331, 1344, 1350, 1352, 1365, 1372, 1375, 1386, 1400, 1404, 1408,
											1430, 1440, 1452, 1456, 1458, 1470, 1485, 1500, 1512, 1521, 1536, 1540, 1560,
											1568, 1573, 1575, 1584, 1600, 1617, 1620, 1625, 1638, 1650, 1664, 1680, 1690,
											1694, 1701, 1715, 1716, 1728, 1750, 1755, 1760, 1764, 1782, 1792, 1800, 1815,
											1820, 1848, 1859, 1872, 1875, 1890, 1911, 1920, 1925, 1936, 1944, 1950, 1960,
											1980, 2000, 2002, 2016, 2025, 2028, 2048, 2058, 2079, 2080, 2100, 2106, 2112,
											2145, 2156, 2160, 2178, 2184, 2187, 2197, 2200, 2205, 2240, 2250, 2268, 2275,
											2288, 2304, 2310, 2340, 2352, 2366, 2376, 2400, 2401, 2420, 2430, 2450, 2457,
											2464, 2475, 2496, 2500, 2520, 2535, 2541, 2548, 2560, 2574, 2592, 2600, 2625,
											2640, 2646, 2662, 2673, 2688, 2695, 2700, 2704, 2730, 2744, 2750, 2772, 2800,
											2808, 2816, 2835, 2860, 2880, 2904, 2912, 2916, 2925, 2940, 2970, 3000, 3003,
											3024, 3025, 3042, 3072, 3080, 3087, 3120, 3125, 3136, 3146, 3150, 3159, 3168,
											3185, 3200, 3234, 3240, 3250, 3267, 3276, 3300, 3328, 3360, 3375, 3380, 3388,
											3402, 3430, 3432, 3456, 3465, 3500, 3510, 3520, 3528, 3549, 3564, 3575, 3584,
											3600, 3630, 3640, 3645, 3675, 3696, 3718, 3744, 3750, 3773, 3780, 3822, 3840,
											3850, 3861, 3872, 3888, 3900, 3920, 3960, 3969, 3993, 4000, 4004, 4032, 4050,
											4056, 4095, 4096 };

					size_t lenSupported = sizeof(supported)/sizeof(supported[0]);
					size_t maxFactoredLength = (supported[lenSupported-1] < Large1DThreshold) ? supported[lenSupported-1] : Large1DThreshold;

					size_t halfPowerLength = (size_t)1 << ( (StockhamGenerator::CeilPo2(fftPlan->length[0]) + 1) / 2 );
					size_t factoredLengthStart =  (halfPowerLength < maxFactoredLength) ? halfPowerLength : maxFactoredLength;

					size_t indexStart = 0;
					while(supported[indexStart] < factoredLengthStart) indexStart++;

					for(size_t i = indexStart; i >= 1; i--)
					{
						if( fftPlan->length[0] % supported[i] == 0 )
						{
							if (Is1DPossible(supported[i], Large1DThreshold))
							{
								clLengths[1] = supported[i];
								break;
							}
						}
					}
				}
				// add some special cases
				/*
				if (fftPlan->length[0] == 10000)
					clLengths[1] = 100;//100 x 100
				if (fftPlan->length[0] == 100000)
					clLengths[1] = 100;//100 x 1,000
				if (fftPlan->length[0] == 10000000)
					clLengths[1] = 1000;//1,000 x 10,000
				if (fftPlan->length[0] == 100000000)
					clLengths[1] = 10000;//10,000 x 10,000
				if (fftPlan->length[0] == 1000000000)
					clLengths[1] = 10000;//10,000 x 100,000
				
				if (fftPlan->length[0] == 3099363912)
					clLengths[1] = 78732;//39366 x 78732
				if (fftPlan->length[0] == 39366)
					clLengths[1] = 81;//81*486
				if (fftPlan->length[0] == 78732)
					clLengths[1] = 162;//162*486
				if (fftPlan->length[0] == 354294)
					clLengths[1] = 243;
				*/
				size_t threshold = 4096;
				if (fftPlan->precision == CLFFT_DOUBLE)
					threshold = 2048;
				if (clfftGetRequestLibNoMemAlloc() &&
					fftPlan->placeness == CLFFT_INPLACE &&
					(fftPlan->inputLayout == fftPlan->outputLayout)
					&& fftPlan->length[0] > threshold)
				{
					//for inplace fft with inplace transpose, the split logic is different
					vector<vector<size_t> > splitNums;
					bool implemented = split1D_for_inplace(fftPlan->length[0], splitNums, fftPlan->precision, threshold);
					if (implemented)
						clLengths[1] = splitNums[0][0];
				}

				clLengths[0] = fftPlan->length[0]/clLengths[1];

                // Start of block where transposes are generated; 1D FFT
				while (1 && (fftPlan->inputLayout != CLFFT_REAL) && (fftPlan->outputLayout != CLFFT_REAL))
				{
					if (fftPlan->length[0] <= Large1DThreshold) break;

					if (fftPlan->inStride[0] != 1 || fftPlan->outStride[0] != 1) break;

					if ( IsPo2(fftPlan->length[0]) &&
						 (fftPlan->length[0] <= 262144/PrecisionWidth(fftPlan->precision)) && (fftPlan->length.size() <= 1) &&
						 (!clfftGetRequestLibNoMemAlloc() || (fftPlan->placeness == CLFFT_OUTOFPLACE)) ) break;

					if ( clLengths[0]<=32 && clLengths[1]<=32) break;


					size_t biggerDim = clLengths[0] > clLengths[1] ? clLengths[0] : clLengths[1];
					size_t smallerDim = biggerDim == clLengths[0] ? clLengths[1] : clLengths[0];
					size_t padding = 0;
					if( (smallerDim % 64 == 0) || (biggerDim % 64 == 0) )
						padding = 64;

					clfftGenerators transGen = Transpose_GCN;
					
					size_t dim_ratio = biggerDim / smallerDim;
					size_t dim_residue = biggerDim % smallerDim;
					//    If this is an in-place transform the
					//    input and output layout, dimensions and strides
					//    *MUST* be the same.
					//
					bool inStrideEqualsOutStride = true;
					for (size_t u = fftPlan->inStride.size(); u-- > 0; ) {
						if (fftPlan->inStride[u] != fftPlan->outStride[u])
						{
							inStrideEqualsOutStride = false;
							break;
						}
					}
					//packed data is required for inplace transpose
					bool isDataPacked = true;
					for (size_t u = 0; u < fftPlan->inStride.size(); u++)
					{
						if (u == 0)
						{
							if (fftPlan->inStride[0] == 1)
								continue;
							else
							{
								isDataPacked = false;
								break;
							}
						}
						else
						{
							size_t packDataSize = 1;
							for (size_t i = 0; i < u; i++)
								packDataSize *= fftPlan->length[i];
							if (fftPlan->inStride[u] == packDataSize)
								continue;
							else
							{
								isDataPacked = false;
								break;
							}
						}
					}
					if (clfftGetRequestLibNoMemAlloc() &&
						dim_residue == 0 &&
						((dim_ratio % 2 == 0) ||
						 (dim_ratio % 3 == 0) ||
						 (dim_ratio % 5 == 0) ||
						 (dim_ratio % 10 == 0)) &&
						 fftPlan->placeness == CLFFT_INPLACE &&
						 (fftPlan->inputLayout == fftPlan->outputLayout) &&
						 (inStrideEqualsOutStride) && (isDataPacked))
					{
						padding = 0;
						fftPlan->allOpsInplace = true;
						transGen = Transpose_NONSQUARE;
						//std::cout << "Transpose_NONSQUARE" << std::endl;
					}

					if( clfftGetRequestLibNoMemAlloc() &&
						(clLengths[0] == clLengths[1]) &&
						fftPlan->placeness == CLFFT_INPLACE )
					{
						padding = 0;
						fftPlan->allOpsInplace = true;
						transGen = Transpose_SQUARE;
					}

					if (fftPlan->tmpBufSize != 0)
						padding = 0;

					if ( (fftPlan->tmpBufSize==0 ) && !fftPlan->allOpsInplace)
					{
						fftPlan->tmpBufSize = (smallerDim + padding) * biggerDim *
							fftPlan->batchsize * fftPlan->ElementSize();

						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSize *= fftPlan->length[index];
						}
					}

					//Transpose
					//Input --> tmp buffer
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan Large1d transpose 1 failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->placeness     = fftPlan->allOpsInplace ? CLFFT_INPLACE : CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->inputLayout   = fftPlan->inputLayout;
					trans1Plan->outputLayout  = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					trans1Plan->inStride[0]   = fftPlan->inStride[0];
					trans1Plan->inStride[1]   = clLengths[0];
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = clLengths[1] + padding;
					trans1Plan->iDist         = fftPlan->iDist;
					trans1Plan->oDist         = clLengths[0] * trans1Plan->outStride[1];
					trans1Plan->gen           = transGen;
					trans1Plan->transflag     = true;

					if (trans1Plan->gen == Transpose_NONSQUARE || trans1Plan->gen == Transpose_SQUARE)// inplace transpose
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							//trans1Plan->length.push_back(fftPlan->length[index]);
							/*
							replacing the line above with the two lines below since:
							fftPlan is still 1D, thus the broken down transpose should be 2D not 3D
							the batchSize for the transpose should increase accordingly.
							the iDist should decrease accordingly. Push back to length will cause a 3D transpose
							*/
							trans1Plan->batchsize = trans1Plan->batchsize * fftPlan->length[index];
							trans1Plan->iDist = trans1Plan->iDist / fftPlan->length[index];

							trans1Plan->inStride.push_back(fftPlan->inStride[index]);
							trans1Plan->outStride.push_back(trans1Plan->oDist);
							trans1Plan->oDist *= fftPlan->length[index];
						}
					}
					else
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							trans1Plan->length.push_back(fftPlan->length[index]);

							trans1Plan->inStride.push_back(fftPlan->inStride[index]);
							trans1Plan->outStride.push_back(trans1Plan->oDist);
							trans1Plan->oDist *= fftPlan->length[index];
						}
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						trans1Plan->hasPreCallback = true;
						trans1Plan->preCallback = fftPlan->preCallback;
						trans1Plan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans1 plan failed" ) );

					//Row transform
					//tmp->output
					//size clLengths[1], batch clLengths[0], with length[0] twiddle factor multiplication
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[1] ),
						_T( "CreateDefaultPlan Large1d column failed" ) );

					FFTPlan* row1Plan	= NULL;
					lockRAII* row1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, row1Plan, row1Lock ), _T( "fftRepo.getPlan failed" ) );

					row1Plan->placeness     = fftPlan->allOpsInplace ? CLFFT_INPLACE : CLFFT_OUTOFPLACE;
					row1Plan->precision     = fftPlan->precision;
					row1Plan->forwardScale  = 1.0f;
					row1Plan->backwardScale = 1.0f;
					row1Plan->tmpBufSize    = 0;
					row1Plan->batchsize     = fftPlan->batchsize;

					row1Plan->gen			= fftPlan->gen;
					row1Plan->envelope		= fftPlan->envelope;

					// twiddling is done in row2
					row1Plan->large1D		= 0;

					row1Plan->length.push_back(clLengths[0]);
					row1Plan->inputLayout   = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					row1Plan->outputLayout  = fftPlan->outputLayout;
					row1Plan->inStride[0]   = 1;
					row1Plan->outStride[0]  = fftPlan->outStride[0];
					row1Plan->inStride.push_back(clLengths[1]+padding);
					row1Plan->outStride.push_back(clLengths[1]);
					row1Plan->iDist         = clLengths[0] * row1Plan->inStride[1];
					row1Plan->oDist         = fftPlan->oDist;

					for (size_t index = 1; index < fftPlan->length.size(); index++)
					{
						row1Plan->length.push_back(fftPlan->length[index]);
						row1Plan->inStride.push_back(row1Plan->iDist);
						row1Plan->iDist *= fftPlan->length[index];
						row1Plan->outStride.push_back(fftPlan->outStride[index]);
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d first row plan failed" ) );

					//Transpose 2
					//Output --> tmp buffer
					clLengths[2] = clLengths[0];
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, &clLengths[1] ),
						_T( "CreateDefaultPlan Large1d transpose 2 failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTY, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->placeness     = fftPlan->allOpsInplace ? CLFFT_INPLACE : CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->inputLayout   = fftPlan->outputLayout;
					trans2Plan->outputLayout  = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					trans2Plan->inStride[0]   = fftPlan->outStride[0];
					trans2Plan->inStride[1]   = clLengths[1];
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = clLengths[0] + padding;
					trans2Plan->iDist         = fftPlan->oDist;
					trans2Plan->oDist         = clLengths[1] * trans2Plan->outStride[1];
                    trans2Plan->gen           = transGen;

					//if (transGen != Transpose_NONSQUARE)//twiddle
						trans2Plan->large1D		  = fftPlan->length[0];

					trans2Plan->transflag     = true;

					if (trans2Plan->gen == Transpose_NONSQUARE || trans2Plan->gen == Transpose_SQUARE)// inplace transpose
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							//trans2Plan->length.push_back(fftPlan->length[index]);
							/*
							replacing the line above with the two lines below since:
							fftPlan is still 1D, thus the broken down transpose should be 2D not 3D
							the batchSize for the transpose should increase accordingly.
							the iDist should decrease accordingly. Push back to length will cause a 3D transpose
							*/
							trans2Plan->batchsize = trans2Plan->batchsize * fftPlan->length[index];
							trans2Plan->iDist = trans2Plan->iDist / fftPlan->length[index];
							trans2Plan->inStride.push_back(fftPlan->outStride[index]);
							trans2Plan->outStride.push_back(trans2Plan->oDist);
							trans2Plan->oDist *= fftPlan->length[index];
						}
					}
					else
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							trans2Plan->length.push_back(fftPlan->length[index]);

							trans2Plan->inStride.push_back(fftPlan->outStride[index]);
							trans2Plan->outStride.push_back(trans2Plan->oDist);
							trans2Plan->oDist *= fftPlan->length[index];
						}
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans2 plan failed" ) );

					//Row transform 2
					//tmp->tmp
					//size clLengths[0], batch clLengths[1]
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &clLengths[0] ),
						_T( "CreateDefaultPlan Large1d second row plan failed" ) );

					FFTPlan* row2Plan	= NULL;
					lockRAII* row2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, row2Plan, row2Lock ), _T( "fftRepo.getPlan failed" ) );

					row2Plan->placeness     = CLFFT_INPLACE;
					row2Plan->precision     = fftPlan->precision;
					row2Plan->forwardScale  = fftPlan->forwardScale;
					row2Plan->backwardScale = fftPlan->backwardScale;
					row2Plan->tmpBufSize    = 0;
					row2Plan->batchsize     = fftPlan->batchsize;

					row2Plan->gen			= fftPlan->gen;
					row2Plan->envelope		= fftPlan->envelope;


					row2Plan->length.push_back(clLengths[1]);
					row2Plan->inputLayout   = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					row2Plan->outputLayout  = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					row2Plan->inStride[0]   = 1;
					row2Plan->outStride[0]  = 1;
					row2Plan->inStride.push_back(clLengths[0] + padding);
					row2Plan->outStride.push_back(clLengths[0] + padding);
					row2Plan->iDist         = clLengths[1] * row2Plan->inStride[1];
					row2Plan->oDist         = clLengths[1] * row2Plan->outStride[1];

					for (size_t index = 1; index < fftPlan->length.size(); index++)
					{
						row2Plan->length.push_back(fftPlan->length[index]);
						row2Plan->inStride.push_back(row2Plan->iDist);
						row2Plan->outStride.push_back(row2Plan->oDist);
						row2Plan->iDist *= fftPlan->length[index];
						row2Plan->oDist *= fftPlan->length[index];
					}
					
					//if (transGen != Transpose_NONSQUARE)//twiddle in transform
					//{
					//	row2Plan->large1D = fftPlan->length[0];
					//	row2Plan->twiddleFront = true;
					//}

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d second row plan failed" ) );

					//Transpose 3
					//tmp --> output
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTZ, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan Large1d transpose 3 failed" ) );

					FFTPlan* trans3Plan	= NULL;
					lockRAII* trans3Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTZ, trans3Plan, trans3Lock ), _T( "fftRepo.getPlan failed" ) );

					trans3Plan->placeness     = fftPlan->allOpsInplace ? CLFFT_INPLACE : CLFFT_OUTOFPLACE;
					trans3Plan->precision     = fftPlan->precision;
					trans3Plan->tmpBufSize    = 0;
					trans3Plan->batchsize     = fftPlan->batchsize;
					trans3Plan->envelope	  = fftPlan->envelope;
					trans3Plan->inputLayout   = fftPlan->allOpsInplace ? fftPlan->inputLayout : CLFFT_COMPLEX_INTERLEAVED;
					trans3Plan->outputLayout  = fftPlan->outputLayout;
					trans3Plan->inStride[0]   = 1;
					trans3Plan->inStride[1]   = clLengths[0] + padding;
					trans3Plan->outStride[0]  = fftPlan->outStride[0];
					trans3Plan->outStride[1]  = clLengths[1];
					trans3Plan->iDist         = clLengths[1] * trans3Plan->inStride[1];
					trans3Plan->oDist         = fftPlan->oDist;
                    trans3Plan->gen           = transGen;
					trans3Plan->transflag     = true;
					trans3Plan->transOutHorizontal = true;


					if (trans3Plan->gen == Transpose_NONSQUARE)// inplace transpose
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							//trans3Plan->length.push_back(fftPlan->length[index]);
							/*
							replacing the line above with the two lines below since:
							fftPlan is still 1D, thus the broken down transpose should be 2D not 3D
							the batchSize for the transpose should increase accordingly.
							the iDist should decrease accordingly. Push back to length will cause a 3D transpose
							*/
							trans3Plan->batchsize = trans3Plan->batchsize * fftPlan->length[index];
							//trans3Plan->iDist = trans3Plan->iDist / fftPlan->length[index];
							//trans3Plan->inStride.push_back(trans3Plan->iDist);
							trans3Plan->inStride.push_back(fftPlan->inStride[index]);
							//trans3Plan->iDist *= fftPlan->length[index];
							trans3Plan->outStride.push_back(fftPlan->outStride[index]);
						}
					}
					else if (trans3Plan->gen == Transpose_SQUARE)
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							trans3Plan->batchsize = trans3Plan->batchsize * fftPlan->length[index];
							//trans3Plan->iDist = trans3Plan->iDist / fftPlan->length[index];
							//trans3Plan->inStride.push_back(trans3Plan->iDist);
							trans3Plan->inStride.push_back(fftPlan->inStride[index]);
							//trans3Plan->iDist *= fftPlan->length[index];
							trans3Plan->outStride.push_back(fftPlan->outStride[index]);
						}
					}
					else
					{
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							trans3Plan->length.push_back(fftPlan->length[index]);

							trans3Plan->inStride.push_back(trans3Plan->iDist);
							trans3Plan->iDist *= fftPlan->length[index];
							trans3Plan->outStride.push_back(fftPlan->outStride[index]);
						}
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						trans3Plan->hasPostCallback = true;
						trans3Plan->postCallbackParam = fftPlan->postCallbackParam;
						trans3Plan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTZ, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans3 plan failed" ) );

					fftPlan->transflag = true;
					fftPlan->baked = true;
					return	CLFFT_SUCCESS;
				}

				size_t length0 = clLengths[0];
				size_t length1 = clLengths[1];


				// For real transforms
				// Special case optimization with 5-step algorithm
				if( (fftPlan->inputLayout == CLFFT_REAL) && IsPo2(fftPlan->length[0])
					&& (fftPlan->length.size() == 1)
					&& (fftPlan->inStride[0] == 1) && (fftPlan->outStride[0] == 1)
					&& (fftPlan->length[0] > 4096) && (fftPlan->length.size() == 1) )
				{

					ARG_CHECK(clLengths[0] <= Large1DThreshold);


					size_t biggerDim = clLengths[0] > clLengths[1] ? clLengths[0] : clLengths[1];
					size_t smallerDim = biggerDim == clLengths[0] ? clLengths[1] : clLengths[0];
					size_t padding = 0;
					if( (smallerDim % 64 == 0) || (biggerDim % 64 == 0) )
						padding = 64;


					if (fftPlan->tmpBufSize==0 )
					{
						size_t Nf = (1 + smallerDim/2) * biggerDim;
						fftPlan->tmpBufSize = (smallerDim + padding) * biggerDim / 2;

						if(fftPlan->tmpBufSize < Nf) 
							fftPlan->tmpBufSize = Nf;

						fftPlan->tmpBufSize *= ( fftPlan->batchsize * fftPlan->ElementSize() );

						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSize *= fftPlan->length[index];
						}
					}

					if (fftPlan->tmpBufSizeRC==0 )
					{
						fftPlan->tmpBufSizeRC = fftPlan->tmpBufSize;
					}

					//Transpose
					//Input --> tmp buffer
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan Large1d transpose 1 failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->placeness     = CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->inputLayout   = fftPlan->inputLayout;
					trans1Plan->outputLayout  = CLFFT_REAL;
					trans1Plan->inStride[0]   = fftPlan->inStride[0];
					trans1Plan->inStride[1]   = clLengths[0];
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = clLengths[1] + padding;
					trans1Plan->iDist         = fftPlan->iDist;
					trans1Plan->oDist         = clLengths[0] * trans1Plan->outStride[1];
					trans1Plan->gen           = Transpose_GCN;
					trans1Plan->transflag     = true;

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						trans1Plan->hasPreCallback = true;
						trans1Plan->preCallback = fftPlan->preCallback;
						trans1Plan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans1 plan failed" ) );

					//Row transform
					//tmp->output
					//size clLengths[1], batch clLengths[0], with length[0] twiddle factor multiplication
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[1] ),
						_T( "CreateDefaultPlan Large1d column failed" ) );

					FFTPlan* row1Plan	= NULL;
					lockRAII* row1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, row1Plan, row1Lock ), _T( "fftRepo.getPlan failed" ) );

					row1Plan->placeness     = CLFFT_OUTOFPLACE;
					row1Plan->precision     = fftPlan->precision;
					row1Plan->forwardScale  = 1.0f;
					row1Plan->backwardScale = 1.0f;
					row1Plan->tmpBufSize    = 0;
					row1Plan->batchsize     = fftPlan->batchsize;

					row1Plan->gen			= fftPlan->gen;
					row1Plan->envelope		= fftPlan->envelope;

					// twiddling is done in row2
					row1Plan->large1D		= 0;

					row1Plan->length.push_back(clLengths[0]);
					row1Plan->inputLayout   = CLFFT_REAL;
					row1Plan->outputLayout  = CLFFT_HERMITIAN_INTERLEAVED;
					row1Plan->inStride[0]   = 1;
					row1Plan->outStride[0]  = 1;
					row1Plan->inStride.push_back(clLengths[1]+padding);
					row1Plan->outStride.push_back(1 + clLengths[1]/2);
					row1Plan->iDist         = clLengths[0] * row1Plan->inStride[1];
					row1Plan->oDist         = clLengths[0] * row1Plan->outStride[1]; 


					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d first row plan failed" ) );

					//Transpose 2
					//Output --> tmp buffer
					clLengths[2] = clLengths[0];
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, &clLengths[1] ),
						_T( "CreateDefaultPlan Large1d transpose 2 failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTY, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->transflag = true;

					size_t transLengths[2];
					transLengths[0] = 1 + clLengths[1]/2;
					transLengths[1] = clLengths[0];
					OPENCL_V(clfftSetPlanLength( fftPlan->planTY, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTY transpose failed" ) );



					trans2Plan->placeness     = CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					trans2Plan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
					trans2Plan->inStride[0]   = 1;
					trans2Plan->inStride[1]   = 1 + clLengths[1]/2;
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = clLengths[0];
					trans2Plan->iDist         = clLengths[0] * trans2Plan->inStride[1];
					trans2Plan->oDist         = (1 + clLengths[1]/2) * trans2Plan->outStride[1];
                    trans2Plan->gen           = Transpose_GCN;
					trans2Plan->transflag     = true;
					trans2Plan->transOutHorizontal = true;

					OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans2 plan failed" ) );

					//Row transform 2
					//tmp->tmp
					//size clLengths[0], batch clLengths[1]
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &clLengths[0] ),
						_T( "CreateDefaultPlan Large1d second row plan failed" ) );

					FFTPlan* row2Plan	= NULL;
					lockRAII* row2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, row2Plan, row2Lock ), _T( "fftRepo.getPlan failed" ) );

					row2Plan->placeness     = CLFFT_OUTOFPLACE;
					row2Plan->precision     = fftPlan->precision;
					row2Plan->forwardScale  = fftPlan->forwardScale;
					row2Plan->backwardScale = fftPlan->backwardScale;
					row2Plan->tmpBufSize    = 0;
					row2Plan->batchsize     = fftPlan->batchsize;

					row2Plan->gen			= fftPlan->gen;
					row2Plan->envelope		= fftPlan->envelope;


					row2Plan->length.push_back(1+clLengths[1]/2);
					row2Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					row2Plan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
					row2Plan->inStride[0]   = 1;
					row2Plan->outStride[0]  = 1;
					row2Plan->inStride.push_back(clLengths[0]);
					row2Plan->outStride.push_back(1 + clLengths[0]/2);
					row2Plan->iDist         = (1 + clLengths[1]/2) * row2Plan->inStride[1];
					row2Plan->oDist         = clLengths[1] * row2Plan->outStride[1];

					row2Plan->large1D		= fftPlan->length[0];
					row2Plan->twiddleFront  = true;

					row2Plan->realSpecial = true;
					row2Plan->realSpecial_Nr = clLengths[1];

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d second row plan failed" ) );

					//Transpose 3
					//tmp --> output
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTZ, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan Large1d transpose 3 failed" ) );

					FFTPlan* trans3Plan	= NULL;
					lockRAII* trans3Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTZ, trans3Plan, trans3Lock ), _T( "fftRepo.getPlan failed" ) );

					trans3Plan->transflag = true;

					transLengths[0] = 1 + clLengths[0]/2;
					transLengths[1] = clLengths[1];
					OPENCL_V(clfftSetPlanLength( fftPlan->planTZ, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTZ transpose failed" ) );

					trans3Plan->placeness     = CLFFT_OUTOFPLACE;
					trans3Plan->precision     = fftPlan->precision;
					trans3Plan->tmpBufSize    = 0;
					trans3Plan->batchsize     = fftPlan->batchsize;
					trans3Plan->envelope	  = fftPlan->envelope;
					trans3Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					if(fftPlan->outputLayout == CLFFT_HERMITIAN_PLANAR)
						trans3Plan->outputLayout  = CLFFT_COMPLEX_PLANAR;
					else
						trans3Plan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
					trans3Plan->inStride[0]   = 1;
					trans3Plan->inStride[1]   = 1 + clLengths[0]/2;
					trans3Plan->outStride[0]  = 1;
					trans3Plan->outStride[1]  = clLengths[1];
					trans3Plan->iDist         = clLengths[1] * trans3Plan->inStride[1];
					trans3Plan->oDist         = fftPlan->oDist;
                    trans3Plan->gen           = Transpose_GCN;
					trans3Plan->transflag     = true;
					trans3Plan->realSpecial	  = true;
					trans3Plan->transOutHorizontal = true;

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						trans3Plan->hasPostCallback = true;
						trans3Plan->postCallbackParam = fftPlan->postCallbackParam;
						trans3Plan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTZ, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan large1d trans3 plan failed" ) );

					fftPlan->transflag = true;
					fftPlan->baked = true;
					return	CLFFT_SUCCESS;
				}
				else if (fftPlan->inputLayout == CLFFT_REAL)
				{
					if (fftPlan->tmpBufSizeRC == 0)
					{
						fftPlan->tmpBufSizeRC = length0 * length1 *
							fftPlan->batchsize * fftPlan->ElementSize();
						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSizeRC *= fftPlan->length[index];
						}
					}

					// column FFT, size clLengths[1], batch clLengths[0], with length[0] twiddle factor multiplication
					// transposed output
					OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[1]),
						_T("CreateDefaultPlan Large1d column failed"));

					FFTPlan* colTPlan = NULL;
					lockRAII* colLock = NULL;
					OPENCL_V(fftRepo.getPlan(fftPlan->planX, colTPlan, colLock), _T("fftRepo.getPlan failed"));

					// current plan is to create intermediate buffer, packed and interleave
					// This is a column FFT, the first elements distance between each FFT is the distance of the first two
					// elements in the original buffer. Like a transpose of the matrix
					// we need to pass clLengths[0] and instride size to kernel, so kernel can tell the difference

					//this part are common for both passes
					colTPlan->placeness = CLFFT_OUTOFPLACE;
					colTPlan->precision = fftPlan->precision;
					colTPlan->forwardScale = 1.0f;
					colTPlan->backwardScale = 1.0f;
					colTPlan->tmpBufSize = 0;
					colTPlan->batchsize = fftPlan->batchsize;

					colTPlan->gen = fftPlan->gen;
					colTPlan->envelope = fftPlan->envelope;

					//Pass large1D flag to confirm we need multiply twiddle factor
					colTPlan->large1D = fftPlan->length[0];
					colTPlan->RCsimple = true;

					colTPlan->length.push_back(clLengths[0]);

					// first Pass
					colTPlan->inputLayout = fftPlan->inputLayout;
					colTPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
					colTPlan->inStride[0] = fftPlan->inStride[0] * clLengths[0];
					colTPlan->outStride[0] = 1;
					colTPlan->iDist = fftPlan->iDist;
					colTPlan->oDist = length0 * length1;//fftPlan->length[0];
					colTPlan->inStride.push_back(fftPlan->inStride[0]);
					colTPlan->outStride.push_back(length1);//clLengths[1]);

					for (size_t index = 1; index < fftPlan->length.size(); index++)
					{
						colTPlan->length.push_back(fftPlan->length[index]);
						colTPlan->inStride.push_back(fftPlan->inStride[index]);
						// tmp buffer is tightly packed
						colTPlan->outStride.push_back(colTPlan->oDist);
						colTPlan->oDist *= fftPlan->length[index];
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						colTPlan->hasPreCallback = true;
						colTPlan->preCallback = fftPlan->preCallback;
						colTPlan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL), _T("BakePlan large1d first column plan failed"));

					//another column FFT, size clLengths[0], batch clLengths[1], output without transpose
					OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planY, fftPlan->context, CLFFT_1D, &clLengths[0]),
						_T("CreateDefaultPlan large1D row failed"));

					FFTPlan* col2Plan = NULL;
					lockRAII* rowLock = NULL;
					OPENCL_V(fftRepo.getPlan(fftPlan->planY, col2Plan, rowLock), _T("fftRepo.getPlan failed"));

					// This is second column fft, intermediate buffer is packed and interleaved
					// we need to pass clLengths[1] and instride size to kernel, so kernel can tell the difference

					col2Plan->precision = fftPlan->precision;
					col2Plan->forwardScale = fftPlan->forwardScale;
					col2Plan->backwardScale = fftPlan->backwardScale;
					col2Plan->tmpBufSize = 0;
					col2Plan->batchsize = fftPlan->batchsize;

					col2Plan->gen = fftPlan->gen;
					col2Plan->envelope = fftPlan->envelope;

					col2Plan->length.push_back(length1);

					col2Plan->inStride[0] = length1;
					col2Plan->inStride.push_back(1);
					col2Plan->iDist = length0 * length1;

					// make sure colTPlan (first column plan) does not recurse, otherwise large twiddle mul
					// cannot be done with this algorithm sequence
					assert(colTPlan->planX == 0);


					col2Plan->placeness = CLFFT_INPLACE;
					col2Plan->inputLayout = CLFFT_COMPLEX_INTERLEAVED;
					col2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;

					col2Plan->outStride[0] = length1;
					col2Plan->outStride.push_back(1);
					col2Plan->oDist = length0 * length1;

					for (size_t index = 1; index < fftPlan->length.size(); index++)
					{
						col2Plan->length.push_back(fftPlan->length[index]);
						col2Plan->inStride.push_back(col2Plan->iDist);
						col2Plan->outStride.push_back(col2Plan->oDist);
						col2Plan->iDist *= fftPlan->length[index];
						col2Plan->oDist *= fftPlan->length[index];
					}


					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d second column plan failed" ) );

					if ( (fftPlan->outputLayout == CLFFT_HERMITIAN_INTERLEAVED) ||
						 (fftPlan->outputLayout == CLFFT_HERMITIAN_PLANAR) )
					{
						// copy plan to get back to hermitian
						OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planRCcopy, fftPlan->context, CLFFT_1D, &fftPlan->length[0]),
							_T("CreateDefaultPlan RC copy failed"));

						FFTPlan* copyPlan = NULL;
						lockRAII* copyLock = NULL;
						OPENCL_V(fftRepo.getPlan(fftPlan->planRCcopy, copyPlan, copyLock), _T("fftRepo.getPlan failed"));

						// This is second column fft, intermediate buffer is packed and interleaved
						// we need to pass clLengths[1] and instride size to kernel, so kernel can tell the difference

						// common part for both passes
						copyPlan->placeness = CLFFT_OUTOFPLACE;
						copyPlan->inputLayout = CLFFT_COMPLEX_INTERLEAVED;
						copyPlan->outputLayout = fftPlan->outputLayout;

						copyPlan->precision = fftPlan->precision;
						copyPlan->forwardScale = 1.0f;
						copyPlan->backwardScale = 1.0f;
						copyPlan->tmpBufSize = 0;
						copyPlan->batchsize = fftPlan->batchsize;

						copyPlan->gen = Copy;
						copyPlan->envelope = fftPlan->envelope;


						copyPlan->inStride[0] = 1;
						copyPlan->iDist = fftPlan->length[0];

						copyPlan->outStride[0] = fftPlan->outStride[0];
						copyPlan->oDist = fftPlan->oDist;

						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							copyPlan->length.push_back(fftPlan->length[index]);
							copyPlan->inStride.push_back(copyPlan->inStride[index - 1] * fftPlan->length[index - 1]);
							copyPlan->iDist *= fftPlan->length[index];
							copyPlan->outStride.push_back(fftPlan->outStride[index]);
						}

						//Set callback data if set on top level plan
						if (fftPlan->hasPostCallback)
						{
							copyPlan->hasPostCallback = true;
							copyPlan->postCallbackParam = fftPlan->postCallbackParam;
							copyPlan->postcallUserData = fftPlan->postcallUserData;
						}

						OPENCL_V(clfftBakePlan(fftPlan->planRCcopy, numQueues, commQueueFFT, NULL, NULL), _T("BakePlan large1d RC copy plan failed"));
					}

				}
				else if(fftPlan->outputLayout == CLFFT_REAL)
				{
					if (fftPlan->tmpBufSizeRC==0 )
					{
						fftPlan->tmpBufSizeRC = length0 * length1 *
							fftPlan->batchsize * fftPlan->ElementSize();
						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSizeRC *= fftPlan->length[index];
						}
					}

					if ((fftPlan->inputLayout == CLFFT_HERMITIAN_INTERLEAVED) ||
						(fftPlan->inputLayout == CLFFT_HERMITIAN_PLANAR))
					{
						// copy plan to from hermitian to full complex
						OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planRCcopy, fftPlan->context, CLFFT_1D, &fftPlan->length[0]),
							_T("CreateDefaultPlan RC copy failed"));

						FFTPlan* copyPlan = NULL;
						lockRAII* copyLock = NULL;
						OPENCL_V(fftRepo.getPlan(fftPlan->planRCcopy, copyPlan, copyLock), _T("fftRepo.getPlan failed"));

						// This is second column fft, intermediate buffer is packed and interleaved
						// we need to pass clLengths[1] and instride size to kernel, so kernel can tell the difference

						// common part for both passes
						copyPlan->placeness = CLFFT_OUTOFPLACE;
						copyPlan->inputLayout = fftPlan->inputLayout;
						copyPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;

						copyPlan->precision = fftPlan->precision;
						copyPlan->forwardScale = 1.0f;
						copyPlan->backwardScale = 1.0f;
						copyPlan->tmpBufSize = 0;
						copyPlan->batchsize = fftPlan->batchsize;

						copyPlan->gen = Copy;
						copyPlan->envelope = fftPlan->envelope;

						copyPlan->inStride[0] = fftPlan->inStride[0];
						copyPlan->iDist = fftPlan->iDist;

						copyPlan->outStride[0] = 1;
						copyPlan->oDist = fftPlan->length[0];

						for (size_t index = 1; index < fftPlan->length.size(); index++)
						{
							copyPlan->length.push_back(fftPlan->length[index]);
							copyPlan->outStride.push_back(copyPlan->outStride[index - 1] * fftPlan->length[index - 1]);
							copyPlan->oDist *= fftPlan->length[index];
							copyPlan->inStride.push_back(fftPlan->inStride[index]);
						}

						//Set callback data if set on top level plan
						if (fftPlan->hasPreCallback)
						{
							copyPlan->hasPreCallback = true;
							copyPlan->preCallback = fftPlan->preCallback;
							copyPlan->precallUserData = fftPlan->precallUserData;
						}

						OPENCL_V(clfftBakePlan(fftPlan->planRCcopy, numQueues, commQueueFFT, NULL, NULL), _T("BakePlan large1d RC copy plan failed"));
					}

					// column FFT, size clLengths[1], batch clLengths[0], with length[0] twiddle factor multiplication
					// transposed output
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[1] ),
						_T( "CreateDefaultPlan Large1d column failed" ) );

					FFTPlan* colTPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, colTPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					// current plan is to create intermediate buffer, packed and interleave
					// This is a column FFT, the first elements distance between each FFT is the distance of the first two
					// elements in the original buffer. Like a transpose of the matrix
					// we need to pass clLengths[0] and instride size to kernel, so kernel can tell the difference

					//this part are common for both passes
					colTPlan->precision     = fftPlan->precision;
					colTPlan->forwardScale  = 1.0f;
					colTPlan->backwardScale = 1.0f;
					colTPlan->tmpBufSize    = 0;
					colTPlan->batchsize     = fftPlan->batchsize;

					colTPlan->gen			= fftPlan->gen;
					colTPlan->envelope		= fftPlan->envelope;

					//Pass large1D flag to confirm we need multiply twiddle factor
					colTPlan->large1D       = fftPlan->length[0];

					colTPlan->length.push_back(clLengths[0]);

					colTPlan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					colTPlan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
					
					colTPlan->inStride[0]  = length0;
					colTPlan->inStride.push_back(1);
					colTPlan->iDist        = length0 * length1;

					colTPlan->outStride[0] = length0;
					colTPlan->outStride.push_back(1);
					colTPlan->oDist         = length0 * length1;

					for (size_t index=1; index < fftPlan->length.size(); index++)
					{
						colTPlan->length.push_back(fftPlan->length[index]);
						colTPlan->inStride.push_back(colTPlan->iDist);
						colTPlan->outStride.push_back(colTPlan->oDist);
						colTPlan->iDist   *= fftPlan->length[index];
						colTPlan->oDist   *= fftPlan->length[index];
					}

					if ((fftPlan->inputLayout == CLFFT_HERMITIAN_INTERLEAVED) ||
						(fftPlan->inputLayout == CLFFT_HERMITIAN_PLANAR))
					{
						colTPlan->placeness = CLFFT_INPLACE;
					}
					else
					{
						colTPlan->placeness = CLFFT_OUTOFPLACE;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d first column plan failed" ) );

					//another column FFT, size clLengths[0], batch clLengths[1], output without transpose
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D,  &clLengths[0] ),
						_T( "CreateDefaultPlan large1D row failed" ) );

					FFTPlan* col2Plan	= NULL;
					lockRAII* rowLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, col2Plan, rowLock ), _T( "fftRepo.getPlan failed" ) );

					// This is second column fft, intermediate buffer is packed and interleaved
					// we need to pass clLengths[1] and instride size to kernel, so kernel can tell the difference

					// common part for both passes
					col2Plan->placeness     = CLFFT_OUTOFPLACE;
					col2Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					col2Plan->outputLayout  = fftPlan->outputLayout;

					col2Plan->precision     = fftPlan->precision;
					col2Plan->forwardScale  = fftPlan->forwardScale;
					col2Plan->backwardScale = fftPlan->backwardScale;
					col2Plan->tmpBufSize    = 0;
					col2Plan->batchsize     = fftPlan->batchsize;

					col2Plan->gen			= fftPlan->gen;
					col2Plan->envelope			= fftPlan->envelope;

					col2Plan->RCsimple = true;
					col2Plan->length.push_back(length1);

					col2Plan->inStride[0]  = 1;
					col2Plan->inStride.push_back(length0);
					col2Plan->iDist        = length0 * length1;

					col2Plan->outStride[0] = length1 * fftPlan->outStride[0];
					col2Plan->outStride.push_back(fftPlan->outStride[0]);
					col2Plan->oDist         = fftPlan->oDist;

					for (size_t index=1; index < fftPlan->length.size(); index++)
					{
						col2Plan->length.push_back(fftPlan->length[index]);
						col2Plan->inStride.push_back(col2Plan->iDist);
						col2Plan->iDist   *= fftPlan->length[index];
						col2Plan->outStride.push_back(fftPlan->outStride[index]);
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						col2Plan->hasPostCallback = true;
						col2Plan->postCallbackParam = fftPlan->postCallbackParam;
						col2Plan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d second column plan failed" ) );
				}
				else
				{

					if( (fftPlan->length[0] > 262144/PrecisionWidth(fftPlan->precision)) && fftPlan->blockCompute )
					{
						assert(fftPlan->length[0] <= 1048576);


						size_t padding = 64;	
						if (fftPlan->tmpBufSize==0 )
						{
							fftPlan->tmpBufSize = (length1 + padding) * length0 *
									fftPlan->batchsize * fftPlan->ElementSize();
							for (size_t index=1; index < fftPlan->length.size(); index++)
							{
								fftPlan->tmpBufSize *= fftPlan->length[index];
							}
						}

						// Algorithm in this case is 
						// T(with pad, out_of_place), R (in_place), C(in_place), Unpad(out_of_place)

						size_t len[3] = { clLengths[1], clLengths[0], 1 };

						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, len ),
						_T( "CreateDefaultPlan Large1d trans1 failed" ) );

						FFTPlan* trans1Plan	= NULL;
						lockRAII* trans1Lock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

						trans1Plan->placeness     = CLFFT_OUTOFPLACE;
						trans1Plan->precision     = fftPlan->precision;
						trans1Plan->tmpBufSize    = 0;
						trans1Plan->batchsize     = fftPlan->batchsize;
						trans1Plan->envelope	  = fftPlan->envelope;
						trans1Plan->inputLayout   = fftPlan->inputLayout;
						trans1Plan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						trans1Plan->inStride[0]   = fftPlan->inStride[0];
						trans1Plan->inStride[1]   = length1;
						trans1Plan->outStride[0]  = 1;
						trans1Plan->outStride[1]  = length0 + padding;
						trans1Plan->iDist         = fftPlan->iDist;
						trans1Plan->oDist         = length1 * trans1Plan->outStride[1];
						trans1Plan->gen           = Transpose_GCN;
						trans1Plan->transflag     = true;

						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							trans1Plan->length.push_back(fftPlan->length[index]);
							trans1Plan->inStride.push_back(fftPlan->inStride[index]);
							trans1Plan->outStride.push_back(trans1Plan->oDist);
							trans1Plan->oDist *= fftPlan->length[index];
						}

						//Set callback data if set on top level plan
						if (fftPlan->hasPreCallback)
						{
							trans1Plan->hasPreCallback = true;
							trans1Plan->preCallback = fftPlan->preCallback;
							trans1Plan->precallUserData = fftPlan->precallUserData;
						}

						OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
							_T( "BakePlan large1d trans1 plan failed" ) );


						// row FFT
						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[0] ),
							_T( "CreateDefaultPlan Large1d column failed" ) );

						FFTPlan* rowPlan	= NULL;
						lockRAII* rowLock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

						assert(fftPlan->large1D == 0);

						rowPlan->placeness     = CLFFT_INPLACE;
						rowPlan->precision     = fftPlan->precision;
						rowPlan->forwardScale  = 1.0f;
						rowPlan->backwardScale = 1.0f;
						rowPlan->tmpBufSize    = 0;
						rowPlan->batchsize     = fftPlan->batchsize;

						rowPlan->gen			= fftPlan->gen;
						rowPlan->envelope		= fftPlan->envelope;

						rowPlan->length.push_back(length1);


						rowPlan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
						rowPlan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						rowPlan->inStride[0]   = 1;
						rowPlan->outStride[0]  = 1;
						rowPlan->inStride.push_back(length0+padding);
						rowPlan->outStride.push_back(length0+padding);
						rowPlan->iDist         = (length0+padding)*length1;
						rowPlan->oDist         = (length0+padding)*length1;

						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							rowPlan->length.push_back(fftPlan->length[index]);
							rowPlan->inStride.push_back(rowPlan->iDist);
							rowPlan->iDist *= fftPlan->length[index];
							rowPlan->outStride.push_back(rowPlan->oDist);
							rowPlan->oDist *= fftPlan->length[index];
						}


						OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d first row plan failed" ) );

						//column FFT
						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D,  &clLengths[1] ),
							_T( "CreateDefaultPlan large1D column failed" ) );

						FFTPlan* col2Plan	= NULL;
						lockRAII* colLock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planY, col2Plan, colLock ), _T( "fftRepo.getPlan failed" ) );

						col2Plan->placeness     = CLFFT_INPLACE;
						col2Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
						col2Plan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						col2Plan->precision     = fftPlan->precision;
						col2Plan->forwardScale  = fftPlan->forwardScale;
						col2Plan->backwardScale = fftPlan->backwardScale;
						col2Plan->tmpBufSize    = 0;
						col2Plan->batchsize     = fftPlan->batchsize;

						col2Plan->gen			= fftPlan->gen;
						col2Plan->envelope		= fftPlan->envelope;

						col2Plan->large1D       = fftPlan->length[0];
						col2Plan->twiddleFront	= true;

						col2Plan->length.push_back(clLengths[0]);



						col2Plan->blockCompute = true;
						col2Plan->blockComputeType = BCT_C2C;

						col2Plan->inStride[0]  = length0+padding;
						col2Plan->outStride[0] = length0+padding;
						col2Plan->iDist        = (length0+padding) * length1;
						col2Plan->oDist        = (length0+padding) * length1;
						col2Plan->inStride.push_back(1);
						col2Plan->outStride.push_back(1);


						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							col2Plan->length.push_back(fftPlan->length[index]);
							col2Plan->inStride.push_back(col2Plan->iDist);
							col2Plan->outStride.push_back(col2Plan->oDist);
							col2Plan->iDist   *= fftPlan->length[index];
							col2Plan->oDist   *= fftPlan->length[index];
						}


						OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d second column plan failed" ) );


						// copy plan to get results back to packed output
						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planCopy, fftPlan->context, CLFFT_1D,  &clLengths[0] ),
							_T( "CreateDefaultPlan Copy failed" ) );

						FFTPlan* copyPlan	= NULL;
						lockRAII* copyLock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planCopy, copyPlan, copyLock ), _T( "fftRepo.getPlan failed" ) );


						copyPlan->placeness     = CLFFT_OUTOFPLACE;
						copyPlan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
						copyPlan->outputLayout  = fftPlan->outputLayout;

						copyPlan->precision     = fftPlan->precision;
						copyPlan->forwardScale  = 1.0f;
						copyPlan->backwardScale = 1.0f;
						copyPlan->tmpBufSize    = 0;
						copyPlan->batchsize     = fftPlan->batchsize;

						copyPlan->gen			= Copy;
						copyPlan->envelope		= fftPlan->envelope;

						copyPlan->length.push_back(length1);

						copyPlan->inStride[0]  = 1;
						copyPlan->inStride.push_back(length0+padding);
						copyPlan->iDist        = length1*(length0+padding);

						copyPlan->outStride[0] = fftPlan->outStride[0];
						copyPlan->outStride.push_back(length0);
						copyPlan->oDist         = fftPlan->oDist;

						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							copyPlan->length.push_back(fftPlan->length[index]);
							copyPlan->inStride.push_back(copyPlan->inStride[index] * copyPlan->length[index]);
							copyPlan->iDist   *= fftPlan->length[index];
							copyPlan->outStride.push_back(fftPlan->outStride[index]);
						}

						OPENCL_V(clfftBakePlan(fftPlan->planCopy, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d copy plan failed" ) );
					}
					else
					{

						if (fftPlan->tmpBufSize==0 )
						{
							fftPlan->tmpBufSize = length0 * length1 *
								fftPlan->batchsize * fftPlan->ElementSize();
							for (size_t index=1; index < fftPlan->length.size(); index++)
							{
								fftPlan->tmpBufSize *= fftPlan->length[index];
							}
						}

						// column FFT, size clLengths[1], batch clLengths[0], with length[0] twiddle factor multiplication
						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &clLengths[1] ),
							_T( "CreateDefaultPlan Large1d column failed" ) );

						FFTPlan* colTPlan	= NULL;
						lockRAII* colLock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planX, colTPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

						assert(fftPlan->large1D == 0);

						// current plan is to create intermediate buffer, packed and interleave
						// This is a column FFT, the first elements distance between each FFT is the distance of the first two
						// elements in the original buffer. Like a transpose of the matrix
						// we need to pass clLengths[0] and instride size to kernel, so kernel can tell the difference

						//this part are common for both passes
						colTPlan->placeness     = CLFFT_OUTOFPLACE;
						colTPlan->precision     = fftPlan->precision;
						colTPlan->forwardScale  = 1.0f;
						colTPlan->backwardScale = 1.0f;
						colTPlan->tmpBufSize    = 0;
						colTPlan->batchsize     = fftPlan->batchsize;

						colTPlan->gen			= fftPlan->gen;
						colTPlan->envelope			= fftPlan->envelope;

						//Pass large1D flag to confirm we need multiply twiddle factor
						colTPlan->large1D       = fftPlan->length[0];

						colTPlan->length.push_back(length0);


						colTPlan->inputLayout   = fftPlan->inputLayout;
						colTPlan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						colTPlan->inStride[0]   = fftPlan->inStride[0] * length0;
						colTPlan->outStride[0]  = length0;
						colTPlan->iDist         = fftPlan->iDist;
						colTPlan->oDist         = length0 * length1;
						colTPlan->inStride.push_back(fftPlan->inStride[0]);
						colTPlan->outStride.push_back(1);

						//Set callback data if set on top level plan
						if (fftPlan->hasPreCallback)
						{
							colTPlan->hasPreCallback = true;
							colTPlan->preCallback = fftPlan->preCallback;
							colTPlan->precallUserData = fftPlan->precallUserData;
						}

						// Enabling block column compute
						if( (colTPlan->inStride[0] == length0) && IsPo2(fftPlan->length[0]) && (fftPlan->length[0] < 524288) )
						{
							colTPlan->blockCompute = true;
							colTPlan->blockComputeType = BCT_C2C;
						}

						for (size_t index=1; index < fftPlan->length.size(); index++)
						{
							colTPlan->length.push_back(fftPlan->length[index]);
							colTPlan->inStride.push_back(fftPlan->inStride[index]);
							// tmp buffer is tightly packed
							colTPlan->outStride.push_back(colTPlan->oDist);
							colTPlan->oDist        *= fftPlan->length[index];
						}


						OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d first column plan failed" ) );

						//another column FFT, size clLengths[0], batch clLengths[1], output with transpose
						OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D,  &clLengths[0] ),
							_T( "CreateDefaultPlan large1D row failed" ) );

						FFTPlan* col2Plan	= NULL;
						lockRAII* rowLock	= NULL;
						OPENCL_V( fftRepo.getPlan( fftPlan->planY, col2Plan, rowLock ), _T( "fftRepo.getPlan failed" ) );

						// This is second column fft, intermediate buffer is packed and interleaved
						// we need to pass clLengths[1] and instride size to kernel, so kernel can tell the difference

						// common part for both passes
						col2Plan->outputLayout  = fftPlan->outputLayout;
						col2Plan->precision     = fftPlan->precision;
						col2Plan->forwardScale  = fftPlan->forwardScale;
						col2Plan->backwardScale = fftPlan->backwardScale;
						col2Plan->tmpBufSize    = 0;
						col2Plan->batchsize     = fftPlan->batchsize;
						col2Plan->oDist         = fftPlan->oDist;

						col2Plan->gen			= fftPlan->gen;
						col2Plan->envelope		= fftPlan->envelope;


						col2Plan->length.push_back(clLengths[1]);

						bool integratedTranposes = true;


						if( colTPlan->blockCompute && (fftPlan->outStride[0] == 1) && clLengths[0] <= 256)
						{
							col2Plan->blockCompute = true;
							col2Plan->blockComputeType = BCT_R2C;

							col2Plan->placeness    = CLFFT_OUTOFPLACE;
							col2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
							col2Plan->inStride[0]  = 1;
							col2Plan->outStride[0] = length1;
							col2Plan->iDist        = length0 * length1;
							col2Plan->inStride.push_back(length0);
							col2Plan->outStride.push_back(1);
						}
						else if( colTPlan->blockCompute && (fftPlan->outStride[0] == 1) )
						{
							integratedTranposes = false;

							col2Plan->placeness    = CLFFT_INPLACE;
							col2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
							col2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							col2Plan->inStride[0]  = 1;
							col2Plan->outStride[0] = 1;
							col2Plan->iDist        = length0 * length1;
							col2Plan->oDist        = length0 * length1;
							col2Plan->inStride.push_back(length0);
							col2Plan->outStride.push_back(length0);
						}
						else
						{
							//first layer, large 1D from tmp buffer to output buffer
							col2Plan->placeness    = CLFFT_OUTOFPLACE;
							col2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
							col2Plan->inStride[0]  = 1;
							col2Plan->outStride[0] = fftPlan->outStride[0] * clLengths[1];
							col2Plan->iDist        = length0 * length1; //fftPlan->length[0];
							col2Plan->inStride.push_back(length0);
							col2Plan->outStride.push_back(fftPlan->outStride[0]);
						}

						if(!integratedTranposes)
						{
							for (size_t index=1; index < fftPlan->length.size(); index++)
							{
								col2Plan->length.push_back(fftPlan->length[index]);
								col2Plan->inStride.push_back(col2Plan->iDist);
								col2Plan->outStride.push_back(col2Plan->oDist);
								col2Plan->iDist        *= fftPlan->length[index];
								col2Plan->oDist        *= fftPlan->length[index];
							}
						}
						else
						{
							for (size_t index=1; index < fftPlan->length.size(); index++)
							{
								col2Plan->length.push_back(fftPlan->length[index]);
								col2Plan->inStride.push_back(col2Plan->iDist);
								col2Plan->outStride.push_back(fftPlan->outStride[index]);
								col2Plan->iDist   *= fftPlan->length[index];
							}
						}

						//Set callback data if set on top level plan
						if (fftPlan->hasPostCallback && integratedTranposes)
						{
							col2Plan->hasPostCallback = true;
							col2Plan->postCallbackParam = fftPlan->postCallbackParam;
							col2Plan->postcallUserData = fftPlan->postcallUserData;
						}

						OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan large1d second column plan failed" ) );

						if(!integratedTranposes)
						{
							//Transpose 
							//tmp --> output
							OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTZ, fftPlan->context, CLFFT_2D, clLengths ),
								_T( "CreateDefaultPlan Large1d transpose failed" ) );

							FFTPlan* trans3Plan	= NULL;
							lockRAII* trans3Lock	= NULL;
							OPENCL_V( fftRepo.getPlan( fftPlan->planTZ, trans3Plan, trans3Lock ), _T( "fftRepo.getPlan failed" ) );

							trans3Plan->placeness     = CLFFT_OUTOFPLACE;
							trans3Plan->precision     = fftPlan->precision;
							trans3Plan->tmpBufSize    = 0;
							trans3Plan->batchsize     = fftPlan->batchsize;
							trans3Plan->envelope	  = fftPlan->envelope;
							trans3Plan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
							trans3Plan->outputLayout  = fftPlan->outputLayout;
							trans3Plan->inStride[0]   = 1;
							trans3Plan->inStride[1]   = clLengths[0];
							trans3Plan->outStride[0]  = fftPlan->outStride[0];
							trans3Plan->outStride[1]  = clLengths[1] * fftPlan->outStride[0];
							trans3Plan->iDist         = fftPlan->length[0];
							trans3Plan->oDist         = fftPlan->oDist;
							trans3Plan->gen           = Transpose_GCN;
							trans3Plan->transflag     = true;

							for (size_t index=1; index < fftPlan->length.size(); index++)
							{
								trans3Plan->length.push_back(fftPlan->length[index]);
								trans3Plan->inStride.push_back(trans3Plan->iDist);
								trans3Plan->iDist *= fftPlan->length[index];
								trans3Plan->outStride.push_back(fftPlan->outStride[index]);
							}

							//Set callback data if set on top level plan
							if (fftPlan->hasPostCallback)
							{
								trans3Plan->hasPostCallback = true;
								trans3Plan->postCallbackParam = fftPlan->postCallbackParam;
								trans3Plan->postcallUserData = fftPlan->postcallUserData;
							}

							OPENCL_V(clfftBakePlan(fftPlan->planTZ, numQueues, commQueueFFT, NULL, NULL ),
								_T( "BakePlan large1d trans plan failed" ) );
						}
					}
				}

				fftPlan->baked = true;
				return	CLFFT_SUCCESS;
			}
		}
		break;
	case CLFFT_2D:
		{

			if (fftPlan->transflag) //Transpose for 2D
			{
                clfftStatus err = CLFFT_SUCCESS;
				if(fftPlan->gen == Transpose_GCN)
					fftPlan->action = new FFTGeneratedTransposeGCNAction(plHandle, fftPlan, *commQueueFFT, err);
				else if (fftPlan->gen == Transpose_SQUARE)
					fftPlan->action = new FFTGeneratedTransposeSquareAction(plHandle, fftPlan, *commQueueFFT, err);
                else if (fftPlan->gen == Transpose_NONSQUARE)
                {
					if(fftPlan->nonSquareKernelType != NON_SQUARE_TRANS_PARENT)
						fftPlan->action = new FFTGeneratedTransposeNonSquareAction(plHandle, fftPlan, *commQueueFFT, err);
					else
					{
						size_t clLengths[] = { 1, 1, 0 };
						clLengths[0] = fftPlan->length[0];
						clLengths[1] = fftPlan->length[1];

						//NON_SQUARE_KERNEL_ORDER currKernelOrder;
						// controling the transpose and swap kernel order
						// if leading dim is larger than the other dim it makes sense to swap and transpose
						if (clLengths[0] > clLengths[1])
						{
							//Twiddling will be done in swap kernel, in regardless of the order
							fftPlan->nonSquareKernelOrder = SWAP_AND_TRANSPOSE;
						}
						else
						{
							if (fftPlan->large1D != 0 && 0)
							{
                                //this is not going to happen anymore
								fftPlan->nonSquareKernelOrder = TRANSPOSE_LEADING_AND_SWAP;
							}
							else
							{
                                //twiddling can be done in swap
								fftPlan->nonSquareKernelOrder = TRANSPOSE_AND_SWAP;
							}
						}

						//std::cout << "currKernelOrder = " << fftPlan->nonSquareKernelOrder << std::endl;
						//ends tranpose kernel order

						//Transpose stage 1 
						OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planTX, fftPlan->context, CLFFT_2D, clLengths),
							_T("CreateDefaultPlan transpose_nsq_stage1 plan failed"));

						FFTPlan* trans1Plan = NULL;
						lockRAII* trans1Lock = NULL;
						OPENCL_V(fftRepo.getPlan(fftPlan->planTX, trans1Plan, trans1Lock), _T("fftRepo.getPlan failed"));

						trans1Plan->placeness = CLFFT_INPLACE;
						trans1Plan->precision = fftPlan->precision;
						trans1Plan->tmpBufSize = 0;
						trans1Plan->batchsize = fftPlan->batchsize;
						trans1Plan->envelope = fftPlan->envelope;
						trans1Plan->inputLayout = fftPlan->inputLayout;
						trans1Plan->outputLayout = fftPlan->outputLayout;
						trans1Plan->inStride[0] = fftPlan->inStride[0];
						trans1Plan->outStride[0] = fftPlan->outStride[0];
						trans1Plan->inStride[1] = fftPlan->inStride[1];
						trans1Plan->outStride[1] = fftPlan->outStride[1];
						trans1Plan->iDist = fftPlan->iDist;
						trans1Plan->oDist = fftPlan->oDist;
						trans1Plan->gen = Transpose_NONSQUARE;
						trans1Plan->nonSquareKernelOrder = fftPlan->nonSquareKernelOrder;
						if(fftPlan->nonSquareKernelOrder == SWAP_AND_TRANSPOSE)
							trans1Plan->nonSquareKernelType = NON_SQUARE_TRANS_SWAP;
						else if (fftPlan->nonSquareKernelOrder == TRANSPOSE_AND_SWAP)
							trans1Plan->nonSquareKernelType = NON_SQUARE_TRANS_TRANSPOSE_BATCHED;
						else if(fftPlan->nonSquareKernelOrder == TRANSPOSE_LEADING_AND_SWAP)
							trans1Plan->nonSquareKernelType = NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING;
						trans1Plan->transflag = true;
                        trans1Plan->large1D = fftPlan->large1D;//twiddling may happen in this kernel

						if (trans1Plan->nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
						{
							//this should be in a function to avoide duplicate code TODO
							//need to treat a non square matrix as a sqaure matrix with bigger batch size
							size_t lengthX = trans1Plan->length[0];
							size_t lengthY = trans1Plan->length[1];

							size_t BatchFactor = (lengthX > lengthY) ? (lengthX / lengthY) : (lengthY / lengthX);
							trans1Plan->transposeMiniBatchSize = BatchFactor;
							trans1Plan->batchsize *= BatchFactor;
							trans1Plan->iDist = trans1Plan->iDist / BatchFactor;
							if (lengthX > lengthY)
							{
								trans1Plan->length[0] = lengthX / BatchFactor;
								trans1Plan->inStride[1] = lengthX / BatchFactor;
							}
							else if (lengthX < lengthY)
							{
								trans1Plan->length[1] = lengthY / BatchFactor;
								trans1Plan->inStride[1] = lengthX;
							}
						}

						for (size_t index = 2; index < fftPlan->length.size(); index++)
						{
							trans1Plan->length.push_back(fftPlan->length[index]);
							trans1Plan->inStride.push_back(fftPlan->inStride[index]);
							trans1Plan->outStride.push_back(fftPlan->outStride[index]);
						}

						if (fftPlan->hasPreCallback)
						{
							trans1Plan->hasPreCallback = true;
							trans1Plan->preCallback = fftPlan->preCallback;
							trans1Plan->precallUserData = fftPlan->precallUserData;
						}


						OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL),
							_T("BakePlan transpose_nsq_stage1 plan failed"));


						//Transpose stage 2 
						OPENCL_V(clfftCreateDefaultPlanInternal(&fftPlan->planTY, fftPlan->context, CLFFT_2D, clLengths),
							_T("CreateDefaultPlan transpose_nsq_stage2 plan failed"));

						FFTPlan* trans2Plan = NULL;
						lockRAII* trans2Lock = NULL;
						OPENCL_V(fftRepo.getPlan(fftPlan->planTY, trans2Plan, trans2Lock), _T("fftRepo.getPlan failed"));

						trans2Plan->placeness = CLFFT_INPLACE;
						trans2Plan->precision = fftPlan->precision;
						trans2Plan->tmpBufSize = 0;
						trans2Plan->batchsize = fftPlan->batchsize;
						trans2Plan->envelope = fftPlan->envelope;
						trans2Plan->inputLayout = fftPlan->inputLayout;
						trans2Plan->outputLayout = fftPlan->outputLayout;
						trans2Plan->inStride[0] = fftPlan->inStride[0];
						trans2Plan->outStride[0] = fftPlan->outStride[0];
						trans2Plan->inStride[1] = fftPlan->inStride[1];
						trans2Plan->outStride[1] = fftPlan->outStride[1];
						trans2Plan->iDist = fftPlan->iDist;
						trans2Plan->oDist = fftPlan->oDist;
						trans2Plan->gen = Transpose_NONSQUARE;
						trans2Plan->nonSquareKernelOrder = fftPlan->nonSquareKernelOrder;
						if (fftPlan->nonSquareKernelOrder == SWAP_AND_TRANSPOSE)
							trans2Plan->nonSquareKernelType = NON_SQUARE_TRANS_TRANSPOSE_BATCHED;
						else if(fftPlan->nonSquareKernelOrder == TRANSPOSE_AND_SWAP)
							trans2Plan->nonSquareKernelType = NON_SQUARE_TRANS_SWAP;
						else if(fftPlan->nonSquareKernelOrder == TRANSPOSE_LEADING_AND_SWAP)
							trans2Plan->nonSquareKernelType = NON_SQUARE_TRANS_SWAP;
						trans2Plan->transflag = true;
						trans2Plan->large1D = fftPlan->large1D;//twiddling may happen in this kernel

						if (trans2Plan->nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
						{
							//need to treat a non square matrix as a sqaure matrix with bigger batch size
							size_t lengthX = trans2Plan->length[0];
							size_t lengthY = trans2Plan->length[1];

							size_t BatchFactor = (lengthX > lengthY) ? (lengthX/lengthY) : (lengthY/lengthX);
							trans2Plan->transposeMiniBatchSize = BatchFactor;
							trans2Plan->batchsize *= BatchFactor;
							trans2Plan->iDist = trans2Plan->iDist / BatchFactor;
							if (lengthX > lengthY)
							{
								trans2Plan->length[0] = lengthX / BatchFactor;
								trans2Plan->inStride[1] = lengthX / BatchFactor;
							}
							else if(lengthX < lengthY)
							{
								trans2Plan->length[1] = lengthY / BatchFactor;
								trans2Plan->inStride[1] = lengthX;
							}
						}

						for (size_t index = 2; index < fftPlan->length.size(); index++)
						{
							trans2Plan->length.push_back(fftPlan->length[index]);
							trans2Plan->inStride.push_back(fftPlan->inStride[index]);
							trans2Plan->outStride.push_back(fftPlan->outStride[index]);
						}

						if (fftPlan->hasPostCallback)
						{
							trans2Plan->hasPostCallback = true;
							trans2Plan->postCallbackParam = fftPlan->postCallbackParam;
							trans2Plan->postcallUserData = fftPlan->postcallUserData;
						}

						OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL),
							_T("BakePlan transpose_nsq_stage2 plan failed"));
					}
                }
				else
					fftPlan->action = new FFTGeneratedTransposeGCNAction(plHandle, fftPlan, *commQueueFFT, err);

                OPENCL_V( err, "FFTGeneratedTransposeXXXAction failed");

				fftPlan->baked		= true;
				return	CLFFT_SUCCESS;
			}

			size_t length0 = fftPlan->length[0];
			size_t length1 = fftPlan->length[1];


			if (fftPlan->length[0] > Large1DThreshold ||
				fftPlan->length[1] > Large1DThreshold)
				fftPlan->large2D = true;

			while (1 && (fftPlan->inputLayout != CLFFT_REAL) && (fftPlan->outputLayout != CLFFT_REAL))
			{
				//break;


                // TODO : Check for a better way to do this.
                bool isnvidia = false;
                for (size_t Idx = 0; !isnvidia && Idx < numQueues; Idx++)
                {
                    cl_command_queue QIdx = commQueueFFT[Idx];
                    cl_device_id Device;
                    clGetCommandQueueInfo(QIdx, CL_QUEUE_DEVICE, sizeof(Device), &Device, NULL);
                    char Vendor[256];
                    clGetDeviceInfo(Device, CL_DEVICE_VENDOR, sizeof(Vendor), &Vendor, NULL);
                    isnvidia |= (strncmp(Vendor, "NVIDIA", 6) == 0);
                }
                // nvidia gpus are failing when doing transpose for 2D FFTs
                if (isnvidia) break;

				if (fftPlan->length.size() != 2) break;
				if (!(IsPo2(fftPlan->length[0])) || !(IsPo2(fftPlan->length[1])))
					break;
				if (fftPlan->length[1] < 32) break;
				//TBD: restrict the use large2D in x!=y case becase we will need two temp buffers
				//     (1) for 2D usage (2) for 1D large usage
				//if (fftPlan->large2D) break;
				//Performance show 512 is the good case with transpose
				//if user want the result to be transposed, then we will.

				if (fftPlan->length[0] < 64) break;
				//x!=y case, we need tmp buffer, currently temp buffer only support interleaved format
				//if (fftPlan->length[0] != fftPlan->length[1] && fftPlan->outputLayout == CLFFT_COMPLEX_PLANAR) break;
				if (fftPlan->inStride[0] != 1 || fftPlan->outStride[0] != 1 ||
					fftPlan->inStride[1] != fftPlan->length[0] || fftPlan->outStride[1] != fftPlan->length[0])
					break;
				//if (fftPlan->placeness != CLFFT_INPLACE || fftPlan->inputLayout != CLFFT_COMPLEX_PLANAR)
				//	break;
				//if (fftPlan->batchsize != 1) break;
				//if (fftPlan->precision != CLFFT_SINGLE) break;

				fftPlan->transflag = true;

				//create row plan,
				// x=y & x!=y, In->In for inplace, In->out for outofplace
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimX ] ),
					_T( "CreateDefaultPlan for planX failed" ) );

				FFTPlan* rowPlan	= NULL;
				lockRAII* rowLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

				rowPlan->inputLayout     = fftPlan->inputLayout;
				rowPlan->outputLayout    = fftPlan->outputLayout;
				rowPlan->placeness       = fftPlan->placeness;
				rowPlan->outStride[0]    = fftPlan->outStride[0];
				rowPlan->outStride.push_back(fftPlan->outStride[1]);
				rowPlan->oDist           = fftPlan->oDist;
				rowPlan->precision       = fftPlan->precision;
				rowPlan->forwardScale    = 1.0f;
				rowPlan->backwardScale   = 1.0f;
				rowPlan->tmpBufSize      = 0;

				rowPlan->gen			 = fftPlan->gen;
				rowPlan->envelope		 = fftPlan->envelope;
				rowPlan->batchsize       = fftPlan->batchsize;
				rowPlan->inStride[0]     = fftPlan->inStride[0];
				rowPlan->length.push_back(fftPlan->length[1]);
				rowPlan->inStride.push_back(fftPlan->inStride[1]);
				rowPlan->iDist           = fftPlan->iDist;
				
				//Set callback data if set on top level plan
				if (fftPlan->hasPreCallback)
				{
					rowPlan->hasPreCallback = true;
					rowPlan->preCallback = fftPlan->preCallback;
					rowPlan->precallUserData = fftPlan->precallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ),
					_T( "BakePlan for planX failed" ) );

				//Create transpose plan for first transpose
				//x=y: inplace. x!=y inplace: in->tmp, outofplace out->tmp
				size_t clLengths[] = { 1, 1, 0 };
				clLengths[0] = fftPlan->length[0];
				clLengths[1] = fftPlan->length[1];

				size_t biggerDim = clLengths[0] > clLengths[1] ? clLengths[0] : clLengths[1];
				size_t smallerDim = biggerDim == clLengths[0] ? clLengths[1] : clLengths[0];
				size_t padding = 0;

				fftPlan->transpose_in_2d_inplace = (clLengths[0]==clLengths[1]) ? true : false;
				if ( (!fftPlan->transpose_in_2d_inplace) && fftPlan->tmpBufSize==0 && fftPlan->length.size()<=2 )
				{
					if ((smallerDim % 64 == 0) || (biggerDim % 64 == 0))
						if(biggerDim > 512)
							padding = 64;

					// we need tmp buffer for x!=y case
					// we assume the tmp buffer is packed interleaved
					fftPlan->tmpBufSize = (smallerDim + padding) * biggerDim *
						fftPlan->batchsize * fftPlan->ElementSize();
				}

				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, clLengths ),
					_T( "CreateDefaultPlan for planT failed" ) );

				FFTPlan* transPlanX	= NULL;
				lockRAII* transLockX	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planTX, transPlanX, transLockX ), _T( "fftRepo.getPlan failed" ) );

				transPlanX->inputLayout     = fftPlan->outputLayout;
				transPlanX->precision       = fftPlan->precision;
				transPlanX->tmpBufSize      = 0;

				transPlanX->envelope		= fftPlan->envelope;
				transPlanX->batchsize       = fftPlan->batchsize;
				transPlanX->inStride[0]     = fftPlan->outStride[0];
				transPlanX->inStride[1]     = fftPlan->outStride[1];
				transPlanX->iDist           = fftPlan->oDist;
				transPlanX->transflag       = true;

				if (!fftPlan->transpose_in_2d_inplace)
				{
					transPlanX->gen = Transpose_GCN;
					transPlanX->outputLayout    = CLFFT_COMPLEX_INTERLEAVED;
					transPlanX->placeness       = CLFFT_OUTOFPLACE;
					transPlanX->outStride[0]    = 1;
					transPlanX->outStride[1]    = clLengths[1] + padding;
					transPlanX->oDist           = clLengths[0] * transPlanX->outStride[1];
				}
				else
				{
					transPlanX->gen = Transpose_SQUARE;
					transPlanX->outputLayout    = fftPlan->outputLayout;
					transPlanX->placeness       = CLFFT_INPLACE;
					transPlanX->outStride[0]    = fftPlan->outStride[0];
					transPlanX->outStride[1]    = fftPlan->outStride[1];
					transPlanX->oDist           = fftPlan->oDist;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
					_T( "BakePlan for planTX failed" ) );

				//create second row plan
				//x!=y: tmp->tmp, x=y case: In->In or Out->Out
				//if Transposed result is a choice x!=y: tmp->In or out
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
					_T( "CreateDefaultPlan for planY failed" ) );

				FFTPlan* colPlan	= NULL;
				lockRAII* colLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

				if (!fftPlan->transpose_in_2d_inplace)
				{
					colPlan->inputLayout     = CLFFT_COMPLEX_INTERLEAVED;
					colPlan->inStride[0]     = 1;
					colPlan->inStride.push_back(clLengths[1] + padding);
					colPlan->iDist           = clLengths[0] * colPlan->inStride[1];

					if (fftPlan->transposed == CLFFT_NOTRANSPOSE)
					{
						colPlan->outputLayout    = CLFFT_COMPLEX_INTERLEAVED;
						colPlan->outStride[0]    = 1;
						colPlan->outStride.push_back(clLengths[1] + padding);
						colPlan->oDist           = clLengths[0] * colPlan->outStride[1];
						colPlan->placeness       = CLFFT_INPLACE;
					}
					else
					{
						colPlan->outputLayout    = fftPlan->outputLayout;
						colPlan->outStride[0]    = fftPlan->outStride[0];
						colPlan->outStride.push_back(clLengths[1] * fftPlan->outStride[0]);
						colPlan->oDist           = fftPlan->oDist;
						colPlan->placeness       = CLFFT_OUTOFPLACE;
					}
				}
				else
				{
					colPlan->inputLayout     = fftPlan->outputLayout;
					colPlan->outputLayout    = fftPlan->outputLayout;
					colPlan->outStride[0]    = fftPlan->outStride[0];
					colPlan->outStride.push_back(fftPlan->outStride[1]);
					colPlan->oDist           = fftPlan->oDist;
					colPlan->inStride[0]     = fftPlan->outStride[0];
					colPlan->inStride.push_back(fftPlan->outStride[1]);
					colPlan->iDist           = fftPlan->oDist;
					colPlan->placeness       = CLFFT_INPLACE;
				}

				colPlan->precision       = fftPlan->precision;
				colPlan->forwardScale    = fftPlan->forwardScale;
				colPlan->backwardScale   = fftPlan->backwardScale;
				colPlan->tmpBufSize      = 0;

				colPlan->gen			 = fftPlan->gen;
				colPlan->envelope		 = fftPlan->envelope;
				colPlan->batchsize       = fftPlan->batchsize;
				colPlan->length.push_back(fftPlan->length[0]);

				OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ),
					_T( "BakePlan for planY failed" ) );

				if (fftPlan->transposed == CLFFT_TRANSPOSED)
				{
					fftPlan->baked = true;
					return	CLFFT_SUCCESS;
				}

				//Create transpose plan for second transpose
				//x!=y case tmp->In or Out, x=y case In->In or Out->out
				size_t clLengthsY[2] = { clLengths[1], clLengths[0] };
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, clLengthsY ),
					_T( "CreateDefaultPlan for planTY failed" ) );

				FFTPlan* transPlanY	= NULL;
				lockRAII* transLockY	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planTY, transPlanY, transLockY ), _T( "fftRepo.getPlan failed" ) );

				if (!fftPlan->transpose_in_2d_inplace)
				{
					transPlanY->gen = Transpose_GCN;
					transPlanY->inputLayout     = CLFFT_COMPLEX_INTERLEAVED;
					transPlanY->placeness       = CLFFT_OUTOFPLACE;
					transPlanY->inStride[0]     = 1;
					transPlanY->inStride[1]     = clLengths[1] + padding;
					transPlanY->iDist           = clLengths[0] * transPlanY->inStride[1];
					transPlanY->transOutHorizontal = true;
				}
				else
				{
					transPlanY->gen = Transpose_SQUARE;
					transPlanY->inputLayout     = fftPlan->outputLayout;
					transPlanY->placeness       = CLFFT_INPLACE;
					transPlanY->inStride[0]     = fftPlan->outStride[0];
					transPlanY->inStride[1]     = fftPlan->outStride[1];
					transPlanY->iDist           = fftPlan->oDist;
				}
				transPlanY->outputLayout    = fftPlan->outputLayout;
				transPlanY->outStride[0]    = fftPlan->outStride[0];
				transPlanY->outStride[1]    = fftPlan->outStride[1];
				transPlanY->oDist           = fftPlan->oDist;
				transPlanY->precision       = fftPlan->precision;
				transPlanY->tmpBufSize      = 0;

				transPlanY->envelope		= fftPlan->envelope;
				transPlanY->batchsize       = fftPlan->batchsize;
				transPlanY->transflag       = true;

				//Set callback data if set on top level plan
				if (fftPlan->hasPostCallback)
				{
					transPlanY->hasPostCallback = true;
					transPlanY->postCallbackParam = fftPlan->postCallbackParam;
					transPlanY->postcallUserData = fftPlan->postcallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
					_T( "BakePlan for planTY failed" ) );

				fftPlan->baked = true;
				return	CLFFT_SUCCESS;
			}

			//check transposed
			if (fftPlan->transposed != CLFFT_NOTRANSPOSE)
				return CLFFT_TRANSPOSED_NOTIMPLEMENTED;


			if(fftPlan->inputLayout == CLFFT_REAL)
			{
				length0 = fftPlan->length[0];
				length1 = fftPlan->length[1];

				size_t Nt = (1 + length0/2);


				// create row plan
				// real to hermitian

				//create row plan
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimX ] ),
					_T( "CreateDefaultPlan for planX failed" ) );

				FFTPlan* rowPlan	= NULL;
				lockRAII* rowLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );


				rowPlan->outputLayout  = fftPlan->outputLayout;
				rowPlan->inputLayout  = fftPlan->inputLayout;
				rowPlan->placeness     = fftPlan->placeness;
				rowPlan->length.push_back(length1);

				rowPlan->inStride[0]  = fftPlan->inStride[0];
				rowPlan->inStride.push_back(fftPlan->inStride[1]);
				rowPlan->iDist         = fftPlan->iDist;

				rowPlan->precision     = fftPlan->precision;
				rowPlan->forwardScale  = 1.0f;
				rowPlan->backwardScale = 1.0f;
				rowPlan->tmpBufSize    = 0;

				rowPlan->gen			= fftPlan->gen;
				rowPlan->envelope		= fftPlan->envelope;

				rowPlan->batchsize    = fftPlan->batchsize;

				rowPlan->outStride[0]  = fftPlan->outStride[0];
				rowPlan->outStride.push_back(fftPlan->outStride[1]);
				rowPlan->oDist         = fftPlan->oDist;

				//this 2d is decomposed from 3d
				for (size_t index=2; index < fftPlan->length.size(); index++)
				{
					rowPlan->length.push_back(fftPlan->length[index]);
					rowPlan->inStride.push_back(fftPlan->inStride[index]);
					rowPlan->outStride.push_back(fftPlan->outStride[index]);
				}

				//Set callback data if set on top level plan
				if (fftPlan->hasPreCallback)
				{
					rowPlan->hasPreCallback = true;
					rowPlan->preCallback = fftPlan->preCallback;
					rowPlan->precallUserData = fftPlan->precallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planX failed" ) );

				if( (rowPlan->inStride[0] == 1) && (rowPlan->outStride[0] == 1) &&
					( ((rowPlan->inStride[1] == Nt*2) && (rowPlan->placeness == CLFFT_INPLACE)) ||
					  ((rowPlan->inStride[1] == length0) && (rowPlan->placeness == CLFFT_OUTOFPLACE)) )
					&& (rowPlan->outStride[1] == Nt) )
				{
					// calc temp buf size
					if (fftPlan->tmpBufSize==0)
					{
						fftPlan->tmpBufSize = Nt * length1 * fftPlan->batchsize * fftPlan->ElementSize();

						for (size_t index=2; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSize *= fftPlan->length[index];
						}
					}

					// create first transpose plan
					
					//Transpose 
					// output --> tmp
					size_t transLengths[2] = { length0, length1 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, transLengths ),
						_T( "CreateDefaultPlan for planTX transpose failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->transflag = true;

					transLengths[0] = Nt;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTX, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTX transpose failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					trans1Plan->placeness     = CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->forwardScale  = 1.0f;
					trans1Plan->backwardScale = 1.0f;

					trans1Plan->inStride[0]   = 1;
					trans1Plan->inStride[1]   = Nt;
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = length1;
					trans1Plan->iDist         = rowPlan->oDist;
					trans1Plan->oDist		  = Nt*length1;
					trans1Plan->transOutHorizontal = true;

					trans1Plan->gen           = Transpose_GCN;


					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						trans1Plan->length.push_back(fftPlan->length[index]);
						trans1Plan->inStride.push_back(rowPlan->outStride[index]);
						trans1Plan->outStride.push_back(trans1Plan->oDist);
						trans1Plan->oDist *= fftPlan->length[index];
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTX failed" ) );


					// Create column plan as a row plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
						_T( "CreateDefaultPlan for planY failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					colPlan->outputLayout  = trans1Plan->outputLayout;
					colPlan->inputLayout   = trans1Plan->outputLayout;
					colPlan->placeness     = CLFFT_INPLACE;
					colPlan->length.push_back(Nt);

					colPlan->inStride[0]  = 1;
					colPlan->inStride.push_back(length1);
					colPlan->iDist         = Nt*length1;

					colPlan->outStride[0]  = 1;
					colPlan->outStride.push_back(length1);
					colPlan->oDist         = Nt*length1;

					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = fftPlan->forwardScale;
					colPlan->backwardScale = fftPlan->backwardScale;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope		= fftPlan->envelope;

					colPlan->batchsize    = fftPlan->batchsize;

					//this 2d is decomposed from 3d
					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->inStride.push_back(colPlan->iDist);
						colPlan->outStride.push_back(colPlan->oDist);
						colPlan->iDist *= fftPlan->length[index];
						colPlan->oDist *= fftPlan->length[index];
					}

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planY failed" ) );

					if (fftPlan->transposed == CLFFT_TRANSPOSED)
					{
						fftPlan->baked = true;
						return	CLFFT_SUCCESS;
					}

					// create second transpose plan
					
					//Transpose 
					//output --> tmp
					size_t trans2Lengths[2] = { length1, length0 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, trans2Lengths ),
						_T( "CreateDefaultPlan for planTY transpose failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTY, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->transflag = true;

					trans2Lengths[1] = Nt;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTY, CLFFT_2D, trans2Lengths ),
						_T( "clfftSetPlanLength for planTY transpose failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans2Plan->outputLayout = CLFFT_COMPLEX_PLANAR;
							trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					default: assert(false);
					}

					trans2Plan->placeness     = CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->forwardScale  = 1.0f;
					trans2Plan->backwardScale = 1.0f;

					trans2Plan->inStride[0]   = 1;
					trans2Plan->inStride[1]   = length1;
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = Nt;
					trans2Plan->iDist         = Nt*length1;
					trans2Plan->oDist		  = fftPlan->oDist;

					trans2Plan->gen           = Transpose_GCN;
					trans2Plan->transflag     = true;

					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						trans2Plan->length.push_back(fftPlan->length[index]);
						trans2Plan->inStride.push_back(trans2Plan->iDist);
						trans2Plan->iDist *= fftPlan->length[index];
						trans2Plan->outStride.push_back(fftPlan->outStride[index]);

					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						trans2Plan->hasPostCallback = true;
						trans2Plan->postCallbackParam = fftPlan->postCallbackParam;
						trans2Plan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTY failed" ) );

				}
				else
				{
					// create col plan
					// complex to complex

					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
						_T( "CreateDefaultPlan for planY failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_PLANAR;
							colPlan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					colPlan->placeness     = CLFFT_INPLACE;
					colPlan->length.push_back(Nt);

					colPlan->outStride[0]  = fftPlan->outStride[1];
					colPlan->outStride.push_back(fftPlan->outStride[0]);
					colPlan->oDist         = fftPlan->oDist;


					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = fftPlan->forwardScale;
					colPlan->backwardScale = fftPlan->backwardScale;
					colPlan->tmpBufSize    = fftPlan->tmpBufSize;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope			= fftPlan->envelope;

					colPlan->batchsize = fftPlan->batchsize;

					colPlan->inStride[0]  = rowPlan->outStride[1];
					colPlan->inStride.push_back(rowPlan->outStride[0]);
					colPlan->iDist         = rowPlan->oDist;

					//this 2d is decomposed from 3d
					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->outStride.push_back(fftPlan->outStride[index]);
						colPlan->inStride.push_back(rowPlan->outStride[index]);
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						colPlan->hasPostCallback = true;
						colPlan->postCallbackParam = fftPlan->postCallbackParam;
						colPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planY failed" ) );
				}

			}
			else if(fftPlan->outputLayout == CLFFT_REAL)
			{
				length0 = fftPlan->length[0];
				length1 = fftPlan->length[1];

				size_t Nt = (1 + length0/2);
				if (fftPlan->tmpBufSize==0)
				{
					fftPlan->tmpBufSize = Nt * length1 * fftPlan->batchsize * fftPlan->ElementSize();
					for (size_t index=2; index < fftPlan->length.size(); index++)
						fftPlan->tmpBufSize *= fftPlan->length[index];
				}

				if ((fftPlan->tmpBufSizeC2R==0) && (fftPlan->placeness == CLFFT_OUTOFPLACE) && (fftPlan->length.size() == 2))
				{
					fftPlan->tmpBufSizeC2R = fftPlan->tmpBufSize;
				}

				if( (fftPlan->inStride[0] == 1) && (fftPlan->outStride[0] == 1) &&
					( ((fftPlan->outStride[1] == Nt*2) && (fftPlan->oDist == Nt*2*length1) && (fftPlan->placeness == CLFFT_INPLACE)) ||
						((fftPlan->outStride[1] == length0) && (fftPlan->oDist == length0*length1) && (fftPlan->placeness == CLFFT_OUTOFPLACE)) )
					&& (fftPlan->inStride[1] == Nt) && (fftPlan->iDist == Nt*length1) )
				{
					// create first transpose plan
					
					//Transpose 
					// input --> tmp
					size_t transLengths[2] = { length0, length1 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, transLengths ),
						_T( "CreateDefaultPlan for planTY transpose failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTY, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->transflag = true;

					transLengths[0] = Nt;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTY, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTY transpose failed" ) );

					switch(fftPlan->inputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					trans1Plan->placeness     = CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->forwardScale  = 1.0f;
					trans1Plan->backwardScale = 1.0f;

					trans1Plan->inStride[0]   = 1;
					trans1Plan->inStride[1]   = Nt;
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = length1;
					trans1Plan->iDist         = fftPlan->iDist;
					trans1Plan->oDist		  = Nt*length1;
					trans1Plan->transOutHorizontal = true;

					trans1Plan->gen           = Transpose_GCN;


					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						trans1Plan->length.push_back(fftPlan->length[index]);
						trans1Plan->inStride.push_back(fftPlan->inStride[index]);
						trans1Plan->outStride.push_back(trans1Plan->oDist);
						trans1Plan->oDist *= fftPlan->length[index];
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						trans1Plan->hasPreCallback = true;
						trans1Plan->preCallback = fftPlan->preCallback;
						trans1Plan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTY failed" ) );

					// create col plan
					// complex to complex

					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
						_T( "CreateDefaultPlan for planY failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					colPlan->length.push_back(Nt);

					colPlan->inStride[0]  = 1;
					colPlan->inStride.push_back(length1);
					colPlan->iDist         = trans1Plan->oDist;

					colPlan->placeness = CLFFT_INPLACE;
					colPlan->inputLayout = CLFFT_COMPLEX_INTERLEAVED;
					colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;

					colPlan->outStride[0]  = colPlan->inStride[0];
					colPlan->outStride.push_back(colPlan->inStride[1]);
					colPlan->oDist         = colPlan->iDist;

					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->inStride.push_back(trans1Plan->outStride[index]);
						colPlan->outStride.push_back(trans1Plan->outStride[index]);
					}


					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = 1.0f;
					colPlan->backwardScale = 1.0f;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope		= fftPlan->envelope;

					colPlan->batchsize = fftPlan->batchsize;

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planY failed" ) );

					// create second transpose plan
					
					//Transpose 
					//tmp --> output
					size_t trans2Lengths[2] = { length1, length0 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, trans2Lengths ),
						_T( "CreateDefaultPlan for planTX transpose failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->transflag = true;

					trans2Lengths[1] = Nt;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTX, CLFFT_2D, trans2Lengths ),
						_T( "clfftSetPlanLength for planTX transpose failed" ) );


					trans2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
					trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;


					trans2Plan->placeness     = CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->forwardScale  = 1.0f;
					trans2Plan->backwardScale = 1.0f;

					trans2Plan->inStride[0]   = 1;
					trans2Plan->inStride[1]   = length1;
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = Nt;
					trans2Plan->iDist         = colPlan->oDist;
					trans2Plan->oDist		  = Nt*length1;

					trans2Plan->gen           = Transpose_GCN;
					trans2Plan->transflag     = true;

					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						trans2Plan->length.push_back(fftPlan->length[index]);
						trans2Plan->inStride.push_back(colPlan->outStride[index]);
						trans2Plan->outStride.push_back(trans2Plan->oDist);
						trans2Plan->oDist *= fftPlan->length[index];

					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTX failed" ) );

					// create row plan
					// hermitian to real

					//create row plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimX ] ),
						_T( "CreateDefaultPlan for planX failed" ) );

					FFTPlan* rowPlan	= NULL;
					lockRAII* rowLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

					rowPlan->outputLayout  = fftPlan->outputLayout;
					rowPlan->inputLayout   = CLFFT_HERMITIAN_INTERLEAVED;

					rowPlan->length.push_back(length1);

					rowPlan->outStride[0]  = fftPlan->outStride[0];
					rowPlan->outStride.push_back(fftPlan->outStride[1]);
					rowPlan->oDist         = fftPlan->oDist;

					rowPlan->inStride[0]  = trans2Plan->outStride[0];
					rowPlan->inStride.push_back(trans2Plan->outStride[1]);
					rowPlan->iDist         = trans2Plan->oDist;

					for (size_t index=2; index < fftPlan->length.size(); index++)
					{
						rowPlan->length.push_back(fftPlan->length[index]);
						rowPlan->inStride.push_back(trans2Plan->outStride[index]);
						rowPlan->outStride.push_back(fftPlan->outStride[index]);
					}

					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						rowPlan->placeness     = CLFFT_INPLACE;
					}
					else
					{
						rowPlan->placeness     = CLFFT_OUTOFPLACE;
					}				


					rowPlan->precision     = fftPlan->precision;
					rowPlan->forwardScale  = fftPlan->forwardScale;
					rowPlan->backwardScale = fftPlan->backwardScale;
					rowPlan->tmpBufSize    = 0;

					rowPlan->gen			= fftPlan->gen;
					rowPlan->envelope		= fftPlan->envelope;

					rowPlan->batchsize    = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						rowPlan->hasPostCallback = true;
						rowPlan->postCallbackParam = fftPlan->postCallbackParam;
						rowPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planX failed" ) );
				}
				else
				{

					// create col plan
					// complex to complex

					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
						_T( "CreateDefaultPlan for planY failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );


					switch(fftPlan->inputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}


					colPlan->length.push_back(Nt);

					colPlan->inStride[0]  = fftPlan->inStride[1];
					colPlan->inStride.push_back(fftPlan->inStride[0]);
					colPlan->iDist         = fftPlan->iDist;


					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						colPlan->placeness = CLFFT_INPLACE;
					}
					else
					{
						if(fftPlan->length.size() > 2)
							colPlan->placeness = CLFFT_INPLACE;
						else
							colPlan->placeness = CLFFT_OUTOFPLACE;
					}

					if(colPlan->placeness == CLFFT_INPLACE)
					{
						colPlan->outStride[0]  = colPlan->inStride[0];
						colPlan->outStride.push_back(colPlan->inStride[1]);
						colPlan->oDist         = colPlan->iDist;

						for (size_t index=2; index < fftPlan->length.size(); index++)
						{
							colPlan->length.push_back(fftPlan->length[index]);
							colPlan->inStride.push_back(fftPlan->inStride[index]);
							colPlan->outStride.push_back(fftPlan->inStride[index]);
						}
					}
					else
					{
						colPlan->outStride[0]  = Nt;
						colPlan->outStride.push_back(1);
						colPlan->oDist         = Nt*length1;

						for (size_t index=2; index < fftPlan->length.size(); index++)
						{
							colPlan->length.push_back(fftPlan->length[index]);
							colPlan->inStride.push_back(fftPlan->inStride[index]);
							colPlan->outStride.push_back(colPlan->oDist);
							colPlan->oDist *= fftPlan->length[index];
						}
					}

					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = 1.0f;
					colPlan->backwardScale = 1.0f;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope			= fftPlan->envelope;

					colPlan->batchsize = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						colPlan->hasPreCallback = true;
						colPlan->preCallback = fftPlan->preCallback;
						colPlan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planY failed" ) );

					// create row plan
					// hermitian to real

					//create row plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimX ] ),
						_T( "CreateDefaultPlan for planX failed" ) );

					FFTPlan* rowPlan	= NULL;
					lockRAII* rowLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

					rowPlan->outputLayout  = fftPlan->outputLayout;
					rowPlan->inputLayout   = CLFFT_HERMITIAN_INTERLEAVED;

					rowPlan->length.push_back(length1);

					rowPlan->outStride[0]  = fftPlan->outStride[0];
					rowPlan->outStride.push_back(fftPlan->outStride[1]);
					rowPlan->oDist         = fftPlan->oDist;

					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						rowPlan->placeness     = CLFFT_INPLACE;

						rowPlan->inStride[0]  = colPlan->outStride[1];
						rowPlan->inStride.push_back(colPlan->outStride[0]);
						rowPlan->iDist         = colPlan->oDist;

						for (size_t index=2; index < fftPlan->length.size(); index++)
						{
							rowPlan->length.push_back(fftPlan->length[index]);
							rowPlan->inStride.push_back(colPlan->outStride[index]);
							rowPlan->outStride.push_back(fftPlan->outStride[index]);
						}
					}
					else
					{
						rowPlan->placeness     = CLFFT_OUTOFPLACE;

						rowPlan->inStride[0]   = 1;
						rowPlan->inStride.push_back(Nt);
						rowPlan->iDist         = Nt*length1;

						for (size_t index=2; index < fftPlan->length.size(); index++)
						{
							rowPlan->length.push_back(fftPlan->length[index]);
							rowPlan->outStride.push_back(fftPlan->outStride[index]);
							rowPlan->inStride.push_back(rowPlan->iDist);						
							rowPlan->iDist *= fftPlan->length[index];
						}
					}
				

					rowPlan->precision     = fftPlan->precision;
					rowPlan->forwardScale  = fftPlan->forwardScale;
					rowPlan->backwardScale = fftPlan->backwardScale;
					rowPlan->tmpBufSize    = 0;

					rowPlan->gen			= fftPlan->gen;
					rowPlan->envelope		= fftPlan->envelope;

					rowPlan->batchsize    = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						rowPlan->hasPostCallback = true;
						rowPlan->postCallbackParam = fftPlan->postCallbackParam;
						rowPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planX failed" ) );
				}
			}
			else
			{
				if (fftPlan->tmpBufSize==0 && fftPlan->length.size()<=2)
				{
					fftPlan->tmpBufSize = length0 * length1 *
						fftPlan->batchsize * fftPlan->ElementSize();
				}

				//create row plan
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimX ] ),
					_T( "CreateDefaultPlan for planX failed" ) );

				FFTPlan* rowPlan	= NULL;
				lockRAII* rowLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

				rowPlan->inputLayout   = fftPlan->inputLayout;
				if (fftPlan->large2D || fftPlan->length.size()>2)
				{
					rowPlan->outputLayout  = fftPlan->outputLayout;
					rowPlan->placeness     = fftPlan->placeness;
					rowPlan->outStride[0]  = fftPlan->outStride[0];
					rowPlan->outStride.push_back(fftPlan->outStride[1]);
					rowPlan->oDist         = fftPlan->oDist;
				}
				else
				{
					rowPlan->outputLayout  = CLFFT_COMPLEX_INTERLEAVED;
					rowPlan->placeness     = CLFFT_OUTOFPLACE;
					rowPlan->outStride[0]  = length1;//1;
					rowPlan->outStride.push_back(1);//length0);
					rowPlan->oDist         = length0 * length1;
				}
				rowPlan->precision     = fftPlan->precision;
				rowPlan->forwardScale  = 1.0f;
				rowPlan->backwardScale = 1.0f;
				rowPlan->tmpBufSize    = fftPlan->tmpBufSize;

				rowPlan->gen			= fftPlan->gen;
				rowPlan->envelope			= fftPlan->envelope;

				// This is the row fft, the first elements distance between the first two FFTs is the distance of the first elements
				// of the first two rows in the original buffer.
				rowPlan->batchsize    = fftPlan->batchsize;
				rowPlan->inStride[0]  = fftPlan->inStride[0];

				//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
				rowPlan->length.push_back(fftPlan->length[1]);
				rowPlan->inStride.push_back(fftPlan->inStride[1]);

				//this 2d is decomposed from 3d
				if (fftPlan->length.size()>2)
				{
					rowPlan->length.push_back(fftPlan->length[2]);
					rowPlan->inStride.push_back(fftPlan->inStride[2]);
					rowPlan->outStride.push_back(fftPlan->outStride[2]);
				}

				rowPlan->iDist    = fftPlan->iDist;

				//Set callback data if set on top level plan
				if (fftPlan->hasPreCallback)
				{
					rowPlan->hasPreCallback = true;
					rowPlan->preCallback = fftPlan->preCallback;
					rowPlan->precallUserData = fftPlan->precallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planX failed" ) );

				//create col plan
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planY, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimY ] ),
					_T( "CreateDefaultPlan for planY failed" ) );

				FFTPlan* colPlan	= NULL;
				lockRAII* colLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planY, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

				if (fftPlan->large2D || fftPlan->length.size()>2)
				{
					colPlan->inputLayout   = fftPlan->outputLayout;
					colPlan->placeness     = CLFFT_INPLACE;
					colPlan->inStride[0]   = fftPlan->outStride[1];
					colPlan->inStride.push_back(fftPlan->outStride[0]);
					colPlan->iDist         = fftPlan->oDist;
				}
				else
				{
					colPlan->inputLayout   = CLFFT_COMPLEX_INTERLEAVED;
					colPlan->placeness     = CLFFT_OUTOFPLACE;
					colPlan->inStride[0]   = 1;//length0;
					colPlan->inStride.push_back(length1);//1);
					colPlan->iDist         = length0 * length1;
				}

				colPlan->outputLayout  = fftPlan->outputLayout;
				colPlan->precision     = fftPlan->precision;
				colPlan->forwardScale  = fftPlan->forwardScale;
				colPlan->backwardScale = fftPlan->backwardScale;
				colPlan->tmpBufSize    = fftPlan->tmpBufSize;

				colPlan->gen			= fftPlan->gen;
				colPlan->envelope			= fftPlan->envelope;

				// This is a column FFT, the first elements distance between each FFT is the distance of the first two
				// elements in the original buffer. Like a transpose of the matrix
				colPlan->batchsize = fftPlan->batchsize;
				colPlan->outStride[0] = fftPlan->outStride[1];

				//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
				colPlan->length.push_back(fftPlan->length[0]);
				colPlan->outStride.push_back(fftPlan->outStride[0]);
				colPlan->oDist    = fftPlan->oDist;

				//this 2d is decomposed from 3d
				if (fftPlan->length.size()>2)
				{
					//assert(fftPlan->large2D);
					colPlan->length.push_back(fftPlan->length[2]);
					colPlan->inStride.push_back(fftPlan->outStride[2]);
					colPlan->outStride.push_back(fftPlan->outStride[2]);
				}

				//Set callback data if set on top level plan
				if (fftPlan->hasPostCallback)
				{
					colPlan->hasPostCallback = true;
					colPlan->postCallbackParam = fftPlan->postCallbackParam;
					colPlan->postcallUserData = fftPlan->postcallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planY, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planY failed" ) );
			}

			fftPlan->baked = true;
			return	CLFFT_SUCCESS;
		}
	case CLFFT_3D:
		{
			if(fftPlan->inputLayout == CLFFT_REAL)
			{

				size_t length0 = fftPlan->length[ DimX ];
				size_t length1 = fftPlan->length[ DimY ];
				size_t length2 = fftPlan->length[ DimZ ];

				size_t Nt = (1 + length0/2);


				//create 2D xy plan
				size_t clLengths[] = { length0, length1, 0 };
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_2D, clLengths ),
					_T( "CreateDefaultPlan 2D planX failed" ) );

				FFTPlan* xyPlan	= NULL;
				lockRAII* rowLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planX, xyPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

				xyPlan->inputLayout   = fftPlan->inputLayout;
				xyPlan->outputLayout  = fftPlan->outputLayout;
				xyPlan->placeness     = fftPlan->placeness;
				xyPlan->precision     = fftPlan->precision;
				xyPlan->forwardScale  = 1.0f;
				xyPlan->backwardScale = 1.0f;
				xyPlan->tmpBufSize    = fftPlan->tmpBufSize;

				xyPlan->gen			 = fftPlan->gen;
				xyPlan->envelope			 = fftPlan->envelope;

				// This is the xy fft, the first elements distance between the first two FFTs is the distance of the first elements
				// of the first two rows in the original buffer.
				xyPlan->batchsize    = fftPlan->batchsize;
				xyPlan->inStride[0]  = fftPlan->inStride[0];
				xyPlan->inStride[1]  = fftPlan->inStride[1];
				xyPlan->outStride[0] = fftPlan->outStride[0];
				xyPlan->outStride[1] = fftPlan->outStride[1];

				//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
				xyPlan->length.push_back(fftPlan->length[2]);
				xyPlan->inStride.push_back(fftPlan->inStride[2]);
				xyPlan->outStride.push_back(fftPlan->outStride[2]);
				xyPlan->iDist    = fftPlan->iDist;
				xyPlan->oDist    = fftPlan->oDist;

				//this 3d is decomposed from 4d
				for (size_t index=3; index < fftPlan->length.size(); index++)
				{
					xyPlan->length.push_back(fftPlan->length[index]);
					xyPlan->inStride.push_back(fftPlan->inStride[index]);
					xyPlan->outStride.push_back(fftPlan->outStride[index]);
				}

				//Set callback data if set on top level plan
				if (fftPlan->hasPreCallback)
				{
					xyPlan->hasPreCallback = true;
					xyPlan->preCallback = fftPlan->preCallback;
					xyPlan->precallUserData = fftPlan->precallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->2D planX failed" ) );

				if( (xyPlan->inStride[0] == 1) && (xyPlan->outStride[0] == 1) &&
					(xyPlan->outStride[2] == Nt*length1) &&
					( ((xyPlan->inStride[2] == Nt*2*length1) && (xyPlan->placeness == CLFFT_INPLACE)) ||
					  ((xyPlan->inStride[2] == length0*length1) && (xyPlan->placeness == CLFFT_OUTOFPLACE)) ) )
				{

					if (fftPlan->tmpBufSize==0)
					{
						fftPlan->tmpBufSize = Nt * length1 * length2 * fftPlan->batchsize * fftPlan->ElementSize();

						for (size_t index=3; index < fftPlan->length.size(); index++)
						{
							fftPlan->tmpBufSize *= fftPlan->length[index];
						}
					}

					// create first transpose plan
					
					//Transpose 
					// output --> tmp
					size_t transLengths[2] = { length0*length1, length2 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, transLengths ),
						_T( "CreateDefaultPlan for planTX transpose failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->transflag = true;

					transLengths[0] = Nt*length1;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTX, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTX transpose failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					trans1Plan->placeness     = CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->forwardScale  = 1.0f;
					trans1Plan->backwardScale = 1.0f;

					trans1Plan->inStride[0]   = 1;
					trans1Plan->inStride[1]   = Nt*length1;
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = length2;
					trans1Plan->iDist         = xyPlan->oDist;
					trans1Plan->oDist		  = Nt*length1*length2;
					trans1Plan->transOutHorizontal = true;

					trans1Plan->gen           = Transpose_GCN;


					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						trans1Plan->length.push_back(fftPlan->length[index]);
						trans1Plan->inStride.push_back(xyPlan->outStride[index]);
						trans1Plan->outStride.push_back(trans1Plan->oDist);
						trans1Plan->oDist *= fftPlan->length[index];
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTX failed" ) );

					// Create column plan as a row plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planZ, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimZ ] ),
						_T( "CreateDefaultPlan for planZ failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planZ, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					colPlan->outputLayout  = trans1Plan->outputLayout;
					colPlan->inputLayout   = trans1Plan->outputLayout;
					colPlan->placeness     = CLFFT_INPLACE;
					colPlan->length.push_back(Nt*length1);

					colPlan->inStride[0]  = 1;
					colPlan->inStride.push_back(length2);
					colPlan->iDist         = Nt*length1*length2;

					colPlan->outStride[0]  = 1;
					colPlan->outStride.push_back(length2);
					colPlan->oDist         = Nt*length1*length2;

					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = fftPlan->forwardScale;
					colPlan->backwardScale = fftPlan->backwardScale;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope		= fftPlan->envelope;

					colPlan->batchsize    = fftPlan->batchsize;

					//this 2d is decomposed from 3d
					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->inStride.push_back(colPlan->iDist);
						colPlan->outStride.push_back(colPlan->oDist);
						colPlan->iDist *= fftPlan->length[index];
						colPlan->oDist *= fftPlan->length[index];
					}

					OPENCL_V(clfftBakePlan(fftPlan->planZ, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planZ failed" ) );

					if (fftPlan->transposed == CLFFT_TRANSPOSED)
					{
						fftPlan->baked = true;
						return	CLFFT_SUCCESS;
					}

					// create second transpose plan
					
					//Transpose 
					//output --> tmp
					size_t trans2Lengths[2] = { length2, length0*length1 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTY, fftPlan->context, CLFFT_2D, trans2Lengths ),
						_T( "CreateDefaultPlan for planTY transpose failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTY, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->transflag = true;

					trans2Lengths[1] = Nt*length1;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTY, CLFFT_2D, trans2Lengths ),
						_T( "clfftSetPlanLength for planTY transpose failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans2Plan->outputLayout = CLFFT_COMPLEX_PLANAR;
							trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					default: assert(false);
					}

					trans2Plan->placeness     = CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->forwardScale  = 1.0f;
					trans2Plan->backwardScale = 1.0f;

					trans2Plan->inStride[0]   = 1;
					trans2Plan->inStride[1]   = length2;
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = Nt*length1;
					trans2Plan->iDist         = Nt*length1*length2;
					trans2Plan->oDist		  = fftPlan->oDist;

					trans2Plan->gen           = Transpose_GCN;
					trans2Plan->transflag     = true;

					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						trans2Plan->length.push_back(fftPlan->length[index]);
						trans2Plan->inStride.push_back(trans2Plan->iDist);
						trans2Plan->iDist *= fftPlan->length[index];
						trans2Plan->outStride.push_back(fftPlan->outStride[index]);
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						trans2Plan->hasPostCallback = true;
						trans2Plan->postCallbackParam = fftPlan->postCallbackParam;
						trans2Plan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTY, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTY failed" ) );


				}
				else
				{

					clLengths[0] = fftPlan->length[ DimZ ];
					clLengths[1] = clLengths[2] = 0;
					//create 1D col plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planZ, fftPlan->context, CLFFT_1D, clLengths ),
						_T( "CreateDefaultPlan for planZ failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planZ, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					switch(fftPlan->outputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_PLANAR;
							colPlan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					colPlan->placeness     = CLFFT_INPLACE;
					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = fftPlan->forwardScale;
					colPlan->backwardScale = fftPlan->backwardScale;
					colPlan->tmpBufSize    = fftPlan->tmpBufSize;

					colPlan->gen			 = fftPlan->gen;
					colPlan->envelope			 = fftPlan->envelope;

					// This is a column FFT, the first elements distance between each FFT is the distance of the first two
					// elements in the original buffer. Like a transpose of the matrix
					colPlan->batchsize = fftPlan->batchsize;
					colPlan->inStride[0] = fftPlan->outStride[2];
					colPlan->outStride[0] = fftPlan->outStride[2];

					//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
					colPlan->length.push_back(1 + fftPlan->length[0]/2);
					colPlan->length.push_back(fftPlan->length[1]);
					colPlan->inStride.push_back(fftPlan->outStride[0]);
					colPlan->inStride.push_back(fftPlan->outStride[1]);
					colPlan->outStride.push_back(fftPlan->outStride[0]);
					colPlan->outStride.push_back(fftPlan->outStride[1]);
					colPlan->iDist    = fftPlan->oDist;
					colPlan->oDist    = fftPlan->oDist;

					//this 3d is decomposed from 4d
					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->inStride.push_back(xyPlan->outStride[index]);
						colPlan->outStride.push_back(fftPlan->outStride[index]);
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						colPlan->hasPostCallback = true;
						colPlan->postCallbackParam = fftPlan->postCallbackParam;
						colPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planZ, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->1D planZ failed" ) );
				}
			}
			else if(fftPlan->outputLayout == CLFFT_REAL)
			{
				size_t length0 = fftPlan->length[ DimX ];
				size_t length1 = fftPlan->length[ DimY ];
				size_t length2 = fftPlan->length[ DimZ ];

				size_t Nt = (1 + length0/2);

				if (fftPlan->tmpBufSize == 0)
				{
					fftPlan->tmpBufSize = Nt * length1 * length2 * fftPlan->batchsize * fftPlan->ElementSize();
					for (size_t index=3; index < fftPlan->length.size(); index++)
						fftPlan->tmpBufSize *= fftPlan->length[index];
				}

				if ((fftPlan->tmpBufSizeC2R==0) && (fftPlan->placeness == CLFFT_OUTOFPLACE))
				{
					fftPlan->tmpBufSizeC2R = fftPlan->tmpBufSize;
				}

				if( (fftPlan->inStride[0] == 1) && (fftPlan->outStride[0] == 1) &&
					( ((fftPlan->outStride[2] == Nt*2*length1) && (fftPlan->oDist == Nt*2*length1*length2) && (fftPlan->placeness == CLFFT_INPLACE)) ||
						((fftPlan->outStride[2] == length0*length1) && (fftPlan->oDist == length0*length1*length2) && (fftPlan->placeness == CLFFT_OUTOFPLACE)) )
					&& (fftPlan->inStride[2] == Nt*length1) && (fftPlan->iDist == Nt*length1*length2))
				{
					// create first transpose plan
					
					//Transpose 
					// input --> tmp
					size_t transLengths[2] = { length0*length1, length2 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTZ, fftPlan->context, CLFFT_2D, transLengths ),
						_T( "CreateDefaultPlan for planTZ transpose failed" ) );

					FFTPlan* trans1Plan	= NULL;
					lockRAII* trans1Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTZ, trans1Plan, trans1Lock ), _T( "fftRepo.getPlan failed" ) );

					trans1Plan->transflag = true;

					transLengths[0] = Nt*length1;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTZ, CLFFT_2D, transLengths ),
						_T( "clfftSetPlanLength for planTZ transpose failed" ) );

					switch(fftPlan->inputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							trans1Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							trans1Plan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					trans1Plan->placeness     = CLFFT_OUTOFPLACE;
					trans1Plan->precision     = fftPlan->precision;
					trans1Plan->tmpBufSize    = 0;
					trans1Plan->batchsize     = fftPlan->batchsize;
					trans1Plan->envelope	  = fftPlan->envelope;
					trans1Plan->forwardScale  = 1.0f;
					trans1Plan->backwardScale = 1.0f;

					trans1Plan->inStride[0]   = 1;
					trans1Plan->inStride[1]   = Nt*length1;
					trans1Plan->outStride[0]  = 1;
					trans1Plan->outStride[1]  = length2;
					trans1Plan->iDist         = fftPlan->iDist;
					trans1Plan->oDist		  = Nt*length1*length2;
					trans1Plan->transOutHorizontal = true;

					trans1Plan->gen           = Transpose_GCN;


					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						trans1Plan->length.push_back(fftPlan->length[index]);
						trans1Plan->inStride.push_back(fftPlan->inStride[index]);
						trans1Plan->outStride.push_back(trans1Plan->oDist);
						trans1Plan->oDist *= fftPlan->length[index];
					}

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						trans1Plan->hasPreCallback = true;
						trans1Plan->preCallback = fftPlan->preCallback;
						trans1Plan->precallUserData = fftPlan->precallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planTZ, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTZ failed" ) );

					// create col plan
					// complex to complex

					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planZ, fftPlan->context, CLFFT_1D, &fftPlan->length[ DimZ ] ),
						_T( "CreateDefaultPlan for planZ failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planZ, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					colPlan->length.push_back(Nt*length1);

					colPlan->inStride[0]  = 1;
					colPlan->inStride.push_back(length2);
					colPlan->iDist        = trans1Plan->oDist;

					colPlan->placeness = CLFFT_INPLACE;
					colPlan->inputLayout = CLFFT_COMPLEX_INTERLEAVED;
					colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;

					colPlan->outStride[0]  = colPlan->inStride[0];
					colPlan->outStride.push_back(colPlan->inStride[1]);
					colPlan->oDist         = colPlan->iDist;

					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						colPlan->length.push_back(fftPlan->length[index]);
						colPlan->inStride.push_back(trans1Plan->outStride[index-1]);
						colPlan->outStride.push_back(trans1Plan->outStride[index-1]);
					}


					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = 1.0f;
					colPlan->backwardScale = 1.0f;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			= fftPlan->gen;
					colPlan->envelope		= fftPlan->envelope;

					colPlan->batchsize = fftPlan->batchsize;

					OPENCL_V(clfftBakePlan(fftPlan->planZ, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planZ failed" ) );

					// create second transpose plan
					
					//Transpose 
					//tmp --> output
					size_t trans2Lengths[2] = { length2, length0*length1 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planTX, fftPlan->context, CLFFT_2D, trans2Lengths ),
						_T( "CreateDefaultPlan for planTX transpose failed" ) );

					FFTPlan* trans2Plan	= NULL;
					lockRAII* trans2Lock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planTX, trans2Plan, trans2Lock ), _T( "fftRepo.getPlan failed" ) );

					trans2Plan->transflag = true;

					trans2Lengths[1] = Nt*length1;
					OPENCL_V(clfftSetPlanLength( fftPlan->planTX, CLFFT_2D, trans2Lengths ),
						_T( "clfftSetPlanLength for planTX transpose failed" ) );


					trans2Plan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
					trans2Plan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;


					trans2Plan->placeness     = CLFFT_OUTOFPLACE;
					trans2Plan->precision     = fftPlan->precision;
					trans2Plan->tmpBufSize    = 0;
					trans2Plan->batchsize     = fftPlan->batchsize;
					trans2Plan->envelope	  = fftPlan->envelope;
					trans2Plan->forwardScale  = 1.0f;
					trans2Plan->backwardScale = 1.0f;

					trans2Plan->inStride[0]   = 1;
					trans2Plan->inStride[1]   = length2;
					trans2Plan->outStride[0]  = 1;
					trans2Plan->outStride[1]  = Nt*length1;
					trans2Plan->iDist         = colPlan->oDist;
					trans2Plan->oDist		  = Nt*length1*length2;

					trans2Plan->gen           = Transpose_GCN;
					trans2Plan->transflag     = true;

					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						trans2Plan->length.push_back(fftPlan->length[index]);
						trans2Plan->inStride.push_back(colPlan->outStride[index-1]);
						trans2Plan->outStride.push_back(trans2Plan->oDist);
						trans2Plan->oDist *= fftPlan->length[index];

					}

					OPENCL_V(clfftBakePlan(fftPlan->planTX, numQueues, commQueueFFT, NULL, NULL ),
						_T( "BakePlan for planTX failed" ) );

					// create row plan
					// hermitian to real

					//create 2D xy plan
					size_t clLengths[] = { length0, length1, 0 };
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan for 2D planX failed" ) );

					FFTPlan* rowPlan	= NULL;
					lockRAII* rowLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, rowPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

					rowPlan->outputLayout  = fftPlan->outputLayout;
					rowPlan->inputLayout   = CLFFT_HERMITIAN_INTERLEAVED;

					rowPlan->length.push_back(length2);

					rowPlan->outStride[0]  = fftPlan->outStride[0];
					rowPlan->outStride[1]  = fftPlan->outStride[1];
					rowPlan->outStride.push_back(fftPlan->outStride[2]);
					rowPlan->oDist         = fftPlan->oDist;

					rowPlan->inStride[0]  = trans2Plan->outStride[0];
					rowPlan->inStride[1]  = Nt;
					rowPlan->inStride.push_back(Nt*length1);
					rowPlan->iDist         = trans2Plan->oDist;

					for (size_t index=3; index < fftPlan->length.size(); index++)
					{
						rowPlan->length.push_back(fftPlan->length[index]);
						rowPlan->inStride.push_back(trans2Plan->outStride[index-1]);
						rowPlan->outStride.push_back(fftPlan->outStride[index]);
					}

					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						rowPlan->placeness     = CLFFT_INPLACE;
					}
					else
					{
						rowPlan->placeness     = CLFFT_OUTOFPLACE;
					}				


					rowPlan->precision     = fftPlan->precision;
					rowPlan->forwardScale  = fftPlan->forwardScale;
					rowPlan->backwardScale = fftPlan->backwardScale;
					rowPlan->tmpBufSize    = 0;

					rowPlan->gen			= fftPlan->gen;
					rowPlan->envelope		= fftPlan->envelope;

					rowPlan->batchsize    = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						rowPlan->hasPostCallback = true;
						rowPlan->postCallbackParam = fftPlan->postCallbackParam;
						rowPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan for planX failed" ) );
				}
				else
				{

					size_t clLengths[] = { 1, 0, 0 };

					clLengths[0] = fftPlan->length[ DimZ ];

					//create 1D col plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planZ, fftPlan->context, CLFFT_1D, clLengths ),
						_T( "CreateDefaultPlan for planZ failed" ) );

					FFTPlan* colPlan	= NULL;
					lockRAII* colLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planZ, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

					switch(fftPlan->inputLayout)
					{
					case CLFFT_HERMITIAN_INTERLEAVED:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_INTERLEAVED;
						}
						break;
					case CLFFT_HERMITIAN_PLANAR:
						{
							colPlan->outputLayout = CLFFT_COMPLEX_INTERLEAVED;
							colPlan->inputLayout  = CLFFT_COMPLEX_PLANAR;
						}
						break;
					default: assert(false);
					}

					colPlan->length.push_back(Nt);
					colPlan->length.push_back(length1);

					colPlan->inStride[0]  = fftPlan->inStride[2];
					colPlan->inStride.push_back(fftPlan->inStride[0]);
					colPlan->inStride.push_back(fftPlan->inStride[1]);
					colPlan->iDist         = fftPlan->iDist;


					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						colPlan->placeness = CLFFT_INPLACE;

						colPlan->outStride[0]  = colPlan->inStride[0];
						colPlan->outStride.push_back(colPlan->inStride[1]);
						colPlan->outStride.push_back(colPlan->inStride[2]);
						colPlan->oDist         = colPlan->iDist;

						for (size_t index=3; index < fftPlan->length.size(); index++)
						{
							colPlan->length.push_back(fftPlan->length[index]);
							colPlan->inStride.push_back(fftPlan->inStride[index]);
							colPlan->outStride.push_back(fftPlan->inStride[index]);
						}
					}
					else
					{
						colPlan->placeness = CLFFT_OUTOFPLACE;

						colPlan->outStride[0]  = Nt*length1;
						colPlan->outStride.push_back(1);
						colPlan->outStride.push_back(Nt);
						colPlan->oDist         = Nt*length1*length2;

						for (size_t index=3; index < fftPlan->length.size(); index++)
						{
							colPlan->length.push_back(fftPlan->length[index]);
							colPlan->inStride.push_back(fftPlan->inStride[index]);
							colPlan->outStride.push_back(colPlan->oDist);
							colPlan->oDist *= fftPlan->length[index];
						}
					}

				
					colPlan->precision     = fftPlan->precision;
					colPlan->forwardScale  = 1.0f;
					colPlan->backwardScale = 1.0f;
					colPlan->tmpBufSize    = 0;

					colPlan->gen			 = fftPlan->gen;
					colPlan->envelope		 = fftPlan->envelope;

					colPlan->batchsize = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPreCallback)
					{
						colPlan->hasPreCallback = true;
						colPlan->preCallback = fftPlan->preCallback;
						colPlan->precallUserData = fftPlan->precallUserData;
					}
				
					OPENCL_V(clfftBakePlan(fftPlan->planZ, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->1D planZ failed" ) );


					clLengths[0] = fftPlan->length[ DimX ];
					clLengths[1] = fftPlan->length[ DimY ];

					//create 2D xy plan
					OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_2D, clLengths ),
						_T( "CreateDefaultPlan 2D planX failed" ) );

					FFTPlan* xyPlan	= NULL;
					lockRAII* rowLock	= NULL;
					OPENCL_V( fftRepo.getPlan( fftPlan->planX, xyPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

					xyPlan->inputLayout   = CLFFT_HERMITIAN_INTERLEAVED;
					xyPlan->outputLayout  = fftPlan->outputLayout;

					xyPlan->length.push_back(length2);
			
					xyPlan->outStride[0]  = fftPlan->outStride[0];
					xyPlan->outStride[1]  = fftPlan->outStride[1];
					xyPlan->outStride.push_back(fftPlan->outStride[2]);
					xyPlan->oDist         = fftPlan->oDist;

					if (fftPlan->placeness == CLFFT_INPLACE)
					{
						xyPlan->placeness     = CLFFT_INPLACE;

						xyPlan->inStride[0]  = colPlan->outStride[1];
						xyPlan->inStride[1]  = colPlan->outStride[2];
						xyPlan->inStride.push_back(colPlan->outStride[0]);
						xyPlan->iDist         = colPlan->oDist;

						for (size_t index=3; index < fftPlan->length.size(); index++)
						{
							xyPlan->length.push_back(fftPlan->length[index]);
							xyPlan->inStride.push_back(colPlan->outStride[index]);
							xyPlan->outStride.push_back(fftPlan->outStride[index]);
						}
					}
					else
					{
						xyPlan->placeness     = CLFFT_OUTOFPLACE;

						xyPlan->inStride[0]   = 1;
						xyPlan->inStride[1]   = Nt;
						xyPlan->inStride.push_back(Nt*length1);
						xyPlan->iDist         = Nt*length1*length2;

						for (size_t index=3; index < fftPlan->length.size(); index++)
						{
							xyPlan->length.push_back(fftPlan->length[index]);
							xyPlan->outStride.push_back(fftPlan->outStride[index]);
							xyPlan->inStride.push_back(xyPlan->iDist);						
							xyPlan->iDist *= fftPlan->length[index];
						}
					}


					xyPlan->precision     = fftPlan->precision;
					xyPlan->forwardScale  = fftPlan->forwardScale;
					xyPlan->backwardScale = fftPlan->backwardScale;
					xyPlan->tmpBufSize    = fftPlan->tmpBufSize;

					xyPlan->gen			 = fftPlan->gen;
					xyPlan->envelope	 = fftPlan->envelope;

					xyPlan->batchsize    = fftPlan->batchsize;

					//Set callback data if set on top level plan
					if (fftPlan->hasPostCallback)
					{
						xyPlan->hasPostCallback = true;
						xyPlan->postCallbackParam = fftPlan->postCallbackParam;
						xyPlan->postcallUserData = fftPlan->postcallUserData;
					}

					OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->2D planX failed" ) );
				}
			}
			else
			{
				if (fftPlan->tmpBufSize==0 && (
					fftPlan->length[0] > Large1DThreshold ||
					fftPlan->length[1] > Large1DThreshold ||
					fftPlan->length[2] > Large1DThreshold
					))
				{
					fftPlan->tmpBufSize = fftPlan->length[0] * fftPlan->length[1] * fftPlan->length[2] *
						fftPlan->batchsize * fftPlan->ElementSize();
				}

				size_t clLengths[] = { 1, 1, 0 };
				clLengths[0] = fftPlan->length[ DimX ];
				clLengths[1] = fftPlan->length[ DimY ];

				//create 2D xy plan
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planX, fftPlan->context, CLFFT_2D, clLengths ),
					_T( "CreateDefaultPlan 2D planX failed" ) );

				FFTPlan* xyPlan	= NULL;
				lockRAII* rowLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planX, xyPlan, rowLock ), _T( "fftRepo.getPlan failed" ) );

				xyPlan->inputLayout   = fftPlan->inputLayout;
				xyPlan->outputLayout  = fftPlan->outputLayout;
				xyPlan->placeness     = fftPlan->placeness;
				xyPlan->precision     = fftPlan->precision;
				xyPlan->forwardScale  = 1.0f;
				xyPlan->backwardScale = 1.0f;
				xyPlan->tmpBufSize    = fftPlan->tmpBufSize;

				xyPlan->gen			 = fftPlan->gen;
				xyPlan->envelope			 = fftPlan->envelope;

				// This is the xy fft, the first elements distance between the first two FFTs is the distance of the first elements
				// of the first two rows in the original buffer.
				xyPlan->batchsize    = fftPlan->batchsize;
				xyPlan->inStride[0]  = fftPlan->inStride[0];
				xyPlan->inStride[1]  = fftPlan->inStride[1];
				xyPlan->outStride[0] = fftPlan->outStride[0];
				xyPlan->outStride[1] = fftPlan->outStride[1];

				//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
				xyPlan->length.push_back(fftPlan->length[2]);
				xyPlan->inStride.push_back(fftPlan->inStride[2]);
				xyPlan->outStride.push_back(fftPlan->outStride[2]);
				xyPlan->iDist    = fftPlan->iDist;
				xyPlan->oDist    = fftPlan->oDist;

				//Set callback data if set on top level plan
				if (fftPlan->hasPreCallback)
				{
					xyPlan->hasPreCallback = true;
					xyPlan->preCallback = fftPlan->preCallback;
					xyPlan->precallUserData = fftPlan->precallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planX, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->2D planX failed" ) );

				clLengths[0] = fftPlan->length[ DimZ ];
				clLengths[1] = clLengths[2] = 0;
				//create 1D col plan
				OPENCL_V(clfftCreateDefaultPlanInternal( &fftPlan->planZ, fftPlan->context, CLFFT_1D, clLengths ),
					_T( "CreateDefaultPlan for planZ failed" ) );

				FFTPlan* colPlan	= NULL;
				lockRAII* colLock	= NULL;
				OPENCL_V( fftRepo.getPlan( fftPlan->planZ, colPlan, colLock ), _T( "fftRepo.getPlan failed" ) );

				colPlan->inputLayout   = fftPlan->outputLayout;
				colPlan->outputLayout  = fftPlan->outputLayout;
				colPlan->placeness     = CLFFT_INPLACE;
				colPlan->precision     = fftPlan->precision;
				colPlan->forwardScale  = fftPlan->forwardScale;
				colPlan->backwardScale = fftPlan->backwardScale;
				colPlan->tmpBufSize    = fftPlan->tmpBufSize;

				colPlan->gen			 = fftPlan->gen;
				colPlan->envelope			 = fftPlan->envelope;

				// This is a column FFT, the first elements distance between each FFT is the distance of the first two
				// elements in the original buffer. Like a transpose of the matrix
				colPlan->batchsize = fftPlan->batchsize;
				colPlan->inStride[0] = fftPlan->outStride[2];
				colPlan->outStride[0] = fftPlan->outStride[2];

				//pass length and other info to kernel, so the kernel knows this is decomposed from higher dimension
				colPlan->length.push_back(fftPlan->length[0]);
				colPlan->length.push_back(fftPlan->length[1]);
				colPlan->inStride.push_back(fftPlan->outStride[0]);
				colPlan->inStride.push_back(fftPlan->outStride[1]);
				colPlan->outStride.push_back(fftPlan->outStride[0]);
				colPlan->outStride.push_back(fftPlan->outStride[1]);
				colPlan->iDist    = fftPlan->oDist;
				colPlan->oDist    = fftPlan->oDist;

				//Set callback data if set on top level plan
				if (fftPlan->hasPostCallback)
				{
					colPlan->hasPostCallback = true;
					colPlan->postCallbackParam = fftPlan->postCallbackParam;
					colPlan->postcallUserData = fftPlan->postcallUserData;
				}

				OPENCL_V(clfftBakePlan(fftPlan->planZ, numQueues, commQueueFFT, NULL, NULL ), _T( "BakePlan 3D->1D planZ failed" ) );
			}

			fftPlan->baked = true;
			return	CLFFT_SUCCESS;
		}
	}

	
	clfftStatus err = selectAction(fftPlan, fftPlan->action, commQueueFFT);

	//	Allocate resources
	OPENCL_V( fftPlan->AllocateBuffers (), _T("AllocateBuffers() failed"));

	fftPlan->ConstructAndEnqueueConstantBuffers( commQueueFFT );

	//	Record that we baked the plan
	fftPlan->baked		= true;

	return	CLFFT_SUCCESS;
}

clfftStatus clfftCopyPlan( clfftPlanHandle* out_plHandle, cl_context new_context, clfftPlanHandle in_plHandle )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* in_fftPlan	= NULL, *out_fftPlan = NULL;
	lockRAII* in_planLock = NULL, *out_planLock = NULL;

	OPENCL_V( fftRepo.getPlan( in_plHandle, in_fftPlan, in_planLock ), _T( "fftRepo.getPlan failed" ) );

	OPENCL_V( clfftCreateDefaultPlan( out_plHandle, new_context, in_fftPlan->dim, &in_fftPlan->length[ 0 ] ),
		_T( "clfftCreateDefaultPlan failed" ) );

	OPENCL_V( fftRepo.getPlan( *out_plHandle, out_fftPlan, out_planLock ), _T( "fftRepo.getPlan failed" ) );

	//	Let other operations complete before attempting to copy the plan
	scopedLock sLock( *in_planLock, _T( "clfftCopyPlan" ) );

	out_fftPlan->baked = false;
	out_fftPlan->gen = in_fftPlan->gen;
	out_fftPlan->envelope = in_fftPlan->envelope;
	out_fftPlan->dim = in_fftPlan->dim;
	out_fftPlan->inputLayout = in_fftPlan->inputLayout;
	out_fftPlan->outputLayout = in_fftPlan->outputLayout;
	out_fftPlan->placeness = in_fftPlan->placeness;
	out_fftPlan->precision = in_fftPlan->precision;
	out_fftPlan->forwardScale = in_fftPlan->forwardScale;
	out_fftPlan->backwardScale = in_fftPlan->backwardScale;
	out_fftPlan->iDist = in_fftPlan->iDist;
	out_fftPlan->oDist = in_fftPlan->oDist;
	out_fftPlan->length = in_fftPlan->length;
	out_fftPlan->inStride = in_fftPlan->inStride;
	out_fftPlan->outStride = in_fftPlan->outStride;
	out_fftPlan->batchsize = in_fftPlan->batchsize;
	out_fftPlan->transposed = in_fftPlan->transposed;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTPlan::ConstructAndEnqueueConstantBuffers( cl_command_queue* commQueueFFT )
{
	//	Construct the constant buffer and call clEnqueueWriteBuffer
	//
	cb_t ConstantBufferParams [CLFFT_CB_SIZE];
	memset (& ConstantBufferParams, 0, sizeof (ConstantBufferParams));

	ConstantBufferParams[0].u = std::max<cl_uint> (1, cl_uint (/*fftPlan->*/batchsize));


	OPENCL_V(clEnqueueWriteBuffer( *commQueueFFT,
		/*fftPlan->*/const_buffer,
		1,		// TODO? non-blocking write?
		0,
		sizeof(ConstantBufferParams),
		&ConstantBufferParams,
		0,
		NULL,
		NULL), _T("clEnqueueWriteBuffer failed") );

	return CLFFT_SUCCESS;
}


clfftStatus	clfftDestroyPlan( clfftPlanHandle* plHandle )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	OPENCL_V( fftRepo.getPlan( *plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );

	//	Recursively destroy subplans, that are used for higher dimensional FFT's
	if( fftPlan->planX )
		clfftDestroyPlan( &fftPlan->planX );
	if( fftPlan->planY )
		clfftDestroyPlan( &fftPlan->planY );
	if( fftPlan->planZ )
		clfftDestroyPlan( &fftPlan->planZ );
	if( fftPlan->planTX )
		clfftDestroyPlan( &fftPlan->planTX );
	if( fftPlan->planTY )
		clfftDestroyPlan( &fftPlan->planTY );
	if( fftPlan->planTZ )
		clfftDestroyPlan( &fftPlan->planTZ );
	if( fftPlan->planRCcopy )
		clfftDestroyPlan( &fftPlan->planRCcopy );
	if( fftPlan->planCopy )
		clfftDestroyPlan( &fftPlan->planCopy );

	fftRepo.deletePlan( plHandle );

	return	CLFFT_SUCCESS;
}

//	This routine will query the OpenCL context for it's devices
//	and their hardware limitations, which we synthesize into a
//	hardware "envelope".
//	We only query the devices the first time we're called after
//	the object's context is set.  On 2nd and subsequent calls,
//	we just return the pointer.
//
clfftStatus FFTPlan::SetEnvelope ()
{

	// TODO  The caller has already acquired the lock on *this
	//	However, we shouldn't depend on it.

	if (0 == envelope.limit_LocalMemSize) do {
		//	First time, query OpenCL for the device info
		//
		memset (&envelope, 0, sizeof(envelope));

		//	Get the size needed for the device list
		//
		size_t deviceListSize = 0;
		OPENCL_V( ::clGetContextInfo( context, CL_CONTEXT_DEVICES, 0, NULL, &deviceListSize ),
			_T("Getting device array size ( ::clGetContextInfo() )" ));
		cl_uint n = cl_uint (deviceListSize / sizeof(cl_device_id));
		if (n == 0) break;

		std::vector< cl_device_id > devices( n+1 );
		//	Get the device list
		//
		OPENCL_V( ::clGetContextInfo( context, CL_CONTEXT_DEVICES, deviceListSize, &devices[ 0 ], NULL ),
			_T("Getting device array ( ::clGetContextInfo() )") );

		//	Get the # of devices
		//
		cl_uint cContextDevices	= 0;

		size_t deviceVersionSize	= 0;
		OPENCL_V( ::clGetDeviceInfo( devices[0], CL_DEVICE_VERSION, 0, NULL, &deviceVersionSize ),
			_T("Getting CL_DEVICE_VERSION Info string size ( ::clGetDeviceInfo() )" ));

		std::vector< char > szDeviceVersion( deviceVersionSize );
		OPENCL_V( ::clGetDeviceInfo( devices[0], CL_DEVICE_VERSION, deviceVersionSize, &szDeviceVersion[ 0 ], NULL ),
			_T("Getting CL_DEVICE_VERSION Platform Info string ( ::clGetDeviceInfo() )" ));

		char openclstr[11]="OpenCL 1.0";

		if (!strncmp((const char*)&szDeviceVersion[ 0 ], openclstr, 10))
		{
			cContextDevices	= 1;
		}
		else
		{
			OPENCL_V( ::clGetContextInfo( context, CL_CONTEXT_NUM_DEVICES, sizeof( cContextDevices ), &cContextDevices, NULL ),
				_T("Getting number of context devices ( ::clGetContextInfo() )" ));
		}

		cContextDevices = std::min<cl_uint> (cContextDevices, n);
		if (0 == cContextDevices)
			break;

		envelope.limit_LocalMemSize  = 32768;
		envelope.limit_WorkGroupSize = 256;
		envelope.limit_Dimensions    = countOf (envelope.limit_Size);
		for (size_t u = 0; u < countOf (envelope.limit_Size); ++u) {
			envelope.limit_Size[u] = 256;
		}

		for( cl_uint i = 0; i < cContextDevices; ++i )
		{
			cl_device_id devId = devices[i];

			cl_ulong memsize = 0;
			unsigned int maxdim = 0;
			size_t temp[countOf (envelope.limit_Size)];
			memset (&temp, 0, sizeof(temp));

			OPENCL_V( ::clGetDeviceInfo( devId, CL_DEVICE_LOCAL_MEM_SIZE, sizeof( cl_ulong ), &memsize, NULL ),
				_T("Getting CL_DEVICE_LOCAL_MEM_SIZE device info ( ::clGetDeviceInfo() )") );
			envelope.limit_LocalMemSize = std::min<size_t> (envelope.limit_LocalMemSize, memsize);

			OPENCL_V( ::clGetDeviceInfo( devId, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof( unsigned int ), &maxdim, NULL ),
				_T("Getting CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS device info ( ::clGetDeviceInfo() )") );
			BUG_CHECK (countOf (envelope.limit_Size) >= maxdim);
			envelope.limit_Dimensions = std::min<size_t> (envelope.limit_Dimensions, maxdim);

			OPENCL_V( ::clGetDeviceInfo( devId, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof( size_t ), &temp[0], NULL ),
				_T("Getting CL_DEVICE_MAX_WORK_GROUP_SIZE device info ( ::clGetDeviceInfo() )") );
			envelope.limit_WorkGroupSize = std::min<size_t> (envelope.limit_WorkGroupSize, temp[0]);

			OPENCL_V( ::clGetDeviceInfo( devId, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof( temp ), &temp[0], NULL ),
				_T("Getting CL_DEVICE_MAX_WORK_ITEM_SIZES device info ( ::clGetDeviceInfo() )") );
			for (size_t u = 0; u < envelope.limit_Dimensions; ++u) {
				BUG_CHECK (temp[u] > 0)
				envelope.limit_Size[u] = std::min<size_t> (envelope.limit_Size[u], temp[u]);
			}
		}

		BUG_CHECK (envelope.limit_LocalMemSize >= 1024)
	} while (0);

	return CLFFT_SUCCESS;
}

clfftStatus FFTPlan::AllocateBuffers ()
{
	cl_int status = CL_SUCCESS;

	assert (NULL == const_buffer);
	ReleaseBuffers ();

	assert(4 == sizeof(int));

	do {
		const_buffer = clCreateBuffer (context,
										CL_MEM_READ_ONLY,
										CLFFT_CB_SIZE * sizeof (int),
										0,
										&status);
		if (CL_SUCCESS != status)
			break;
	} while (0);

	return	(clfftStatus) status;
}

clfftStatus FFTPlan::ReleaseBuffers ()
{
	clfftStatus result = CLFFT_SUCCESS;
	clfftStatus tmp;

	if( NULL != const_buffer )
	{
		tmp = static_cast< clfftStatus >( clReleaseMemObject( const_buffer ) );
		const_buffer = NULL;
		if( CLFFT_SUCCESS == result )
			result = tmp;
	}

	if( (NULL != intBuffer) && libCreatedIntBuffer )
	{
		tmp = static_cast< clfftStatus >( clReleaseMemObject( intBuffer ) );
		intBuffer = NULL;
		if( CLFFT_SUCCESS == result )
			result = tmp;
	}

	if( NULL != intBufferRC )
	{
		tmp = static_cast< clfftStatus >( clReleaseMemObject( intBufferRC ) );
		intBufferRC = NULL;
		if( CLFFT_SUCCESS == result )
			result = tmp;
	}
	
	if( NULL != intBufferC2R )
	{
		tmp = static_cast< clfftStatus >( clReleaseMemObject( intBufferC2R ) );
		intBufferC2R = NULL;
		if( CLFFT_SUCCESS == result )
			result = tmp;
	}

	return	result;
}



clfftStatus FFTPlan::GetMax1DLength (size_t *longest ) const
{
	switch(gen)
	{
	case Stockham:		return GetMax1DLengthStockham(longest);
    case Transpose_GCN:			*longest = 4096; return CLFFT_SUCCESS;
    case Transpose_SQUARE:		*longest = 4096; return CLFFT_SUCCESS;
	case Transpose_NONSQUARE:	*longest = 4096; return CLFFT_SUCCESS;
    case Copy:					*longest = 4096; return CLFFT_SUCCESS;
	default:			assert(false); return CLFFT_NOTIMPLEMENTED;
	}
}

clfftStatus FFTPlan::GetEnvelope (const FFTEnvelope ** ppEnvelope) const
{
	if( &envelope == NULL )
    { 
        assert( false );
        return CLFFT_NOTIMPLEMENTED;
    }

	*ppEnvelope = &envelope;
	return CLFFT_SUCCESS;
}

size_t FFTPlan::ElementSize() const
{
	return ( ((precision == CLFFT_DOUBLE) || (precision == CLFFT_DOUBLE_FAST)) ? sizeof( std::complex<double> ) : sizeof( std::complex<float> ) );
}

