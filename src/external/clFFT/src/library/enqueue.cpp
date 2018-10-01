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

#include "stdafx.h"
#include <math.h>
#include "private.h"
#include "repo.h"
#include "plan.h"
#include "generator.stockham.h"
#include "../include/convenienceFunctions.h"

#include "action.h"
#include "fft_binary_lookup.h"

#define FFT_CACHE_DEBUG 0



FFTCopyAction::FFTCopyAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTAction(plan, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit constructor
        return;
    }

    err = CLFFT_SUCCESS;
}

FFTTransposeGCNAction::FFTTransposeGCNAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTAction(plan, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit constructor
        return;
    }

    err = CLFFT_SUCCESS;
}

FFTTransposeSquareAction::FFTTransposeSquareAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTAction(plan, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit constructor
        return;
    }

    err = CLFFT_SUCCESS;
}

FFTTransposeNonSquareAction::FFTTransposeNonSquareAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTAction(plan, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit constructor
        return;
    }

    err = CLFFT_SUCCESS;
}

FFTStockhamAction::FFTStockhamAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTAction(plan, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit constructor
        return;
    }

    err = CLFFT_SUCCESS;
}



FFTAction::FFTAction(FFTPlan * fftPlan, clfftStatus & err)
    : plan(fftPlan)
{
    err = CLFFT_SUCCESS;
}

clfftStatus FFTAction::selectBufferArguments(FFTPlan * fftPlan,
                                             cl_mem* clInputBuffers,
                                             cl_mem* clOutputBuffers,
                                             std::vector< cl_mem > &inputBuff,
                                             std::vector< cl_mem > &outputBuff)
{
    
    // 1d with normal length will fall into the below category
    // add: 2d transpose kernel will fall into here too.
    inputBuff.reserve( 2 );
    outputBuff.reserve( 2 );

    //	Decode the relevant properties from the plan paramter to figure out how many input/output buffers we have
    switch( fftPlan->inputLayout )
    {
    case CLFFT_COMPLEX_INTERLEAVED:
    {
        switch( fftPlan->outputLayout )
        {
        case CLFFT_COMPLEX_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_COMPLEX_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                //	Invalid to be an inplace transform, and go from 1 to 2 buffers
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_REAL:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        default:
        {
            //	Don't recognize output layout
            return CLFFT_INVALID_ARG_VALUE;
        }
        }

        break;
    }
    case CLFFT_COMPLEX_PLANAR:
    {
        switch( fftPlan->outputLayout )
        {
        case CLFFT_COMPLEX_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_COMPLEX_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_REAL:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        default:
        {
            //	Don't recognize output layout
            return CLFFT_INVALID_ARG_VALUE;
        }
        }

        break;
    }
    case CLFFT_HERMITIAN_INTERLEAVED:
    {
        switch( fftPlan->outputLayout )
        {
        case CLFFT_COMPLEX_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_COMPLEX_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_INTERLEAVED:
        {
            return CLFFT_INVALID_ARG_VALUE;
        }
        case CLFFT_HERMITIAN_PLANAR:
        {
            return CLFFT_INVALID_ARG_VALUE;
        }
        case CLFFT_REAL:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        default:
        {
            //	Don't recognize output layout
            return CLFFT_INVALID_ARG_VALUE;
        }
        }

        break;
    }
    case CLFFT_HERMITIAN_PLANAR:
    {
        switch( fftPlan->outputLayout )
        {
        case CLFFT_COMPLEX_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_COMPLEX_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_INTERLEAVED:
        {
            return CLFFT_INVALID_ARG_VALUE;
        }
        case CLFFT_HERMITIAN_PLANAR:
        {
            return CLFFT_INVALID_ARG_VALUE;
        }
        case CLFFT_REAL:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                inputBuff.push_back( clInputBuffers[ 1 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        default:
        {
            //	Don't recognize output layout
            return CLFFT_INVALID_ARG_VALUE;
        }
        }

        break;
    }
    case CLFFT_REAL:
    {
        switch( fftPlan->outputLayout )
        {
        case CLFFT_COMPLEX_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_COMPLEX_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_INTERLEAVED:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 0 ] );
            }

            break;
        }
        case CLFFT_HERMITIAN_PLANAR:
        {
            if( fftPlan->placeness == CLFFT_INPLACE )
            {
                return CLFFT_INVALID_ARG_VALUE;
            }
            else
            {
                inputBuff.push_back( clInputBuffers[ 0 ] );

                outputBuff.push_back( clOutputBuffers[ 0 ] );
                outputBuff.push_back( clOutputBuffers[ 1 ] );
            }

            break;
        }
        default:
        {
			if(fftPlan->transflag)
			{
				if( fftPlan->placeness == CLFFT_INPLACE )
				{
					return CLFFT_INVALID_ARG_VALUE;
				}
				else
				{
					inputBuff.push_back( clInputBuffers[ 0 ] );
					outputBuff.push_back( clOutputBuffers[ 0 ] );
				}
			}
			else
			{
				//	Don't recognize output layout
				return CLFFT_INVALID_ARG_VALUE;
			}
        }
        }

        break;
    }
    default:
    {
        //	Don't recognize output layout
        return CLFFT_INVALID_ARG_VALUE;
    }
    }

    return CLFFT_SUCCESS;
}


clfftStatus FFTAction::enqueue(clfftPlanHandle plHandle,
                               clfftDirection dir,
                               cl_uint numQueuesAndEvents,
                               cl_command_queue* commQueues,
                               cl_uint numWaitEvents,
                               const cl_event* waitEvents,
                               cl_event* outEvents,
                               cl_mem* clInputBuffers,
                               cl_mem* clOutputBuffers)
{
    FFTRepo & fftRepo = FFTRepo::getInstance();

    std::vector< cl_mem > inputBuff;
    std::vector< cl_mem > outputBuff;


    clfftStatus status = selectBufferArguments(this->plan,
                                               clInputBuffers, clOutputBuffers,
                                               inputBuff, outputBuff);

    if (status != CLFFT_SUCCESS)
    {
        return status;
    }
    
    //	TODO:  In the case of length == 1, FFT is a trivial NOP, but we still need to apply the forward and backwards tranforms
    //	TODO:  Are map lookups expensive to call here?  We can cache a pointer to the cl_program/cl_kernel in the plan

    //	Translate the user plan into the structure that we use to map plans to clPrograms

    cl_program	prog;
    cl_kernel	kern;
	lockRAII* kernelLock;
    OPENCL_V( fftRepo.getclProgram( this->getGenerator(), this->getSignatureData(), prog, this->plan->bakeDevice, this->plan->context ), _T( "fftRepo.getclProgram failed" ) );
    OPENCL_V( fftRepo.getclKernel( prog, dir, kern, kernelLock), _T( "fftRepo.getclKernels failed" ) );

	scopedLock sLock(*kernelLock, _T("FFTAction::enqueue"));

    cl_uint uarg = 0;
    if (!this->plan->transflag && !(this->plan->gen == Copy))
    {
        //	::clSetKernelArg() is not thread safe, according to the openCL spec for the same cl_kernel object
        //	TODO:  Need to verify that two different plans (which would get through our lock above) with exactly the same
        //	parameters would NOT share the same cl_kernel objects

        /* constant buffer */
        OPENCL_V( clSetKernelArg( kern, uarg++, sizeof( cl_mem ), (void*)&this->plan->const_buffer ), _T( "clSetKernelArg failed" ) );
    }

    //	Input buffer(s)
    //	Input may be 1 buffer  (CLFFT_COMPLEX_INTERLEAVED)
    //	          or 2 buffers (CLFFT_COMPLEX_PLANAR)

    for (size_t i = 0; i < inputBuff.size(); ++i)
    {
        OPENCL_V( clSetKernelArg( kern, uarg++, sizeof( cl_mem ), (void*)&inputBuff[i] ), _T( "clSetKernelArg failed" ) );
    }
    //	Output buffer(s)
    //	Output may be 0 buffers (CLFFT_INPLACE)
    //	           or 1 buffer  (CLFFT_COMPLEX_INTERLEAVED)
    //	           or 2 buffers (CLFFT_COMPLEX_PLANAR)
    for (size_t o = 0; o < outputBuff.size(); ++o)
    {
        OPENCL_V( clSetKernelArg( kern, uarg++, sizeof( cl_mem ), (void*)&outputBuff[o] ), _T( "clSetKernelArg failed" ) );
    }

	//If callback function is set for the plan, pass the appropriate aruments
	if (this->plan->hasPreCallback || this->plan->hasPostCallback)
	{
	if (this->plan->hasPreCallback)
	{
		OPENCL_V( clSetKernelArg( kern, uarg++, sizeof( cl_mem ), (void*)&this->plan->precallUserData ), _T( "clSetKernelArg failed" ) );
		}

		//If post-callback function is set for the plan, pass the appropriate aruments
		if (this->plan->hasPostCallback)
		{
			OPENCL_V( clSetKernelArg( kern, uarg++, sizeof( cl_mem ), (void*)&this->plan->postcallUserData ), _T( "clSetKernelArg failed" ) );
		}

		//Pass LDS size arument if set
		if ((this->plan->hasPreCallback && this->plan->preCallback.localMemSize > 0) || 
			(this->plan->hasPostCallback && this->plan->postCallbackParam.localMemSize > 0))
		{
			int localmemSize = 0;
			if (this->plan->hasPreCallback && this->plan->preCallback.localMemSize > 0)
				localmemSize = this->plan->preCallback.localMemSize;
			if (this->plan->hasPostCallback && this->plan->postCallbackParam.localMemSize > 0)
				localmemSize += this->plan->postCallbackParam.localMemSize;

			OPENCL_V( clSetKernelArg( kern, uarg++, localmemSize, NULL ), _T( "clSetKernelArg failed" ) );
		}
	}

    std::vector< size_t > gWorkSize;
    std::vector< size_t > lWorkSize;
    clfftStatus result = this->getWorkSizes (gWorkSize, lWorkSize);
	//std::cout << "work sizes are " << gWorkSize[0] << ", " << lWorkSize[0] << std::endl;
	/*
	std::cout << "work sizes are ";
	for (auto itor = gWorkSize.begin(); itor != gWorkSize.end(); itor++)
		std::cout << *itor << " ";
	std::cout << ", ";
	for (auto itor = lWorkSize.begin(); itor != lWorkSize.end(); itor++)
		std::cout << *itor << " ";
	std::cout << std::endl;
	*/
    // TODO:  if getWorkSizes returns CLFFT_INVALID_GLOBAL_WORK_SIZE, that means
    // that this multidimensional input data array is too large to be transformed
    // with a single call to clEnqueueNDRangeKernel.  For now, we will just return
    // the error code back up the call stack.
    // The *correct* course of action would be to split the work into mutliple
    // calls to clEnqueueNDRangeKernel.
    if (CLFFT_INVALID_GLOBAL_WORK_SIZE == result)
    {
        OPENCL_V( result, _T("Work size too large for clEnqueNDRangeKernel()"));
    }
    else
    {
        OPENCL_V( result, _T("FFTAction::getWorkSizes failed"));
    }
    BUG_CHECK (gWorkSize.size() == lWorkSize.size());


    cl_int call_status = clEnqueueNDRangeKernel( *commQueues, kern, static_cast< cl_uint >( gWorkSize.size( ) ),
                                            NULL, &gWorkSize[ 0 ],  &lWorkSize[ 0 ], numWaitEvents, waitEvents, outEvents );
    OPENCL_V( call_status, _T( "clEnqueueNDRangeKernel failed" ) );

    if( fftRepo.pStatTimer )
    {
        fftRepo.pStatTimer->AddSample( plHandle, this->plan, kern, numQueuesAndEvents, outEvents, gWorkSize, lWorkSize );
    }

    return CLFFT_SUCCESS;
}



//	Read the kernels that this plan uses from file, and store into the plan
clfftStatus FFTAction::writeKernel( const clfftPlanHandle plHandle, const clfftGenerators gen, const FFTKernelSignatureHeader* data, const cl_context& context, const cl_device_id &device )
{
    FFTRepo& fftRepo	= FFTRepo::getInstance( );

    std::string kernelPath = getKernelName(gen, plHandle, true);

    //	Logic to write string contents out to file
    tofstreamRAII< std::ofstream, std::string > kernelFile( kernelPath.c_str( ) );
    if( !kernelFile.get( ) )
    {
        std::cerr << "Failed to open kernel file for writing: " << kernelPath.c_str( ) << std::endl;
        return CLFFT_FILE_CREATE_FAILURE;
    }

    std::string kernel;
    OPENCL_V( fftRepo.getProgramCode( gen, data, kernel, device, context ), _T( "fftRepo.getProgramCode failed." ) );

    kernelFile.get( ) << kernel << std::endl;

    return	CLFFT_SUCCESS;
}


// **************** TODO TODO TODO ***********************
// Making compileKernels function take in command queue parameter so we can build for 1 particular device only;
// this may not be desirable for persistent plans, where we may have to compile for all devices in the context;
// make changes appropriately before enabling persistent plans and then remove this comment

//	Compile the kernels that this plan uses, and store into the plan
clfftStatus FFTAction::compileKernels( const cl_command_queue commQueueFFT, const clfftPlanHandle plHandle, FFTPlan* fftPlan )
{
    cl_int status = 0;
    size_t deviceListSize = 0;

    FFTRepo& fftRepo	= FFTRepo::getInstance( );

    // create a cl program executable for the device associated with command queue
    // Get the device
    cl_device_id &q_device = fftPlan->bakeDevice;

    cl_program program;
    if( fftRepo.getclProgram( this->getGenerator(), this->getSignatureData(), program, q_device, fftPlan->context ) == CLFFT_INVALID_PROGRAM )
    {
        FFTBinaryLookup lookup (this->getGenerator(), plHandle, fftPlan->context, q_device);

        lookup.variantRaw(this->getSignatureData(), this->getSignatureData()->datasize);

        if (lookup.found())
        {
#if FFT_CACHE_DEBUG
            // debug message in debug mode to ensure that the cache is used
            fprintf(stderr, "Kernel loaded from cache\n");
#endif

            program = lookup.getProgram();
        }
        else
        {
#if FFT_CACHE_DEBUG
            fprintf(stderr, "Kernel built from source\n");
#endif

            //	If the user wishes us to write the kernels out to disk, we do so
            if( fftRepo.setupData.debugFlags & CLFFT_DUMP_PROGRAMS )
            {
				OPENCL_V( writeKernel( plHandle, this->getGenerator(), this->getSignatureData(), fftPlan->context, fftPlan->bakeDevice ), _T( "writeKernel failed." ) );
            }

            std::string programCode;
            OPENCL_V( fftRepo.getProgramCode( this->getGenerator(), this->getSignatureData(), programCode, q_device, fftPlan->context  ), _T( "fftRepo.getProgramCode failed." ) );

            const char* source = programCode.c_str();
            program = clCreateProgramWithSource( fftPlan->context, 1, &source, NULL, &status );
            OPENCL_V( status, _T( "clCreateProgramWithSource failed." ) );

            // create a cl program executable for the device associated with command queue

#if defined(DEBUGGING)
            status = clBuildProgram( program, 1, &q_device, "-g -cl-opt-disable", NULL, NULL); // good for debugging kernels

// if you have trouble creating smbols that GDB can pick up to set a breakpoint after kernels are loaded into memory
// this can be used to stop execution to allow you to set a breakpoint in a kernel after kernel symbols are in memory.
#ifdef DEBUG_BREAK_GDB
            __debugbreak();
#endif
#else
            status = clBuildProgram( program, 1, &q_device, "", NULL, NULL);
#endif
            if( status != CL_SUCCESS )
            {
                if( status == CL_BUILD_PROGRAM_FAILURE )
                {
                    size_t buildLogSize = 0;
                    OPENCL_V( clGetProgramBuildInfo( program, q_device, CL_PROGRAM_BUILD_LOG, 0, NULL, &buildLogSize ),
                              _T( "clGetProgramBuildInfo failed" ) );

                    std::vector< char > buildLog( buildLogSize );
                    ::memset( &buildLog[ 0 ], 0x0, buildLogSize );

                    OPENCL_V( clGetProgramBuildInfo( program, q_device, CL_PROGRAM_BUILD_LOG, buildLogSize, &buildLog[ 0 ], NULL ),
                              _T( "clGetProgramBuildInfo failed" ) );

                    std::cerr << "\n\t\t\tBUILD LOG\n";
                    std::cerr << "************************************************\n";
                    std::cerr << &buildLog[ 0 ] << std::endl;
                    std::cerr << "************************************************\n";
                }

                OPENCL_V( status, _T( "clBuildProgram failed" ) );
            }

            lookup.setProgram(program, source);
            lookup.populateCache();
        }

        fftRepo.setclProgram( this->getGenerator(), this->getSignatureData(), program, q_device, fftPlan->context );


        // For real transforms we compile either forward or backward kernel
        bool buildFwdKernel = buildForwardKernel();
        bool buildBwdKernel = buildBackwardKernel();

        // get a kernel object handle for a kernel with the given name
        cl_kernel kernel;
        if( buildFwdKernel )
        {
			lockRAII *kernelLock;
            if( fftRepo.getclKernel( program, CLFFT_FORWARD, kernel, kernelLock) == CLFFT_INVALID_KERNEL )
            {
                std::string entryPoint;
                OPENCL_V( fftRepo.getProgramEntryPoint( this->getGenerator(), this->getSignatureData(), CLFFT_FORWARD, entryPoint, q_device, fftPlan->context ), _T( "fftRepo.getProgramEntryPoint failed." ) );

                kernel = clCreateKernel( program, entryPoint.c_str( ), &status );
                OPENCL_V( status, _T( "clCreateKernel failed" ) );

                fftRepo.setclKernel( program, CLFFT_FORWARD, kernel );
            }
        }

        if( buildBwdKernel )
        {
			lockRAII *kernelLock;
            if( fftRepo.getclKernel( program, CLFFT_BACKWARD, kernel, kernelLock ) == CLFFT_INVALID_KERNEL )
            {
                std::string entryPoint;
                OPENCL_V( fftRepo.getProgramEntryPoint( this->getGenerator(), this->getSignatureData(), CLFFT_BACKWARD, entryPoint, q_device, fftPlan->context ), _T( "fftRepo.getProgramEntryPoint failed." ) );

                kernel = clCreateKernel( program, entryPoint.c_str( ), &status );
                OPENCL_V( status, _T( "clCreateKernel failed" ) );

                fftRepo.setclKernel( program, CLFFT_BACKWARD, kernel );
            }
        }
    }

    return	CLFFT_SUCCESS;
}


