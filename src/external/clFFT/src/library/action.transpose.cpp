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


// action.transpose.nonsquare.cpp provides the entry points of "baking"
// nonsquare inplace transpose kernels called in plan.cpp.
// the actual kernel string generation is provided by generator.transpose.cpp

#include "stdafx.h"

#include <math.h>
#include <iomanip>
#include "generator.transpose.h"
#include "action.transpose.h"
#include "generator.stockham.h"

#include "action.h"

FFTGeneratedTransposeNonSquareAction::FFTGeneratedTransposeNonSquareAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTTransposeNonSquareAction(plHandle, plan, queue, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTTransposeNonSquareAction() failed, exit
        fprintf(stderr, "FFTTransposeNonSquareAction() failed!\n");
        return;
    }

    // Initialize the FFTAction::FFTKernelGenKeyParams member
    err = this->initParams();

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedTransposeNonSquareAction::initParams() failed!\n");
        return;
    }

    FFTRepo &fftRepo = FFTRepo::getInstance();

    err = this->generateKernel(fftRepo, queue);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedTransposeNonSquareAction::generateKernel failed\n");
        return;
    }

    err = compileKernels(queue, plHandle, plan);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedTransposeNonSquareAction::compileKernels failed\n");
        return;
    }

    err = CLFFT_SUCCESS;
}


bool FFTGeneratedTransposeNonSquareAction::buildForwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool r2c_transform = (inputLayout == CLFFT_REAL);
    bool c2r_transform = (outputLayout == CLFFT_REAL);
    bool real_transform = (r2c_transform || c2r_transform);

    return (!real_transform) || r2c_transform;
}

bool FFTGeneratedTransposeNonSquareAction::buildBackwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool r2c_transform = (inputLayout == CLFFT_REAL);
    bool c2r_transform = (outputLayout == CLFFT_REAL);
    bool real_transform = (r2c_transform || c2r_transform);

    return (!real_transform) || c2r_transform;
}

// These strings represent the names that are used as strKernel parameters
const std::string pmRealIn("pmRealIn");
const std::string pmImagIn("pmImagIn");
const std::string pmRealOut("pmRealOut");
const std::string pmImagOut("pmImagOut");
const std::string pmComplexIn("pmComplexIn");
const std::string pmComplexOut("pmComplexOut");

clfftStatus FFTGeneratedTransposeNonSquareAction::initParams()
{

    this->signature.fft_precision = this->plan->precision;
    this->signature.fft_placeness = this->plan->placeness;
    this->signature.fft_inputLayout = this->plan->inputLayout;
    this->signature.fft_outputLayout = this->plan->outputLayout;
    this->signature.fft_3StepTwiddle = false;
    this->signature.nonSquareKernelType = this->plan->nonSquareKernelType;

    this->signature.fft_realSpecial = this->plan->realSpecial;

    this->signature.transOutHorizontal = this->plan->transOutHorizontal;	// using the twiddle front flag to specify horizontal write
                                                                            // we do this so as to reuse flags in FFTKernelGenKeyParams
                                                                            // and to avoid making a new one 

    ARG_CHECK(this->plan->inStride.size() == this->plan->outStride.size());

    if (CLFFT_INPLACE == this->signature.fft_placeness)
    {
        //	If this is an in-place transform the
        //	input and output layout
        //	*MUST* be the same.
        //
        ARG_CHECK(this->signature.fft_inputLayout == this->signature.fft_outputLayout)

    /*        for (size_t u = this->plan->inStride.size(); u-- > 0; )
            {
                ARG_CHECK(this->plan->inStride[u] == this->plan->outStride[u]);
            }*/
    }

    this->signature.fft_DataDim = this->plan->length.size() + 1;

    int i = 0;
    for (i = 0; i < (this->signature.fft_DataDim - 1); i++)
    {
        this->signature.fft_N[i] = this->plan->length[i];
        this->signature.fft_inStride[i] = this->plan->inStride[i];
        this->signature.fft_outStride[i] = this->plan->outStride[i];

    }
    this->signature.fft_inStride[i] = this->plan->iDist;
    this->signature.fft_outStride[i] = this->plan->oDist;

    if (this->plan->large1D != 0) {
        ARG_CHECK(this->signature.fft_N[0] != 0)
            //ToDo:ENABLE ASSERT
       //     ARG_CHECK((this->plan->large1D % this->signature.fft_N[0]) == 0)
            this->signature.fft_3StepTwiddle = true;
        //ToDo:ENABLE ASSERT
       // ARG_CHECK(this->plan->large1D == (this->signature.fft_N[1] * this->signature.fft_N[0]));
    }

    //	Query the devices in this context for their local memory sizes
    //	How we generate a kernel depends on the *minimum* LDS size for all devices.
    //
    const FFTEnvelope * pEnvelope = NULL;
    OPENCL_V(this->plan->GetEnvelope(&pEnvelope), _T("GetEnvelope failed"));
    BUG_CHECK(NULL != pEnvelope);

    // TODO:  Since I am going with a 2D workgroup size now, I need a better check than this 1D use
    // Check:  CL_DEVICE_MAX_WORK_GROUP_SIZE/CL_KERNEL_WORK_GROUP_SIZE
    // CL_DEVICE_MAX_WORK_ITEM_SIZES
    this->signature.fft_R = 1; // Dont think i'll use
    this->signature.fft_SIMD = pEnvelope->limit_WorkGroupSize; // Use devices maximum workgroup size

                                                               //Set callback if specified
    if (this->plan->hasPreCallback)
    {
        this->signature.fft_hasPreCallback = true;
        this->signature.fft_preCallback = this->plan->preCallback;
    }
	if (this->plan->hasPostCallback)
	{
		this->signature.fft_hasPostCallback = true;
		this->signature.fft_postCallback = this->plan->postCallbackParam;
	}
	this->signature.limit_LocalMemSize = this->plan->envelope.limit_LocalMemSize;

	this->signature.transposeMiniBatchSize = this->plan->transposeMiniBatchSize;
	this->signature.nonSquareKernelOrder = this->plan->nonSquareKernelOrder;
	this->signature.transposeBatchSize = this->plan->batchsize;

    return CLFFT_SUCCESS;
}


static const size_t lwSize = 256;
static const size_t reShapeFactor = 2;


//	OpenCL does not take unicode strings as input, so this routine returns only ASCII strings
//	Feed this generator the FFTPlan, and it returns the generated program as a string
clfftStatus FFTGeneratedTransposeNonSquareAction::generateKernel(FFTRepo& fftRepo, const cl_command_queue commQueueFFT)
{


    std::string programCode;
	std::string kernelFuncName;//applied to swap kernel for now
    if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING)
    {
		//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by transpose kernel
		if (this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0)
		{
			assert(!this->signature.fft_hasPostCallback);

			bool validLDSSize = false;
			size_t requestedCallbackLDS = 0;

			requestedCallbackLDS = this->signature.fft_preCallback.localMemSize;
			
			validLDSSize = ((2 * this->plan->ElementSize() * 16 * reShapeFactor * 16 * reShapeFactor) + requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;
		
			if(!validLDSSize)
			{
				fprintf(stderr, "Requested local memory size not available\n");
				return CLFFT_INVALID_ARG_VALUE;
			}
		}
        OPENCL_V(clfft_transpose_generator::genTransposeKernelLeadingDimensionBatched(this->signature, programCode, lwSize, reShapeFactor), _T("genTransposeKernel() failed!"));
    }
	else if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
	{
		//pre call back check
		//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by transpose kernel
		if (this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0)
		{
			assert(!this->signature.fft_hasPostCallback);

			bool validLDSSize = false;
			size_t requestedCallbackLDS = 0;

			requestedCallbackLDS = this->signature.fft_preCallback.localMemSize;

			validLDSSize = ((2 * this->plan->ElementSize() * 16 * reShapeFactor * 16 * reShapeFactor) + requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;

			if (!validLDSSize)
			{
				fprintf(stderr, "Requested local memory size not available\n");
				return CLFFT_INVALID_ARG_VALUE;
			}
		}
		OPENCL_V(clfft_transpose_generator::genTransposeKernelBatched(this->signature, programCode, lwSize, reShapeFactor), _T("genTransposeKernel() failed!"));
	}
    else
    {
		//pre-callback is possible in swap kernel now
		if (this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0)
		{
			assert(!this->signature.fft_hasPostCallback);

			bool validLDSSize = false;
			size_t requestedCallbackLDS = 0;

			requestedCallbackLDS = this->signature.fft_preCallback.localMemSize;
			//LDS usage of swap lines is exactly 2 lines
			size_t lineSize = (this->signature.fft_N[0]) < (this->signature.fft_N[1]) ? this->signature.fft_N[0] : this->signature.fft_N[1];
			validLDSSize = ((2 * this->plan->ElementSize() * lineSize) + requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;

			if (!validLDSSize)
			{
				fprintf(stderr, "Requested local memory size not available\n");
				return CLFFT_INVALID_ARG_VALUE;
			}
		}
		//here we should decide generate what kind of swap kernel. 1:2 and 1:3 probably need different swap kernels
		/*
		if (this->signature.fft_N[0] == 2 * this->signature.fft_N[1] || 2 * this->signature.fft_N[0] == this->signature.fft_N[1])
		{
			OPENCL_V(clfft_transpose_generator::genSwapKernel(this->signature, programCode, kernelFuncName, lwSize, reShapeFactor), _T("genSwapKernel() failed!"));
		}
		else
		{
			OPENCL_V(clfft_transpose_generator::genSwapKernelGeneral(this->signature, programCode, kernelFuncName, lwSize, reShapeFactor), _T("genSwapKernel() failed!"));
		}
		*/
		//general swap kernel takes care of all ratio
		OPENCL_V(clfft_transpose_generator::genSwapKernelGeneral(this->signature, programCode, kernelFuncName, lwSize, reShapeFactor), _T("genSwapKernel() failed!"));
    }
	//std::cout << programCode << std::endl;
    cl_int status = CL_SUCCESS;
    cl_device_id Device = NULL;
    status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_DEVICE, sizeof(cl_device_id), &Device, NULL);
    OPENCL_V(status, _T("clGetCommandQueueInfo failed"));

    cl_context QueueContext = NULL;
    status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_CONTEXT, sizeof(cl_context), &QueueContext, NULL);
    OPENCL_V(status, _T("clGetCommandQueueInfo failed"));


    OPENCL_V(fftRepo.setProgramCode(Transpose_NONSQUARE, this->getSignatureData(), programCode, Device, QueueContext), _T("fftRepo.setclString() failed!"));
    if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING)
    {
        // Note:  See genFunctionPrototype( )
        if (this->signature.fft_3StepTwiddle)
        {
            OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), "transpose_nonsquare_tw_fwd", "transpose_nonsquare_tw_back", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
        }
        else
        {
            OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), "transpose_nonsquare", "transpose_nonsquare", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
        }
    }
	else if(this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
	{
        //for non square we do twiddling in swap kernel
        /*
		if (this->signature.fft_3StepTwiddle && (this->signature.transposeMiniBatchSize == 1))
		{
			OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), "transpose_square_tw_fwd", "transpose_square_tw_back", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
		}
		else
		{
			OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), "transpose_square", "transpose_square", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
		}
        */
        OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), "transpose_square", "transpose_square", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
	}
    else
    {
        if (this->signature.fft_3StepTwiddle)//if miniBatchSize > 1 twiddling is done in swap kernel
        {
            std::string kernelFwdFuncName = kernelFuncName + "_tw_fwd";
            std::string kernelBwdFuncName = kernelFuncName + "_tw_back";
            OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), kernelFwdFuncName.c_str(), kernelBwdFuncName.c_str(), Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
        }
        else
            OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_NONSQUARE, this->getSignatureData(), kernelFuncName.c_str(), kernelFuncName.c_str(), Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
    }
    return CLFFT_SUCCESS;
}


clfftStatus FFTGeneratedTransposeNonSquareAction::getWorkSizes(std::vector< size_t >& globalWS, std::vector< size_t >& localWS)
{

    size_t wg_slice;
    size_t smaller_dim = (this->signature.fft_N[0] < this->signature.fft_N[1]) ? this->signature.fft_N[0] : this->signature.fft_N[1];
	size_t bigger_dim = (this->signature.fft_N[0] >= this->signature.fft_N[1]) ? this->signature.fft_N[0] : this->signature.fft_N[1];
	size_t dim_ratio = bigger_dim / smaller_dim;
    size_t global_item_size;

    if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING)
    {
        if (smaller_dim % (16 * reShapeFactor) == 0)
            wg_slice = smaller_dim / 16 / reShapeFactor;
        else
            wg_slice = (smaller_dim / (16 * reShapeFactor)) + 1;

        global_item_size = wg_slice*(wg_slice + 1) / 2 * 16 * 16 * this->plan->batchsize;

        for (int i = 2; i < this->signature.fft_DataDim - 1; i++)
        {
            global_item_size *= this->signature.fft_N[i];
        }

        /*Push the data required for the transpose kernels*/
        globalWS.clear();
		if(this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED_LEADING)
			globalWS.push_back(global_item_size * dim_ratio);
		else if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
			globalWS.push_back(global_item_size);


        localWS.clear();
        localWS.push_back(lwSize);
    }
	else if (this->signature.nonSquareKernelType == NON_SQUARE_TRANS_TRANSPOSE_BATCHED)
	{
		if (smaller_dim % (16 * reShapeFactor) == 0)
			wg_slice = smaller_dim / 16 / reShapeFactor;
		else
			wg_slice = (smaller_dim / (16 * reShapeFactor)) + 1;

		global_item_size = wg_slice*(wg_slice + 1) / 2 * 16 * 16 * this->plan->batchsize;

		for (int i = 2; i < this->plan->length.size(); i++)
		{
			global_item_size *= this->plan->length[i];
		}

		/*Push the data required for the transpose kernels*/
		globalWS.clear();
		globalWS.push_back(global_item_size);


		localWS.clear();
		localWS.push_back(lwSize);
	}
    else
    {
        /*Now calculate the data for the swap kernels */
		// general swap kernel takes care of all ratio. need clean up here
		if(dim_ratio == 2 && 0){
			//1:2 ratio
			size_t input_elm_size_in_bytes;
			switch (this->signature.fft_precision)
			{
			case CLFFT_SINGLE:
			case CLFFT_SINGLE_FAST:
				input_elm_size_in_bytes = 4;
				break;
			case CLFFT_DOUBLE:
			case CLFFT_DOUBLE_FAST:
				input_elm_size_in_bytes = 8;
				break;
			default:
				return CLFFT_TRANSPOSED_NOTIMPLEMENTED;
			}

			switch (this->signature.fft_outputLayout)
			{
			case CLFFT_COMPLEX_INTERLEAVED:
			case CLFFT_COMPLEX_PLANAR:
				input_elm_size_in_bytes *= 2;
				break;
			case CLFFT_REAL:
				break;
			default:
				return CLFFT_TRANSPOSED_NOTIMPLEMENTED;
			}
			size_t max_elements_loaded = AVAIL_MEM_SIZE / input_elm_size_in_bytes;
			size_t num_elements_loaded;
			size_t local_work_size_swap, num_grps_pro_row;

			if ((max_elements_loaded >> 1) > smaller_dim)
			{
				local_work_size_swap = (smaller_dim < 256) ? smaller_dim : 256;
				num_elements_loaded = smaller_dim;
				num_grps_pro_row = 1;
			}
			else
			{
				num_grps_pro_row = (smaller_dim << 1) / max_elements_loaded;
				num_elements_loaded = max_elements_loaded >> 1;
				local_work_size_swap = (num_elements_loaded < 256) ? num_elements_loaded : 256;
			}
			size_t num_reduced_row;
			size_t num_reduced_col;

			if (this->signature.fft_N[1] == smaller_dim)
			{
				num_reduced_row = smaller_dim;
				num_reduced_col = 2;
			}
			else
			{
				num_reduced_row = 2;
				num_reduced_col = smaller_dim;
			}

			size_t *cycle_map = new size_t[num_reduced_row * num_reduced_col * 2];
			/* The memory required by cycle_map cannot exceed 2 times row*col by design*/
			clfft_transpose_generator::get_cycles(cycle_map, num_reduced_row, num_reduced_col);

			global_item_size = local_work_size_swap * num_grps_pro_row * cycle_map[0] * this->plan->batchsize;

			for (int i = 2; i < this->signature.fft_DataDim - 1; i++)
			{
				global_item_size *= this->signature.fft_N[i];
			}
			delete[] cycle_map;

			globalWS.push_back(global_item_size);
			localWS.push_back(local_work_size_swap);
		}
		else
		{
			//if (dim_ratio == 2 || dim_ratio == 3 || dim_ratio == 5 || dim_ratio == 10)
			if (dim_ratio % 2 == 0 || dim_ratio % 3 == 0 || dim_ratio % 5 == 0 || dim_ratio % 10 == 0)
			{
				size_t local_work_size_swap = 256;
				std::vector<std::vector<size_t> > permutationTable;
				clfft_transpose_generator::permutation_calculation(dim_ratio, smaller_dim, permutationTable);
				size_t global_item_size;
				if(this->plan->large1D && (dim_ratio > 1))
					global_item_size = (permutationTable.size() + 2) * local_work_size_swap * this->plan->batchsize;
				else
					global_item_size = (permutationTable.size() + 2) * local_work_size_swap * this->plan->batchsize;
				//for (int i = 2; i < this->plan->length.size(); i++)
				//	global_item_size *= this->plan->length[i];
				size_t LDS_per_WG = smaller_dim;
				while (LDS_per_WG > 1024)//avoiding using too much lds memory. the biggest LDS memory we will allocate would be 1024*sizeof(float2/double2)*2
				{
					if (LDS_per_WG % 2 == 0)
					{
						LDS_per_WG /= 2;
						continue;
					}
					if (LDS_per_WG % 3 == 0)
					{
						LDS_per_WG /= 3;
						continue;
					}
					if (LDS_per_WG % 5 == 0)
					{
						LDS_per_WG /= 5;
						continue;
					}
					return CLFFT_NOTIMPLEMENTED;
				}

				size_t WG_per_line = smaller_dim / LDS_per_WG;
				global_item_size *= WG_per_line;
				globalWS.push_back(global_item_size);
				localWS.push_back(local_work_size_swap);
			}
			else
				return CLFFT_NOTIMPLEMENTED;
		}
    }
    return CLFFT_SUCCESS;
}

FFTGeneratedTransposeSquareAction::FFTGeneratedTransposeSquareAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
	: FFTTransposeSquareAction(plHandle, plan, queue, err)
{
	if (err != CLFFT_SUCCESS)
	{
		// FFTTransposeSquareAction() failed, exit
		fprintf(stderr, "FFTTransposeSquareAction() failed!\n");
		return;
	}

	// Initialize the FFTAction::FFTKernelGenKeyParams member
	err = this->initParams();

	if (err != CLFFT_SUCCESS)
	{
		fprintf(stderr, "FFTGeneratedTransposeSquareAction::initParams() failed!\n");
		return;
	}

	FFTRepo &fftRepo = FFTRepo::getInstance();

	err = this->generateKernel(fftRepo, queue);

	if (err != CLFFT_SUCCESS)
	{
		fprintf(stderr, "FFTGeneratedTransposeSquareAction::generateKernel failed\n");
		return;
	}

	err = compileKernels(queue, plHandle, plan);

	if (err != CLFFT_SUCCESS)
	{
		fprintf(stderr, "FFTGeneratedTransposeSquareAction::compileKernels failed\n");
		return;
	}

	err = CLFFT_SUCCESS;
}


bool FFTGeneratedTransposeSquareAction::buildForwardKernel()
{
	clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
	clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

	bool r2c_transform = (inputLayout == CLFFT_REAL);
	bool c2r_transform = (outputLayout == CLFFT_REAL);
	bool real_transform = (r2c_transform || c2r_transform);

	return (!real_transform) || r2c_transform;
}

bool FFTGeneratedTransposeSquareAction::buildBackwardKernel()
{
	clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
	clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

	bool r2c_transform = (inputLayout == CLFFT_REAL);
	bool c2r_transform = (outputLayout == CLFFT_REAL);
	bool real_transform = (r2c_transform || c2r_transform);

	return (!real_transform) || c2r_transform;
}

/*sqaure action*/
clfftStatus FFTGeneratedTransposeSquareAction::initParams()
{

	this->signature.fft_precision = this->plan->precision;
	this->signature.fft_placeness = this->plan->placeness;
	this->signature.fft_inputLayout = this->plan->inputLayout;
	this->signature.fft_outputLayout = this->plan->outputLayout;
	this->signature.fft_3StepTwiddle = false;

	this->signature.fft_realSpecial = this->plan->realSpecial;

	this->signature.transOutHorizontal = this->plan->transOutHorizontal;	// using the twiddle front flag to specify horizontal write
																			// we do this so as to reuse flags in FFTKernelGenKeyParams
																			// and to avoid making a new one 

	ARG_CHECK(this->plan->inStride.size() == this->plan->outStride.size());

	if (CLFFT_INPLACE == this->signature.fft_placeness)
	{
		//	If this is an in-place transform the
		//	input and output layout, dimensions and strides
		//	*MUST* be the same.
		//
		ARG_CHECK(this->signature.fft_inputLayout == this->signature.fft_outputLayout)

			for (size_t u = this->plan->inStride.size(); u-- > 0; )
			{
				ARG_CHECK(this->plan->inStride[u] == this->plan->outStride[u]);
			}
	}

	this->signature.fft_DataDim = this->plan->length.size() + 1;
	int i = 0;
	for (i = 0; i < (this->signature.fft_DataDim - 1); i++)
	{
		this->signature.fft_N[i] = this->plan->length[i];
		this->signature.fft_inStride[i] = this->plan->inStride[i];
		this->signature.fft_outStride[i] = this->plan->outStride[i];

	}
	this->signature.fft_inStride[i] = this->plan->iDist;
	this->signature.fft_outStride[i] = this->plan->oDist;

	if (this->plan->large1D != 0) {
		ARG_CHECK(this->signature.fft_N[0] != 0)
			ARG_CHECK((this->plan->large1D % this->signature.fft_N[0]) == 0)
			this->signature.fft_3StepTwiddle = true;
		ARG_CHECK(this->plan->large1D == (this->signature.fft_N[1] * this->signature.fft_N[0]));
	}

	//	Query the devices in this context for their local memory sizes
	//	How we generate a kernel depends on the *minimum* LDS size for all devices.
	//
	const FFTEnvelope * pEnvelope = NULL;
	OPENCL_V(this->plan->GetEnvelope(&pEnvelope), _T("GetEnvelope failed"));
	BUG_CHECK(NULL != pEnvelope);

	// TODO:  Since I am going with a 2D workgroup size now, I need a better check than this 1D use
	// Check:  CL_DEVICE_MAX_WORK_GROUP_SIZE/CL_KERNEL_WORK_GROUP_SIZE
	// CL_DEVICE_MAX_WORK_ITEM_SIZES
	this->signature.fft_R = 1; // Dont think i'll use
	this->signature.fft_SIMD = pEnvelope->limit_WorkGroupSize; // Use devices maximum workgroup size

															   //Set callback if specified
	if (this->plan->hasPreCallback)
	{
		this->signature.fft_hasPreCallback = true;
		this->signature.fft_preCallback = this->plan->preCallback;
	}
	if (this->plan->hasPostCallback)
	{
		this->signature.fft_hasPostCallback = true;
		this->signature.fft_postCallback = this->plan->postCallbackParam;
	}
	this->signature.limit_LocalMemSize = this->plan->envelope.limit_LocalMemSize;

	this->signature.transposeMiniBatchSize = this->plan->transposeMiniBatchSize;
	this->signature.transposeBatchSize = this->plan->batchsize;

	return CLFFT_SUCCESS;
}


//	OpenCL does not take unicode strings as input, so this routine returns only ASCII strings
//	Feed this generator the FFTPlan, and it returns the generated program as a string
clfftStatus FFTGeneratedTransposeSquareAction::generateKernel(FFTRepo& fftRepo, const cl_command_queue commQueueFFT)
{
	//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by main FFT kernel
	if ((this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0) ||
		(this->signature.fft_hasPostCallback && this->signature.fft_postCallback.localMemSize > 0))
	{
		assert(!(this->signature.fft_hasPreCallback && this->signature.fft_hasPostCallback));

		bool validLDSSize = false;
		size_t requestedCallbackLDS = 0;

		if (this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0)
			requestedCallbackLDS = this->signature.fft_preCallback.localMemSize;
		else if (this->signature.fft_hasPostCallback && this->signature.fft_postCallback.localMemSize > 0)
			requestedCallbackLDS = this->signature.fft_postCallback.localMemSize;

		validLDSSize = ((2 * this->plan->ElementSize() * 16 * reShapeFactor * 16 * reShapeFactor) + requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;

		if (!validLDSSize)
		{
			fprintf(stderr, "Requested local memory size not available\n");
			return CLFFT_INVALID_ARG_VALUE;
		}
	}

	std::string programCode;
	OPENCL_V(clfft_transpose_generator::genTransposeKernelBatched(this->signature, programCode, lwSize, reShapeFactor), _T("GenerateTransposeKernel() failed!"));

	cl_int status = CL_SUCCESS;
	cl_device_id Device = NULL;
	status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_DEVICE, sizeof(cl_device_id), &Device, NULL);
	OPENCL_V(status, _T("clGetCommandQueueInfo failed"));

	cl_context QueueContext = NULL;
	status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_CONTEXT, sizeof(cl_context), &QueueContext, NULL);
	OPENCL_V(status, _T("clGetCommandQueueInfo failed"));


	OPENCL_V(fftRepo.setProgramCode(Transpose_SQUARE, this->getSignatureData(), programCode, Device, QueueContext), _T("fftRepo.setclString() failed!"));

	// Note:  See genFunctionPrototype( )
	if (this->signature.fft_3StepTwiddle)
	{
		OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_SQUARE, this->getSignatureData(), "transpose_square_tw_fwd", "transpose_square_tw_back", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
	}
	else
	{
		OPENCL_V(fftRepo.setProgramEntryPoints(Transpose_SQUARE, this->getSignatureData(), "transpose_square", "transpose_square", Device, QueueContext), _T("fftRepo.setProgramEntryPoint() failed!"));
	}

	return CLFFT_SUCCESS;
}


clfftStatus FFTGeneratedTransposeSquareAction::getWorkSizes(std::vector< size_t >& globalWS, std::vector< size_t >& localWS)
{

	size_t wg_slice;
	if (this->signature.fft_N[0] % (16 * reShapeFactor) == 0)
		wg_slice = this->signature.fft_N[0] / 16 / reShapeFactor;
	else
		wg_slice = (this->signature.fft_N[0] / (16 * reShapeFactor)) + 1;

	size_t global_item_size = wg_slice*(wg_slice + 1) / 2 * 16 * 16 * this->plan->batchsize;

	for (int i = 2; i < this->signature.fft_DataDim - 1; i++)
	{
		global_item_size *= this->signature.fft_N[i];
	}

	globalWS.clear();
	globalWS.push_back(global_item_size);

	localWS.clear();
	localWS.push_back(lwSize);

	return CLFFT_SUCCESS;
}
