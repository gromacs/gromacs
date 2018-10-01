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
#include <list>
#include "generator.stockham.h"
#include "action.h"

FFTGeneratedCopyAction::FFTGeneratedCopyAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTCopyAction(plHandle, plan, queue, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTCopyAction() failed, exit
        fprintf(stderr, "FFTCopyAction() failed!\n");
        return;
    }

    // Initialize the FFTAction::FFTKernelGenKeyParams member
    err = this->initParams();

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedCopyAction::initParams() failed!\n");
        return;
    }

    FFTRepo &fftRepo = FFTRepo::getInstance();

    err = this->generateKernel(fftRepo, queue);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedCopyAction::generateKernel failed\n");
        return;
    }

    err = compileKernels( queue, plHandle, plan);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedCopyAction::compileKernels failed\n");
        return;
    }

    err = CLFFT_SUCCESS;
}

bool FFTGeneratedCopyAction::buildForwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool r2c_transform = (inputLayout == CLFFT_REAL);
	bool h2c = (inputLayout == CLFFT_HERMITIAN_PLANAR) || (inputLayout == CLFFT_HERMITIAN_INTERLEAVED);
    bool c2h = (outputLayout == CLFFT_HERMITIAN_PLANAR) || (outputLayout == CLFFT_HERMITIAN_INTERLEAVED);

    return (r2c_transform || c2h) || (!(h2c || c2h));
}

bool FFTGeneratedCopyAction::buildBackwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool c2r_transform = (outputLayout == CLFFT_REAL);
    bool h2c = (inputLayout == CLFFT_HERMITIAN_PLANAR) || (inputLayout == CLFFT_HERMITIAN_INTERLEAVED);
	bool c2h = (outputLayout == CLFFT_HERMITIAN_PLANAR) || (outputLayout == CLFFT_HERMITIAN_INTERLEAVED);

    return (c2r_transform || h2c) || (!(h2c || c2h));
}

using namespace StockhamGenerator;

namespace CopyGenerator
{
    // Copy kernel
    template <Precision PR>
    class CopyKernel
    {
        size_t N;
		size_t Nt;
		const FFTGeneratedCopyAction::Signature & params;
		bool h2c, c2h;
		bool general;

		inline std::string OffsetCalc(const std::string &off, bool input = true)
		{
			std::string str;

			const size_t *pStride = input ? params.fft_inStride : params.fft_outStride;

			str += "\t"; str += off; str += " = ";
			std::string nextBatch = "batch";
			for(size_t i=(params.fft_DataDim - 1); i>1; i--)
			{
				size_t currentLength = 1;
				for(int j=1; j<i; j++) currentLength *= params.fft_N[j];

				str += "("; str += nextBatch; str += "/"; str += SztToStr(currentLength);
				str += ")*"; str += SztToStr(pStride[i]); str += " + ";

				nextBatch = "(" + nextBatch + "%" + SztToStr(currentLength) + ")";
			}

			str += nextBatch; str += "*"; str += SztToStr(pStride[1]); str += ";\n";

			return str;
		}

    public:
        CopyKernel( const FFTGeneratedCopyAction::Signature &paramsVal) :
					params(paramsVal)

        {
			N = params.fft_N[0];
			Nt = 1 + N/2;

			h2c = (	(params.fft_inputLayout == CLFFT_HERMITIAN_PLANAR) ||
					(params.fft_inputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;
			c2h = (	(params.fft_outputLayout == CLFFT_HERMITIAN_PLANAR) ||
					(params.fft_outputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;

			general = !(h2c || c2h);

			// We only do out-of-place copies at this point
			assert(params.fft_placeness == CLFFT_OUTOFPLACE);
		}

        void GenerateKernel(std::string &str)
		{
			std::string rType  = RegBaseType<PR>(1);
			std::string r2Type  = RegBaseType<PR>(2);

			bool inIlvd; // Input is interleaved format
			bool outIlvd; // Output is interleaved format
			inIlvd  = (	(params.fft_inputLayout == CLFFT_COMPLEX_INTERLEAVED) ||
						(params.fft_inputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;
			outIlvd = (	(params.fft_outputLayout == CLFFT_COMPLEX_INTERLEAVED) ||
						(params.fft_outputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;



			// Pragma
			str += ClPragma<PR>();

			std::string sfx = FloatSuffix<PR>();

			//If pre-callback is set for the plan
			if (params.fft_hasPreCallback && h2c)
			{
				//Insert callback function code at the beginning 
				str += params.fft_preCallback.funcstring;
				str += "\n\n";
			}

			//if postcallback is set
			if (params.fft_hasPostCallback)
			{
				//Insert callback function code at the beginning 
				str += params.fft_postCallback.funcstring;
				str += "\n\n";
			}

			// Copy kernel begin
			str += "__kernel void ";

			// Function name
			if(general)
					str += "copy_general";
			else
			{
				if(h2c)	str += "copy_h2c";
				else	str += "copy_c2h";
			}

			str += "(";

			if(inIlvd)
			{
				str += "__global const "; str += r2Type; str += " * restrict gbIn, ";
			}
			else
			{
				str += "__global const "; str += rType; str += " * restrict gbInRe, ";
				str += "__global const "; str += rType; str += " * restrict gbInIm, ";
			}

			if(outIlvd)
			{
				str += "__global "; str += r2Type; str += " * restrict gbOut";
			}
			else
			{
				str += "__global "; str += rType; str += " * restrict gbOutRe, ";
				str += "__global "; str += rType; str += " * restrict gbOutIm";
			}

			if (params.fft_hasPreCallback && h2c)
			{
				assert(!params.fft_hasPostCallback);

				str += ", __global void* pre_userdata";
				if (params.fft_preCallback.localMemSize > 0)
				{
					str += ", __local void* localmem";
				}
			}

			if (params.fft_hasPostCallback)
			{
				assert(!params.fft_hasPreCallback);

				str += ", __global void* post_userdata";
				if (params.fft_postCallback.localMemSize > 0)
				{
					str += ", __local void* localmem";
				}
			}

			str += ")\n";

			str += "{\n";

			// Initialize
			if(general)
			{
				str += "\tuint me = get_local_id(0);\n\t";
				str += "uint batch = get_group_id(0);\n\t";
			}
			else
			{
				str += "\tuint me = get_global_id(0);\n\t";
			}

			// Declare memory pointers
			str += "\n\t";
			str += "uint iOffset;\n\t";
			str += "uint oOffset;\n\t";

			if (!(params.fft_hasPreCallback && h2c))
			{
				// input
				if(inIlvd)
				{
					str += "__global "; str += r2Type; str += " *lwbIn;\n\t";
				}
				else
				{
					str += "__global "; str += rType; str += " *lwbInRe;\n\t";
					str += "__global "; str += rType; str += " *lwbInIm;\n\t";
				}
			}

			// output
			if(outIlvd)
			{
				if (!params.fft_hasPostCallback)
				{
					str += "__global "; str += r2Type; str += " *lwbOut;\n";
				}
				if(h2c)
				{
					str += "\t";
					str += "__global "; str += r2Type; str += " *lwbOut2;\n\n";
				}
			}
			else
			{
				if (!params.fft_hasPostCallback)
				{
					str += "__global "; str += rType; str += " *lwbOutRe;\n\t";
					str += "__global "; str += rType; str += " *lwbOutIm;\n";
				}
				if(h2c)
				{
					str += "\t";
					str += "__global "; str += rType; str += " *lwbOutRe2;\n\t";
					str += "__global "; str += rType; str += " *lwbOutIm2;\n\n";
				}
			}
			
			// Setup registers
			str += "\t"; str += RegBaseType<PR>(2); str += " R;\n\n";
			
			size_t NtRounded64 = DivRoundingUp<size_t>(Nt,64) * 64;

			if(!general)
			{
				// Setup variables
				str += "\tuint batch, meg, mel, mel2;\n\t";
				str += "batch = me/"; str += SztToStr(NtRounded64); str += ";\n\t";
				str += "meg = me%"; str += SztToStr(NtRounded64); str += ";\n\t";
				str += "mel = me%"; str += SztToStr(Nt); str += ";\n\t";
				str += "mel2 = ("; str += SztToStr(N); str += " - mel)%"; str += SztToStr(N); str += ";\n\n";
			}


			// Setup memory pointers
			str += OffsetCalc("iOffset", true);
			str += OffsetCalc("oOffset", false);

			// offset strings
			std::string inF, inF2, outF, outF2;
			if(general)
			{
				inF = inF2 = outF = outF2 = "";
			}
			else
			{
				inF   = " + (mel*";  inF   += SztToStr(params.fft_inStride[0]);  inF   += ")";
				inF2  = " + (mel2*"; inF2  += SztToStr(params.fft_inStride[0]);  inF2  += ")";
				outF  = " + (mel*";  outF  += SztToStr(params.fft_outStride[0]); outF  += ")";
				outF2 = " + (mel2*"; outF2 += SztToStr(params.fft_outStride[0]); outF2 += ")";
			}

			str += "\n\t";

			if (!(params.fft_hasPreCallback && h2c))
			{
				// inputs
				if(inIlvd)
				{
					str += "lwbIn = gbIn + iOffset"; str += inF; str += ";\n\t";
				}
				else
				{
					str += "lwbInRe = gbInRe + iOffset"; str += inF; str += ";\n\t";
					str += "lwbInIm = gbInIm + iOffset"; str += inF; str += ";\n\t";
				}
			}

			// outputs
			if(outIlvd)
			{
				if (!params.fft_hasPostCallback)
				{
					str += "lwbOut = gbOut + oOffset"; str += outF; str += ";\n";
				}
				if(h2c)
				{
					str += "\t";
					str += "lwbOut2 = gbOut + oOffset"; str += outF2; str += ";\n";
				}
			}
			else
			{
				if (!params.fft_hasPostCallback)
				{
					str += "lwbOutRe = gbOutRe + oOffset"; str += outF; str += ";\n\t";
					str += "lwbOutIm = gbOutIm + oOffset"; str += outF; str += ";\n";
				}
				if(h2c)
				{
					str += "\t";
					str += "lwbOutRe2 = gbOutRe + oOffset"; str += outF2; str += ";\n\t";
					str += "lwbOutIm2 = gbOutIm + oOffset"; str += outF2; str += ";\n";
				}
			}

			str += "\n\t";

			// Do the copy
			if(general)
			{
				str += "for(uint t=0; t<"; str += SztToStr(N/64); str += "; t++)\n\t{\n\t\t";
				
				if(inIlvd)
				{
					str += "R = lwbIn[me + t*64];\n\t\t";
				}
				else
				{
					str += "R.x = lwbInRe[me + t*64];\n\t\t";
					str += "R.y = lwbInIm[me + t*64];\n\t\t";
				}
				
				if(outIlvd)
				{
					str += "lwbOut[me + t*64] = R;\n";
				}
				else
				{
					str += "lwbOutRe[me + t*64] = R.x;\n\t\t";
					str += "lwbOutIm[me + t*64] = R.y;\n";
				}

				str += "\t}\n\n";
			}
			else
			{
				str += "if(meg < "; str += SztToStr(Nt); str += ")\n\t{\n\t";
				if(c2h)
				{	
					if(inIlvd)
					{
						str += "R = lwbIn[0];\n\t";
					}
					else
					{
						str += "R.x = lwbInRe[0];\n\t";
						str += "R.y = lwbInIm[0];\n\t";
					}
				
					if(outIlvd)
					{
						if (params.fft_hasPostCallback)
						{
							str += params.fft_postCallback.funcname; str += "(gbOut, oOffset"; str += outF;
							str += ", post_userdata, R";
							if (params.fft_postCallback.localMemSize > 0)
							{
								str += ", localmem";
							}
							str += ");\n\n";
						}
						else
						{
						str += "lwbOut[0] = R;\n\n";
					}
					}
					else
					{
						if (params.fft_hasPostCallback)
						{
							str += params.fft_postCallback.funcname; str += "(gbOutRe, gbOutIm, oOffset"; str += outF;
							str += ", post_userdata, R.x, R.y";

							if (params.fft_postCallback.localMemSize > 0)
							{
								str += ", localmem";
							}
							str += ");\n\t";
						}
					else
					{
						str += "lwbOutRe[0] = R.x;\n\t";
						str += "lwbOutIm[0] = R.y;\n\t";
					}
				}
				}
				else
				{
					if (params.fft_hasPreCallback)
					{
						if(inIlvd)
						{
							str += "R = "; str += params.fft_preCallback.funcname; str += "( gbIn, (iOffset"; str += inF; str += "), pre_userdata"; 
						}
						else
						{
							str += "R = "; str += params.fft_preCallback.funcname; str += "( gbInRe, gbInIm, (iOffset"; str += inF; str += "), pre_userdata";
						}
						if (params.fft_preCallback.localMemSize > 0)
						{
							str += ", localmem";
						}
						str += ");\n\t\t";
					}
					else
					{
						if(inIlvd)
						{
							str += "R = lwbIn[0];\n\t";
						}
						else
						{
							str += "R.x = lwbInRe[0];\n\t";
							str += "R.y = lwbInIm[0];\n\t";
						}
					}

					if(outIlvd)
					{
						str += "lwbOut[0] = R;\n\t";
						str += "R.y = -R.y;\n\t";
						str += "lwbOut2[0] = R;\n\t";
					}
					else
					{
						str += "lwbOutRe[0] = R.x;\n\t";
						str += "lwbOutIm[0] = R.y;\n\t";
						str += "R.y = -R.y;\n\t";
						str += "lwbOutRe2[0] = R.x;\n\t";
						str += "lwbOutIm2[0] = R.y;\n\t";
					}
				}
				str += "}\n\n";
			}

			str += "}\n";
		}
    };
};


clfftStatus FFTGeneratedCopyAction::initParams ()
{

    //    Query the devices in this context for their local memory sizes
    //    How we generate a kernel depends on the *minimum* LDS size for all devices.
    //
    const FFTEnvelope * pEnvelope = NULL;
    OPENCL_V(this->plan->GetEnvelope (& pEnvelope), _T("GetEnvelope failed"));
    BUG_CHECK (NULL != pEnvelope);

    this->signature.fft_precision    = this->plan->precision;
    this->signature.fft_placeness    = this->plan->placeness;
    this->signature.fft_inputLayout  = this->plan->inputLayout;
	this->signature.fft_MaxWorkGroupSize = this->plan->envelope.limit_WorkGroupSize;

    ARG_CHECK (this->plan->inStride.size() == this->plan->outStride.size())

    this->signature.fft_outputLayout = this->plan->outputLayout;

	this->signature.fft_DataDim = this->plan->length.size() + 1;
	int i = 0;
	for(i = 0; i < (this->signature.fft_DataDim - 1); i++)
	{
        this->signature.fft_N[i]         = this->plan->length[i];
        this->signature.fft_inStride[i]  = this->plan->inStride[i];
        this->signature.fft_outStride[i] = this->plan->outStride[i];

	}
    this->signature.fft_inStride[i]  = this->plan->iDist;
    this->signature.fft_outStride[i] = this->plan->oDist;

    this->signature.fft_fwdScale  = this->plan->forwardScale;
    this->signature.fft_backScale = this->plan->backwardScale;

	//Set callback if specified
	if (this->plan->hasPreCallback)
	{
		assert(!this->plan->hasPostCallback);

		this->signature.fft_hasPreCallback = true;
		this->signature.fft_preCallback = this->plan->preCallback;

		//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by main FFT kernel
		if (this->plan->preCallback.localMemSize > this->plan->envelope.limit_LocalMemSize)
		{
			fprintf(stderr, "Requested local memory size not available\n");
			return CLFFT_INVALID_ARG_VALUE;
		}
	}

	if (this->plan->hasPostCallback)
	{
		assert(!this->plan->hasPreCallback);

		this->signature.fft_hasPostCallback = true;
		this->signature.fft_postCallback = this->plan->postCallbackParam;

		//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by main FFT kernel
		//Assumes copy kernel does not have both pre and post callback 
		if (this->plan->postCallbackParam.localMemSize > this->plan->envelope.limit_LocalMemSize)
		{
			fprintf(stderr, "Requested local memory size not available\n");
			return CLFFT_INVALID_ARG_VALUE;
		}
	}
	this->signature.limit_LocalMemSize = this->plan->envelope.limit_LocalMemSize;

    return CLFFT_SUCCESS;
}

clfftStatus FFTGeneratedCopyAction::getWorkSizes (std::vector<size_t> & globalWS, std::vector<size_t> & localWS)
{
	bool h2c, c2h;
	h2c = (	(this->signature.fft_inputLayout == CLFFT_HERMITIAN_PLANAR) || (this->signature.fft_inputLayout == CLFFT_HERMITIAN_INTERLEAVED) );
	c2h = (	(this->signature.fft_outputLayout == CLFFT_HERMITIAN_PLANAR) || (this->signature.fft_outputLayout == CLFFT_HERMITIAN_INTERLEAVED) );

	bool general = !(h2c || c2h);

	size_t count = this->plan->batchsize;

	switch(this->signature.fft_DataDim)
	{
	case 5: assert(false);
	case 4: count *= this->signature.fft_N[2];
	case 3: count *= this->signature.fft_N[1];
	case 2:
			{
				if(general)
				{
					count *= 64;
				}
				else
				{
					count *= (DivRoundingUp<size_t>((1 + this->signature.fft_N[0]/2), 64) * 64); 
				}
			}
			break;
	case 1: assert(false);
	}

	globalWS.push_back( count );
    localWS.push_back( 64 );

    return    CLFFT_SUCCESS;
}


using namespace CopyGenerator;

clfftStatus FFTGeneratedCopyAction::generateKernel(FFTRepo& fftRepo, const cl_command_queue commQueueFFT )
{

  bool h2c, c2h;
  h2c = ( (this->signature.fft_inputLayout == CLFFT_HERMITIAN_PLANAR) || (this->signature.fft_inputLayout == CLFFT_HERMITIAN_INTERLEAVED) );
  c2h = ( (this->signature.fft_outputLayout == CLFFT_HERMITIAN_PLANAR) || (this->signature.fft_outputLayout == CLFFT_HERMITIAN_INTERLEAVED) );
  
  bool general = !(h2c || c2h);

  std::string programCode;
  Precision pr = (this->signature.fft_precision == CLFFT_SINGLE) ? P_SINGLE : P_DOUBLE;
  switch(pr)
  {
  case P_SINGLE:
    {
      CopyKernel<P_SINGLE> kernel(this->signature);
      kernel.GenerateKernel(programCode);
    } break;
  case P_DOUBLE:
    {
      CopyKernel<P_DOUBLE> kernel(this->signature);
      kernel.GenerateKernel(programCode);
    } break;
  }

	cl_int status = CL_SUCCESS;
	cl_device_id Device = NULL;
	status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_DEVICE, sizeof(cl_device_id), &Device, NULL);
	OPENCL_V( status, _T( "clGetCommandQueueInfo failed" ) );

    cl_context QueueContext = NULL;
    status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_CONTEXT, sizeof(cl_context), &QueueContext, NULL);
    OPENCL_V( status, _T( "clGetCommandQueueInfo failed" ) );

  OPENCL_V( fftRepo.setProgramCode( this->getGenerator(), this->getSignatureData(), programCode, Device, QueueContext ), _T( "fftRepo.setclString() failed!" ) );

  if(general)
  {
  OPENCL_V( fftRepo.setProgramEntryPoints( this->getGenerator(), this->getSignatureData(), "copy_general", "copy_general", Device, QueueContext ), _T( "fftRepo.setProgramEntryPoint() failed!" ) );
  }
  else
  {
  OPENCL_V( fftRepo.setProgramEntryPoints( this->getGenerator(), this->getSignatureData(), "copy_c2h", "copy_h2c", Device, QueueContext ), _T( "fftRepo.setProgramEntryPoint() failed!" ) );
  }

  return CLFFT_SUCCESS;
}
