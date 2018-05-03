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
#include "generator.stockham.h"
#include <list>
#include "action.h"


FFTGeneratedStockhamAction::FFTGeneratedStockhamAction(clfftPlanHandle plHandle, FFTPlan * plan, cl_command_queue queue, clfftStatus & err)
    : FFTStockhamAction(plHandle, plan, queue, err)
{
    if (err != CLFFT_SUCCESS)
    {
        // FFTAction() failed, exit
        fprintf(stderr, "FFTStockhamAction() failed!\n");
        return;
    }

    // Initialize the FFTAction::FFTKernelGenKeyParams member
    err = this->initParams();

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedStockhamAction::initParams() failed!\n");
        return;
    }
	
    FFTRepo &fftRepo = FFTRepo::getInstance();

    err = this->generateKernel(fftRepo, queue);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedStockhamAction::generateKernel failed\n");
        return;
    }

	err = compileKernels( queue, plHandle, plan);

    if (err != CLFFT_SUCCESS)
    {
        fprintf(stderr, "FFTGeneratedStockhamAction::compileKernels failed\n");
        return;
    }

    err = CLFFT_SUCCESS;
}

bool FFTGeneratedStockhamAction::buildForwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool r2c_transform = (inputLayout == CLFFT_REAL);
    bool c2r_transform = (outputLayout == CLFFT_REAL);
    bool real_transform = (r2c_transform || c2r_transform);

    return (!real_transform) || r2c_transform;
}

bool FFTGeneratedStockhamAction::buildBackwardKernel()
{
    clfftLayout inputLayout = this->getSignatureData()->fft_inputLayout;
    clfftLayout outputLayout = this->getSignatureData()->fft_outputLayout;

    bool r2c_transform = (inputLayout == CLFFT_REAL);
    bool c2r_transform = (outputLayout == CLFFT_REAL);
    bool real_transform = (r2c_transform || c2r_transform);

    return (!real_transform) || c2r_transform;
}


// FFT Stockham Autosort Method
//
//   Each pass does one digit reverse in essence. Hence by the time all passes are done, complete
//   digit reversal is done and output FFT is in correct order. Intermediate FFTs are stored in natural order,
//   which is not the case with basic Cooley-Tukey algorithm. Natural order in intermediate data makes it
//   convenient for stitching together passes with different radices.
//
//  Basic FFT algorithm:
//
//        Pass loop
//        {
//            Outer loop
//            {
//                Inner loop
//                {
//                }
//            }
//        }
//
//  The sweeps of the outer and inner loop resemble matrix indexing, this matrix changes shape with every pass as noted below
//
//   FFT pass diagram (radix 2)
//
//                k            k+R                                    k
//            * * * * * * * * * * * * * * * *                     * * * * * * * *
//            *   |             |           *                     *   |         *
//            *   |             |           *                     *   |         *
//            *   |             |           * LS        -->       *   |         *
//            *   |             |           *                     *   |         *
//            *   |             |           *                     *   |         *
//            * * * * * * * * * * * * * * * *                     *   |         *
//                         RS                                     *   |         * L
//                                                                *   |         *
//                                                                *   |         *
//                                                                *   |         *
//                                                                *   |         *
//                                                                *   |         *
//                                                                *   |         *
//                                                                *   |         *
//                                                                * * * * * * * *
//                                                                       R
//
//
//    With every pass, the matrix doubles in height and halves in length
//
//
//  N = 2^T = Length of FFT
//  q = pass loop index
//  k = outer loop index = (0 ... R-1)
//  j = inner loop index = (0 ... LS-1)
//
//  Tables shows how values change as we go through the passes
//
//    q | LS   |  R   |  L  | RS
//   ___|______|______|_____|___
//    0 |  1   | N/2  |  2  | N
//    1 |  2   | N/4  |  4  | N/2
//    2 |  4   | N/8  |  8  | N/4
//    . |  .   | .    |  .  | .
//  T-1 |  N/2 | 1    |  N  | 2
//
//
//   Data Read Order
//     Radix 2: k*LS + j, (k+R)*LS + j
//     Radix 3: k*LS + j, (k+R)*LS + j, (k+2R)*LS + j
//     Radix 4: k*LS + j, (k+R)*LS + j, (k+2R)*LS + j, (k+3R)*LS + j
//     Radix 5: k*LS + j, (k+R)*LS + j, (k+2R)*LS + j, (k+3R)*LS + j, (k+4R)*LS + j
//
//   Data Write Order
//       Radix 2: k*L + j, k*L + j + LS
//       Radix 3: k*L + j, k*L + j + LS, k*L + j + 2*LS
//       Radix 4: k*L + j, k*L + j + LS, k*L + j + 2*LS, k*L + j + 3*LS
//       Radix 5: k*L + j, k*L + j + LS, k*L + j + 2*LS, k*L + j + 3*LS, k*L + j + 4*LS
//

namespace StockhamGenerator
{
	// Experimnetal Start =========================================
	// Kernel Generator Parameterization ==========================

	// Uncomment this directive to activate parameter reads from file
//#define PARMETERS_TO_BE_READ

	// Parameters to read
	struct ParamRead
	{
		size_t	workGroupSize;
		size_t	numTransformsPerWg;
		std::vector<size_t> radices;
		bool	halfLds;
	};

	// File format

	// WorkGroupSize:
	// TransformsPerWorkGroup:
	// Radices:
	// LdsUse:

	void ReadParameterFile(ParamRead &readParam)
	{
		const char *fileName = "parameters.txt";
		std::ifstream file(fileName);

		if(!file.is_open())
		{
			std::cout << "File: " << fileName << " could not be opened, exiting ...." << std::endl;
			exit(-1);
		}

		std::string strWgs = "WorkGroupSize:";
		std::string strNtw = "TransformsPerWorkGroup:";
		std::string strRad = "Radices:";
		std::string strLds = "LdsUse:";
		std::string numbers = "0123456789";

		std::string line;
		while(std::getline(file, line))
		{

			size_t pos;

			pos = line.find(strWgs);
			if(pos != std::string::npos)
			{
				line.erase(pos, strWgs.length());
				size_t numStart	= line.find_first_of(numbers);
				size_t numEnd	= line.find_first_not_of(numbers, numStart);
				std::string val = line.substr(numStart, numEnd-numStart);
				readParam.workGroupSize = strtol(val.c_str(), NULL, 10);
				continue;
			}

			pos = line.find(strNtw);
			if(pos != std::string::npos)
			{
				line.erase(pos, strNtw.length());
				size_t numStart	= line.find_first_of(numbers);
				size_t numEnd	= line.find_first_not_of(numbers, numStart);
				std::string val = line.substr(numStart, numEnd-numStart);
				readParam.numTransformsPerWg = strtol(val.c_str(), NULL, 10);
				continue;
			}

			pos = line.find(strRad);
			if(pos != std::string::npos)
			{
				line.erase(pos, strRad.length());
				while(std::string::npos != line.find_first_of(numbers))
				{
					size_t numStart	= line.find_first_of(numbers);
					size_t numEnd	= line.find_first_not_of(numbers, numStart);
					std::string val = line.substr(numStart, numEnd-numStart);
					readParam.radices.push_back(strtol(val.c_str(), NULL, 10));
					line.erase(0, numEnd);
				}
				continue;
			}
		}

		//std::cout << std::endl;
		//std::cout << "File Parameters" << std::endl;
		//std::cout << strWgs << " " << readParam.workGroupSize << std::endl;
		//std::cout << strNtw << " " << readParam.numTransformsPerWg << std::endl;
		//std::cout << strRad << " "; for(size_t i=0; i<readParam.radices.size(); i++) std::cout << readParam.radices[i] << " "; std::cout << std::endl;
		//std::cout << std::endl;
	}

	// Uncomment this directive to experiment with kernels read from file
//#define KERNEL_INTERJECT

	void ReadKernelFromFile(std::string &str)
	{
		const char *fileName = "fft_kernel.cl";
		std::ifstream file(fileName);

		if(!file.is_open())
		{
			std::cout << "File: " << fileName << " could not be opened, exiting ...." << std::endl;
			exit(-1);
		}

		str.clear();

		std::string line;
		while(std::getline(file, line))
		{
			str += line;
			str += '\n';
		}
	}

	// Experimental End ===========================================

#define RADIX_TABLE_COMMON 	{     2048,           256,             1,         4,     8, 8, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{      512,            64,             1,         3,     8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{      256,            64,             1,         4,     4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{       64,            64,             4,         3,     4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{       32,            64,            16,         2,     8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{       16,            64,            16,         2,     4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{        4,            64,            32,         2,     2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },	\
							{        2,            64,            64,         1,     2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },



    template <Precision PR>
	class KernelCoreSpecs
	{
		struct SpecRecord
		{
			size_t length;
			size_t workGroupSize;
			size_t numTransforms;
			size_t numPasses;
			size_t radices[12]; // Setting upper limit of number of passes to 12
		};

		typedef typename std::map<size_t, SpecRecord> SpecTable;
		SpecTable specTable;

	public:
		KernelCoreSpecs()
		{
			switch(PR)
			{
			case P_SINGLE:
				{
					SpecRecord specRecord[] = {

					RADIX_TABLE_COMMON

					//  Length, WorkGroupSize, NumTransforms, NumPasses,  Radices
					{     4096,           256,             1,         4,     8, 8, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0 },
					{     1024,           128,             1,         4,     8, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0 },
					{      128,            64,             4,         3,     8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
					{        8,            64,            32,         2,     4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },

					};

					size_t tableLength = sizeof(specRecord)/sizeof(specRecord[0]);
					for(size_t i=0; i<tableLength; i++) specTable[specRecord[i].length] = specRecord[i];

				} break;

			case P_DOUBLE:
				{
					SpecRecord specRecord[] = {

					RADIX_TABLE_COMMON

					//  Length, WorkGroupSize, NumTransforms, NumPasses,  Radices
					{     1024,           128,             1,         4,     8, 8, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0 },
					{      128,            64,             4,         3,     8, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
					{        8,            64,            16,         3,     2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },

					};

					size_t tableLength = sizeof(specRecord)/sizeof(specRecord[0]);
					for(size_t i=0; i<tableLength; i++) specTable[specRecord[i].length] = specRecord[i];
				} break;

			default:
				assert(false);
			}
		}

		void GetRadices(size_t length, size_t &numPasses, const size_t * &pRadices) const
		{
			pRadices = NULL;
			numPasses = 0;

			typename SpecTable::const_iterator it = specTable.find(length);
			if(it != specTable.end())
			{
				pRadices = it->second.radices;
				numPasses = it->second.numPasses;
			}
		}

		void GetWGSAndNT(size_t length, size_t &workGroupSize, size_t &numTransforms) const
		{
			workGroupSize = 0;
			numTransforms = 0;

			typename SpecTable::const_iterator it = specTable.find(length);
			if(it != specTable.end())
			{
				workGroupSize = it->second.workGroupSize;
				numTransforms = it->second.numTransforms;
			}
		}
	};



	// Given the length of 1d fft, this function determines the appropriate work group size
	// and the number of transforms per work group
	// TODO for optimizations - experiment with different possibilities for work group sizes and num transforms for improving performance
	void DetermineSizes(const size_t &MAX_WGS, const size_t &length, size_t &workGroupSize, size_t &numTrans, Precision &pr)
	{
		assert(MAX_WGS >= 64);

		if(length == 1) // special case
		{
			workGroupSize = 64;
			numTrans = 64;
			return;
		}

		size_t baseRadix[] = {13,11,7,5,3,2}; // list only supported primes
		size_t baseRadixSize = sizeof(baseRadix)/sizeof(baseRadix[0]);

		size_t l = length;
		std::map<size_t, size_t> primeFactorsExpanded;
		for(size_t r=0; r<baseRadixSize; r++)
		{
			size_t rad = baseRadix[r];
			size_t e = 1;
			while(!(l%rad))
			{
				l /= rad;
				e *= rad;
			}

			primeFactorsExpanded[rad] = e;
		}

		assert(l == 1); // Makes sure the number is composed of only supported primes

		if		(primeFactorsExpanded[2] == length)	// Length is pure power of 2
		{
			if		(length >= 1024)	{ workGroupSize = (MAX_WGS >= 256) ? 256 : MAX_WGS; numTrans = 1; }
			else if (length == 512)		{ workGroupSize = 64; numTrans = 1; }
			else if	(length >= 16)		{ workGroupSize = 64;  numTrans = 256/length; }
			else						{ workGroupSize = 64;  numTrans = 128/length; }
		}
		else if	(primeFactorsExpanded[3] == length) // Length is pure power of 3
		{
			workGroupSize = (MAX_WGS >= 256) ? 243 : 27;
			numTrans = length >= 3*workGroupSize ? 1 : (3*workGroupSize)/length;
		}
		else if	(primeFactorsExpanded[5] == length) // Length is pure power of 5
		{
			workGroupSize = (MAX_WGS >= 128) ? 125 : 25;
			numTrans = length >= 5*workGroupSize ? 1 : (5*workGroupSize)/length;
		}
		else if	(primeFactorsExpanded[7] == length) // Length is pure power of 7
		{
			workGroupSize = 49;
			numTrans = length >= 7*workGroupSize ? 1 : (7*workGroupSize)/length;
		}
		else if (primeFactorsExpanded[11] == length) // Length is pure power of 11
		{
			workGroupSize = 121;
			numTrans = length >= 11 * workGroupSize ? 1 : (11 * workGroupSize) / length;
		}
		else if (primeFactorsExpanded[13] == length) // Length is pure power of 13
		{
			workGroupSize = 169;
			numTrans = length >= 13 * workGroupSize ? 1 : (13 * workGroupSize) / length;
		}
		else
		{
			size_t leastNumPerWI = 1; // least number of elements in one work item
			size_t maxWorkGroupSize = MAX_WGS; // maximum work group size desired


			if        (primeFactorsExpanded[2] * primeFactorsExpanded[3] == length) { 
				if (length % 12 == 0) { 
					leastNumPerWI = 12; maxWorkGroupSize = 128;
				} else { 
					leastNumPerWI =  6; maxWorkGroupSize = 256;
				}
			} else if (primeFactorsExpanded[2] * primeFactorsExpanded[5] == length) { 
				if (length % 20 == 0) { 
					leastNumPerWI = 20; maxWorkGroupSize = 64;
				} else { 
					leastNumPerWI = 10; maxWorkGroupSize = 128;
				}
			} else if (primeFactorsExpanded[2] * primeFactorsExpanded[7] == length) { 
					leastNumPerWI = 14; maxWorkGroupSize = 64;
			} else if (primeFactorsExpanded[3] * primeFactorsExpanded[5] == length) { 
				    leastNumPerWI = 15; maxWorkGroupSize = 128;
			} else if (primeFactorsExpanded[3] * primeFactorsExpanded[7] == length) { 
				    leastNumPerWI = 21; maxWorkGroupSize = 128;
			} else if (primeFactorsExpanded[5] * primeFactorsExpanded[7] == length) { 
				    leastNumPerWI = 35; maxWorkGroupSize = 64;
			} else if (primeFactorsExpanded[2] * primeFactorsExpanded[3] * primeFactorsExpanded[5] == length) { 
				    leastNumPerWI = 30; maxWorkGroupSize = 64;
			} else if (primeFactorsExpanded[2] * primeFactorsExpanded[3] * primeFactorsExpanded[7] == length) { 
				    leastNumPerWI = 42; maxWorkGroupSize = 60;
			} else if (primeFactorsExpanded[2] * primeFactorsExpanded[5] * primeFactorsExpanded[7] == length) { 
				    leastNumPerWI = 70; maxWorkGroupSize = 36;
			} else if (primeFactorsExpanded[3] * primeFactorsExpanded[5] * primeFactorsExpanded[7] == length) { 
				    leastNumPerWI =105; maxWorkGroupSize = 24;
			}
			else if (primeFactorsExpanded[2] * primeFactorsExpanded[11] == length) {
				leastNumPerWI = 22; maxWorkGroupSize = 128;
			}
			else if (primeFactorsExpanded[2] * primeFactorsExpanded[13] == length) {
				leastNumPerWI = 26; maxWorkGroupSize = 128;
			}
			else {
				    leastNumPerWI =210; maxWorkGroupSize = 12;
			}
			if (pr==P_DOUBLE)
			{
				//leastNumPerWI /= 2; 
				maxWorkGroupSize /= 2;
			}
			

			if (maxWorkGroupSize > MAX_WGS)
				maxWorkGroupSize = MAX_WGS;
			assert (leastNumPerWI > 0 && length % leastNumPerWI == 0);

			for (size_t lnpi = leastNumPerWI; lnpi <= length; lnpi += leastNumPerWI) {
				if (length % lnpi != 0) continue;

				if (length / lnpi <= MAX_WGS) {
					leastNumPerWI = lnpi;
					break;
				}
			}

			numTrans = maxWorkGroupSize / (length / leastNumPerWI);
			numTrans = numTrans < 1 ? 1 : numTrans;
			workGroupSize = numTrans * (length / leastNumPerWI);
		}

		assert(workGroupSize <= MAX_WGS);
	}

	// Twiddle factors table
    class TwiddleTable
    {
        size_t N; // length
		double *wc, *ws; // cosine, sine arrays

	public:
		TwiddleTable(size_t length) : N(length)
		{
			// Allocate memory for the tables
			// We compute twiddle factors in double precision for both P_SINGLE and P_DOUBLE
			wc = new double[N];
			ws = new double[N];
		}

		~TwiddleTable()
		{
			// Free
			delete[] wc;
			delete[] ws;
		}

		template <Precision PR>
		void GenerateTwiddleTable(const std::vector<size_t> &radices, std::string &twStr)
		{
			const double TWO_PI = -6.283185307179586476925286766559;

			// Make sure the radices vector sums up to N
			size_t sz = 1;
			for(std::vector<size_t>::const_iterator i = radices.begin();
				i != radices.end(); i++)
			{
				sz *= (*i);
			}
			assert(sz == N);

			// Generate the table
			size_t L = 1;
			size_t nt = 0;
			for(std::vector<size_t>::const_iterator i = radices.begin();
				i != radices.end(); i++)
			{
				size_t radix = *i;

				L *= radix;

				// Twiddle factors
				for(size_t k=0; k<(L/radix); k++)
				{
					double theta = TWO_PI * ((double)k)/((double)L);

					for(size_t j=1; j<radix; j++)
					{
						double c = cos(((double)j) * theta);
						double s = sin(((double)j) * theta);

						//if (fabs(c) < 1.0E-12)	c = 0.0;
						//if (fabs(s) < 1.0E-12)	s = 0.0;

						wc[nt]   = c;
						ws[nt++] = s;
					}
				}
			}

			std::string sfx = FloatSuffix<PR>();

			// Stringize the table
			std::stringstream ss;
			ss.imbue(std::locale("C"));
			ss.precision(34);
			ss << std::scientific;
			for(size_t i = 0; i < (N-1); i++)
			{
				ss << "("; ss << RegBaseType<PR>(2); ss << ")(";
				ss << wc[i] << sfx << ", ";
				ss << ws[i] << sfx << "),\n";
			}
			twStr += ss.str();
		}
	};

    // A pass inside an FFT kernel
    template <Precision PR>
    class Pass
    {
		size_t position;					// Position in the kernel

		size_t algL;						// 'L' value from fft algorithm
		size_t algLS;						// 'LS' value
		size_t algR;						// 'R' value

		size_t length;						// Length of FFT
        size_t radix;						// Base radix
		size_t cnPerWI;						// Complex numbers per work-item

		size_t workGroupSize;				// size of the workgroup = (length / cnPerWI)
											// this number is essentially number of work-items needed to compute 1 transform
											// this number will be different from the kernel class workGroupSize if there
											// are multiple transforms per workgroup

		size_t numButterfly;				// Number of basic FFT butterflies = (cnPerWI / radix)
		size_t numB1, numB2, numB4;			// number of different types of butterflies

		bool r2c;							// real to complex transform
		bool c2r;							// complex to real transform
		bool rcFull;
		bool rcSimple;

		bool realSpecial;

		bool enableGrouping;				
		bool linearRegs;					// scalar registers (non-vectorized registers) to be used
		bool halfLds;						// only half the LDS of a complex length need to be used
		Pass<PR> *nextPass;

		//callback members
		bool fft_doPreCallback;
		clfftCallbackParam fft_preCallback;

		bool fft_doPostCallback;
		clfftCallbackParam fft_postCallback;

		inline void RegBase(size_t regC, std::string &str) const
		{
			str += "B";
			str += SztToStr(regC);
		}

		inline void RegBaseAndCount(size_t num, std::string &str) const
		{
			str += "C";
			str += SztToStr(num);
		}

		inline void RegBaseAndCountAndPos(const std::string &RealImag, size_t radPos, std::string &str) const
		{
			str += RealImag;
			str += SztToStr(radPos);
		}

		void RegIndex(size_t regC, size_t num, const std::string &RealImag, size_t radPos, std::string &str) const
		{
			RegBase(regC, str);
			RegBaseAndCount(num, str);
			RegBaseAndCountAndPos(RealImag, radPos, str);
		}

		void DeclareRegs(const std::string &regType, size_t regC, size_t numB, std::string &passStr) const
		{
			std::string regBase;
			RegBase(regC, regBase);

			if(linearRegs)
			{
				assert(regC == 1);
				assert(numB == numButterfly);
			}

			for(size_t i=0; i<numB; i++)
			{
				passStr += "\n\t";
				passStr += regType;
				passStr += " ";

				std::string regBaseCount = regBase;
				RegBaseAndCount(i, regBaseCount);

				for(size_t r=0; ; r++)
				{
					if(linearRegs)
					{
						std::string regIndex = "R";
						RegBaseAndCountAndPos("", i*radix + r, regIndex);

						passStr += regIndex;
					}
					else
					{
						std::string regRealIndex(regBaseCount), regImagIndex(regBaseCount);

						RegBaseAndCountAndPos("R", r, regRealIndex); // real
						RegBaseAndCountAndPos("I", r, regImagIndex); // imaginary

						passStr += regRealIndex; passStr += ", ";
						passStr += regImagIndex;
					}

					if(r == radix-1)
					{
						passStr += ";";
						break;
					}
					else
					{
						passStr += ", ";
					}
				}
			}
		}

		inline std::string IterRegArgs() const
		{
			std::string str = "";

			if(linearRegs)
			{
				std::string regType = RegBaseType<PR>(2);

				for(size_t i=0; i<cnPerWI; i++)
				{
					if(i != 0) str += ", ";
					str += regType; str += " *R";
					str += SztToStr(i);
				}
			}

			return str;
		}

#define SR_READ			1
#define SR_TWMUL		2
#define SR_TWMUL_3STEP	3
#define SR_WRITE		4

#define SR_COMP_REAL 0 // real
#define SR_COMP_IMAG 1 // imag
#define SR_COMP_BOTH 2 // real & imag

		// SweepRegs is to iterate through the registers to do the three basic operations:
		// reading, twiddle multiplication, writing
		void SweepRegs(	size_t flag, bool fwd, bool interleaved, size_t stride, size_t component,
						double scale, bool frontTwiddle,
						const std::string &bufferRe, const std::string &bufferIm, const std::string &offset,
						size_t regC, size_t numB, size_t numPrev, std::string &passStr, bool initZero = false, bool isPrecallVector = false, bool oddt = false) const
		{
			assert( (flag == SR_READ )			||
					(flag == SR_TWMUL)			||
					(flag == SR_TWMUL_3STEP)	||
					(flag == SR_WRITE) );

			const std::string twTable = TwTableName();
			const std::string tw3StepFunc = TwTableLargeFunc();

			// component: 0 - real, 1 - imaginary, 2 - both
			size_t cStart, cEnd;
			switch(component)
			{
			case SR_COMP_REAL:	cStart = 0; cEnd = 1; break;
			case SR_COMP_IMAG:	cStart = 1; cEnd = 2; break;
			case SR_COMP_BOTH:	cStart = 0; cEnd = 2; break;
			default:	assert(false);
			}

			// Read/Write logic:
			// The double loop inside pass loop of FFT algorithm is mapped into the
			// workGroupSize work items with each work item handling cnPerWI numbers

			// Read logic:
			// Reads for any pass appear the same with the stockham algorithm when mapped to
			// the work items. The buffer is divided into (L/radix) sized blocks and the
			// values are read in linear order inside each block.

			// Vector reads are possible if we have unit strides
			// since read pattern remains the same for all passes and they are contiguous
			// Writes are not contiguous

			// TODO : twiddle multiplies can be combined with read
			// TODO : twiddle factors can be reordered in the table to do vector reads of them

			// Write logic:
			// outer loop index k and the inner loop index j map to 'me' as follows:
			// In one work-item (1 'me'), there are 'numButterfly' fft butterflies. They
			// are indexed as numButterfly*me + butterflyIndex, where butterflyIndex's range is
			// 0 ... numButterfly-1. The total number of butterflies needed is covered over all
			// the work-items. So essentially the double loop k,j is flattened to fit this linearly
			// increasing 'me'.
			// j = (numButterfly*me + butterflyIndex)%LS
			// k = (numButterfly*me + butterflyIndex)/LS


			std::string twType = RegBaseType<PR>(2);
			std::string rType  = RegBaseType<PR>(1);

			size_t butterflyIndex = numPrev;
			std::string bufOffset;

			std::string regBase;
			RegBase(regC, regBase);

			// special write back to global memory with float4 grouping, writing 2 complex numbers at once
			if( numB && (numB%2 == 0) && (regC == 1) && (stride == 1) && (numButterfly%2 == 0) && (algLS%2 == 0) && (flag == SR_WRITE) &&
				(nextPass == NULL) && interleaved && (component == SR_COMP_BOTH) && linearRegs && enableGrouping && !fft_doPostCallback )
			{
				assert((numButterfly * workGroupSize) == algLS);
				assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data

				passStr += "\n\t";
				passStr += "__global "; passStr += RegBaseType<PR>(4);
				passStr += " *buff4g = (__global "; passStr += RegBaseType<PR>(4); passStr += " *)";
				passStr += bufferRe; passStr += ";\n\t"; // Assuming 'outOffset' is 0, so not adding it here

				for(size_t r=0; r<radix; r++) // setting the radix loop outside to facilitate grouped writing
				{
					butterflyIndex = numPrev;

					for(size_t i=0; i<(numB/2); i++)
					{
						std::string regIndexA = "(*R";
						std::string regIndexB = "(*R";

						RegBaseAndCountAndPos("", (2*i + 0)*radix + r, regIndexA); regIndexA += ")";
						RegBaseAndCountAndPos("", (2*i + 1)*radix + r, regIndexB); regIndexB += ")";

						passStr += "\n\t";
						passStr += "buff4g"; passStr += "[ ";
						passStr += SztToStr(numButterfly/2); passStr += "*me + "; passStr += SztToStr(butterflyIndex);
						passStr += " + ";
						passStr += SztToStr(r*(algLS/2)); passStr += " ]";
						passStr += " = "; passStr += "("; passStr += RegBaseType<PR>(4); passStr += ")(";
						passStr += regIndexA; passStr += ".x, ";
						passStr += regIndexA; passStr += ".y, ";
						passStr += regIndexB; passStr += ".x, ";
						passStr += regIndexB; passStr += ".y) ";
						if(scale != 1.0f) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); }
						passStr += ";";

						butterflyIndex++;
					}
				}

				return;
			}

			size_t hid = 0;
			bool swapElement = false;
			size_t tIter = numB * radix;

			// block to rearrange reads of adjacent memory locations together
			if(linearRegs && (flag == SR_READ))
			{
				for(size_t r=0; r<radix; r++)
				{
					for(size_t i=0; i<numB; i++)
					{
						for(size_t c=cStart; c<cEnd; c++) // component loop: 0 - real, 1 - imaginary
						{
							swapElement = (fft_doPreCallback && c2r && component == SR_COMP_REAL); //reset at start of loop

							std::string tail;
							std::string regIndex;
							std::string regIndexC;
							regIndex = "(*R";
							std::string buffer;

							// Read real & imag at once
							if(interleaved && (component == SR_COMP_BOTH))
							{
								assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data
								buffer = bufferRe;
								RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ")";
								tail = ";";
							}
							else
							{
								if(c == 0)
								{
									RegBaseAndCountAndPos("", i*radix + r, regIndex); 
									
									hid = (i * radix + r) / ( tIter > 1 ? (tIter / 2) : 1 );
									swapElement = swapElement && hid != 0;
									swapElement = (oddt && ((i * radix + r) >= (tIter - 1))) ? false : swapElement;  //for c2r odd size don't swap for last register
									if (swapElement)
									{
										regIndexC = regIndex; regIndexC += ").y";
									}
									
									regIndex += ").x";
									buffer = bufferRe;
									tail = interleaved ? ".x;" : ";";
								}
								else
								{
									RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").y";
									buffer = bufferIm;
									tail = interleaved ? ".y;" : ";";
								}
							}

							//get offset 
							bufOffset.clear();
							bufOffset += offset; bufOffset += " + ( "; bufOffset += SztToStr(numPrev); bufOffset += " + ";
							bufOffset += "me*"; bufOffset += SztToStr(numButterfly); bufOffset += " + ";
							bufOffset += SztToStr(i); bufOffset += " + ";
							bufOffset += SztToStr(r*length/radix); bufOffset += " )*";
							bufOffset += SztToStr(stride);

							//If precallback is set invoke callback function
							//Invoke callback only once in Planar data layout (i.e.c==0)
							if (fft_doPreCallback && c == 0 && component == SR_COMP_BOTH)
							{
								passStr += "\n\t";
								passStr += "retPrecallback = "; passStr += fft_preCallback.funcname; passStr += "("; 
								if(interleaved)
								{
									passStr += buffer; passStr += ", ";
								}
								else
								{
									passStr += bufferRe; passStr += ", "; passStr += bufferIm; passStr += ", ";
								}
								passStr += bufOffset; passStr += ", pre_userdata";
								if (fft_preCallback.localMemSize > 0)
								{
									passStr += ", localmem";
								}
								passStr += ");";
							}

							if (swapElement)
							{
								passStr += "\n\t";
								passStr += regIndexC; passStr += " = "; passStr += regIndex; passStr += ";";
							}

							passStr += "\n\t";
							passStr += regIndex;
							passStr += " = ";

							if (initZero)
							{
								if (interleaved && (component == SR_COMP_BOTH))
									passStr += "(fvect2)(0, 0);";
								else
									passStr += "0;";
							}
							else
							{
								//Use the return value from precallback if set
								if (fft_doPreCallback && (component == SR_COMP_BOTH || r2c))
								{
									if (component == SR_COMP_BOTH)
									{
										passStr += "retPrecallback";
										passStr += interleaved ? tail : (c == 0) ? ".x;" : ".y;";
									}
									else if (r2c)
									{
										passStr += fft_preCallback.funcname; passStr += "("; passStr += buffer; passStr += ", ";
										passStr += bufOffset; passStr += ", pre_userdata";

										if (fft_preCallback.localMemSize > 0)
										{
											passStr += ", localmem";
										}
										passStr += ");";
									}
								}
								else
								{
									passStr += buffer;
									passStr += "["; passStr += bufOffset; passStr += "]"; passStr += tail;
								}
							}

							// Since we read real & imag at once, we break the loop
							if(interleaved && (component == SR_COMP_BOTH) )
								break;
						}
					}
				}
				return;
			}

			// block to rearrange writes of adjacent memory locations together
			if(linearRegs && (flag == SR_WRITE) && (nextPass == NULL))
			{
				for(size_t r=0; r<radix; r++)
				{
					butterflyIndex = numPrev;

					for(size_t i=0; i<numB; i++)
					{
						if(realSpecial && (nextPass == NULL) && (r > (radix/2)))
							break;

						if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i != 0))
							break;

						if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i == 0))
							passStr += "\n\t}\n\tif( rw && !me)\n\t{";

						std::string regIndexC0;
						for(size_t c=cStart; c<cEnd; c++) // component loop: 0 - real, 1 - imaginary
						{
							std::string tail;
							std::string regIndex;
							regIndex = "(*R";
							std::string buffer;

							// Write real & imag at once
							if(interleaved && (component == SR_COMP_BOTH))
							{
								assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data
								buffer = bufferRe;
								RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ")";
								tail = "";
							}
							else
							{
								if(c == 0)
								{
									RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").x";
									buffer = bufferRe;
									tail = interleaved ? ".x" : "";
								}
								else
								{
									RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").y";
									buffer = bufferIm;
									tail = interleaved ? ".y" : "";
								}
							}

							bufOffset.clear();
							bufOffset += offset; bufOffset += " + ( "; 
							if( (numButterfly * workGroupSize) > algLS )
							{
								bufOffset += "(("; bufOffset += SztToStr(numButterfly);
								bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex); bufOffset += ")/";
								bufOffset += SztToStr(algLS); bufOffset += ")*"; bufOffset += SztToStr(algL); bufOffset += " + (";
								bufOffset += SztToStr(numButterfly); bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex);
								bufOffset += ")%"; bufOffset += SztToStr(algLS); bufOffset += " + ";
							}
							else
							{
								bufOffset += SztToStr(numButterfly); bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex);
								bufOffset += " + ";
							}
							bufOffset += SztToStr(r*algLS); bufOffset += " )*"; bufOffset += SztToStr(stride);

							if(scale != 1.0f) { regIndex += " * "; regIndex += FloatToStr(scale); regIndex += FloatSuffix<PR>(); }
							if (c == cStart)	regIndexC0 = regIndex;

							if (fft_doPostCallback && !r2c)
							{
								if (interleaved || c == (cEnd - 1))
								{
							passStr += "\n\t";
									passStr += fft_postCallback.funcname; passStr += "(";

									if (interleaved || (c2r && bufferRe.compare(bufferIm) == 0))
							{
										passStr += buffer;
							}
							else
							{
										passStr += bufferRe; passStr += ", "; passStr += bufferIm;
							}
									passStr += ", ";
									passStr += bufOffset; passStr += ", post_userdata, ("; passStr += regIndexC0; passStr += ")";
									if (!(interleaved || (c2r && bufferRe.compare(bufferIm) == 0))) { passStr += ", ("; passStr += regIndex; passStr += ")"; }

									if (fft_postCallback.localMemSize > 0)
									{
										passStr += ", post_localmem";
									}
									passStr += ");";
								}
							}
							else
							{	
								passStr += "\n\t";
								passStr += buffer; passStr += "["; passStr += bufOffset; passStr += "]";
								passStr += tail; passStr += " = "; passStr += regIndex; passStr += ";";
							}

							// Since we write real & imag at once, we break the loop
							if(interleaved && (component == SR_COMP_BOTH))
								break;
						}

						if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i == 0))
							passStr += "\n\t}\n\tif(rw)\n\t{";

						butterflyIndex++;
					}
				}

				return;
			}
			
			
			for(size_t i=0; i<numB; i++)
			{
				std::string regBaseCount = regBase;
				RegBaseAndCount(i, regBaseCount);

				if(flag == SR_READ) // read operation
				{
					// the 'r' (radix index) loop is placed outer to the
					// 'v' (vector index) loop to make possible vectorized reads

					for(size_t r=0; r<radix; r++)
					{
						for(size_t c=cStart; c<cEnd; c++) // component loop: 0 - real, 1 - imaginary
						{
							std::string tail;
							std::string regIndex;
							std::string regIndexC;
							regIndex = linearRegs ? "(*R" : regBaseCount;
							std::string buffer;

							// Read real & imag at once
							if(interleaved && (component == SR_COMP_BOTH) && linearRegs)
							{
								assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data
								buffer = bufferRe;
								RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ")";
								tail = ";";
							}
							else
							{
								if(c == 0)
								{
									if(linearRegs) 
									{ 
										RegBaseAndCountAndPos("", i*radix + r, regIndex); 
									
										hid = (i * radix + r) / (numB * radix / 2);
										if (fft_doPreCallback && c2r && component == SR_COMP_REAL && hid != 0)
										{
											regIndexC = regIndex; regIndexC += ").y";
										}
										regIndex += ").x"; 
									}
									else		   { RegBaseAndCountAndPos("R", r, regIndex); }
									buffer = bufferRe;
									tail = interleaved ? ".x;" : ";";
								}
								else
								{
									if(linearRegs) { RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").y"; }
									else		   { RegBaseAndCountAndPos("I", r, regIndex); }
									buffer = bufferIm;
									tail = interleaved ? ".y;" : ";";
								}
							}

							for(size_t v=0; v<regC; v++) // TODO: vectorize the reads; instead of reading individually for consecutive reads of vector elements
							{
								std::string regIndexSub(regIndex);
								if(regC != 1)
								{
									regIndexSub += ".s";
									regIndexSub += SztToStr(v);
								}

								//get offset 
								bufOffset.clear();
								bufOffset += offset; bufOffset += " + ( "; bufOffset += SztToStr(numPrev); bufOffset += " + ";
								bufOffset += "me*"; bufOffset += SztToStr(numButterfly); bufOffset += " + ";
								bufOffset += SztToStr(i*regC + v); bufOffset += " + ";
								bufOffset += SztToStr(r*length/radix); bufOffset += " )*";
								bufOffset += SztToStr(stride);

								//If precallback is set invoke callback function
								//Invoke callback only once in Planar data layout (i.e.c==0)
								if (fft_doPreCallback && c == 0 && component == SR_COMP_BOTH)
								{
									passStr += "\n\t";
									passStr += "retPrecallback"; 
									
									if (isPrecallVector)
									{
										passStr += "["; passStr += SztToStr(v); passStr += "]"; 
									}

									passStr += " = "; passStr += fft_preCallback.funcname; passStr += "("; 
									if(interleaved)
									{
										passStr += buffer; passStr += ", ";
									}
									else
									{
										passStr += bufferRe; passStr += ", "; passStr += bufferIm; passStr += ", ";
									}
									passStr += bufOffset; passStr += ", pre_userdata";
									if (fft_preCallback.localMemSize > 0)
									{
										passStr += ", localmem";
									}
									passStr += ");";
								}

								if (fft_doPreCallback && c2r && component == SR_COMP_REAL && hid != 0)
								{
									passStr += "\n\t";
									passStr += regIndexC; passStr += " = "; passStr += regIndexSub; passStr += ";";
								}

								passStr += "\n\t";
								passStr += regIndexSub;
								passStr += " = "; 

								//Use the return value from precallback if set
								if (fft_doPreCallback && (component == SR_COMP_BOTH || r2c))
								{
									if (component == SR_COMP_BOTH)
									{
										passStr += "retPrecallback"; 
									
										if (isPrecallVector)
										{
											passStr += "["; passStr += SztToStr(v); passStr += "]"; 
										} 
										passStr += interleaved ? tail : (c == 0) ? ".x;" : ".y;";
									}
									else if (r2c)
									{
										passStr += fft_preCallback.funcname; passStr += "("; passStr += buffer; passStr += ", ";
										passStr += bufOffset; passStr += ", pre_userdata";

										if (fft_preCallback.localMemSize > 0)
										{
											passStr += ", localmem";
										}
										passStr += ");";
									}
								}
								else
								{
									passStr += buffer;
									passStr += "["; passStr += bufOffset; passStr += "]"; passStr += tail;
								}
							}

							// Since we read real & imag at once, we break the loop
							if(interleaved && (component == SR_COMP_BOTH) && linearRegs)
								break;
						}
					}
				}
				else if( (flag == SR_TWMUL) || (flag == SR_TWMUL_3STEP) ) // twiddle multiplies and writes require that 'r' loop be innermost
				{
					for(size_t v=0; v<regC; v++)
					{
						for(size_t r=0; r<radix; r++)
						{

							std::string regRealIndex, regImagIndex;
							regRealIndex = linearRegs ? "(*R" : regBaseCount;
							regImagIndex = linearRegs ? "(*R" : regBaseCount;

							if(linearRegs)
							{
								RegBaseAndCountAndPos("", i*radix + r, regRealIndex); regRealIndex += ").x";
								RegBaseAndCountAndPos("", i*radix + r, regImagIndex); regImagIndex += ").y";
							}
							else
							{
								RegBaseAndCountAndPos("R", r, regRealIndex);
								RegBaseAndCountAndPos("I", r, regImagIndex);
							}

							if(regC != 1)
							{
								regRealIndex += ".s"; regRealIndex += SztToStr(v);
								regImagIndex += ".s"; regImagIndex += SztToStr(v);
							}


							if(flag == SR_TWMUL) // twiddle multiply operation
							{
								if(r == 0) // no twiddle muls needed
									continue;

								passStr += "\n\t{\n\t\t"; passStr += twType; passStr += " W = ";
								passStr += twTable; passStr += "["; passStr += SztToStr(algLS-1); passStr += " + ";
								passStr += SztToStr(radix-1); passStr += "*(("; passStr += SztToStr(numButterfly);
								passStr += "*me + "; passStr += SztToStr(butterflyIndex); passStr += ")%";
								passStr += SztToStr(algLS); passStr += ") + "; passStr += SztToStr(r-1);
								passStr += "];\n\t\t";
							}
							else	// 3-step twiddle
							{
								passStr += "\n\t{\n\t\t"; passStr += twType; passStr += " W = ";
								passStr += tw3StepFunc; passStr += "( ";

								if(frontTwiddle)
								{
									assert(linearRegs);
									passStr += "("; passStr += "me*"; passStr += SztToStr(numButterfly);
									passStr += " + "; passStr += SztToStr(i); passStr += " + ";
									passStr += SztToStr(r*length/radix); passStr += ") * b";
								}
								else
								{
									passStr += "(("; passStr += SztToStr(numButterfly); passStr += "*me + ";
									passStr += SztToStr(butterflyIndex);
									passStr += ")%"; passStr += SztToStr(algLS); passStr += " + ";
									passStr += SztToStr(r*algLS); passStr += ") * b";
								}

								passStr += " );\n\t\t";
							}

							passStr += rType; passStr += " TR, TI;\n\t\t";

							if(realSpecial && (flag == SR_TWMUL_3STEP))
							{
								if(fwd)
								{
									passStr += "if(t==0)\n\t\t{\n\t\t";

									passStr += "TR = (W.x * "; passStr += regRealIndex; passStr += ") - (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = (W.y * "; passStr += regRealIndex; passStr += ") + (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";

									passStr += "}\n\t\telse\n\t\t{\n\t\t";
									
									passStr += "TR = (W.x * "; passStr += regRealIndex; passStr += ") + (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = (W.y * "; passStr += regRealIndex; passStr += ") - (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									
									passStr += "}\n\t\t";
								}
								else
								{
									passStr += "if(t==0)\n\t\t{\n\t\t";

									passStr += "TR = (W.x * "; passStr += regRealIndex; passStr += ") + (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = (W.y * "; passStr += regRealIndex; passStr += ") - (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";

									passStr += "}\n\t\telse\n\t\t{\n\t\t";

									passStr += "TR = (W.x * "; passStr += regRealIndex; passStr += ") - (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = (W.y * "; passStr += regRealIndex; passStr += ") + (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";

									passStr += "}\n\t\t";
								}
							}
							else
							{
								if(fwd)
								{
									passStr += "TR = (W.x * "; passStr += regRealIndex; passStr += ") - (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = (W.y * "; passStr += regRealIndex; passStr += ") + (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
								}
								else
								{
									passStr += "TR =  (W.x * "; passStr += regRealIndex; passStr += ") + (W.y * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
									passStr += "TI = -(W.y * "; passStr += regRealIndex; passStr += ") + (W.x * ";
									passStr += regImagIndex; passStr += ");\n\t\t";
								}
							}

							passStr += regRealIndex; passStr += " = TR;\n\t\t";
							passStr += regImagIndex; passStr += " = TI;\n\t}\n";

						}

						butterflyIndex++;
					}
				}
				else // write operation
				{
					for(size_t v=0; v<regC; v++)
					{
						for(size_t r=0; r<radix; r++)
						{
							if(realSpecial && (nextPass == NULL) && (r > (radix/2)))
								break;

							if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i != 0))
								break;

							if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i == 0))
								passStr += "\n\t}\n\tif( rw && !me)\n\t{";

							std::string regIndexC0;
								
							for(size_t c=cStart; c<cEnd; c++) // component loop: 0 - real, 1 - imaginary
							{
								std::string tail;
								std::string regIndex;
								regIndex = linearRegs ? "(*R" : regBaseCount;
								std::string buffer;

								// Write real & imag at once
								if(interleaved && (component == SR_COMP_BOTH) && linearRegs)
								{
									assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data
									buffer = bufferRe;
									RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ")";
									tail = "";
								}
								else
								{
									if(c == 0)
									{
										if(linearRegs) { RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").x"; }
										else		   { RegBaseAndCountAndPos("R", r, regIndex); }
										buffer = bufferRe;
										tail = interleaved ? ".x" : "";
									}
									else
									{
										if(linearRegs) { RegBaseAndCountAndPos("", i*radix + r, regIndex); regIndex += ").y"; }
										else		   { RegBaseAndCountAndPos("I", r, regIndex); }
										buffer = bufferIm;
										tail = interleaved ? ".y" : "";
									}
								}

								if(regC != 1)
								{
									regIndex += ".s";
									regIndex += SztToStr(v);
								}

								passStr += "\n\t";

								if(scale != 1.0f) { regIndex += " * "; regIndex += FloatToStr(scale); regIndex += FloatSuffix<PR>(); }
								if (c == 0) regIndexC0 += regIndex;

								bufOffset.clear();
								bufOffset += offset; bufOffset += " + ( ";
								if( (numButterfly * workGroupSize) > algLS )
								{
									bufOffset += "(("; bufOffset += SztToStr(numButterfly);
									bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex); bufOffset += ")/";
									bufOffset += SztToStr(algLS); bufOffset += ")*"; bufOffset += SztToStr(algL); bufOffset += " + (";
									bufOffset += SztToStr(numButterfly); bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex);
									bufOffset += ")%"; bufOffset += SztToStr(algLS); bufOffset += " + ";
								}
								else
								{
									bufOffset += SztToStr(numButterfly); bufOffset += "*me + "; bufOffset += SztToStr(butterflyIndex);
									bufOffset += " + ";
								}

								bufOffset += SztToStr(r*algLS); bufOffset += " )*"; bufOffset += SztToStr(stride); 
								
								if (fft_doPostCallback)
								{
									if(interleaved && (component == SR_COMP_BOTH))
									{
										if (c == (cEnd - 1))
										{
											passStr += "tempC.x = "; passStr += regIndexC0; passStr += ";\n\t";
											passStr += "tempC.y = "; passStr += regIndex; passStr += ";\n\t";

											passStr += fft_postCallback.funcname; passStr += "("; 
											passStr += buffer; passStr += ", (";
											passStr += bufOffset; passStr += "), post_userdata, tempC";
											if (fft_postCallback.localMemSize > 0)
											{
												passStr += ", post_localmem";
											}
											passStr += ");";
										}
									}
									else if (c == (cEnd - 1))
									{
										passStr += fft_postCallback.funcname; passStr += "("; 
										passStr += bufferRe; passStr += ", "; passStr += bufferIm; passStr += ", (";
										passStr += bufOffset; passStr += "), post_userdata, ("; 
										passStr += regIndexC0; passStr += "), ("; passStr += regIndex; passStr += ")";
										if (fft_postCallback.localMemSize > 0)
										{
											passStr += ", post_localmem";
										}
										passStr += ");";
									}
								}
								else
								{
									passStr += buffer; passStr += "["; passStr += bufOffset; passStr += "]";
								passStr += tail; passStr += " = "; passStr += regIndex;
								passStr += ";";
								}

								// Since we write real & imag at once, we break the loop
								if(interleaved && (component == SR_COMP_BOTH) && linearRegs)
									break;
							}

							if(realSpecial && (nextPass == NULL) && (r == radix/2) && (i == 0))
								passStr += "\n\t}\n\tif(rw)\n\t{";

						}

						butterflyIndex++;
					}

				}
			}

			assert(butterflyIndex <= numButterfly);
		}


		// Special SweepRegs function to carry out some R-C/C-R specific operations
		void SweepRegsRC(	size_t flag, bool fwd, bool interleaved, size_t stride, size_t component,
							double scale, bool setZero, bool batch2, bool oddt,
							const std::string &bufferRe, const std::string &bufferIm, const std::string &offset,
							std::string &passStr) const
		{
			assert( (flag == SR_READ ) ||
					(flag == SR_WRITE) );


			// component: 0 - real, 1 - imaginary, 2 - both
			size_t cStart, cEnd;
			switch(component)
			{
			case SR_COMP_REAL:	cStart = 0; cEnd = 1; break;
			case SR_COMP_IMAG:	cStart = 1; cEnd = 2; break;
			case SR_COMP_BOTH:	cStart = 0; cEnd = 2; break;
			default:	assert(false);
			}

			std::string rType  = RegBaseType<PR>(1);

			assert(r2c || c2r);
			assert(linearRegs);
			bool singlePass = ((position == 0) && (nextPass == NULL));

			size_t numCR = numButterfly * radix;
			if(!(numCR%2)) assert(!oddt);

			size_t rStart = 0;
			size_t rEnd = numCR;

			bool oddp = ((numCR%2) && (numCR > 1) && !setZero);
			if(oddp)
			{
				if(oddt)	{ rStart = numCR-1; rEnd = numCR+1; }
				else		{ rStart = 0;		rEnd = numCR-1; }
			}

			if(!oddp) assert(!oddt);

			for(size_t r=rStart; r<rEnd; r++)
			{
				std::string val1StrExt;

				for(size_t c=cStart; c<cEnd; c++) // component loop: 0 - real, 1 - imaginary
				{
					if(flag == SR_READ) // read operation
					{
						std::string tail, tail2;
						std::string regIndex = "(*R";
						std::string buffer;

						RegBaseAndCountAndPos("", r, regIndex); 
						if (fft_doPreCallback && c2r)
						{
							 regIndex += ")";
							 if (interleaved)
							 {
								buffer = (c == 0) ? bufferRe : bufferIm;
							 }
							 else
							 {
								 buffer += bufferRe; buffer += ", "; buffer += bufferIm;
							 }
						}
						else
						{
							if(c == 0)
							{
								regIndex += ").x";
								buffer = bufferRe;
								tail  = interleaved ? ".x;" : ";";
								tail2 = interleaved ? ".y;" : ";";
							}
							else
							{
								regIndex += ").y";
								buffer = bufferIm;
								tail  = interleaved ? ".y;" : ";";
								tail2 = interleaved ? ".x;" : ";";
							}
						}

						size_t bid = numCR/2;
						bid = bid ? bid : 1;
						size_t cid, lid;

						if(oddt)
						{
							cid = r%2;
							lid = 1 + (numCR/2);
						}
						else
						{
							cid = r/bid;
							lid = 1 + r%bid;
						}

						std::string oddpadd = oddp ? " (me/2) + " : " ";

						std::string idxStr, idxStrRev;
						if((length <= 2) || ((length & (length - 1)) != 0))
						{
							idxStr += SztToStr(bid); idxStr += "*me +"; idxStr += oddpadd; idxStr += SztToStr(lid);
						}
						else
						{
							idxStr += "me + "; idxStr += SztToStr(1 + length*(r%bid)/numCR); idxStr += oddpadd;
						}
						idxStrRev += SztToStr(length); idxStrRev += " - ("; idxStrRev += idxStr; idxStrRev += " )";

						bool act = ( fwd || ((cid == 0) && (!batch2)) || ((cid != 0) && batch2) );
						if(act)
						{
							passStr += "\n\t";
							passStr += regIndex;
							passStr += " = ";
						}

						if(setZero)
						{
							if(act) passStr += "0;";
						}
						else
						{
							if(act)
							{
								if (fft_doPreCallback)
								{
									passStr += fft_preCallback.funcname; passStr += "(";
									passStr += buffer; passStr += ", ";  
								}
								else
								{
									passStr += buffer;
									passStr += "["; 
								}
								passStr += offset; passStr += " + ( ";
							}

							if(fwd)
							{
								if(cid == 0)	passStr += idxStr;
								else			passStr += idxStrRev;
							}
							else
							{
								if(cid == 0)	{ if(!batch2) passStr += idxStr; }
								else			{ if(batch2)  passStr += idxStr; }
							}

							if(act)
							{
								passStr += " )*"; passStr += SztToStr(stride); 
								
								if (fft_doPreCallback)
								{
									passStr += ", pre_userdata";
									passStr += (fft_preCallback.localMemSize > 0) ? ", localmem);" : ");";
								}
								else
								{
									passStr += "]";

									if(fwd) { passStr += tail; }
									else	{ if(!batch2) passStr += tail; else passStr += tail2; }
								}
							}
						}
					}
					else // write operation
					{
						std::string tail;
						std::string regIndex = "(*R";
						std::string regIndexPair = "(*R";
						std::string buffer;

						// Write real & imag at once
						if(interleaved && (component == SR_COMP_BOTH))
						{
							assert(bufferRe.compare(bufferIm) == 0); // Make sure Real & Imag buffer strings are same for interleaved data
							buffer = bufferRe;
						}
						else
						{
							if(c == 0)
							{
								buffer = bufferRe;
								tail = interleaved ? ".x" : "";
							}
							else
							{
								buffer = bufferIm;
								tail = interleaved ? ".y" : "";
							}
						}


						size_t bid, cid, lid;
						if(singlePass && fwd)
						{
							bid = 1 + radix/2;
							lid = r;
							cid = r/bid;

							RegBaseAndCountAndPos("", r, regIndex); regIndex += ")";
							RegBaseAndCountAndPos("", (radix - r)%radix , regIndexPair); regIndexPair += ")";
						}
						else
						{
							bid = numCR/2;

							if(oddt)
							{
								cid = r%2;
								lid = 1 + (numCR/2);

								RegBaseAndCountAndPos("", r, regIndex); regIndex += ")";
								RegBaseAndCountAndPos("", r + 1, regIndexPair); regIndexPair += ")";
							}
							else
							{
								cid = r/bid;
								lid = 1 + r%bid;

								RegBaseAndCountAndPos("", r, regIndex); regIndex += ")";
								RegBaseAndCountAndPos("", r + bid, regIndexPair); regIndexPair += ")";
							}
						}


						if(!cid)
						{
							std::string oddpadd = oddp ? " (me/2) + " : " ";

							std::string sclStr = "";
							if(scale != 1.0f) { sclStr += " * "; sclStr += FloatToStr(scale); sclStr += FloatSuffix<PR>(); }

							if(fwd)
							{
								std::string idxStr, idxStrRev;
								if((length <= 2) || ((length & (length - 1)) != 0))
								{
								idxStr += SztToStr(length/(2*workGroupSize)); idxStr += "*me +"; idxStr += oddpadd; idxStr += SztToStr(lid);
								}
								else
								{
								idxStr += "me + "; idxStr += SztToStr(1 + length*(r%bid)/numCR); idxStr += oddpadd;
								}
								idxStrRev += SztToStr(length); idxStrRev += " - ("; idxStrRev += idxStr; idxStrRev += " )";

								std::string val1Str, val2Str;

								if (fft_doPostCallback && !rcFull)
								{
									if (interleaved)
									{
										val1Str += "\n\t";
										val1Str += fft_postCallback.funcname; val1Str += "("; val1Str += buffer; val1Str += ", ";
										val1Str += offset; val1Str += " + ( "; val1Str += idxStr; val1Str += " )*"; val1Str += SztToStr(stride);
										val1Str += ", post_userdata, ";
									}
									else if (c == 0)
									{
										val1StrExt += "\n\t";
										val1StrExt += fft_postCallback.funcname; val1StrExt += "("; val1StrExt += bufferRe; val1StrExt += ", ";
										val1StrExt += bufferIm; val1StrExt += ", "; val1StrExt += offset; val1StrExt += " + ( "; val1StrExt += idxStr; 
										val1StrExt += " )*"; val1StrExt += SztToStr(stride); val1StrExt += ", post_userdata, ";
									}									
								}
								else
								{
								val1Str += "\n\t";
								val1Str += buffer; val1Str += "["; val1Str += offset; val1Str += " + ( ";
								val1Str += idxStr; val1Str += " )*"; val1Str += SztToStr(stride); val1Str += "]";
								val1Str += tail; val1Str += " = ";
								}

								val2Str += "\n\t";
								val2Str += buffer; val2Str += "["; val2Str += offset; val2Str += " + ( ";
								val2Str += idxStrRev; val2Str += " )*"; val2Str += SztToStr(stride); val2Str += "]";
								val2Str += tail; val2Str += " = ";

								std::string real1, imag1, real2, imag2;

								real1 +=  "("; real1 += regIndex; real1 += ".x + "; real1 += regIndexPair; real1 += ".x)*0.5";
								imag1 +=  "("; imag1 += regIndex; imag1 += ".y - "; imag1 += regIndexPair; imag1 += ".y)*0.5";
								real2 +=  "("; real2 += regIndex; real2 += ".y + "; real2 += regIndexPair; real2 += ".y)*0.5";
								imag2 += "(-"; imag2 += regIndex; imag2 += ".x + "; imag2 += regIndexPair; imag2 += ".x)*0.5";

								if(interleaved && (component == SR_COMP_BOTH))
								{
									val1Str += "("; val1Str += RegBaseType<PR>(2); val1Str += ")( ";
									val2Str += "("; val2Str += RegBaseType<PR>(2); val2Str += ")( ";

									if(!batch2) { val1Str += real1; val1Str += ", "; val1Str += "+"; val1Str += imag1;
												  val2Str += real1; val2Str += ", "; val2Str += "-"; val2Str += imag1; }
									else		{ val1Str += real2; val1Str += ", "; val1Str += "+"; val1Str += imag2;
												  val2Str += real2; val2Str += ", "; val2Str += "-"; val2Str += imag2; }

									val1Str += " )";
									val2Str += " )";
								}
								else
								{
									val1Str += " (";
									val2Str += " (";
									if(c == 0)
									{
										if(!batch2) { val1Str += real1;
													  val2Str += real1; }
										else		{ val1Str += real2;
													  val2Str += real2; }
									}
									else
									{
										if(!batch2) { val1Str += "+"; val1Str += imag1;
													  val2Str += "-"; val2Str += imag1; }
										else		{ val1Str += "+"; val1Str += imag2;
													  val2Str += "-"; val2Str += imag2; }
									}
									val1Str += " )";
									val2Str += " )";
								}

								val1Str += sclStr;
								val2Str += sclStr;

								if (fft_doPostCallback && !rcFull) 
								{
									if (!interleaved) 
									{
										val1StrExt += val1Str;
										val1Str.clear();

										if(c == 0) val1StrExt += ", ";
										else	val1Str += val1StrExt;
									}

									if (interleaved || c == (cEnd - 1))
									{
										if (fft_postCallback.localMemSize > 0)	val1Str += ", localmem";
										val1Str += ");";
									}
								}
								else
								{
									val1Str += ";";
								}

												passStr += val1Str; 
								if(rcFull)	{	passStr += val2Str; passStr += ";"; }
							}
							else
							{
								std::string idxStr, idxStrRev;
								if((length <= 2) || ((length & (length - 1)) != 0))
								{
								idxStr += SztToStr(bid); idxStr += "*me +"; idxStr += oddpadd; idxStr += SztToStr(lid);
								}
								else
								{								
								idxStr += "me + "; idxStr += SztToStr(1 + length*(r%bid)/numCR); idxStr += oddpadd;
								}
								idxStrRev += SztToStr(length); idxStrRev += " - ("; idxStrRev += idxStr; idxStrRev += " )";

								passStr += "\n\t";
								passStr += buffer; passStr += "["; passStr += offset; passStr += " + ( ";

								if(!batch2)	passStr += idxStr;
								else		passStr += idxStrRev;

								passStr += " )*"; passStr += SztToStr(stride); passStr += "]";
								passStr += tail; passStr += " = ";

								passStr += "( ";
								if(c == 0)
								{
									regIndex += ".x"; 
									regIndexPair += fft_doPreCallback ? ".y" : ".x";

									if(!batch2)	{ passStr += regIndex; passStr += " - "; passStr += regIndexPair; }
									else		{ passStr += regIndex; passStr += " + "; passStr += regIndexPair; }
								}
								else
								{
									regIndex += ".y"; regIndexPair += (fft_doPreCallback && oddt) ? ".x" : ".y";

									if(!batch2)	{					passStr += regIndex; passStr += " + "; passStr += regIndexPair; }
									else		{ passStr += " - "; passStr += regIndex; passStr += " + "; passStr += regIndexPair; }
								}
								passStr += " )";
								passStr += sclStr;
								passStr += ";";
							}



							// Since we write real & imag at once, we break the loop
							if(interleaved && (component == SR_COMP_BOTH))
								break;
						}
					}
				}
			}

		}


		void CallButterfly(const std::string &bflyName, size_t regC, size_t numB, std::string &passStr) const
		{
			std::string regBase;
			RegBase(regC, regBase);

			for(size_t i=0; i<numB; i++)
			{
				std::string regBaseCount = regBase;
				RegBaseAndCount(i, regBaseCount);

				passStr += "\n\t";
				passStr += bflyName;
				passStr += "(";

				for(size_t r=0; ; r++)
				{
					if(linearRegs)
					{
						std::string regIndex = "R";
						RegBaseAndCountAndPos("", i*radix + r, regIndex);

						passStr += regIndex;
					}
					else
					{
						std::string regRealIndex(regBaseCount);
						std::string regImagIndex(regBaseCount);
						RegBaseAndCountAndPos("R", r, regRealIndex);
						RegBaseAndCountAndPos("I", r, regImagIndex);

						passStr += "&"; passStr += regRealIndex; passStr += ", ";
						passStr += "&"; passStr += regImagIndex;
					}

					if(r == radix-1)
					{
						passStr += ");";
						break;
					}
					else
					{
						passStr += ", ";
					}
				}
			}
		}

    public:
		Pass(	size_t positionVal, size_t lengthVal, size_t radixVal, size_t cnPerWIVal,
				size_t L, size_t LS, size_t R, bool linearRegsVal, bool halfLdsVal,
				bool r2cVal, bool c2rVal, bool rcFullVal, bool rcSimpleVal, bool realSpecialVal) :
			position(positionVal), length(lengthVal), radix(radixVal), cnPerWI(cnPerWIVal),
			algL(L), algLS(LS), algR(R), linearRegs(linearRegsVal), halfLds(halfLdsVal),
			r2c(r2cVal), c2r(c2rVal), rcFull(rcFullVal), rcSimple(rcSimpleVal), realSpecial(realSpecialVal),
			enableGrouping(true),
			numB1(0), numB2(0), numB4(0),
			nextPass(NULL), fft_doPreCallback(false), fft_doPostCallback(false)
		{
			assert(radix <= length);
			assert(length%radix == 0);

			numButterfly = cnPerWI/radix;
			workGroupSize = length/cnPerWI;

			// Total number of butterflies (over all work-tems) must be divisible by LS
			assert( ((numButterfly*workGroupSize)%algLS) == 0 );

			// All butterflies in one work-item should always be part of no more than 1 FFT transform.
			// In other words, there should not be more than 1 FFT transform per work-item.
			assert(cnPerWI <= length);

			// Calculate the different types of Butterflies needed
			if(linearRegs || r2c || c2r)
			{
				numB1 = numButterfly;
			}
			else
			{
				numB4 = numButterfly/4;
				numB2 = (numButterfly%4)/2; // can be 0 or 1
				numB1 = (numButterfly%2); // can be 0 or 1

				assert(numButterfly == (numB4*4 + numB2*2 + numB1));
			}

			// if only half LDS can be used, we need the passes to share registers
			// and hence they need to be linear registers
			if(halfLds) assert(linearRegs);
		}

		size_t GetNumB1() const { return numB1; }
		size_t GetNumB2() const { return numB2; }
		size_t GetNumB4() const { return numB4; }

		size_t GetPosition() const { return position; }
		size_t GetRadix() const { return radix; }

		void SetNextPass(Pass<PR> *np) { nextPass = np; }
		void SetGrouping(bool grp) { enableGrouping = grp; }

		void SetPrecallback(bool hasPrecallback, clfftCallbackParam precallbackParam) 
		{ 
			fft_doPreCallback = hasPrecallback;
			fft_preCallback = precallbackParam;
		}

		void SetPostcallback(bool hasPostcallback, clfftCallbackParam postcallbackParam) 
		{ 
			fft_doPostCallback = hasPostcallback;
			fft_postCallback = postcallbackParam;
		}

		void GeneratePass(	bool fwd, std::string &passStr, bool fft_3StepTwiddle, bool twiddleFront,
							bool inInterleaved, bool outInterleaved,
							bool inReal, bool outReal,
							size_t inStride, size_t outStride, double scale,
							bool gIn = false, bool gOut = false) const
		{
			const std::string bufferInRe  = (inReal || inInterleaved) ?   "bufIn"  : "bufInRe";
			const std::string bufferInIm  = (inReal || inInterleaved) ?   "bufIn"  : "bufInIm";
			const std::string bufferOutRe = (outReal || outInterleaved) ? "bufOut" : "bufOutRe";
			const std::string bufferOutIm = (outReal || outInterleaved) ? "bufOut" : "bufOutIm";

			const std::string bufferInRe2  = (inReal || inInterleaved) ?   "bufIn2"  : "bufInRe2";
			const std::string bufferInIm2  = (inReal || inInterleaved) ?   "bufIn2"  : "bufInIm2";
			const std::string bufferOutRe2 = (outReal || outInterleaved) ? "bufOut2" : "bufOutRe2";
			const std::string bufferOutIm2 = (outReal || outInterleaved) ? "bufOut2" : "bufOutIm2";

			// for real transforms we use only B1 butteflies (regC = 1)
			if(r2c || c2r)
			{
				assert(numB1 == numButterfly);
				assert(linearRegs);
			}

			// Check if it is single pass transform
			bool singlePass = ((position == 0) && (nextPass == NULL));
			if(singlePass) assert(numButterfly == 1); // for single pass transforms, there can be only 1 butterfly per transform
			if(singlePass) assert(workGroupSize == 1);

			// Register types
			std::string regB1Type = RegBaseType<PR>(1);
			std::string regB2Type = RegBaseType<PR>(2);
			std::string regB4Type = RegBaseType<PR>(4);

			//Function attribute
			passStr += "__attribute__((always_inline)) void\n";

			//Function name
			passStr += PassName(position, fwd);

			// Function arguments
			passStr += "(";
			passStr += "uint rw, uint b, ";
			if(realSpecial) passStr += "uint t, ";
			passStr += "uint me, uint inOffset, uint outOffset, ";

			if(r2c || c2r)
			{
				assert(halfLds);

				if(gIn)
				{
					if(inInterleaved)
					{
										passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
						if(!rcSimple) {	passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferInRe2; passStr += ", "; }
					}
					else if(inReal)
					{
										passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
						if(!rcSimple) {	passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe2; passStr += ", "; }
					}
					else
					{
										passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
						if(!rcSimple) {	passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe2; passStr += ", "; }
										passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInIm;  passStr += ", ";
						if(!rcSimple) {	passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInIm2; passStr += ", "; }
					}
				}
				else
				{
					passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe; passStr += ", ";
					passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferInIm; passStr += ", ";
				}

				if(gOut)
				{
					if(outInterleaved)
					{
															passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferOutRe;
						if(!rcSimple) { passStr += ", ";	passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferOutRe2; }
					}
					else if(outReal)
					{
															passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe;
						if(!rcSimple) { passStr += ", ";	passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe2; }
					}
					else
					{
															passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe;  passStr += ", ";
						if(!rcSimple) {						passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe2; passStr += ", "; }
															passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutIm;
						if(!rcSimple) { passStr += ", ";	passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutIm2; }
					}
				}
				else
				{
					passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe; passStr += ", ";
					passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferOutIm;
				}
			}
			else
			{
				if(gIn)
				{
					if(inInterleaved)
					{
						passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
					}
					else
					{
						passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
						passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferInIm;  passStr += ", ";
					}
				}
				else
				{
					if(inInterleaved)
					{
						passStr += "__local "; passStr += regB2Type; passStr += " *"; passStr += bufferInRe;  passStr += ", ";
					}
					else
					{
						passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferInRe; passStr += ", ";
						passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferInIm; passStr += ", ";
					}
				}


				if(gOut)
				{
					if(outInterleaved)
					{
						passStr += "__global "; passStr += regB2Type; passStr += " *"; passStr += bufferOutRe;
					}
					else
					{
						passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe;  passStr += ", ";
						passStr += "__global "; passStr += regB1Type; passStr += " *"; passStr += bufferOutIm;
					}
				}
				else
				{
					if(outInterleaved)
					{
						passStr += "__local "; passStr += regB2Type; passStr += " *"; passStr += bufferOutRe;
					}
					else
					{
						passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferOutRe; passStr += ", ";
						passStr += "__local "; passStr += regB1Type; passStr += " *"; passStr += bufferOutIm;
					}
				}
			}

			// Register arguments
			if(linearRegs)
			{
				passStr += ", "; passStr += IterRegArgs();
			}

			if (fft_doPreCallback || fft_doPostCallback)
			{
				//Include pre-callback parameters if pre-callback is set
			if (fft_doPreCallback )
			{
				if ((r2c && !rcSimple) || c2r)
				{
					passStr += ", uint inOffset2";
				}

					passStr += ", __global void* pre_userdata";
				}

				//Include post-callback parameters if post-callback is set
				if (fft_doPostCallback )
				{
					if (r2c || (c2r && !rcSimple))
					{
						passStr += ", uint outOffset2";
					}
					passStr += ", __global void* post_userdata";
				}

				if (fft_doPreCallback && fft_preCallback.localMemSize > 0)
				{
					passStr += ", __local void* localmem";
				}
				if (fft_doPostCallback && fft_postCallback.localMemSize > 0)
				{
					passStr += ", __local void* post_localmem";
				}
			}

			passStr += ")\n{\n";

			// Register Declarations
			if(!linearRegs)
			{
				DeclareRegs(regB1Type, 1, numB1, passStr);
				DeclareRegs(regB2Type, 2, numB2, passStr);
				DeclareRegs(regB4Type, 4, numB4, passStr);
			}

			// odd cnPerWI processing
			bool oddp = false;
			oddp = ((cnPerWI%2) && (length > 1) && (!singlePass));

			// additional register for odd
			if( !rcSimple && oddp && ((r2c && (nextPass == NULL)) || (c2r && (position == 0))) )
			{
				passStr += "\n\t";
				passStr += "uint brv = 0;\n\t";
				passStr += "\n\t";
				passStr += regB2Type; passStr += " R"; passStr += SztToStr(cnPerWI); passStr += "[1];\n\t";
				passStr += "(*R"; passStr += SztToStr(cnPerWI); passStr += ").x = 0; ";
				passStr += "(*R"; passStr += SztToStr(cnPerWI); passStr += ").y = 0;\n";
			}

			// Special private memory for c-r 1 pass transforms
			if( !rcSimple && (c2r && (position == 0)) && singlePass )
			{
				assert(radix == length);

				passStr += "\n\t";
				passStr += regB1Type;
				passStr += " mpvt["; passStr += SztToStr(length); passStr += "];\n";
			}

			passStr += "\n";

			// Read into registers
			if(r2c)
			{
				if(position == 0)
				{
					passStr += "\n\tif(rw)\n\t{";
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 1, numB1, 0, passStr);
					passStr += "\n\t}\n";

					if(rcSimple)
					{
						passStr += "\n";
						SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, true, true, false, bufferInRe2, bufferInIm2, "inOffset", passStr);
						passStr += "\n";
					}
					else
					{
						passStr += "\n\tif(rw > 1)\n\t{";
						if (fft_doPreCallback)
						{
							SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, bufferInRe2, bufferInIm2, "inOffset2", 1, numB1, 0, passStr);
						}
						else
						{
							SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, bufferInRe2, bufferInIm2, "inOffset", 1, numB1, 0, passStr);
						}
						passStr += "\n\t}\n";

						passStr += "\telse\n\t{";
						SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, true, true, false, bufferInRe2, bufferInIm2, "inOffset", passStr);
						passStr += "\n\t}\n";
					}
				}
			}
			else if(c2r && !rcSimple)
			{
				if(position == 0)
				{
					std::string processBufRe = bufferOutRe;
					std::string processBufIm = bufferOutIm;
					std::string processBufOffset = "outOffset";
					size_t processBufStride = outStride;

					if(singlePass)
					{
						processBufRe = "mpvt";
						processBufIm = "mpvt";
						processBufOffset = "0";
						processBufStride = 1;
					}

					passStr += "\n\tif(rw && !me)\n\t{\n\t";
					passStr += processBufRe; passStr += "["; passStr += processBufOffset; passStr += "] = ";
					
					if (fft_doPreCallback)
					{
						passStr += fft_preCallback.funcname; passStr += "("; passStr += bufferInRe; 
						if (!inInterleaved) { passStr += ", "; passStr += bufferInIm; }
						passStr += ", inOffset, pre_userdata";
						passStr += fft_preCallback.localMemSize > 0 ? ", localmem)" : ")";
					}
					else
					{
						passStr += bufferInRe; passStr+= "[inOffset]";
					}

					if(inInterleaved || fft_doPreCallback) passStr += ".x;\n\t}"; else passStr += ";\n\t}";

					if(length > 1)
					{
						passStr += "\n\n\tif(rw)\n\t{";
						if (fft_doPreCallback && !inInterleaved)
						{
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, false, bufferInRe, bufferInIm, "inOffset", passStr);
						}
						else
						{
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, false, bufferInRe, bufferInRe, "inOffset", passStr);
						}
						passStr += "\n\t}\n";

						passStr += "\n\tif(rw > 1)\n\t{";
						if (fft_doPreCallback)
						{
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, true, false, bufferInRe2, bufferInIm2, "inOffset2", passStr);
						}
						else
						{
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, true, false, bufferInIm2, bufferInIm2, "inOffset", passStr);
						}

						passStr += "\n\t}\n\telse\n\t{";
						SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, true, true, false, bufferInIm2, bufferInIm2, "inOffset", passStr);
						passStr += "\n\t}\n";

						if(oddp)
						{
							passStr += "\n\tif(rw && (me%2))\n\t{";
							if (fft_doPreCallback)
							{
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, true, bufferInRe, bufferInIm, "inOffset", passStr);
							}
							else
							{
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, true, bufferInRe, bufferInRe, "inOffset", passStr);
							}
							passStr += "\n\t}";
							passStr += "\n\tif((rw > 1) && (me%2))\n\t{";
							if (fft_doPreCallback)
							{
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, true, true, bufferInRe2, bufferInIm2, "inOffset2", passStr);
							}
							else
							{
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, true, true, bufferInIm2, bufferInIm2, "inOffset", passStr);
							}
							passStr += "\n\t}\n";
						}


						SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_REAL, 1.0f, false, true, false, processBufRe, processBufIm, processBufOffset, passStr);
						if(oddp)
						{
							passStr += "\n\tif(me%2)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_REAL, 1.0f, false, true, true, processBufRe, processBufIm, processBufOffset, passStr);
							passStr += "\n\t}\n";
						}
						SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_REAL, 1.0f, false, false, false, processBufRe, processBufIm, processBufOffset, passStr);
						if(oddp)
						{
							passStr += "\n\tif(me%2)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_REAL, 1.0f, false, false, true, processBufRe, processBufIm, processBufOffset, passStr);
							passStr += "\n\t}\n";
						}
					}

					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
					SweepRegs(SR_READ, fwd, outInterleaved, processBufStride, SR_COMP_REAL, 1.0f, false, processBufRe, processBufIm, processBufOffset, 1, numB1, 0, passStr, false, false, oddp);
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";


					passStr += "\n\tif((rw > 1) && !me)\n\t{\n\t";
					passStr += processBufIm; passStr += "["; passStr += processBufOffset; passStr += "] = ";
					
					if (fft_doPreCallback)
					{
						passStr += fft_preCallback.funcname; passStr += "("; passStr += bufferInRe2; 
						if (!inInterleaved) { passStr += ", "; passStr += bufferInIm2; }
						passStr += ", inOffset2, pre_userdata";
						passStr += fft_preCallback.localMemSize > 0 ? ", localmem)" : ")";
					}
					else
					{
						passStr += bufferInRe2; passStr+= "[inOffset]";
					}
					if(inInterleaved || fft_doPreCallback) passStr += ".x;\n\t}"; else passStr += ";\n\t}";
					passStr += "\n\tif((rw == 1) && !me)\n\t{\n\t"; passStr += processBufIm; passStr += "["; passStr += processBufOffset; passStr += "] = 0;\n\t}";


					if(length > 1)
					{
						if (!fft_doPreCallback)
						{
							passStr += "\n\n\tif(rw)\n\t{";
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, false, false, bufferInIm, bufferInIm, "inOffset", passStr);
							passStr += "\n\t}\n";

							passStr += "\n\tif(rw > 1)\n\t{";
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, true, false, bufferInRe2, bufferInRe2, "inOffset", passStr);
							passStr += "\n\t}\n\telse\n\t{";
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, true, true, false, bufferInRe2, bufferInRe2, "inOffset", passStr);
							passStr += "\n\t}";
						
						
							if(oddp)
							{
								passStr += "\n\tif(rw && (me%2))\n\t{";
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, false, true, bufferInIm, bufferInIm, "inOffset", passStr);
								passStr += "\n\t}";
								passStr += "\n\tif((rw > 1) && (me%2))\n\t{";
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, true, true, bufferInRe2, bufferInRe2, "inOffset", passStr);
								passStr += "\n\t}";
							}
						}
						passStr += "\n";

						SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_IMAG, 1.0f, false, true, false, processBufRe, processBufIm, processBufOffset, passStr);
						if(oddp)
						{
							passStr += "\n\tif(me%2)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_IMAG, 1.0f, false, true, true, processBufRe, processBufIm, processBufOffset, passStr);
							passStr += "\n\t}\n";
						}
						SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_IMAG, 1.0f, false, false, false, processBufRe, processBufIm, processBufOffset, passStr);
						if(oddp)
						{
							passStr += "\n\tif(me%2)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, processBufStride, SR_COMP_IMAG, 1.0f, false, false, true, processBufRe, processBufIm, processBufOffset, passStr);
							passStr += "\n\t}\n";
						}
					}

					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
					SweepRegs(SR_READ, fwd, outInterleaved, processBufStride, SR_COMP_IMAG, 1.0f, false, processBufRe, processBufIm, processBufOffset, 1, numB1, 0, passStr);
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
				}
			}
			else
			{
				if( (!halfLds) || (halfLds && (position == 0)) )
				{
					bool isPrecallVector = false;
					//If precallback is set
					if (fft_doPreCallback)
					{
						passStr += "\n\t"; passStr += regB2Type; passStr += " retPrecallback"; 
						
						if (numB4 > 0 || numB2 > 0)
						{
							passStr += "["; 
							passStr += (numB4 > 0) ? "4" : (numB2 > 0) ? "2" : "1"; 
							passStr += "]";

							isPrecallVector = true;
						}
						passStr += ";";
					}
					passStr += "\n\tif(rw)\n\t{";
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 1, numB1, 0, passStr, false, isPrecallVector);
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 2, numB2, numB1, passStr, false, isPrecallVector);
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 4, numB4, 2*numB2 + numB1, passStr, false, isPrecallVector);
					passStr += "\n\t}\n";
					passStr += "\n\telse\n\t{";
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 1, numB1, 0, passStr, true, isPrecallVector);
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 2, numB2, numB1, passStr, true, isPrecallVector);
					SweepRegs(SR_READ, fwd, inInterleaved, inStride, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 4, numB4, 2 * numB2 + numB1, passStr, true, isPrecallVector);
					passStr += "\n\t}\n";
				}
			}

			passStr += "\n";

			// 3-step twiddle multiplies done in the front
			bool tw3Done = false;
			if(fft_3StepTwiddle && twiddleFront)
			{
				tw3Done = true;
				if(linearRegs)
				{
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, true, bufferInRe, bufferInIm, "", 1, numB1, 0, passStr);
				}
				else
				{
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, true, bufferInRe, bufferInIm, "", 1, numB1, 0, passStr);
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, true, bufferInRe, bufferInIm, "", 2, numB2, numB1, passStr);
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, true, bufferInRe, bufferInIm, "", 4, numB4, 2*numB2 + numB1, passStr);
				}
			}

			passStr += "\n";

			// Twiddle multiply
			if( (position > 0) && (radix > 1) )
			{
				SweepRegs(SR_TWMUL, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 1, numB1, 0, passStr);
				SweepRegs(SR_TWMUL, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 2, numB2, numB1, passStr);
				SweepRegs(SR_TWMUL, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 4, numB4, 2*numB2 + numB1, passStr);
			}

			// Butterfly calls
			if(radix > 1)
			{
				if(numB1) CallButterfly(ButterflyName(radix, 1, fwd), 1, numB1, passStr);
				if(numB2) CallButterfly(ButterflyName(radix, 2, fwd), 2, numB2, passStr);
				if(numB4) CallButterfly(ButterflyName(radix, 4, fwd), 4, numB4, passStr);
			}


			if(!halfLds) passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
			passStr += "\n\n";

			// 3-step twiddle multiplies
			if(fft_3StepTwiddle && !tw3Done)
			{
				assert(nextPass == NULL);
				if(linearRegs)
				{
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 1, numB1, 0, passStr);
				}
				else
				{
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 1, numB1, 0, passStr);
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 2, numB2, numB1, passStr);
					SweepRegs(SR_TWMUL_3STEP, fwd, false, 1, SR_COMP_BOTH, 1.0f, false, bufferInRe, bufferInIm, "", 4, numB4, 2*numB2 + numB1, passStr);
				}
			}

			// Write back from registers
			if(halfLds)
			{
				// In this case, we have to write & again read back for the next pass since we are
				// using only half the lds. Number of barriers will increase at the cost of halving the lds.

				if(nextPass == NULL) // last pass
				{
					if(r2c && !rcSimple)
					{
						if(!singlePass)
						{
							SweepRegs(SR_WRITE, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 1, numB1, 0, passStr);
							passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, false, bufferInRe, bufferInIm, "inOffset", passStr);
							if(oddp)
							{
								passStr += "\n\tif(me%2)\n\t{";
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_REAL, 1.0f, false, false, true, bufferInRe, bufferInIm, "inOffset", passStr);
								passStr += "\n\t}\n";
							}

							passStr += "\n\tif(rw && !me)\n\t{\n\t";
							if(outInterleaved)
							{
								if (fft_doPostCallback)
								{
									passStr += fft_postCallback.funcname; passStr += "(bufOut, outOffset, post_userdata, ";
									passStr += "("; passStr += RegBaseType<PR>(2); passStr += ") ( ("; passStr += bufferInRe; passStr += "[inOffset]";
									if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ") , 0 )";
									if (fft_postCallback.localMemSize > 0)
									{
										passStr += ", localmem";
									}
									passStr += ");\n\t}";
								}
								else
								{
								passStr += bufferOutRe; passStr+= "[outOffset].x = "; passStr += bufferInRe; passStr += "[inOffset]";
								if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ";\n\t";
								passStr += bufferOutIm; passStr+= "[outOffset].y = "; passStr += "0;\n\t}";
							}
							}
							else
							{
								if (fft_doPostCallback)
								{
									passStr += fft_postCallback.funcname; passStr += "("; passStr += bufferOutRe; passStr += ", "; passStr += bufferOutIm;
									passStr += ", outOffset, post_userdata, "; passStr += bufferInRe; passStr += "[inOffset]";
									if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ", 0";
									if (fft_postCallback.localMemSize > 0)
									{
										passStr += ", localmem";
									}
									passStr += ");\n\t}";
								}
							else
							{
								passStr += bufferOutRe; passStr+= "[outOffset] = ";   passStr += bufferInRe; passStr += "[inOffset]";
								if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ";\n\t";
								passStr += bufferOutIm; passStr+= "[outOffset] = ";   passStr += "0;\n\t}";
							}
							}
							passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";


							SweepRegs(SR_WRITE, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, bufferInRe, bufferInIm, "inOffset", 1, numB1, 0, passStr);
							passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
							SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, false, false, bufferInRe, bufferInIm, "inOffset", passStr);
							if(oddp)
							{
								passStr += "\n\tif(me%2)\n\t{";
								SweepRegsRC(SR_READ, fwd, inInterleaved, inStride, SR_COMP_IMAG, 1.0f, false, false, true, bufferInRe, bufferInIm, "inOffset", passStr);
								passStr += "\n\t}\n";
							}

							passStr += "\n\tif((rw > 1) && !me)\n\t{\n\t";
							if(outInterleaved)
							{
								if (fft_doPostCallback)
								{
									passStr += fft_postCallback.funcname; passStr += "(bufOut2, outOffset2, post_userdata, ";
									passStr += "("; passStr += RegBaseType<PR>(2); passStr += ") ( ("; passStr += bufferInIm; passStr += "[inOffset]";
									if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ") , 0 )";
									if (fft_postCallback.localMemSize > 0)
									{
										passStr += ", localmem";
									}
									passStr += ");\n\t}";
								}
								else
								{
								passStr += bufferOutRe2; passStr+= "[outOffset].x = "; passStr += bufferInIm; passStr += "[inOffset]";
								if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ";\n\t";
								passStr += bufferOutIm2; passStr+= "[outOffset].y = "; passStr += "0;\n\t}";
							}
							}
							else
							{
								if (fft_doPostCallback)
								{
									passStr += fft_postCallback.funcname; passStr += "("; passStr += bufferOutRe2; passStr += ", "; passStr += bufferOutIm2;
									passStr+= ", outOffset2, post_userdata, "; passStr += bufferInIm; passStr += "[inOffset]";
									if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ", 0";
									if (fft_postCallback.localMemSize > 0)
									{
										passStr += ", localmem";
									}
									passStr += ");\n\t}";
								}
							else
							{
								passStr += bufferOutRe2; passStr+= "[outOffset] = ";   passStr += bufferInIm; passStr += "[inOffset]";
								if(scale != 1.0) { passStr += " * "; passStr += FloatToStr(scale); passStr += FloatSuffix<PR>(); } passStr += ";\n\t";
								passStr += bufferOutIm2; passStr+= "[outOffset] = ";   passStr += "0;\n\t}";
							}
							}
							passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
						}


						passStr += "\n\n\tif(rw)\n\t{";
						SweepRegsRC(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, false, false, bufferOutRe, bufferOutIm, "outOffset", passStr);
						passStr += "\n\t}\n";
						if(oddp)
						{
							passStr += "\n\n\tbrv = ((rw != 0) & (me%2 == 1));\n\t";
							passStr += "if(brv)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, false, true, bufferOutRe, bufferOutIm, "outOffset", passStr);
							passStr += "\n\t}\n";
						}

						passStr += "\n\n\tif(rw > 1)\n\t{";
						
						std::string outOffset;
						outOffset += "outOffset";
						if (fft_doPostCallback) outOffset += "2";

						SweepRegsRC(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, true, false, bufferOutRe2, bufferOutIm2, outOffset, passStr);
						passStr += "\n\t}\n";
						if(oddp)
						{
							passStr += "\n\n\tbrv = ((rw > 1) & (me%2 == 1));\n\t";
							passStr += "if(brv)\n\t{";
							SweepRegsRC(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, true, true, bufferOutRe2, bufferOutIm2, outOffset, passStr);
							passStr += "\n\t}\n";
						}

					}
					else if(c2r)
					{
						passStr += "\n\tif(rw)\n\t{";
						SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_REAL, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, numB1, 0, passStr);
						passStr += "\n\t}\n";

						if(!rcSimple)
						{
							std::string outOffset;
							outOffset += "outOffset";
							if (fft_doPostCallback) outOffset += "2";

							passStr += "\n\tif(rw > 1)\n\t{";
							SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_IMAG, scale, false, bufferOutRe2, bufferOutIm2, outOffset, 1, numB1, 0, passStr);
							passStr += "\n\t}\n";
						}
					}
					else
					{
						passStr += "\n\tif(rw)\n\t{";
						SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, numB1, 0, passStr);
						passStr += "\n\t}\n";
					}
				}
				else
				{
					passStr += "\n\tif(rw)\n\t{";
					SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_REAL, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, numB1, 0, passStr);
					passStr += "\n\t}\n";
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
					passStr += "\n\tif(rw)\n\t{";
					nextPass->SweepRegs(SR_READ, fwd, outInterleaved, outStride, SR_COMP_REAL, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, nextPass->GetNumB1(), 0, passStr);
					passStr += "\n\t}\n";
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
					passStr += "\n\tif(rw)\n\t{";
					SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_IMAG, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, numB1, 0, passStr);
					passStr += "\n\t}\n";
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
					passStr += "\n\tif(rw)\n\t{";
					nextPass->SweepRegs(SR_READ, fwd, outInterleaved, outStride, SR_COMP_IMAG, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, nextPass->GetNumB1(), 0, passStr);
					passStr += "\n\t}\n";
					passStr += "\n\n\tbarrier(CLK_LOCAL_MEM_FENCE);\n";
				}
			}
			else
			{
				if (fft_doPostCallback && outInterleaved)
				{
					passStr += "\n\t"; passStr += regB2Type; passStr += " tempC;";
				}
				passStr += "\n\tif(rw)\n\t{";
				SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, bufferOutRe, bufferOutIm, "outOffset", 1, numB1, 0, passStr);
				SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, bufferOutRe, bufferOutIm, "outOffset", 2, numB2, numB1, passStr);
				SweepRegs(SR_WRITE, fwd, outInterleaved, outStride, SR_COMP_BOTH, scale, false, bufferOutRe, bufferOutIm, "outOffset", 4, numB4, 2*numB2 + numB1, passStr);
				passStr += "\n\t}\n";
			}


			passStr += "\n}\n\n";
		}
    };

    // FFT kernel
    template <Precision PR>
    class Kernel
    {
        size_t length;							// Length of FFT
        size_t workGroupSize;					// Work group size
		size_t cnPerWI;							// complex numbers per work-item

		size_t numTrans;						// Number of transforms per work-group
		size_t workGroupSizePerTrans;			// Work group subdivision per transform
		size_t numPasses;						// Number of FFT passes
        std::vector<size_t> radices;			// Base radix at each pass
        std::vector<Pass<PR> > passes;			// Array of pass objects

		bool halfLds;							// LDS used to store one component (either real or imaginary) at a time
												// for passing intermediate data between the passes, if this is set
												// then each pass-function should accept same set of registers

		bool linearRegs;						// scalar registers

		// Future optimization ideas
		// bool limitRegs;							// TODO: Incrementally write to LDS, thereby using same set of registers for more than 1 butterflies
		// bool combineReadTwMul;					// TODO: Combine reading into registers and Twiddle multiply

		bool r2c2r;								// real to complex or complex to real transform
		bool r2c, c2r;
		bool rcFull;
		bool rcSimple;

		bool blockCompute;						// When we have to compute FFT in blocks (either read or write is along columns)
		BlockComputeType blockComputeType;
		size_t blockWidth, blockWGS, blockLDS;

		bool realSpecial;

		const FFTKernelGenKeyParams params;		// key params


		inline std::string IterRegs(const std::string &pfx, bool initComma = true)
		{
			std::string str = "";

			if(linearRegs)
			{
				if(initComma) str += ", ";

				for(size_t i=0; i<cnPerWI; i++)
				{
					if(i != 0) str += ", ";
					str += pfx; str += "R";
					str += SztToStr(i);
				}
			}

			return str;
		}

		inline bool IsGroupedReadWritePossible()
		{
			bool possible = true;
			const size_t *iStride, *oStride;

			if(r2c2r)
				return false;

			if(realSpecial)
				return false;

			if(params.fft_placeness == CLFFT_INPLACE)
			{
				iStride = oStride = params.fft_inStride;
			}
			else
			{
				iStride = params.fft_inStride;
				oStride = params.fft_outStride;
			}

			for(size_t i=1; i < params.fft_DataDim; i++)
			{
				if(iStride[i] % 2) { possible = false; break; }
				if(oStride[i] % 2) { possible = false; break; }
			}

			return possible;
		}

		inline std::string OffsetCalcBlock(const std::string &off, bool input = true)
		{
			std::string str;

			const size_t *pStride = input ? params.fft_inStride : params.fft_outStride;

			str += "\t"; str += off; str += " = ";
			std::string nextBatch = "batch";
			for(size_t i=(params.fft_DataDim - 1); i>2; i--)
			{
				size_t currentLength = 1;
				for(int j=2; j<i; j++) currentLength *= params.fft_N[j];
				currentLength *= (params.fft_N[1]/blockWidth);

				str += "("; str += nextBatch; str += "/"; str += SztToStr(currentLength);
				str += ")*"; str += SztToStr(pStride[i]); str += " + ";

				nextBatch = "(" + nextBatch + "%" + SztToStr(currentLength) + ")";
			}

			str += "("; str += nextBatch; str += "/"; str += SztToStr(params.fft_N[1]/blockWidth);
			str += ")*"; str += SztToStr(pStride[2]); str += " + ("; str += nextBatch;
			str += "%"; str += SztToStr(params.fft_N[1]/blockWidth); str += ")*";
			if( (input && (blockComputeType == BCT_R2C)) || (!input && (blockComputeType == BCT_C2R)) )
				str += SztToStr(blockWidth*length);
			else
				str += SztToStr(blockWidth);
			str += ";\n";

			return str;
		}

		inline std::string OffsetCalc(const std::string &off, bool input = true, bool rc_second_index = false)
		{
			std::string str;

			const size_t *pStride = input ? params.fft_inStride : params.fft_outStride;

			std::string batch;
			if(r2c2r && !rcSimple)
			{
				batch += "(batch*"; batch += SztToStr(2*numTrans);
				if(rc_second_index) batch += " + 1";
				else				batch += " + 0";

				if(numTrans != 1)	{ batch += " + 2*(me/"; batch += SztToStr(workGroupSizePerTrans); batch += "))"; }
				else				{ batch += ")"; }
			}
			else
			{
				if(numTrans == 1)	{	batch += "batch"; }
				else				{	batch += "(batch*"; batch += SztToStr(numTrans);
										batch += " + (me/"; batch += SztToStr(workGroupSizePerTrans); batch += "))"; }
			}

			str += "\t"; str += off; str += " = ";
			std::string nextBatch = batch;
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
        Kernel( const FFTKernelGenKeyParams &paramsVal) :
			params(paramsVal), r2c2r(false)

        {
			length = params.fft_N[0];
			workGroupSize = params.fft_SIMD;
			numTrans = (workGroupSize * params.fft_R) / length;

			r2c = false;
			c2r = false;
			// Check if it is R2C or C2R transform
			if(params.fft_inputLayout == CLFFT_REAL)  r2c = true;
			if(params.fft_outputLayout == CLFFT_REAL) c2r = true;
			r2c2r = (r2c || c2r);

			if(r2c)
			{
				rcFull = (	(params.fft_outputLayout == CLFFT_COMPLEX_INTERLEAVED) ||
							(params.fft_outputLayout == CLFFT_COMPLEX_PLANAR) ) ? true : false;
			}
			if(c2r)
			{
				rcFull = (	(params.fft_inputLayout  == CLFFT_COMPLEX_INTERLEAVED) ||
							(params.fft_inputLayout  == CLFFT_COMPLEX_PLANAR) ) ? true : false;
			}

			rcSimple = params.fft_RCsimple;

			halfLds = true;
			linearRegs = true;

			realSpecial = params.fft_realSpecial;

			blockCompute = params.blockCompute;
			blockComputeType = params.blockComputeType;
			// Make sure we can utilize all Lds if we are going to
			// use blocked columns to compute FFTs
			if(blockCompute)
			{
				assert(length <= 256);  // 256 parameter comes from prototype experiments
										// largest length at which block column possible given 32KB LDS limit
										// if LDS limit is different this number need to be changed appropriately
				halfLds = false;
				linearRegs = true;
			}

			assert( ((length*numTrans)%workGroupSize) == 0 );
			cnPerWI = (numTrans * length) / workGroupSize;
			workGroupSizePerTrans = workGroupSize/numTrans;

			// !!!! IMPORTANT !!!! Keep these assertions unchanged, algorithm depend on these to be true
			assert( (cnPerWI * workGroupSize) == (numTrans * length) );
			assert( cnPerWI <= length ); // Don't do more than 1 fft per work-item

			// Breakdown into passes

			size_t LS = 1;
			size_t L;
			size_t R = length;
			size_t pid = 0;

			// See if we can get radices from the lookup table
			const size_t *pRadices = NULL;
			size_t nPasses;
			KernelCoreSpecs<PR> kcs;
			kcs.GetRadices(length, nPasses, pRadices);
			if((params.fft_MaxWorkGroupSize >= 256) && (pRadices != NULL))
			{
				for(size_t i=0; i<nPasses; i++)
				{
					size_t rad = pRadices[i];
					L = LS * rad;
					R /= rad;

					radices.push_back(rad);
					passes.push_back(Pass<PR>(i, length, rad, cnPerWI, L, LS, R, linearRegs, halfLds, r2c, c2r, rcFull, rcSimple, realSpecial));
					
					//Pass precallback information to Pass object if its the first pass. 
					//This will be used in single kernel transforms
					if (params.fft_hasPreCallback && i == 0 && !params.blockCompute)
					{
						passes[0].SetPrecallback(params.fft_hasPreCallback, params.fft_preCallback);
					}

					//Pass post-callback information to Pass object if its the last pass. 
					//This will be used in single kernel transforms
					if (params.fft_hasPostCallback && i == (nPasses - 1) && !params.blockCompute)
					{
						passes[i].SetPostcallback(params.fft_hasPostCallback, params.fft_postCallback);
					}

					LS *= rad;
				}
				assert(R == 1); // this has to be true for correct radix composition of the length
				numPasses = nPasses;
			}
			else
			{
				// Possible radices
				size_t cRad[] = {13,11,10,8,7,6,5,4,3,2,1}; // Must be in descending order
				size_t cRadSize = (sizeof(cRad)/sizeof(cRad[0]));

				// Generate the radix and pass objects
				while(true)
				{
					size_t rad;

					assert(cRadSize >= 1);

					// Picks the radices in descending order (biggest radix first)
					for(size_t r=0; r<cRadSize; r++)
					{
						rad = cRad[r];

						if((rad > cnPerWI) || (cnPerWI%rad))
							continue;

						if(!(R % rad))
							break;
					}

					assert((cnPerWI%rad) == 0);

					L = LS * rad;
					R /= rad;

					radices.push_back(rad);
					passes.push_back(Pass<PR>(pid, length, rad, cnPerWI, L, LS, R, linearRegs, halfLds, r2c, c2r, rcFull, rcSimple, realSpecial));

					//Pass precallback information to Pass object if its the first pass. 
					//This will be used in single kernel transforms
					if (pid == 0 && params.fft_hasPreCallback)
					{
						passes[0].SetPrecallback(params.fft_hasPreCallback, params.fft_preCallback);
					}

					pid++;
					LS *= rad;

					assert(R >= 1);
					if(R == 1)
						break;
				}
				numPasses = pid;

				//Pass post-callback information to Pass object if its the last pass. 
				//This will be used in single kernel transforms
				if (params.fft_hasPostCallback)
				{
					passes[numPasses - 1].SetPostcallback(params.fft_hasPostCallback, params.fft_postCallback);
				}
			}

			assert(numPasses == passes.size());
			assert(numPasses == radices.size());

#ifdef PARMETERS_TO_BE_READ

			ParamRead pr;
			ReadParameterFile(pr);

			radices.clear();
			passes.clear();

			radices = pr.radices;
			numPasses = radices.size();

			LS = 1;
			R = length;
			for(size_t i=0; i<numPasses; i++)
			{
				size_t rad = radices[i];
				L = LS * rad;
				R /= rad;

				passes.push_back(Pass<PR>(i, length, rad, cnPerWI, L, LS, R, linearRegs));

				LS *= rad;
			}
			assert(R == 1);
#endif

			// Grouping read/writes ok?
			bool grp = IsGroupedReadWritePossible();
			for(size_t i=0; i < numPasses; i++)
				passes[i].SetGrouping(grp);

			// Store the next pass-object pointers
			if(numPasses > 1)
				for(size_t i=0; i < (numPasses - 1); i++)
					passes[i].SetNextPass(&passes[i+1]);


			if(blockCompute)
			{
				blockWidth = BlockSizes::BlockWidth(length);
				blockWGS = BlockSizes::BlockWorkGroupSize(length);
				blockLDS = BlockSizes::BlockLdsSize(length);
			}
			else
			{
				blockWidth = blockWGS = blockLDS = 0;
			}
		}

		class BlockSizes
		{
		public:
			enum ValType
			{
				BS_VT_WGS,
				BS_VT_BWD,
				BS_VT_LDS,
			};

			static size_t BlockLdsSize(size_t N) { return GetValue(N, BS_VT_LDS); }
			static size_t BlockWidth(size_t N) { return GetValue(N, BS_VT_BWD); }
			static size_t BlockWorkGroupSize(size_t N) { return GetValue(N, BS_VT_WGS); }

		private:

			static size_t GetValue(size_t N, ValType vt)
			{
				size_t wgs; // preferred work group size
				size_t bwd; // block width to be used
				size_t lds; // LDS size to be used for the block


				KernelCoreSpecs<PR> kcs;
				size_t t_wgs, t_nt;
				kcs.GetWGSAndNT(N, t_wgs, t_nt);

				switch(N)
				{
				case 256:	bwd = 8/PrecisionWidth<PR>();   wgs = (bwd > t_nt) ? 256 : t_wgs; break;
				case 128:	bwd = 8/PrecisionWidth<PR>();   wgs = (bwd > t_nt) ? 128 : t_wgs; break;
				case 64:	bwd = 16/PrecisionWidth<PR>();  wgs = (bwd > t_nt) ? 128 : t_wgs; break;
				case 32:	bwd = 32/PrecisionWidth<PR>();  wgs = (bwd > t_nt) ? 64  : t_wgs; break;
				case 16:	bwd = 64/PrecisionWidth<PR>();  wgs = (bwd > t_nt) ? 64  : t_wgs; break;
				case 8:		bwd = 128/PrecisionWidth<PR>(); wgs = (bwd > t_nt) ? 64  : t_wgs; break;
				default:	assert(false);
				}

				// block width cannot be less than numTrans, math in other parts of code depend on this assumption
				assert(bwd >= t_nt);

				lds = N*bwd;

				switch(vt)
				{
				case BS_VT_WGS: return wgs;
				case BS_VT_BWD: return bwd;
				case BS_VT_LDS: return lds;
				default: assert(false); return 0;
				}
			}
		};

        void GenerateKernel(std::string &str, cl_device_id Dev_ID)
		{
			std::string twType = RegBaseType<PR>(2);
			std::string rType  = RegBaseType<PR>(1);
			std::string r2Type  = RegBaseType<PR>(2);

			bool inInterleaved;	 // Input is interleaved format
			bool outInterleaved; // Output is interleaved format
			inInterleaved  = (	(params.fft_inputLayout == CLFFT_COMPLEX_INTERLEAVED) ||
								(params.fft_inputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;
			outInterleaved = (	(params.fft_outputLayout == CLFFT_COMPLEX_INTERLEAVED) ||
								(params.fft_outputLayout == CLFFT_HERMITIAN_INTERLEAVED) ) ? true : false;

			// use interleaved LDS when halfLds constraint absent
			bool ldsInterleaved = inInterleaved || outInterleaved;
			ldsInterleaved = halfLds ? false : ldsInterleaved;
			ldsInterleaved = blockCompute ? true : ldsInterleaved;

			bool inReal;  // Input is real format
			bool outReal; // Output is real format
			inReal  = (params.fft_inputLayout == CLFFT_REAL) ? true : false;
			outReal = (params.fft_outputLayout == CLFFT_REAL) ? true : false;

			size_t large1D = 0;
			if(params.fft_realSpecial)
				large1D = params.fft_N[0] * params.fft_realSpecial_Nr;
			else
				large1D = params.fft_N[0] * params.fft_N[1];

			// Pragma
			str += ClPragma<PR>();

			// Twiddle table
			if(length > 1)
			{
				TwiddleTable twTable(length);

				str += "\n__constant ";
				str += twType; str += " ";
				str += TwTableName();
				str += "["; str += SztToStr(length-1); str += "] = {\n";
				twTable.GenerateTwiddleTable<PR>(radices, str);
				str += "};\n\n";
			}
			str += "\n";

			// twiddle factors for 1d-large 3-step algorithm
			if(params.fft_3StepTwiddle)
			{
				TwiddleTableLarge twLarge(large1D);
				twLarge.GenerateTwiddleTable<PR>(str);
			}

			std::string sfx = FloatSuffix<PR>();

			// Base type
			str += "#define fptype "; str += RegBaseType<PR>(1); str += "\n\n";

			// Vector type
			str += "#define fvect2 "; str += RegBaseType<PR>(2); str += "\n\n";

			bool cReg = linearRegs ? true : false;

			// Generate butterflies for all unique radices
			std::list<size_t> uradices;
			for(std::vector<size_t>::const_iterator r = radices.begin(); r != radices.end(); r++)
				uradices.push_back(*r);

			uradices.sort();
			uradices.unique();


			//constants
			if (length%8 == 0)
			{
				str += "#define C8Q  0.70710678118654752440084436210485"; str += sfx; str += "\n";
			}

			if (length % 5 == 0)
			{
				str += "#define C5QA 0.30901699437494742410229341718282"; str += sfx; str += "\n";
				str += "#define C5QB 0.95105651629515357211643933337938"; str += sfx; str += "\n";
				str += "#define C5QC 0.50000000000000000000000000000000"; str += sfx; str += "\n";
				str += "#define C5QD 0.58778525229247312916870595463907"; str += sfx; str += "\n";
				str += "#define C5QE 0.80901699437494742410229341718282"; str += sfx; str += "\n";
			}

			if (length % 3 == 0)
			{
				str += "#define C3QA 0.50000000000000000000000000000000"; str += sfx; str += "\n";
				str += "#define C3QB 0.86602540378443864676372317075294"; str += sfx; str += "\n";
			}

			if (length % 7 == 0)
			{
				str += "#define C7Q1 -1.16666666666666651863693004997913"; str += sfx; str += "\n";
				str += "#define C7Q2  0.79015646852540022404554065360571"; str += sfx; str += "\n";
				str += "#define C7Q3  0.05585426728964774240049351305970"; str += sfx; str += "\n";
				str += "#define C7Q4  0.73430220123575240531721419756650"; str += sfx; str += "\n";
				str += "#define C7Q5  0.44095855184409837868031445395900"; str += sfx; str += "\n";
				str += "#define C7Q6  0.34087293062393136944265847887436"; str += sfx; str += "\n";
				str += "#define C7Q7 -0.53396936033772524066165487965918"; str += sfx; str += "\n";
				str += "#define C7Q8  0.87484229096165666561546458979137"; str += sfx; str += "\n";
			}

			if (length % 11 == 0)
			{
				str += "#define b11_0 0.9898214418809327"; str += sfx; str += "\n";
				str += "#define b11_1 0.9594929736144973"; str += sfx; str += "\n";
				str += "#define b11_2 0.9189859472289947"; str += sfx; str += "\n";
				str += "#define b11_3 0.8767688310025893"; str += sfx; str += "\n";
				str += "#define b11_4 0.8308300260037728"; str += sfx; str += "\n";
				str += "#define b11_5 0.7784344533346518"; str += sfx; str += "\n";
				str += "#define b11_6 0.7153703234534297"; str += sfx; str += "\n";
				str += "#define b11_7 0.6343562706824244"; str += sfx; str += "\n";
				str += "#define b11_8 0.3425847256816375"; str += sfx; str += "\n";
				str += "#define b11_9 0.5211085581132027"; str += sfx; str += "\n";
			}

			if (length % 13 == 0)
			{
				str += "#define b13_0  0.9682872443619840"; str += sfx; str += "\n";
				str += "#define b13_1  0.9578059925946651"; str += sfx; str += "\n";
				str += "#define b13_2  0.8755023024091479"; str += sfx; str += "\n";
				str += "#define b13_3  0.8660254037844386"; str += sfx; str += "\n";
				str += "#define b13_4  0.8595425350987748"; str += sfx; str += "\n";
				str += "#define b13_5  0.8534800018598239"; str += sfx; str += "\n";
				str += "#define b13_6  0.7693388175729806"; str += sfx; str += "\n";
				str += "#define b13_7  0.6865583707817543"; str += sfx; str += "\n";
				str += "#define b13_8  0.6122646503767565"; str += sfx; str += "\n";
				str += "#define b13_9  0.6004772719326652"; str += sfx; str += "\n";
				str += "#define b13_10 0.5817047785105157"; str += sfx; str += "\n";
				str += "#define b13_11 0.5751407294740031"; str += sfx; str += "\n";
				str += "#define b13_12 0.5220263851612750"; str += sfx; str += "\n";
				str += "#define b13_13 0.5200285718888646"; str += sfx; str += "\n";
				str += "#define b13_14 0.5165207806234897"; str += sfx; str += "\n";
				str += "#define b13_15 0.5149187780863157"; str += sfx; str += "\n";
				str += "#define b13_16 0.5035370328637666"; str += sfx; str += "\n";
				str += "#define b13_17 0.5000000000000000"; str += sfx; str += "\n";
				str += "#define b13_18 0.3027756377319946"; str += sfx; str += "\n";
				str += "#define b13_19 0.3014792600477098"; str += sfx; str += "\n";
				str += "#define b13_20 0.3004626062886657"; str += sfx; str += "\n";
				str += "#define b13_21 0.2517685164318833"; str += sfx; str += "\n";
				str += "#define b13_22 0.2261094450357824"; str += sfx; str += "\n";
				str += "#define b13_23 0.0833333333333333"; str += sfx; str += "\n";
				str += "#define b13_24 0.0386329546443481"; str += sfx; str += "\n";
			}

			str += "\n";


			//If pre-callback is set for the plan
			std::string callbackstr;
			if (params.fft_hasPreCallback)
			{
				//Insert pre-callback function code at the beginning 
				callbackstr += params.fft_preCallback.funcstring;
				callbackstr += "\n\n";

				str += callbackstr;
			}

			//If post-callback is set for the plan
			if (params.fft_hasPostCallback)
			{
				//Insert post-callback function code 
				str += params.fft_postCallback.funcstring;
				str += "\n\n";
			}

			typename std::vector< Pass<PR> >::const_iterator p;
			if(length > 1)
			{
				for(std::list<size_t>::const_iterator r = uradices.begin(); r != uradices.end(); r++)
				{
					size_t rad = *r;
					p = passes.begin();
					while(p->GetRadix() != rad) p++;

					for(size_t d=0; d<2; d++)
					{
						bool fwd = d ? false : true;

						if(p->GetNumB1()) { Butterfly<PR> bfly(rad, 1, fwd, cReg); bfly.GenerateButterfly(str); str += "\n"; }
						if(p->GetNumB2()) { Butterfly<PR> bfly(rad, 2, fwd, cReg); bfly.GenerateButterfly(str); str += "\n"; }
						if(p->GetNumB4()) { Butterfly<PR> bfly(rad, 4, fwd, cReg); bfly.GenerateButterfly(str); str += "\n"; }
					}
				}
			}

			// Generate passes
			for(size_t d=0; d<2; d++)
			{
				bool fwd;

				if(r2c2r)
				{
					fwd = r2c;
				}
				else
				{
					fwd = d ? false : true;
				}

				double scale = fwd ? params.fft_fwdScale : params.fft_backScale;

				for(p = passes.begin(); p != passes.end(); p++)
				{
					double s = 1.0;
					size_t ins = 1, outs = 1;
					bool gIn = false, gOut = false;
					bool inIlvd = false, outIlvd = false;
					bool inRl = false, outRl = false;
					bool tw3Step = false;


					if(p == passes.begin() && params.fft_twiddleFront ) { tw3Step = params.fft_3StepTwiddle; }
					if((p+1) == passes.end())	{ s = scale; if(!params.fft_twiddleFront) tw3Step = params.fft_3StepTwiddle; }

					if(blockCompute && !r2c2r)
					{
						inIlvd = ldsInterleaved;
						outIlvd = ldsInterleaved;
					}
					else
					{
						if(p == passes.begin())		{ inIlvd  = inInterleaved;  inRl  = inReal;  gIn  = true; ins  = params.fft_inStride[0];  }
						if((p+1) == passes.end())	{ outIlvd = outInterleaved; outRl = outReal; gOut = true; outs = params.fft_outStride[0]; }

						if(p != passes.begin())		{ inIlvd = ldsInterleaved; }
						if((p+1) != passes.end())	{ outIlvd = ldsInterleaved; }
					}

					p->GeneratePass(fwd, str, tw3Step, params.fft_twiddleFront, inIlvd, outIlvd, inRl, outRl, ins, outs, s, gIn, gOut);
				}

				// if real transform we do only 1 direction
				if(r2c2r)
					break;
			}



			// TODO : address this kludge
			str += " typedef union  { uint u; int i; } cb_t;\n\n";

			for(size_t d=0; d<2; d++)
			{
				bool fwd;

				if(r2c2r)
				{
					fwd = inReal ? true : false;
				}
				else
				{
					fwd = d ? false : true;
				}

				// FFT kernel begin
				// Function attribute
				str += "__kernel __attribute__((reqd_work_group_size (";
				if(blockCompute)	str += SztToStr(blockWGS);
				else				str += SztToStr(workGroupSize);
				str += ",1,1)))\nvoid ";

				// Function name
				if(fwd) str += "fft_fwd";
				else	str += "fft_back";
				str += "(";

        // TODO : address this kludge
        size_t SizeParam_ret = 0;
        clGetDeviceInfo(Dev_ID, CL_DEVICE_VENDOR, 0, NULL, &SizeParam_ret);
        char* nameVendor = new char[SizeParam_ret];
        clGetDeviceInfo(Dev_ID, CL_DEVICE_VENDOR, SizeParam_ret, nameVendor, NULL);

        //nv compiler doesn't support __constant kernel argument
        if (strncmp(nameVendor, "NVIDIA",6)!=0)
          str += "__constant cb_t *cb __attribute__((max_constant_size(32))), ";
        else
          str += "__global cb_t *cb, ";

        delete [] nameVendor;

		//If plan has pre/post callback
		callbackstr.clear();
		bool hasCallback = params.fft_hasPreCallback || params.fft_hasPostCallback;

		if (hasCallback)
		{
		if (params.fft_hasPreCallback)
		{
				callbackstr += ", __global void* pre_userdata";
			}
			if (params.fft_hasPostCallback)
			{
				callbackstr += ", __global void* post_userdata";
			}

			if (params.fft_preCallback.localMemSize > 0 || params.fft_postCallback.localMemSize > 0)
			{
				callbackstr += ", __local void* localmem";
			}
		}

				// Function attributes
				if(params.fft_placeness == CLFFT_INPLACE)
				{
					if(r2c2r)
					{
						if(outInterleaved)
						{
							str += "__global "; str += r2Type; str += " * restrict gb";
						}
						else
						{
							str += "__global "; str += rType; str += " * restrict gb";
						}

						//If plan has callback
						if (hasCallback)
						{
							str += callbackstr;
						}

						str += ")\n";
					}
					else
					{
						assert(inInterleaved == outInterleaved);
						assert(params.fft_inStride[1] == params.fft_outStride[1]);
						assert(params.fft_inStride[0] == params.fft_outStride[0]);

						if(inInterleaved)
						{
							str += "__global "; str += r2Type; str += " * restrict gb";

							//If plan has callback
							if (hasCallback)
							{
								str += callbackstr;
							}
							
							str += ")\n";
						}
						else
						{
							str += "__global "; str += rType; str += " * restrict gbRe, ";
							str += "__global "; str += rType; str += " * restrict gbIm";

							//If plan has callback
							if (hasCallback)
							{
								str += callbackstr;
							}

							str += ")\n";
						}
					}
				}
				else
				{
					if(r2c2r)
					{
						if(inInterleaved)
						{
							str += "__global "; str += r2Type; str += " * restrict gbIn, ";
						}
						else if(inReal)
						{
							str += "__global "; str += rType; str += " * restrict gbIn, ";
						}
						else
						{
							str += "__global const "; str += rType; str += " * restrict gbInRe, ";
							str += "__global const "; str += rType; str += " * restrict gbInIm, ";
						}

						if(outInterleaved)
						{
							str += "__global "; str += r2Type; str += " * restrict gbOut";
						}
						else if(outReal)
						{
							str += "__global "; str += rType; str += " * restrict gbOut";
						}
						else
						{
							str += "__global "; str += rType; str += " * restrict gbOutRe, ";
							str += "__global "; str += rType; str += " * restrict gbOutIm";
						}

						//If plan has callback
						if (hasCallback)
						{
							str += callbackstr;
						}

						str += ")\n";
					}
					else
					{
						if(inInterleaved)
						{
							str += "__global const "; str += r2Type; str += " * restrict gbIn, ";
						}
						else
						{
							str += "__global const "; str += rType; str += " * restrict gbInRe, ";
							str += "__global const "; str += rType; str += " * restrict gbInIm, ";
						}

						if(outInterleaved)
						{
							str += "__global "; str += r2Type; str += " * restrict gbOut";
						}
						else
						{
							str += "__global "; str += rType; str += " * restrict gbOutRe, ";
							str += "__global "; str += rType; str += " * restrict gbOutIm";
						}

						//If plan has callback
						if (hasCallback)
						{
							str += callbackstr;
						}

						str += ")\n";
					}
				}

				str += "{\n";

				// Initialize
				str += "\t";
				str += "uint me = get_local_id(0);\n\t";
				str += "uint batch = get_group_id(0);";
				str += "\n";



				// Allocate LDS
				if(blockCompute)
				{
					str += "\n\t"; str += "__local "; str += r2Type; str += " lds[";
					str += SztToStr(blockLDS); str += "];\n";
				}
				else
				{
					size_t ldsSize = halfLds ? length*numTrans : 2*length*numTrans;
					ldsSize = ldsInterleaved ? ldsSize/2 : ldsSize;

					if(numPasses > 1)
					{
						str += "\n\t";
						str += "__local "; str += ldsInterleaved ? r2Type: rType; str += " lds[";
						str += SztToStr(ldsSize); str += "];\n";
					}
				}

				// Declare memory pointers
				str += "\n\t";
				if(r2c2r)
				{
					str += "uint iOffset;\n\t";
					str += "uint oOffset;\n\n\t";
					if(!rcSimple)
					{
						str += "uint iOffset2;\n\t";
						str += "uint oOffset2;\n\n\t";
					}

					if (!params.fft_hasPreCallback)	
					{ 
						if(inInterleaved)
						{
							if(!rcSimple)	{	str += "__global "; str += r2Type; str += " *lwbIn2;\n\t"; }
												str += "__global "; str += r2Type; str += " *lwbIn;\n\t";  
						}
						else if(inReal)
						{
							if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbIn2;\n\t"; }
												str += "__global "; str += rType; str += " *lwbIn;\n\t";

						}
						else
						{
							if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbInRe2;\n\t"; }
							if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbInIm2;\n\t"; }
							  
												str += "__global "; str += rType; str += " *lwbInRe;\n\t"; 
												str += "__global "; str += rType; str += " *lwbInIm;\n\t"; 
							
						}
					}

					if(outInterleaved)
					{
						if (!params.fft_hasPostCallback)	
						{ 
						if(!rcSimple)	{	str += "__global "; str += r2Type; str += " *lwbOut2;\n\t"; }
												str += "__global "; str += r2Type; str += " *lwbOut;\n";
						}
					}
					else if(outReal)
					{
						if (!params.fft_hasPostCallback)	
						{
						if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbOut2;\n\t"; }
												str += "__global "; str += rType; str += " *lwbOut;\n";
						}
					}
					else
					{
						if (!params.fft_hasPostCallback)	
						{
						if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbOutRe2;\n\t"; }
						if(!rcSimple)	{	str += "__global "; str += rType; str += " *lwbOutIm2;\n\t"; }
											str += "__global "; str += rType; str += " *lwbOutRe;\n\t";
												str += "__global "; str += rType; str += " *lwbOutIm;\n";
						}
					}
					str += "\n";
				}
				else
				{
					if(params.fft_placeness == CLFFT_INPLACE)
					{
						str += "uint ioOffset;\n\t";

						//Skip if callback is set 
						if (!params.fft_hasPreCallback || !params.fft_hasPostCallback)
						{
							if(inInterleaved)
							{
								str += "__global "; str += r2Type; str += " *lwb;\n";
							}
							else
							{
								str += "__global "; str += rType; str += " *lwbRe;\n\t";
								str += "__global "; str += rType; str += " *lwbIm;\n";
							}
						}
						str += "\n";
					}
					else
					{
						str += "uint iOffset;\n\t";
						str += "uint oOffset;\n\t";

						//Skip if precallback is set 
						if (!(params.fft_hasPreCallback))
						{
							if(inInterleaved)
							{
								str += "__global "; str += r2Type; str += " *lwbIn;\n\t";
							}
							else
							{
								str += "__global "; str += rType; str += " *lwbInRe;\n\t";
								str += "__global "; str += rType; str += " *lwbInIm;\n\t";
							}
						}

						//Skip if postcallback is set 
						if (!params.fft_hasPostCallback)
						{
						if(outInterleaved)
						{
								str += "__global "; str += r2Type; str += " *lwbOut;\n";
						}
						else
						{
							str += "__global "; str += rType; str += " *lwbOutRe;\n\t";
								str += "__global "; str += rType; str += " *lwbOutIm;\n";
							}
						}
						str += "\n";
					}
				}

				// Setup registers if needed
				if(linearRegs)
				{
					str += "\t"; str += RegBaseType<PR>(2);
					str += " "; str += IterRegs("", false);
					str += ";\n\n";
				}

				// Calculate total transform count
				std::string totalBatch = "(";
				size_t i = 0;
				while(i < (params.fft_DataDim - 2))
				{
					totalBatch += SztToStr(params.fft_N[i+1]); totalBatch += " * ";
					i++;
				}
				totalBatch += "cb[0].u)";

				// Conditional read-write ('rw') for arbitrary batch number
				if(r2c2r && !rcSimple)
				{
					str += "\tuint this = "; str += totalBatch; str += " - batch*";
					str +=  SztToStr(2*numTrans); str += ";\n";
					str += "\tuint rw = (me < ((this+1)/2)*"; str += SztToStr(workGroupSizePerTrans);
					str += ") ? (this - 2*(me/"; str += SztToStr(workGroupSizePerTrans); str += ")) : 0;\n\n";
				}
				else
				{
					if( (numTrans > 1) && !blockCompute )
					{
						str += "\tuint rw = (me < ("; str += totalBatch;
						str += " - batch*"; str += SztToStr(numTrans); str += ")*";
						str += SztToStr(workGroupSizePerTrans); str += ") ? 1 : 0;\n\n";
					}
					else
					{
						str += "\tuint rw = 1;\n\n";
					}
				}

				// Transform index for 3-step twiddles
				if(params.fft_3StepTwiddle && !blockCompute)
				{
					if(numTrans == 1)
					{
						str += "\tuint b = batch%";
					}
					else
					{
						str += "\tuint b = (batch*"; str += SztToStr(numTrans); str += " + (me/";
						str += SztToStr(workGroupSizePerTrans); str += "))%";
					}

					str += SztToStr(params.fft_N[1]); str += ";\n\n";

					if(params.fft_realSpecial)
					{
						str += "\tuint bt = b;\n\n";
					}
				}
				else
				{
					str += "\tuint b = 0;\n\n";
				}

				// Setup memory pointers
				if(r2c2r)
				{
					str += OffsetCalc("iOffset", true);
					str += OffsetCalc("oOffset", false);
					if(!rcSimple) { str += OffsetCalc("iOffset2",  true, true); }
					if(!rcSimple) { str += OffsetCalc("oOffset2", false, true); }

					str += "\n\t";
					if(params.fft_placeness == CLFFT_INPLACE)
					{
						if(!params.fft_hasPreCallback)
						{
							if(inInterleaved)
							{
								if(!rcSimple) {	str += "lwbIn2 = (__global "; str += r2Type; str += " *)gb + iOffset2;\n\t"; }
												str += "lwbIn  = (__global "; str += r2Type; str += " *)gb + iOffset;\n\t"; 
							}
							else
							{
								if(!rcSimple) {	str += "lwbIn2 = (__global "; str += rType; str += " *)gb + iOffset2;\n\t"; }
												str += "lwbIn  = (__global "; str += rType; str += " *)gb + iOffset;\n\t";

							}
						}

						if(!params.fft_hasPostCallback)
						{
						if(!rcSimple) {	str += "lwbOut2 = gb + oOffset2;\n\t"; }
											str += "lwbOut = gb + oOffset;\n";
						}
						str += "\n";
					}
					else
					{
						if (!params.fft_hasPreCallback)
						{
							if(inInterleaved || inReal)
							{
								if(!rcSimple) {	str += "lwbIn2 = gbIn + iOffset2;\n\t"; }
								str += "lwbIn = gbIn + iOffset;\n\t"; 
							}
							else
							{
								if(!rcSimple) {	str += "lwbInRe2 = gbInRe + iOffset2;\n\t"; }
								if(!rcSimple) {	str += "lwbInIm2 = gbInIm + iOffset2;\n\t"; }
												str += "lwbInRe = gbInRe + iOffset;\n\t";
												str += "lwbInIm = gbInIm + iOffset;\n\t";
							}
						}

						if (!params.fft_hasPostCallback)
						{
						if(outInterleaved || outReal)
						{
							if(!rcSimple) {	str += "lwbOut2 = gbOut + oOffset2;\n\t"; }
												str += "lwbOut = gbOut + oOffset;\n";
						}
						else
						{
								
							if(!rcSimple) {	str += "lwbOutRe2 = gbOutRe + oOffset2;\n\t"; }
							if(!rcSimple) {	str += "lwbOutIm2 = gbOutIm + oOffset2;\n\t"; }
											str += "lwbOutRe = gbOutRe + oOffset;\n\t";
												str += "lwbOutIm = gbOutIm + oOffset;\n";
						}
					}
						str += "\n";
					}
				}
				else
				{
					if(params.fft_placeness == CLFFT_INPLACE)
					{
						if(blockCompute)
							str += OffsetCalcBlock("ioOffset", true);
						else
							str += OffsetCalc("ioOffset", true);

						str += "\t";

						//Skip if callback is set 
						if (!params.fft_hasPreCallback || !params.fft_hasPostCallback)
						{
							if(inInterleaved)
							{
								str += "lwb = gb + ioOffset;\n";
							}
							else
							{
								str += "lwbRe = gbRe + ioOffset;\n\t";
								str += "lwbIm = gbIm + ioOffset;\n";
							}
						}
						str += "\n";
					}
					else
					{
						if(blockCompute)
						{
							str += OffsetCalcBlock("iOffset", true);
							str += OffsetCalcBlock("oOffset", false);
						}
						else
						{
							str += OffsetCalc("iOffset", true);
							str += OffsetCalc("oOffset", false);
						}

						str += "\t";

						//Skip if precallback is set 
						if (!(params.fft_hasPreCallback))
						{
							if(inInterleaved)
							{
								str += "lwbIn = gbIn + iOffset;\n\t";
							}
							else
							{
								str += "lwbInRe = gbInRe + iOffset;\n\t";
								str += "lwbInIm = gbInIm + iOffset;\n\t";
							}
						}

						//Skip if postcallback is set 
						if (!params.fft_hasPostCallback)
						{
						if(outInterleaved)
						{
								str += "lwbOut = gbOut + oOffset;\n";
						}
						else
						{
							str += "lwbOutRe = gbOutRe + oOffset;\n\t";
								str += "lwbOutIm = gbOutIm + oOffset;\n";
							}
						}
						str += "\n";
					}
				}

				std::string inOffset;
				std::string outOffset;
				if (params.fft_placeness == CLFFT_INPLACE && !r2c2r)
				{
					inOffset += "ioOffset";
					outOffset += "ioOffset";
				}
				else
				{
					inOffset += "iOffset";
					outOffset += "oOffset";
				}
				
				// Read data into LDS for blocked access
				if(blockCompute)
				{

					size_t loopCount = (length * blockWidth)/blockWGS;
					
					if ((blockComputeType == BCT_C2C) && params.fft_hasPreCallback)
					{
						str += "\n\t"; str += r2Type; str += " retCallback;";
					}

					str += "\n\tfor(uint t=0; t<"; str += SztToStr(loopCount);
					str += "; t++)\n\t{\n";

					//get offset 
					std::string bufOffset;

					for(size_t c=0; c<2; c++)
					{
						std::string comp = "";
						std::string readBuf = (params.fft_placeness == CLFFT_INPLACE) ? "lwb" : "lwbIn";
						if(!inInterleaved) comp = c ? ".y" : ".x";
						if(!inInterleaved)
							readBuf = (params.fft_placeness == CLFFT_INPLACE) ? (c ? "lwbIm" : "lwbRe") : (c ? "lwbInIm" : "lwbInRe");

						if( (blockComputeType == BCT_C2C) || (blockComputeType == BCT_C2R) )
						{
							bufOffset.clear();
							bufOffset += "(me%"; bufOffset += SztToStr(blockWidth); bufOffset += ") + ";
							bufOffset += "(me/"; bufOffset+= SztToStr(blockWidth); bufOffset+= ")*"; bufOffset += SztToStr(params.fft_inStride[0]);
							bufOffset += " + t*"; bufOffset += SztToStr(params.fft_inStride[0]*blockWGS/blockWidth);

							if ((blockComputeType == BCT_C2C) && params.fft_hasPreCallback)
							{
								if (c == 0)
								{
									str += "\t\tretCallback = "; str += params.fft_preCallback.funcname; str += "(";
								
									if(inInterleaved)
									{
										str += (params.fft_placeness == CLFFT_INPLACE) ? "gb, " : "gbIn, ";
									}
									else
									{
										str += (params.fft_placeness == CLFFT_INPLACE) ? "gbRe, gbIm, " : "gbInRe, gbInIm, ";
									}

									str += inOffset; str += " + "; str += bufOffset; str += ", pre_userdata";
									str += (params.fft_preCallback.localMemSize > 0) ? str += ", localmem);\n" : ");\n";
								}

								str += "\t\tR0"; str+= comp; str+= " = retCallback"; str+= comp; str += ";\n";
							}
							else
							{
								str += "\t\tR0"; str+= comp; str+= " = "; 
								str += readBuf; str += "[";	str += bufOffset; str += "];\n";
							}
						}
						else
						{
							str += "\t\tR0"; str+= comp; str+= " = "; str += readBuf; str += "[me + t*"; str += SztToStr(blockWGS); str += "];\n";
						}


						if(inInterleaved) break;
					}

					if( (blockComputeType == BCT_C2C) || (blockComputeType == BCT_C2R) )
					{
						str += "\t\tlds[t*"; str += SztToStr(blockWGS/blockWidth); str += " + ";
						str += "(me%"; str+= SztToStr(blockWidth); str+= ")*"; str += SztToStr(length); str += " + ";
						str += "(me/"; str+= SztToStr(blockWidth); str+= ")] = R0;"; str +="\n";
					}
					else
					{
						str += "\t\tlds[t*"; str += SztToStr(blockWGS); str += " + me] = R0;"; str +="\n";
					}

					str += "\t}\n\n";		
					str += "\tbarrier(CLK_LOCAL_MEM_FENCE);\n\n";
				}


				// Set rw and 'me' per transform
				// rw string also contains 'b'
				std::string rw, me;

				if(r2c2r && !rcSimple)	rw = "rw, b, ";
				else					rw = ((numTrans > 1) || realSpecial) ? "rw, b, " : "1, b, ";

				if(numTrans > 1)	{ me += "me%"; me += SztToStr(workGroupSizePerTrans); me += ", "; }
				else				{ me += "me, "; }

				if(blockCompute) { me = "me%"; me += SztToStr(workGroupSizePerTrans); me += ", "; }

				// Buffer strings
				std::string inBuf, outBuf;
				if(r2c2r)
				{
					if(rcSimple)
					{
						if(inInterleaved || inReal)		inBuf  = params.fft_hasPreCallback ?  "gbIn, " : "lwbIn, ";
						else							inBuf  = "lwbInRe, lwbInIm, ";
						if(outInterleaved || outReal)	outBuf = params.fft_hasPostCallback ? "gbOut" : "lwbOut";
						else							outBuf = "lwbOutRe, lwbOutIm";
					}
					else
					{
						if(inInterleaved || inReal)		
						{
							if (!params.fft_hasPreCallback)
							{
								inBuf  = "lwbIn, lwbIn2, ";
							}
							else
							{
								if (params.fft_placeness == CLFFT_INPLACE) 
								{
									inBuf = "(__global "; inBuf += r2c ? rType : r2Type; inBuf += "*) gb, ";
									inBuf += "(__global "; inBuf += r2c ? rType : r2Type; inBuf += "*) gb, ";
								}
								else
								{
									inBuf  = "gbIn, gbIn, " ;
								}
							}	
						}
						else							inBuf  = (params.fft_hasPreCallback) ? "gbInRe, gbInRe, gbInIm, gbInIm, " : "lwbInRe, lwbInRe2, lwbInIm, lwbInIm2, ";

						if(outInterleaved || outReal)	outBuf = params.fft_hasPostCallback ? ((params.fft_placeness == CLFFT_INPLACE) ? "gb, gb" : "gbOut, gbOut") : "lwbOut, lwbOut2";
						else							outBuf = params.fft_hasPostCallback ? "gbOutRe, gbOutRe, gbOutIm, gbOutIm" : "lwbOutRe, lwbOutRe2, lwbOutIm, lwbOutIm2";
					}
				}
				else
				{
					if(params.fft_placeness == CLFFT_INPLACE)
					{
						if(inInterleaved)	
						{
							inBuf = params.fft_hasPreCallback ? "gb, " : "lwb, ";
							outBuf = params.fft_hasPostCallback ? "gb" : "lwb"; 
						}
						else	
						{ 
							inBuf = params.fft_hasPreCallback ? "gbRe, gbIm, " : "lwbRe, lwbIm, ";
							outBuf = params.fft_hasPostCallback ? "gbRe, gbIm" : "lwbRe, lwbIm"; 
						}
					}
					else
					{
						if(inInterleaved)	inBuf  = params.fft_hasPreCallback ? "gbIn, " : "lwbIn, ";
						else				inBuf  = params.fft_hasPreCallback ? "gbInRe, gbInIm, " : "lwbInRe, lwbInIm, ";
						if(outInterleaved)	outBuf = params.fft_hasPostCallback ? "gbOut" : "lwbOut";
						else				outBuf = params.fft_hasPostCallback ? "gbOutRe, gbOutIm" : "lwbOutRe, lwbOutIm";
					}
				}


				if(blockCompute)
				{
					str += "\n\tfor(uint t=0; t<"; str += SztToStr(blockWidth/(blockWGS/workGroupSizePerTrans));
					str += "; t++)\n\t{\n\n";

					inBuf = "lds, ";
					outBuf = "lds";

					if(params.fft_3StepTwiddle)
					{
						str += "\t\tb = (batch%"; str += SztToStr(params.fft_N[1]/blockWidth); str += ")*";
						str += SztToStr(blockWidth); str += " + t*"; str += SztToStr(blockWGS/workGroupSizePerTrans);
						str += " + (me/"; str += SztToStr(workGroupSizePerTrans); str += ");\n\n";
					}
				}

				if(realSpecial)
				{
					str += "\n\tfor(uint t=0; t<2; t++)\n\t{\n\n";
				}

				// Call passes
				if(numPasses == 1)
				{
					str += "\t";
					str += PassName(0, fwd);
					str += "("; str += rw; str += me;
					
					str += (params.fft_hasPreCallback) ? inOffset : "0";
					
					if (params.fft_hasPostCallback)
					{
						str += ", "; str += outOffset; str += ", ";
					}
					else
					{
					str += ", 0, ";
					}

					str += inBuf; str += outBuf;
					str += IterRegs("&");

					//If callback is set
					if (hasCallback)
					{
						//if pre-calback set 
					if (params.fft_hasPreCallback)
					{
							str += (r2c2r && !rcSimple) ?  ", iOffset2, pre_userdata" : ", pre_userdata";
						}

						//if post-calback set 
						if (params.fft_hasPostCallback)
						{
							if ((r2c || c2r) && !rcSimple) { str += ", "; str += outOffset; str += "2"; }

							str += ", post_userdata";
						}

						if (params.fft_preCallback.localMemSize > 0)
						{
							str += ", localmem";
						}
						if (params.fft_postCallback.localMemSize > 0)
						{
							//if precallback localmem also requested, send the localmem with the right offset
							if (params.fft_hasPreCallback && params.fft_preCallback.localMemSize > 0)
							{
								str += ", ((__local char *)localmem + "; str += SztToStr(params.fft_preCallback.localMemSize); str += ")";
							}
							else
							{
								str += ", localmem";
							}
						}
					}

					str += ");\n";
				}
				else
				{
					for(typename std::vector<Pass<PR> >::const_iterator p = passes.begin(); p != passes.end(); p++)
					{
						std::string exTab = "";
						if(blockCompute || realSpecial) exTab = "\t";

						str += exTab;
						str += "\t";
						str += PassName(p->GetPosition(), fwd);
						str += "(";

						std::string ldsOff;
						if(blockCompute)
						{
							ldsOff += "t*"; ldsOff += SztToStr(length*(blockWGS/workGroupSizePerTrans)); ldsOff += " + (me/";
							ldsOff += SztToStr(workGroupSizePerTrans); ldsOff += ")*"; ldsOff += SztToStr(length);
						}
						else
						{
							if(numTrans > 1)
							{
								ldsOff += "(me/"; ldsOff += SztToStr(workGroupSizePerTrans);
								ldsOff += ")*"; ldsOff += SztToStr(length);
							}
							else
							{
								ldsOff += "0";
							}
						}

						std::string ldsArgs;
						if(halfLds) { ldsArgs += "lds, lds"; }
						else		{	if(ldsInterleaved) { ldsArgs += "lds"; }
										else { ldsArgs += "lds, lds + "; ldsArgs += SztToStr(length*numTrans); } }

						str += rw;
						if(params.fft_realSpecial) str += "t, ";
						str += me;
						if(p == passes.begin()) // beginning pass
						{
							if (blockCompute)
							{
								str += ldsOff;
							}
							else
							{
								str += (params.fft_hasPreCallback) ? inOffset : "0";
							}
							str += ", ";
							str += ldsOff;
							str += ", ";
							str += inBuf;
							str += ldsArgs; str += IterRegs("&"); 
							
							//if precalback set, append additional arguments
							if (!blockCompute && params.fft_hasPreCallback)
							{
								str += (r2c2r && !rcSimple) ?  ", iOffset2, pre_userdata" : ", pre_userdata";

								if (params.fft_preCallback.localMemSize > 0)
								{
									str += ", localmem";
								}
							}

							str += ");\n";
							if(!halfLds) { str += exTab; str += "\tbarrier(CLK_LOCAL_MEM_FENCE);\n"; }
						}
						else if((p+1) == passes.end()) // ending pass
						{
							str += ldsOff;
							str += ", ";
							if (blockCompute)
							{
								str += ldsOff;
							}
							else
							{
								str += (params.fft_hasPostCallback) ? outOffset : "0";
							}
							str += ", ";
							str += ldsArgs; str += ", ";
							str += outBuf;

							str += IterRegs("&"); 
							
							if (!blockCompute && params.fft_hasPostCallback)
							{
								if ((c2r || r2c) && !rcSimple) { str += ", "; str += outOffset; str += "2"; }

								str += ", post_userdata";

								if (params.fft_postCallback.localMemSize > 0)
								{
									//if precallback localmem also requested, send the localmem with the right offset
									if (params.fft_hasPreCallback && params.fft_preCallback.localMemSize > 0)
									{
										str += ", ((__local char *)localmem + "; str += SztToStr(params.fft_preCallback.localMemSize); str += ")";
									}
									else
									{
										str += ", localmem";
									}
								}
							}
							str += ");\n";
							
							if (!halfLds) { str += exTab; str += "\tbarrier(CLK_LOCAL_MEM_FENCE);\n"; }
						}
						else // intermediate pass
						{
							str += ldsOff;
							str += ", ";
							str += ldsOff;
							str += ", ";
							str += ldsArgs; str += ", ";
							str += ldsArgs; str += IterRegs("&"); str += ");\n";
							if(!halfLds) { str += exTab; str += "\tbarrier(CLK_LOCAL_MEM_FENCE);\n"; }
						}
					}
				}

				if(realSpecial)
				{
					size_t Nt = 1 + length/2;
					str += 	"\n\t\tif( (bt == 0) || (2*bt == ";
					str += SztToStr(params.fft_realSpecial_Nr); str += ") ) { rw = 0; }\n";

					str += "\t\tlwbOut += ("; str += SztToStr(params.fft_realSpecial_Nr);
					str += " - 2*bt)*"; str += SztToStr(Nt); str += ";\n";
					str += "\t\tb = "; str += SztToStr(params.fft_realSpecial_Nr);
					str += " - b;\n\n";
				}

				if(blockCompute || realSpecial)
				{
					str += "\n\t}\n\n";
				}



				// Write data from LDS for blocked access
				if(blockCompute)
				{

					size_t loopCount = (length * blockWidth)/blockWGS;
					
					str += "\tbarrier(CLK_LOCAL_MEM_FENCE);\n\n";
					str += "\n\tfor(uint t=0; t<"; str += SztToStr(loopCount);
					str += "; t++)\n\t{\n";

					if( (blockComputeType == BCT_C2C) || (blockComputeType == BCT_R2C) )
					{
						str += "\t\tR0 = lds[t*"; str += SztToStr(blockWGS/blockWidth); str += " + ";
						str += "(me%"; str+= SztToStr(blockWidth); str+= ")*"; str += SztToStr(length); str += " + ";
						str += "(me/"; str+= SztToStr(blockWidth); str+= ")];"; str +="\n";
					}
					else
					{
						str += "\t\tR0 = lds[t*"; str += SztToStr(blockWGS); str += " + me];"; str +="\n";
					}

					for(size_t c=0; c<2; c++)
					{
						std::string comp = "";
						std::string writeBuf = (params.fft_placeness == CLFFT_INPLACE) ? "lwb" : "lwbOut";
						if(!outInterleaved) comp = c ? ".y" : ".x";
						if(!outInterleaved)
							writeBuf = (params.fft_placeness == CLFFT_INPLACE) ? (c ? "lwbIm" : "lwbRe") : (c ? "lwbOutIm" : "lwbOutRe");

						if( (blockComputeType == BCT_C2C) || (blockComputeType == BCT_R2C) )
						{
							if (blockComputeType == BCT_R2C && params.fft_hasPostCallback)
							{
								if (outInterleaved)
									writeBuf = (params.fft_placeness == CLFFT_INPLACE) ? "gb" : "gbOut";
								else
									writeBuf = (params.fft_placeness == CLFFT_INPLACE) ? "gbRe, gbIm" : "gbOutRe, gbOutIm";
								
								str += "\t\t"; str += params.fft_postCallback.funcname; str += "("; str += writeBuf; str += ", (";
								str += outOffset; str += " + (me%"; str+= SztToStr(blockWidth); str += ") + ";
								str += "(me/"; str+= SztToStr(blockWidth); str+= ")*"; str += SztToStr(params.fft_outStride[0]);
								str += " + t*"; str += SztToStr(params.fft_outStride[0]*blockWGS/blockWidth); 
								str += "), post_userdata, R0"; 
								if (!outInterleaved) str += ".x, R0.y";

								if (params.fft_postCallback.localMemSize > 0)
								{
									if (params.fft_hasPreCallback && params.fft_preCallback.localMemSize > 0)
									{
										str += ", (char *)(localmem + "; str += SztToStr(params.fft_preCallback.localMemSize); str += ")";
									}
									else
									{
										str += ", localmem";
									}
								}
								str += ");\n";

								//in the planar case, break from for loop since both real and imag components are handled
								//together in post-callback
								if (!outInterleaved) break;
							}
							else
							{
							str += "\t\t"; str += writeBuf; str += "[(me%"; str+= SztToStr(blockWidth); str += ") + ";
							str += "(me/"; str+= SztToStr(blockWidth); str+= ")*"; str += SztToStr(params.fft_outStride[0]);
							str += " + t*"; str += SztToStr(params.fft_outStride[0]*blockWGS/blockWidth); str += "] = R0"; str+= comp; str += ";\n";
						}
						}
						else
						{
							str += "\t\t"; str += writeBuf; str += "[me + t*"; str += SztToStr(blockWGS); str += "] = R0"; str+= comp; str += ";\n";
						}

						if(outInterleaved) break;
					}
					
					str += "\t}\n\n";		
				}
				


				str += "}\n\n";

				if(r2c2r)
					break;
			}
		}
    };
};

using namespace StockhamGenerator;

clfftStatus FFTGeneratedStockhamAction::initParams ()
{

    //    Query the devices in this context for their local memory sizes
    //    How we generate a kernel depends on the *minimum* LDS size for all devices.
    //
    const FFTEnvelope * pEnvelope = NULL;
    OPENCL_V(this->plan->GetEnvelope (& pEnvelope), _T("GetEnvelope failed"));
    BUG_CHECK (NULL != pEnvelope);

    // Remainder: params was properly cleared by its constructor
    //            clearing it again would destroy datasize and id!!
    this->signature.fft_precision    = this->plan->precision;
    this->signature.fft_placeness    = this->plan->placeness;
    this->signature.fft_inputLayout  = this->plan->inputLayout;
	this->signature.fft_MaxWorkGroupSize = this->plan->envelope.limit_WorkGroupSize;

    ARG_CHECK(this->plan->length.size()    > 0);
	ARG_CHECK(this->plan->inStride.size()  > 0);
    ARG_CHECK(this->plan->outStride.size() > 0);

    ARG_CHECK (this->plan->inStride.size() == this->plan->outStride.size())

	bool real_transform = ((this->plan->inputLayout == CLFFT_REAL) || (this->plan->outputLayout == CLFFT_REAL));

    if ( (CLFFT_INPLACE == this->plan->placeness) && (!real_transform) ) {
        //    If this is an in-place transform the
        //    input and output layout, dimensions and strides
        //    *MUST* be the same.
        //
        ARG_CHECK (this->plan->inputLayout == this->plan->outputLayout)
        this->signature.fft_outputLayout = this->plan->inputLayout;
        for (size_t u = this->plan->inStride.size(); u-- > 0; ) {
            ARG_CHECK (this->plan->inStride[u] == this->plan->outStride[u]);
        }
    } else {
        this->signature.fft_outputLayout = this->plan->outputLayout;
    }

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


	this->signature.fft_RCsimple = this->plan->RCsimple;

	this->signature.fft_realSpecial = this->plan->realSpecial;
	this->signature.fft_realSpecial_Nr = this->plan->realSpecial_Nr;

	this->signature.blockCompute = this->plan->blockCompute;
	this->signature.blockComputeType = this->plan->blockComputeType;

	this->signature.fft_twiddleFront = this->plan->twiddleFront;

	size_t wgs, nt;
#ifdef PARMETERS_TO_BE_READ
	ParamRead pr;
	ReadParameterFile(pr);
	wgs = pr.workGroupSize;
	nt = pr.numTransformsPerWg;
#else
	size_t t_wgs, t_nt;
	Precision pr = (this->signature.fft_precision == CLFFT_SINGLE) ? P_SINGLE : P_DOUBLE;
	switch(pr)
	{
	case P_SINGLE:
		{
			KernelCoreSpecs<P_SINGLE> kcs;
			kcs.GetWGSAndNT(this->signature.fft_N[0], t_wgs, t_nt);
			if(this->signature.blockCompute)
			{
				this->signature.blockSIMD = Kernel<P_SINGLE>::BlockSizes::BlockWorkGroupSize(this->signature.fft_N[0]);
				this->signature.blockLDS  = Kernel<P_SINGLE>::BlockSizes::BlockLdsSize(this->signature.fft_N[0]);
			}
		} break;
	case P_DOUBLE:
		{
			KernelCoreSpecs<P_DOUBLE> kcs;
			kcs.GetWGSAndNT(this->signature.fft_N[0], t_wgs, t_nt);
			if(this->signature.blockCompute)
			{
				this->signature.blockSIMD = Kernel<P_DOUBLE>::BlockSizes::BlockWorkGroupSize(this->signature.fft_N[0]);
				this->signature.blockLDS  = Kernel<P_DOUBLE>::BlockSizes::BlockLdsSize(this->signature.fft_N[0]);
			}
		} break;
	}

	if((t_wgs != 0) && (t_nt != 0) && (this->plan->envelope.limit_WorkGroupSize >= 256))
	{
		wgs = t_wgs;
		nt = t_nt;
	}
	else
		DetermineSizes(this->plan->envelope.limit_WorkGroupSize, this->signature.fft_N[0], wgs, nt, pr);
#endif

	assert((nt * this->signature.fft_N[0]) >= wgs);
	assert((nt * this->signature.fft_N[0])%wgs == 0);

	this->signature.fft_R = (nt * this->signature.fft_N[0])/wgs;
	this->signature.fft_SIMD = wgs;

	//Set pre-callback if specified
	if (this->plan->hasPreCallback)
	{
		this->signature.fft_hasPreCallback = true;
		this->signature.fft_preCallback = this->plan->preCallback;
	}

	//Set post-callback if specified
	if (this->plan->hasPostCallback)
	{
		this->signature.fft_hasPostCallback = true;
		this->signature.fft_postCallback = this->plan->postCallbackParam;
	}
	this->signature.limit_LocalMemSize = this->plan->envelope.limit_LocalMemSize;

    if (this->plan->large1D != 0) {
        ARG_CHECK (this->signature.fft_N[0] != 0)
        ARG_CHECK ((this->plan->large1D % this->signature.fft_N[0]) == 0)
        this->signature.fft_3StepTwiddle = true;
		if(!(this->plan->realSpecial))
			ARG_CHECK ( this->plan->large1D  == (this->signature.fft_N[1] * this->signature.fft_N[0]) );
    }

    this->signature.fft_fwdScale  = this->plan->forwardScale;
    this->signature.fft_backScale = this->plan->backwardScale;

    return CLFFT_SUCCESS;
}

clfftStatus FFTGeneratedStockhamAction::getWorkSizes (std::vector<size_t> & globalWS, std::vector<size_t> & localWS)
{
    //    How many complex numbers in the input mutl-dimensional array?
    //
    unsigned long long count = 1;
    for (unsigned u = 0; u < this->plan->length.size(); ++u) {
        count *= std::max<size_t> (1, this->plan->length[ u ]);
    }
    count *= this->plan->batchsize;


	if(this->signature.blockCompute)
	{
		count = DivRoundingUp<unsigned long long> (count, this->signature.blockLDS); 
		count = count * this->signature.blockSIMD; 

		globalWS.push_back( static_cast< size_t >( count ) );
		localWS.push_back( this->signature.blockSIMD );

		return    CLFFT_SUCCESS;
	}

    count = DivRoundingUp<unsigned long long> (count, this->signature.fft_R);      // count of WorkItems
    count = DivRoundingUp<unsigned long long> (count, this->signature.fft_SIMD);   // count of WorkGroups

	// for real transforms we only need half the work groups since we do twice the work in 1 work group
	if( !(this->signature.fft_RCsimple) && ((this->signature.fft_inputLayout == CLFFT_REAL) || (this->signature.fft_outputLayout == CLFFT_REAL)) )
		count = DivRoundingUp<unsigned long long> (count, 2);

    count = std::max<unsigned long long> (count, 1) * this->signature.fft_SIMD;
        // .. count of WorkItems, rounded up to next multiple of fft_SIMD.

	// 1 dimension work group size
	globalWS.push_back( static_cast< size_t >( count ) );

    localWS.push_back( this->signature.fft_SIMD );

    return    CLFFT_SUCCESS;
}

clfftStatus FFTPlan::GetMax1DLengthStockham (size_t * longest) const
{
	// TODO  The caller has already acquired the lock on *this
	//	However, we shouldn't depend on it.

	//	Query the devices in this context for their local memory sizes
	//	How large a kernel we can generate depends on the *minimum* LDS
	//	size for all devices.
	//
	const FFTEnvelope * pEnvelope = NULL;
	OPENCL_V(this->GetEnvelope (& pEnvelope), _T("GetEnvelope failed"));
	BUG_CHECK (NULL != pEnvelope);

	ARG_CHECK (NULL != longest)
	size_t LdsperElement = this->ElementSize();
	size_t result = pEnvelope->limit_LocalMemSize /
		(1 * LdsperElement);
	result = FloorPo2 (result);
	*longest = result;
	return CLFFT_SUCCESS;
}

clfftStatus FFTGeneratedStockhamAction::generateKernel(FFTRepo& fftRepo, const cl_command_queue commQueueFFT )
{
    cl_int status = CL_SUCCESS;
    cl_device_id Device = NULL;
    status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_DEVICE, sizeof(cl_device_id), &Device, NULL);
    OPENCL_V( status, _T( "clGetCommandQueueInfo failed" ) );

    cl_context QueueContext = NULL;
    status = clGetCommandQueueInfo(commQueueFFT, CL_QUEUE_CONTEXT, sizeof(cl_context), &QueueContext, NULL);
    OPENCL_V( status, _T( "clGetCommandQueueInfo failed" ) );

	std::string programCode;
	Precision pr = (this->signature.fft_precision == CLFFT_SINGLE) ? P_SINGLE : P_DOUBLE;
	switch(pr)
	{
	case P_SINGLE:
		{
			Kernel<P_SINGLE> kernel(this->signature);
			kernel.GenerateKernel(programCode, Device);
		} break;
	case P_DOUBLE:
		{
			Kernel<P_DOUBLE> kernel(this->signature);
			kernel.GenerateKernel(programCode, Device);
		} break;
	}

	//Requested local memory size by callback must not exceed the device LDS limits after factoring the LDS size required by main FFT kernel
	if ((this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0) || 
		(this->signature.fft_hasPostCallback && this->signature.fft_postCallback.localMemSize > 0))
	{
		bool validLDSSize = false;
		size_t requestedCallbackLDS = 0;

		if (this->signature.fft_hasPreCallback && this->signature.fft_preCallback.localMemSize > 0)
			requestedCallbackLDS = this->signature.fft_preCallback.localMemSize;
		if (this->signature.fft_hasPostCallback && this->signature.fft_postCallback.localMemSize > 0)
			requestedCallbackLDS += this->signature.fft_postCallback.localMemSize;

		if (this->plan->blockCompute)
		{
			validLDSSize = ((this->signature.blockLDS * this->plan->ElementSize()) +  requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;
		}
		else
		{
			size_t length = this->signature.fft_N[0];
			size_t workGroupSize = this->signature.fft_SIMD;
			size_t numTrans = (workGroupSize * this->signature.fft_R) / length;

			//TODO - Need to abstract this out. Repeating the same compute as in GenerateKernel. 
			// Set half lds only for power-of-2 problem sizes & interleaved data
			bool halfLds = ( (this->signature.fft_inputLayout == CLFFT_COMPLEX_INTERLEAVED) &&
						(this->signature.fft_outputLayout == CLFFT_COMPLEX_INTERLEAVED) ) ? true : false;
			halfLds = halfLds ? ((length & (length-1)) ? false : true) : false;

			// Set half lds for real transforms
			halfLds = ( (this->signature.fft_inputLayout == CLFFT_REAL) &&
						(this->signature.fft_outputLayout == CLFFT_REAL) ) ? true : halfLds;

			size_t ldsSize = halfLds ? length*numTrans : 2*length*numTrans;
			size_t elementSize = ((this->signature.fft_precision == CLFFT_DOUBLE) || (this->signature.fft_precision == CLFFT_DOUBLE_FAST)) ? sizeof(double) : sizeof(float);

			validLDSSize = ((ldsSize * elementSize) + requestedCallbackLDS) < this->plan->envelope.limit_LocalMemSize;
		}

		if(!validLDSSize)
		{
			fprintf(stderr, "Requested local memory size not available\n");
			return CLFFT_INVALID_ARG_VALUE;
		}
	}

#ifdef KERNEL_INTERJECT
	ReadKernelFromFile(programCode);
#endif

    OPENCL_V( fftRepo.setProgramCode( this->getGenerator(), this->getSignatureData(), programCode, Device, QueueContext ), _T( "fftRepo.setclString() failed!" ) );
    OPENCL_V( fftRepo.setProgramEntryPoints( this->getGenerator(), this->getSignatureData(), "fft_fwd", "fft_back", Device, QueueContext ), _T( "fftRepo.setProgramEntryPoint() failed!" ) );

    return CLFFT_SUCCESS;
}
