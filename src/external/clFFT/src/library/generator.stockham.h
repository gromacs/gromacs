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

#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif

#pragma once
#if !defined( AMD_CLFFT_generator_stockham_H )
#define AMD_CLFFT_generator_stockham_H
#include <stdio.h>
#include "private.h"
#include "repo.h"
#include "plan.h"

typedef union {
	cl_float f;
	cl_uint  u;
	cl_int   i;
} cb_t;

namespace StockhamGenerator
{
	// Precision
	enum Precision
	{
		P_SINGLE,
		P_DOUBLE,
	};

	template <Precision PR>
	inline size_t PrecisionWidth()
	{
		switch(PR)
		{
		case P_SINGLE:	return 1;
		case P_DOUBLE:	return 2;
		default:		assert(false); return 1;
		}
	}

	template <Precision PR>
	inline std::string ClPragma()
	{
		switch(PR)
		{
		case P_SINGLE:	return "";
		case P_DOUBLE:	return	"\n#ifdef cl_khr_fp64\n"
								"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n"
								"#else\n"
								"#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n"
								"#endif\n\n";
		default:		assert(false); return "";
		}
	}

	// Convert unsigned integers to string
	inline std::string SztToStr(size_t i)
	{
		std::stringstream ss;
		ss << i;
		return ss.str();
	}

	inline std::string FloatToStr(double f)
	{
		std::stringstream ss;
		ss.imbue(std::locale("C"));
		ss.precision(16);
		ss << std::scientific << f;
		return ss.str();
	}


	//	Find the smallest power of 2 that is >= n; return its power of 2 factor
	//	e.g., CeilPo2 (7) returns 3 : (2^3 >= 7)
	inline size_t CeilPo2 (size_t n)
	{
		size_t v = 1, t = 0;
		while(v < n)
		{
			v <<= 1;
			t++;
		}

		return t;
	}

	inline size_t FloorPo2 (size_t n)
	//	return the largest power of 2 that is <= n.
	//	e.g., FloorPo2 (7) returns 4.
	// *** TODO use x86 BSR instruction, using compiler intrinsics.
	{
		size_t tmp;
		while (0 != (tmp = n & (n-1)))
			n = tmp;
		return n;
	}

	typedef std::pair<std::string,std::string> stringpair;
	inline stringpair ComplexMul(const char *type, const char * a, const char * b, bool forward = true)
	{
		stringpair result;
		result.first = "(";
		result.first += type;
		result.first += ") ((";
		result.first += a;
		result.first += ".x * ";
		result.first += b;
		result.first += (forward ? ".x - " : ".x + ");
		result.first += a;
		result.first += ".y * ";
		result.first += b;
		result.first += ".y),";
		result.second = "(";
		result.second += a;
		result.second += ".y * ";
		result.second += b;
		result.second += (forward ? ".x + " : ".x - ");
		result.second += a;
		result.second += ".x * ";
		result.second += b;
		result.second += ".y))";
		return result;
	}


	// Register data base types
	template <Precision PR>
	inline std::string RegBaseType(size_t count)
	{
		switch(PR)
		{
		case P_SINGLE:
			switch(count)
			{
			case 1: return "float";
			case 2: return "float2";
			case 4: return "float4";
			default: assert(false); return "";
			}
			break;
		case P_DOUBLE:
			switch(count)
			{
			case 1: return "double";
			case 2: return "double2";
			case 4: return "double4";
			default: assert(false); return "";
			}
			break;
		default:
			assert(false); return "";
		}
	}

	template <Precision PR>
	inline std::string FloatSuffix()
	{
		// Suffix for constants
		std::string sfx;
		switch(PR)
		{
		case P_SINGLE: sfx = "f"; break;
		case P_DOUBLE: sfx = "";  break;
		default: assert(false);
		}

		return sfx;
	}

	inline std::string ButterflyName(size_t radix, size_t count, bool fwd)
	{
		std::string str;
		if(fwd) str += "Fwd";
		else	str += "Inv";
		str += "Rad"; str += SztToStr(radix);
		str += "B"; str += SztToStr(count);
		return str;
	}

	inline std::string PassName(size_t pos, bool fwd)
	{
		std::string str;
		if(fwd) str += "Fwd";
		else	str += "Inv";
		str += "Pass"; str += SztToStr(pos);
		return str;
	}

	inline std::string TwTableName()
	{
		return "twiddles";
	}

	inline std::string TwTableLargeName()
	{
		return "twiddle_dee";
	}

	inline std::string TwTableLargeFunc()
	{
		return "TW3step";
	}

	// Twiddle factors table for large N
	// used in 3-step algorithm
    class TwiddleTableLarge
    {
        size_t N; // length
		size_t X, Y;
		size_t tableSize;
		double *wc, *ws; // cosine, sine arrays

	public:
		TwiddleTableLarge(size_t length) : N(length)
		{
			X = size_t(1) << ARBITRARY::TWIDDLE_DEE;
			Y = DivRoundingUp<size_t> (CeilPo2(N), ARBITRARY::TWIDDLE_DEE);
			tableSize = X * Y;

			// Allocate memory for the tables
			wc = new double[tableSize];
			ws = new double[tableSize];
		}

		~TwiddleTableLarge()
		{
			// Free
			delete[] wc;
			delete[] ws;
		}

		template <Precision PR>
		void GenerateTwiddleTable(std::string &twStr)
		{
			const double TWO_PI = -6.283185307179586476925286766559;

			// Generate the table
			size_t nt = 0;
			double phi = TWO_PI / double (N);
			for (size_t iY = 0; iY < Y; ++iY)
			{
				size_t i = size_t(1) << (iY * ARBITRARY::TWIDDLE_DEE);
				for (size_t iX = 0; iX < X; ++iX)
				{
					size_t j = i * iX;

					double c = cos(phi * (double)j);
					double s = sin(phi * (double)j);

					//if (fabs(c) < 1.0E-12)	c = 0.0;
					//if (fabs(s) < 1.0E-12)	s = 0.0;

					wc[nt]   = c;
					ws[nt++] = s;
				}
			}

			std::string sfx = FloatSuffix<PR>();

			// Stringize the table
			std::stringstream ss;
			ss.imbue(std::locale("C"));
			ss.precision(34);
			ss << std::scientific;
			nt = 0;

			ss << "\n __constant ";
			ss << RegBaseType<PR>(2);
			ss << " " << TwTableLargeName();
			ss << "[" << Y << "][" << X << "] = {\n";
			for (size_t iY = 0; iY < Y; ++iY)
			{
				ss << "{ ";
				for (size_t iX = 0; iX < X; ++iX)
				{
					ss << "("; ss << RegBaseType<PR>(2); ss << ")(";
					ss << wc[nt] << sfx << ", ";
					ss << ws[nt++] << sfx << "),\n";
				}
				ss << " },\n";
			}
			ss << "};\n\n";

			// Twiddle calc function
			ss << "__attribute__((always_inline)) ";
			ss << RegBaseType<PR>(2);
			ss << "\n" << TwTableLargeFunc() << "(size_t u)\n{\n";

			ss << "\t" "size_t j = u & " << unsigned(X-1) << ";\n";
			ss << "\t" ; ss << RegBaseType<PR>(2); ss << " result = ";
			ss << TwTableLargeName();
			ss << "[0][j];\n";

			for (size_t iY = 1; iY < Y; ++iY)
			{
				std::string phasor = TwTableLargeName();
				phasor += "[";
				phasor += SztToStr(iY);
				phasor += "][j]";

				stringpair product = ComplexMul((RegBaseType<PR>(2)).c_str(), "result", phasor.c_str());

				ss << "\t" "u >>= " << unsigned (ARBITRARY::TWIDDLE_DEE) << ";\n";
				ss << "\t" "j = u & " << unsigned(X-1) << ";\n";
				ss << "\t" "result = " << product.first << "\n";
				ss << "\t" "\t" << product.second <<";\n";
			}
			ss << "\t" "return result;\n}\n\n";

			twStr += ss.str();
		}
	};

	// FFT butterfly
    template <Precision PR>
    class Butterfly
    {
		size_t radix;		// Base radix
        size_t count;       // Number of basic butterflies, valid values: 1,2,4
		bool fwd;			// FFT direction
		bool cReg;			// registers are complex numbers, .x (real), .y(imag)

		size_t BitReverse (size_t n, size_t N) const
		{
			return (N < 2) ? n : (BitReverse (n >> 1, N >> 1) | ((n & 1) != 0 ? (N >> 1) : 0));
		}

		void GenerateButterflyStr(std::string &bflyStr) const
		{
			std::string regType = cReg ? RegBaseType<PR>(2) : RegBaseType<PR>(count);

			// Function attribute
			bflyStr += "__attribute__((always_inline)) void \n";

			// Function name
			bflyStr += ButterflyName(radix, count, fwd);

			// Function Arguments
			bflyStr += "(";
			for(size_t i=0;;i++)
			{
				if(cReg)
				{
					bflyStr += regType; bflyStr += " *R";
					if(radix & (radix-1))	bflyStr += SztToStr(i);
					else					bflyStr += SztToStr(BitReverse(i,radix));
				}
				else
				{
					bflyStr += regType; bflyStr += " *R"; bflyStr += SztToStr(i); bflyStr += ", ";	// real arguments
					bflyStr += regType; bflyStr += " *I"; bflyStr += SztToStr(i);					// imaginary arguments
				}

				if(i == radix-1)
				{
					bflyStr += ")";
					break;
				}
				else
				{
					bflyStr += ", ";
				}
			}

			bflyStr += "\n{\n\n";


			// Temporary variables
			// Allocate temporary variables if we are not using complex registers (cReg = 0) or if cReg is true, then
			// allocate temporary variables only for non power-of-2 radices
			if (!( (radix == 7 && cReg) || (radix == 11 && cReg) || (radix == 13 && cReg) ))
			{
			if( (radix & (radix-1)) || (!cReg) )
			{
				bflyStr += "\t";
				if(cReg)
					bflyStr += RegBaseType<PR>(1);
				else
					bflyStr += regType;

				for(size_t i=0;;i++)
				{
					bflyStr += " TR"; bflyStr += SztToStr(i); bflyStr += ",";	// real arguments
					bflyStr += " TI"; bflyStr += SztToStr(i);					// imaginary arguments

					if(i == radix-1)
					{
						bflyStr += ";";
						break;
					}
					else
					{
						bflyStr += ",";
					}
				}
			}
			else
			{
				bflyStr += "\t";
				bflyStr += RegBaseType<PR>(2);
				bflyStr += " T;";
				}
			}


			bflyStr += "\n\n\t";

			// Butterfly for different radices
			switch(radix)
			{
			case 2:
				{
					if(cReg)
					{
						bflyStr +=
						"(*R1) = (*R0) - (*R1);\n\t"
						"(*R0) = 2.0f * (*R0) - (*R1);\n\t";
					}
					else
					{
						bflyStr +=
						"TR0 = (*R0) + (*R1);\n\t"
						"TI0 = (*I0) + (*I1);\n\t"
						"TR1 = (*R0) - (*R1);\n\t"
						"TI1 = (*I0) - (*I1);\n\t";
					}

				} break;
			case 3:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R1).x + (*R2).x;\n\t"
							"TR1 = ((*R0).x - C3QA*((*R1).x + (*R2).x)) + C3QB*((*R1).y - (*R2).y);\n\t"
							"TR2 = ((*R0).x - C3QA*((*R1).x + (*R2).x)) - C3QB*((*R1).y - (*R2).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R1).y + (*R2).y;\n\t"
							"TI1 = ((*R0).y - C3QA*((*R1).y + (*R2).y)) - C3QB*((*R1).x - (*R2).x);\n\t"
							"TI2 = ((*R0).y - C3QA*((*R1).y + (*R2).y)) + C3QB*((*R1).x - (*R2).x);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R1 + *R2;\n\t"
							"TR1 = (*R0 - C3QA*(*R1 + *R2)) + C3QB*(*I1 - *I2);\n\t"
							"TR2 = (*R0 - C3QA*(*R1 + *R2)) - C3QB*(*I1 - *I2);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I1 + *I2;\n\t"
							"TI1 = (*I0 - C3QA*(*I1 + *I2)) - C3QB*(*R1 - *R2);\n\t"
							"TI2 = (*I0 - C3QA*(*I1 + *I2)) + C3QB*(*R1 - *R2);\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R1).x + (*R2).x;\n\t"
							"TR1 = ((*R0).x - C3QA*((*R1).x + (*R2).x)) - C3QB*((*R1).y - (*R2).y);\n\t"
							"TR2 = ((*R0).x - C3QA*((*R1).x + (*R2).x)) + C3QB*((*R1).y - (*R2).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R1).y + (*R2).y;\n\t"
							"TI1 = ((*R0).y - C3QA*((*R1).y + (*R2).y)) + C3QB*((*R1).x - (*R2).x);\n\t"
							"TI2 = ((*R0).y - C3QA*((*R1).y + (*R2).y)) - C3QB*((*R1).x - (*R2).x);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R1 + *R2;\n\t"
							"TR1 = (*R0 - C3QA*(*R1 + *R2)) - C3QB*(*I1 - *I2);\n\t"
							"TR2 = (*R0 - C3QA*(*R1 + *R2)) + C3QB*(*I1 - *I2);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I1 + *I2;\n\t"
							"TI1 = (*I0 - C3QA*(*I1 + *I2)) + C3QB*(*R1 - *R2);\n\t"
							"TI2 = (*I0 - C3QA*(*I1 + *I2)) - C3QB*(*R1 - *R2);\n\t";
						}
					}
				} break;
			case 4:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"(*R1) = (*R0) - (*R1);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R1);\n\t"
							"(*R3) = (*R2) - (*R3);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R3);\n\t"
							"\n\t"
							"(*R2) = (*R0) - (*R2);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R2);\n\t"
							"(*R3) = (*R1) + (fvect2)(-(*R3).y, (*R3).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R3);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = (*R0) + (*R2) + (*R1) + (*R3);\n\t"
							"TR1 = (*R0) - (*R2) + (*I1) - (*I3);\n\t"
							"TR2 = (*R0) + (*R2) - (*R1) - (*R3);\n\t"
							"TR3 = (*R0) - (*R2) - (*I1) + (*I3);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*I0) + (*I2) + (*I1) + (*I3);\n\t"
							"TI1 = (*I0) - (*I2) - (*R1) + (*R3);\n\t"
							"TI2 = (*I0) + (*I2) - (*I1) - (*I3);\n\t"
							"TI3 = (*I0) - (*I2) + (*R1) - (*R3);\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"(*R1) = (*R0) - (*R1);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R1);\n\t"
							"(*R3) = (*R2) - (*R3);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R3);\n\t"
							"\n\t"
							"(*R2) = (*R0) - (*R2);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R2);\n\t"
							"(*R3) = (*R1) + (fvect2)((*R3).y, -(*R3).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R3);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = (*R0) + (*R2) + (*R1) + (*R3);\n\t"
							"TR1 = (*R0) - (*R2) - (*I1) + (*I3);\n\t"
							"TR2 = (*R0) + (*R2) - (*R1) - (*R3);\n\t"
							"TR3 = (*R0) - (*R2) + (*I1) - (*I3);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*I0) + (*I2) + (*I1) + (*I3);\n\t"
							"TI1 = (*I0) - (*I2) + (*R1) - (*R3);\n\t"
							"TI2 = (*I0) + (*I2) - (*I1) - (*I3);\n\t"
							"TI3 = (*I0) - (*I2) - (*R1) + (*R3);\n\t";
						}
					}
				} break;
			case 5:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R1).x + (*R2).x + (*R3).x + (*R4).x;\n\t"
							"TR1 = ((*R0).x - C5QC*((*R2).x + (*R3).x)) + C5QB*((*R1).y - (*R4).y) + C5QD*((*R2).y - (*R3).y) + C5QA*(((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));\n\t"
							"TR4 = ((*R0).x - C5QC*((*R2).x + (*R3).x)) - C5QB*((*R1).y - (*R4).y) - C5QD*((*R2).y - (*R3).y) + C5QA*(((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));\n\t"
							"TR2 = ((*R0).x - C5QC*((*R1).x + (*R4).x)) - C5QB*((*R2).y - (*R3).y) + C5QD*((*R1).y - (*R4).y) + C5QA*(((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));\n\t"
							"TR3 = ((*R0).x - C5QC*((*R1).x + (*R4).x)) + C5QB*((*R2).y - (*R3).y) - C5QD*((*R1).y - (*R4).y) + C5QA*(((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R1).y + (*R2).y + (*R3).y + (*R4).y;\n\t"
							"TI1 = ((*R0).y - C5QC*((*R2).y + (*R3).y)) - C5QB*((*R1).x - (*R4).x) - C5QD*((*R2).x - (*R3).x) + C5QA*(((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));\n\t"
							"TI4 = ((*R0).y - C5QC*((*R2).y + (*R3).y)) + C5QB*((*R1).x - (*R4).x) + C5QD*((*R2).x - (*R3).x) + C5QA*(((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));\n\t"
							"TI2 = ((*R0).y - C5QC*((*R1).y + (*R4).y)) + C5QB*((*R2).x - (*R3).x) - C5QD*((*R1).x - (*R4).x) + C5QA*(((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));\n\t"
							"TI3 = ((*R0).y - C5QC*((*R1).y + (*R4).y)) - C5QB*((*R2).x - (*R3).x) + C5QD*((*R1).x - (*R4).x) + C5QA*(((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R1 + *R2 + *R3 + *R4;\n\t"
							"TR1 = (*R0 - C5QC*(*R2 + *R3)) + C5QB*(*I1 - *I4) + C5QD*(*I2 - *I3) + C5QA*((*R1 - *R2) + (*R4 - *R3));\n\t"
							"TR4 = (*R0 - C5QC*(*R2 + *R3)) - C5QB*(*I1 - *I4) - C5QD*(*I2 - *I3) + C5QA*((*R1 - *R2) + (*R4 - *R3));\n\t"
							"TR2 = (*R0 - C5QC*(*R1 + *R4)) - C5QB*(*I2 - *I3) + C5QD*(*I1 - *I4) + C5QA*((*R2 - *R1) + (*R3 - *R4));\n\t"
							"TR3 = (*R0 - C5QC*(*R1 + *R4)) + C5QB*(*I2 - *I3) - C5QD*(*I1 - *I4) + C5QA*((*R2 - *R1) + (*R3 - *R4));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I1 + *I2 + *I3 + *I4;\n\t"
							"TI1 = (*I0 - C5QC*(*I2 + *I3)) - C5QB*(*R1 - *R4) - C5QD*(*R2 - *R3) + C5QA*((*I1 - *I2) + (*I4 - *I3));\n\t"
							"TI4 = (*I0 - C5QC*(*I2 + *I3)) + C5QB*(*R1 - *R4) + C5QD*(*R2 - *R3) + C5QA*((*I1 - *I2) + (*I4 - *I3));\n\t"
							"TI2 = (*I0 - C5QC*(*I1 + *I4)) + C5QB*(*R2 - *R3) - C5QD*(*R1 - *R4) + C5QA*((*I2 - *I1) + (*I3 - *I4));\n\t"
							"TI3 = (*I0 - C5QC*(*I1 + *I4)) - C5QB*(*R2 - *R3) + C5QD*(*R1 - *R4) + C5QA*((*I2 - *I1) + (*I3 - *I4));\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R1).x + (*R2).x + (*R3).x + (*R4).x;\n\t"
							"TR1 = ((*R0).x - C5QC*((*R2).x + (*R3).x)) - C5QB*((*R1).y - (*R4).y) - C5QD*((*R2).y - (*R3).y) + C5QA*(((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));\n\t"
							"TR4 = ((*R0).x - C5QC*((*R2).x + (*R3).x)) + C5QB*((*R1).y - (*R4).y) + C5QD*((*R2).y - (*R3).y) + C5QA*(((*R1).x - (*R2).x) + ((*R4).x - (*R3).x));\n\t"
							"TR2 = ((*R0).x - C5QC*((*R1).x + (*R4).x)) + C5QB*((*R2).y - (*R3).y) - C5QD*((*R1).y - (*R4).y) + C5QA*(((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));\n\t"
							"TR3 = ((*R0).x - C5QC*((*R1).x + (*R4).x)) - C5QB*((*R2).y - (*R3).y) + C5QD*((*R1).y - (*R4).y) + C5QA*(((*R2).x - (*R1).x) + ((*R3).x - (*R4).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R1).y + (*R2).y + (*R3).y + (*R4).y;\n\t"
							"TI1 = ((*R0).y - C5QC*((*R2).y + (*R3).y)) + C5QB*((*R1).x - (*R4).x) + C5QD*((*R2).x - (*R3).x) + C5QA*(((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));\n\t"
							"TI4 = ((*R0).y - C5QC*((*R2).y + (*R3).y)) - C5QB*((*R1).x - (*R4).x) - C5QD*((*R2).x - (*R3).x) + C5QA*(((*R1).y - (*R2).y) + ((*R4).y - (*R3).y));\n\t"
							"TI2 = ((*R0).y - C5QC*((*R1).y + (*R4).y)) - C5QB*((*R2).x - (*R3).x) + C5QD*((*R1).x - (*R4).x) + C5QA*(((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));\n\t"
							"TI3 = ((*R0).y - C5QC*((*R1).y + (*R4).y)) + C5QB*((*R2).x - (*R3).x) - C5QD*((*R1).x - (*R4).x) + C5QA*(((*R2).y - (*R1).y) + ((*R3).y - (*R4).y));\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R1 + *R2 + *R3 + *R4;\n\t"
							"TR1 = (*R0 - C5QC*(*R2 + *R3)) - C5QB*(*I1 - *I4) - C5QD*(*I2 - *I3) + C5QA*((*R1 - *R2) + (*R4 - *R3));\n\t"
							"TR4 = (*R0 - C5QC*(*R2 + *R3)) + C5QB*(*I1 - *I4) + C5QD*(*I2 - *I3) + C5QA*((*R1 - *R2) + (*R4 - *R3));\n\t"
							"TR2 = (*R0 - C5QC*(*R1 + *R4)) + C5QB*(*I2 - *I3) - C5QD*(*I1 - *I4) + C5QA*((*R2 - *R1) + (*R3 - *R4));\n\t"
							"TR3 = (*R0 - C5QC*(*R1 + *R4)) - C5QB*(*I2 - *I3) + C5QD*(*I1 - *I4) + C5QA*((*R2 - *R1) + (*R3 - *R4));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I1 + *I2 + *I3 + *I4;\n\t"
							"TI1 = (*I0 - C5QC*(*I2 + *I3)) + C5QB*(*R1 - *R4) + C5QD*(*R2 - *R3) + C5QA*((*I1 - *I2) + (*I4 - *I3));\n\t"
							"TI4 = (*I0 - C5QC*(*I2 + *I3)) - C5QB*(*R1 - *R4) - C5QD*(*R2 - *R3) + C5QA*((*I1 - *I2) + (*I4 - *I3));\n\t"
							"TI2 = (*I0 - C5QC*(*I1 + *I4)) - C5QB*(*R2 - *R3) + C5QD*(*R1 - *R4) + C5QA*((*I2 - *I1) + (*I3 - *I4));\n\t"
							"TI3 = (*I0 - C5QC*(*I1 + *I4)) + C5QB*(*R2 - *R3) - C5QD*(*R1 - *R4) + C5QA*((*I2 - *I1) + (*I3 - *I4));\n\t";
						}
					}
				} break;
			case 6:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R2).x + (*R4).x;\n\t"
							"TR2 = ((*R0).x - C3QA*((*R2).x + (*R4).x)) + C3QB*((*R2).y - (*R4).y);\n\t"
							"TR4 = ((*R0).x - C3QA*((*R2).x + (*R4).x)) - C3QB*((*R2).y - (*R4).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R2).y + (*R4).y;\n\t"
							"TI2 = ((*R0).y - C3QA*((*R2).y + (*R4).y)) - C3QB*((*R2).x - (*R4).x);\n\t"
							"TI4 = ((*R0).y - C3QA*((*R2).y + (*R4).y)) + C3QB*((*R2).x - (*R4).x);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = (*R1).x + (*R3).x + (*R5).x;\n\t"
							"TR3 = ((*R1).x - C3QA*((*R3).x + (*R5).x)) + C3QB*((*R3).y - (*R5).y);\n\t"
							"TR5 = ((*R1).x - C3QA*((*R3).x + (*R5).x)) - C3QB*((*R3).y - (*R5).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = (*R1).y + (*R3).y + (*R5).y;\n\t"
							"TI3 = ((*R1).y - C3QA*((*R3).y + (*R5).y)) - C3QB*((*R3).x - (*R5).x);\n\t"
							"TI5 = ((*R1).y - C3QA*((*R3).y + (*R5).y)) + C3QB*((*R3).x - (*R5).x);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).x = TR0 + TR1;\n\t"
							"(*R1).x = TR2 + ( C3QA*TR3 + C3QB*TI3);\n\t"
							"(*R2).x = TR4 + (-C3QA*TR5 + C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).y = TI0 + TI1;\n\t"
							"(*R1).y = TI2 + (-C3QB*TR3 + C3QA*TI3);\n\t"
							"(*R2).y = TI4 + (-C3QB*TR5 - C3QA*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3).x = TR0 - TR1;\n\t"
							"(*R4).x = TR2 - ( C3QA*TR3 + C3QB*TI3);\n\t"
							"(*R5).x = TR4 - (-C3QA*TR5 + C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3).y = TI0 - TI1;\n\t"
							"(*R4).y = TI2 - (-C3QB*TR3 + C3QA*TI3);\n\t"
							"(*R5).y = TI4 - (-C3QB*TR5 - C3QA*TI5);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R2 + *R4;\n\t"
							"TR2 = (*R0 - C3QA*(*R2 + *R4)) + C3QB*(*I2 - *I4);\n\t"
							"TR4 = (*R0 - C3QA*(*R2 + *R4)) - C3QB*(*I2 - *I4);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I2 + *I4;\n\t"
							"TI2 = (*I0 - C3QA*(*I2 + *I4)) - C3QB*(*R2 - *R4);\n\t"
							"TI4 = (*I0 - C3QA*(*I2 + *I4)) + C3QB*(*R2 - *R4);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = *R1 + *R3 + *R5;\n\t"
							"TR3 = (*R1 - C3QA*(*R3 + *R5)) + C3QB*(*I3 - *I5);\n\t"
							"TR5 = (*R1 - C3QA*(*R3 + *R5)) - C3QB*(*I3 - *I5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = *I1 + *I3 + *I5;\n\t"
							"TI3 = (*I1 - C3QA*(*I3 + *I5)) - C3QB*(*R3 - *R5);\n\t"
							"TI5 = (*I1 - C3QA*(*I3 + *I5)) + C3QB*(*R3 - *R5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0) = TR0 + TR1;\n\t"
							"(*R1) = TR2 + ( C3QA*TR3 + C3QB*TI3);\n\t"
							"(*R2) = TR4 + (-C3QA*TR5 + C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I0) = TI0 + TI1;\n\t"
							"(*I1) = TI2 + (-C3QB*TR3 + C3QA*TI3);\n\t"
							"(*I2) = TI4 + (-C3QB*TR5 - C3QA*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3) = TR0 - TR1;\n\t"
							"(*R4) = TR2 - ( C3QA*TR3 + C3QB*TI3);\n\t"
							"(*R5) = TR4 - (-C3QA*TR5 + C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I3) = TI0 - TI1;\n\t"
							"(*I4) = TI2 - (-C3QB*TR3 + C3QA*TI3);\n\t"
							"(*I5) = TI4 - (-C3QB*TR5 - C3QA*TI5);\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R2).x + (*R4).x;\n\t"
							"TR2 = ((*R0).x - C3QA*((*R2).x + (*R4).x)) - C3QB*((*R2).y - (*R4).y);\n\t"
							"TR4 = ((*R0).x - C3QA*((*R2).x + (*R4).x)) + C3QB*((*R2).y - (*R4).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R2).y + (*R4).y;\n\t"
							"TI2 = ((*R0).y - C3QA*((*R2).y + (*R4).y)) + C3QB*((*R2).x - (*R4).x);\n\t"
							"TI4 = ((*R0).y - C3QA*((*R2).y + (*R4).y)) - C3QB*((*R2).x - (*R4).x);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = (*R1).x + (*R3).x + (*R5).x;\n\t"
							"TR3 = ((*R1).x - C3QA*((*R3).x + (*R5).x)) - C3QB*((*R3).y - (*R5).y);\n\t"
							"TR5 = ((*R1).x - C3QA*((*R3).x + (*R5).x)) + C3QB*((*R3).y - (*R5).y);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = (*R1).y + (*R3).y + (*R5).y;\n\t"
							"TI3 = ((*R1).y - C3QA*((*R3).y + (*R5).y)) + C3QB*((*R3).x - (*R5).x);\n\t"
							"TI5 = ((*R1).y - C3QA*((*R3).y + (*R5).y)) - C3QB*((*R3).x - (*R5).x);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).x = TR0 + TR1;\n\t"
							"(*R1).x = TR2 + ( C3QA*TR3 - C3QB*TI3);\n\t"
							"(*R2).x = TR4 + (-C3QA*TR5 - C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).y = TI0 + TI1;\n\t"
							"(*R1).y = TI2 + ( C3QB*TR3 + C3QA*TI3);\n\t"
							"(*R2).y = TI4 + ( C3QB*TR5 - C3QA*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3).x = TR0 - TR1;\n\t"
							"(*R4).x = TR2 - ( C3QA*TR3 - C3QB*TI3);\n\t"
							"(*R5).x = TR4 - (-C3QA*TR5 - C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3).y = TI0 - TI1;\n\t"
							"(*R4).y = TI2 - ( C3QB*TR3 + C3QA*TI3);\n\t"
							"(*R5).y = TI4 - ( C3QB*TR5 - C3QA*TI5);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R2 + *R4;\n\t"
							"TR2 = (*R0 - C3QA*(*R2 + *R4)) - C3QB*(*I2 - *I4);\n\t"
							"TR4 = (*R0 - C3QA*(*R2 + *R4)) + C3QB*(*I2 - *I4);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I2 + *I4;\n\t"
							"TI2 = (*I0 - C3QA*(*I2 + *I4)) + C3QB*(*R2 - *R4);\n\t"
							"TI4 = (*I0 - C3QA*(*I2 + *I4)) - C3QB*(*R2 - *R4);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = *R1 + *R3 + *R5;\n\t"
							"TR3 = (*R1 - C3QA*(*R3 + *R5)) - C3QB*(*I3 - *I5);\n\t"
							"TR5 = (*R1 - C3QA*(*R3 + *R5)) + C3QB*(*I3 - *I5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = *I1 + *I3 + *I5;\n\t"
							"TI3 = (*I1 - C3QA*(*I3 + *I5)) + C3QB*(*R3 - *R5);\n\t"
							"TI5 = (*I1 - C3QA*(*I3 + *I5)) - C3QB*(*R3 - *R5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0) = TR0 + TR1;\n\t"
							"(*R1) = TR2 + ( C3QA*TR3 - C3QB*TI3);\n\t"
							"(*R2) = TR4 + (-C3QA*TR5 - C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I0) = TI0 + TI1;\n\t"
							"(*I1) = TI2 + ( C3QB*TR3 + C3QA*TI3);\n\t"
							"(*I2) = TI4 + ( C3QB*TR5 - C3QA*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R3) = TR0 - TR1;\n\t"
							"(*R4) = TR2 - ( C3QA*TR3 - C3QB*TI3);\n\t"
							"(*R5) = TR4 - (-C3QA*TR5 - C3QB*TI5);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I3) = TI0 - TI1;\n\t"
							"(*I4) = TI2 - ( C3QB*TR3 + C3QA*TI3);\n\t"
							"(*I5) = TI4 - ( C3QB*TR5 - C3QA*TI5);\n\t";
						}
					}
				} break;
			case 7:
				{
					static const char *C7SFR = "\
					/*FFT7 Forward Real */ \n\
					\n\
						pr0 = *R1 + *R6; \n\
						pi0 = *I1 + *I6; \n\
						pr1 = *R1 - *R6; \n\
						pi1 = *I1 - *I6; \n\
						pr2 = *R2 + *R5; \n\
						pi2 = *I2 + *I5; \n\
						pr3 = *R2 - *R5; \n\
						pi3 = *I2 - *I5; \n\
						pr4 = *R4 + *R3; \n\
						pi4 = *I4 + *I3; \n\
						pr5 = *R4 - *R3; \n\
						pi5 = *I4 - *I3; \n\
					\n\
						pr6 = pr2 + pr0; \n\
						pi6 = pi2 + pi0; \n\
						qr4 = pr2 - pr0; \n\
						qi4 = pi2 - pi0; \n\
						qr2 = pr0 - pr4; \n\
						qi2 = pi0 - pi4; \n\
						qr3 = pr4 - pr2; \n\
						qi3 = pi4 - pi2; \n\
						pr7 = pr5 + pr3; \n\
						pi7 = pi5 + pi3; \n\
						qr7 = pr5 - pr3; \n\
						qi7 = pi5 - pi3; \n\
						qr6 = pr1 - pr5; \n\
						qi6 = pi1 - pi5; \n\
						qr8 = pr3 - pr1; \n\
						qi8 = pi3 - pi1; \n\
						qr1 = pr6 + pr4; \n\
						qi1 = pi6 + pi4; \n\
						qr5 = pr7 + pr1; \n\
						qi5 = pi7 + pi1; \n\
						qr0 = *R0 + qr1; \n\
						qi0 = *I0 + qi1; \n\
					\n\
						qr1 *= C7Q1; \n\
						qi1 *= C7Q1; \n\
						qr2 *= C7Q2; \n\
						qi2 *= C7Q2; \n\
						qr3 *= C7Q3; \n\
						qi3 *= C7Q3; \n\
						qr4 *= C7Q4; \n\
						qi4 *= C7Q4; \n\
					\n\
						qr5 *= (C7Q5); \n\
						qi5 *= (C7Q5); \n\
						qr6 *= (C7Q6); \n\
						qi6 *= (C7Q6); \n\
						qr7 *= (C7Q7); \n\
						qi7 *= (C7Q7); \n\
						qr8 *= (C7Q8); \n\
						qi8 *= (C7Q8); \n\
					\n\
						pr0 =  qr0 + qr1; \n\
						pi0 =  qi0 + qi1; \n\
						pr1 =  qr2 + qr3; \n\
						pi1 =  qi2 + qi3; \n\
						pr2 =  qr4 - qr3; \n\
						pi2 =  qi4 - qi3; \n\
						pr3 = -qr2 - qr4; \n\
						pi3 = -qi2 - qi4; \n\
						pr4 =  qr6 + qr7; \n\
						pi4 =  qi6 + qi7; \n\
						pr5 =  qr8 - qr7; \n\
						pi5 =  qi8 - qi7; \n\
						pr6 = -qr8 - qr6; \n\
						pi6 = -qi8 - qi6; \n\
						pr7 =  pr0 + pr1; \n\
						pi7 =  pi0 + pi1; \n\
						pr8 =  pr0 + pr2; \n\
						pi8 =  pi0 + pi2; \n\
						pr9 =  pr0 + pr3; \n\
						pi9 =  pi0 + pi3; \n\
						qr6 =  pr4 + qr5; \n\
						qi6 =  pi4 + qi5; \n\
						qr7 =  pr5 + qr5; \n\
						qi7 =  pi5 + qi5; \n\
						qr8 =  pr6 + qr5; \n\
						qi8 =  pi6 + qi5; \n\
					\n\
						TR0 = qr0; TI0 = qi0; \n\
						TR1 = pr7 + qi6; \n\
						TI1 = pi7 - qr6; \n\
						TR2 = pr9 + qi8; \n\
						TI2 = pi9 - qr8; \n\
						TR3 = pr8 - qi7; \n\
						TI3 = pi8 + qr7; \n\
						TR4 = pr8 + qi7; \n\
						TI4 = pi8 - qr7; \n\
						TR5 = pr9 - qi8; \n\
						TI5 = pi9 + qr8; \n\
						TR6 = pr7 - qi6; \n\
						TI6 = pi7 + qr6; \n\
					";

					static const char *C7SBR = "\
					/*FFT7 Backward Real */ \n\
					\n\
						pr0 = *R1 + *R6; \n\
						pi0 = *I1 + *I6; \n\
						pr1 = *R1 - *R6; \n\
						pi1 = *I1 - *I6; \n\
						pr2 = *R2 + *R5; \n\
						pi2 = *I2 + *I5; \n\
						pr3 = *R2 - *R5; \n\
						pi3 = *I2 - *I5; \n\
						pr4 = *R4 + *R3; \n\
						pi4 = *I4 + *I3; \n\
						pr5 = *R4 - *R3; \n\
						pi5 = *I4 - *I3; \n\
					\n\
						pr6 = pr2 + pr0; \n\
						pi6 = pi2 + pi0; \n\
						qr4 = pr2 - pr0; \n\
						qi4 = pi2 - pi0; \n\
						qr2 = pr0 - pr4; \n\
						qi2 = pi0 - pi4; \n\
						qr3 = pr4 - pr2; \n\
						qi3 = pi4 - pi2; \n\
						pr7 = pr5 + pr3; \n\
						pi7 = pi5 + pi3; \n\
						qr7 = pr5 - pr3; \n\
						qi7 = pi5 - pi3; \n\
						qr6 = pr1 - pr5; \n\
						qi6 = pi1 - pi5; \n\
						qr8 = pr3 - pr1; \n\
						qi8 = pi3 - pi1; \n\
						qr1 = pr6 + pr4; \n\
						qi1 = pi6 + pi4; \n\
						qr5 = pr7 + pr1; \n\
						qi5 = pi7 + pi1; \n\
						qr0 = *R0 + qr1; \n\
						qi0 = *I0 + qi1; \n\
					\n\
						qr1 *= C7Q1; \n\
						qi1 *= C7Q1; \n\
						qr2 *= C7Q2; \n\
						qi2 *= C7Q2; \n\
						qr3 *= C7Q3; \n\
						qi3 *= C7Q3; \n\
						qr4 *= C7Q4; \n\
						qi4 *= C7Q4; \n\
					\n\
						qr5 *= -(C7Q5); \n\
						qi5 *= -(C7Q5); \n\
						qr6 *= -(C7Q6); \n\
						qi6 *= -(C7Q6); \n\
						qr7 *= -(C7Q7); \n\
						qi7 *= -(C7Q7); \n\
						qr8 *= -(C7Q8); \n\
						qi8 *= -(C7Q8); \n\
					\n\
						pr0 =  qr0 + qr1; \n\
						pi0 =  qi0 + qi1; \n\
						pr1 =  qr2 + qr3; \n\
						pi1 =  qi2 + qi3; \n\
						pr2 =  qr4 - qr3; \n\
						pi2 =  qi4 - qi3; \n\
						pr3 = -qr2 - qr4; \n\
						pi3 = -qi2 - qi4; \n\
						pr4 =  qr6 + qr7; \n\
						pi4 =  qi6 + qi7; \n\
						pr5 =  qr8 - qr7; \n\
						pi5 =  qi8 - qi7; \n\
						pr6 = -qr8 - qr6; \n\
						pi6 = -qi8 - qi6; \n\
						pr7 =  pr0 + pr1; \n\
						pi7 =  pi0 + pi1; \n\
						pr8 =  pr0 + pr2; \n\
						pi8 =  pi0 + pi2; \n\
						pr9 =  pr0 + pr3; \n\
						pi9 =  pi0 + pi3; \n\
						qr6 =  pr4 + qr5; \n\
						qi6 =  pi4 + qi5; \n\
						qr7 =  pr5 + qr5; \n\
						qi7 =  pi5 + qi5; \n\
						qr8 =  pr6 + qr5; \n\
						qi8 =  pi6 + qi5; \n\
					\n\
						TR0 = qr0; TI0 = qi0; \n\
						TR1 = pr7 + qi6; \n\
						TI1 = pi7 - qr6; \n\
						TR2 = pr9 + qi8; \n\
						TI2 = pi9 - qr8; \n\
						TR3 = pr8 - qi7; \n\
						TI3 = pi8 + qr7; \n\
						TR4 = pr8 + qi7; \n\
						TI4 = pi8 - qr7; \n\
						TR5 = pr9 - qi8; \n\
						TI5 = pi9 + qr8; \n\
						TR6 = pr7 - qi6; \n\
						TI6 = pi7 + qr6; \n\
					";

					static const char *C7SFC = "\
					/*FFT7 Forward Complex */ \n\
					\n\
						p0 = *R1 + *R6; \n\
						p1 = *R1 - *R6; \n\
						p2 = *R2 + *R5; \n\
						p3 = *R2 - *R5; \n\
						p4 = *R4 + *R3; \n\
						p5 = *R4 - *R3; \n\
					\n\
						p6 = p2 + p0; \n\
						q4 = p2 - p0; \n\
						q2 = p0 - p4; \n\
						q3 = p4 - p2; \n\
						p7 = p5 + p3; \n\
						q7 = p5 - p3; \n\
						q6 = p1 - p5; \n\
						q8 = p3 - p1; \n\
						q1 = p6 + p4; \n\
						q5 = p7 + p1; \n\
						q0 = *R0 + q1; \n\
					\n\
						q1 *= C7Q1; \n\
						q2 *= C7Q2; \n\
						q3 *= C7Q3; \n\
						q4 *= C7Q4; \n\
					\n\
						q5 *= (C7Q5); \n\
						q6 *= (C7Q6); \n\
						q7 *= (C7Q7); \n\
						q8 *= (C7Q8); \n\
					\n\
						p0 = q0 + q1; \n\
						p1 = q2 + q3; \n\
						p2 = q4 - q3; \n\
						p3 = -q2 - q4; \n\
						p4 = q6 + q7; \n\
						p5 = q8 - q7; \n\
						p6 = -q8 - q6; \n\
						p7 = p0 + p1; \n\
						p8 = p0 + p2; \n\
						p9 = p0 + p3; \n\
						q6 = p4 + q5; \n\
						q7 = p5 + q5; \n\
						q8 = p6 + q5; \n\
					\n\
						*R0 = q0; \n\
						(*R1).x = p7.x + q6.y; \n\
						(*R1).y = p7.y - q6.x; \n\
						(*R2).x = p9.x + q8.y; \n\
						(*R2).y = p9.y - q8.x; \n\
						(*R3).x = p8.x - q7.y; \n\
						(*R3).y = p8.y + q7.x; \n\
						(*R4).x = p8.x + q7.y; \n\
						(*R4).y = p8.y - q7.x; \n\
						(*R5).x = p9.x - q8.y; \n\
						(*R5).y = p9.y + q8.x; \n\
						(*R6).x = p7.x - q6.y; \n\
						(*R6).y = p7.y + q6.x; \n\
					";

					static const char *C7SBC = "\
					/*FFT7 Backward Complex */ \n\
					\n\
						p0 = *R1 + *R6; \n\
						p1 = *R1 - *R6; \n\
						p2 = *R2 + *R5; \n\
						p3 = *R2 - *R5; \n\
						p4 = *R4 + *R3; \n\
						p5 = *R4 - *R3; \n\
					\n\
						p6 = p2 + p0; \n\
						q4 = p2 - p0; \n\
						q2 = p0 - p4; \n\
						q3 = p4 - p2; \n\
						p7 = p5 + p3; \n\
						q7 = p5 - p3; \n\
						q6 = p1 - p5; \n\
						q8 = p3 - p1; \n\
						q1 = p6 + p4; \n\
						q5 = p7 + p1; \n\
						q0 = *R0 + q1; \n\
					\n\
						q1 *= C7Q1; \n\
						q2 *= C7Q2; \n\
						q3 *= C7Q3; \n\
						q4 *= C7Q4; \n\
					\n\
						q5 *= -(C7Q5); \n\
						q6 *= -(C7Q6); \n\
						q7 *= -(C7Q7); \n\
						q8 *= -(C7Q8); \n\
					\n\
						p0 = q0 + q1; \n\
						p1 = q2 + q3; \n\
						p2 = q4 - q3; \n\
						p3 = -q2 - q4; \n\
						p4 = q6 + q7; \n\
						p5 = q8 - q7; \n\
						p6 = -q8 - q6; \n\
						p7 = p0 + p1; \n\
						p8 = p0 + p2; \n\
						p9 = p0 + p3; \n\
						q6 = p4 + q5; \n\
						q7 = p5 + q5; \n\
						q8 = p6 + q5; \n\
					\n\
						*R0 = q0; \n\
						(*R1).x = p7.x + q6.y; \n\
						(*R1).y = p7.y - q6.x; \n\
						(*R2).x = p9.x + q8.y; \n\
						(*R2).y = p9.y - q8.x; \n\
						(*R3).x = p8.x - q7.y; \n\
						(*R3).y = p8.y + q7.x; \n\
						(*R4).x = p8.x + q7.y; \n\
						(*R4).y = p8.y - q7.x; \n\
						(*R5).x = p9.x - q8.y; \n\
						(*R5).y = p9.y + q8.x; \n\
						(*R6).x = p7.x - q6.y; \n\
						(*R6).y = p7.y + q6.x; \n\
					";



					if (!cReg) {
						for (size_t i = 0; i < 10; i++)
							bflyStr += regType + " pr" + SztToStr(i) + ", pi" + SztToStr(i) + ";\n\t";
						for (size_t i = 0; i < 9; i++)
							bflyStr += regType + " qr" + SztToStr(i) + ", qi" + SztToStr(i) + ";\n\t";

						if (fwd)
							bflyStr += C7SFR;
						else
							bflyStr += C7SBR;
					} else {
						for (size_t i = 0; i < 10; i++)
							bflyStr += regType + " p" + SztToStr(i) + ";\n\t";
						for (size_t i = 0; i < 9; i++)
							bflyStr += regType + " q" + SztToStr(i) + ";\n\t";
						if (fwd)
							bflyStr += C7SFC;
						else
							bflyStr += C7SBC;
					}
				}
				break;

			case 8:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"(*R1) = (*R0) - (*R1);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R1);\n\t"
							"(*R3) = (*R2) - (*R3);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R3);\n\t"
							"(*R5) = (*R4) - (*R5);\n\t"
							"(*R4) = 2.0f * (*R4) - (*R5);\n\t"
							"(*R7) = (*R6) - (*R7);\n\t"
							"(*R6) = 2.0f * (*R6) - (*R7);\n\t"
							"\n\t"
							"(*R2) = (*R0) - (*R2);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R2);\n\t"
							"(*R3) = (*R1) + (fvect2)(-(*R3).y, (*R3).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R3);\n\t"
							"(*R6) = (*R4) - (*R6);\n\t"
							"(*R4) = 2.0f * (*R4) - (*R6);\n\t"
							"(*R7) = (*R5) + (fvect2)(-(*R7).y, (*R7).x);\n\t"
							"(*R5) = 2.0f * (*R5) - (*R7);\n\t"
							"\n\t"
							"(*R4) = (*R0) - (*R4);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R4);\n\t"
							"(*R5) = ((*R1) - C8Q * (*R5)) - C8Q * (fvect2)((*R5).y, -(*R5).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R5);\n\t"
							"(*R6) = (*R2) + (fvect2)(-(*R6).y, (*R6).x);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R6);\n\t"
							"(*R7) = ((*R3) + C8Q * (*R7)) - C8Q * (fvect2)((*R7).y, -(*R7).x);\n\t"
							"(*R3) = 2.0f * (*R3) - (*R7);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = (*R0) + (*R4) + (*R2) + (*R6) +     (*R1)             +     (*R3)             +     (*R5)             +     (*R7)            ;\n\t"
							"TR1 = (*R0) - (*R4) + (*I2) - (*I6) + C8Q*(*R1) + C8Q*(*I1) - C8Q*(*R3) + C8Q*(*I3) - C8Q*(*R5) - C8Q*(*I5) + C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TR2 = (*R0) + (*R4) - (*R2) - (*R6)             +     (*I1)             -     (*I3)             +     (*I5)             -     (*I7);\n\t"
							"TR3 = (*R0) - (*R4) - (*I2) + (*I6) - C8Q*(*R1) + C8Q*(*I1) + C8Q*(*R3) + C8Q*(*I3) + C8Q*(*R5) - C8Q*(*I5) - C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TR4 = (*R0) + (*R4) + (*R2) + (*R6) -     (*R1)             -     (*R3)             -     (*R5)             -     (*R7)            ;\n\t"
							"TR5 = (*R0) - (*R4) + (*I2) - (*I6) - C8Q*(*R1) - C8Q*(*I1) + C8Q*(*R3) - C8Q*(*I3) + C8Q*(*R5) + C8Q*(*I5) - C8Q*(*R7) + C8Q*(*I7);\n\t"
							"TR6 = (*R0) + (*R4) - (*R2) - (*R6)             -    (*I1)              +     (*I3)             -     (*I5)             +     (*I7);\n\t"
							"TR7 = (*R0) - (*R4) - (*I2) + (*I6) + C8Q*(*R1) - C8Q*(*I1) - C8Q*(*R3) - C8Q*(*I3) - C8Q*(*R5) + C8Q*(*I5) + C8Q*(*R7) + C8Q*(*I7);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*I0) + (*I4) + (*I2) + (*I6)             +     (*I1)             +     (*I3)             +     (*I5)             +     (*I7);\n\t"
							"TI1 = (*I0) - (*I4) - (*R2) + (*R6) - C8Q*(*R1) + C8Q*(*I1) - C8Q*(*R3) - C8Q*(*I3) + C8Q*(*R5) - C8Q*(*I5) + C8Q*(*R7) + C8Q*(*I7);\n\t"
							"TI2 = (*I0) + (*I4) - (*I2) - (*I6) -     (*R1)             +     (*R3)             -     (*R5)             +     (*R7)            ;\n\t"
							"TI3 = (*I0) - (*I4) + (*R2) - (*R6) - C8Q*(*R1) - C8Q*(*I1) - C8Q*(*R3) + C8Q*(*I3) + C8Q*(*R5) + C8Q*(*I5) + C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TI4 = (*I0) + (*I4) + (*I2) + (*I6)             -    (*I1)              -     (*I3)             -     (*I5)             -     (*I7);\n\t"
							"TI5 = (*I0) - (*I4) - (*R2) + (*R6) + C8Q*(*R1) - C8Q*(*I1) + C8Q*(*R3) + C8Q*(*I3) - C8Q*(*R5) + C8Q*(*I5) - C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TI6 = (*I0) + (*I4) - (*I2) - (*I6) +     (*R1)             -     (*R3)             +     (*R5)             -     (*R7)            ;\n\t"
							"TI7 = (*I0) - (*I4) + (*R2) - (*R6) + C8Q*(*R1) + C8Q*(*I1) + C8Q*(*R3) - C8Q*(*I3) - C8Q*(*R5) - C8Q*(*I5) - C8Q*(*R7) + C8Q*(*I7);\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"(*R1) = (*R0) - (*R1);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R1);\n\t"
							"(*R3) = (*R2) - (*R3);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R3);\n\t"
							"(*R5) = (*R4) - (*R5);\n\t"
							"(*R4) = 2.0f * (*R4) - (*R5);\n\t"
							"(*R7) = (*R6) - (*R7);\n\t"
							"(*R6) = 2.0f * (*R6) - (*R7);\n\t"
							"\n\t"
							"(*R2) = (*R0) - (*R2);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R2);\n\t"
							"(*R3) = (*R1) + (fvect2)((*R3).y, -(*R3).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R3);\n\t"
							"(*R6) = (*R4) - (*R6);\n\t"
							"(*R4) = 2.0f * (*R4) - (*R6);\n\t"
							"(*R7) = (*R5) + (fvect2)((*R7).y, -(*R7).x);\n\t"
							"(*R5) = 2.0f * (*R5) - (*R7);\n\t"
							"\n\t"
							"(*R4) = (*R0) - (*R4);\n\t"
							"(*R0) = 2.0f * (*R0) - (*R4);\n\t"
							"(*R5) = ((*R1) - C8Q * (*R5)) + C8Q * (fvect2)((*R5).y, -(*R5).x);\n\t"
							"(*R1) = 2.0f * (*R1) - (*R5);\n\t"
							"(*R6) = (*R2) + (fvect2)((*R6).y, -(*R6).x);\n\t"
							"(*R2) = 2.0f * (*R2) - (*R6);\n\t"
							"(*R7) = ((*R3) + C8Q * (*R7)) + C8Q * (fvect2)((*R7).y, -(*R7).x);\n\t"
							"(*R3) = 2.0f * (*R3) - (*R7);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = (*R0) + (*R4) + (*R2) + (*R6) +     (*R1)             +     (*R3)             +     (*R5)             +     (*R7)            ;\n\t"
							"TR1 = (*R0) - (*R4) - (*I2) + (*I6) + C8Q*(*R1) - C8Q*(*I1) - C8Q*(*R3) - C8Q*(*I3) - C8Q*(*R5) + C8Q*(*I5) + C8Q*(*R7) + C8Q*(*I7);\n\t"
							"TR2 = (*R0) + (*R4) - (*R2) - (*R6)             -     (*I1)             +     (*I3)             -     (*I5)             +     (*I7);\n\t"
							"TR3 = (*R0) - (*R4) + (*I2) - (*I6) - C8Q*(*R1) - C8Q*(*I1) + C8Q*(*R3) - C8Q*(*I3) + C8Q*(*R5) + C8Q*(*I5) - C8Q*(*R7) + C8Q*(*I7);\n\t"
							"TR4 = (*R0) + (*R4) + (*R2) + (*R6) -     (*R1)             -    (*R3)              -     (*R5)             -     (*R7)            ;\n\t"
							"TR5 = (*R0) - (*R4) - (*I2) + (*I6) - C8Q*(*R1) + C8Q*(*I1) + C8Q*(*R3) + C8Q*(*I3) + C8Q*(*R5) - C8Q*(*I5) - C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TR6 = (*R0) + (*R4) - (*R2) - (*R6)             +     (*I1)             -     (*I3)             +     (*I5)             -     (*I7);\n\t"
							"TR7 = (*R0) - (*R4) + (*I2) - (*I6) + C8Q*(*R1) + C8Q*(*I1) - C8Q*(*R3) + C8Q*(*I3) - C8Q*(*R5) - C8Q*(*I5) + C8Q*(*R7) - C8Q*(*I7);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*I0) + (*I4) + (*I2) + (*I6)             +     (*I1)             +    (*I3)              +     (*I5)             +     (*I7);\n\t"
							"TI1 = (*I0) - (*I4) + (*R2) - (*R6) + C8Q*(*R1) + C8Q*(*I1) + C8Q*(*R3) - C8Q*(*I3) - C8Q*(*R5) - C8Q*(*I5) - C8Q*(*R7) + C8Q*(*I7);\n\t"
							"TI2 = (*I0) + (*I4) - (*I2) - (*I6) +     (*R1)             -     (*R3)             +     (*R5)             -     (*R7)            ;\n\t"
							"TI3 = (*I0) - (*I4) - (*R2) + (*R6) + C8Q*(*R1) - C8Q*(*I1) + C8Q*(*R3) + C8Q*(*I3) - C8Q*(*R5) + C8Q*(*I5) - C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TI4 = (*I0) + (*I4) + (*I2) + (*I6)             -     (*I1)             -     (*I3)             -     (*I5)             -     (*I7);\n\t"
							"TI5 = (*I0) - (*I4) + (*R2) - (*R6) - C8Q*(*R1) - C8Q*(*I1) - C8Q*(*R3) + C8Q*(*I3) + C8Q*(*R5) + C8Q*(*I5) + C8Q*(*R7) - C8Q*(*I7);\n\t"
							"TI6 = (*I0) + (*I4) - (*I2) - (*I6) -     (*R1)             +     (*R3)             -     (*R5)             +     (*R7)            ;\n\t"
							"TI7 = (*I0) - (*I4) - (*R2) + (*R6) - C8Q*(*R1) + C8Q*(*I1) - C8Q*(*R3) - C8Q*(*I3) + C8Q*(*R5) - C8Q*(*I5) + C8Q*(*R7) + C8Q*(*I7);\n\t";
						}
					}
				} break;
			case 10:
				{
					if(fwd)
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R2).x + (*R4).x + (*R6).x + (*R8).x;\n\t"
							"TR2 = ((*R0).x - C5QC*((*R4).x + (*R6).x)) + C5QB*((*R2).y - (*R8).y) + C5QD*((*R4).y - (*R6).y) + C5QA*(((*R2).x - (*R4).x) + ((*R8).x - (*R6).x));\n\t"
							"TR8 = ((*R0).x - C5QC*((*R4).x + (*R6).x)) - C5QB*((*R2).y - (*R8).y) - C5QD*((*R4).y - (*R6).y) + C5QA*(((*R2).x - (*R4).x) + ((*R8).x - (*R6).x));\n\t"
							"TR4 = ((*R0).x - C5QC*((*R2).x + (*R8).x)) - C5QB*((*R4).y - (*R6).y) + C5QD*((*R2).y - (*R8).y) + C5QA*(((*R4).x - (*R2).x) + ((*R6).x - (*R8).x));\n\t"
							"TR6 = ((*R0).x - C5QC*((*R2).x + (*R8).x)) + C5QB*((*R4).y - (*R6).y) - C5QD*((*R2).y - (*R8).y) + C5QA*(((*R4).x - (*R2).x) + ((*R6).x - (*R8).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R2).y + (*R4).y + (*R6).y + (*R8).y;\n\t"
							"TI2 = ((*R0).y - C5QC*((*R4).y + (*R6).y)) - C5QB*((*R2).x - (*R8).x) - C5QD*((*R4).x - (*R6).x) + C5QA*(((*R2).y - (*R4).y) + ((*R8).y - (*R6).y));\n\t"
							"TI8 = ((*R0).y - C5QC*((*R4).y + (*R6).y)) + C5QB*((*R2).x - (*R8).x) + C5QD*((*R4).x - (*R6).x) + C5QA*(((*R2).y - (*R4).y) + ((*R8).y - (*R6).y));\n\t"
							"TI4 = ((*R0).y - C5QC*((*R2).y + (*R8).y)) + C5QB*((*R4).x - (*R6).x) - C5QD*((*R2).x - (*R8).x) + C5QA*(((*R4).y - (*R2).y) + ((*R6).y - (*R8).y));\n\t"
							"TI6 = ((*R0).y - C5QC*((*R2).y + (*R8).y)) - C5QB*((*R4).x - (*R6).x) + C5QD*((*R2).x - (*R8).x) + C5QA*(((*R4).y - (*R2).y) + ((*R6).y - (*R8).y));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = (*R1).x + (*R3).x + (*R5).x + (*R7).x + (*R9).x;\n\t"
							"TR3 = ((*R1).x - C5QC*((*R5).x + (*R7).x)) + C5QB*((*R3).y - (*R9).y) + C5QD*((*R5).y - (*R7).y) + C5QA*(((*R3).x - (*R5).x) + ((*R9).x - (*R7).x));\n\t"
							"TR9 = ((*R1).x - C5QC*((*R5).x + (*R7).x)) - C5QB*((*R3).y - (*R9).y) - C5QD*((*R5).y - (*R7).y) + C5QA*(((*R3).x - (*R5).x) + ((*R9).x - (*R7).x));\n\t"
							"TR5 = ((*R1).x - C5QC*((*R3).x + (*R9).x)) - C5QB*((*R5).y - (*R7).y) + C5QD*((*R3).y - (*R9).y) + C5QA*(((*R5).x - (*R3).x) + ((*R7).x - (*R9).x));\n\t"
							"TR7 = ((*R1).x - C5QC*((*R3).x + (*R9).x)) + C5QB*((*R5).y - (*R7).y) - C5QD*((*R3).y - (*R9).y) + C5QA*(((*R5).x - (*R3).x) + ((*R7).x - (*R9).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = (*R1).y + (*R3).y + (*R5).y + (*R7).y + (*R9).y;\n\t"
							"TI3 = ((*R1).y - C5QC*((*R5).y + (*R7).y)) - C5QB*((*R3).x - (*R9).x) - C5QD*((*R5).x - (*R7).x) + C5QA*(((*R3).y - (*R5).y) + ((*R9).y - (*R7).y));\n\t"
							"TI9 = ((*R1).y - C5QC*((*R5).y + (*R7).y)) + C5QB*((*R3).x - (*R9).x) + C5QD*((*R5).x - (*R7).x) + C5QA*(((*R3).y - (*R5).y) + ((*R9).y - (*R7).y));\n\t"
							"TI5 = ((*R1).y - C5QC*((*R3).y + (*R9).y)) + C5QB*((*R5).x - (*R7).x) - C5QD*((*R3).x - (*R9).x) + C5QA*(((*R5).y - (*R3).y) + ((*R7).y - (*R9).y));\n\t"
							"TI7 = ((*R1).y - C5QC*((*R3).y + (*R9).y)) - C5QB*((*R5).x - (*R7).x) + C5QD*((*R3).x - (*R9).x) + C5QA*(((*R5).y - (*R3).y) + ((*R7).y - (*R9).y));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).x = TR0 + TR1;\n\t"
							"(*R1).x = TR2 + ( C5QE*TR3 + C5QD*TI3);\n\t"
							"(*R2).x = TR4 + ( C5QA*TR5 + C5QB*TI5);\n\t"
							"(*R3).x = TR6 + (-C5QA*TR7 + C5QB*TI7);\n\t"
							"(*R4).x = TR8 + (-C5QE*TR9 + C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).y = TI0 + TI1;\n\t"
							"(*R1).y = TI2 + (-C5QD*TR3 + C5QE*TI3);\n\t"
							"(*R2).y = TI4 + (-C5QB*TR5 + C5QA*TI5);\n\t"
							"(*R3).y = TI6 + (-C5QB*TR7 - C5QA*TI7);\n\t"
							"(*R4).y = TI8 + (-C5QD*TR9 - C5QE*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5).x = TR0 - TR1;\n\t"
							"(*R6).x = TR2 - ( C5QE*TR3 + C5QD*TI3);\n\t"
							"(*R7).x = TR4 - ( C5QA*TR5 + C5QB*TI5);\n\t"
							"(*R8).x = TR6 - (-C5QA*TR7 + C5QB*TI7);\n\t"
							"(*R9).x = TR8 - (-C5QE*TR9 + C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5).y = TI0 - TI1;\n\t"
							"(*R6).y = TI2 - (-C5QD*TR3 + C5QE*TI3);\n\t"
							"(*R7).y = TI4 - (-C5QB*TR5 + C5QA*TI5);\n\t"
							"(*R8).y = TI6 - (-C5QB*TR7 - C5QA*TI7);\n\t"
							"(*R9).y = TI8 - (-C5QD*TR9 - C5QE*TI9);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R2 + *R4 + *R6 + *R8;\n\t"
							"TR2 = (*R0 - C5QC*(*R4 + *R6)) + C5QB*(*I2 - *I8) + C5QD*(*I4 - *I6) + C5QA*((*R2 - *R4) + (*R8 - *R6));\n\t"
							"TR8 = (*R0 - C5QC*(*R4 + *R6)) - C5QB*(*I2 - *I8) - C5QD*(*I4 - *I6) + C5QA*((*R2 - *R4) + (*R8 - *R6));\n\t"
							"TR4 = (*R0 - C5QC*(*R2 + *R8)) - C5QB*(*I4 - *I6) + C5QD*(*I2 - *I8) + C5QA*((*R4 - *R2) + (*R6 - *R8));\n\t"
							"TR6 = (*R0 - C5QC*(*R2 + *R8)) + C5QB*(*I4 - *I6) - C5QD*(*I2 - *I8) + C5QA*((*R4 - *R2) + (*R6 - *R8));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I2 + *I4 + *I6 + *I8;\n\t"
							"TI2 = (*I0 - C5QC*(*I4 + *I6)) - C5QB*(*R2 - *R8) - C5QD*(*R4 - *R6) + C5QA*((*I2 - *I4) + (*I8 - *I6));\n\t"
							"TI8 = (*I0 - C5QC*(*I4 + *I6)) + C5QB*(*R2 - *R8) + C5QD*(*R4 - *R6) + C5QA*((*I2 - *I4) + (*I8 - *I6));\n\t"
							"TI4 = (*I0 - C5QC*(*I2 + *I8)) + C5QB*(*R4 - *R6) - C5QD*(*R2 - *R8) + C5QA*((*I4 - *I2) + (*I6 - *I8));\n\t"
							"TI6 = (*I0 - C5QC*(*I2 + *I8)) - C5QB*(*R4 - *R6) + C5QD*(*R2 - *R8) + C5QA*((*I4 - *I2) + (*I6 - *I8));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = *R1 + *R3 + *R5 + *R7 + *R9;\n\t"
							"TR3 = (*R1 - C5QC*(*R5 + *R7)) + C5QB*(*I3 - *I9) + C5QD*(*I5 - *I7) + C5QA*((*R3 - *R5) + (*R9 - *R7));\n\t"
							"TR9 = (*R1 - C5QC*(*R5 + *R7)) - C5QB*(*I3 - *I9) - C5QD*(*I5 - *I7) + C5QA*((*R3 - *R5) + (*R9 - *R7));\n\t"
							"TR5 = (*R1 - C5QC*(*R3 + *R9)) - C5QB*(*I5 - *I7) + C5QD*(*I3 - *I9) + C5QA*((*R5 - *R3) + (*R7 - *R9));\n\t"
							"TR7 = (*R1 - C5QC*(*R3 + *R9)) + C5QB*(*I5 - *I7) - C5QD*(*I3 - *I9) + C5QA*((*R5 - *R3) + (*R7 - *R9));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = *I1 + *I3 + *I5 + *I7 + *I9;\n\t"
							"TI3 = (*I1 - C5QC*(*I5 + *I7)) - C5QB*(*R3 - *R9) - C5QD*(*R5 - *R7) + C5QA*((*I3 - *I5) + (*I9 - *I7));\n\t"
							"TI9 = (*I1 - C5QC*(*I5 + *I7)) + C5QB*(*R3 - *R9) + C5QD*(*R5 - *R7) + C5QA*((*I3 - *I5) + (*I9 - *I7));\n\t"
							"TI5 = (*I1 - C5QC*(*I3 + *I9)) + C5QB*(*R5 - *R7) - C5QD*(*R3 - *R9) + C5QA*((*I5 - *I3) + (*I7 - *I9));\n\t"
							"TI7 = (*I1 - C5QC*(*I3 + *I9)) - C5QB*(*R5 - *R7) + C5QD*(*R3 - *R9) + C5QA*((*I5 - *I3) + (*I7 - *I9));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0) = TR0 + TR1;\n\t"
							"(*R1) = TR2 + ( C5QE*TR3 + C5QD*TI3);\n\t"
							"(*R2) = TR4 + ( C5QA*TR5 + C5QB*TI5);\n\t"
							"(*R3) = TR6 + (-C5QA*TR7 + C5QB*TI7);\n\t"
							"(*R4) = TR8 + (-C5QE*TR9 + C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I0) = TI0 + TI1;\n\t"
							"(*I1) = TI2 + (-C5QD*TR3 + C5QE*TI3);\n\t"
							"(*I2) = TI4 + (-C5QB*TR5 + C5QA*TI5);\n\t"
							"(*I3) = TI6 + (-C5QB*TR7 - C5QA*TI7);\n\t"
							"(*I4) = TI8 + (-C5QD*TR9 - C5QE*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5) = TR0 - TR1;\n\t"
							"(*R6) = TR2 - ( C5QE*TR3 + C5QD*TI3);\n\t"
							"(*R7) = TR4 - ( C5QA*TR5 + C5QB*TI5);\n\t"
							"(*R8) = TR6 - (-C5QA*TR7 + C5QB*TI7);\n\t"
							"(*R9) = TR8 - (-C5QE*TR9 + C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I5) = TI0 - TI1;\n\t"
							"(*I6) = TI2 - (-C5QD*TR3 + C5QE*TI3);\n\t"
							"(*I7) = TI4 - (-C5QB*TR5 + C5QA*TI5);\n\t"
							"(*I8) = TI6 - (-C5QB*TR7 - C5QA*TI7);\n\t"
							"(*I9) = TI8 - (-C5QD*TR9 - C5QE*TI9);\n\t";
						}
					}
					else
					{
						if(cReg)
						{
							bflyStr +=
							"TR0 = (*R0).x + (*R2).x + (*R4).x + (*R6).x + (*R8).x;\n\t"
							"TR2 = ((*R0).x - C5QC*((*R4).x + (*R6).x)) - C5QB*((*R2).y - (*R8).y) - C5QD*((*R4).y - (*R6).y) + C5QA*(((*R2).x - (*R4).x) + ((*R8).x - (*R6).x));\n\t"
							"TR8 = ((*R0).x - C5QC*((*R4).x + (*R6).x)) + C5QB*((*R2).y - (*R8).y) + C5QD*((*R4).y - (*R6).y) + C5QA*(((*R2).x - (*R4).x) + ((*R8).x - (*R6).x));\n\t"
							"TR4 = ((*R0).x - C5QC*((*R2).x + (*R8).x)) + C5QB*((*R4).y - (*R6).y) - C5QD*((*R2).y - (*R8).y) + C5QA*(((*R4).x - (*R2).x) + ((*R6).x - (*R8).x));\n\t"
							"TR6 = ((*R0).x - C5QC*((*R2).x + (*R8).x)) - C5QB*((*R4).y - (*R6).y) + C5QD*((*R2).y - (*R8).y) + C5QA*(((*R4).x - (*R2).x) + ((*R6).x - (*R8).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = (*R0).y + (*R2).y + (*R4).y + (*R6).y + (*R8).y;\n\t"
							"TI2 = ((*R0).y - C5QC*((*R4).y + (*R6).y)) + C5QB*((*R2).x - (*R8).x) + C5QD*((*R4).x - (*R6).x) + C5QA*(((*R2).y - (*R4).y) + ((*R8).y - (*R6).y));\n\t"
							"TI8 = ((*R0).y - C5QC*((*R4).y + (*R6).y)) - C5QB*((*R2).x - (*R8).x) - C5QD*((*R4).x - (*R6).x) + C5QA*(((*R2).y - (*R4).y) + ((*R8).y - (*R6).y));\n\t"
							"TI4 = ((*R0).y - C5QC*((*R2).y + (*R8).y)) - C5QB*((*R4).x - (*R6).x) + C5QD*((*R2).x - (*R8).x) + C5QA*(((*R4).y - (*R2).y) + ((*R6).y - (*R8).y));\n\t"
							"TI6 = ((*R0).y - C5QC*((*R2).y + (*R8).y)) + C5QB*((*R4).x - (*R6).x) - C5QD*((*R2).x - (*R8).x) + C5QA*(((*R4).y - (*R2).y) + ((*R6).y - (*R8).y));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = (*R1).x + (*R3).x + (*R5).x + (*R7).x + (*R9).x;\n\t"
							"TR3 = ((*R1).x - C5QC*((*R5).x + (*R7).x)) - C5QB*((*R3).y - (*R9).y) - C5QD*((*R5).y - (*R7).y) + C5QA*(((*R3).x - (*R5).x) + ((*R9).x - (*R7).x));\n\t"
							"TR9 = ((*R1).x - C5QC*((*R5).x + (*R7).x)) + C5QB*((*R3).y - (*R9).y) + C5QD*((*R5).y - (*R7).y) + C5QA*(((*R3).x - (*R5).x) + ((*R9).x - (*R7).x));\n\t"
							"TR5 = ((*R1).x - C5QC*((*R3).x + (*R9).x)) + C5QB*((*R5).y - (*R7).y) - C5QD*((*R3).y - (*R9).y) + C5QA*(((*R5).x - (*R3).x) + ((*R7).x - (*R9).x));\n\t"
							"TR7 = ((*R1).x - C5QC*((*R3).x + (*R9).x)) - C5QB*((*R5).y - (*R7).y) + C5QD*((*R3).y - (*R9).y) + C5QA*(((*R5).x - (*R3).x) + ((*R7).x - (*R9).x));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = (*R1).y + (*R3).y + (*R5).y + (*R7).y + (*R9).y;\n\t"
							"TI3 = ((*R1).y - C5QC*((*R5).y + (*R7).y)) + C5QB*((*R3).x - (*R9).x) + C5QD*((*R5).x - (*R7).x) + C5QA*(((*R3).y - (*R5).y) + ((*R9).y - (*R7).y));\n\t"
							"TI9 = ((*R1).y - C5QC*((*R5).y + (*R7).y)) - C5QB*((*R3).x - (*R9).x) - C5QD*((*R5).x - (*R7).x) + C5QA*(((*R3).y - (*R5).y) + ((*R9).y - (*R7).y));\n\t"
							"TI5 = ((*R1).y - C5QC*((*R3).y + (*R9).y)) - C5QB*((*R5).x - (*R7).x) + C5QD*((*R3).x - (*R9).x) + C5QA*(((*R5).y - (*R3).y) + ((*R7).y - (*R9).y));\n\t"
							"TI7 = ((*R1).y - C5QC*((*R3).y + (*R9).y)) + C5QB*((*R5).x - (*R7).x) - C5QD*((*R3).x - (*R9).x) + C5QA*(((*R5).y - (*R3).y) + ((*R7).y - (*R9).y));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).x = TR0 + TR1;\n\t"
							"(*R1).x = TR2 + ( C5QE*TR3 - C5QD*TI3);\n\t"
							"(*R2).x = TR4 + ( C5QA*TR5 - C5QB*TI5);\n\t"
							"(*R3).x = TR6 + (-C5QA*TR7 - C5QB*TI7);\n\t"
							"(*R4).x = TR8 + (-C5QE*TR9 - C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0).y = TI0 + TI1;\n\t"
							"(*R1).y = TI2 + ( C5QD*TR3 + C5QE*TI3);\n\t"
							"(*R2).y = TI4 + ( C5QB*TR5 + C5QA*TI5);\n\t"
							"(*R3).y = TI6 + ( C5QB*TR7 - C5QA*TI7);\n\t"
							"(*R4).y = TI8 + ( C5QD*TR9 - C5QE*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5).x = TR0 - TR1;\n\t"
							"(*R6).x = TR2 - ( C5QE*TR3 - C5QD*TI3);\n\t"
							"(*R7).x = TR4 - ( C5QA*TR5 - C5QB*TI5);\n\t"
							"(*R8).x = TR6 - (-C5QA*TR7 - C5QB*TI7);\n\t"
							"(*R9).x = TR8 - (-C5QE*TR9 - C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5).y = TI0 - TI1;\n\t"
							"(*R6).y = TI2 - ( C5QD*TR3 + C5QE*TI3);\n\t"
							"(*R7).y = TI4 - ( C5QB*TR5 + C5QA*TI5);\n\t"
							"(*R8).y = TI6 - ( C5QB*TR7 - C5QA*TI7);\n\t"
							"(*R9).y = TI8 - ( C5QD*TR9 - C5QE*TI9);\n\t";
						}
						else
						{
							bflyStr +=
							"TR0 = *R0 + *R2 + *R4 + *R6 + *R8;\n\t"
							"TR2 = (*R0 - C5QC*(*R4 + *R6)) - C5QB*(*I2 - *I8) - C5QD*(*I4 - *I6) + C5QA*((*R2 - *R4) + (*R8 - *R6));\n\t"
							"TR8 = (*R0 - C5QC*(*R4 + *R6)) + C5QB*(*I2 - *I8) + C5QD*(*I4 - *I6) + C5QA*((*R2 - *R4) + (*R8 - *R6));\n\t"
							"TR4 = (*R0 - C5QC*(*R2 + *R8)) + C5QB*(*I4 - *I6) - C5QD*(*I2 - *I8) + C5QA*((*R4 - *R2) + (*R6 - *R8));\n\t"
							"TR6 = (*R0 - C5QC*(*R2 + *R8)) - C5QB*(*I4 - *I6) + C5QD*(*I2 - *I8) + C5QA*((*R4 - *R2) + (*R6 - *R8));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI0 = *I0 + *I2 + *I4 + *I6 + *I8;\n\t"
							"TI2 = (*I0 - C5QC*(*I4 + *I6)) + C5QB*(*R2 - *R8) + C5QD*(*R4 - *R6) + C5QA*((*I2 - *I4) + (*I8 - *I6));\n\t"
							"TI8 = (*I0 - C5QC*(*I4 + *I6)) - C5QB*(*R2 - *R8) - C5QD*(*R4 - *R6) + C5QA*((*I2 - *I4) + (*I8 - *I6));\n\t"
							"TI4 = (*I0 - C5QC*(*I2 + *I8)) - C5QB*(*R4 - *R6) + C5QD*(*R2 - *R8) + C5QA*((*I4 - *I2) + (*I6 - *I8));\n\t"
							"TI6 = (*I0 - C5QC*(*I2 + *I8)) + C5QB*(*R4 - *R6) - C5QD*(*R2 - *R8) + C5QA*((*I4 - *I2) + (*I6 - *I8));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TR1 = *R1 + *R3 + *R5 + *R7 + *R9;\n\t"
							"TR3 = (*R1 - C5QC*(*R5 + *R7)) - C5QB*(*I3 - *I9) - C5QD*(*I5 - *I7) + C5QA*((*R3 - *R5) + (*R9 - *R7));\n\t"
							"TR9 = (*R1 - C5QC*(*R5 + *R7)) + C5QB*(*I3 - *I9) + C5QD*(*I5 - *I7) + C5QA*((*R3 - *R5) + (*R9 - *R7));\n\t"
							"TR5 = (*R1 - C5QC*(*R3 + *R9)) + C5QB*(*I5 - *I7) - C5QD*(*I3 - *I9) + C5QA*((*R5 - *R3) + (*R7 - *R9));\n\t"
							"TR7 = (*R1 - C5QC*(*R3 + *R9)) - C5QB*(*I5 - *I7) + C5QD*(*I3 - *I9) + C5QA*((*R5 - *R3) + (*R7 - *R9));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"TI1 = *I1 + *I3 + *I5 + *I7 + *I9;\n\t"
							"TI3 = (*I1 - C5QC*(*I5 + *I7)) + C5QB*(*R3 - *R9) + C5QD*(*R5 - *R7) + C5QA*((*I3 - *I5) + (*I9 - *I7));\n\t"
							"TI9 = (*I1 - C5QC*(*I5 + *I7)) - C5QB*(*R3 - *R9) - C5QD*(*R5 - *R7) + C5QA*((*I3 - *I5) + (*I9 - *I7));\n\t"
							"TI5 = (*I1 - C5QC*(*I3 + *I9)) - C5QB*(*R5 - *R7) + C5QD*(*R3 - *R9) + C5QA*((*I5 - *I3) + (*I7 - *I9));\n\t"
							"TI7 = (*I1 - C5QC*(*I3 + *I9)) + C5QB*(*R5 - *R7) - C5QD*(*R3 - *R9) + C5QA*((*I5 - *I3) + (*I7 - *I9));\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R0) = TR0 + TR1;\n\t"
							"(*R1) = TR2 + ( C5QE*TR3 - C5QD*TI3);\n\t"
							"(*R2) = TR4 + ( C5QA*TR5 - C5QB*TI5);\n\t"
							"(*R3) = TR6 + (-C5QA*TR7 - C5QB*TI7);\n\t"
							"(*R4) = TR8 + (-C5QE*TR9 - C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I0) = TI0 + TI1;\n\t"
							"(*I1) = TI2 + ( C5QD*TR3 + C5QE*TI3);\n\t"
							"(*I2) = TI4 + ( C5QB*TR5 + C5QA*TI5);\n\t"
							"(*I3) = TI6 + ( C5QB*TR7 - C5QA*TI7);\n\t"
							"(*I4) = TI8 + ( C5QD*TR9 - C5QE*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*R5) = TR0 - TR1;\n\t"
							"(*R6) = TR2 - ( C5QE*TR3 - C5QD*TI3);\n\t"
							"(*R7) = TR4 - ( C5QA*TR5 - C5QB*TI5);\n\t"
							"(*R8) = TR6 - (-C5QA*TR7 - C5QB*TI7);\n\t"
							"(*R9) = TR8 - (-C5QE*TR9 - C5QD*TI9);\n\t";

							bflyStr += "\n\t";

							bflyStr +=
							"(*I5) = TI0 - TI1;\n\t"
							"(*I6) = TI2 - ( C5QD*TR3 + C5QE*TI3);\n\t"
							"(*I7) = TI4 - ( C5QB*TR5 + C5QA*TI5);\n\t"
							"(*I8) = TI6 - ( C5QB*TR7 - C5QA*TI7);\n\t"
							"(*I9) = TI8 - ( C5QD*TR9 - C5QE*TI9);\n\t";
						}
					}
				} break;
			case 11:
				{
					static const char *radix11str = " \
						fptype p0, p1, p2, p3, p4, p5, p6, p7, p8, p9; \n\
						p0 = ((*R1).x - (*R10).x)*dir; \n\
						p1 = (*R1).x + (*R10).x; \n\
						p2 = ((*R5).x - (*R6).x)*dir; \n\
						p3 = (*R5).x + (*R6).x; \n\
						p4 = ((*R2).x - (*R9).x)*dir; \n\
						p5 = (*R2).x + (*R9).x; \n\
						p6 = ((*R3).x - (*R8).x)*dir; \n\
						p7 = (*R3).x + (*R8).x; \n\
						p8 = (*R4).x + (*R7).x; \n\
						p9 = ((*R4).x - (*R7).x)*dir; \n\
						\n\
						fptype r0, r1, r2, r3, r4, r5, r6, r7, r8, r9; \n\
						r0 = p4 - p0 * b11_9; \n\
						r1 = p0 + p2 * b11_9; \n\
						r2 = p2 + p6 * b11_9; \n\
						r3 = p6 + p9 * b11_9; \n\
						r4 = p9 - p4 * b11_9; \n\
						r5 = p7 - p1 * b11_8; \n\
						r6 = p5 - p7 * b11_8; \n\
						r7 = p1 - p8 * b11_8; \n\
						r8 = p3 - p5 * b11_8; \n\
						r9 = p8 - p3 * b11_8; \n\
						\n\
						fptype s0, s1, s2, s3, s4, s5, s6, s7, s8, s9; \n\
						s0 = p6 - r0 * b11_6; \n\
						s1 = p9 + r1 * b11_6; \n\
						s2 = p4 - r2 * b11_6; \n\
						s3 = p0 + r3 * b11_6; \n\
						s4 = p2 + r4 * b11_6; \n\
						s5 = p3 - r5 * b11_7; \n\
						s6 = p8 - r6 * b11_7; \n\
						s7 = p5 - r7 * b11_7; \n\
						s8 = p1 - r8 * b11_7; \n\
						s9 = p7 - r9 * b11_7; \n\
						\n\
						fptype p10, p11, p12, p13, p14, p15, p16, p17, p18, p19; \n\
						p10 = ((*R10).y - (*R1).y)*dir; \n\
						p11 = (*R1).y + (*R10).y; \n\
						p12 = ((*R9).y - (*R2).y)*dir; \n\
						p13 = (*R2).y + (*R9).y; \n\
						p14 = ((*R8).y - (*R3).y)*dir; \n\
						p15 = (*R3).y + (*R8).y; \n\
						p16 = ((*R7).y - (*R4).y)*dir; \n\
						p17 = (*R4).y + (*R7).y; \n\
						p18 = ((*R6).y - (*R5).y)*dir; \n\
						p19 = (*R5).y + (*R6).y; \n\
						\n\
						fptype r10, r11, r12, r13, r14, r15, r16, r17, r18, r19; \n\
						r10 = p12 - p10 * b11_9; \n\
						r11 = p16 - p12 * b11_9; \n\
						r12 = p18 + p14 * b11_9; \n\
						r13 = p14 + p16 * b11_9; \n\
						r14 = p10 + p18 * b11_9; \n\
						r15 = p15 - p11 * b11_8; \n\
						r16 = p19 - p13 * b11_8; \n\
						r17 = p13 - p15 * b11_8; \n\
						r18 = p11 - p17 * b11_8; \n\
						r19 = p17 - p19 * b11_8; \n\
						\n\
						fptype s10, s11, s12, s13, s14, s15, s16, s17, s18, s19; \n\
						s10 = p14 - r10 * b11_6; \n\
						s11 = p18 + r11 * b11_6; \n\
						s12 = p12 - r12 * b11_6; \n\
						s13 = p10 + r13 * b11_6; \n\
						s14 = p16 + r14 * b11_6; \n\
						s15 = p19 - r15 * b11_7; \n\
						s16 = p11 - r16 * b11_7; \n\
						s17 = p17 - r17 * b11_7; \n\
						s18 = p13 - r18 * b11_7; \n\
						s19 = p15 - r19 * b11_7; \n\
						\n\
						fptype v0, v1, v2, v3, v4, v5, v6, v7, v8, v9; \n\
						fptype v10, v11, v12, v13, v14, v15, v16, v17, v18, v19; \n\
						v0 = p9 - s0 * b11_4; \n\
						v1 = p4 + s1 * b11_4; \n\
						v2 = p0 + s2 * b11_4; \n\
						v3 = p2 - s3 * b11_4; \n\
						v4 = p6 - s4 * b11_4; \n\
						v5 = p8 - s5 * b11_5; \n\
						v6 = p1 - s6 * b11_5; \n\
						v7 = p3 - s7 * b11_5; \n\
						v8 = p7 - s8 * b11_5; \n\
						v9 = p5 - s9 * b11_5; \n\
						v10 = p16 - s10 * b11_4; \n\
						v11 = p14 - s11 * b11_4; \n\
						v12 = p10 + s12 * b11_4; \n\
						v13 = p18 - s13 * b11_4; \n\
						v14 = p12 + s14 * b11_4; \n\
						v15 = p17 - s15 * b11_5; \n\
						v16 = p15 - s16 * b11_5; \n\
						v17 = p11 - s17 * b11_5; \n\
						v18 = p19 - s18 * b11_5; \n\
						v19 = p13 - s19 * b11_5; \n\
						\n\
						fptype w0, w1, w2, w3, w4, w5, w6, w7, w8, w9; \n\
						fptype w10, w11, w12, w13, w14, w15, w16, w17, w18, w19; \n\
						w0 = p2 - v0 * b11_2; \n\
						w1 = p6 + v1 * b11_2; \n\
						w2 = p9 - v2 * b11_2; \n\
						w3 = p4 + v3 * b11_2; \n\
						w4 = p0 - v4 * b11_2; \n\
						w5 = p5 - v5 * b11_3; \n\
						w6 = p3 - v6 * b11_3; \n\
						w7 = p7 - v7 * b11_3; \n\
						w8 = p8 - v8 * b11_3; \n\
						w9 = p1 - v9 * b11_3; \n\
						w10 = p18 - v10 * b11_2; \n\
						w11 = p10 - v11 * b11_2; \n\
						w12 = p16 - v12 * b11_2; \n\
						w13 = p12 + v13 * b11_2; \n\
						w14 = p14 + v14 * b11_2; \n\
						w15 = p13 - v15 * b11_3; \n\
						w16 = p17 - v16 * b11_3; \n\
						w17 = p19 - v17 * b11_3; \n\
						w18 = p15 - v18 * b11_3; \n\
						w19 = p11 - v19 * b11_3; \n\
						\n\
						fptype z0, z1, z2, z3, z4, z5, z6, z7, z8, z9; \n\
						z0 = (*R0).x - w5 * b11_1; \n\
						z1 = (*R0).x - w6 * b11_1; \n\
						z2 = (*R0).x - w7 * b11_1; \n\
						z3 = (*R0).x - w8 * b11_1; \n\
						z4 = (*R0).x - w9 * b11_1; \n\
						z5 = (*R0).y - w15 * b11_1; \n\
						z6 = (*R0).y - w16 * b11_1; \n\
						z7 = (*R0).y - w17 * b11_1; \n\
						z8 = (*R0).y - w18 * b11_1; \n\
						z9 = (*R0).y - w19 * b11_1; \n\
						\n\
						(*R0).x = (*R0).x + p1 + p3 + p5 + p7 + p8; \n\
						(*R0).y = (*R0).y + p11 + p13 + p15 + p17 + p19; \n\
						(*R1).x = z1 + w14* b11_0; \n\
						(*R1).y = z7 + w1* b11_0; \n\
						(*R2).x = z2 - w12* b11_0; \n\
						(*R2).y = z8 - w2* b11_0; \n\
						(*R3).x = z0 + w11* b11_0; \n\
						(*R3).y = z5 + w4* b11_0; \n\
						(*R4).x = z3 - w13* b11_0; \n\
						(*R4).y = z6 - w3* b11_0; \n\
						(*R5).x = z4 + w10* b11_0; \n\
						(*R5).y = z9 + w0* b11_0; \n\
						(*R6).x = z4 - w10* b11_0; \n\
						(*R6).y = z9 - w0* b11_0; \n\
						(*R7).x = z3 + w13* b11_0; \n\
						(*R7).y = z6 + w3* b11_0; \n\
						(*R8).x = z0 - w11* b11_0; \n\
						(*R8).y = z5 - w4* b11_0; \n\
						(*R9).x = z2 + w12* b11_0; \n\
						(*R9).y = z8 + w2* b11_0; \n\
						(*R10).x = z1 - w14* b11_0; \n\
						(*R10).y = z7 - w1* b11_0; \n";

					if (fwd)
					{
						bflyStr += "fptype dir = -1;\n\n";
					}
					else
					{
						bflyStr += "fptype dir = 1;\n\n";
					}

					bflyStr += radix11str;

				} break;
			case 13:
				{

					static const char *radix13str = " \
						fptype p0, p1, p2, p3, p4, p5, p6, p7, p8, p9;\n\
						p0 = (*R7).x - (*R2).x;\n\
						p1 = (*R7).x + (*R2).x;\n\
						p2 = (*R8).x - (*R5).x;\n\
						p3 = (*R8).x + (*R5).x;\n\
						p4 = (*R9).x - (*R3).x;\n\
						p5 = (*R3).x + (*R9).x;\n\
						p6 = (*R10).x + (*R4).x;\n\
						p7 = (*R10).x - (*R4).x;\n\
						p8 = (*R11).x + (*R6).x;\n\
						p9 = (*R11).x - (*R6).x;\n\
						\n\
						fptype p10, p11, p12, p13, p14, p15, p16, p17, p18, p19;\n\
						p10 = (*R12).x + p6;\n\
						p11 = (*R1).x + p5;\n\
						p12 = p8 - p1;\n\
						p13 = p8 + p1;\n\
						p14 = p9 + p0;\n\
						p15 = p9 - p0;\n\
						p16 = p7 - p4;\n\
						p17 = p4 + p7;\n\
						p18 = p11 + p10;\n\
						p19 = p11 - p10;\n\
						\n\
						fptype s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11;\n\
						s0 = p3 + p13;\n\
						s1 = p2 + p14;\n\
						s2 = p16 - p15;\n\
						s3 = p16 + p15;\n\
						s4 = -(*R12).x + p6 * b13_17;\n\
						s5 =   (*R1).x - p5 * b13_17;\n\
						s6 = s5 - s4;\n\
						s7 = s5 + s4;\n\
						s8 = p18 + s0;\n\
						s9 = p18 - s0;\n\
						fptype c2 = p3 - p13 * b13_17;\n\
						s10 = s6 - c2;\n\
						s11 = s6 + c2;\n\
						\n\
						fptype r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11;\n\
						r0 = (*R7).y + (*R2).y;\n\
						r1 = (*R7).y - (*R2).y;\n\
						r2 = (*R8).y + (*R5).y;\n\
						r3 = (*R8).y - (*R5).y;\n\
						r4 = (*R9).y - (*R3).y;\n\
						r5 = (*R3).y + (*R9).y;\n\
						r6 = (*R10).y + (*R4).y;\n\
						r7 = (*R10).y - (*R4).y;\n\
						r8 = (*R11).y - (*R6).y;\n\
						r9 = (*R11).y + (*R6).y;\n\
						r10 = (*R12).y + r6;\n\
						r11 = (*R1).y + r5;\n\
						\n\
						fptype m0, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10;\n\
						fptype m11, m12, m13, m14, m15, m16, m17, m18, m19, m20;\n\
						m0 = r4 + r7;\n\
						m1 = r7 - r4;\n\
						m2 = r8 - r1;\n\
						m3 = r8 + r1;\n\
						m4 = r9 + r0;\n\
						m5 = r9 - r0;\n\
						m6 = r11 + r10;\n\
						m7 = r11 - r10;\n\
						m8 = m1 - m2;\n\
						m9 = m1 + m2;\n\
						m10 = r3 + m3;\n\
						m11 = r2 + m4;\n\
						m12 = m6 - m11;\n\
						m13 = m6 + m11;\n\
						\n\
						m14 =  (*R1).y - r5 * b13_17;\n\
						m15 = -(*R12).y + r6 * b13_17;\n\
						m16 =  r2      - m4 * b13_17;\n\
						\n\
						m17 = m14 + m15;\n\
						m18 = m14 - m15;\n\
						m19 = m18 + m16;\n\
						m20 = m18 - m16;\n\
						\n\
						fptype c0, c1, c3, c4, c5, c6, c7, c8, c9;\n\
						fptype c10, c11, c12, c13, c14, c15, c16, c17, c18, c19;\n\
						fptype c20, c21, c22, c23, c24;\n\
						c0  =  s7 - p12 * b13_3;\n\
						c1  =  s7 + p12 * b13_3;\n\
						c3  =  p2 - p14 * b13_17;\n\
						c4  =  s1 - p19 * b13_18;\n\
						c5  = p19 + s1 * b13_18;\n\
						c6  = s10 - s2 * b13_15;\n\
						c7  = s11 - s3 * b13_22;\n\
						c8  = (*R0).x - s8 * b13_23;\n\
						c9  =  s2 + s10 * b13_7;\n\
						c10 =  s3 + s11 * b13_19;\n\
						c11 =  r3 - m3 * b13_17;\n\
						c12 = m17 - m5 * b13_3;\n\
						c13 = m17 + m5 * b13_3;\n\
						c14 = m10 - m7 * b13_18;\n\
						c15 = m20 - m8 * b13_15;\n\
						c16 = m19 - m9 * b13_22;\n\
						c17 =  m7 + m10 * b13_18;\n\
						c18 = (*R0).y- m13 * b13_23;\n\
						c19 =  m9 + m19 * b13_19;\n\
						c20 =  m8 + m20 * b13_7;\n\
						c21 =  c3 + p17 * b13_3;\n\
						c22 =  c3 - p17 * b13_3;\n\
						c23 = c11 + m0 * b13_3;\n\
						c24 = c11 - m0 * b13_3;\n\
						\n\
						fptype d0, d1, d2, d3, d4, d5, d6, d7, d8, d9;\n\
						fptype d10, d11, d12, d13, d14, d15, d16, d17, d18, d19;\n\
						d0  = c22 +  c0 * b13_8;\n\
						d1  =  c0 - c22 * b13_8;\n\
						d2  = c21 +  c1 * b13_24;\n\
						d3  =  c1 - c21 * b13_24;\n\
						d4  =  s9 -  c6 * b13_4;\n\
						d5  =  c6 +  s9 * b13_10;\n\
						d6  =  c7 +  c9 * b13_6;\n\
						d7  =  c7 -  c9 * b13_6;\n\
						d8  =  c8 - c10 * b13_21;\n\
						d9  =  c8 + c10 * b13_16;\n\
						d10 = c24 + c12 * b13_8;\n\
						d11 = c12 - c24 * b13_8;\n\
						d12 = c23 + c13 * b13_24;\n\
						d13 = c13 - c23 * b13_24;\n\
						d14 = m12 - c15 * b13_4;\n\
						d15 = c15 + m12 * b13_10;\n\
						d16 = c18 + c19 * b13_16;\n\
						d17 = c18 - c19 * b13_21;\n\
						d18 = c16 - c20 * b13_6;\n\
						d19 = c16 + c20 * b13_6;\n\
						\n\
						fptype e0, e1, e2, e3, e4, e5, e6, e7, e8, e9;\n\
						fptype e10, e11, e12, e13, e14, e15;\n\
						e0  = d2  +  d0 * b13_5;\n\
						e1  = d2  -  d0 * b13_5;\n\
						e2  = d3  -  d1 * b13_5;\n\
						e3  = d3  +  d1 * b13_5;\n\
						e4  = d8  -  d4 * b13_20;\n\
						e5  = d8  +  d4 * b13_20;\n\
						e6  = d9  +  d5 * b13_14;\n\
						e7  = d9  -  d5 * b13_14;\n\
						e8  = d12 + d10 * b13_5;\n\
						e9  = d12 - d10 * b13_5;\n\
						e10 = d13 - d11 * b13_5;\n\
						e11 = d13 + d11 * b13_5;\n\
						e12 = d16 + d15 * b13_14;\n\
						e13 = d16 - d15 * b13_14;\n\
						e14 = d17 + d14 * b13_20;\n\
						e15 = d17 - d14 * b13_20;\n\
						\n\
						fptype f0, f1, f2, f3, f4, f5, f6, f7, f8, f9;\n\
						fptype f10, f11, f12, f13, f14, f15, f16, f17, f18, f19;\n\
						fptype f20, f21, f22, f23;\n\
						f0  = c17 - e10 * b13_12;\n\
						f1  = e10 + c17 * b13_1;\n\
						f2  = e9  + c14 * b13_1;\n\
						f3  = c14 -  e9 * b13_12;\n\
						f4  = e11 + dir * d7 * b13_0;\n\
						f5  = e11 - dir * d7 * b13_0;\n\
						f6  = e5  + dir * f3 * b13_11;\n\
						f7  = e5  - dir * f3 * b13_11;\n\
						f8  = e4  + dir * e8 * b13_13;\n\
						f9  = e4  - dir * e8 * b13_13;\n\
						f10 = f0  - dir * d6 * b13_2;\n\
						f11 = f0  + dir * d6 * b13_2;\n\
						f12 = e1  +  c4 * b13_1;\n\
						f13 = c4  -  e1 * b13_12;\n\
						f14 = c5  -  e2 * b13_12;\n\
						f15 = e2  +  c5 * b13_1;\n\
						f16 = f14 + dir * d19 * b13_2;\n\
						f17 = f14 - dir * d19 * b13_2;\n\
						f18 = e15 - dir *  e0 * b13_13;\n\
						f19 = e15 + dir *  e0 * b13_13;\n\
						f20 = e14 - dir * f13 * b13_11;\n\
						f21 = e14 + dir * f13 * b13_11;\n\
						f22 = e3  - dir * d18 * b13_0;\n\
						f23 = e3  + dir * d18 * b13_0;\n\
						\n\
						(*R0).x  = (*R0).x + s8;\n\
						(*R0).y  = (*R0).y + m13;\n\
						(*R1).x  =  e6 +  f2 * dir * b13_9 ;\n\
						(*R1).y  = e12 - f12 * dir * b13_9 ;\n\
						(*R2).x  =  f9 - f10 * dir * b13_11;\n\
						(*R2).y  = f19 + f16 * dir * b13_11;\n\
						(*R3).x  =  f6 -  f5 * dir * b13_13;\n\
						(*R3).y  = f20 + f23 * dir * b13_13;\n\
						(*R4).x  =  f7 -  f4 * dir * b13_13;\n\
						(*R4).y  = f21 + f22 * dir * b13_13;\n\
						(*R5).x  =  e7 -  f1 * dir * b13_9 ;\n\
						(*R5).y  = e13 + f15 * dir * b13_9 ;\n\
						(*R6).x  =  f8 - f11 * dir * b13_11;\n\
						(*R6).y  = f18 + f17 * dir * b13_11;\n\
						(*R7).x  =  f9 + f10 * dir * b13_11;\n\
						(*R7).y  = f19 - f16 * dir * b13_11;\n\
						(*R8).x  =  e7 +  f1 * dir * b13_9 ;\n\
						(*R8).y  = e13 - f15 * dir * b13_9 ;\n\
						(*R9).x  =  f6 +  f5 * dir * b13_13;\n\
						(*R9).y  = f20 - f23 * dir * b13_13;\n\
						(*R10).x =  f7 +  f4 * dir * b13_13;\n\
						(*R10).y = f21 - f22 * dir * b13_13;\n\
						(*R11).x =  f8 + f11 * dir * b13_11;\n\
						(*R11).y = f18 - f17 * dir * b13_11;\n\
						(*R12).x =  e6 -  f2 * dir * b13_9 ;\n\
						(*R12).y = e12 + f12 * dir * b13_9 ;\n";

						if (fwd)
						{
							bflyStr += "fptype dir = -1;\n\n";
						}
						else
						{
							bflyStr += "fptype dir = 1;\n\n";
						}

						bflyStr += radix13str;

				} break;

			default:
				assert(false);
			}

			bflyStr += "\n\t";

			// Assign results
			if( (radix & (radix-1)) || (!cReg) )
			{
				if( (radix != 10) && (radix != 6) )
				{
				for(size_t i=0; i<radix;i++)
				{
					if(cReg)
					{
						if ( (radix != 7) && (radix != 11) && (radix != 13) )
						{
						bflyStr += "((*R"; bflyStr += SztToStr(i); bflyStr += ").x) = TR"; bflyStr += SztToStr(i); bflyStr += "; ";
						bflyStr += "((*R"; bflyStr += SztToStr(i); bflyStr += ").y) = TI"; bflyStr += SztToStr(i); bflyStr += ";\n\t";
						}
					}
					else
					{
						bflyStr += "(*R"; bflyStr += SztToStr(i); bflyStr += ") = TR"; bflyStr += SztToStr(i); bflyStr += "; ";
						bflyStr += "(*I"; bflyStr += SztToStr(i); bflyStr += ") = TI"; bflyStr += SztToStr(i); bflyStr += ";\n\t";
					}
				}
				}
			}
			else
			{
				for(size_t i=0; i<radix;i++)
				{
					size_t j = BitReverse(i, radix);

					if(i < j)
					{
						bflyStr += "T = (*R"; bflyStr += SztToStr(i); bflyStr += "); (*R";
						bflyStr += SztToStr(i); bflyStr += ") = (*R"; bflyStr += SztToStr(j); bflyStr += "); (*R";
						bflyStr += SztToStr(j); bflyStr += ") = T;\n\t";
					}
				}
			}

			bflyStr += "\n}\n";
		}

	public:
		Butterfly(size_t radixVal, size_t countVal, bool fwdVal, bool cRegVal) : radix(radixVal), count(countVal), fwd(fwdVal), cReg(cRegVal) {}

		void GenerateButterfly(std::string &bflyStr) const
		{
			assert(count <= 4);
			if(count > 0)
				GenerateButterflyStr(bflyStr);
		}
    };

};

#endif

