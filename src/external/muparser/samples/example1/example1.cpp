/*

	 _____  __ _____________ _______  ______ ___________
	/     \|  |  \____ \__  \\_  __ \/  ___// __ \_  __ \
   |  Y Y  \  |  /  |_> > __ \|  | \/\___ \\  ___/|  | \/
   |__|_|  /____/|   __(____  /__|  /____  >\___  >__|
		 \/      |__|       \/           \/     \/
   Copyright (C) 2022 Ingo Berg

	Redistribution and use in source and binary forms, with or without modification, are permitted
	provided that the following conditions are met:

	  * Redistributions of source code must retain the above copyright notice, this list of
		conditions and the following disclaimer.
	  * Redistributions in binary form must reproduce the above copyright notice, this list of
		conditions and the following disclaimer in the documentation and/or other materials provided
		with the distribution.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
	IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
	FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
	IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
	OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <locale>
#include <limits>
#include <ios> 
#include <iomanip>
#include <numeric>

#include "muParserTest.h"
#include "muParser.h"

using namespace std;
using namespace mu;


// Forward declarations
void CalcBulk();

// Operator callback functions
static value_type Mega(value_type a_fVal) { return a_fVal * 1e6; }
static value_type Milli(value_type a_fVal) { return a_fVal / (value_type)1e3; }
static value_type Rnd(value_type v) { return v * std::rand() / (value_type)(RAND_MAX + 1.0); }
static value_type Not(value_type v) { return v == 0; }
static value_type Add(value_type v1, value_type v2) { return v1 + v2; }
static value_type Mul(value_type v1, value_type v2) { return v1 * v2; }
static value_type Arg2Of2(value_type /* v1 */, value_type v2) { return v2; }
static value_type Arg1Of2(value_type v1, value_type /* v2 */) { return v1; }


static value_type ThrowAnException(value_type)
{
	throw std::runtime_error("This function does throw an exception.");
}


static value_type BulkFun1(int nBulkIdx, int nThreadIdx, value_type v1)
{
	// Note: I'm just doing something with all three parameters to shut 
	// compiler warnings up!
	return (value_type)nBulkIdx + nThreadIdx + v1;
}


static value_type Ping()
{
	mu::console() << "ping\n";
	return 0;
}


static value_type StrFun0(const char_type* szMsg)
{
	if (szMsg)
		mu::console() << szMsg << std::endl;

	return 999;
}


static value_type StrFun2(const char_type* v1, value_type v2, value_type v3)
{
	mu::console() << v1 << std::endl;
	return v2 + v3;
}


static value_type Debug(mu::value_type v1, mu::value_type v2)
{
	ParserBase::EnableDebugDump(v1 != 0, v2 != 0);
	mu::console() << _T("Bytecode dumping ") << ((v1 != 0) ? _T("active") : _T("inactive")) << _T("\n");
	return 1;
}


// Factory function for creating new parser variables
// This could as well be a function performing database queries.
static value_type* AddVariable(const char_type* a_szName, void* a_pUserData)
{
	// I don't want dynamic allocation here, so i used this static buffer
	// If you want dynamic allocation you must allocate all variables dynamically
	// in order to delete them later on. Or you find other ways to keep track of 
	// variables that have been created implicitely.
	static value_type afValBuf[100];
	static int iVal = -1;

	++iVal;

	mu::console()
		<< _T("Generating new variable \"")
		<< a_szName << std::dec << _T("\" (slots left: ")
		<< 99 - iVal << _T(")")
		<< _T(" User data pointer is:")
		<< std::hex << a_pUserData << endl;

	afValBuf[iVal] = 0;

	if (iVal >= 99)
		throw mu::ParserError(_T("Variable buffer overflow."));
	else
		return &afValBuf[iVal];
}


int IsBinValue(const char_type* a_szExpr, int* a_iPos, value_type* a_fVal)
{
	if (a_szExpr[0] != 0 && a_szExpr[1] != 'b')
		return 0;

	unsigned iVal = 0;
	unsigned iBits = sizeof(iVal) * 8;
	unsigned i = 0;

	for (i = 0; (a_szExpr[i + 2] == '0' || a_szExpr[i + 2] == '1') && i < iBits; ++i)
		iVal |= (int)(a_szExpr[i + 2] == '1') << ((iBits - 1) - i);

	if (i == 0)
		return 0;

	if (i == iBits)
		throw mu::Parser::exception_type(_T("Binary to integer conversion error (overflow)."));

	*a_fVal = (unsigned)(iVal >> (iBits - i));
	*a_iPos += i + 2;

	return 1;
}


static int IsHexValue(const char_type* a_szExpr, int* a_iPos, value_type* a_fVal)
{
	if (a_szExpr[1] == 0 || (a_szExpr[0] != '0' || a_szExpr[1] != 'x'))
		return 0;

	unsigned iVal(0);

	// New code based on streams for UNICODE compliance:
	stringstream_type::pos_type nPos(0);
	stringstream_type ss(a_szExpr + 2);
	ss >> std::hex >> iVal;
	nPos = ss.tellg();

	if (nPos == (stringstream_type::pos_type)0)
		return 1;

	*a_iPos += (int)(2 + nPos);
	*a_fVal = (value_type)iVal;

	return 1;
}


static void Splash()
{
	mu::console() << _T("\n");
	mu::console() << _T(R"(   _____  __ _____________ ________  _____ ____________  )") << _T("\n");
	mu::console() << _T(R"(  /     \|  |  \____ \__   \\_  __ \/ ___// __  \_  __ \ )") << _T("\n");
	mu::console() << _T(R"( |  Y Y  \  |  /  |_> > ___ \|  | \/\___\\  ___/ |  | \/ )") << _T("\n");
	mu::console() << _T(R"( |__|_|  /____/|   __(____  /___|  /___  >\___  >|__|    )") << _T("\n");
	mu::console() << _T(R"(       \/      |__|       \/           \/     \/        )") << _T("\n");
	mu::console() << _T("  Version ") << Parser().GetVersion(pviFULL) << _T("\n");
	mu::console() << _T("  (C) 2022 Ingo Berg\n");
	mu::console() << _T("\n");
	mu::console() << _T("-----------------------------------------------------------\n");

#if defined(__clang__)
	// Note: CLANG also identifies as GCC 4.2.1
	mu::console() << _T("  Compiled with CLANG Version ") << __clang_major__ << _T(".") << __clang_minor__ << _T(".") << __clang_patchlevel__ << _T("\n");
#elif defined (__GNUC__)
	mu::console() << _T("  Compiled with GCC Version ") << __GNUC__ << _T(".") << __GNUC_MINOR__ << _T(".") << __GNUC_PATCHLEVEL__ << _T("\n");
#elif defined(_MSC_VER)
	mu::console() << _T("  Compiled with MSVC Version ") << _MSC_VER << _T("\n");
#endif

	mu::console() << _T("  IEEE 754 (IEC 559) is ") << ((std::numeric_limits<double>::is_iec559) ? "Available" : " NOT AVAILABLE") << _T("\n");
	mu::console() << _T("  ") << sizeof(void*) * 8 << _T("-bit build\n");
}


static value_type SelfTest()
{
	mu::console() << _T("-----------------------------------------------------------\n");
	mu::console() << _T("Running unit tests:\n\n");

	// Skip the self test if the value type is set to an integer type.
	if (mu::TypeInfo<mu::value_type>::IsInteger())
	{
		mu::console() << _T("  Test skipped: integer data type are not compatible with the unit test!\n\n");
	}
	else
	{
		mu::Test::ParserTester pt;
		pt.Run();
	}

	return 0;
}


static value_type Help()
{
	mu::console() << _T("-----------------------------------------------------------\n");
	mu::console() << _T("Commands:\n\n");
	mu::console() << _T("  list var     - list parser variables\n");
	mu::console() << _T("  list exprvar - list expression variables\n");
	mu::console() << _T("  list const   - list all numeric parser constants\n");
	mu::console() << _T("  opt on       - enable optimizer (default)\n");
	mu::console() << _T("  opt off      - disable optimizer\n");
	mu::console() << _T("  locale de    - switch to german locale\n");
	mu::console() << _T("  locale en    - switch to english locale\n");
	mu::console() << _T("  locale reset - reset locale\n");
	mu::console() << _T("  test bulk    - test bulk mode\n");
	mu::console() << _T("  quit         - exits the parser\n");
	mu::console() << _T("\nConstants:\n\n");
	mu::console() << _T("  \"_e\"   2.718281828459045235360287\n");
	mu::console() << _T("  \"_pi\"  3.141592653589793238462643\n");
	mu::console() << _T("-----------------------------------------------------------\n");
	return 0;
}


static void ListVar(const mu::ParserBase& parser)
{
	// Query the used variables (must be done after calc)
	mu::varmap_type variables = parser.GetVar();
	if (!variables.size())
		return;

	cout << "\nParser variables:\n";
	cout << "-----------------\n";
	cout << "Number: " << (int)variables.size() << "\n";
	varmap_type::const_iterator item = variables.begin();
	for (; item != variables.end(); ++item)
		mu::console() << _T("Name: ") << item->first << _T("   Address: [0x") << item->second << _T("]\n");
}


static void ListConst(const mu::ParserBase& parser)
{
	mu::console() << _T("\nParser constants:\n");
	mu::console() << _T("-----------------\n");

	mu::valmap_type cmap = parser.GetConst();
	if (!cmap.size())
	{
		mu::console() << _T("Expression does not contain constants\n");
	}
	else
	{
		valmap_type::const_iterator item = cmap.begin();
		for (; item != cmap.end(); ++item)
			mu::console() << _T("  ") << item->first << _T(" =  ") << item->second << _T("\n");
	}
}


static void ListExprVar(const mu::ParserBase& parser)
{
	string_type sExpr = parser.GetExpr();
	if (sExpr.length() == 0)
	{
		cout << _T("Expression string is empty\n");
		return;
	}

	// Query the used variables (must be done after calc)
	mu::console() << _T("\nExpression variables:\n");
	mu::console() << _T("---------------------\n");
	mu::console() << _T("Expression: ") << parser.GetExpr() << _T("\n");

	varmap_type variables = parser.GetUsedVar();
	if (!variables.size())
	{
		mu::console() << _T("Expression does not contain variables\n");
	}
	else
	{
		mu::console() << _T("Number: ") << (int)variables.size() << _T("\n");
		mu::varmap_type::const_iterator item = variables.begin();
		for (; item != variables.end(); ++item)
			mu::console() << _T("Name: ") << item->first << _T("   Address: [0x") << item->second << _T("]\n");
	}
}


/** \brief Check for external keywords.
*/
static int CheckKeywords(const mu::char_type* a_szLine, mu::Parser& a_Parser)
{
	string_type sLine(a_szLine);

	if (sLine == _T("quit"))
	{
		return -1;
	}
	else if (sLine == _T("list var"))
	{
		ListVar(a_Parser);
		return 1;
	}
	else if (sLine == _T("opt on"))
	{
		a_Parser.EnableOptimizer(true);
		mu::console() << _T("Optimizer enabled\n");
		return 1;
	}
	else if (sLine == _T("opt off"))
	{
		a_Parser.EnableOptimizer(false);
		mu::console() << _T("Optimizer disabled\n");
		return 1;
	}
	else if (sLine == _T("list const"))
	{
		ListConst(a_Parser);
		return 1;
	}
	else if (sLine == _T("list exprvar"))
	{
		ListExprVar(a_Parser);
		return 1;
	}
	else if (sLine == _T("locale de"))
	{
		mu::console() << _T("Setting german locale: ArgSep=';' DecSep=',' ThousandsSep='.'\n");
		a_Parser.SetArgSep(';');
		a_Parser.SetDecSep(',');
		a_Parser.SetThousandsSep('.');
		return 1;
	}
	else if (sLine == _T("locale en"))
	{
		mu::console() << _T("Setting english locale: ArgSep=',' DecSep='.' ThousandsSep=''\n");
		a_Parser.SetArgSep(',');
		a_Parser.SetDecSep('.');
		a_Parser.SetThousandsSep();
		return 1;
	}
	else if (sLine == _T("locale reset"))
	{
		mu::console() << _T("Resetting locale\n");
		a_Parser.ResetLocale();
		return 1;
	}
	else if (sLine == _T("test bulk"))
	{
		mu::console() << _T("Testing bulk mode\n");
		CalcBulk();
		return 1;
	}
	else if (sLine == _T("dbg"))
	{
		string_type dbg = _T("((\"\")), 7");
		a_Parser.SetExpr(dbg);
		mu::console() << dbg;

		int stackSize;
		double* v = a_Parser.Eval(stackSize);
		mu::console() << "=" <<  *v << std::endl;
		return 1;
	}

	return 0;
}


void CalcBulk()
{
	const int nBulkSize = 200;
	value_type* x = new value_type[nBulkSize];
	value_type* y = new value_type[nBulkSize];
	value_type* result = new value_type[nBulkSize];

	try
	{
		for (int i = 0; i < nBulkSize; ++i)
		{
			x[i] = i;
			y[i] = (value_type)i / 10;
		}
		mu::Parser  parser;
		parser.DefineVar(_T("x"), x);
		parser.DefineVar(_T("y"), y);
		parser.DefineFun(_T("fun1"), BulkFun1);
		parser.SetExpr(_T("fun1(0)+x+y"));
		parser.Eval(result, nBulkSize);

		for (int i = 0; i < nBulkSize; ++i)
		{
			mu::console() << _T("Eqn. ") << i << _T(": x=") << x[i] << _T("; y=") << y[i] << _T("; result=") << result[i] << _T("\n");
		}
	}
	catch (...)
	{
		delete[] x;
		delete[] y;
		delete[] result;
		throw;
	}

	delete[] x;
	delete[] y;
	delete[] result;
}


static void Calc()
{
	mu::Parser  parser;

	// Add some variables
	value_type  vVarVal[] = { 1, 2 }; // Values of the parser variables
	parser.DefineVar(_T("a"), &vVarVal[0]);  // Assign Variable names and bind them to the C++ variables
	parser.DefineVar(_T("b"), &vVarVal[1]);
	parser.DefineVar(_T("ft"), &vVarVal[1]);
	parser.DefineStrConst(_T("sVar1"), _T("Sample string 1"));
	parser.DefineStrConst(_T("sVar2"), _T("Sample string 2"));
	parser.AddValIdent(IsHexValue);
	parser.AddValIdent(IsBinValue);

	// Add user defined unary operators
	parser.DefinePostfixOprt(_T("M"), Mega);
	parser.DefinePostfixOprt(_T("m"), Milli);
	parser.DefineInfixOprt(_T("!"), Not);
	parser.DefineFun(_T("strfun0"), StrFun0);
	parser.DefineFun(_T("strfun2"), StrFun2);
	parser.DefineFun(_T("ping"), Ping);
	parser.DefineFun(_T("rnd"), Rnd, false);     // Add an unoptimizeable function
	parser.DefineFun(_T("throw"), ThrowAnException);

	parser.DefineOprt(_T("add"), Add, 0);
	parser.DefineOprt(_T("mul"), Mul, 1);

	// These are service and debug functions
	parser.DefineFun(_T("debug"), Debug);
	parser.DefineFun(_T("selftest"), SelfTest);
	parser.DefineFun(_T("help"), Help);
	parser.DefineFun(_T("arg2of2"), Arg2Of2);
	parser.DefineFun(_T("arg1of2"), Arg1Of2, false);

	parser.DefinePostfixOprt(_T("{ft}"), Milli);
	parser.DefinePostfixOprt(_T("ft"), Milli);

	// Define the variable factory
	parser.SetVarFactory(AddVariable, &parser);

	for (;;)
	{
		try
		{
			string_type sLine;
			std::getline(mu::console_in(), sLine);

			switch (CheckKeywords(sLine.c_str(), parser))
			{
			case  0: break;
			case  1: continue;
			case -1: return;
			}

			if (!sLine.length())
				continue;

			parser.SetExpr(sLine);
			mu::console() << std::setprecision(12);

			// There are multiple ways to retrieve the result...
			// 1.) If you know there is only a single return value or in case you only need the last 
			//     result of an expression consisting of comma separated subexpressions you can 
			//     simply use: 
			mu::console() << _T("ans=") << parser.Eval() << _T("\n");

			// 2.) As an alternative you can also retrieve multiple return values using this API:
			int nNum = parser.GetNumResults();
			if (nNum > 1)
			{
				mu::console() << _T("Multiple return values detected! Complete list:\n");

				// this is the hard way if you need to retrieve multiple subexpression
				// results
				value_type* v = parser.Eval(nNum);
				mu::console() << std::setprecision(12);
				for (int i = 0; i < nNum; ++i)
				{
					mu::console() << v[i] << _T("\n");
				}
			}
		}
		catch (mu::Parser::exception_type& e)
		{
			mu::console() << _T("\nError:\n");
			mu::console() << _T("------\n");
			mu::console() << _T("Message:     ") << e.GetMsg() << _T("\n");
			mu::console() << _T("Expression:  \"") << e.GetExpr() << _T("\"\n");
			mu::console() << _T("Token:       \"") << e.GetToken() << _T("\"\n");
			mu::console() << _T("Position:    ") << (int)e.GetPos() << _T("\n");
			mu::console() << _T("Errc:        ") << std::dec << e.GetCode() << _T("\n");
		}
	} // while running
}


int main(int, char**)
{
	Splash();
	SelfTest();
	Help();

	mu::console() << _T("Enter an expression or a command:\n");

	try
	{
		Calc();
	}
	catch (Parser::exception_type& e)
	{
		// Only erros raised during the initialization will end up here
		// formula related errors are treated in Calc()
		console() << _T("Initialization error:  ") << e.GetMsg() << endl;
		console() << _T("aborting...") << endl;
		string_type sBuf;
		console_in() >> sBuf;
	}
	catch (std::exception& /*exc*/)
	{
		// there is no unicode compliant way to query exc.what()
		// i'll leave it for this example.
		console() << _T("aborting...\n");
	}

	return 0;
}
