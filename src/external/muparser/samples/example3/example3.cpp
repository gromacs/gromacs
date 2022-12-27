/*

     _____  __ _____________ _______  ______ ___________
    /     \|  |  \____ \__  \\_  __ \/  ___// __ \_  __ \
   |  Y Y  \  |  /  |_> > __ \|  | \/\___ \\  ___/|  | \/
   |__|_|  /____/|   __(____  /__|  /____  >\___  >__|
         \/      |__|       \/           \/     \/
   Copyright (C) 2004 - 2022 Ingo Berg

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

// Small example using the cmake imported target. Include file and link library
// should work automagically.

#include <muParser.h>
#include <muParserDef.h>

int main()
{
   mu::Parser  parser;

   mu::value_type  values[] = { 1, 2 };
   parser.DefineVar("a", &values[0]);
   parser.DefineVar("b", &values[1]);

   std::string expr = "a + b";
   parser.SetExpr("a + b");
   mu::value_type ans = parser.Eval();
   std::cout << expr << " == " << ans << "\n";

   return (ans == 3.0) ? 0 : -1;
}
