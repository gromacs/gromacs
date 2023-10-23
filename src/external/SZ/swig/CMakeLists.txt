#[[
A set of python bindings for SZ 

Developed by Robert Underwood while he was at Clemson University
This material is based upon work supported by the National Science 
Foundation under Grant No. 1633608.

Copyright © 2019 Robert Underwood
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Robert Underwood ''AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Robert Underwood BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation
are those of the authors and should not be interpreted as representing
official policies, either expressed or implied, of Robert Underwood.
]]

cmake_minimum_required(VERSION 3.14)
find_package(SWIG)
if(SWIG_FOUND AND BUILD_SHARED_LIBS)
  include(UseSWIG)


  find_package(Python COMPONENTS Interpreter Development NumPy)
  if(NOT SZ_PYTHON_SITELIB)
    set(SZ_PYTHON_SITELIB ${Python_SITELIB} CACHE PATH "path to install python libraries to")
  endif()

  set_property(SOURCE pysz.i PROPERTY CPLUSPLUS ON)
  set_property(SOURCE pysz.i PROPERTY USE_TARGET_INCLUDE_DIRECTORIES ON)

  swig_add_library(
    pysz
    LANGUAGE python
    SOURCES pysz.i pysz.cc
    )

  target_link_libraries(pysz PUBLIC SZ Python::Python Python::NumPy)
  target_include_directories(pysz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  get_property(swig_generated_module TARGET pysz PROPERTY SWIG_SUPPORT_FILES)
  install(TARGETS pysz DESTINATION ${SZ_PYTHON_SITELIB})
  install(FILES ${swig_generated_module} DESTINATION ${SZ_PYTHON_SITELIB})

endif()
