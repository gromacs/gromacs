# GROMACS-owned build for a trimmed subset of the RCSB mmCIF Dictionary Suite.
#
# Upstream: https://github.com/rcsb/cpp-dict-pack  (mmcif-dict-suite v2.500).
# Upstream ships a heavyweight build (flex/bison code generation plus network
# download and compilation of the mmCIF dictionaries). GROMACS only needs the
# CIF/PDBx *file* parser and writer, so this file builds a single static library
# `mmcif` from just the parser/table/cif modules, using pre-generated flex/bison
# parser sources committed under generated/. As a result the GROMACS build needs
# no new tools (no flex, bison, python) and performs no network access.
#
# See UPDATING.md for how to refresh the vendored copy / regenerate the parsers.

set(MMCIF_MODULES "${CMAKE_CURRENT_SOURCE_DIR}/modules")
set(MMCIF_GENERATED "${CMAKE_CURRENT_SOURCE_DIR}/generated")

# Explicit source lists, mirroring the object modules that upstream combines
# into `mmciflib-all` (cmake/commonparser.cmake), minus the dict-object-file
# module and command-line tools, which the CIF read/write path does not use.
set(MMCIF_SOURCES
    # cpp-common
    ${MMCIF_MODULES}/cpp-common/src/RcsbPlatform.C
    ${MMCIF_MODULES}/cpp-common/src/RcsbFile.C
    ${MMCIF_MODULES}/cpp-common/src/BlockIO.C
    ${MMCIF_MODULES}/cpp-common/src/CifString.C
    ${MMCIF_MODULES}/cpp-common/src/Serializer.C
    ${MMCIF_MODULES}/cpp-common/src/GenString.C
    ${MMCIF_MODULES}/cpp-common/src/GenCont.C
    ${MMCIF_MODULES}/cpp-common/src/Exceptions.C
    ${MMCIF_MODULES}/cpp-common/src/DataInfo.C
    ${MMCIF_MODULES}/cpp-common/src/mapped_vector.C
    ${MMCIF_MODULES}/cpp-common/src/mapped_ptr_vector.C
    # cpp-tables
    ${MMCIF_MODULES}/cpp-tables/src/ISTable.C
    ${MMCIF_MODULES}/cpp-tables/src/ITTable.C
    ${MMCIF_MODULES}/cpp-tables/src/TTable.C
    ${MMCIF_MODULES}/cpp-tables/src/TableFile.C
    # cc-regex (C sources; engine.c and split.c are #included, not compiled directly)
    ${MMCIF_MODULES}/cc-regex/src/regcomp.c
    ${MMCIF_MODULES}/cc-regex/src/regexec.c
    ${MMCIF_MODULES}/cc-regex/src/regerror.c
    ${MMCIF_MODULES}/cc-regex/src/regfree.c
    # cpp-cif-file
    ${MMCIF_MODULES}/cpp-cif-file/src/CifFile.C
    ${MMCIF_MODULES}/cpp-cif-file/src/DicFile.C
    ${MMCIF_MODULES}/cpp-cif-file/src/ParentChild.C
    ${MMCIF_MODULES}/cpp-cif-file/src/CifParentChild.C
    ${MMCIF_MODULES}/cpp-cif-file/src/CifDataInfo.C
    ${MMCIF_MODULES}/cpp-cif-file/src/CifExcept.C
    # cpp-cif-file-util
    ${MMCIF_MODULES}/cpp-cif-file-util/src/CifFileUtil.C
    ${MMCIF_MODULES}/cpp-cif-file-util/src/CifCorrector.C
    # cpp-cif-parser (hand-written base classes)
    ${MMCIF_MODULES}/cpp-cif-parser/src/CifFileReadDef.C
    ${MMCIF_MODULES}/cpp-cif-parser/src/CifScannerBase.C
    ${MMCIF_MODULES}/cpp-cif-parser/src/CifParserBase.C
    ${MMCIF_MODULES}/cpp-cif-parser/src/DICScannerBase.C
    ${MMCIF_MODULES}/cpp-cif-parser/src/DICParserBase.C
    # pre-generated flex/bison parser sources (compiled as C)
    ${MMCIF_GENERATED}/CifParser.c
    ${MMCIF_GENERATED}/CifScanner.c
    ${MMCIF_GENERATED}/DICParser.c
    ${MMCIF_GENERATED}/DICScanner.c
    )

add_library(mmcif STATIC ${MMCIF_SOURCES})

# Headers reference each other by flat name across modules, so every module
# include directory must be on the search path; generated/ holds the bison
# headers (CifParser.h, DICParser.h). Exposed PUBLIC + SYSTEM so the single
# consumer (fileio/cifio.cpp) can include CifFile.h etc. without warnings.
target_include_directories(mmcif SYSTEM PUBLIC
    ${MMCIF_MODULES}/cpp-common/include
    # ISTable.h / TableFile.h #include the mapped_*.C template bodies from here
    ${MMCIF_MODULES}/cpp-common/src
    ${MMCIF_MODULES}/cpp-tables/include
    ${MMCIF_MODULES}/cc-regex/include
    ${MMCIF_MODULES}/cpp-cif-file/include
    ${MMCIF_MODULES}/cpp-cif-file-util/include
    ${MMCIF_MODULES}/cpp-cif-parser/include
    ${MMCIF_GENERATED}
    )

# This is third-party code: keep it position-independent for a shared libgromacs
# and silence its (many) warnings so it does not pollute the GROMACS build.
set_target_properties(mmcif PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|Intel|IntelLLVM")
    target_compile_options(mmcif PRIVATE -w)
elseif(MSVC)
    target_compile_options(mmcif PRIVATE /w)
endif()

# ---------------------------------------------------------------------------
# WORKAROUNDS for defects in the vendored RCSB sources.
#
# The two settings below work around problems that live in the upstream RCSB
# code, not in GROMACS. They are TEMPORARY: the intent is to submit fixes to
# the upstream project (https://github.com/rcsb/cpp-dict-pack); if/when those
# are accepted and this vendored copy is updated, remove these workarounds.
# See UPDATING.md ("Local workarounds") for details.
#
# 1. GenString.C/.h and GenCont.C use std::bind2nd / std::not1, deprecated in
#    C++11 and removed from the standard library in C++17. libstdc++ still
#    offers them as an extension, but libc++ does not once the TU is compiled
#    as C++17+, so a C++20 GROMACS build fails to compile this library. Build
#    just this vendored target as C++14, where the binders still exist
#    everywhere; it is a standalone static library, so its language standard is
#    independent of the rest of GROMACS.
set_target_properties(mmcif PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED ON)
#
# 2. The Clang static analyzer (GMX_CLANG_ANALYZER) reports harmless dead stores
#    in this third-party code, which fails the analyzer CI job. Those stores are
#    annotated in the sources with [[clang::suppress]] (grep for it under
#    modules/); a compiler-flag bypass is not possible because the analyzer
#    build compiles with g++ -- which rejects Clang-only analyzer flags -- while
#    only the analysis pass is Clang, and both share the compile command.
# ---------------------------------------------------------------------------

# Keep the target out of IDE clutter alongside the other bundled externals.
set_target_properties(mmcif PROPERTIES FOLDER "External")
