#!/usr/bin/env bash
#
# Regenerate the committed flex/bison parser sources under generated/ from the
# grammar sources in modules/cpp-cif-parser/src/. Run this ONLY when bumping the
# vendored RCSB mmCIF suite to a new upstream version (see UPDATING.md). A normal
# GROMACS build never invokes flex or bison — it compiles the committed output.
#
# The flex/bison flags mirror upstream cmake/commonparser.cmake so the generated
# code keeps the cifparser_/dicparser_ symbol prefixes the base classes expect.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src="$here/modules/cpp-cif-parser/src"
gen="$here/generated"

command -v flex  >/dev/null || { echo "error: flex not found";  exit 1; }
command -v bison >/dev/null || { echo "error: bison not found"; exit 1; }

mkdir -p "$gen"

echo "Regenerating CIF parser/scanner..."
bison -d -v -l -p cifparser_ -o "$gen/CifParser.c"  "$src/CifParser.y"
flex  -Cfr -L -Pcifparser_    -o "$gen/CifScanner.c" "$src/CifScanner.l"

echo "Regenerating DIC parser/scanner..."
bison -d -v -l -p dicparser_ -o "$gen/DICParser.c"  "$src/DICParser.y"
flex  -Cfr -L -Pdicparser_    -o "$gen/DICScanner.c" "$src/DICScanner.l"

# bison writes CifParser.h / DICParser.h next to the .c with -d; drop the .output
rm -f "$gen"/*.output
echo "Done. Review 'git diff generated/' before committing."
