#!/usr/bin/perl -w
#
# Usage example:
#
#  perl -0777 -i admin/refactoring/remove-cpp-guards.pl $(git grep -l "extern \"C" -- src/gromacs)
#
BEGIN { $/ = undef; $\ = undef; }
LINE: while (defined($_ = readline ARGV)) {
    s/\#ifdef __cplusplus\nextern "C"\s\{\n\#endif\n\n?//gm;
    s/\#ifdef 0\n\}\n\#endif\n\n?//gm;
    s/\#ifdef __cplusplus\n\}.*\n\#endif\n\n?//gm;
    s/\#if GMX_CXX11_COMPILATION\n(    static_assert.*\n?.*;)\n#endif\n\n?/$1\n/gm;
}
continue {
    die "-p destination: $!\n" unless print $_;
}
