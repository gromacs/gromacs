int main()
{
/* This detects 2.8, 2.9 and 3.0 versions for both C and C++ clang. It
 * detects the version of the LLVM back end, and not (for example) the
 * Apple clang version number (which might be 4.1 or some number based
 * on its "compatibility with gcc 4.2.1," even though the LLVM back
 * end is 3.0!)
 */
#if ((__clang_major__ == 3) && (__clang_minor__ == 0)) ||
    ((__clang_major__ == 2) && ((__clang_minor__ == 8) ||
                                (__clang_minor__ == 9)))
    return 0;
#else
    trigger syntax error
#endif
}
