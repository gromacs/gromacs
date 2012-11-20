int main()
{
/* This works for both C and C++ clang, and detects the
 * LLVM back end, not the Apple clang version number (which
 * might be 4.2.1 for LLVM back end 3.0!)
 */
#if (__clang_major__ == 3) && (__clang_minor__ == 0)
    return 0;
#else
    trigger syntax error
#endif
}
