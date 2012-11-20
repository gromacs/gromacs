int main()
{
/* This detects 3.0 versions for both C and C++ clang. It detects the
 * version of the LLVM back end, and not (for example) the Apple clang
 * version number (which might be 4.1 or some number based on its
 * "compatibility with gcc 4.2.1," even though the LLVM back end is
 * 3.0!).
 *
 * If/when we have time or user complaints, we can maybe ban earlier
 * versions of clang, but we don't actually know there's a problem
 * with them at the time of this commit.
 */
#if (__clang_major__ == 3) && (__clang_minor__ == 0)
    return 0;
#else
#error clang version information not found
#endif
}
