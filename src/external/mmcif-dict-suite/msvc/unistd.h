/*
 * MSVC has no <unistd.h>. The bundled mmCIF suite includes it only for open(),
 * read(), write(), lseek(), close() and isatty(), which the Windows C runtime
 * declares in <io.h>, and for the permission bits Serializer.C passes to open().
 *
 * Without O_BINARY, Serializer opens its binary database files in text mode on
 * Windows. GROMACS reads and writes only text mmCIF files, never those.
 */
#ifndef GMX_MMCIF_MSVC_UNISTD_H
#define GMX_MMCIF_MSVC_UNISTD_H

#include <io.h>
#include <sys/stat.h>

#ifndef S_IRUSR
#    define S_IRUSR _S_IREAD
#endif
#ifndef S_IWUSR
#    define S_IWUSR _S_IWRITE
#endif

#endif
