

#ifdef HAVE_TMPI_CONFIG_H
#include "tmpi_config.h"
#endif


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef THREAD_WINDOWS
    #ifdef __MINGW32__
       #define _WIN32_WINNT 0x0601 /* Windows 7*/
    #endif
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

#ifndef THREAD_WINDOWS

/* We don't have specific NUMA aware allocators: */

void *tMPI_Malloc_local(size_t size)
{
    return malloc(size);
}

void *tMPI_Calloc_local(size_t nmemb, size_t size)
{
    return calloc(nmemb, size);
}

void *tMPI_Realloc_local(void *ptr, size_t size)
{
    return realloc(ptr, size);
}

int tMPI_Free_numa(void *ptr)
{
    free(ptr);
    return 0; /* we don't detect errors here */
}

#else

#define TMPI_NUMA_MALLOC

/*
    Windows NUMA memory allocation support.

    NUMA support is implemented to maximize the chance that memory access
    patterns remain Local to the NUMA node.  This avoids penalties accessing
    "remote" memory.
    An important assumption here is that code paths which allocate and
    reallocate heap blocks are likely to be accessing that allocated memory
    on the same NUMA node.
    Testing has shown the above criteria to be met, yielding gains of > 15%
    when on Windows with NUMA hardware.

    The high level approach is:
    1. Use a separate heap per NUMA node.  This reduces heap contention, steers
       allocations to the local NUMA node, and avoids re-use of freed heap
       blocks across (remote) NUMA nodes.
    2. Allocate each heap locally to each NUMA node, such that heap control
       structures are on the NUMA node accessing the heap.
    3. During realloc operations, transfer the new block to the local NUMA
       node, if appropriate.  This is a rare case when thread affinity and
       access patterns are correct.
    4. Use GetProcAddress() to obtain function pointers to functions that are
       operating system version dependent, to allow maximum binary
       compatibility.

    Scott Field (sfield@microsoft.com)      Jan-2011
 */

#include <windows.h>

/*
    __declspec(align()) may not be supported by all compilers, so define the
    size of the structure manually to force alignment
    note that HeapAlloc() already returns aligned blocks.
    typedef __declspec(align(32)) struct ...
 */
typedef struct {
    DWORD  dwMagic;
    HANDLE hHeap;               /* 8 */
    size_t cbAllocationSize;    /* 8 */
    ULONG  ProcessorNumber;     /* processor number at time of allocation
                                   (development/testing) */
    USHORT NodeNumber;          /* NUMA node number at time of allocation
                                   (development/testing) */
} HEAPHEADER, *PHEAPHEADER;

#define HEAP_NUMA_MAGIC     0x05AF0777
#define HEAPHEADER_SIZE     (32)

#ifdef C_ASSERT
/* fail compile if size of HEAPHEADER exceeds pre-defined value */
C_ASSERT(sizeof(HEAPHEADER) <= HEAPHEADER_SIZE);
#endif

/* function prototypes and variables to support obtaining function
   addresses dynamically -- supports down-level operating systems */

typedef BOOL (WINAPI *func_GetNumaHighestNodeNumber_t)( PULONG HighestNodeNumber );
typedef BOOL (WINAPI *func_GetNumaProcessorNodeEx_t)( PPROCESSOR_NUMBER Processor, PUSHORT NodeNumber );
typedef VOID (WINAPI *func_GetCurrentProcessorNumberEx_t)( PPROCESSOR_NUMBER ProcNumber );

func_GetNumaHighestNodeNumber_t             smalloc_GetNumaHighestNodeNumber;         /* WinXP SP2, WinXP64, WinSrv 2003 */
func_GetNumaProcessorNodeEx_t               smalloc_GetNumaProcessorNodeEx;           /* Windows 7, WinSrv 2008R2 */
func_GetCurrentProcessorNumberEx_t          smalloc_GetCurrentProcessorNumberEx;      /* Windows 7, WinSrv 2008R2 */

#define NUMA_STATUS_UNKNOWN     (0)
#define NUMA_STATUS_NOT_NUMA    (1)
#define NUMA_STATUS_NUMA        (2)

DWORD   g_dwTlsHeap;    /* TLS slot used for preferred heap handle */
HANDLE *g_hHeap;        /* array of heap handles */
ULONG   g_ulNumaStatus; /* 0 = unknown, 1 = not NUMA, 2 = NUMA */

VOID
InitNumaHeapSupport(
        VOID
        )
{
    HMODULE hModKernel32;   /* module handle to kernel32.dll -- we already
                               reference it, so it's already loaded */
    ULONG   ulNumaHighestNodeNumber;

    /* grab the addresses for the NUMA functions.
       It's fine if there is a race condition reaching this routine */
    hModKernel32 = GetModuleHandleA("kernel32.dll");
    if (hModKernel32 == NULL)
    {
        g_ulNumaStatus = NUMA_STATUS_NOT_NUMA;
        return;
    }

    smalloc_GetNumaHighestNodeNumber    = (func_GetNumaHighestNodeNumber_t)GetProcAddress( hModKernel32, "GetNumaHighestNodeNumber" );
    smalloc_GetCurrentProcessorNumberEx = (func_GetCurrentProcessorNumberEx_t)GetProcAddress( hModKernel32, "GetCurrentProcessorNumberEx" );
    smalloc_GetNumaProcessorNodeEx      = (func_GetNumaProcessorNodeEx_t)GetProcAddress( hModKernel32, "GetNumaProcessorNodeEx" );

    if ( (smalloc_GetNumaHighestNodeNumber == NULL) ||
         (smalloc_GetCurrentProcessorNumberEx == NULL) ||
         (smalloc_GetNumaProcessorNodeEx == NULL) )
    {
        g_ulNumaStatus = NUMA_STATUS_NOT_NUMA;
        return;
    }

    /* determine how many NUMA nodes are present */

    if (!smalloc_GetNumaHighestNodeNumber(&ulNumaHighestNodeNumber) ||
        (ulNumaHighestNodeNumber == 0) )
    {
        g_ulNumaStatus = NUMA_STATUS_NOT_NUMA;
        return;
    }

    /* handle deferred creation of TLS slot.
       note: this could be moved to one-time init path.
       failures here result in assuming the system is not NUMA.
     */

    if (g_dwTlsHeap == 0)
    {
        DWORD dwTlsHeap = TlsAlloc();
        DWORD dwPriorValue;

        if (dwTlsHeap == TLS_OUT_OF_INDEXES)
        {
            g_ulNumaStatus = NUMA_STATUS_NOT_NUMA;
            return;
        }

        dwPriorValue = (DWORD)InterlockedCompareExchange(
                    (LONG volatile *)&g_dwTlsHeap,
                    (LONG) dwTlsHeap,
                    0
                    );

        if (dwPriorValue != 0)
        {
            TlsFree( dwTlsHeap );
        }
    }

    /* handle deferred creation of heap handle array.
       note: this could be moved to one-time init path.
     */

    if (g_hHeap == NULL)
    {
        HANDLE *hHeapNew;
        HANDLE *hPriorValue;

        /* allocate an array to contain a heap handle for each NUMA node */
        hHeapNew = (HANDLE*)HeapAlloc(
                    GetProcessHeap(),
                    HEAP_ZERO_MEMORY,
                    sizeof(HANDLE) * (ulNumaHighestNodeNumber+1)
                    );

        if (hHeapNew == NULL)
        {
            g_ulNumaStatus = NUMA_STATUS_NOT_NUMA;
            return;
        }

#if defined(WIN64) || defined( _WIN64 )
        hPriorValue = (HANDLE *)InterlockedCompareExchange64(
                    (LONGLONG volatile *)&g_hHeap,
                    (LONGLONG) hHeapNew,
                    0
                    );
#else
        hPriorValue = (HANDLE *)InterlockedCompareExchange(
                    (LONG volatile *)&g_hHeap,
                    (LONG) hHeapNew,
                    0
                    );
#endif

        if (hPriorValue != NULL)
        {
            HeapFree(GetProcessHeap(), 0, hHeapNew);
        }
    }

    /* indicate system is NUMA */
    g_ulNumaStatus = NUMA_STATUS_NUMA;

    return;
}

HANDLE
ReturnHeapHandle(
        VOID
        )
{
    HANDLE           hHeap;                     /* preferred heap handle to
                                                   return to caller */
    PROCESSOR_NUMBER CurrentProcessorNumber;    /* processor number associated
                                                   with calling thread */
    USHORT           CurrentNumaNodeNumber;     /* NUMA node number assocaited
                                                   with calling thread */

    /* determine NUMA status of system. */

    if (g_ulNumaStatus == NUMA_STATUS_UNKNOWN)
    {
        InitNumaHeapSupport();
        if (g_ulNumaStatus == NUMA_STATUS_NOT_NUMA)
        {
            return GetProcessHeap();
        }
    }
    else if (g_ulNumaStatus == NUMA_STATUS_NOT_NUMA)
    {
        /* not NUMA, return the process heap handle */
        return GetProcessHeap();
    }


    /* return the preferred heap handle from the TLS slot, if set.
       This is the commonly taken path. */

    hHeap = (HANDLE)TlsGetValue( g_dwTlsHeap );

    if (hHeap != NULL)
    {
        return hHeap;
    }


    /* preferred heap handle not yet set.
       determine the numa node we're executing on, and create a heap which
       is assigned to this node.
       one (soft) assumption that is made here is that thread affinity has
       been set such that threads do not move between NUMA nodes.
     */

    smalloc_GetCurrentProcessorNumberEx(&CurrentProcessorNumber);

    if (!smalloc_GetNumaProcessorNodeEx(&CurrentProcessorNumber, &CurrentNumaNodeNumber))
    {
        /* GetNumaProcessorNodeEx() can fail on WOW64/32bit if invoked
           against processor numbers > 32.
           this should never be reached for 64bit builds.
         */
        CurrentNumaNodeNumber = 0;
    }


    /* check if the NUMA node array slot already contains a heap */
    /* CurrentNumaNodeNumber cannot execeed count of heaps, as NUMA nodes
       cannot be added. */

    hHeap = g_hHeap[ CurrentNumaNodeNumber ];

    if (hHeap == NULL)
    {
        HANDLE hHeapPrior = NULL;
        ULONG  ulOption   = 2; /* HEAP_LFH */

        /* create a heap for this numa node
           defer creating the heap - while running on each node - to ensure
           the heap control structures get created on the local NUMA node.
         */

        hHeap = HeapCreate(0, 0, 0);

        if (hHeap == NULL)
        {
            /* just return the process heap.  We'll try to create a heap
               again next time */
            return GetProcessHeap();
        }

        /* make the new heap a low-fragmentation heap */

        HeapSetInformation(
                hHeap,
                0,          /* HeapCompatibilityInformation */
                &ulOption,
                sizeof(ulOption)
                );

        /* set the array slot entry for this NUMA node to contain the newly
           allocated heap */

        hHeapPrior = (HANDLE)InterlockedCompareExchangePointer(&(g_hHeap[CurrentNumaNodeNumber]), hHeap, NULL);
        if (hHeapPrior != NULL)
        {
            HeapDestroy( hHeap );
            hHeap = hHeapPrior;
        }
    }

    /* we reached here since there was no heap assigned to the TLS slot.
       Assign it. */
    TlsSetValue(g_dwTlsHeap, hHeap);

    return hHeap;
}




void *tMPI_Malloc_local(size_t size)
{
    HANDLE         hHeap;
    unsigned char *ptr;
    size_t         new_size;
    HEAPHEADER    *phdr;

    hHeap = ReturnHeapHandle();

    new_size = size + HEAPHEADER_SIZE;

    ptr = (unsigned char *)HeapAlloc( hHeap, 0, new_size );

    if (ptr == NULL)
    {
        return NULL;
    }

    phdr = (HEAPHEADER*)ptr;

    phdr->dwMagic          = HEAP_NUMA_MAGIC;
    phdr->hHeap            = hHeap;     /* track the heap handle for realloc
                                           and free */
    phdr->cbAllocationSize = new_size;  /* track the allocation size for
                                           realloc and debugging */

    return ( ptr + HEAPHEADER_SIZE );
}


void *tMPI_Calloc_local(size_t nelem, size_t elsize)
{
    void  *ptr;
    size_t size = nelem * elsize;

    ptr = tMPI_Malloc_local(size);

    if (ptr != NULL)
    {
        memset(ptr, 0, size);
    }

    return ptr;
}


void *tMPI_Realloc_local(void *ptr, size_t size)
{
    HANDLE         hHeap;
    HEAPHEADER    *phdr;
    unsigned char *new_ptr;
    size_t         new_size;

    /* calculate the allocation address and check for presence of the hint
       which indicates this was allocated by our allocator.
     */

    phdr = (HEAPHEADER*)((unsigned char*)ptr - HEAPHEADER_SIZE);

    if (phdr->dwMagic != HEAP_NUMA_MAGIC)
    {
        /* TODO: call tMPI_Error() */
        /*gmx_fatal(errno,__FILE__,__LINE__,
           "Invalid Heap header during realloc. %p", ptr);*/

        return NULL;
    }

    /* calculate size of new/realloc'd block.
     */

    new_size = size + HEAPHEADER_SIZE;

    /* if the NUMA Node changed between the initial allocation and the
       reallocation, copy the memory to an allocation on the new local node:
       we assume the new realloc'd block is more likely to be manipulated by
       the current thread which is calling realloc.
       the simple way to detect this condition is to compare the preferred heap
       handle with the heap handle stored in the current memory block.
     */

    hHeap = ReturnHeapHandle();

    if (hHeap != phdr->hHeap)
    {
        new_ptr = HeapAlloc( hHeap, 0, new_size );

        /* if the new allocation succeeded, copy the buffer and free the
           original buffer.
         */

        if (new_ptr != NULL)
        {
            size_t copy_size;

            /* the realloc can be larger or smaller than the original
               allocation size.
             */

            if (new_size > phdr->cbAllocationSize)
            {
                copy_size = phdr->cbAllocationSize;
            }
            else
            {
                copy_size = new_size;
            }

            /* copy the current memory block contents into the newly allocated
               buffer, and then free the original buffer.
             */

            memcpy( new_ptr, phdr, copy_size );

            HeapFree( phdr->hHeap, 0, phdr );
        }

    }
    else
    {

        /* NodeNumber of existing allocation matches current node number.
           realloc from the heap associated with the existing allocation.
         */

        hHeap = phdr->hHeap;

        new_ptr = HeapReAlloc(
                    hHeap,
                    0,
                    phdr,
                    new_size
                    );
    }

    if (new_ptr == NULL)
    {
        return NULL;
    }

    phdr                   = (HEAPHEADER*)new_ptr;
    phdr->cbAllocationSize = new_size;  /* update allocation size to match
                                           realloc */
    phdr->hHeap            = hHeap;

    return ( new_ptr + HEAPHEADER_SIZE );
}


int tMPI_Free_numa(void *ptr)
{
    HEAPHEADER *phdr;

    /* caller doesn't call us on ptr == NULL case, so we don't need to
       check that here. */

    phdr = (HEAPHEADER*)((unsigned char*)ptr - HEAPHEADER_SIZE);

    /* this check should happen in __try / __except block until the
       mis-matched free/sfree calls are fixed, but this is here primarilly
       for debugging and catching mis-matched memory alloc and free references.
     */

    if (phdr->dwMagic != HEAP_NUMA_MAGIC)
    {
        /* ptr is leaked here, rather than faulting in the allocator.
           this is done in order to track mis-matched alloc/free calls.
         */
        return 1;
    }

    phdr->dwMagic = 0;

    HeapFree( phdr->hHeap, 0, phdr );

    return 0;
}

#endif /* NUMA allocation functions for (_WIN32 || _WIN64) */
