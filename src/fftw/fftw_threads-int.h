/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef FFTW_THREADS_INT_H
#define FFTW_THREADS_INT_H

#include <fftw-int.h>
#include <fftw_threads-int.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/************************* Thread Glue *************************/

/* Adding support for a new shared memory thread API should be easy.  You
   simply do the following things (look at the POSIX and Solaris
   threads code for examples):

   * Invent a symbol of the form FFTW_USING_FOO_THREADS to denote
     the use of your thread API, and add an 
              #elif defined(FFTW_USING_FOO_THREADS)
      before the #else clause below.  This is where you will put
      your thread definitions.  In this #elif, insert the following:

      -- #include any header files needed to use the thread API.

      -- Typedef fftw_thread_function to be a function pointer
         of the type used as a argument to fftw_thread_spawn
	 (i.e. the entry function for a thread).

      -- Define fftw_thread_id, via a typedef, to be the type
         that is used for thread identifiers.

      -- #define fftw_thread_spawn(tid_ptr, proc, data) to
         call whatever function to spawn a new thread.  The
         new thread should call proc(data) as its starting point,
         and tid_ptr is a pointer to a fftw_thread_id that
         is set to an identifier for the thread.  You can also
         define this as a subroutine (put it in fftw_threads.c)
	 if it is too complicated for a macro.  The prototype should
	 be:

	 void fftw_thread_spawn(fftw_thread_id *tid_ptr,
	                        fftw_thread_function proc,
				void *data);

      -- #define fftw_thread_wait(tid) to block until the thread
         whose identifier is tid has terminated.  You can also
         define this as a subroutine (put it in fftw_threads.c) if
	 it is too complicated for a macro.  The prototype should be:
	
	 void fftw_thread_wait(fftw_thread_id tid);

   * If you need to perform any initialization before using threads,
     put your initialization code in the fftw_threads_init() function
     in fftw_threads.c, bracketed by the appropriate #ifdef of course.

   * Also, of course, you should modify fftw/config.h to #define
     FFTW_USING_FOO_THREADS, or better yet modify config.h.in
     and configure.in so that autoconf can automatically detect
     your threads library.

   * Finally, if you do implement support for a new threads API, be
     sure to let us know at fftw@theory.lcs.mit.edu so that we can
     distribute your code to others!

*/

/************************** Solaris Threads ****************************/
      
#if defined(FFTW_USING_SOLARIS_THREADS)

/* Solaris threads glue.  Tested. */

/* link with -lthread */

#include <thread.h>

/* Thread entry point: */
typedef void * (*fftw_thread_function) (void *);

typedef thread_t fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) \
     thr_create(0,0,proc,data,THR_BOUND,tid_ptr)

#define fftw_thread_wait(tid) thr_join(tid,0,0)

/************************** BeOS Threads ****************************/
      
#elif defined(FFTW_USING_BEOS_THREADS)

/* BeOS threads glue.  Tested for DR8.2. */

#include <OS.h>

/* Thread entry point: */
typedef thread_entry fftw_thread_function;

typedef thread_id fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) { \
     *(tid_ptr) = spawn_thread(proc,"FFTW",B_NORMAL_PRIORITY,data); \
     resume_thread(*(tid_ptr)); \
}

/* wait_for_thread requires that we pass a valid pointer as the
   second argument, even if we're not interested in the result. */
#define fftw_thread_wait(tid) {long exit_val;wait_for_thread(tid, &exit_val);}

/************************** MacOS Threads ****************************/

#elif defined(FFTW_USING_MACOS_THREADS)

/* MacOS MP threads glue. Experimental, untested! I do not have an
   MP MacOS system available to me...I just read the documentation.
   There is actually a good chance that this will work (since the
   code below is so short), but I make no guarantees.  Consider
   it to be a starting point for your own implementation. 

   I also had to insert some code in fftw_threads.c. 

   MacOS X will have real SMP support, thank goodness. */

/* Using this code in the MacOS: (See the README file for general
   documenation on the FFTW threads code.)  To use this code, you have
   to do two things.  First of all, you have to #define the symbol
   FFTW_USING_MACOS_THREADS. This can be done at the top of this file
   or perhaps in your compiler options.  Second, you have to weak-link
   your project to the MP library.

   In your code, you should check at run-time with MPLibraryIsLoaded()
   to see if the MP library is available. If it is not, it is
   still safe to call the fftw_threads routines...in this case,
   however, you must always pass 1 for the nthreads parameter!
   (Otherwise, you will probably want to pass the value of
   MPProcessors() for the nthreads parameter.) */

#include <MP.h>

typedef TaskProc fftw_thread_function;

typedef MPQueueID fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) { \
     MPTaskID task; \
     MPCreateQueue(tid_ptr); \
     MPCreateTask(proc,data,kMPUseDefaultStackSize,*(tid_ptr),0,0, \
		  kMPNormalTaskOptions,&task); \
}

#define fftw_thread_wait(tid) { \
     void *param1,*param2,*param3; \
     MPWaitOnQueue(tid,&param1,&param2,&param3,kDurationForever); \
     MPDeleteQueue(tid); \
}

/************************** Win32 Threads ****************************/

#elif defined(FFTW_USING_WIN32_THREADS)

/* Win32 threads glue.  We have not tested this code!  (I just implemented
   it by looking at a Win32 threads manual.)  Users have reported that this
   code works under NT using Microsoft compilers.
   
   To use it, you should #define the symbol FFTW_USING_WIN32_THREADS. */

#include <windows.h>

typedef LPTHREAD_START_ROUTINE fftw_thread_function;
typedef HANDLE fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) { \
     DWORD thrid; \
     *(tid_ptr) = CreateThread((LPSECURITY_ATTRIBUTES) NULL, 0, \
			       (fftw_thread_function) proc, (LPVOID) data, \
			       0, &thrid); \
}

#define fftw_thread_wait(tid) { \
     WaitForSingleObject(tid, INFINITE); \
     CloseHandle(tid); \
}

/************************** Mach cthreads ****************************/

#elif defined(FFTW_USING_MACH_THREADS)

#ifdef HAVE_MACH_CTHREADS_H
#include <mach/cthreads.h>
#elif defined(HAVE_CTHREADS_H)
#include <cthreads.h>
#elif defined(HAVE_CTHREAD_H)
#include <cthread.h>
#endif

typedef cthread_fn_t fftw_thread_function;

typedef cthread_t fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) \
     *(tid_ptr) = cthread_fork(proc, (any_t) (data))

#define fftw_thread_wait(tid) cthread_join(tid)

/************************** POSIX Threads ****************************/

#else /* use the default, POSIX threads: */

/* POSIX threads glue.  Tested. */

#ifndef FFTW_USING_POSIX_THREADS
#define FFTW_USING_POSIX_THREADS
#endif

/* link with -lpthread */

#include <pthread.h>

/* Thread entry point: */
typedef void * (*fftw_thread_function) (void *);

extern pthread_attr_t *fftw_pthread_attributes_p;

typedef pthread_t fftw_thread_id;

#define fftw_thread_spawn(tid_ptr, proc, data) \
     pthread_create(tid_ptr,fftw_pthread_attributes_p,proc,data)

#define fftw_thread_wait(tid) pthread_join(tid,0)

#endif

/************************ Function prototypes ***********************/

#ifndef __GNUC__
#     if HAVE_ALLOCA_H
#          include <alloca.h>
#     else
#          ifndef alloca /* predefined by HP cc +Olibcalls */
extern void *alloca(size_t n);
#          endif
#     endif
#endif

/* Use alloca instead of malloc, if possible, in the hope that alloca
   will be faster at allocating the small blocks we need.  (In principle,
   alloca just needs to increment the stack pointer.) */
#ifdef HAVE_ALLOCA
#     define ALLOCA(n) alloca(n)
#     define ALLOCA_CLEANUP(p) ;
#else
#     define ALLOCA(n) fftw_malloc(n)
#     define ALLOCA_CLEANUP(p) fftw_free(p);
#endif

typedef struct {
     int min, max, thread_num;
     void *data;
} fftw_loop_data;

typedef void *(*fftw_loop_function) (fftw_loop_data *);

extern void fftw_thread_spawn_loop(int loopmax, int nthreads,
				   fftw_loop_function proc, void *data);

extern void fftw_executor_many_inplace_threads(int n, fftw_complex *in,
					       fftw_complex *work,
					       fftw_plan_node *p,
					       int istride,
					       int howmany, int idist,
					       int nthreads);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_THREADS_INT_H */
