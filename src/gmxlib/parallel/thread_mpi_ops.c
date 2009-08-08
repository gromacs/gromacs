
#ifdef THREAD_MPI_OPS

/* cpp wizardry follows... 

This file is #included directly from thread_mpi.c, and constructs 
MPI_Reduce operators. 

What this does is create the min, max, sum, prod, etc. functions for a given
datatype (pre-defined as TYPE, with identifier name TYPENM) and puts pointers 
to these functions in an array called oplist_TYPENM. 

gmx_thread_mpi_reduce.c includes this file once for each type used by MPI, 
and thus builds up a set of arrays of function pointers, that then get used
in the mpi_datatype_ structure. This way, each operation/datatype entry
that makes sense can be extracted easily. Note that we don't (yet) support
user-defined ops */

#define FNAMEr(tp,fn) tMPI_##tp##_##fn
#define FNAME(tp,fn) FNAMEr(tp,fn)

/* macros to define functions and prototypes based on a name and an operation */
#define FNr(tp,fname,fn) \
static void tMPI_##tp##_##fname  (void *dest, void *src_a, void *src_b, \
                                  int count) \
{ \
    /*printf("in function %s, count=%d\n", __FUNCTION__, count);*/\
    TYPE *a=(TYPE*)src_a; \
    TYPE *b=(TYPE*)src_b; \
    TYPE *d=(TYPE*)dest; \
    int i; \
    for(i=0;i<count;i++) \
        d[i]=(TYPE)(fn(a[i],b[i])); \
}  

#define FN(tp,fname,fn) FNr(tp,fname,fn)

#define OPFNr(tp,fname,operator)  \
static void tMPI_##tp##_##fname  (void *dest, void *src_a, void *src_b, \
                                  int count) \
{ \
    /*printf("in function %s, count=%d\n", __FUNCTION__, count);*/\
    TYPE *a=(TYPE*)src_a; \
    TYPE *b=(TYPE*)src_b; \
    TYPE *d=(TYPE*)dest; \
    int i; \
    for(i=0;i<count;i++) \
        d[i]=(TYPE)(a[i] operator b[i]); \
}  

#define OPFN(tp,fname,operator) OPFNr(tp,fname,operator)


/* these are the function prototypes + definitions: */
#define MAX(a, b)  ( (a) > (b) ) ? (a) : (b)
FN(TYPENM,max,MAX)
#undef MAX
#define MIN(a, b)  ( (a) < (b) ) ? (a) : (b)
FN(TYPENM,min,MIN)
#undef MIN
OPFN(TYPENM,sum,+)
OPFN(TYPENM,prod,*)
#if INTTYPE!=0
OPFN(TYPENM,land,&&)
OPFN(TYPENM,band,&)
OPFN(TYPENM,lor,||)
OPFN(TYPENM,bor,|)
OPFN(TYPENM,bxor,^)
#define XOR(a, b)  ( (!a) ^ (!b) ) 
FN(TYPENM,lxor,XOR)
#undef XOR
#endif

#define OPARRAYr(tp) oplist_##tp
#define OPARRAY(tp) OPARRAYr(tp)

tMPI_Op_fn OPARRAY(TYPENM)[] = {
    FNAME(TYPENM,min),
    FNAME(TYPENM,max),
    FNAME(TYPENM,sum),
    FNAME(TYPENM,prod),
#if INTTYPE
    FNAME(TYPENM,land),
    FNAME(TYPENM,band),
    FNAME(TYPENM,lor),
    FNAME(TYPENM,bor),
    FNAME(TYPENM,lxor),
    FNAME(TYPENM,bxor)
#else
    0,
    0,
    0,
    0,
    0,
    0
#endif
};


#undef FNAME
#undef FNAMEr
#undef OPARRAYr
#undef OPARRAY
#undef FN
#undef FNr
#undef OPFN
#undef OPFNr

#undef TYPE
#undef TYPENM
#undef INTTYPE

#endif
