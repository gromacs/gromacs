set(CONTRIB_PROGRAMS 
     #add here any programs you want to compile

)

mark_as_advanced(CLEAR GMX_CPU_ACCELERATION GMX_MPI GMX_FFT_LIBRARY
    GMX_QMMM_PROGRAM GMX_THREAD_MPI GMX_DOUBLE)

foreach(PROG ${CONTRIB_PROGRAMS})
        add_executable(${PROG} ${PROG}.c ${NGMX_COMMON_SOURCE})
	set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG)
