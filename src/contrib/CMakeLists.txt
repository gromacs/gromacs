set(CONTRIB_PROGRAMS 
     #add here any programs you want to compile

)

foreach(PROG ${CONTRIB_PROGRAMS})
        add_executable(${PROG} ${PROG}.c ${NGMX_COMMON_SOURCE})
	set_target_properties(${PROG} PROPERTIES OUTPUT_NAME "${PROG}${GMX_BINARY_SUFFIX}")
endforeach(PROG)
