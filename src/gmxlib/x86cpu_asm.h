segment .text
	global gmxcpuid		;  issues the cpuid instruction with supplied args
gmxcpuid:	
	push ebp
	mov  ebp,esp
	push edi
	push ebx
	push ecx
	push edx
	mov  eax, [ebp+8]	
	cpuid
	mov  edi, [ebp+12]
	mov  [edi],eax
	mov  edi, [ebp+16]
	mov  [edi],ebx
	mov  edi, [ebp+20]
	mov  [edi],ecx
	mov  edi, [ebp+24]
	mov  [edi],edx
	pop edx
	pop ecx
	pop ebx
	pop edi
	mov esp, ebp
	pop ebp
	ret
	
	global checksse		;  tries to issue a simple SSE instruction
checksse:	
	xorps xmm0,xmm0
	ret
	
	global check3dnow       ;  tries to issue a simple 3DNOW instruction 	
check3dnow:	
	pfmul mm0,mm0
	ret
	
