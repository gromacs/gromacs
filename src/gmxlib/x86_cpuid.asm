	;; this file must be processed with a version
	;; of nasm that supports the extended 3dnow instructions.
	;; you can find a binary of such a version on the
	;; gromacs homepage.
segment .text
	global x86_cpuid		;  issues the cpuid instruction with supplied args
x86_cpuid:	
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
