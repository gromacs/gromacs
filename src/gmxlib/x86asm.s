	;; this file must be processed with a version
	;; of nasm that supports the extended 3dnow instructions.
	;; you can find a binary of such a version on the
	;; gromacs homepage.
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
	emms
	xorps xmm0,xmm0
	emms
	ret
	
	global check3dnow       ;  tries to issue a simple 3DNOW instruction 	
check3dnow:	
	femms
	pfmul mm0,mm0
	femms
	ret
	

segment .text
	
	global vecrecip_3dnow
vecrecip_3dnow
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
        shr ecx, 2 
        jecxz .tail
        emms	
.mainloop:	
        movq mm0,[eax]
	add eax, byte 8
        pfrcp mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax, byte 8
        pfrcp mm2,mm0
	pswapd mm0,mm0
        pfrcp mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrcp mm6,mm4
	pswapd mm4,mm4
	pfrsqit1 mm0,mm1
	punpckldq mm5,mm6	
	pfrcpit2 mm0,mm1
        movq [ebx],mm0
	pfrsqit1 mm4,mm5
	add ebx, byte 8
	pfrcpit2 mm4,mm5	
        movq [ebx],mm4
	add ebx, byte 8	
        dec ecx
        jecxz .tail
        jmp short .mainloop
.tail:
        mov ecx,edx
        and ecx,3
        jecxz .end
.tailloop:	
        movd mm0,[eax]
	add eax, byte 4
        pfrcp mm1,mm0
        pfrsqit1 mm0,mm1
        pfrcpit2 mm0,mm1
        movd [ebx],mm0	
	add ebx, byte 4
	dec ecx
	jecxz .end
	jmp short .tailloop
.end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret

		
segment .text

	global vecinvsqrt_3dnow
vecinvsqrt_3dnow
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
        shr ecx, 2 
        jecxz .tail
        emms	
.mainloop:	
        movq mm0,[eax]
	add eax, byte 8
        pfrsqrt mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax, byte 8
        pfrsqrt mm2,mm0
	pswapd mm0,mm0
        pfrsqrt mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrsqrt mm6,mm4
	movq mm3,mm1
	pswapd mm4,mm4
	pfmul mm1,mm1
	punpckldq mm5,mm6	
	pfrsqit1 mm1,mm0
	movq mm7,mm5	
	pfrcpit2 mm1,mm3
	pfmul mm5,mm5
        movq [ebx],mm1
	pfrsqit1 mm5,mm4
	add ebx, byte 8
	pfrcpit2 mm5,mm7	
        movq [ebx],mm5
	add ebx, byte 8	
        dec ecx
        jecxz .tail
        jmp short .mainloop
.tail:
        mov ecx,edx
        and ecx,3
        jecxz .end
.tailloop:	
        movd mm0,[eax]
	add eax, byte 4
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0
        pfrcpit2 mm1,mm2
        movd [ebx],mm1		
	add ebx, byte 4
	dec ecx
	jecxz .end
	jmp short .tailloop
.end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

segment .data	
sse_three
	dd 3.0
	dd 3.0
	dd 3.0
	dd 3.0
sse_half
	dd -0.5
	dd -0.5
	dd -0.5
	dd -0.5
sse_two
	dd 2.0
	dd 2.0
	dd 2.0
	dd 2.0

			
segment .text
align 16
	global vecinvsqrt_sse
vecinvsqrt_sse
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
	movups xmm6,[sse_three]
	movups xmm7,[sse_half]
        shr ecx, 3
        jecxz .iter4
        emms	
.loop8:	
	movaps xmm0,[eax]
	add eax, byte 16
	rsqrtps xmm1,xmm0
	movaps xmm2,[eax]
	add eax, byte 16
	rsqrtps xmm3,xmm2
	mulps xmm0,xmm1
        mulps xmm2,xmm3
	mulps xmm0,xmm1
        mulps xmm2,xmm3
	subps xmm0,xmm6
	subps xmm2,xmm6
	mulps xmm0,xmm1
	mulps xmm2,xmm3
	mulps xmm0,xmm7
	mulps xmm2,xmm7
	movaps [ebx],xmm0
	add ebx, byte 16
	movaps [ebx],xmm2
	add ebx, byte 16
        dec ecx
        jecxz .iter4
        jmp .loop8
.iter4:
        mov ecx,edx
        and ecx,4
        jecxz .iter2
	movaps xmm0,[eax]
	add eax, byte 16
	rsqrtps xmm1,xmm0
	mulps xmm0,xmm1
	mulps xmm0,xmm1
	subps xmm0,xmm6
	mulps xmm0,xmm1
	mulps xmm0,xmm7
	movaps [ebx],xmm0
	add ebx, byte 16        
.iter2:
        mov ecx,edx
        and ecx,2
        jecxz .iter1
	movlps xmm0,[eax]
	add eax, byte 8
	rsqrtps xmm1,xmm0
	mulps xmm0,xmm1
	mulps xmm0,xmm1
	subps xmm0,xmm6
	mulps xmm0,xmm1
	mulps xmm0,xmm7
	movlps [ebx],xmm0
	add ebx, byte 8     
.iter1:
        mov ecx,edx
        and ecx,1
        jecxz .end
	movss xmm0,[eax]
	rsqrtss xmm1,xmm0
	mulss xmm0,xmm1
	mulss xmm0,xmm1
	subss xmm0,xmm6
	mulss xmm0,xmm1
	mulss xmm0,xmm7
	movss [ebx],xmm0        
.end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	
	global vecrecip_sse
vecrecip_sse
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
        mov edx, ecx
	movups xmm6,[sse_two]
        shr ecx, 3
        jecxz .iter4
        emms	
.loop8:	
	movaps xmm0,[eax]
	add eax, byte 16
	rcpps xmm1,xmm0
	movaps xmm3,[eax]
	add eax, byte 16
	rcpps xmm4,xmm3
	movaps xmm2,xmm6
	mulps xmm0,xmm1
	movaps xmm5,xmm6	
	subps xmm2,xmm0
	mulps xmm3,xmm4
	mulps xmm2,xmm1	
	subps xmm5,xmm3	
	movaps [ebx],xmm2
	mulps xmm5,xmm4
	add ebx, byte 16
	movaps [ebx],xmm5
	add ebx, byte 16
        dec ecx
        jecxz .iter4
        jmp .loop8
.iter4:
        mov ecx,edx
        and ecx,4
        jecxz .iter2
	movaps xmm0,[eax]
	add eax, byte 16
	rcpps xmm1,xmm0
	movaps xmm2,xmm6
	mulps xmm0,xmm1		
	subps xmm2,xmm0
	mulps xmm2,xmm1
	movaps [ebx],xmm2
	add ebx, byte 16        
.iter2:
        mov ecx,edx
        and ecx,2
        jecxz .iter1
	movlps xmm0,[eax]
	add eax, byte 8
	rcpps xmm1,xmm0
	movaps xmm2,xmm6
	mulps xmm0,xmm1		
	subps xmm2,xmm0
	mulps xmm2,xmm1
	movlps [ebx],xmm2
	add ebx, byte 8     
.iter1:
        mov ecx,edx
        and ecx,1
        jecxz .end
	movss xmm0,[eax]
	rcpss xmm1,xmm0
	movss xmm2,xmm6
	mulss xmm0,xmm1		
	subss xmm2,xmm0
	mulss xmm2,xmm1
	movss [ebx],xmm2        
.end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	
