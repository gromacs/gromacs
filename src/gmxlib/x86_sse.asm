;;
;;                This source code is part of
;;
;;                 G   R   O   M   A   C   S
;;
;;          GROningen MAchine for Chemical Simulations
;;
;;                        VERSION 3.0
;;
;; Copyright (c) 1991-2001
;; BIOSON Research Institute, Dept. of Biophysical Chemistry
;; University of Groningen, The Netherlands
;;
;; This program is free software; you can redistribute it and/or
;; modify it under the terms of the GNU General Public License
;; as published by the Free Software Foundation; either version 2
;; of the License, or (at your option) any later version.
;;
;; If you want to redistribute modifications, please consider that
;; scientific software is very special. Version control is crucial -
;; bugs must be traceable. We will be happy to consider code for
;; inclusion in the official distribution, but derived work must not
;; be called official GROMACS. Details are found in the README & COPYING
;; files - if they are missing, get the official version at www.gromacs.org.
;;
;; To help us fund GROMACS development, we humbly ask that you cite
;; the papers on the package - you can find them in the top README file.
;;
;; Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
;;
;; And Hey:
;; GROup of MAchos and Cynical Suckers

; NASM macro set to make interfacing to 32-bit programs easier -*- nasm -*-
%imacro proc 1                  ; begin a procedure definition
%push proc
          global %1
%1:       push ebp
          mov ebp,esp
%assign %$arg 8
%define %$procname %1
%endmacro



%imacro arg 0-1 4               ; used with the argument name as a label
%00       equ %$arg
%assign %$arg %1+%$arg
%endmacro



%imacro endproc 0
%ifnctx proc
%error Mismatched `endproc'/`proc'

%else
          leave
          ret
__end_%$procname:               ; useful for calculating function size

%pop
%endif
%endmacro


segment .data

sse_minushalf
	dd -0.5
	dd -0.5
	dd -0.5
	dd -0.5
sse_half
	dd 0.5
	dd 0.5
	dd 0.5
	dd 0.5
sse_two
	dd 2.0
	dd 2.0
	dd 2.0
	dd 2.0
sse_three
	dd 3.0
	dd 3.0
	dd 3.0
	dd 3.0
sse_six
	dd 6.0
	dd 6.0
	dd 6.0
	dd 6.0
sse_twelve
	dd 12.0
	dd 12.0
	dd 12.0
	dd 12.0


segment .text

	global checksse		;  tries to issue a simple SSE instruction
checksse:
	emms
	xorps xmm0,xmm0
	emms
	ret

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
	
	
proc inl0100_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg		
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.dx          equ    48
.dy          equ    64
.dz          equ    80
.two         equ    96		
.c6          equ   112
.c12         equ   128
.six         equ   144
.twelve      equ   160		 
.vnbtot      equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.ntia	     equ   280	
.innerjjnr   equ   284
.innerk      equ   288
.salign	     equ   292								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 296		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm1, [sse_two]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movaps [esp + .two], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vnbtot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	


	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4


	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:			
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 296
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl0110_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg		
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.dx          equ    48
.dy          equ    64
.dz          equ    80
.two         equ    96		
.c6          equ   112
.c12         equ   128
.six         equ   144
.twelve      equ   160		 
.vnbtot      equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.shX	     equ   280
.shY         equ   284
.shZ         equ   288
.ntia	     equ   292	
.innerjjnr0  equ   296
.innerjjnr   equ   300
.innerk0     equ   304
.innerk      equ   308
.salign	     equ   312						
.nsvdwc      equ   316
.nscoul      equ   320
.nsvdw       equ   324
.solnr	     equ   328		
	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 332		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm1, [sse_two]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movaps [esp + .two], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movlps xmm0, [eax + ebx*4]	; getting the shiftvector
	movss xmm1, [eax + ebx*4 + 8] 
	movlps [esp + .shX], xmm0
	movss [esp + .shZ], xmm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	add   [ebp + %$nsatoms], dword 12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + .nsvdwc], edx
	mov   [esp + .nscoul], eax
	mov   [esp + .nsvdw], ecx

	; clear vnbtot 
	xorps xmm4, xmm4
	movaps [esp + .vnbtot], xmm4
	mov   [esp + .solnr],  ebx
		
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]

	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testvdw
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	; quad-unroll innerloop here
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdwc
	jmp   .checksingle_vdwc
.dopair_vdwc:	

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	


	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4


	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdwc:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdwc
	jmp    .updateouterdata_vdwc
.dosingle_vdwc:			
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdwc:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testvdw
	jmp .mno_vdwc
.testvdw
	mov  ebx,  [esp + .nscoul]
	add  [esp + .solnr], dword ebx

	mov  ecx, [esp + .nsvdw]
	cmp  ecx, byte 0
	jnz  .mno_vdw
	jmp  .last_mno
.mno_vdw:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdw
	jmp   .checksingle_vdw
.dopair_vdw:	

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	


	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4


	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdw:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdw
	jmp    .updateouterdata_vdw
.dosingle_vdw:			
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdw:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5
	
	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 332
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl0300_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale      arg
%$VFtab         arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.dx          equ    48
.dy          equ    64
.dz          equ    80
.two	     equ    96
.tabscale    equ   112
.c6          equ   128
.c12         equ   144
.fs          equ   160
.vnbtot      equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.ntia	     equ   280	
.innerjjnr   equ   284
.innerk      equ   288
.salign	     equ   292								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 296		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear tot potential and i forces
	xorps xmm4, xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	; dispersion
	movlps xmm5, [esi + eax*4 + 0]
	movlps xmm7, [esi + ecx*4 + 0]
	movhps xmm5, [esi + ebx*4 + 0]
	movhps xmm7, [esi + edx*4 + 0] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	; dispersion
	movlps xmm5, [esi + ecx*4 + 0]
	movhps xmm5, [esi + edx*4 + 0]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	; dispersion
	movlps xmm4, [esi + ebx*4 + 0]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 296
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl0310_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg		
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
%$tabscale      arg
%$VFtab         arg
%$nsatoms       arg
	;; stack offsets for local variables
        ;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.dx          equ    48
.dy          equ    64
.dz          equ    80
.two         equ    96   
.tabscale    equ   112
.c6          equ   128
.c12         equ   144
.fs          equ   160
.vnbtot      equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.shX	     equ   280
.shY         equ   284
.shZ         equ   288
.ntia	     equ   292	
.innerjjnr0  equ   296
.innerjjnr   equ   300
.innerk0     equ   304
.innerk      equ   308
.salign	     equ   312						
.nsvdwc      equ   316
.nscoul      equ   320
.nsvdw       equ   324
.solnr	     equ   328
        push eax      
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 332		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three], xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movlps xmm0, [eax + ebx*4]	; getting the shiftvector
	movss xmm1, [eax + ebx*4 + 8] 
	movlps [esp + .shX], xmm0
	movss [esp + .shZ], xmm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	add   [ebp + %$nsatoms], dword 12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + .nsvdwc], edx
	mov   [esp + .nscoul], eax
	mov   [esp + .nsvdw], ecx

	; clear vnbtot 
	xorps xmm4, xmm4
	movaps [esp + .vnbtot], xmm4
	mov   [esp + .solnr],  ebx
		
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]

	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testvdw
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	; dispersion
	movlps xmm5, [esi + eax*4 + 0]
	movlps xmm7, [esi + ecx*4 + 0]
	movhps xmm5, [esi + ebx*4 + 0]
	movhps xmm7, [esi + edx*4 + 0] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdwc
	jmp   .checksingle_vdwc
.dopair_vdwc:	
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	; dispersion
	movlps xmm5, [esi + ecx*4 + 0]
	movhps xmm5, [esi + edx*4 + 0]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 

	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdwc:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdwc
	jmp    .updateouterdata_vdwc
.dosingle_vdwc:
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	; dispersion
	movlps xmm4, [esi + ebx*4 + 0]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdwc:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testvdw
	jmp .mno_vdwc
.testvdw
	mov  ebx,  [esp + .nscoul]
	add  [esp + .solnr], dword ebx

	mov  ecx, [esp + .nsvdw]
	cmp  ecx, byte 0
	jnz  .mno_vdw
	jmp  .last_mno
.mno_vdw:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	; dispersion
	movlps xmm5, [esi + eax*4 + 0]
	movlps xmm7, [esi + ecx*4 + 0]
	movhps xmm5, [esi + ebx*4 + 0]
	movhps xmm7, [esi + edx*4 + 0] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdw
	jmp   .checksingle_vdw
.dopair_vdw:	
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	; dispersion
	movlps xmm5, [esi + ecx*4 + 0]
	movhps xmm5, [esi + edx*4 + 0]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 

	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdw:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdw
	jmp    .updateouterdata_vdw
.dosingle_vdw:
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 3

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	; dispersion
	movlps xmm4, [esi + ebx*4 + 0]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5		
.updateouterdata_vdw:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5
	
	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw	
.last_mno:	

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 332
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl1000_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg		
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.vctot       equ   112
.fix         equ   128
.fiy         equ   144
.fiz         equ   160
.half        equ   176
.three       equ   192
.is3         equ   208
.ii3         equ   212
.innerjjnr   equ   216
.innerk      equ   220		
.salign	     equ   224							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 228		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unrolled innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm5, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b	      
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm5
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]	;x1 y1 - -
	movlps xmm5, [esi + ecx*4]	;x3 y3 - -
	movss xmm2, [esi + eax*4 + 8]	;z1 -  - -
	movss xmm6, [esi + ecx*4 + 8]   ;z3 -  - -

	movhps xmm4, [esi + ebx*4]	;x1 y1 x2 y2
	movhps xmm5, [esi + edx*4]	;x3 y3 x4 y4

	movss xmm0, [esi + ebx*4 + 8]	;z2 - - -
	movss xmm1, [esi + edx*4 + 8]	;z4 - - -

	shufps xmm2, xmm0, 0b		;z1 z1 z2 z2
	shufps xmm6, xmm1, 0b		;z3 z3 z4 z4
	
	movaps xmm0, xmm4		;x1 y1 x2 y2	
	movaps xmm1, xmm4		;x1 y1 x2 y2

	shufps xmm2, xmm6, 10001000b	;z1 z2 z3 z4
	
	shufps xmm0, xmm5, 10001000b	;x1 x2 x3 x4
	shufps xmm1, xmm5, 11011101b	;y1 y2 y3 y4		

	mov    edi, [ebp + %$faction]

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's 
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b
	
	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	
.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	
	
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	mov    edi, [ebp + %$faction]
	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec dword ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 228
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl1010_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg		
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.vctot       equ   112
.fix         equ   128
.fiy         equ   144
.fiz         equ   160
.half        equ   176
.three       equ   192
.is3         equ   208
.ii3         equ   212
.shX	     equ   216
.shY         equ   220
.shZ         equ   224
.ntia	     equ   228	
.innerjjnr0  equ   232
.innerk0     equ   236
.innerjjnr   equ   240
.innerk      equ   244		
.salign	     equ   248						
.nscoul      equ   252
.solnr	     equ   256		
	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 260		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	add   [ebp + %$nsatoms], dword 8

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 
	movss [esp + .shX], xmm0
	movss [esp + .shY], xmm1
	movss [esp + .shZ], xmm2

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	mov   ecx, [eax]
	add   [ebp + %$nsatoms], dword 12
	mov   [esp + .nscoul], ecx	

	; clear vctot 
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	mov   [esp + .solnr], ebx

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]
	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nscoul]
	cmp   ecx, dword 0
	jnz   .mno_coul
	jmp   .last_mno
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner

.unroll_coul_loop:	
	;; quad-unrolled innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm5, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b	      
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm5
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	mov    edi, [ebp + %$faction]

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5


	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_coul
	jmp   .checksingle_coul
.dopair_coul:	
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's 
	movss xmm3, [edi + eax*4]
	movss xmm4, [edi + eax*4 + 4]
	movss xmm5, [edi + eax*4 + 8]
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	subps  xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b
	
	movss xmm3, [edi + ebx*4]
	movss xmm4, [edi + ebx*4 + 4]
	movss xmm5, [edi + ebx*4 + 8]
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	subps  xmm5, xmm2
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5

.checksingle_coul:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_coul
	jmp    .updateouterdata_coul
.dosingle_coul:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	
	
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	mov    edi, [ebp + %$faction]
	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
        subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_coul:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .last_mno
	jmp .mno_coul
	
.last_mno:	
	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 260
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl1020_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.vctot       equ   352
.fixO        equ   368
.fiyO        equ   384
.fizO        equ   400
.fixH1       equ   416
.fiyH1       equ   432
.fizH1       equ   448
.fixH2       equ   464
.fiyH2       equ   480
.fizH2       equ   496
.fjx	     equ   512
.fjy         equ   528
.fjz         equ   544
.half        equ   560
.three       equ   576
.is3         equ   592
.ii3         equ   596
.innerjjnr   equ   600
.innerk      equ   604
.salign	     equ   608								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 612		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4	

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	; start with rsqO - seed in xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm7, xmm4	; rinvO in xmm7
	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm6, xmm4	; rinvH1 in xmm6
	; rsqH2 - seed in xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm5, xmm4	; rinvH2 in xmm5

	; do O interactions
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	; xmm7=rinv, xmm4=rinvsq
	mulps  xmm7, [esp + .qqO]	;xmm7=vcoul
	
	mulps  xmm4, xmm7	; total fsO in xmm4

	addps  xmm7, [esp + .vctot]
	
	movaps [esp + .vctot], xmm7

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H1 interactions
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	; xmm6=rinv, xmm4=rinvsq
	mulps  xmm6, [esp + .qqH]	;xmm6=vcoul
	mulps  xmm4, xmm6		; total fsH1 in xmm4
	
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	movaps [esp + .vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H2 interactions
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	; xmm5=rinv, xmm4=rinvsq
	mulps  xmm5, [esp + .qqH]	;xmm5=vcoul
	mulps  xmm4, xmm5		; total fsH1 in xmm4
	
	addps  xmm5, [esp + .vctot]

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	movaps [esp + .vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm3, [esp + .qqO]

	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	movaps [esp + .vctot], xmm3

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 	
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 612
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl1030_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use	
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.vctot       equ   768		
.fixO        equ   784
.fiyO        equ   800
.fizO        equ   816
.fixH1       equ   832
.fiyH1       equ   848
.fizH1       equ   864
.fixH2       equ   880
.fiyH2       equ   896
.fizH2       equ   912
.fjxO	     equ   928
.fjyO        equ   944
.fjzO        equ   960
.fjxH1	     equ   976
.fjyH1       equ   992
.fjzH1       equ  1008
.fjxH2	     equ  1024
.fjyH2       equ  1040
.fjzH2       equ  1056
.half        equ  1072
.three       equ  1088
.rsqOO       equ  1104
.rsqOH1      equ  1120
.rsqOH2      equ  1136
.rsqH1O      equ  1152
.rsqH1H1     equ  1168
.rsqH1H2     equ  1184
.rsqH2O      equ  1200
.rsqH2H1     equ  1216
.rsqH2H2     equ  1232
.rinvOO      equ  1248
.rinvOH1     equ  1264
.rinvOH2     equ  1280
.rinvH1O     equ  1296
.rinvH1H1    equ  1312
.rinvH1H2    equ  1328
.rinvH2O     equ  1344
.rinvH2H1    equ  1360
.rinvH2H2    equ  1376
.is3         equ  1392
.ii3         equ  1396
.innerjjnr   equ  1400
.innerk      equ  1404
.salign	     equ  1408							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1412		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	mulps  xmm7, [esp + .qqOO]
	mulps  xmm0, xmm7	
	addps  xmm7, [esp + .vctot] 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsOH1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsOH2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsH1O 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH1H1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsOH2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsH2O 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH2H1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH2H2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps [esp + .vctot], xmm7
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]
		
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	; xmm0=rinvsq
	;; fetch charges to xmm4 (temporary)
	movss   xmm4, [esp + .qqOO]

	movhps  xmm4, [esp + .qqOH]

	mulps   xmm3, xmm4	; xmm3=vcoul
	mulps   xmm0, xmm3	;  total fscal
	addps   xmm3, [esp + .vctot]
	movaps  [esp + .vctot], xmm3	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.

	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	;; assemble charges in xmm6
	xorps   xmm6, xmm6
	; do coulomb interaction
	movaps  xmm0, xmm3
	movss   xmm6, [esp + .qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, [esp + .qqHH]
	mulps   xmm0, xmm0	;  rinvsq
	mulps   xmm4, xmm4	;  rinvsq
	mulps   xmm3, xmm6	;  vcoul
	mulps   xmm7, xmm6	;  vcoul
	movaps  xmm2, xmm3
	addps   xmm2, xmm7	;  total vcoul
	mulps   xmm0, xmm3	;  fscal
	
	addps   xmm2, [esp + .vctot]
	mulps   xmm7, xmm4	;  fscal
	movaps  [esp + .vctot], xmm2
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do forces H2 - j water
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 	

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1412
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc







proc inl1100_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96	
.c6          equ   112
.c12         equ   128
.six         equ   144
.twelve      equ   160		 
.vctot       equ   176
.vnbtot      equ   192
.fix         equ   208
.fiy         equ   224
.fiz         equ   240
.half        equ   256
.three       equ   272
.is3         equ   288
.ii3         equ   292
.ntia	     equ   296	
.innerjjnr   equ   300
.innerk      equ   304
.salign	     equ   308								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, dword 312		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00001100b 
	shufps xmm3, xmm3, 01011000b ; xmm3(0,1) has the charges.

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, dword 312
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl2100_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg	
%$crf		arg	
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96	
.c6          equ   112
.c12         equ   128
.six         equ   144
.twelve      equ   160		 
.vctot       equ   176
.vnbtot      equ   192
.fix         equ   208
.fiy         equ   224
.fiz         equ   240
.half        equ   256
.three       equ   272
.two         equ   288
.krf	     equ   304	 
.crf	     equ   320	 
.is3         equ   336
.ii3         equ   340
.ntia	     equ   344
.innerjjnr   equ   348
.innerk      equ   352
.salign	     equ   356								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, dword 360		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	shufps xmm6, xmm6, 0b
	movaps [esp + .krf], xmm5
	movaps [esp + .crf], xmm6

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4
	
	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq
	movaps xmm1, xmm4
	subps  xmm6, [esp + .crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm6, xmm3	; xmm6=vcoul=qq*(rinv+krsq)
	mulps  xmm7, [esp + .two]
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6
	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00001100b 
	shufps xmm3, xmm3, 01011000b ; xmm3(0,1) has the charges.

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq
	movaps xmm1, xmm4
	subps  xmm6, [esp + .crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm6, xmm3	; xmm6=vcoul=qq*(rinv+krsq-crf)
	mulps  xmm7, [esp + .two]	
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0	
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq
	movaps xmm1, xmm4
	subps  xmm6, [esp + .crf]	
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm6, xmm3	; xmm6=vcoul
	mulps  xmm7, [esp + .two]
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, dword 360
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl2000_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg	
%$crf           arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96	
.vctot       equ   112
.fix         equ   128
.fiy         equ   144
.fiz         equ   160
.half        equ   176
.three       equ   192
.two         equ   208
.krf	     equ   224	 
.crf	     equ   240	 
.is3         equ   256
.ii3         equ   260
.innerjjnr   equ   264
.innerk      equ   268
.salign	     equ   272								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, dword 276		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	movaps [esp + .krf], xmm5
	shufps xmm6, xmm6, 0b
	movaps [esp + .crf], xmm6

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4
	
	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq

	subps  xmm6, [esp + .crf] ;  xmm6=rinv+krsq-crf

	mulps  xmm6, xmm3	; xmm6=vcoul=qq*(rinv+krsq)
	mulps  xmm7, [esp + .two]

	subps  xmm0, xmm7
	mulps  xmm3, xmm0	
	mulps  xmm4, xmm3	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00001100b 
	shufps xmm3, xmm3, 01011000b ; xmm3(0,1) has the charges.	

	mov edi, [ebp + %$pos]	
				
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]

	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq

	subps  xmm6, [esp + .crf] ;  xmm6=rinv+krsq-crf

	mulps  xmm6, xmm3	; xmm6=vcoul=qq*(rinv+krsq-crf)
	mulps  xmm7, [esp + .two]	

	subps  xmm0, xmm7
	mulps  xmm3, xmm0	

	mulps  xmm4, xmm3	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge		
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	movaps xmm7, [esp + .krf]
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	mulps  xmm7, xmm4	;  xmm7=krsq
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;  xmm6=rinv+krsq

	subps  xmm6, [esp + .crf] ;  xmm6=rinv+krsq-crf

	mulps  xmm6, xmm3	; xmm6=vcoul
	mulps  xmm7, [esp + .two]

	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	mulps  xmm4, xmm3	; xmm4=total fscal
	addps  xmm6, [esp + .vctot]
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm6

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, dword 276
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc





proc inl1110_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96	
.c6          equ   112
.c12         equ   128
.two         equ   144
.six         equ   160
.twelve      equ   176		 
.vctot       equ   192
.vnbtot      equ   208
.fix         equ   224
.fiy         equ   240
.fiz         equ   256
.half        equ   272
.three       equ   288
.is3         equ   304
.ii3         equ   308
.shX	     equ   312
.shY         equ   316
.shZ         equ   320
.ntia	     equ   324	
.innerjjnr0  equ   328
.innerk0     equ   332
.innerjjnr   equ   336
.innerk      equ   340
.salign	     equ   344							
.nsvdwc      equ   348
.nscoul      equ   352
.nsvdw       equ   356
.solnr	     equ   360		

	push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 364		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movups xmm3, [sse_six]
	movups xmm4, [sse_twelve]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three], xmm2
	movaps [esp + .six],  xmm3
	movaps [esp + .twelve], xmm4

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movlps xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 8] 
	movlps [esp + .shX], xmm0
	movss [esp + .shZ], xmm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	add   [ebp + %$nsatoms], dword 12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + .nsvdwc], edx
	mov   [esp + .nscoul], eax
	mov   [esp + .nsvdw], ecx
		
	;; clear potential
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	mov   [esp + .solnr],  ebx

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]
	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]
	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdwc
	jmp   .checksingle_vdwc
.dopair_vdwc:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdwc:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdwc
	jmp    .updateouterdata_vdwc
.dosingle_vdwc:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdwc:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5


	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul:
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner

.unroll_coul_loop:	
	;; quad-unrolled innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm5, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b	      
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm5
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	mov    edi, [ebp + %$faction]

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5


	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_coul
	jmp   .checksingle_coul
.dopair_coul:	
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + .iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's 
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b
	
	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_coul:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_coul
	jmp    .updateouterdata_coul
.dosingle_coul:			
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	; xmm3(0) has the charge	
	
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + .iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	mov    edi, [ebp + %$faction]
	movaps xmm5, [esp + .vctot]
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm4, xmm3	; xmm4=fscal
	addps  xmm5, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

.updateouterdata_coul:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw:
	mov  ecx, [esp + .nsvdw]
	cmp  ecx, byte 0
	jnz  .mno_vdw
	jmp  .last_mno
.mno_vdw:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdw
	jmp   .checksingle_vdw
.dopair_vdw:	

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	


	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4


	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdw:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdw
	jmp    .updateouterdata_vdw
.dosingle_vdw:			
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdw:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5
	
	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
.last_mno:	
	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 364
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl1120_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.c6          equ   352
.c12         equ   368
.six         equ   384
.twelve      equ   400		 
.vctot       equ   416
.vnbtot      equ   432
.fixO        equ   448
.fiyO        equ   464
.fizO        equ   480
.fixH1       equ   496
.fiyH1       equ   512
.fizH1       equ   528
.fixH2       equ   544
.fiyH2       equ   560
.fizH2       equ   576
.fjx	     equ   592
.fjy         equ   608
.fjz         equ   624
.half        equ   640
.three       equ   656
.is3         equ   672
.ii3         equ   676
.ntia	     equ   680	
.innerjjnr   equ   684
.innerk      equ   688
.salign	     equ   692								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 696		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	; start with rsqO - seed in xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm7, xmm4	; rinvO in xmm7
	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm6, xmm4	; rinvH1 in xmm6
	; rsqH2 - seed in xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm5, xmm4	; rinvH2 in xmm5

	; do O interactions
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	; xmm7=rinv, xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm7, [esp + .qqO]	;xmm7=vcoul
	
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	; vnb=vnb12-vnb6		
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm7	
	mulps  xmm4, xmm2	; total fsO in xmm4

	addps  xmm7, [esp + .vctot]
	
	movaps [esp + .vnbtot], xmm3
	movaps [esp + .vctot], xmm7

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H1 interactions
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	; xmm6=rinv, xmm4=rinvsq
	mulps  xmm6, [esp + .qqH]	;xmm6=vcoul
	mulps  xmm4, xmm6		; total fsH1 in xmm4
	
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	movaps [esp + .vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H2 interactions
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	; xmm5=rinv, xmm4=rinvsq
	mulps  xmm5, [esp + .qqH]	;xmm5=vcoul
	mulps  xmm4, xmm5		; total fsH1 in xmm4
	
	addps  xmm5, [esp + .vctot]

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	movaps [esp + .vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	xorps xmm6, xmm6
	mov esi, [ebp + %$type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + %$nbfp]
	shl ebx, 1	
	add ebx, [esp + .ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 11111100b
	shufps xmm7, xmm7, 11111101b
	movaps [esp + .c6], xmm6
	movaps [esp + .c12], xmm7

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	movaps xmm3, [esp + .qqO]
	mulss  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm3, xmm0	; xmm3=vcoul
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulss  xmm1, [esp + .six]
	mulss  xmm2, [esp + .twelve]
	subss  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm3, [esp + .vctot]

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	movaps [esp + .vctot], xmm3
	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 696
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl1130_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.c6          equ   768
.c12         equ   784
.six         equ   800
.twelve      equ   816		 
.vctot       equ   832
.vnbtot      equ   848
.fixO        equ   864
.fiyO        equ   880
.fizO        equ   896
.fixH1       equ   912
.fiyH1       equ   928
.fizH1       equ   944
.fixH2       equ   960
.fiyH2       equ   976
.fizH2       equ   992
.fjxO	     equ  1008
.fjyO        equ  1024
.fjzO        equ  1040
.fjxH1	     equ  1056
.fjyH1       equ  1072
.fjzH1       equ  1088
.fjxH2	     equ  1104
.fjyH2       equ  1120
.fjzH2       equ  1136
.half        equ  1152
.three       equ  1168
.rsqOO       equ  1184
.rsqOH1      equ  1200
.rsqOH2      equ  1216
.rsqH1O      equ  1232
.rsqH1H1     equ  1248
.rsqH1H2     equ  1264
.rsqH2O      equ  1280
.rsqH2H1     equ  1296
.rsqH2H2     equ  1312
.rinvOO      equ  1328
.rinvOH1     equ  1344
.rinvOH2     equ  1360
.rinvH1O     equ  1376
.rinvH1H1    equ  1392
.rinvH1H2    equ  1408
.rinvH2O     equ  1424
.rinvH2H1    equ  1440
.rinvH2H2    equ  1456
.is3         equ  1472
.ii3         equ  1476
.innerjjnr   equ  1480
.innerk      equ  1484
.salign	     equ  1488							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1492		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 01010101b
	movaps [esp + .c6], xmm0
	movaps [esp + .c12], xmm1

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm7, xmm0
	mulps  xmm0, xmm0
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	; xmm1=rinvsix
	mulps  xmm7, [esp + .qqOO]
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=rinvtwelve
	mulps  xmm1, [esp + .c6]	
	mulps  xmm2, [esp + .c12]	
	movaps xmm3, xmm2
	subps  xmm3, xmm1	; xmm3=vnb12-vnb6
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	movaps [esp + .vnbtot], xmm3
	subps  xmm2, xmm1
	addps  xmm2, xmm7
	addps  xmm7, [esp + .vctot]
	mulps  xmm0, xmm2	
 
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsOH1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsOH2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsH1O 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH1H1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsOH2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqOH]
	mulps xmm0, xmm1	; fsH2O 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH2H1 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm1, xmm0
	mulps xmm0, xmm0
	mulps xmm1, [esp + .qqHH]
	mulps xmm0, xmm1	; fsH2H2 
	addps xmm7, xmm1	; add to local vctot.
	movaps xmm1, xmm0
	movaps [esp + .vctot], xmm7
	movaps xmm2, xmm0
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]
		
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	; xmm0=rinvsq
	;; fetch charges to xmm4 (temporary)
	movss   xmm4, [esp + .qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, [esp + .qqOH]
	mulss   xmm1, xmm0
	mulps   xmm3, xmm4	; xmm3=vcoul
	mulss   xmm1, xmm0	;  xmm1(0)=rinvsix
	movaps  xmm2, xmm1	;  zero everything else in xmm2
	mulss   xmm2, xmm2	;  xmm2=rinvtwelve

	mulss   xmm1, [esp + .c6]
	mulss   xmm2, [esp + .c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;  vnbtot=vnb12-vnb6
	addps   xmm4, [esp + .vnbtot]
	mulss   xmm1, [esp + .six]
	mulss   xmm2, [esp + .twelve]	
	movaps  [esp + .vnbtot], xmm4
	subss   xmm2, xmm1	; fsD+fsR
	addps   xmm2, xmm3	; fsC+fsD+fsR

	addps   xmm3, [esp + .vctot]
	mulps   xmm0, xmm2	;  total fscal
	movaps  [esp + .vctot], xmm3	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.

	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	;; assemble charges in xmm6
	xorps   xmm6, xmm6
	; do coulomb interaction
	movaps  xmm0, xmm3
	movss   xmm6, [esp + .qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, [esp + .qqHH]
	mulps   xmm0, xmm0	;  rinvsq
	mulps   xmm4, xmm4	;  rinvsq
	mulps   xmm3, xmm6	;  vcoul
	mulps   xmm7, xmm6	;  vcoul
	movaps  xmm2, xmm3
	addps   xmm2, xmm7	;  total vcoul
	mulps   xmm0, xmm3	;  fscal
	
	addps   xmm2, [esp + .vctot]
	mulps   xmm7, xmm4	;  fscal
	movaps  [esp + .vctot], xmm2
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do forces H2 - j water
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1492
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl2120_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg	
%$crf		arg	
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.c6          equ   352
.c12         equ   368
.six         equ   384
.twelve      equ   400		 
.vctot       equ   416
.vnbtot      equ   432
.fixO        equ   448
.fiyO        equ   464
.fizO        equ   480
.fixH1       equ   496
.fiyH1       equ   512
.fizH1       equ   528
.fixH2       equ   544
.fiyH2       equ   560
.fizH2       equ   576
.fjx	     equ   592
.fjy         equ   608
.fjz         equ   624
.half        equ   640
.three       equ   656
.two	     equ   672
.krf	     equ   688
.crf	     equ   704
.krsqO       equ   720
.krsqH1      equ   736
.krsqH2	     equ   752	 		
.is3         equ   768
.ii3         equ   772
.ntia	     equ   776	
.innerjjnr   equ   780
.innerk      equ   784
.salign	     equ   788								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 792		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]

	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	shufps xmm6, xmm6, 0b
	movaps [esp + .krf], xmm5
	movaps [esp + .crf], xmm6
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + .krf]	
	mulps  xmm1, [esp + .krf]	
	mulps  xmm2, [esp + .krf]	

	movaps [esp + .krsqH2], xmm0
	movaps [esp + .krsqH1], xmm1
	movaps [esp + .krsqO], xmm2
	
	; start with rsqO - seed in xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm7, xmm4	; rinvO in xmm7
	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm6, xmm4	; rinvH1 in xmm6
	; rsqH2 - seed in xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm5, xmm4	; rinvH2 in xmm5

	; do O interactions
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	; xmm7=rinv, xmm4=rinvsq
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	; vnb=vnb12-vnb6		
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1	;  nb part of fs

	movaps xmm0, xmm7
	movaps xmm1, [esp + .krsqO]
	addps  xmm0, xmm1
	mulps  xmm1, [esp + .two]
	subps  xmm0, [esp + .crf] ;  xmm0=rinv+krsq-crf
	subps  xmm7, xmm1
	mulps  xmm0, [esp + .qqO]
	mulps  xmm7, [esp + .qqO]
	addps  xmm2, xmm7

	mulps  xmm4, xmm2	; total fsO in xmm4

	addps  xmm0, [esp + .vctot]
	movaps [esp + .vnbtot], xmm3
	movaps [esp + .vctot], xmm0

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H1 interactions
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	; xmm6=rinv, xmm4=rinvsq
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + .krsqH1]
	addps   xmm6, xmm0	; xmm6=rinv+krsq
	mulps   xmm0, [esp + .two]
	subps   xmm6, [esp + .crf]
	subps   xmm7, xmm0	; xmm7=rinv-2*krsq
	mulps   xmm6, [esp + .qqH] ;  vcoul
	mulps   xmm7, [esp + .qqH]
	mulps  xmm4, xmm7		; total fsH1 in xmm4
	
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	movaps [esp + .vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H2 interactions
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	; xmm5=rinv, xmm4=rinvsq
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + .krsqH2]
	addps   xmm5, xmm0	; xmm5=rinv+krsq
	mulps   xmm0, [esp + .two]
	subps   xmm5, [esp + .crf]
	subps   xmm7, xmm0	; xmm7=rinv-2*krsq
	mulps   xmm5, [esp + .qqH] ;  vcoul
	mulps   xmm7, [esp + .qqH]
	mulps  xmm4, xmm7		; total fsH2 in xmm4
	
	addps  xmm5, [esp + .vctot]

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	movaps [esp + .vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	xorps xmm6, xmm6
	mov esi, [ebp + %$type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + %$nbfp]
	shl ebx, 1	
	add ebx, [esp + .ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 11111100b
	shufps xmm7, xmm7, 11111101b
	movaps [esp + .c6], xmm6
	movaps [esp + .c12], xmm7

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	movaps xmm0, xmm4
	mulps xmm0, [esp + .krf]
	movaps [esp + .krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	mulss  xmm1, xmm4	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulss  xmm1, [esp + .six]
	mulss  xmm2, [esp + .twelve]
	subss  xmm2, xmm1

	movaps xmm1, xmm0	; xmm1=rinv
	movaps xmm3, [esp + .krsqO]
	addps  xmm0, xmm3	; xmm0=rinv+krsq
	mulps  xmm3, [esp + .two]
	subps  xmm0, [esp + .crf] ;  xmm0=rinv+krsq-crf
	subps  xmm1, xmm3	; xmm1=rinv-2*krsq
	mulps  xmm0, [esp + .qqO]	; xmm0=vcoul
	mulps  xmm1, [esp + .qqO] 	; xmm1=coul part of fs

	addps xmm2, xmm1	;  total fs
	
	mulps  xmm4, xmm2	; xmm4=total fscal
	addps  xmm0, [esp + .vctot]
	movaps [esp + .vctot], xmm0
	
	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 792
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl2130_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg
%$crf		arg
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.c6          equ   768
.c12         equ   784
.six         equ   800
.twelve      equ   816		 
.vctot       equ   832
.vnbtot      equ   848
.fixO        equ   864
.fiyO        equ   880
.fizO        equ   896
.fixH1       equ   912
.fiyH1       equ   928
.fizH1       equ   944
.fixH2       equ   960
.fiyH2       equ   976
.fizH2       equ   992
.fjxO	     equ  1008
.fjyO        equ  1024
.fjzO        equ  1040
.fjxH1	     equ  1056
.fjyH1       equ  1072
.fjzH1       equ  1088
.fjxH2	     equ  1104
.fjyH2       equ  1120
.fjzH2       equ  1136
.half        equ  1152
.three       equ  1168
.rsqOO       equ  1184
.rsqOH1      equ  1200
.rsqOH2      equ  1216
.rsqH1O      equ  1232
.rsqH1H1     equ  1248
.rsqH1H2     equ  1264
.rsqH2O      equ  1280
.rsqH2H1     equ  1296
.rsqH2H2     equ  1312
.rinvOO      equ  1328
.rinvOH1     equ  1344
.rinvOH2     equ  1360
.rinvH1O     equ  1376
.rinvH1H1    equ  1392
.rinvH1H2    equ  1408
.rinvH2O     equ  1424
.rinvH2H1    equ  1440
.rinvH2H2    equ  1456
.two         equ  1472
.krf	     equ  1488	
.crf	     equ  1504
.is3         equ  1520
.ii3         equ  1524
.innerjjnr   equ  1528
.innerk      equ  1532
.salign	     equ  1536							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1540		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm2, [sse_six]
	movups xmm3, [sse_twelve]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .six],  xmm2
	movaps [esp + .twelve], xmm3
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	shufps xmm6, xmm6, 0b
	movaps [esp + .krf], xmm5
	movaps [esp + .crf], xmm6
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 01010101b
	movaps [esp + .c6], xmm0
	movaps [esp + .c12], xmm1

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]
	mulps  xmm0, xmm0
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	; xmm1=rinvsix
	mulps  xmm5, [esp + .rsqOO] ;  xmm5=krsq
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;  xmm6=rinv+krsq
	subps  xmm6, [esp + .crf]
	
	mulps  xmm6, [esp + .qqOO] ;  xmm6=voul=qq*(rinv+krsq-crf)
	mulps xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOO] ; xmm7 = coul part of fscal
	
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=rinvtwelve
	mulps  xmm1, [esp + .c6]	
	mulps  xmm2, [esp + .c12]	
	movaps xmm3, xmm2
	subps  xmm3, xmm1	; xmm3=vnb12-vnb6
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	movaps [esp + .vnbtot], xmm3
	subps  xmm2, xmm1
	addps  xmm2, xmm7
	addps  xmm6, [esp + .vctot] ;  local vctot summation variable
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqOH1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	mulps  xmm0, xmm0
	subps  xmm4, [esp + .crf]
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH1 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqOH2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	mulps xmm0, xmm0
	subps  xmm4, [esp + .crf]
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1O] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	mulps xmm0, xmm0
	subps  xmm4, [esp + .crf]
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1H1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1H2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	mulps xmm0, xmm0
	subps  xmm4, [esp + .crf]
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2O] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2H1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2H2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [esp + .vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]
		
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.

	movaps xmm6, xmm0
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + .krf] ; xmm6=krsq
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	addps   xmm6, xmm3	;  xmm6=rinv+krsq
	mulps   xmm7, [esp + .two]
	subps  xmm6, [esp + .crf]	;  xmm6=rinv+krsq-crf
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	; xmm3=rinv-2*krsq
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	; xmm0=rinvsq
	;; fetch charges to xmm4 (temporary)
	movss   xmm4, [esp + .qqOO]
	movss   xmm1, xmm0
	movhps  xmm4, [esp + .qqOH]
	mulss   xmm1, xmm0

	mulps xmm6, xmm4	;  vcoul 
	mulps xmm3, xmm4	;  coul part of fs
	
	mulss   xmm1, xmm0	;  xmm1(0)=rinvsix
	movaps  xmm2, xmm1	;  zero everything else in xmm2
	mulss   xmm2, xmm2	;  xmm2=rinvtwelve

	mulss   xmm1, [esp + .c6]
	mulss   xmm2, [esp + .c12]
	movaps  xmm4, xmm2
	subss   xmm4, xmm1	;  vnbtot=vnb12-vnb6
	addps   xmm4, [esp + .vnbtot]
	mulss   xmm1, [esp + .six]
	mulss   xmm2, [esp + .twelve]	
	movaps  [esp + .vnbtot], xmm4
	subss   xmm2, xmm1	; fsD+fsR
	addps   xmm2, xmm3	; fsC+fsD+fsR

	addps   xmm6, [esp + .vctot]
	mulps   xmm0, xmm2	;  total fscal
	movaps  [esp + .vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	mulps xmm0, [esp + .krf] ;  krsq
	mulps xmm4, [esp + .krf] ;  krsq 


	;; assemble charges in xmm6
	xorps   xmm6, xmm6
	movss   xmm6, [esp + .qqOH]
	movhps  xmm6, [esp + .qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	; krsq+rinv
	addps   xmm4, xmm7	; krsq+rinv
	subps xmm0, [esp + .crf]
	subps xmm4, [esp + .crf]
	mulps   xmm1, [esp + .two]
	mulps   xmm5, [esp + .two]
	mulps   xmm0, xmm6	;  vcoul
	mulps   xmm4, xmm6	;  vcoul
	addps   xmm4, xmm0		
	addps   xmm4, [esp + .vctot]
	movaps  [esp + .vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;  fscal
	mulps   xmm7, xmm4	;  fscal
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do forces H2 - j water
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	
 
	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1540
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	

proc inl2020_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg	
%$crf		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.vctot       equ   352
.fixO        equ   384
.fiyO        equ   400
.fizO        equ   416
.fixH1       equ   432
.fiyH1       equ   448
.fizH1       equ   464
.fixH2       equ   480
.fiyH2       equ   496
.fizH2       equ   512
.fjx	     equ   528
.fjy         equ   544
.fjz         equ   560
.half        equ   576
.three       equ   592
.two	     equ   608
.krf	     equ   624
.crf	     equ   640
.krsqO       equ   656
.krsqH1      equ   672
.krsqH2	     equ   688	 		
.is3         equ   704
.ii3         equ   708
.innerjjnr   equ   712
.innerk      equ   716
.salign	     equ   720								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 724		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]

	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	shufps xmm6, xmm6, 0b
	movaps [esp + .krf], xmm5
	movaps [esp + .crf], xmm6
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
			
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + .krf]	
	mulps  xmm1, [esp + .krf]	
	mulps  xmm2, [esp + .krf]	

	movaps [esp + .krsqH2], xmm0
	movaps [esp + .krsqH1], xmm1
	movaps [esp + .krsqO], xmm2
	
	; start with rsqO - seed in xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm7, xmm4	; rinvO in xmm7
	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm6, xmm4	; rinvH1 in xmm6
	; rsqH2 - seed in xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  xmm5, xmm4	; rinvH2 in xmm5

	; do O interactions
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	; xmm7=rinv, xmm4=rinvsq

	movaps xmm0, xmm7
	movaps xmm1, [esp + .krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [esp + .crf] ;  xmm0=rinv+krsq-crf
	mulps  xmm1, [esp + .two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + .qqO]
	mulps  xmm7, [esp + .qqO]

	mulps  xmm4, xmm7	; total fsO in xmm4

	addps  xmm0, [esp + .vctot]
	movaps [esp + .vctot], xmm0

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H1 interactions
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	; xmm6=rinv, xmm4=rinvsq
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + .krsqH1]
	addps   xmm6, xmm0	; xmm6=rinv+krsq
	subps   xmm6, [esp + .crf] ;  xmm6=rinv+krsq-crf
	mulps   xmm0, [esp + .two]
	subps   xmm7, xmm0	; xmm7=rinv-2*krsq
	mulps   xmm6, [esp + .qqH] ;  vcoul
	mulps   xmm7, [esp + .qqH]
	mulps  xmm4, xmm7		; total fsH1 in xmm4
	
	addps  xmm6, [esp + .vctot]

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	movaps [esp + .vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; H2 interactions
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	; xmm5=rinv, xmm4=rinvsq
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + .krsqH2]
	addps   xmm5, xmm0	; xmm6=rinv+krsq
	subps   xmm5, [esp + .crf] ;  xmm5=rinv+krsq-crf
	mulps   xmm0, [esp + .two]
	subps   xmm7, xmm0	; xmm7=rinv-2*krsq
	mulps   xmm5, [esp + .qqH] ;  vcoul
	mulps   xmm7, [esp + .qqH]
	mulps  xmm4, xmm7		; total fsH2 in xmm4
	
	addps  xmm5, [esp + .vctot]

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	movaps [esp + .vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	movaps xmm0, xmm4
	mulps xmm0, [esp + .krf]
	movaps [esp + .krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	; xmm4=rinvsq

	movaps xmm1, xmm0	; xmm1=rinv
	movaps xmm3, [esp + .krsqO]
	addps  xmm0, xmm3	; xmm0=rinv+krsq
	subps  xmm0, [esp + .crf] ;  xmm0=rinv+krsq-crf
	mulps  xmm3, [esp + .two]
	subps  xmm1, xmm3	; xmm1=rinv-2*krsq
	mulps  xmm0, [esp + .qqO]	; xmm0=vcoul
	mulps  xmm1, [esp + .qqO] 	; xmm1=coul part of fs

	
	mulps  xmm4, xmm1	; xmm4=total fscal
	addps  xmm0, [esp + .vctot]
	movaps [esp + .vctot], xmm0
	
	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 724
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	
proc inl2030_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$krf		arg
%$crf		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.vctot       equ   768
.fixO        equ   784
.fiyO        equ   800
.fizO        equ   816
.fixH1       equ   832
.fiyH1       equ   848
.fizH1       equ   864
.fixH2       equ   880
.fiyH2       equ   896
.fizH2       equ   912
.fjxO	     equ   928
.fjyO        equ   944
.fjzO        equ   960
.fjxH1	     equ   976
.fjyH1       equ   992
.fjzH1       equ  1008
.fjxH2	     equ  1024
.fjyH2       equ  1040
.fjzH2       equ  1056
.half        equ  1072
.three       equ  1088
.rsqOO       equ  1104
.rsqOH1      equ  1120
.rsqOH2      equ  1136
.rsqH1O      equ  1152
.rsqH1H1     equ  1168
.rsqH1H2     equ  1184
.rsqH2O      equ  1200
.rsqH2H1     equ  1216
.rsqH2H2     equ  1232
.rinvOO      equ  1248
.rinvOH1     equ  1264
.rinvOH2     equ  1280
.rinvH1O     equ  1296
.rinvH1H1    equ  1312
.rinvH1H2    equ  1328
.rinvH2O     equ  1344
.rinvH2H1    equ  1360
.rinvH2H2    equ  1376
.two         equ  1392
.krf	     equ  1408	
.crf	     equ  1424
.is3         equ  1440
.ii3         equ  1444
.innerjjnr   equ  1448
.innerk      equ  1452
.salign	     equ  1456							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1460		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_three]
	movups xmm4, [sse_two]
	movss xmm5, [ebp + %$krf]
	movss xmm6, [ebp + %$crf]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .three], xmm1
	movaps [esp + .two], xmm4
	shufps xmm5, xmm5, 0b
	shufps xmm6, xmm6, 0b
	movaps [esp + .krf], xmm5
	movaps [esp + .crf], xmm6
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5
	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]
	mulps  xmm0, xmm0	;  xmm0=rinvsq

	mulps  xmm5, [esp + .rsqOO] ;  xmm5=krsq
	movaps xmm6, xmm5
	addps  xmm6, xmm7	;  xmm6=rinv+krsq
	subps  xmm6, [esp + .crf]
	mulps  xmm6, [esp + .qqOO] ;  xmm6=voul=qq*(rinv+krsq-crf)
	mulps xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOO] ; xmm7 = coul part of fscal
	
	addps  xmm6, [esp + .vctot] ;  local vctot summation variable
	mulps  xmm0, xmm7
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqOH1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps  xmm0, xmm0
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH1 
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqOH2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1O] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1H1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH1H2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2O] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqOH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqOH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2H1] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm7, xmm0	;  xmm7=rinv
	movaps xmm5, [esp + .krf]	
	movaps xmm1, xmm0
	mulps  xmm5, [esp + .rsqH2H2] ;  xmm5=krsq
	movaps xmm4, xmm5
	addps  xmm4, xmm7	; xmm4=rinv+krsq
	subps  xmm4, [esp + .crf]
	mulps xmm0, xmm0
	mulps  xmm4, [esp + .qqHH] ;  xmm4=voul=qq*(rinv+krsq-crf)
	mulps  xmm5, [esp + .two]
	subps  xmm7, xmm5	; xmm7=rinv-2*krsq
	mulps  xmm7, [esp + .qqHH] ; xmm7 = coul part of fscal
	addps  xmm6, xmm4	; add to local vctot.
	mulps xmm0, xmm7	; fsOH2
	movaps xmm1, xmm0
	movaps xmm2, xmm0

	movaps xmm1, xmm0
	movaps [esp + .vctot], xmm6
	movaps xmm2, xmm0
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]
		
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.

	movaps xmm6, xmm0
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	mulps   xmm6, [esp + .krf] ; xmm6=krsq
	movaps  xmm2, xmm1
	movaps  xmm7, xmm6         ; xmm7=krsq
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water


	
	addps   xmm6, xmm3	;  xmm6=rinv+krsq
	mulps   xmm7, [esp + .two]
	subps   xmm6, [esp + .crf] ; xmm6=rinv+krsq-crf
	
	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	subps   xmm3, xmm7	; xmm3=rinv-2*krsq
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	; xmm0=rinvsq
	;; fetch charges to xmm4 (temporary)
	movss   xmm4, [esp + .qqOO]
	movhps  xmm4, [esp + .qqOH]

	mulps xmm6, xmm4	;  vcoul 
	mulps xmm3, xmm4	;  coul part of fs


	addps   xmm6, [esp + .vctot]
	mulps   xmm0, xmm3	;  total fscal
	movaps  [esp + .vctot], xmm6	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	mulps xmm0, [esp + .krf] ;  krsq
	mulps xmm4, [esp + .krf] ;  krsq 

	;; assemble charges in xmm6
	xorps   xmm6, xmm6
	movss   xmm6, [esp + .qqOH]
	movhps  xmm6, [esp + .qqHH]
	movaps  xmm1, xmm0
	movaps  xmm5, xmm4
	addps   xmm0, xmm3	; krsq+rinv
	addps   xmm4, xmm7	; krsq+rinv
	subps   xmm0, [esp + .crf]
	subps   xmm4, [esp + .crf]
	mulps   xmm1, [esp + .two]
	mulps   xmm5, [esp + .two]
	mulps   xmm0, xmm6	;  vcoul
	mulps   xmm4, xmm6	;  vcoul
	addps   xmm4, xmm0		
	addps   xmm4, [esp + .vctot]
	movaps  [esp + .vctot], xmm4
	movaps  xmm0, xmm3
	movaps  xmm4, xmm7
	mulps   xmm3, xmm3
	mulps   xmm7, xmm7
	subps   xmm0, xmm1
	subps   xmm4, xmm5
	mulps   xmm0, xmm6
	mulps   xmm4, xmm6
	mulps   xmm0, xmm3	;  fscal
	mulps   xmm7, xmm4	;  fscal
	
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do forces H2 - j water
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	
 
	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 	
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1460
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	
		

proc inl3000_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$tabscale      arg
%$VFtab         arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.tabscale    equ   128
.qq          equ   144	
.fs          equ   160
.vctot       equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.innerjjnr   equ   280
.innerk      equ   284
.salign	     equ   288								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 292		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	mulps  xmm3, xmm2

	movaps [esp + .qq], xmm3	
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
		
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov edi, [ebp + %$pos]	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 292
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl3010_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg		
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg
%$tabscale      arg
%$VFtab         arg
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.tabscale    equ   128
.qq          equ   144	
.fs          equ   160
.vctot       equ   176
.fix         equ   192
.fiy         equ   208
.fiz         equ   224
.half        equ   240
.three       equ   256
.is3         equ   272
.ii3         equ   276
.shX	     equ   280
.shY         equ   284
.shZ         equ   288
.ntia	     equ   292	
.innerjjnr0  equ   296
.innerk0     equ   300
.innerjjnr   equ   304
.innerk      equ   308
.salign	     equ   312							
.nscoul      equ   316
.solnr	     equ   320		
	
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 324		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	add   [ebp + %$nsatoms], dword 8

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 
	movss [esp + .shX], xmm0
	movss [esp + .shY], xmm1
	movss [esp + .shZ], xmm2

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	mov   ecx, [eax]
	add   [ebp + %$nsatoms], dword 12
	mov   [esp + .nscoul], ecx	

	; clear vctot 
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	mov   [esp + .solnr], ebx

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]
	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nscoul]
	cmp   ecx, dword 0
	jnz   .mno_coul
	jmp   .last_mno
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner

.unroll_coul_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	mulps  xmm3, xmm2

	movaps [esp + .qq], xmm3	
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
		
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_coul
	jmp   .checksingle_coul
.dopair_coul:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov edi, [ebp + %$pos]	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_coul:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_coul
	jmp    .updateouterdata_coul
.dosingle_coul:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_coul:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .last_mno
	jmp .mno_coul
	
.last_mno:	
	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 324
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl3020_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$tabscale	arg	
%$VFtab		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.rinvO       equ   352
.rinvH1      equ   368
.rinvH2	     equ   384		
.rO          equ   400
.rH1         equ   416
.rH2         equ   432
.tabscale    equ   448	
.two         equ   464
.vctot       equ   480
.fixO        equ   496
.fiyO        equ   512
.fizO        equ   528
.fixH1       equ   544
.fiyH1       equ   560
.fizH1       equ   576
.fixH2       equ   592
.fiyH2       equ   608
.fizH2       equ   624
.fjx	     equ   640
.fjy         equ   656
.fjz         equ   672
.half        equ   688
.three       equ   704
.is3         equ   720
.ii3         equ   724
.innerjjnr   equ   728
.innerk      equ   732
.salign	     equ   736								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 740		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp +%$tabscale]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b 
	movaps [esp + .tabscale], xmm3
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4	

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	; start with rsqO - seed to xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvO], xmm4	; rinvO in xmm4
	mulps   xmm7, xmm4
	movaps  [esp + .rO], xmm7	

	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH1], xmm4	; rinvH1 in xmm4
	mulps   xmm6, xmm4
	movaps  [esp + .rH1], xmm6

	; rsqH2 - seed to xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH2], xmm4	; rinvH2 in xmm4
	mulps   xmm5, xmm4
	movaps  [esp + .rH2], xmm5

	; do O interactions
	;; rO is still in xmm7.
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3

	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd mm0, eax   
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq

        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5 
	xorps  xmm4, xmm4

	mulps  xmm0, [esp + .tabscale]
	mulps  xmm0, [esp + .rinvO]	
	subps  xmm4, xmm0

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4	;  tx in xmm0-xmm2

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	;;  Done with O interactions - now H1!
	movaps xmm7, [esp + .rH1]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm7 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH1]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps xmm4, xmm7

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; Done with H1, finally we do H2 interactions
	movaps xmm7, [esp + .rH2]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH2]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

        movd eax, mm0   
        movd ebx, mm1
        movd ecx, mm2
        movd edx, mm3
	
	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	movaps [esp + .rinvO], xmm0
	
	mulps xmm4, [esp + .tabscale]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm7, mm7
        movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	; xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
	
        movd mm0, eax   
        movd mm1, ecx
        movd mm2, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5

	xorps xmm4, xmm4
	mulps  xmm0, [esp + .tabscale]
	mulps  xmm0, [esp + .rinvO]	
	subps  xmm4, xmm0
		
        movd eax, mm0   
        movd ecx, mm1
        movd edx, mm2	
		
	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ; xmm0-xmm2 now contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 740
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc
	

	
proc inl3030_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$tabscale	arg	
%$VFtab		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.two         equ   768
.tabscale    equ   784
.vctot       equ   800
.fixO        equ   816
.fiyO        equ   832
.fizO        equ   848
.fixH1       equ   864
.fiyH1       equ   880
.fizH1       equ   896
.fixH2       equ   912
.fiyH2       equ   928
.fizH2       equ   944
.fjxO	     equ   960
.fjyO        equ   976
.fjzO        equ   992
.fjxH1	     equ  1008
.fjyH1       equ  1024
.fjzH1       equ  1040
.fjxH2	     equ  1056
.fjyH2       equ  1072
.fjzH2       equ  1088
.half        equ  1104
.three       equ  1120
.rsqOO       equ  1136
.rsqOH1      equ  1152
.rsqOH2      equ  1168
.rsqH1O      equ  1184
.rsqH1H1     equ  1200
.rsqH1H2     equ  1216
.rsqH2O      equ  1232
.rsqH2H1     equ  1248
.rsqH2H2     equ  1264
.rinvOO      equ  1280
.rinvOH1     equ  1296
.rinvOH2     equ  1312
.rinvH1O     equ  1328
.rinvH1H1    equ  1344
.rinvH1H2    equ  1360
.rinvH2O     equ  1376
.rinvH2H1    equ  1392
.rinvH2H2    equ  1408	
.is3         equ  1424
.ii3         equ  1428
.innerjjnr   equ  1432
.innerk      equ  1436
.salign	     equ  1440							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1444		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp +%$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three], xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale],  xmm3

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5		

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
		
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOO] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
		
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2
	
        movd mm0, eax
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.
        ; increment vcoul - then we can get rid of mm5.
        ;; update vctot
        addps  xmm5, [esp + .vctot]
	xorps  xmm2, xmm2
        movaps [esp + .vctot], xmm5
	mulps  xmm3, [esp + .tabscale]
	
	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2
	
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]

	movd eax, mm0
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3
	
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	; xmm1=r
	movaps  xmm0, xmm3	; xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2
	
        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

	mov esi, [ebp + %$VFtab]
	
        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOO]
	movhps  xmm3, [esp + .qqOH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
	
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm2, xmm2
	mulps  xmm3, [esp + .tabscale]

	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.

	;;  start with H1, save H2 data
	movaps [esp + .rsqH2O], xmm4
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	;;  start with H1, save H2 data
	movaps [esp + .rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;  xmm1=r
	movaps xmm0, xmm3	;  xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2

        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do table for H2 - j water interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, [esp + .rsqH2O]
	mulps  xmm1, xmm0	; xmm0=rinv, xmm1=r
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2

        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1444
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc
	



proc inl3100_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale      arg
%$VFtab         arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.six	     equ   128
.twelve	     equ   144
.tabscale    equ   160
.qq          equ   176	
.c6          equ   192
.c12         equ   208
.fs          equ   224
.vctot       equ   240
.vnbtot      equ   256
.fix         equ   272
.fiy         equ   288
.fiz         equ   304
.half        equ   320
.three       equ   336
.is3         equ   352
.ii3         equ   356
.ntia	     equ   360	
.innerjjnr   equ   364
.innerk      equ   368
.salign	     equ   372
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 376		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movups xmm3, [sse_six]
	movups xmm4, [sse_twelve]
	movss xmm5, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three],  xmm2
	movaps [esp + .six],  xmm3
	movaps [esp + .twelve],  xmm4
	shufps xmm5, xmm5, 0b
	movaps [esp + .tabscale], xmm5

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [esp + .qq], xmm3
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]
        mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 376
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl3110_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale      arg
%$VFtab         arg
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.tabscale    equ   128
.qq          equ   144	
.c6          equ   160
.c12         equ   176
.six	     equ   192
.twelve      equ   208
.fs          equ   224
.vctot       equ   240
.vnbtot      equ   256
.fix         equ   272
.fiy         equ   288
.fiz         equ   304
.half        equ   320
.three       equ   336
.is3         equ   352
.ii3         equ   356
.shX	     equ   360
.shY         equ   364
.shZ         equ   368
.ntia	     equ   372	
.innerjjnr0  equ   376
.innerk0     equ   380	
.innerjjnr   equ   384
.innerk      equ   388
.salign	     equ   392							
.nsvdwc      equ   396
.nscoul      equ   400
.nsvdw       equ   404
.solnr	     equ   408		
	push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 412		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movups xmm3, [sse_six]
	movups xmm4, [sse_twelve]
	movss xmm5, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three], xmm2
	movaps [esp + .six],  xmm3
	movaps [esp + .twelve], xmm4
	shufps xmm5, xmm5, 0b
	movaps [esp + .tabscale], xmm5

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movlps xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 8] 
	movlps [esp + .shX], xmm0
	movss [esp + .shZ], xmm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	add   [ebp + %$nsatoms], dword 12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + .nsvdwc], edx
	mov   [esp + .nscoul], eax
	mov   [esp + .nsvdw], ecx
		
	;; clear potential
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	mov   [esp + .solnr],  ebx

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]
	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]
	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [esp + .qq], xmm3
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdwc
	jmp   .checksingle_vdwc
.dopair_vdwc:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	mov    edi, [ebp + %$faction]
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdwc:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdwc
	jmp    .updateouterdata_vdwc
.dosingle_vdwc:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
						
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; L-J
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	; xmm4=rinvsq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [esp + .vctot], xmm5 

	mulps  xmm6, xmm4	; xmm6=rinvsix
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm6, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm7, [esp + .vnbtot]
	addps  xmm7, xmm4
	mulps  xmm4, [esp + .twelve]
	subps  xmm7, xmm6
	mulps  xmm3, [esp + .tabscale]
	mulps  xmm6, [esp + .six]
	movaps [esp + .vnbtot], xmm7
	subps  xmm4, xmm6
	mulps  xmm4, xmm0
	subps  xmm4, xmm3
	mulps  xmm4, xmm0

	mov edi, [ebp +%$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdwc:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5


	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul:
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner

.unroll_coul_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	mulps  xmm3, xmm2

	movaps [esp + .qq], xmm3	
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_coul
	jmp   .checksingle_coul
.dopair_coul:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov edi, [ebp + %$pos]	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_coul:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_coul
	jmp    .updateouterdata_coul
.dosingle_coul:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_coul:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw:
	mov  ecx, [esp + .nsvdw]
	cmp  ecx, byte 0
	jnz  .mno_vdw
	jmp  .last_mno
.mno_vdw:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdw
	jmp   .checksingle_vdw
.dopair_vdw:	
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdw:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdw
	jmp    .updateouterdata_vdw
.dosingle_vdw:
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rcpps xmm5, xmm4
	; 1/x lookup seed in xmm5
	movaps xmm0, [esp + .two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;  xmm0=rinvsq
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve

	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;  vnb=vnb12-vnb6
	addps  xmm5, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	; xmm4=total fscal
	
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movaps [esp + .vnbtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	mov edi, [ebp +%$faction]

	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5
.updateouterdata_vdw:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5
	
	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
.last_mno:	
	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 412
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl3120_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
%$tabscale	arg	
%$VFtab		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.rinvO       equ   352
.rinvH1      equ   368
.rinvH2	     equ   384		
.rO          equ   400
.rH1         equ   416
.rH2         equ   432
.tabscale    equ   448	
.two         equ   464
.c6          equ   480
.c12         equ   496
.six         equ   512
.twelve      equ   528
.vctot       equ   544
.vnbtot      equ   560
.fixO        equ   576
.fiyO        equ   592
.fizO        equ   608
.fixH1       equ   624
.fiyH1       equ   640
.fizH1       equ   656
.fixH2       equ   672
.fiyH2       equ   688
.fizH2       equ   704
.fjx	     equ   720
.fjy         equ   736
.fjz         equ   752
.half        equ   768
.three       equ   784
.is3         equ   800
.ii3         equ   804
.ntia	     equ   808	
.innerjjnr   equ   812
.innerk      equ   816
.salign	     equ   820								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 824		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movups xmm3, [sse_six]
	movups xmm4, [sse_twelve]
	movss xmm5, [ebp +%$tabscale]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three],  xmm2
	movaps [esp + .six],  xmm3
	movaps [esp + .twelve],  xmm4
	shufps xmm5, xmm5, 0b 
	movaps [esp + .tabscale], xmm5
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	; start with rsqO - seed to xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvO], xmm4	; rinvO in xmm4
	mulps   xmm7, xmm4
	movaps  [esp + .rO], xmm7	

	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH1], xmm4	; rinvH1 in xmm4
	mulps   xmm6, xmm4
	movaps  [esp + .rH1], xmm6

	; rsqH2 - seed to xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH2], xmm4	; rinvH2 in xmm4
	mulps   xmm5, xmm4
	movaps  [esp + .rH2], xmm5

	; do O interactions
	;; rO is still in xmm7.
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3

	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd mm0, eax   
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq

	; do nontable L-J
	movaps xmm2, [esp + .rinvO]
	mulps  xmm2, xmm2

        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5 

	movaps xmm1, xmm2
	mulps  xmm1, xmm1
	mulps  xmm1, xmm2	; xmm1=rinvsix
	movaps xmm4, xmm1
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm1, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm3, xmm4
	subps  xmm3, xmm1	; xmm3=vnb12-vnb6
	mulps  xmm1, [esp + .six]
	mulps  xmm4, [esp + .twelve]
	subps  xmm4, xmm1
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm4, [esp + .rinvO]
	mulps  xmm0, [esp + .tabscale]
	subps  xmm4, xmm0
	movaps [esp + .vnbtot], xmm3
	mulps  xmm4, [esp + .rinvO]	

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4	;  tx in xmm0-xmm2

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	;;  Done with O interactions - now H1!
	movaps xmm7, [esp + .rH1]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm7 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH1]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps xmm4, xmm7

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; Done with H1, finally we do H2 interactions
	movaps xmm7, [esp + .rH2]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH2]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

        movd eax, mm0   
        movd ebx, mm1
        movd ecx, mm2
        movd edx, mm3
	
	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	xorps xmm6, xmm6
	mov esi, [ebp + %$type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + %$nbfp]
	shl ebx, 1	
	add ebx, [esp + .ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 11111100b
	shufps xmm7, xmm7, 11111101b
	movaps [esp + .c6], xmm6
	movaps [esp + .c12], xmm7

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	movaps [esp + .rinvO], xmm0
	
	mulps xmm4, [esp + .tabscale]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm7, mm7
        movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	; xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
	
        movd mm0, eax   
        movd mm1, ecx
        movd mm2, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5

	; do nontable L-J
	movaps xmm2, [esp + .rinvO]
	mulps  xmm2, xmm2
	movaps xmm1, xmm2
	mulps  xmm1, xmm1
	mulps  xmm1, xmm2	; xmm1=rinvsix
	movaps xmm4, xmm1
	mulps  xmm4, xmm4	; xmm4=rinvtwelve
	mulps  xmm1, [esp + .c6]
	mulps  xmm4, [esp + .c12]
	movaps xmm3, xmm4
	subps  xmm3, xmm1	; xmm3=vnb12-vnb6
	mulps  xmm1, [esp + .six]
	mulps  xmm4, [esp + .twelve]
	subps  xmm4, xmm1
	addps  xmm3, [esp + .vnbtot]
	mulps  xmm4, [esp + .rinvO]
	mulps  xmm0, [esp + .tabscale]
	subps  xmm4, xmm0
	movaps [esp + .vnbtot], xmm3
	mulps  xmm4, [esp + .rinvO]	
		
        movd eax, mm0   
        movd ecx, mm1
        movd edx, mm2	
		
	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ; xmm0-xmm2 now contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 824
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc
	

	
proc inl3130_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale	arg	
%$VFtab		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.two         equ   768
.tabscale    equ   784
.c6          equ   800
.c12	     equ   816		 
.six         equ   832
.twelve      equ   848		 
.vctot       equ   864
.vnbtot      equ   880
.fixO        equ   896
.fiyO        equ   912
.fizO        equ   928
.fixH1       equ   944
.fiyH1       equ   960
.fizH1       equ   976
.fixH2       equ   992
.fiyH2       equ  1008
.fizH2       equ  1024
.fjxO	     equ  1040
.fjyO        equ  1056
.fjzO        equ  1072
.fjxH1	     equ  1088
.fjyH1       equ  1104
.fjzH1       equ  1120
.fjxH2	     equ  1136
.fjyH2       equ  1152
.fjzH2       equ  1168
.half        equ  1184
.three       equ  1200
.rsqOO       equ  1216
.rsqOH1      equ  1232
.rsqOH2      equ  1248
.rsqH1O      equ  1264
.rsqH1H1     equ  1280
.rsqH1H2     equ  1296
.rsqH2O      equ  1312
.rsqH2H1     equ  1328
.rsqH2H2     equ  1344
.rinvOO      equ  1360
.rinvOH1     equ  1376
.rinvOH2     equ  1392
.rinvH1O     equ  1408
.rinvH1H1    equ  1424
.rinvH1H2    equ  1440
.rinvH2O     equ  1456
.rinvH2H1    equ  1472
.rinvH2H2    equ  1488
.fstmp	     equ  1504	
.is3         equ  1520
.ii3         equ  1524
.innerjjnr   equ  1528
.innerk      equ  1532
.salign	     equ  1536							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1540		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movups xmm3, [sse_six]
	movups xmm4, [sse_twelve]
	movss xmm5, [ebp +%$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three], xmm2
	movaps [esp + .six], xmm3
	movaps [esp + .twelve], xmm4
	shufps xmm5, xmm5, 0b
	movaps [esp + .tabscale],  xmm5

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 01010101b
	movaps [esp + .c6], xmm0
	movaps [esp + .c12], xmm1

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
		
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOO] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
		
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2
	
        movd mm0, eax
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.
        ; increment vcoul - then we can get rid of mm5.
        ;; update vctot
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	mulps  xmm3, [esp + .tabscale]
	
	;;  start doing lj
	movaps xmm2, xmm0
	mulps  xmm2, xmm2
	movaps xmm1, xmm2
	mulps  xmm1, xmm2
	mulps  xmm1, xmm2	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;  xmm2=rinvtwelve
	mulps  xmm1, [esp + .c6]
	mulps  xmm2, [esp + .c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [esp + .vnbtot]
	mulps  xmm1, [esp + .six]
	mulps  xmm2, [esp + .twelve]
	movaps [esp + .vnbtot], xmm4
	subps  xmm2, xmm1
	mulps  xmm2, xmm0

	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2
	
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2

	pslld   mm6, 2
	pslld   mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]

	movd eax, mm0
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3
	
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	; xmm1=r
	movaps  xmm0, xmm3	; xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2
	
        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

	mov esi, [ebp + %$VFtab]
	
        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOO]
	movhps  xmm3, [esp + .qqOH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
	
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5

	mulps  xmm3, [esp + .tabscale]
	
	;;  start doing lj
	xorps  xmm2, xmm2
	movss  xmm2, xmm0
	mulss  xmm2, xmm2
	movaps xmm1, xmm2
	mulss  xmm1, xmm2
	mulss  xmm1, xmm2	;  xmm1=rinvsix
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;  xmm2=rinvtwelve
	mulss  xmm1, [esp + .c6]
	mulss  xmm2, [esp + .c12]
	movaps xmm4, xmm2
	subss  xmm4, xmm1
	addps  xmm4, [esp + .vnbtot]
	mulss  xmm1, [esp + .six]
	mulss  xmm2, [esp + .twelve]
	movaps [esp + .vnbtot], xmm4
	subss  xmm2, xmm1
	mulss  xmm2, xmm0

	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.

	;;  start with H1, save H2 data
	movaps [esp + .rsqH2O], xmm4
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	;;  start with H1, save H2 data
	movaps [esp + .rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;  xmm1=r
	movaps xmm0, xmm3	;  xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2

        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do table for H2 - j water interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, [esp + .rsqH2O]
	mulps  xmm1, xmm0	; xmm0=rinv, xmm1=r
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
	pslld   mm6, 2
	pslld   mm7, 2

        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1540
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	


proc inl3300_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale      arg
%$VFtab         arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.tabscale    equ   128
.qq          equ   144	
.c6          equ   160
.c12         equ   176
.fs          equ   192
.vctot       equ   208
.vnbtot      equ   224
.fix         equ   240
.fiy         equ   256
.fiz         equ   272
.half        equ   288
.three       equ   304
.is3         equ   320
.ii3         equ   324
.ntia	     equ   328	
.innerjjnr   equ   332
.innerk      equ   336
.salign	     equ   340								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 344		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [esp + .qq], xmm3
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]
		
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 32]
	movlps xmm7, [esi + ecx*4 + 32]
	movhps xmm5, [esi + ebx*4 + 32]
	movhps xmm7, [esi + edx*4 + 32] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 40]
	movlps xmm3, [esi + ecx*4 + 40]
	movhps xmm7, [esi + ebx*4 + 40]
	movhps xmm3, [esi + edx*4 + 40] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair
	jmp   .checksingle
.dopair:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6
	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 32]
	movhps xmm5, [esi + edx*4 + 32] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + ecx*4 + 40]
	movhps xmm7, [esi + edx*4 + 40] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle
	jmp    .updateouterdata
.dosingle:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	lea  ebx, [ebx + ebx*2]
						
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 32]
	movlps xmm6, [esi + ebx*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 344
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc





proc inl3310_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale      arg
%$VFtab         arg
%$nsatoms       arg			
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ix	     equ     0
.iy	     equ    16
.iz          equ    32
.iq          equ    48
.dx          equ    64
.dy          equ    80
.dz          equ    96
.two	     equ   112
.tabscale    equ   128
.qq          equ   144	
.c6          equ   160
.c12         equ   176
.fs          equ   192
.vctot       equ   208
.vnbtot      equ   224
.fix         equ   240
.fiy         equ   256
.fiz         equ   272
.half        equ   288
.three       equ   304
.is3         equ   320
.ii3         equ   324
.shX	     equ   328
.shY         equ   332
.shZ         equ   336
.ntia	     equ   340	
.innerjjnr0  equ   344
.innerk0     equ   348	
.innerjjnr   equ   352
.innerk      equ   356
.salign	     equ   360							
.nsvdwc      equ   364
.nscoul      equ   368
.nsvdw       equ   372
.solnr	     equ   376		
	push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 380		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp + %$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two], xmm1
	movaps [esp + .three], xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale], xmm3

	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movlps xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 8] 
	movlps [esp + .shX], xmm0
	movss [esp + .shZ], xmm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	add   [ebp + %$nsatoms], dword 12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + .nsvdwc], edx
	mov   [esp + .nscoul], eax
	mov   [esp + .nsvdw], ecx
		
	;; clear potential
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	mov   [esp + .solnr],  ebx

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr0], eax     ;  pointer to jjnr[nj0]
	mov   [esp + .innerk0], edx        ;  number of innerloop atoms

	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]
	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4
	
	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [esp + .qq], xmm3
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]
		
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 32]
	movlps xmm7, [esi + ecx*4 + 32]
	movhps xmm5, [esi + ebx*4 + 32]
	movhps xmm7, [esi + edx*4 + 32] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 40]
	movlps xmm3, [esi + ecx*4 + 40]
	movhps xmm7, [esi + ebx*4 + 40]
	movhps xmm3, [esi + edx*4 + 40] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdwc
	jmp   .checksingle_vdwc
.dopair_vdwc:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6
	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 32]
	movhps xmm5, [esi + edx*4 + 32] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 

	movlps xmm7, [esi + ecx*4 + 40]
	movhps xmm7, [esi + edx*4 + 40] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdwc:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdwc
	jmp    .updateouterdata_vdwc
.dosingle_vdwc:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	lea  ebx, [ebx + ebx*2]
						
	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	; put scalar force on stack temporarily...
	movaps [esp + .fs], xmm3

	; dispersion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6
	addps  xmm7, [esp + .fs] ; add to fscal

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 32]
	movlps xmm6, [esi + ebx*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdwc:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5


	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul:
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [ebp + %$facel]
	shufps xmm3, xmm3, 0b
	
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + .iq], xmm3
	
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   [esp + .ii3], ebx
	
	; clear  i forces
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner

.unroll_coul_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + .iq]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	mulps  xmm3, xmm2

	movaps [esp + .qq], xmm3	
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]
		
	movlps xmm5, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm5, [esi + ebx*4]
	movhps xmm7, [esi + edx*4] ; got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; coulomb table ready, in xmm4-xmm7 	
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_coul
	jmp   .checksingle_coul
.dopair_coul:	
	mov esi, [ebp + %$charge]

        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 00000000b 
	shufps xmm3, xmm3, 00001000b ; xmm3(0,1) has the charges.

	mulps  xmm3, [esp + .iq]
	movlhps xmm3, xmm7
	movaps [esp + .qq], xmm3

	mov edi, [ebp + %$pos]	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	lea   ecx, [ecx + ecx*2]
	lea   edx, [edx + edx*2]

	movlps xmm5, [esi + ecx*4]
	movhps xmm5, [esi + edx*4] ; got half coulomb table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b 
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps  xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_coul:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_coul
	jmp    .updateouterdata_coul
.dosingle_coul:
	mov esi, [ebp + %$charge]
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6
	movss xmm6, [esi + eax*4]	; xmm6(0) has the charge	
	mulps  xmm6, [esp + .iq]
	movaps [esp + .qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6
	
	lea   ebx, [ebx + ebx*2]

	movlps xmm4, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7

	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	movaps xmm3, [esp + .qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
	mulps  xmm5, xmm3 ; vcoul=qq*VV
	mulps  xmm3, xmm7 ; fijC=FF*qq
	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	addps  xmm5, [esp + .vctot]
	movaps [esp + .vctot], xmm5 

	xorps xmm4, xmm4

	mulps xmm3, [esp + .tabscale]
	mulps xmm3, xmm0
	subps  xmm4, xmm3
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_coul:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw:
	mov  ecx, [esp + .nsvdw]
	cmp  ecx, byte 0
	jnz  .mno_vdw
	jmp  .last_mno
.mno_vdw:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]

        mov   edx, [ebp + %$type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [ebp + %$ntype]
        shl   edx, 1
        mov   [esp + .ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	movss xmm0, [esp + .shX]
	movss xmm1, [esp + .shY]
	movss xmm2, [esp + .shZ]

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]
	
	xorps xmm4, xmm4
	movaps [esp + .fix], xmm4
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm4

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b

	movaps [esp + .ix], xmm0
	movaps [esp + .iy], xmm1
	movaps [esp + .iz], xmm2

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 4
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6
	
	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ix-iz to xmm4-xmm6
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	; mm6/mm7 contain lu indices
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  esi, [ebp + %$VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   eax, [eax + eax*2] 
	lea   ebx, [ebx + ebx*2] 
	lea   ecx, [ecx + ecx*2] 
	lea   edx, [edx + edx*2] 

	; dispersion
	movlps xmm5, [esi + eax*4 + 0]
	movlps xmm7, [esi + ecx*4 + 0]
	movhps xmm5, [esi + ebx*4 + 0]
	movhps xmm7, [esi + edx*4 + 0] ; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 
	
	movlps xmm7, [esi + eax*4 + 8]
	movlps xmm3, [esi + ecx*4 + 8]
	movhps xmm7, [esi + ebx*4 + 8]
	movhps xmm3, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + eax*4 + 16]
	movlps xmm7, [esi + ecx*4 + 16]
	movhps xmm5, [esi + ebx*4 + 16]
	movhps xmm7, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + eax*4 + 24]
	movlps xmm3, [esi + ecx*4 + 24]
	movhps xmm7, [esi + ebx*4 + 24]
	movhps xmm3, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    edi, [ebp + %$faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; the fj's - start by accumulating x & y forces from memory
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 10001000b
	shufps xmm4, xmm6, 11011101b			      

	; now xmm3-xmm5 contains fjx, fjy, fjz
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	; unpack them back so we can store them - first x & y in xmm3/xmm4

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	; xmm6(l)=x & y for j1, (h) for j2
	; xmm3(l)=x & y for j3, (h) for j4
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;;  and the z forces
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:
	;;  check if at least two particles remain
	add   [esp + .innerk], dword 4
	mov   edx, [esp + .innerk]
	and   edx, 10b
	jnz   .dopair_vdw
	jmp   .checksingle_vdw
.dopair_vdw:	
        mov   ecx, [esp + .innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   [esp + .innerjjnr], dword 8	
	xorps xmm7, xmm7

	mov esi, [ebp + %$type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + %$pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 1000b 	
	shufps xmm6, xmm6, 1101b
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	; move coordinates to xmm0-xmm2
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0b
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 10001000b
	
	shufps xmm0, xmm0, 10001000b
	shufps xmm1, xmm1, 11011101b
			
	mov    edi, [ebp + %$faction]
	; move ix-iz to xmm4-xmm6
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	lea   ecx, [ecx + ecx*2] 
	lea   edx, [edx + edx*2] 

	; dispersion
	movlps xmm5, [esi + ecx*4 + 0]
	movhps xmm5, [esi + edx*4 + 0]; got half dispersion table
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 10001000b
	shufps xmm5, xmm5, 11011101b 
	
	movlps xmm7, [esi + ecx*4 + 8]
	movhps xmm7, [esi + edx*4 + 8] ; other half of dispersion table
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 10001000b
	shufps xmm7, xmm7, 11011101b
	; dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm5, [esi + ecx*4 + 16]
	movhps xmm5, [esi + edx*4 + 16] ; got half repulsion table
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 10001000b
	shufps xmm5, xmm7, 11011101b 

	movlps xmm7, [esi + ecx*4 + 24]
	movhps xmm7, [esi + edx*4 + 24] ; other half of repulsion table
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 10001000b
	shufps xmm7, xmm3, 11011101b 
	; table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update the fj's
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 11100001b
	shufps  xmm1, xmm1, 11100001b
	shufps  xmm2, xmm2, 11100001b

	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	

.checksingle_vdw:				
	mov   edx, [esp + .innerk]
	and   edx, 1b
	jnz    .dosingle_vdw
	jmp    .updateouterdata_vdw
.dosingle_vdw:
	mov edi, [ebp + %$pos]
	mov   ecx, [esp + .innerjjnr]
	mov   eax, [ecx]	
	xorps  xmm6, xmm6

	mov esi, [ebp + %$type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + %$nbfp]
	shl ecx, 1
	add ecx, [esp + .ntia]
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 11111100b	
	shufps xmm6, xmm6, 11111101b	
			
	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	; move coordinates to xmm0-xmm2
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	 
	
	movaps xmm4, [esp + .ix]
	movaps xmm5, [esp + .iy]
	movaps xmm6, [esp + .iz]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dx], xmm4
	movaps [esp + .dy], xmm5
	movaps [esp + .dz], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	; rsq in xmm4

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv

	mulps xmm4, xmm0	; xmm4=r
	mulps xmm4, [esp + .tabscale]

	cvttps2pi mm6, xmm4     ; mm6 contain lu indices
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;xmm1=eps
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;xmm2=eps2

	pslld mm6, 2

	mov  esi, [ebp + %$VFtab]
	movd ebx, mm6

	lea   ebx, [ebx + ebx*2] 	

	; dispersion
	movlps xmm4, [esi + ebx*4 + 0]
	movlps xmm6, [esi + ebx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV

	movaps xmm4, [esp + .c6]
	mulps  xmm7, xmm4	 ; fijD
	mulps  xmm5, xmm4	 ; vnb6

	; put scalar force on stack. Update vnbtot directly.
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .fs], xmm7
	movaps [esp + .vnbtot], xmm5

	; repulsion
	movlps xmm4, [esi + ebx*4 + 16]
	movlps xmm6, [esi + ebx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1b
	shufps xmm7, xmm7, 1b
	; table ready in xmm4-xmm7
	
	mulps  xmm6, xmm1	; xmm6=Geps
	mulps  xmm7, xmm2	; xmm7=Heps2
	addps  xmm5, xmm6
	addps  xmm5, xmm7	; xmm5=Fp	
	mulps  xmm7, [esp + .two]	; two*Heps2
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ; xmm7=FF
	mulps  xmm5, xmm1 ; xmm5=eps*Fp
	addps  xmm5, xmm4 ; xmm5=VV
 	
	movaps xmm4, [esp + .c12]
	mulps  xmm7, xmm4 ; fijR
	mulps  xmm5, xmm4 ; vnb12
	addps  xmm7, [esp + .fs] 
	
	addps  xmm5, [esp + .vnbtot]
	movaps [esp + .vnbtot], xmm5
	xorps  xmm4, xmm4

	mulps xmm7, [esp + .tabscale]
	mulps xmm7, xmm0
	subps  xmm4, xmm7
	mov    edi, [ebp + %$faction]

	movaps xmm0, [esp + .dx]
	movaps xmm1, [esp + .dy]
	movaps xmm2, [esp + .dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	; xmm0-xmm2 contains tx-tz (partial force)
	; now update f_i 
	movaps xmm3, [esp + .fix]
	movaps xmm4, [esp + .fiy]
	movaps xmm5, [esp + .fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + .fix], xmm3
	movaps [esp + .fiy], xmm4
	movaps [esp + .fiz], xmm5
	; update fj
	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.updateouterdata_vdw:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fix]
	movaps xmm1, [esp + .fiy]
	movaps xmm2, [esp + .fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; increment fshift force
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5
	
	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
.last_mno:	
	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 380
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl3320_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg	
%$tabscale	arg	
%$VFtab		arg	
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.iqO         equ   144 
.iqH         equ   160 
.dxO         equ   176
.dyO         equ   192
.dzO         equ   208	
.dxH1        equ   224
.dyH1        equ   240
.dzH1        equ   256	
.dxH2        equ   272
.dyH2        equ   288
.dzH2        equ   304	
.qqO         equ   320
.qqH         equ   336
.rinvO       equ   352
.rinvH1      equ   368
.rinvH2	     equ   384		
.rO          equ   400
.rH1         equ   416
.rH2         equ   432
.tabscale    equ   448	
.two         equ   464
.c6          equ   480
.c12         equ   496
.vctot       equ   512
.vnbtot      equ   528
.fixO        equ   544
.fiyO        equ   560
.fizO        equ   576
.fixH1       equ   592
.fiyH1       equ   608
.fizH1       equ   624
.fixH2       equ   640
.fiyH2       equ   656
.fizH2       equ   672
.fjx	     equ   688
.fjy         equ   704
.fjz         equ   720
.half        equ   736
.three       equ   752
.is3         equ   768
.ii3         equ   772
.ntia	     equ   776	
.innerjjnr   equ   780
.innerk      equ   784
.salign	     equ   788								
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 792		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp +%$tabscale]
	
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three],  xmm2
	shufps xmm3, xmm3, 0b 
	movaps [esp + .tabscale], xmm3
	
	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	movss xmm5, [ebp + %$facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	movaps [esp + .iqO], xmm3
	movaps [esp + .iqH], xmm4
	
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5
	
	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .odd_inner
.unroll_loop:
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 

	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$charge]        ; base of charge[]
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 00000000b 
	shufps xmm4, xmm7, 00000000b 
	shufps xmm3, xmm4, 10001000b ;  all charges in xmm3
	movaps xmm4, xmm3	     ;  and in xmm4
	mulps  xmm3, [esp + .iqO]
	mulps  xmm4, [esp + .iqH]

	movd  mm0, eax		;  use mmx registers as temp. storage
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + .qqO], xmm3
	movaps  [esp + .qqH], xmm4
	
	mov esi, [ebp + %$type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + %$nbfp]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + .ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 10001000b
	shufps xmm6, xmm7, 11011101b
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + .c6], xmm4
	movaps [esp + .c12], xmm6

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	

	; move four coordinates to xmm0-xmm2	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0b
	shufps xmm6, xmm1, 0b
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 10001000b
	
	shufps xmm0, xmm5, 10001000b
	shufps xmm1, xmm5, 11011101b		

	; move ixO-izO to xmm4-xmm6
	movaps xmm4, [esp + .ixO]
	movaps xmm5, [esp + .iyO]
	movaps xmm6, [esp + .izO]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxO], xmm4
	movaps [esp + .dyO], xmm5
	movaps [esp + .dzO], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	; rsqO in xmm7

	; move ixH1-izH1 to xmm4-xmm6
	movaps xmm4, [esp + .ixH1]
	movaps xmm5, [esp + .iyH1]
	movaps xmm6, [esp + .izH1]

	; calc dr
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	; store dr
	movaps [esp + .dxH1], xmm4
	movaps [esp + .dyH1], xmm5
	movaps [esp + .dzH1], xmm6
	; square it
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	; rsqH1 in xmm6

	; move ixH2-izH2 to xmm3-xmm5
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]

	; calc dr
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	; store dr
	movaps [esp + .dxH2], xmm3
	movaps [esp + .dyH2], xmm4
	movaps [esp + .dzH2], xmm5
	; square it
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	; rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7

	; start with rsqO - seed to xmm2	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm7	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvO], xmm4	; rinvO in xmm4
	mulps   xmm7, xmm4
	movaps  [esp + .rO], xmm7	

	; rsqH1 - seed in xmm2	
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm6	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH1], xmm4	; rinvH1 in xmm4
	mulps   xmm6, xmm4
	movaps  [esp + .rH1], xmm6

	; rsqH2 - seed to xmm2	
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + .three]
	mulps   xmm2, xmm5	; rsq*lu*lu
	subps   xmm4, xmm2	; 3.0-rsq*lu*lu
	mulps   xmm4, xmm3	; lu*(3-rsq*lu*lu)
	mulps   xmm4, [esp + .half]
	movaps  [esp + .rinvH2], xmm4	; rinvH2 in xmm4
	mulps   xmm5, xmm4
	movaps  [esp + .rH2], xmm5

	; do O interactions
	;; rO is still in xmm7.
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3

	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd mm0, eax   
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5 

        ; dispersion
        movlps xmm5, [esi + eax*4 + 16]
        movlps xmm7, [esi + ecx*4 + 16]
        movhps xmm5, [esi + ebx*4 + 16]
        movhps xmm7, [esi + edx*4 + 16] ; got half dispersion table
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 
        
        movlps xmm7, [esi + eax*4 + 24]
        movlps xmm3, [esi + ecx*4 + 24]
        movhps xmm7, [esi + ebx*4 + 24]
        movhps xmm3, [esi + edx*4 + 24] ; other half of dispersion table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; dispersion table ready, in xmm4-xmm7  
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV

        movaps xmm4, [esp + .c6]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb6
        addps  xmm0, xmm7 ; add to fscal

        ; Update vnbtot directly.
        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .vnbtot], xmm5

        ; repulsion
        movlps xmm5, [esi + eax*4 + 32]
        movlps xmm7, [esi + ecx*4 + 32]
        movhps xmm5, [esi + ebx*4 + 32]
        movhps xmm7, [esi + edx*4 + 32] ; got half repulsion table
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 40]
        movlps xmm3, [esi + ecx*4 + 40]
        movhps xmm7, [esi + ebx*4 + 40]
        movhps xmm3, [esi + edx*4 + 40] ; other half of repulsion table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; repulsion table ready, in xmm4-xmm7  	
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV

        movaps xmm4, [esp + .c12]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb12
        addps  xmm7, xmm0 ; add to fscal
        addps  xmm5, [esp + .vnbtot] ;  total nonbonded potential in xmm5.
	xorps xmm4, xmm4
	
	mulps  xmm7, [esp + .rinvO] ;  total fscal now in xmm7

	mulps  xmm7, [esp + .tabscale]
        movaps [esp + .vnbtot], xmm5
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4	;  tx in xmm0-xmm2

	; update O forces
	movaps xmm3, [esp + .fixO]
	movaps xmm4, [esp + .fiyO]
	movaps xmm7, [esp + .fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixO], xmm3
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm7
	; update j forces with water O
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	;;  Done with O interactions - now H1!
	movaps xmm7, [esp + .rH1]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm7 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH1]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps xmm4, xmm7

	movaps xmm0, [esp + .dxH1]
	movaps xmm1, [esp + .dyH1]
	movaps xmm2, [esp + .dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	; update H1 forces
	movaps xmm3, [esp + .fixH1]
	movaps xmm4, [esp + .fiyH1]
	movaps xmm7, [esp + .fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH1], xmm3
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm7
	; update j forces with water H1
	addps  xmm0, [esp + .fjx]
	addps  xmm1, [esp + .fjy]
	addps  xmm2, [esp + .fjz]
	movaps [esp + .fjx], xmm0
	movaps [esp + .fjy], xmm1
	movaps [esp + .fjz], xmm2

	; Done with H1, finally we do H2 interactions
	movaps xmm7, [esp + .rH2]
	mulps   xmm7, [esp + .tabscale]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ; mm6/mm7 contain lu indices
	
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm4, mm7
        movlhps xmm3, xmm4
	
        subps xmm7, xmm3
	movaps xmm1, xmm7	;  xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
		
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm7, xmm0 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul
	xorps  xmm4, xmm4
        addps  xmm5, [esp + .vctot]
	mulps  xmm7, [esp + .rinvH2]
        movaps [esp + .vctot], xmm5 
	mulps  xmm7, [esp + .tabscale]
	subps  xmm4, xmm7

	movaps xmm0, [esp + .dxH2]
	movaps xmm1, [esp + .dyH2]
	movaps xmm2, [esp + .dzH2]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

        movd eax, mm0   
        movd ebx, mm1
        movd ecx, mm2
        movd edx, mm3
	
	; update H2 forces
	movaps xmm3, [esp + .fixH2]
	movaps xmm4, [esp + .fiyH2]
	movaps xmm7, [esp + .fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + .fixH2], xmm3
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm7

	mov edi, [ebp +%$faction]
	; update j forces
	addps xmm0, [esp + .fjx]
	addps xmm1, [esp + .fjy]
	addps xmm2, [esp + .fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 10001000b
	shufps xmm4, xmm7, 11011101b			      
	; xmm3 has fjx, xmm4 has fjy.
	subps xmm3, xmm0
	subps xmm4, xmm1
	; unpack the back for storing.
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	; finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 11100101b
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 11101010b
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 11111111b
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .odd_inner
	jmp   .unroll_loop
.odd_inner:	
	add   [esp + .innerk], dword 4
	jnz   .odd_loop
	jmp   .updateouterdata
.odd_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + .iqO]
	mov esi, [ebp + %$charge] 
	movhps xmm4, [esp + .iqH]     
	movss xmm3, [esi + eax*4]	; charge in xmm3
	shufps xmm3, xmm3, 0b
	mulps xmm3, xmm4
	movaps [esp + .qqO], xmm3	; use oxygen qq for storage.

	xorps xmm6, xmm6
	mov esi, [ebp + %$type]
	mov ebx, [esi + eax*4]
	mov esi, [ebp + %$nbfp]
	shl ebx, 1	
	add ebx, [esp + .ntia]
	movlps xmm6, [esi + ebx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 11111100b
	shufps xmm7, xmm7, 11111101b
	movaps [esp + .c6], xmm6
	movaps [esp + .c12], xmm7

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  
	
	; move j coords to xmm0-xmm2
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	
	movss xmm3, [esp + .ixO]
	movss xmm4, [esp + .iyO]
	movss xmm5, [esp + .izO]
		
	movlps xmm6, [esp + .ixH1]
	movlps xmm7, [esp + .ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + .iyH1]
	movlps xmm7, [esp + .iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + .izH1]
	movlps xmm7, [esp + .izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + .dxO], xmm3
	movaps [esp + .dyO], xmm4
	movaps [esp + .dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	; rsq in xmm4.

	rsqrtps xmm5, xmm4
	; lookup seed in xmm5
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + .three]
	mulps xmm5, xmm4	;rsq*lu*lu			
	movaps xmm0, [esp + .half]
	subps xmm1, xmm5	; 3.0-rsq*lu*lu
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	; xmm0=rinv
	mulps xmm4, xmm0	; xmm4=r
	movaps [esp + .rinvO], xmm0
	
	mulps xmm4, [esp + .tabscale]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm7, mm7
        movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	; xmm1=eps
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
	
        movd mm0, eax   
        movd mm1, ecx
        movd mm2, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]
	
        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7     
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp       
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm0, [esp + .qqO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm0 ; vcoul=qq*VV
        mulps  xmm0, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and xmm0 fijC.
        ; increment vcoul - then we can get rid of mm5.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	
        ; dispersion
        movlps xmm5, [esi + eax*4 + 16]	;  half table 
        movaps xmm4, xmm5
        shufps xmm4, xmm4, 11111100b
        shufps xmm5, xmm5, 11111101b 
        
        movlps xmm7, [esi + eax*4 + 24] ; other half of dispersion table
        movaps xmm6, xmm7
        shufps xmm6, xmm6, 11111100b
        shufps xmm7, xmm7, 11111101b 
        ; dispersion table ready, in xmm4-xmm7  
        mulss  xmm6, xmm1       ; xmm6=Geps
        mulss  xmm7, xmm2       ; xmm7=Heps2
        addss  xmm5, xmm6
        addss  xmm5, xmm7       ; xmm5=Fp       
        mulss  xmm7, [esp + .two]       ; two*Heps2
        addss  xmm7, xmm6
        addss  xmm7, xmm5 ; xmm7=FF
        mulss  xmm5, xmm1 ; xmm5=eps*Fp
        addss  xmm5, xmm4 ; xmm5=VV

        movaps xmm4, [esp + .c6]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb6
        addps  xmm0, xmm7 ; add to fscal

        ; Update vnbtot directly.
        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .vnbtot], xmm5

        ; repulsion
        movlps xmm5, [esi + eax*4 + 32] ; got half repulsion table
        movaps xmm4, xmm5
        shufps xmm4, xmm4, 10001000b
        shufps xmm5, xmm5, 11011101b 

        movlps xmm7, [esi + eax*4 + 40] ; other half of repulsion table
        movaps xmm6, xmm7
        shufps xmm6, xmm6, 10001000b
        shufps xmm7, xmm7, 11011101b 
        ; repulsion table ready, in xmm4-xmm7  	
        mulss  xmm6, xmm1       ; xmm6=Geps
        mulss  xmm7, xmm2       ; xmm7=Heps2
        addss  xmm5, xmm6
        addss  xmm5, xmm7       ; xmm5=Fp       
        mulss  xmm7, [esp + .two]       ; two*Heps2
        addss  xmm7, xmm6
        addss  xmm7, xmm5 ; xmm7=FF
        mulss  xmm5, xmm1 ; xmm5=eps*Fp
        addss  xmm5, xmm4 ; xmm5=VV

        movaps xmm4, [esp + .c12]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb12
        addps  xmm7, xmm0 ; add to fscal
        addps  xmm5, [esp + .vnbtot] ;  total nonbonded potential in xmm5.

	xorps  xmm4, xmm4
        movd eax, mm0   
        movd ecx, mm1
        movd edx, mm2	
		
	mulps  xmm7, [esp + .rinvO] ;  total fscal now in xmm7
        movaps [esp + .vnbtot], xmm5
	mulps  xmm7, [esp + .tabscale]
	subps xmm4, xmm7

	movaps xmm0, [esp + .dxO]
	movaps xmm1, [esp + .dyO]
	movaps xmm2, [esp + .dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ; xmm0-xmm2 now contains tx-tz (partial force)
	movss  xmm3, [esp + .fixO]	
	movss  xmm4, [esp + .fiyO]	
	movss  xmm5, [esp + .fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + .fixO], xmm3	
	movss  [esp + .fiyO], xmm4	
	movss  [esp + .fizO], xmm5	; updated the O force. now do the H's
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 11100110b	; shift right
	shufps xmm4, xmm4, 11100110b
	shufps xmm5, xmm5, 11100110b
	addss  xmm3, [esp + .fixH1]
	addss  xmm4, [esp + .fiyH1]
	addss  xmm5, [esp + .fizH1]
	movss  [esp + .fixH1], xmm3	
	movss  [esp + .fiyH1], xmm4	
	movss  [esp + .fizH1], xmm5	; updated the H1 force. 

	mov edi, [ebp + %$faction]
	shufps xmm3, xmm3, 11100111b	; shift right
	shufps xmm4, xmm4, 11100111b
	shufps xmm5, xmm5, 11100111b
	addss  xmm3, [esp + .fixH2]
	addss  xmm4, [esp + .fiyH2]
	addss  xmm5, [esp + .fizH2]
	movss  [esp + .fixH2], xmm3	
	movss  [esp + .fiyH2], xmm4	
	movss  [esp + .fizH2], xmm5	; updated the H2 force. 

	; the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1.
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	; x,y sum in xmm0

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1b 
	addss    xmm2, xmm1    ; z sum in xmm2
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .odd_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO]
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	mov   edx, [ebp + %$gid]  
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4	

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		
        
	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 792
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc
	

	
proc inl3330_sse
%$nri		arg
%$iinr		arg
%$jindex	arg
%$jjnr		arg
%$shift		arg
%$shiftvec	arg
%$fshift	arg
%$gid		arg
%$pos		arg		
%$faction	arg
%$charge	arg
%$facel		arg
%$Vc		arg			
%$type		arg
%$ntype  	arg
%$nbfp		arg	
%$Vnb		arg
%$tabscale	arg	
%$VFtab		arg
	;; stack offsets for local variables
	;; bottom of stack is cache-aligned for sse use
.ixO	     equ     0
.iyO	     equ    16
.izO         equ    32
.ixH1	     equ    48
.iyH1	     equ    64
.izH1        equ    80
.ixH2	     equ    96
.iyH2	     equ   112
.izH2        equ   128
.jxO	     equ   144
.jyO	     equ   160
.jzO         equ   176
.jxH1	     equ   192
.jyH1	     equ   208
.jzH1        equ   224
.jxH2	     equ   240
.jyH2	     equ   256
.jzH2        equ   272
.dxOO        equ   288
.dyOO        equ   304
.dzOO        equ   320	
.dxOH1       equ   336
.dyOH1       equ   352
.dzOH1       equ   368	
.dxOH2       equ   384
.dyOH2       equ   400
.dzOH2       equ   416	
.dxH1O       equ   432
.dyH1O       equ   448
.dzH1O       equ   464	
.dxH1H1      equ   480
.dyH1H1      equ   496
.dzH1H1      equ   512	
.dxH1H2      equ   528
.dyH1H2      equ   544
.dzH1H2      equ   560	
.dxH2O       equ   576
.dyH2O       equ   592
.dzH2O       equ   608	
.dxH2H1      equ   624
.dyH2H1      equ   640
.dzH2H1      equ   656	
.dxH2H2      equ   672
.dyH2H2      equ   688
.dzH2H2      equ   704
.qqOO        equ   720
.qqOH        equ   736
.qqHH        equ   752
.two         equ   768
.tabscale    equ   784
.c6          equ   800
.c12	     equ   816		 
.vctot       equ   832
.vnbtot      equ   848
.fixO        equ   864
.fiyO        equ   880
.fizO        equ   896
.fixH1       equ   912
.fiyH1       equ   928
.fizH1       equ   944
.fixH2       equ   960
.fiyH2       equ   976
.fizH2       equ   992
.fjxO	     equ  1008
.fjyO        equ  1024
.fjzO        equ  1040
.fjxH1	     equ  1056
.fjyH1       equ  1072
.fjzH1       equ  1088
.fjxH2	     equ  1104
.fjyH2       equ  1120
.fjzH2       equ  1136
.half        equ  1152
.three       equ  1168
.rsqOO       equ  1184
.rsqOH1      equ  1200
.rsqOH2      equ  1216
.rsqH1O      equ  1232
.rsqH1H1     equ  1248
.rsqH1H2     equ  1264
.rsqH2O      equ  1280
.rsqH2H1     equ  1296
.rsqH2H2     equ  1312
.rinvOO      equ  1328
.rinvOH1     equ  1344
.rinvOH2     equ  1360
.rinvH1O     equ  1376
.rinvH1H1    equ  1392
.rinvH1H2    equ  1408
.rinvH2O     equ  1424
.rinvH2H1    equ  1440
.rinvH2H2    equ  1456
.fstmp	     equ  1472	
.is3         equ  1488
.ii3         equ  1492
.innerjjnr   equ  1496
.innerk      equ  1500
.salign	     equ  1504							
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 1508		; local stack space
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + .salign], eax

	emms

	movups xmm0, [sse_half]
	movups xmm1, [sse_two]
	movups xmm2, [sse_three]
	movss xmm3, [ebp +%$tabscale]
	movaps [esp + .half],  xmm0
	movaps [esp + .two],  xmm1
	movaps [esp + .three], xmm2
	shufps xmm3, xmm3, 0b
	movaps [esp + .tabscale],  xmm3

	;; assume we have at least one i particle - start directly
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, xmm3	
	movss xmm5, [edx + ebx*4 + 4]	
	movss xmm6, [ebp + %$facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .qqOO], xmm3
	movaps [esp + .qqOH], xmm4
	movaps [esp + .qqHH], xmm5
		
	xorps xmm0, xmm0
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movlps xmm0, [eax + edx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 01010101b
	movaps [esp + .c6], xmm0
	movaps [esp + .c12], xmm1

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixO], xmm3
	movaps [esp + .iyO], xmm4
	movaps [esp + .izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0b
	shufps xmm1, xmm1, 0b
	shufps xmm2, xmm2, 0b
	shufps xmm3, xmm3, 0b
	shufps xmm4, xmm4, 0b
	shufps xmm5, xmm5, 0b
	movaps [esp + .ixH1], xmm0
	movaps [esp + .iyH1], xmm1
	movaps [esp + .izH1], xmm2
	movaps [esp + .ixH2], xmm3
	movaps [esp + .iyH2], xmm4
	movaps [esp + .izH2], xmm5

	; clear vctot and i forces
	xorps xmm4, xmm4
	movaps [esp + .vctot], xmm4
	movaps [esp + .vnbtot], xmm4
	movaps [esp + .fixO], xmm4
	movaps [esp + .fiyO], xmm4
	movaps [esp + .fizO], xmm4
	movaps [esp + .fixH1], xmm4
	movaps [esp + .fiyH1], xmm4
	movaps [esp + .fizH1], xmm4
	movaps [esp + .fixH2], xmm4
	movaps [esp + .fiyH2], xmm4
	movaps [esp + .fizH2], xmm4
	
	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
	sub   edx, dword 4
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .single_check
.unroll_loop:	
	;; quad-unroll innerloop here.
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 

	mov   eax, [edx]	
	mov   ebx, [edx + 4] 
	mov   ecx, [edx + 8]
	mov   edx, [edx + 12]             ; eax-edx=jnr1-4 
	
	add   [esp + .innerjjnr], dword 16 ; advance pointer (unrolled 4) 

	mov esi, [ebp + %$pos]        ; base of pos[]

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]         ;  replace jnr with j3
	lea   edx, [edx + edx*2]	
	
	; move j coordinates to local temp. variables
	movlps xmm2, [esi + eax*4]
	movlps xmm3, [esi + eax*4 + 12]
	movlps xmm4, [esi + eax*4 + 24]

	movlps xmm5, [esi + ebx*4]
	movlps xmm6, [esi + ebx*4 + 12]
	movlps xmm7, [esi + ebx*4 + 24]

	movhps xmm2, [esi + ecx*4]
	movhps xmm3, [esi + ecx*4 + 12]
	movhps xmm4, [esi + ecx*4 + 24]

	movhps xmm5, [esi + edx*4]
	movhps xmm6, [esi + edx*4 + 12]
	movhps xmm7, [esi + edx*4 + 24]

	;; current state:	
	;;  xmm2= jxOa  jyOa  jxOc  jyOc
	;;  xmm3= jxH1a jyH1a jxH1c jyH1c
	;;  xmm4= jxH2a jyH2a jxH2c jyH2c
	;;  xmm5= jxOb  jyOb  jxOd  jyOd
	;;  xmm6= jxH1b jyH1b jxH1d jyH1d
	;;  xmm7= jxH2b jyH2b jxH2d jyH2d
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	; xmm0= jxOa  jxOb  jyOa  jyOb
	unpcklps xmm1, xmm6	; xmm1= jxH1a jxH1b jyH1a jyH1b
	unpckhps xmm2, xmm5	; xmm2= jxOc  jxOd  jyOc  jyOd
	unpckhps xmm3, xmm6	; xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	; xmm4= jxH2a jxH2b jyH2a jyH2b		
	unpckhps xmm5, xmm7	; xmm5= jxH2c jxH2d jyH2c jyH2d
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	; xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [esp + .jxO], xmm0
	movhlps  xmm2, xmm6	; xmm2= jyOa  jyOb  jyOc  jyOd
	movaps [esp + .jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [esp + .jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [esp + .jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [esp + .jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [esp + .jyH2], xmm5

	movss  xmm0, [esi + eax*4 + 8]
	movss  xmm1, [esi + eax*4 + 20]
	movss  xmm2, [esi + eax*4 + 32]

	movss  xmm3, [esi + ecx*4 + 8]
	movss  xmm4, [esi + ecx*4 + 20]
	movss  xmm5, [esi + ecx*4 + 32]

	movhps xmm0, [esi + ebx*4 + 4]
	movhps xmm1, [esi + ebx*4 + 16]
	movhps xmm2, [esi + ebx*4 + 28]
	
	movhps xmm3, [esi + edx*4 + 4]
	movhps xmm4, [esi + edx*4 + 16]
	movhps xmm5, [esi + edx*4 + 28]
	
	shufps xmm0, xmm3, 11001100b
	shufps xmm1, xmm4, 11001100b
	shufps xmm2, xmm5, 11001100b
	movaps [esp + .jzO],  xmm0
	movaps [esp + .jzH1],  xmm1
	movaps [esp + .jzH2],  xmm2

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixO]
	movaps xmm4, [esp + .iyO]
	movaps xmm5, [esp + .izO]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxOH1], xmm3
	movaps [esp + .dyOH1], xmm4
	movaps [esp + .dzOH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOO], xmm0
	movaps [esp + .rsqOH1], xmm3

	movaps xmm0, [esp + .ixO]
	movaps xmm1, [esp + .iyO]
	movaps xmm2, [esp + .izO]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	subps  xmm3, [esp + .jxO]
	subps  xmm4, [esp + .jyO]
	subps  xmm5, [esp + .jzO]
	movaps [esp + .dxOH2], xmm0
	movaps [esp + .dyOH2], xmm1
	movaps [esp + .dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1O], xmm3
	movaps [esp + .dyH1O], xmm4
	movaps [esp + .dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqOH2], xmm0
	movaps [esp + .rsqH1O], xmm3

	movaps xmm0, [esp + .ixH1]
	movaps xmm1, [esp + .iyH1]
	movaps xmm2, [esp + .izH1]
	movaps xmm3, [esp + .ixH1]
	movaps xmm4, [esp + .iyH1]
	movaps xmm5, [esp + .izH1]
	subps  xmm0, [esp + .jxH1]
	subps  xmm1, [esp + .jyH1]
	subps  xmm2, [esp + .jzH1]
	subps  xmm3, [esp + .jxH2]
	subps  xmm4, [esp + .jyH2]
	subps  xmm5, [esp + .jzH2]
	movaps [esp + .dxH1H1], xmm0
	movaps [esp + .dyH1H1], xmm1
	movaps [esp + .dzH1H1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH1H2], xmm3
	movaps [esp + .dyH1H2], xmm4
	movaps [esp + .dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [esp + .rsqH1H1], xmm0
	movaps [esp + .rsqH1H2], xmm3

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	movaps xmm3, [esp + .ixH2]
	movaps xmm4, [esp + .iyH2]
	movaps xmm5, [esp + .izH2]
	subps  xmm0, [esp + .jxO]
	subps  xmm1, [esp + .jyO]
	subps  xmm2, [esp + .jzO]
	subps  xmm3, [esp + .jxH1]
	subps  xmm4, [esp + .jyH1]
	subps  xmm5, [esp + .jzH1]
	movaps [esp + .dxH2O], xmm0
	movaps [esp + .dyH2O], xmm1
	movaps [esp + .dzH2O], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [esp + .dxH2H1], xmm3
	movaps [esp + .dyH2H1], xmm4
	movaps [esp + .dzH2H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [esp + .rsqH2O], xmm0
	movaps [esp + .rsqH2H1], xmm4

	movaps xmm0, [esp + .ixH2]
	movaps xmm1, [esp + .iyH2]
	movaps xmm2, [esp + .izH2]
	subps  xmm0, [esp + .jxH2]
	subps  xmm1, [esp + .jyH2]
	subps  xmm2, [esp + .jzH2]
	movaps [esp + .dxH2H2], xmm0
	movaps [esp + .dyH2H2], xmm1
	movaps [esp + .dzH2H2], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [esp + .rsqH2H2], xmm0
		
	; start doing invsqrt. use rsq values in xmm0, xmm4
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinvH2H2
	mulps   xmm7, [esp + .half] ; rinvH2H1
	movaps  [esp + .rinvH2H2], xmm3
	movaps  [esp + .rinvH2H1], xmm7
		
	rsqrtps xmm1, [esp + .rsqOO]
	rsqrtps xmm5, [esp + .rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOO]
	mulps   xmm5, [esp + .rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOO], xmm3
	movaps  [esp + .rinvOH1], xmm7
	
	rsqrtps xmm1, [esp + .rsqOH2]
	rsqrtps xmm5, [esp + .rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqOH2]
	mulps   xmm5, [esp + .rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvOH2], xmm3
	movaps  [esp + .rinvH1O], xmm7
	
	rsqrtps xmm1, [esp + .rsqH1H1]
	rsqrtps xmm5, [esp + .rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, [esp + .rsqH1H1]
	mulps   xmm5, [esp + .rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] 
	mulps   xmm7, [esp + .half]
	movaps  [esp + .rinvH1H1], xmm3
	movaps  [esp + .rinvH1H2], xmm7
	
	rsqrtps xmm1, [esp + .rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, [esp + .rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [esp + .half] 
	movaps  [esp + .rinvH2O], xmm3

	;; start with OO interaction.
	movaps xmm0, [esp + .rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOO] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
		
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
	
        movd mm0, eax
        movd mm1, ebx
        movd mm2, ecx
        movd mm3, edx

        mov  esi, [ebp + %$VFtab]
        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOO]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.
        ; increment vcoul - then we can get rid of mm5.
        ;; update vctot
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5 

        ; put scalar force on stack temporarily...
        movaps [esp + .fstmp], xmm3

        ; dispersion
        movlps xmm5, [esi + eax*4 + 16]
        movlps xmm7, [esi + ecx*4 + 16]
        movhps xmm5, [esi + ebx*4 + 16]
        movhps xmm7, [esi + edx*4 + 16] ; got half dispersion table
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 24]
        movlps xmm3, [esi + ecx*4 + 24]
        movhps xmm7, [esi + ebx*4 + 24]
        movhps xmm3, [esi + edx*4 + 24] ; other half of dispersion table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; dispersion table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV

        movaps xmm4, [esp + .c6]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb6
        addps  xmm7, [esp + .fstmp] ; add to fscal

        ; put scalar force on stack. Update vnbtot directly.
        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .fstmp], xmm7
        movaps [esp + .vnbtot], xmm5

        ; repulsion
        movlps xmm5, [esi + eax*4 + 32]
        movlps xmm7, [esi + ecx*4 + 32]
        movhps xmm5, [esi + ebx*4 + 32]
        movhps xmm7, [esi + edx*4 + 32] ; got half repulsion table
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 40]
        movlps xmm3, [esi + ecx*4 + 40]
        movhps xmm7, [esi + ebx*4 + 40]
        movhps xmm3, [esi + edx*4 + 40] ; other half of repulsion table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
 
        movaps xmm4, [esp + .c12]
        mulps  xmm7, xmm4 ; fijR
        mulps  xmm5, xmm4 ; vnb12
        addps  xmm7, [esp + .fstmp] 

        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .vnbtot], xmm5
        xorps  xmm1, xmm1

        mulps xmm7, [esp + .tabscale]
        mulps xmm7, xmm0
        subps  xmm1, xmm7

	movaps xmm0, xmm1
	movaps xmm2, xmm1		

	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOO]
	mulps xmm1, [esp + .dyOO]
	mulps xmm2, [esp + .dzOO]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H1 interaction
	movaps xmm0, [esp + .rinvOH1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH1]
	mulps xmm1, [esp + .dyOH1]
	mulps xmm2, [esp + .dzOH1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; O-H2 interaction
	movaps xmm0, [esp + .rinvOH2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqOH2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	xorps xmm3, xmm3
	movaps xmm4, xmm3
	movaps xmm5, xmm3
	mulps xmm0, [esp + .dxOH2]
	mulps xmm1, [esp + .dyOH2]
	mulps xmm2, [esp + .dzOH2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixO]
	addps xmm1, [esp + .fiyO]
	addps xmm2, [esp + .fizO]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixO], xmm0
	movaps [esp + .fiyO], xmm1
	movaps [esp + .fizO], xmm2

	; H1-O interaction
	movaps xmm0, [esp + .rinvH1O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH1O]
	mulps xmm1, [esp + .dyH1O]
	mulps xmm2, [esp + .dzH1O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H1 interaction
	movaps xmm0, [esp + .rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH1H1]
	mulps xmm1, [esp + .dyH1H1]
	mulps xmm2, [esp + .dzH1H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H1-H2 interaction
	movaps xmm0, [esp + .rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH1H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH1H2]
	mulps xmm1, [esp + .dyH1H2]
	mulps xmm2, [esp + .dzH1H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH1]
	addps xmm1, [esp + .fiyH1]
	addps xmm2, [esp + .fizH1]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH1], xmm0
	movaps [esp + .fiyH1], xmm1
	movaps [esp + .fizH1], xmm2

	; H2-O interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2O] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqOH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1

	movaps xmm3, [esp + .fjxO]
	movaps xmm4, [esp + .fjyO]
	movaps xmm5, [esp + .fjzO]
	mulps xmm0, [esp + .dxH2O]
	mulps xmm1, [esp + .dyH2O]
	mulps xmm2, [esp + .dzH2O]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxO], xmm3
	movaps [esp + .fjyO], xmm4
	movaps [esp + .fjzO], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H1 interaction
	movaps xmm0, [esp + .rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H1] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH1]
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	mulps xmm0, [esp + .dxH2H1]
	mulps xmm1, [esp + .dyH2H1]
	mulps xmm2, [esp + .dzH2H1]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH1], xmm3
	movaps [esp + .fjyH1], xmm4
	movaps [esp + .fjzH1], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	; H2-H2 interaction
	movaps xmm0, [esp + .rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [esp + .rsqH2H2] ; xmm1=r
	mulps  xmm1, [esp + .tabscale]	
	movhlps xmm2, xmm1
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ;xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2

        movd eax, mm6
        psrlq mm6, 32
        movd ecx, mm7
        psrlq mm7, 32
        movd ebx, mm6
        movd edx, mm7

        lea   eax, [eax + eax*2]
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]

        movlps xmm5, [esi + eax*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm5, [esi + ebx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 

        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + eax*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm7, [esi + ebx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 

        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2
        movaps xmm3, [esp + .qqHH]
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point mm5 contains vcoul and mm3 fijC.

        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
	xorps  xmm1, xmm1
	mulps  xmm3,  [esp + .tabscale]
	mulps  xmm3, xmm0
	subps  xmm1, xmm3

	movaps xmm0, xmm1
	movaps xmm2, xmm1
	
	movaps xmm3, [esp + .fjxH2]
	movaps xmm4, [esp + .fjyH2]
	movaps xmm5, [esp + .fjzH2]
	mulps xmm0, [esp + .dxH2H2]
	mulps xmm1, [esp + .dyH2H2]
	mulps xmm2, [esp + .dzH2H2]
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	addps xmm0, [esp + .fixH2]
	addps xmm1, [esp + .fiyH2]
	addps xmm2, [esp + .fizH2]
	movaps [esp + .fjxH2], xmm3
	movaps [esp + .fjyH2], xmm4
	movaps [esp + .fjzH2], xmm5
	movaps [esp + .fixH2], xmm0
	movaps [esp + .fiyH2], xmm1
	movaps [esp + .fizH2], xmm2

	mov edi, [ebp +%$faction]

	movd eax, mm0
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3
	
	; Did all interactions - now update j forces.
	; 4 j waters with three atoms each - first do a & b j particles
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpcklps xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjxOb  fjyOb
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOb  fjyOb
	unpcklps xmm1, xmm2	   ; xmm1= fjzOa  fjxH1a fjzOb  fjxH1b
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpcklps xmm4, xmm5	   ; xmm4= fjyH1a fjzH1a fjyH1b fjzH1b
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1b fjzH1b
	unpcklps xmm5, xmm6	   ; xmm5= fjxH2a fjyH2a fjxH2b fjyH2b
	movlhps  xmm0, xmm1    	   ; xmm0= fjxOa  fjyOa  fjzOa  fjxH1a  
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOb  fjyOb  fjzOb  fjxH1b
	movlhps  xmm4, xmm5   	   ; xmm4= fjyH1a fjzH1a fjxH2a fjyH2a
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1b fjzH1b fjxH2b fjyH2b
	movups   xmm1, [edi + eax*4]
	movups   xmm2, [edi + eax*4 + 16]
	movups   xmm5, [edi + ebx*4]
	movups   xmm6, [edi + ebx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + eax*4 + 32]
	movss    xmm3, [edi + ebx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm7, xmm7, 1b
	
	movups   [edi + eax*4],     xmm1
	movups   [edi + eax*4 + 16],xmm2
	movups   [edi + ebx*4],     xmm5
	movups   [edi + ebx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + eax*4 + 32], xmm0
	movss    [edi + ebx*4 + 32], xmm3	

	;; then do the second pair (c & d)
	movaps xmm0, [esp + .fjxO] ; xmm0= fjxOa  fjxOb  fjxOc  fjxOd
	movaps xmm1, [esp + .fjyO] ; xmm1= fjyOa  fjyOb  fjyOc  fjyOd 
	unpckhps xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjxOd  fjyOd
	movaps xmm1, [esp + .fjzO]
	movaps xmm2, [esp + .fjxH1]
	movhlps  xmm3, xmm0	   ; xmm3= fjxOd  fjyOd
	unpckhps xmm1, xmm2	   ; xmm1= fjzOc  fjxH1c fjzOd  fjxH1d
	movaps xmm4, [esp + .fjyH1]
	movaps xmm5, [esp + .fjzH1]
	unpckhps xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjyH1d fjzH1d	
	movaps xmm5, [esp + .fjxH2]
	movaps xmm6, [esp + .fjyH2]
	movhlps  xmm7, xmm4	   ; xmm7= fjyH1d fjzH1d	 
	unpckhps xmm5, xmm6	   ; xmm5= fjxH2c fjyH2c fjxH2d fjyH2d
	movlhps  xmm0, xmm1	   ; xmm0= fjxOc  fjyOc  fjzOc  fjxH1c 
	shufps   xmm3, xmm1, 11100100b
                                   ; xmm3= fjxOd  fjyOd  fjzOd  fjxH1d
	movlhps  xmm4, xmm5	   ; xmm4= fjyH1c fjzH1c fjxH2c fjyH2c 
	shufps   xmm7, xmm5, 11100100b
                                   ; xmm7= fjyH1d fjzH1d fjxH2d fjyH2d
	movups   xmm1, [edi + ecx*4]
	movups   xmm2, [edi + ecx*4 + 16]
	movups   xmm5, [edi + edx*4]
	movups   xmm6, [edi + edx*4 + 16]
	addps    xmm1, xmm0
	addps    xmm2, xmm4
	addps    xmm5, xmm3
	addps    xmm6, xmm7
	movss    xmm0, [edi + ecx*4 + 32]
	movss    xmm3, [edi + edx*4 + 32]
	
	movaps   xmm4, [esp + .fjzH2]
	movaps   xmm7, xmm4
	shufps   xmm4, xmm4, 10b
	shufps   xmm7, xmm7, 11b
	movups   [edi + ecx*4],     xmm1
	movups   [edi + ecx*4 + 16],xmm2
	movups   [edi + edx*4],     xmm5
	movups   [edi + edx*4 + 16],xmm6	
	addss    xmm0, xmm4
	addss    xmm3, xmm7
	movss    [edi + ecx*4 + 32], xmm0
	movss    [edi + edx*4 + 32], xmm3	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 4
	jl    .single_check
	jmp   .unroll_loop
.single_check:
	add   [esp + .innerk], dword 4
	jnz   .single_loop
	jmp   .updateouterdata
.single_loop:
	mov   edx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [edx]	
	add   [esp + .innerjjnr], dword 4	

	mov esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]  

	; fetch j coordinates
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + eax*4 + 4]
	movss xmm5, [esi + eax*4 + 8]

	movlps xmm6, [esi + eax*4 + 12]
	movhps xmm6, [esi + eax*4 + 24]	; xmm6=jxH1 jyH1 jxH2 jyH2
	;;  fetch both z coords in one go, to positions 0 and 3 in xmm7
	movups xmm7, [esi + eax*4 + 20] ; xmm7=jzH1 jxH2 jyH2 jzH2
	shufps xmm6, xmm6, 11011000b    ;  xmm6=jxH1 jxH2 jyH1 jyH2
	movlhps xmm3, xmm6      	; xmm3= jxO   0  jxH1 jxH2 
	movaps  xmm0, [esp + .ixO]     
	movaps  xmm1, [esp + .iyO]
	movaps  xmm2, [esp + .izO]	
	shufps  xmm4, xmm6, 11100100b ;  xmm4= jyO   0   jyH1 jyH2
	shufps xmm5, xmm7, 11000100b  ;  xmm5= jzO   0   jzH1 jzH2  
	;;  store all j coordinates in jO
	movaps [esp + .jxO], xmm3
	movaps [esp + .jyO], xmm4
	movaps [esp + .jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	movaps [esp + .dxOO], xmm0
	movaps [esp + .dyOO], xmm1
	movaps [esp + .dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;  have rsq in xmm0.
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [esp + .three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [esp + .half] ; rinv iO - j water

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	; xmm1=r
	movaps  xmm0, xmm3	; xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

	mov esi, [ebp + %$VFtab]
	
        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]
	
        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOO]
	movhps  xmm3, [esp + .qqOH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
	
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5
        ; put scalar force on stack temporarily...
        movaps [esp + .fstmp], xmm3

        ; dispersion
	movss  xmm4, [esi + ebx*4 + 16]	
	movss  xmm5, [esi + ebx*4 + 20]	
	movss  xmm6, [esi + ebx*4 + 24]	
	movss  xmm7, [esi + ebx*4 + 28]
        ; dispersion table ready, in xmm4-xmm7 
        mulss  xmm6, xmm1       ; xmm6=Geps
        mulss  xmm7, xmm2       ; xmm7=Heps2
        addss  xmm5, xmm6
        addss  xmm5, xmm7       ; xmm5=Fp
        mulss  xmm7, [esp + .two]       ; two*Heps2
        addss  xmm7, xmm6
        addss  xmm7, xmm5 ; xmm7=FF
        mulss  xmm5, xmm1 ; xmm5=eps*Fp
        addss  xmm5, xmm4 ; xmm5=VV
	xorps  xmm4, xmm4
        movss  xmm4, [esp + .c6]
        mulps  xmm7, xmm4        ; fijD
        mulps  xmm5, xmm4        ; vnb6
        addps  xmm7, [esp + .fstmp] ; add to fscal

        ; put scalar force on stack. Update vnbtot directly.
        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .fstmp], xmm7
        movaps [esp + .vnbtot], xmm5

        ; repulsion
	movss  xmm4, [esi + ebx*4 + 32]	
	movss  xmm5, [esi + ebx*4 + 36]	
	movss  xmm6, [esi + ebx*4 + 40]	
	movss  xmm7, [esi + ebx*4 + 44]
        ; table ready, in xmm4-xmm7 
        mulss  xmm6, xmm1       ; xmm6=Geps
        mulss  xmm7, xmm2       ; xmm7=Heps2
        addss  xmm5, xmm6
        addss  xmm5, xmm7       ; xmm5=Fp
        mulss  xmm7, [esp + .two]       ; two*Heps2
        addss  xmm7, xmm6
        addss  xmm7, xmm5 ; xmm7=FF
        mulss  xmm5, xmm1 ; xmm5=eps*Fp
        addss  xmm5, xmm4 ; xmm5=VV

	xorps  xmm4, xmm4
        movss  xmm4, [esp + .c12]
        mulps  xmm7, xmm4 ; fijR
        mulps  xmm5, xmm4 ; vnb12
        addps  xmm7, [esp + .fstmp] 

        addps  xmm5, [esp + .vnbtot]
        movaps [esp + .vnbtot], xmm5
        xorps  xmm1, xmm1

        mulps xmm7, [esp + .tabscale]
        mulps xmm7, xmm0
        subps  xmm1, xmm7

	movaps xmm0, xmm1
	movaps xmm2, xmm1		

	mulps   xmm0, [esp + .dxOO]
	mulps   xmm1, [esp + .dyOO]
	mulps   xmm2, [esp + .dzOO]
	;; initial update for j forces
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixO]
	addps   xmm1, [esp + .fiyO]
	addps   xmm2, [esp + .fizO]
	movaps  [esp + .fixO], xmm0
	movaps  [esp + .fiyO], xmm1
	movaps  [esp + .fizO], xmm2

	
	;;  done with i O. Now do i H1 & H2 simultaneously. first get i particle coords:
	movaps  xmm0, [esp + .ixH1]
	movaps  xmm1, [esp + .iyH1]
	movaps  xmm2, [esp + .izH1]	
	movaps  xmm3, [esp + .ixH2] 
	movaps  xmm4, [esp + .iyH2] 
	movaps  xmm5, [esp + .izH2] 
	subps   xmm0, [esp + .jxO]
	subps   xmm1, [esp + .jyO]
	subps   xmm2, [esp + .jzO]
	subps   xmm3, [esp + .jxO]
	subps   xmm4, [esp + .jyO]
	subps   xmm5, [esp + .jzO]
	movaps [esp + .dxH1O], xmm0
	movaps [esp + .dyH1O], xmm1
	movaps [esp + .dzH1O], xmm2
	movaps [esp + .dxH2O], xmm3
	movaps [esp + .dyH2O], xmm4
	movaps [esp + .dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;  have rsqH1 in xmm0.
	addps xmm4, xmm5	;  have rsqH2 in xmm4.

	;;  start with H1, save H2 data
	movaps [esp + .rsqH2O], xmm4
	
	;;  do invsqrt
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [esp + .three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [esp + .half] ; rinv H1 - j water
	mulps   xmm7, [esp + .half] ; rinv H2 - j water

	;;  start with H1, save H2 data
	movaps [esp + .rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;  xmm1=r
	movaps xmm0, xmm3	;  xmm0=rinv
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]
	
        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [esp + .dxH1O]
	mulps   xmm1, [esp + .dyH1O]
	mulps   xmm2, [esp + .dzH1O]
	;;  update forces H1 - j water
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH1]
	addps   xmm1, [esp + .fiyH1]
	addps   xmm2, [esp + .fizH1]
	movaps  [esp + .fixH1], xmm0
	movaps  [esp + .fiyH1], xmm1
	movaps  [esp + .fizH1], xmm2
	;; do table for H2 - j water interaction
	movaps xmm0, [esp + .rinvH2O]
	movaps xmm1, [esp + .rsqH2O]
	mulps  xmm1, xmm0	; xmm0=rinv, xmm1=r
	mulps  xmm1, [esp + .tabscale]
	
	movhlps xmm2, xmm1	
        cvttps2pi mm6, xmm1
        cvttps2pi mm7, xmm2     ; mm6/mm7 contain lu indices
        cvtpi2ps xmm3, mm6
        cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;  xmm1=eps
        movaps xmm2, xmm1
        mulps  xmm2, xmm2       ; xmm2=eps2
        pslld mm6, 2
        pslld mm7, 2
        movd ebx, mm6
        movd ecx, mm7
        psrlq mm7, 32
        movd edx, mm7		; table indices in ebx,ecx,edx

        lea   ebx, [ebx + ebx*2]
        lea   ecx, [ecx + ecx*2]
        lea   edx, [edx + edx*2]
	
        movlps xmm5, [esi + ebx*4]
        movlps xmm7, [esi + ecx*4]
        movhps xmm7, [esi + edx*4] ; got half coulomb table 
        movaps xmm4, xmm5
        shufps xmm4, xmm7, 10001000b
        shufps xmm5, xmm7, 11011101b 

        movlps xmm7, [esi + ebx*4 + 8]
        movlps xmm3, [esi + ecx*4 + 8]
        movhps xmm3, [esi + edx*4 + 8] ; other half of coulomb table
        movaps xmm6, xmm7
        shufps xmm6, xmm3, 10001000b
        shufps xmm7, xmm3, 11011101b 
        ; coulomb table ready, in xmm4-xmm7 
        mulps  xmm6, xmm1       ; xmm6=Geps
        mulps  xmm7, xmm2       ; xmm7=Heps2
        addps  xmm5, xmm6
        addps  xmm5, xmm7       ; xmm5=Fp
        mulps  xmm7, [esp + .two]       ; two*Heps2

	xorps  xmm3, xmm3
	;; fetch charges to xmm3 (temporary)
	movss   xmm3, [esp + .qqOH]
	movhps  xmm3, [esp + .qqHH]
		
        addps  xmm7, xmm6
        addps  xmm7, xmm5 ; xmm7=FF
        mulps  xmm5, xmm1 ; xmm5=eps*Fp
        addps  xmm5, xmm4 ; xmm5=VV
        mulps  xmm5, xmm3 ; vcoul=qq*VV
        mulps  xmm3, xmm7 ; fijC=FF*qq
        ; at this point xmm5 contains vcoul and xmm3 fijC.
        addps  xmm5, [esp + .vctot]
        movaps [esp + .vctot], xmm5	

        xorps  xmm1, xmm1

        mulps xmm3, [esp + .tabscale]
        mulps xmm3, xmm0
        subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [esp + .dxH2O]
	mulps   xmm1, [esp + .dyH2O]
	mulps   xmm2, [esp + .dzH2O]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	subps   xmm3, xmm0
	subps   xmm4, xmm1
	subps   xmm5, xmm2
	mov     esi, [ebp + %$faction]
	movaps  [esp + .fjxO], xmm3
	movaps  [esp + .fjyO], xmm4
	movaps  [esp + .fjzO], xmm5
	addps   xmm0, [esp + .fixH2]
	addps   xmm1, [esp + .fiyH2]
	addps   xmm2, [esp + .fizH2]
	movaps  [esp + .fixH2], xmm0
	movaps  [esp + .fiyH2], xmm1
	movaps  [esp + .fizH2], xmm2

	;; update j water forces from local variables
	movlps  xmm0, [esi + eax*4]
	movlps  xmm1, [esi + eax*4 + 12]
	movhps  xmm1, [esi + eax*4 + 24]
	movaps  xmm3, [esp + .fjxO]
	movaps  xmm4, [esp + .fjyO]
	movaps  xmm5, [esp + .fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 10b
	shufps  xmm7, xmm7, 11b
	addss   xmm5, [esi + eax*4 + 8]
	addss   xmm6, [esi + eax*4 + 20]
	addss   xmm7, [esi + eax*4 + 32]
	movss   [esi + eax*4 + 8], xmm5
	movss   [esi + eax*4 + 20], xmm6
	movss   [esi + eax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [esi + eax*4], xmm0 
	movlps  [esi + eax*4 + 12], xmm1 
	movhps  [esi + eax*4 + 24], xmm1 
	
	dec   dword [esp + .innerk]
	jz    .updateouterdata
	jmp   .single_loop
.updateouterdata:
	mov   ecx, [esp + .ii3]
	mov   edi, [ebp + %$faction]
	mov   esi, [ebp + %$fshift]
	mov   edx, [esp + .is3]

	; accumulate Oi forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixO]
	movaps xmm1, [esp + .fiyO] 
	movaps xmm2, [esp + .fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	; accumulate force in xmm6/xmm7 for fshift
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 1000b	

	; accumulate H1i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH1]
	movaps xmm1, [esp + .fiyH1]
	movaps xmm2, [esp + .fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; accumulate H2i forces in xmm0, xmm1, xmm2
	movaps xmm0, [esp + .fixH2]
	movaps xmm1, [esp + .fiyH2]
	movaps xmm2, [esp + .fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ; summan ligger i 1/2 i xmm0-xmm2

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1b
	shufps xmm4, xmm4, 1b
	shufps xmm5, xmm5, 1b
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	; xmm0-xmm2 has single force in pos0.

	; increment i force
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;accumulate force in xmm6/xmm7 for fshift
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 1000b	
	addps   xmm6, xmm0

	; increment fshift force
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	; get group index for i particle
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	; accumulate total potential energy and update it.
	movaps xmm7, [esp + .vctot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vc]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	; accumulate total lj energy and update it.
	movaps xmm7, [esp + .vnbtot]
	; accumulate
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	; pos 0-1 in xmm7 have the sum now
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1b
	addss  xmm7, xmm6		

	; add earlier value from mem.
	mov   eax, [ebp + %$Vnb]
	addss xmm7, [eax + edx*4] 
	; move back to mem.
	movss [eax + edx*4], xmm7 
	
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	emms
	mov eax, [esp + .salign]
	add esp, eax
	add esp, 1508
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


