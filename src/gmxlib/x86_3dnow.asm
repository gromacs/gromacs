

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

	;;	This file contains a subset of the gromacs innerloops
	;;      manually written in assembly to optimize performance
	;;      on AMD extended 3DNow-enabled processors like Athlon 
	;; 	and later generations.
	;;	Erik Lindahl, 2000, erik@theophys.kth.se
	
segment .data
mm_two
	dd 2.0
	dd 2.0
mm_six
	dd 6.0
	dd 6.0
mm_twelve
	dd 12.0
	dd 12.0


segment .text


	global check3dnow       ;  tries to issue a simple 3DNOW instruction 	
check3dnow:	
	femms
	pfmul mm0,mm0
	femms
	ret
	
	
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
	pfrcpit1 mm0,mm1
	punpckldq mm5,mm6	
	pfrcpit2 mm0,mm1
        movq [ebx],mm0
	pfrcpit1 mm4,mm5
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
        pfrcpit1 mm0,mm1
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
	


proc inl0100_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.vnbtot      equ    20           ; repeated (64bit) to fill 3dnow reg
.c6          equ    28           ; repeated (64bit) to fill 3dnow reg
.c12         equ    36           ; repeated (64bit) to fill 3dnow reg
.six         equ    44           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    52           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    60
.innerjjnr   equ    64
.innerk      equ    68		
.fix         equ    72
.fiy         equ    76
.fiz	     equ    80
.dx1	     equ    84
.dy1	     equ    88
.dz1	     equ    92
.dx2	     equ    96
.dy2	     equ   100
.dz2	     equ   104						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + .six], mm0
	movq  [esp + .twelve], mm1
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

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
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best	

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrcp mm0, mm4	                 ; lookup reciprocal seed 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
                  	        	; amd 3dnow N-R iteration to get full precision.
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	;; mm4 now contains invsq,
	;; do potential and fscal
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	;  mm4=invsq 
	;;  calculate potentials and scalar force
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



		
		
proc inl0110_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.vnbtot      equ    32           ; repeated (64bit) to fill 3dnow reg
.c6          equ    40           ; repeated (64bit) to fill 3dnow reg
.c12         equ    48           ; repeated (64bit) to fill 3dnow reg
.six         equ    56           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    64           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    72	
.innerjjnr0  equ    76
.innerk0     equ    80		
.innerjjnr   equ    84
.innerk      equ    88	
.fix         equ    92
.fiy         equ    96
.fiz	     equ    100
.dx1	     equ    104
.dy1	     equ    108
.dz1	     equ    112
.dx2	     equ    116
.dy2	     equ    120
.dz2	     equ    124								
.nsvdwc      equ    128
.nscoul      equ    132
.nsvdw       equ    136
.solnr	     equ    140		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 144		;  local stack space
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + .six],    mm0
	movq  [esp + .twelve], mm1
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

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
	pxor  mm7,mm7
	movq  [esp + .vnbtot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	
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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best	

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrcp mm0, mm4	                 ; lookup reciprocal seed 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
                  	        	; amd 3dnow N-R iteration to get full precision.
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	;; mm4 now contains invsq,
	;; do potential and fscal
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdwc_inner
	jmp  .updateouterdata_vdwc		
.single_vdwc_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	;  mm4=invsq 
	;;  calculate potentials and scalar force
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdwc:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best	

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrcp mm0, mm4	                 ; lookup reciprocal seed 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
                  	        	; amd 3dnow N-R iteration to get full precision.
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	;; mm4 now contains invsq,
	;; do potential and fscal
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdw_inner
	jmp  .updateouterdata_vdw		
.single_vdw_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	;  mm4=invsq 
	;;  calculate potentials and scalar force
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdw:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 144
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl0300_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.vnbtot      equ    20           ; repeated (64bit) to fill 3dnow reg
.c6          equ    28           ; repeated (64bit) to fill 3dnow reg
.c12         equ    36           ; repeated (64bit) to fill 3dnow reg
.two         equ    44           ; repeated (64bit) to fill 3dnow reg
.n1          equ    52           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    60           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    68
.innerjjnr   equ    72
.innerk      equ    76		
.fix         equ    80
.fiy         equ    84
.fiz	     equ    88
.dx1	     equ    92
.dy1	     equ    96
.dz1	     equ    100
.dx2	     equ    104
.dy2	     equ    108
.dz2	     equ    112						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 116		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

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
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	; dispersion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 3
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 116
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


			
	
proc inl0310_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.vnbtot      equ    32           ; repeated (64bit) to fill 3dnow reg
.c6          equ    40           ; repeated (64bit) to fill 3dnow reg
.c12         equ    48           ; repeated (64bit) to fill 3dnow reg
.two         equ    56           ; repeated (64bit) to fill 3dnow reg
.n1          equ    64           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    72           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    80	
.innerjjnr0  equ    84
.innerk0     equ    88		
.innerjjnr   equ    92
.innerk      equ    96	
.fix         equ    100
.fiy         equ    104
.fiz	     equ    108
.dx1	     equ    112
.dy1	     equ    116
.dz1	     equ    120
.dx2	     equ    124
.dy2	     equ    128
.dz2	     equ    132								
.nsvdwc      equ    136
.nscoul      equ    140
.nsvdw       equ    144
.solnr	     equ    148		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 152		;  local stack space
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

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
	pxor  mm7,mm7
	movq  [esp + .vnbtot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	
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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	; dispersion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
		
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdwc_inner
	jmp  .updateouterdata_vdwc		
.single_vdwc_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 3
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdwc:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	; dispersion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	shl ecx, 3
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdw_inner
	jmp  .updateouterdata_vdw		
.single_vdw_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 3
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdw:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 152
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl1000_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.iq          equ    20		; repeated (64bit) to fill 3dnow reg
.vctot       equ    28           ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    36
.innerk      equ    40		
.fix         equ    44
.fiy         equ    48
.fiz	     equ    52
.dx1	     equ    56
.dy1	     equ    60
.dz1	     equ    64
.dx2	     equ    68
.dy2	     equ    72
.dz2	     equ    76									
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 80		;  80 bytes local stack space
	femms
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
	movd mm7, [ecx + ebx*4]  	 ; charge[jnr2] 
	punpckldq mm3,mm7	         ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	prefetchw [esp + .dx1]	;  prefetch i forces to cache
	
	pfmul mm3,mm1		;  6 has both vcoul
	pfmul mm0,mm3		;  0 has both fscal 

	;; update vctot

	pfadd mm3, [esp + .vctot]      ; add the earlier value 
	movq [esp + .vctot], mm3       ;  store the sum
	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]
	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm6, [esp + .iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;  mm6=qq
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + .dx1], mm0
	pfmul mm0,mm0
	movd  [esp + .dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;  mm0=rsq
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar
	;;  update vctot
	movq mm5, [esp + .vctot]
	pfadd mm5, mm6
	movq [esp + .vctot], mm5
	;;  spread fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dx1]
	movd mm1,  [esp + .dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;; update i particle force
	movq mm2,  [esp + .fix]
	movd mm3,  [esp + .fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fix], mm2
	movd [esp + .fiz], mm3
	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 80
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl1010_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.iq          equ    32		; repeated (64bit) to fill 3dnow reg
.vctot       equ    40           ; repeated (64bit) to fill 3dnow reg
.innerjjnr0  equ    48
.innerk0     equ    52		
.innerjjnr   equ    56
.innerk      equ    60		
.fix         equ    64
.fiy         equ    68
.fiz	     equ    72
.dx1	     equ    76
.dy1	     equ    80
.dz1	     equ    84
.dx2	     equ    88
.dy2	     equ    92
.dz2	     equ    96									
.nscoul      equ    100
.solnr	     equ    104		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		;  local stack space
	femms
	;; assume we have at least one i particle - start directly	
	add   [ebp + %$nsatoms], dword 8

.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	mov   ecx, [eax]
	add   [ebp + %$nsatoms], dword 12
	mov   [esp + .nscoul], ecx
		
	;; clear potential
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]

	mov   ecx, [esp + .nscoul]
	cmp   ecx, dword 0
	jnz   .mno_coul
	jmp   .last_mno
.mno_coul:				
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner
.unroll_coul_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
	movd mm7, [ecx + ebx*4]  	 ; charge[jnr2] 
	punpckldq mm3,mm7	         ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	prefetchw [esp + .dx1]	;  prefetch i forces to cache
	
	pfmul mm3,mm1		;  6 has both vcoul
	pfmul mm0,mm3		;  0 has both fscal 

	;; update vctot

	pfadd mm3, [esp + .vctot]      ; add the earlier value 
	movq [esp + .vctot], mm3       ;  store the sum
	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]
	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_coul_inner
	jmp  .updateouterdata_coul		
.single_coul_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm6, [esp + .iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;  mm6=qq
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + .dx1], mm0
	pfmul mm0,mm0
	movd  [esp + .dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;  mm0=rsq
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar
	;;  update vctot
	movq mm5, [esp + .vctot]
	pfadd mm5, mm6
	movq [esp + .vctot], mm5
	;;  spread fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dx1]
	movd mm1,  [esp + .dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;; update i particle force
	movq mm2,  [esp + .fix]
	movd mm3,  [esp + .fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fix], mm2
	movd [esp + .fiz], mm3
	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;;  done!
.updateouterdata_coul:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .last_mno
	jmp .mno_coul
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

			
proc inl1020_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.iqO         equ    44		; repeated (64bit) to fill 3dnow reg
.iqH         equ    52		; repeated (64bit) to fill 3dnow reg
.vctot       equ    60          ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    68
.innerk      equ    72		
.fixO        equ    76 
.fiyO        equ    80
.fizO        equ    84
.fixH        equ    88          ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    96          ; repeated (64bit) to fill 3dnow reg
.fizH        equ    104         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    112
.dyO	     equ    116
.dzO	     equ    120
.dxH	     equ    124         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    132         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    140         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 148		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0]
	pfmul mm2, mm1		
	movq  [esp + .iqO], mm2	        ;  iqO = facel*charge[ii]
	
	movd  mm2, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1]
	pfmul mm2, mm1
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iqH], mm2	        ;  iqH = facel*charge[ii0+1]
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	
					
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .fixO],   mm7
	movd  [esp + .fizO],   mm7
	movq  [esp + .fixH],   mm7
	movq  [esp + .fiyH],   mm7
	movq  [esp + .fizH],   mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 
	;prefetch [ecx + 16]	   ; prefetch data - trial and error says 16 is best

	mov ecx, [ebp + %$charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + .iqO]
	pfmul mm7, [esp + .iqH]	;  mm6=qqO, mm7=qqH
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2
	
	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 148
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl1030_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.qqOO        equ    44		; repeated (64bit) to fill 3dnow reg
.qqOH        equ    52		; repeated (64bit) to fill 3dnow reg
.qqHH        equ    60     	; repeated (64bit) to fill 3dnow reg
.vctot       equ    68          ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    76
.innerk      equ    80		
.fixO        equ    84 
.fiyO        equ    88
.fizO        equ    92
.fixH        equ    96          ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    104         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    112         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    120
.dyO	     equ    124
.dzO	     equ    128
.dxH	     equ    132         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    140         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    148         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 156		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]	;  mm1=facel
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0] (O)
	movd  mm3, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1] (H)
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			; mm4=qqOO*facel
	pfmul mm5, mm3			; mm5=qqOH*facel
	pfmul mm6, mm3			; mm6=qqHH*facel
	punpckldq mm5,mm5	        ;  spread to both halves
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .qqOO], mm4
	movq  [esp + .qqOH], mm5
	movq  [esp + .qqHH], mm6
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	

	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .fixO],  mm7
	movq  [esp + .fizO],  mm7
	movq  [esp + .fixH],  mm7
	movq  [esp + .fiyH],  mm7
	movq  [esp + .fizH],  mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 

	movd  mm6, [esp + .qqOO]
	movq  mm7, [esp + .qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	; interactions with j H1.
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + .qqOH]
	movq mm7, [esp + .qqHH]
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	; interactions with j H2
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + .qqOH]
	movq mm7, [esp + .qqHH]

	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4	

	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop	
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2

	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 156
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl1100_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.iq    	     equ    20		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    28           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    36           ; repeated (64bit) to fill 3dnow reg
.c6          equ    44           ; repeated (64bit) to fill 3dnow reg
.c12         equ    52           ; repeated (64bit) to fill 3dnow reg
.six         equ    60           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    68           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    76
.innerjjnr   equ    80
.innerk      equ    84		
.fix         equ    88
.fiy         equ    92
.fiz	     equ    96
.dx1	     equ    100
.dy1	     equ    104
.dz1	     equ    108
.dx2	     equ    112
.dy2	     equ    116
.dz2	     equ    120						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 124		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + .six],    mm0
	movq  [esp + .twelve], mm1
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm3, mm1		;  mm3 has vcoul for both interactions
	movq  mm7, mm3	        ;  use mm7 for sum to make fscal

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;; update vctot
	pfadd mm3, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm3       ;  store the sum       

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5


	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	movq  mm1, mm0
	pfmul mm0, mm0		;  mm0=invsq
	;;  calculate potentials and scalar force
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm3, mm1		;  mm3 has vcoul for both interactions
	movq  mm7, mm3	        ;  use mm7 for sum to make fscal

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        ; mm0 is total fscal now

	;;  update vctot
	pfadd mm3, [esp + .vctot]
	movq [esp + .vctot], mm3

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 124
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	



proc inl1110_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.iq          equ    32		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    40           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    48           ; repeated (64bit) to fill 3dnow reg
.c6          equ    56           ; repeated (64bit) to fill 3dnow reg
.c12         equ    64           ; repeated (64bit) to fill 3dnow reg
.six         equ    72           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    80           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    88	
.innerjjnr0  equ    92
.innerk0     equ    96		
.innerjjnr   equ    100
.innerk      equ    104	
.fix         equ    108
.fiy         equ    112
.fiz	     equ    116
.dx1	     equ    120
.dy1	     equ    124
.dz1	     equ    128
.dx2	     equ    132
.dy2	     equ    136
.dz2	     equ    140								
.nsvdwc      equ    144
.nscoul      equ    148
.nsvdw       equ    152
.solnr	     equ    156		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 160		;  local stack space
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + .six],    mm0
	movq  [esp + .twelve], mm1
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

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
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	
	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx	

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm3, mm1		;  mm3 has vcoul for both interactions
	movq  mm7, mm3	        ;  use mm7 for sum to make fscal

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;; update vctot
	pfadd mm3, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm3       ;  store the sum       

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdwc_inner
	jmp  .updateouterdata_vdwc		
.single_vdwc_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5


	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	movq  mm1, mm0
	pfmul mm0, mm0		;  mm0=invsq
	;;  calculate potentials and scalar force
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm3, mm1		;  mm3 has vcoul for both interactions
	movq  mm7, mm3	        ;  use mm7 for sum to make fscal

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7        ; mm0 is total fscal now

	;;  update vctot
	pfadd mm3, [esp + .vctot]
	movq [esp + .vctot], mm3

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdwc:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner
.unroll_coul_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
	movd mm7, [ecx + ebx*4]  	 ; charge[jnr2] 
	punpckldq mm3,mm7	         ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	prefetchw [esp + .dx1]	;  prefetch i forces to cache
	
	pfmul mm3,mm1		;  6 has both vcoul
	pfmul mm0,mm3		;  0 has both fscal 

	;; update vctot

	pfadd mm3, [esp + .vctot]      ; add the earlier value 
	movq [esp + .vctot], mm3       ;  store the sum
	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]
	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_coul_inner
	jmp  .updateouterdata_coul		
.single_coul_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm6, [esp + .iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;  mm6=qq
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + .dx1], mm0
	pfmul mm0,mm0
	movd  [esp + .dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;  mm0=rsq
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar
	;;  update vctot
	movq mm5, [esp + .vctot]
	pfadd mm5, mm6
	movq [esp + .vctot], mm5
	;;  spread fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dx1]
	movd mm1,  [esp + .dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;; update i particle force
	movq mm2,  [esp + .fix]
	movd mm3,  [esp + .fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fix], mm2
	movd [esp + .fiz], mm3
	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;;  done!
.updateouterdata_coul:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw
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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best	

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;; mm0 now contains invsq, and mm1 invsqrt
	;; do potential and fscal
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdw_inner
	jmp  .updateouterdata_vdw		
.single_vdw_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5


	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	movq  mm1, mm0
	pfmul mm0, mm0		;  mm0=invsq
	;;  calculate potentials and scalar force
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7        ; mm0 is total fscal now

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdw:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 160
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



proc inl1120_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.iqO         equ    44		; repeated (64bit) to fill 3dnow reg
.iqH         equ    52		; repeated (64bit) to fill 3dnow reg
.vctot       equ    60          ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    68          ; repeated (64bit) to fill 3dnow reg
.c6          equ    76          ; repeated (64bit) to fill 3dnow reg
.c12         equ    84          ; repeated (64bit) to fill 3dnow reg
.six         equ    92          ; repeated (64bit) to fill 3dnow reg
.twelve      equ    100         ; repeated (64bit) to fill 3dnow reg
.ntia        equ    108         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    116
.innerk      equ    120	
.fixO        equ    124
.fiyO        equ    128
.fizO        equ    132
.fixH        equ    136         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    144         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    152         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    160
.dyO	     equ    164
.dzO	     equ    168
.dxH	     equ    172         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    180         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    188         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0]
	pfmul mm2, mm1		
	movq  [esp + .iqO], mm2	        ;  iqO = facel*charge[ii]
	
	movd  mm2, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1]
	pfmul mm2, mm1
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iqH], mm2	        ;  iqH = facel*charge[ii0+1]

	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx
	
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + .six],    mm3
	movq  [esp + .twelve], mm4  
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	
					
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],   mm7
	movd  [esp + .fizO],   mm7
	movq  [esp + .fixH],   mm7
	movq  [esp + .fiyH],   mm7
	movq  [esp + .fizH],   mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 
	;prefetch [ecx + 16]	   ; prefetch data - trial and error says 16 is best

	mov ecx, [ebp + %$charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + .iqO]
	pfmul mm7, [esp + .iqH]	;  mm6=qqO, mm7=qqH

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr]
	mov ecx, [ebp + %$nbfp]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [ecx + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [ecx + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5	
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		;  mm0=rinvsix
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm2=rintwelve
	
	;;  calculate potential and scalar force
	pfmul mm6, mm1		;  mm6=vcoul
	movq  mm1, mm6		;  use mm1 for fscal sum

	;; LJ for the oxygen
	pfmul mm0, [esp + .c6]	 
	pfmul mm2, [esp + .c12]	 

	;; calc nb potential
	movq mm5, mm2
	pfsub mm5, mm0

	;; calc nb force
	pfmul mm0, [esp + .six]
	pfmul mm2, [esp + .twelve]
	
	;; increment scalar force
	pfsub mm1, mm0
	pfadd mm1, mm2
	pfmul mm4, mm1		;  total scalar force on oxygen.
	
	;;  update nb potential
	pfadd mm5, [esp + .vnbtot]
	movq [esp + .vnbtot], mm5
	
	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2
	
	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  same for Vnb.
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	

proc inl1130_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.qqOO        equ    44		; repeated (64bit) to fill 3dnow reg
.qqOH        equ    52		; repeated (64bit) to fill 3dnow reg
.qqHH        equ    60     	; repeated (64bit) to fill 3dnow reg
.c6          equ    68     	; repeated (64bit) to fill 3dnow reg
.c12         equ    76     	; repeated (64bit) to fill 3dnow reg
.six         equ    84     	; repeated (64bit) to fill 3dnow reg
.twelve	     equ    92     	; repeated (64bit) to fill 3dnow reg
.vctot       equ    100         ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    108         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    116
.innerk      equ    120	
.fixO        equ    124
.fiyO        equ    128
.fizO        equ    132
.fixH        equ    136         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    144         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    152         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    160
.dyO	     equ    164
.dzO	     equ    168
.dxH	     equ    172         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    180         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    188         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]	;  mm1=facel
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0] (O)
	movd  mm3, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1] (H)
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			; mm4=qqOO*facel
	pfmul mm5, mm3			; mm5=qqOH*facel
	pfmul mm6, mm3			; mm6=qqHH*facel
	punpckldq mm5,mm5	        ;  spread to both halves
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .qqOO], mm4
	movq  [esp + .qqOH], mm5
	movq  [esp + .qqHH], mm6
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + .c6], mm0
	movq  [esp + .c12], mm1
	movq  mm2, [mm_six]
	movq  mm3, [mm_twelve]
	movq  [esp + .six], mm2
	movq  [esp + .twelve], mm3
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	

	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],  mm7
	movq  [esp + .fizO],  mm7
	movq  [esp + .fixH],  mm7
	movq  [esp + .fiyH],  mm7
	movq  [esp + .fizH],  mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 

	movd  mm6, [esp + .qqOO]
	movq  mm7, [esp + .qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt OO
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq OO

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + .c6]
	pfmul mm0, [esp + .c12]
	movq mm5, mm0
	pfsub mm5, mm2		; vnb

	pfmul mm2, [esp + .six]
	pfmul mm0, [esp + .twelve]

	pfsub mm0, mm2
	
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfadd mm0, mm6
	pfmul mm4, mm0		;  mm4=fscalar

	;;  update nb potential
	pfadd mm5, [esp + .vnbtot]
	movq [esp + .vnbtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	; interactions with j H1.
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + .qqOH]
	movq mm7, [esp + .qqHH]
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	; interactions with j H2
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + .qqOH]
	movq mm7, [esp + .qqHH]

	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	movq  mm4, mm1
	pfmul mm4, mm4		;  mm4=invsq
	;;  calculate potential and scalar force
	pfmul mm6, mm1		; mm6=vcoul
	pfmul mm4, mm6		;  mm4=fscalar

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;  seeds are in mm5 now, and rsq in mm3.

	movq mm2, mm5
	pfmul mm5,mm5
        pfrsqit1 mm5,mm3				
        pfrcpit2 mm5,mm2	;  mm5=invsqrt
	movq mm3,mm5
	pfmul mm3,mm3		; mm3=invsq
	pfmul mm7, mm5		;  mm7=vcoul
	pfmul mm3, mm7		;  mm3=fscal for the two H's.

	;;  update vctot
	pfadd mm7, mm6
	pfadd mm7, [esp + .vctot]
	movq [esp + .vctot], mm7
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm4,mm4
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4	

	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop	
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2

	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnbtot[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl3000_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.iq    	     equ    20		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    28           ; repeated (64bit) to fill 3dnow reg
.two         equ    36           ; repeated (64bit) to fill 3dnow reg
.n1          equ    44           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    52           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    60
.innerjjnr   equ    64
.innerk      equ    68		
.fix         equ    72
.fiy         equ    76
.fiz	     equ    80
.dx1	     equ    84
.dy1	     equ    88
.dz1	     equ    92
.dx2	     equ    96
.dy2	     equ    100
.dz2	     equ    104						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 108		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + .tabscale]	
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 108
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc



	
proc inl3010_3dnow
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
%$VFtab		arg
%$nsatoms       arg			
	;; stack offsets for local variables
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.iq          equ    32		; repeated (64bit) to fill 3dnow reg
.vctot       equ    40          ; repeated (64bit) to fill 3dnow reg
.two         equ    48		; repeated (64bit) to fill 3dnow reg
.n1          equ    56		; repeated (64bit) to fill 3dnow reg
.tabscale    equ    64		; repeated (64bit) to fill 3dnow reg			
.innerjjnr0  equ    72
.innerk0     equ    76		
.innerjjnr   equ    80
.innerk      equ    84		
.fix         equ    88
.fiy         equ    92
.fiz	     equ    96
.dx1	     equ    100
.dy1	     equ    104
.dz1	     equ    108
.dx2	     equ    112
.dy2	     equ    116
.dz2	     equ    120								
.nscoul      equ    124
.solnr	     equ    128		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 132		;  local stack space
	femms
	
	add   [ebp + %$nsatoms], dword 8
	movq  mm2, [mm_two]
	movq  [esp + .two], mm2
	movd  mm3, [ebp + %$tabscale]
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3
	
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   eax, [ebp + %$nsatoms]
	mov   ecx, [eax]
	add   [ebp + %$nsatoms], dword 12
	mov   [esp + .nscoul], ecx
		
	;; clear potential
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	mov   ecx, [esp + .nscoul]
	cmp   ecx, byte 0
	jnz  .mno_coul
	jmp  .last_mno
.mno_coul:				
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner
.unroll_coul_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + .tabscale]	
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_coul_inner
	jmp  .updateouterdata_coul		
.single_coul_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_coul:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .last_mno
	jmp .mno_coul
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 132
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	


proc inl3020_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.iqO         equ    44		; repeated (64bit) to fill 3dnow reg
.iqH         equ    52		; repeated (64bit) to fill 3dnow reg
.qqO         equ    60		; repeated (64bit) to fill 3dnow reg
.qqH         equ    68		; repeated (64bit) to fill 3dnow reg
.vctot       equ    76          ; repeated (64bit) to fill 3dnow reg
.two         equ    84          ; repeated (64bit) to fill 3dnow reg
.n1          equ    92          ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    100          ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    108
.innerk      equ    112	
.fixO        equ    116
.fiyO        equ    120
.fizO        equ    124
.fixH        equ    128         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    136         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    144         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    152
.dyO	     equ    156
.dzO	     equ    160
.dxH	     equ    164         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    172         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    180         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    188        	; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 196		;   local stack space
	femms

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0]
	pfmul mm2, mm1		
	movq  [esp + .iqO], mm2	        ;  iqO = facel*charge[ii]
	
	movd  mm2, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1]
	pfmul mm2, mm1
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iqH], mm2	        ;  iqH = facel*charge[ii0+1]

	movq  mm3, [mm_two]
	movd  mm4, [ebp + %$tabscale]
	punpckldq mm4,mm4	        ;  spread to both halves
	movq  [esp + .two],    mm3
	movq  [esp + .tabscale], mm4	      
	;; assume we have at least one i particle - start directly	 
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	
					
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .fixO],   mm7
	movd  [esp + .fizO],   mm7
	movq  [esp + .fixH],   mm7
	movq  [esp + .fiyH],   mm7
	movq  [esp + .fizH],   mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:	
	;; a single j particle iteration.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 
	;prefetch [ecx + 16]	   ; prefetch data - trial and error says 16 is best

	mov ecx, [ebp + %$charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + .iqO]
	pfmul mm7, [esp + .iqH]	;  mm6=qqO, mm7=qqH
	movd [esp + .qqO], mm6
	movq [esp + .qqH], mm7
		
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt

	pfmul mm0, mm1		;  mm0=r

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqO]	; fijC=qq*FF
	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	movq mm3, mm7	

	; change sign of fscal and multiply with rinv
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now		
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; now do the two hydrogens.
	movq mm0, [esp + .tmprsqH] ; mm0=rsqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqH]	; fijC=qq*FF

	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 + 8], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3, mm3	
	;; mm1 is fzH1
	;; mm3 is fzH2
	
	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 196
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	

proc inl3030_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.qqOO        equ    44		; repeated (64bit) to fill 3dnow reg
.qqOH        equ    52		; repeated (64bit) to fill 3dnow reg
.qqHH        equ    60     	; repeated (64bit) to fill 3dnow reg
.two         equ    68     	; repeated (64bit) to fill 3dnow reg
.n1	     equ    76     	; repeated (64bit) to fill 3dnow reg
.tabscale    equ    84     	; repeated (64bit) to fill 3dnow reg
.vctot       equ    92          ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    100
.innerk      equ    104	
.fixO        equ    108
.fiyO        equ    112
.fizO        equ    116
.fixH        equ    120         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    128         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    136         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    144
.dyO	     equ    148
.dzO	     equ    152
.dxH	     equ    156         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    164         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    172         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    180         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 188		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]	;  mm1=facel
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0] (O)
	movd  mm3, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1] (H)
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			; mm4=qqOO*facel
	pfmul mm5, mm3			; mm5=qqOH*facel
	pfmul mm6, mm3			; mm6=qqHH*facel
	punpckldq mm5,mm5	        ;  spread to both halves
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .qqOO], mm4
	movq  [esp + .qqOH], mm5
	movq  [esp + .qqHH], mm6
	movq  mm2, [mm_two]
	movq  [esp + .two], mm2
	movd  mm3, [ebp + %$tabscale]
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	

	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .fixO],  mm7
	movq  [esp + .fizO],  mm7
	movq  [esp + .fixH],  mm7
	movq  [esp + .fiyH],  mm7
	movq  [esp + .fizH],  mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt OO
	pfmul mm0, mm1		;  mm0=rsq OO

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOO]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	movq mm3, mm7

	; change sign of fscal and multiply with rinv
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; time for hydrogens!

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	; interactions with j H1.

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1		;  mm0=rsq 
	
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, force is moved to mm3.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now		

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	; interactions with j H2
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4	

	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop	
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2

	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 188
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc




proc inl3100_3dnow
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
.is3         equ     0 
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.iq    	     equ    20		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    28           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    36           ; repeated (64bit) to fill 3dnow reg
.c6          equ    44           ; repeated (64bit) to fill 3dnow reg
.c12         equ    52           ; repeated (64bit) to fill 3dnow reg
.six         equ    60           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    68           ; repeated (64bit) to fill 3dnow reg
.two         equ    76           ; repeated (64bit) to fill 3dnow reg
.n1          equ    84           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    92           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    100
.innerjjnr   equ    104
.innerk      equ    108	
.fix         equ    112
.fiy         equ    116
.fiz	     equ    120
.dx1	     equ    124
.dy1	     equ    128
.dz1	     equ    132
.dx2	     equ    136
.dy2	     equ    140
.dz2	     equ    144						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 148		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_two]
	movq  mm1, [mm_six]
	movq  mm2, [mm_twelve]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	movq  [esp + .six],    mm1
	movq  [esp + .twelve],    mm2
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	movq mm1, mm0
	pfmul mm1,mm1 	; mm1=invsq
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	; mm2=rinvsix
	movq  mm1,mm2
	pfmul mm1,mm1	; mm1=rinvtwelve
	
	pfmul mm3, [esp + .tabscale]
	
	pfmul mm1, [esp + .c12]

	pfmul mm2, [esp + .c6]

	movq mm4, mm1
	pfsub mm4, mm2	; mm4 = vnb12-vnb6

	pfmul mm2, [esp + .six]
	pfmul mm1, [esp + .twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	; mm1=	(12*vnb12-6*vnb6)*rinv11

	pfsub mm1, mm3

	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm4, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm4       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5


	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq
	
	; at this point mm5 contains vcoul and mm3 fijC.

	movq mm1, mm0
	pfmul mm1,mm1 	; mm1=invsq
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	; mm2=rinvsix
	movq  mm1,mm2
	pfmul mm1,mm1	; mm1=rinvtwelve
	
	pfmul mm3, [esp + .tabscale]
	
	pfmul mm1, [esp + .c12]

	pfmul mm2, [esp + .c6]

	movq mm4, mm1
	pfsub mm4, mm2	; mm4 = vnb12-vnb6

	pfmul mm2, [esp + .six]
	pfmul mm1, [esp + .twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	; mm1=	(12*vnb12-6*vnb6)*rinv11

	pfsub mm1, mm3

	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm4, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm4       ;  store the sum       

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix] 
	pfadd mm7, [esp + .fiz] 
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
 
	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 148
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc







proc inl3110_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.iq          equ    32		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    40           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    48           ; repeated (64bit) to fill 3dnow reg
.c6          equ    56           ; repeated (64bit) to fill 3dnow reg
.c12         equ    64           ; repeated (64bit) to fill 3dnow reg
.six         equ    72           ; repeated (64bit) to fill 3dnow reg
.twelve      equ    80           ; repeated (64bit) to fill 3dnow reg
.two         equ    88           ; repeated (64bit) to fill 3dnow reg
.n1          equ    96           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    104           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    112
.innerjjnr0  equ    116
.innerk0     equ    120	
.innerjjnr   equ    124
.innerk      equ    128	
.fix         equ    132
.fiy         equ    136
.fiz	     equ    140
.dx1	     equ    144
.dy1	     equ    148
.dz1	     equ    152
.dx2	     equ    156
.dy2	     equ    160
.dz2	     equ    164								
.nsvdwc      equ    168
.nscoul      equ    172
.nsvdw       equ    176
.solnr	     equ    180		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 184		;  local stack space
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + .six],    mm0
	movq  [esp + .twelve], mm1
	movq  mm2, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm2
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

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
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	
	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx	

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	movq mm1, mm0
	pfmul mm1,mm1 	; mm1=invsq
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	; mm2=rinvsix
	movq  mm1,mm2
	pfmul mm1,mm1	; mm1=rinvtwelve
	
	pfmul mm3, [esp + .tabscale]
	
	pfmul mm1, [esp + .c12]

	pfmul mm2, [esp + .c6]

	movq mm4, mm1
	pfsub mm4, mm2	; mm4 = vnb12-vnb6

	pfmul mm2, [esp + .six]
	pfmul mm1, [esp + .twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	; mm1=	(12*vnb12-6*vnb6)*rinv11

	pfsub mm1, mm3

	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm4, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm4       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7	
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdwc_inner
	jmp  .updateouterdata_vdwc		
.single_vdwc_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5


	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	movq mm1, mm0
	pfmul mm1,mm1 	; mm1=invsq
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	; mm2=rinvsix
	movq  mm1,mm2
	pfmul mm1,mm1	; mm1=rinvtwelve
	
	pfmul mm3, [esp + .tabscale]
	
	pfmul mm1, [esp + .c12]

	pfmul mm2, [esp + .c6]

	movq mm4, mm1
	pfsub mm4, mm2	; mm4 = vnb12-vnb6

	pfmul mm2, [esp + .six]
	pfmul mm1, [esp + .twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	; mm1=	(12*vnb12-6*vnb6)*rinv11

	pfsub mm1, mm3

	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm4, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm4       ;  store the sum       

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdwc:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner
.unroll_coul_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + .tabscale]	
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_coul_inner
	jmp  .updateouterdata_coul		
.single_coul_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_coul:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw
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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best	

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrcp mm0, mm4	                 ; lookup reciprocal seed 
        pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
                  	        	; amd 3dnow N-R iteration to get full precision.
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	
	;; mm4 now contains invsq,
	;; do potential and fscal
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdw_inner
	jmp  .updateouterdata_vdw		
.single_vdw_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm4=rsq
	
        pfrcp mm0,mm4
        pfrcpit1 mm4,mm0				
        pfrcpit2 mm4,mm0	;  mm4=invsq 
	;;  calculate potentials and scalar force
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	; mm4=rinvsix
	movq  mm5, mm4	
	pfmul mm5, mm5	                ; mm5=rinvtwelve

	pfmul mm5, [esp + .c12]
	pfmul mm4, [esp + .c6]	
	movq mm6, mm5	; mm6 is vnb12-vnb6
	pfsub mm6, mm4

	pfmul mm4, [esp + .six]

	pfmul mm5, [esp + .twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5        ; mm0 is total fscal now

	;; update vnbtot
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdw:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 184
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


	

proc inl3120_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.iqO         equ    44		; repeated (64bit) to fill 3dnow reg
.iqH         equ    52		; repeated (64bit) to fill 3dnow reg
.qqO         equ    60		; repeated (64bit) to fill 3dnow reg
.qqH         equ    68		; repeated (64bit) to fill 3dnow reg
.vctot       equ    76          ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    84          ; repeated (64bit) to fill 3dnow reg
.c6          equ    92          ; repeated (64bit) to fill 3dnow reg
.c12         equ    100         ; repeated (64bit) to fill 3dnow reg
.six         equ    108         ; repeated (64bit) to fill 3dnow reg
.twelve      equ    116         ; repeated (64bit) to fill 3dnow reg
.two         equ    124         ; repeated (64bit) to fill 3dnow reg
.n1          equ    132         ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    140         ; repeated (64bit) to fill 3dnow reg
.ntia        equ    148         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    156
.innerk      equ    160	
.fixO        equ    164
.fiyO        equ    168
.fizO        equ    172
.fixH        equ    176         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    184         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    192         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    200
.dyO	     equ    204
.dzO	     equ    208
.dxH	     equ    212         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    220         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    228         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    236        	; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 244		;   local stack space
	femms

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0]
	pfmul mm2, mm1
	movq  [esp + .iqO], mm2	        ;  iqO = facel*charge[ii]
	
	movd  mm2, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1]
	pfmul mm2, mm1
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iqH], mm2	        ;  iqH = facel*charge[ii0+1]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [mm_six]
	movq  mm5, [mm_twelve]
	movq  mm6, [ebp + %$tabscale]
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .two], mm3
	movq  [esp + .six], mm4
	movq  [esp + .twelve], mm5
	movq  [esp + .tabscale], mm6	      
 	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	
					
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],   mm7
	movd  [esp + .fizO],   mm7
	movq  [esp + .fixH],   mm7
	movq  [esp + .fiyH],   mm7
	movq  [esp + .fizH],   mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 
	;prefetch [ecx + 16]	   ; prefetch data - trial and error says 16 is best

	mov ecx, [ebp + %$charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + .iqO]
	pfmul mm7, [esp + .iqH]	;  mm6=qqO, mm7=qqH
	movd [esp + .qqO], mm6
	movq [esp + .qqH], mm7

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr]
	mov ecx, [ebp + %$nbfp]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [ecx + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [ecx + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt

	pfmul mm0, mm1		;  mm0=r

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqO]	; fijC=qq*FF
	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5

	movq mm3, mm7
	pfmul mm3, [esp + .tabscale]
	
	; nontabulated LJ - mm1 is invsqrt. - keep mm1!
	movq mm0, mm1
	pfmul mm0, mm0		;  mm0 is invsq
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		; mm2 = rinvsix
	movq mm4, mm2
	pfmul mm4, mm4		;  mm4=rinvtwelve

	pfmul mm4, [esp + .c12]
	pfmul mm2, [esp + .c6]
	movq mm5, mm4
	pfsub mm5, mm2		; mm5=vnb12-vnb6

	pfmul mm2, [esp + .six]
	pfmul mm4, [esp + .twelve]
	pfsub mm4, mm2
	pfmul mm4, mm1        ; mm4=(12*vnb12-6*vnb6)*rinv11

	pfsubr mm3, mm4 
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	
	
	;; update vnbtot 
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; now do the two hydrogens.
	movq mm0, [esp + .tmprsqH] ; mm0=rsqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7
	pfmul mm4, [esp + .tabscale]	
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2
	
	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  same for Vnb.
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 244
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc





proc inl3130_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.qqOO        equ    44		; repeated (64bit) to fill 3dnow reg
.qqOH        equ    52		; repeated (64bit) to fill 3dnow reg
.qqHH        equ    60     	; repeated (64bit) to fill 3dnow reg
.c6          equ    68     	; repeated (64bit) to fill 3dnow reg
.c12         equ    76     	; repeated (64bit) to fill 3dnow reg
.six         equ    84     	; repeated (64bit) to fill 3dnow reg
.twelve      equ    92     	; repeated (64bit) to fill 3dnow reg
.two         equ    100     	; repeated (64bit) to fill 3dnow reg
.n1	     equ    108    	; repeated (64bit) to fill 3dnow reg
.tabscale    equ    116    	; repeated (64bit) to fill 3dnow reg
.vctot       equ    124         ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    132         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    140
.innerk      equ    144	
.fixO        equ    148
.fiyO        equ    152
.fizO        equ    156
.fixH        equ    160         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    168         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    176         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    184
.dyO	     equ    188
.dzO	     equ    192
.dxH	     equ    200         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    208         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    216         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    224         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 232		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]	;  mm1=facel
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0] (O)
	movd  mm3, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1] (H)
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			; mm4=qqOO*facel
	pfmul mm5, mm3			; mm5=qqOH*facel
	pfmul mm6, mm3			; mm6=qqHH*facel
	punpckldq mm5,mm5	        ;  spread to both halves
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .qqOO], mm4
	movq  [esp + .qqOH], mm5
	movq  [esp + .qqHH], mm6
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + .c6], mm0
	movq  [esp + .c12], mm1
	movq  mm2, [mm_two]
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + .two], mm2
	movq  [esp + .six], mm3
	movq  [esp + .twelve], mm4
	movd  mm5, [ebp + %$tabscale]
	punpckldq mm5,mm5
	movq  [esp + .tabscale], mm5
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	

	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],  mm7
	movq  [esp + .fizO],  mm7
	movq  [esp + .fixH],  mm7
	movq  [esp + .fiyH],  mm7
	movq  [esp + .fizH],  mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt OO
	pfmul mm0, mm1		;  mm0=rsq OO

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOO]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	movq mm3, mm7
	pfmul mm3, [esp + .tabscale]
	
	movq mm5, mm1
	pfmul mm5,mm5
	movq mm4, mm5
	pfmul mm4,mm5
	pfmul mm4,mm5
	movq mm5, mm4
	pfmul mm5,mm5	; mm4=rinvsix, mm5=rinvtwelve

	pfmul mm4, [esp + .c6]
	pfmul mm5, [esp + .c12]
	movq mm6,mm5
	pfsub mm6,mm4

	pfmul mm4, [esp + .six]
	pfmul mm5, [esp + .twelve]
	pfsub mm5,mm4
	pfmul mm5, mm1
	pfsubr mm3, mm5

 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	
	
	;; update vnbtot 
	pfadd mm6, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm6       ;  store the sum       
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; time for hydrogens!

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	; interactions with j H1.

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1		;  mm0=rsq 
	
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, force is moved to mm3.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now		

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	; interactions with j H2
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4	

	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop	
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2

	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnbtot[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 232
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl3300_3dnow
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
.is3         equ     0
.ii3         equ     4
.ix          equ     8
.iy          equ    12
.iz          equ    16
.iq    	     equ    20		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    28           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    36           ; repeated (64bit) to fill 3dnow reg
.c6          equ    44           ; repeated (64bit) to fill 3dnow reg
.c12         equ    52           ; repeated (64bit) to fill 3dnow reg
.two         equ    60           ; repeated (64bit) to fill 3dnow reg
.n1          equ    68           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    76           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    84
.innerjjnr   equ    88
.innerk      equ    92		
.fix         equ    96
.fiy         equ    100
.fiz	     equ    104
.dx1	     equ    108
.dy1	     equ    112
.dz1	     equ    116
.dx2	     equ    120
.dy2	     equ    124
.dz2	     equ    128						
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 132		;   local stack space
	femms
	; move data to local stack 
	movq  mm0, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly	
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm0, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm3, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx
	pfadd mm1, mm3
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	
				
	;; clear total potential and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fix],    mm7
	movd  [esp + .fiz],    mm7

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
	sub   edx, dword 2
	mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_loop
	jmp   .finish_inner
.unroll_loop:
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; dispersion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_inner
	jmp   .unroll_loop
.finish_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_inner
	jmp  .updateouterdata		
.single_inner:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]

	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 132
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc





proc inl3310_3dnow
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
.is3         equ     0
.ii3         equ     4
.shX	     equ     8
.shY         equ    12 
.shZ	     equ    16	
.ix          equ    20
.iy          equ    24
.iz          equ    28	
.iq          equ    32		 ; repeated (64bit) to fill 3dnow reg
.vctot       equ    40           ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    48           ; repeated (64bit) to fill 3dnow reg
.c6          equ    56           ; repeated (64bit) to fill 3dnow reg
.c12         equ    64           ; repeated (64bit) to fill 3dnow reg
.two         equ    72           ; repeated (64bit) to fill 3dnow reg
.n1          equ    80           ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    88           ; repeated (64bit) to fill 3dnow reg
.ntia	     equ    96	
.innerjjnr0  equ    100
.innerk0     equ    104	
.innerjjnr   equ    108
.innerk      equ    112	
.fix         equ    116
.fiy         equ    120
.fiz	     equ    124
.dx1	     equ    128
.dy1	     equ    132
.dz1	     equ    136
.dx2	     equ    140
.dy2	     equ    144
.dz2	     equ    148								
.nsvdwc      equ    152
.nscoul      equ    156
.nsvdw       equ    160
.solnr	     equ    164		
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 168		;  local stack space
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + %$tabscale]
	movq  [esp + .two],    mm0
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3	
	;; assume we have at least one i particle - start directly		
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;  move shX/shY to mm0 and shZ to mm1.
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + .shX], mm0
	movd  [esp + .shZ], mm1

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
	pxor  mm7,mm7
	movq  [esp + .vctot],  mm7
	movq  [esp + .vnbtot], mm7
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
	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]
	
	mov   ecx, [esp + .nsvdwc]
	cmp   ecx, dword 0
	jnz   .mno_vdwc
	jmp   .testcoul
.mno_vdwc:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	mov   edx, [ebp + %$type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + %$ntype]
	shl   edx, 1
	mov   [esp + .ntia], edx	

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdwc_loop
	jmp   .finish_vdwc_inner
.unroll_vdwc_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; dispersion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdwc_inner
	jmp   .unroll_vdwc_loop
.finish_vdwc_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdwc_inner
	jmp  .updateouterdata_vdwc		
.single_vdwc_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdwc:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdwc]
	jz  .testcoul
	jmp .mno_vdwc
.testcoul
	mov  ecx, [esp + .nscoul]
	cmp  ecx, byte 0
	jnz  .mno_coul
	jmp  .testvdw
.mno_coul:
	mov   ebx,  [esp + .solnr]
	inc   dword [esp + .solnr]
	mov   edx, [ebp + %$charge]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii]
	pfmul mm2, [ebp + %$facel]
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iq], mm2	        ;  iq =facel*charge[ii]

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	mov   [esp + .ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_coul_loop
	jmp   .finish_coul_inner
.unroll_coul_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$charge]        ; base of charge[]
	movq mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]	         ; charge[jnr1]
        punpckldq mm3, [ecx + ebx*4]     ; move charge 2 to high part of mm3 
	pfmul mm3,mm5		         ;  mm3 now has qq for both particles

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; coulomb table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + .tabscale]	
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_coul_inner
	jmp   .unroll_coul_loop
.finish_coul_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_coul_inner
	jmp  .updateouterdata_coul		
.single_coul_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov ecx, [ebp + %$charge]
	movd mm5, [esp + .iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;  mm3=qq

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2

	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, mm3	; vcoul=qq*VV
	pfmul mm3, mm7	; fijC=FF*qq

	; at this point mm5 contains vcoul and mm3 fijC.
	; increment vcoul - then we can get rid of mm5.
	;; update vctot
	pfadd mm5, [esp + .vctot]      ; add the earlier value
	movq [esp + .vctot], mm5       ;  store the sum       
	
	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_coul:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nscoul]
	jz  .testvdw
	jmp .mno_coul
.testvdw
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
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + .shX]
	pfadd mm1, [esp + .shZ]
	movq  [esp + .ix], mm0	
	movd  [esp + .iz], mm1	

	;; clear forces.
	pxor  mm7,mm7
	movq  [esp + .fix],   mm7
	movd  [esp + .fiz],   mm7

	mov   ecx, [esp + .innerjjnr0]
	mov   [esp + .innerjjnr], ecx
	mov   edx, [esp + .innerk0]
        sub   edx, dword 2
        mov   [esp + .innerk], edx        ;  number of innerloop atoms
	jge   .unroll_vdw_loop
	jmp   .finish_vdw_inner
.unroll_vdw_loop:	
	;; paired innerloop here.
	mov   ecx, [esp + .innerjjnr]     ; pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]             ; eax/ebx=jnr 
	add   [esp + .innerjjnr], dword 8 ; advance pointer (unrolled 2) 
	prefetch [ecx + 16]	         ; prefetch data - trial and error says 16 is best
	
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	mov ecx, [ecx + ebx*4]           ; type [jnr2]

	mov esi, [ebp + %$nbfp]		 ; base of nbfp 
	shl edx, 1
	shl ecx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	add ecx, [esp + .ntia]

	movq mm5, [esi + edx*4]		; mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		; mm7 = 2nd c6 / c12	
	movq mm6,mm5			
	punpckldq mm5,mm7		; mm5 = 1st c6 / 2nd c6  
	punpckhdq mm6,mm7		; mm6 = 1st c12 / 2nd c12
	movq [esp + .c6], mm5
	movq [esp + .c12], mm6

	lea   eax, [eax + eax*2]         ;  replace jnr with j3
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + %$pos]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]	 	
	movq  mm4, [esi + eax*4]         ;  fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		         ;  dr = ir - jr 
	pfsubr mm5,mm1
	movq  [esp + .dx1], mm4	         ;  store dr
	movd  [esp + .dz1], mm5
	pfmul mm4,mm4	                 ;  square dx,dy,dz		         
	pfmul mm5,mm5		
	pfacc mm4, mm5                   ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm4, mm5		         ;  first rsq in lower mm4

	movq  mm6, [esi + ebx*4]         ;  fetch second j coordinates 
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	                 ;  dr = ir - jr 
	pfsubr mm7,mm1
	movq  [esp + .dx2], mm6	         ;  store dr
	movd  [esp + .dz2], mm7
	pfmul mm6,mm6	                 ;  square dx,dy,dz
	pfmul mm7,mm7
	pfacc mm6, mm7		         ;  accumulate to get dx*dx+dy*dy+dz*dz
	pfacc mm6, mm7	                 ;  second rsq in lower mm6

        pfrsqrt mm0, mm4	         ; lookup inverse square root seed 
        pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;  now 4 has rsq and 0 the seed for both pairs.
        movq mm2,mm0	        	; amd 3dnow N-R iteration to get full precision.
	pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.
	;; do potential and fscal
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.
	
	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	; dispersion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	prefetchw [esp + .dx1]	;  prefetch i forces to cache

	;;  spread fscalar to both positions
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;; calc vector force
	prefetchw [edi + eax*4]	; prefetch the 1st faction to cache
	movq mm2,  [esp + .dx1]	; fetch dr
	movd mm3,  [esp + .dz1]

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	prefetchw [edi + ebx*4]	; prefetch the 2nd faction to cache
	pfmul mm2, mm0		; mult by fs 
	pfmul mm3, mm0

	movq mm4,  [esp + .dx2] 	; fetch dr
	movd mm5,  [esp + .dz2]
	pfmul mm4, mm1   	; mult by fs 
	pfmul mm5, mm1
	;; update i forces

	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j forces

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;; should we do one more iteration?
	sub   [esp + .innerk], dword 2
	jl    .finish_vdw_inner
	jmp   .unroll_vdw_loop
.finish_vdw_inner:	
	and [esp + .innerk], dword 1
	jnz  .single_vdw_inner
	jmp  .updateouterdata_vdw		
.single_vdw_inner:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset

	mov esi, [ebp + %$nbfp]
	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr1]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [esi + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [esi + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5

	mov   esi, [ebp + %$pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + .ix]
	movd  mm1, [esp + .iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + .dx1], mm4
	pfmul mm4,mm4
	movd  [esp + .dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;  mm0=rsq
	
        pfrsqrt mm0,mm4
        movq mm2,mm0
        pfmul mm0,mm0
        pfrsqit1 mm0,mm4				
        pfrcpit2 mm0,mm2	;  mm1=invsqrt
	pfmul mm4, mm0
	movq mm1, mm4
	;; mm0 is invsqrt, and mm1 r.

	;;  calculate potentials and scalar force
	pfmul mm1, [esp + .tabscale]	; mm1=rt
	pf2iw mm4,mm1
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4                   ; now mm1 is eps and mm4 n0.

	movq mm2,mm1
	pfmul mm2,mm2	; mm1 is eps, mm2 is eps2
	
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	movq mm3, mm7	; add to fscal

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm1  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of mm3
        pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + .tabscale]
 	pfmul mm0, mm1        ; mm0 is total fscal now	

	;; update vnbtot
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	;;  spread fscalar to both positions
	punpckldq mm0,mm0
	;;  calc vectorial force
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm2,  [esp + .dx1]
	movd mm3,  [esp + .dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;; update i particle force
	movq mm0,  [esp + .fix]
	movd mm1,  [esp + .fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + .fix], mm0
	movd [esp + .fiz], mm1
	;; update j particle force
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;;  done!
.updateouterdata_vdw:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment i force
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fix]
	pfadd mm7, [esp + .fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;;  loop back to mno.
	dec dword [esp + .nsvdw]
	jz  .last_mno
	jmp .mno_vdw
	
.last_mno:	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]
	;; finish if last
	mov   ecx, [ebp + %$nri]
	dec ecx
	jecxz .end
	;;  not last, iterate once more!
	mov [ebp + %$nri], ecx
	jmp .outer
.end:
	femms
	add esp, 168
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc


proc inl3320_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.iqO         equ    44		; repeated (64bit) to fill 3dnow reg
.iqH         equ    52		; repeated (64bit) to fill 3dnow reg
.qqO         equ    60		; repeated (64bit) to fill 3dnow reg
.qqH         equ    68		; repeated (64bit) to fill 3dnow reg
.vctot       equ    76          ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    84          ; repeated (64bit) to fill 3dnow reg
.c6          equ    92          ; repeated (64bit) to fill 3dnow reg
.c12         equ    100          ; repeated (64bit) to fill 3dnow reg
.two         equ    108         ; repeated (64bit) to fill 3dnow reg
.n1          equ    116         ; repeated (64bit) to fill 3dnow reg
.tabscale    equ    124         ; repeated (64bit) to fill 3dnow reg
.ntia        equ    132         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    140
.innerk      equ    144	
.fixO        equ    148
.fiyO        equ    152
.fizO        equ    156
.fixH        equ    160         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    168         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    176         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    184
.dyO	     equ    188
.dzO	     equ    192
.dxH	     equ    196         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    204         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    212         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    220        	; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 228		;   local stack space
	femms

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0]
	pfmul mm2, mm1		
	movq  [esp + .iqO], mm2	        ;  iqO = facel*charge[ii]
	
	movd  mm2, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1]
	pfmul mm2, mm1
	punpckldq mm2,mm2	        ;  spread to both halves
	movq  [esp + .iqH], mm2	        ;  iqH = facel*charge[ii0+1]

	mov   edx, [ebp + %$type] 	
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1		        
	imul  ecx, [ebp + %$ntype]      ; ecx = ntia = 2*ntype*type[ii0]
	mov   [esp + .ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [ebp + %$tabscale]
	punpckldq mm4,mm4	        ;  spread to both halves
	movq  [esp + .two],    mm3
	movq  [esp + .tabscale], mm4	      
	;; assume we have at least one i particle - start directly	 
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]
	
	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	
					
	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],   mm7
	movd  [esp + .fizO],   mm7
	movq  [esp + .fixH],   mm7
	movq  [esp + .fiyH],   mm7
	movq  [esp + .fizH],   mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:	
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 
	;prefetch [ecx + 16]	   ; prefetch data - trial and error says 16 is best

	mov ecx, [ebp + %$charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + .iqO]
	pfmul mm7, [esp + .iqH]	;  mm6=qqO, mm7=qqH
	movd [esp + .qqO], mm6
	movq [esp + .qqH], mm7

	mov ecx, [ebp + %$type]
	mov edx, [ecx + eax*4]        	 ; type [jnr]
	mov ecx, [ebp + %$nbfp]
	shl edx, 1
	add edx, [esp + .ntia]	         ; tja = ntia + 2*type
	movd mm5, [ecx + edx*4]		; mm5 = 1st c6 		
	movq [esp + .c6], mm5
	movd mm5, [ecx + edx*4 + 4]	; mm5 = 1st c12 		
	movq [esp + .c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3
	
        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt

	pfmul mm0, mm1		;  mm0=r

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqO]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	movq mm3, mm7
	
	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal  

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of fscal and multiply with rinv
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	
	
	;; update vnbtot 
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; now do the two hydrogens.
	movq mm0, [esp + .tmprsqH] ; mm0=rsqH

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	prefetchw [edi + eax*4]	; prefetch faction to cache
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2
	
	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential
	
	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  same for Vnb.
	
	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnb[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 228
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc

	

proc inl3330_3dnow
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
.is3         equ     0
.ii3         equ     4
.ixO         equ     8
.iyO         equ    12
.izO         equ    16	
.ixH         equ    20          ; repeated (64bit) to fill 3dnow reg
.iyH         equ    28          ; repeated (64bit) to fill 3dnow reg
.izH         equ    36          ; repeated (64bit) to fill 3dnow reg
.qqOO        equ    44		; repeated (64bit) to fill 3dnow reg
.qqOH        equ    52		; repeated (64bit) to fill 3dnow reg
.qqHH        equ    60     	; repeated (64bit) to fill 3dnow reg
.c6          equ    68     	; repeated (64bit) to fill 3dnow reg
.c12         equ    76     	; repeated (64bit) to fill 3dnow reg
.two         equ    84     	; repeated (64bit) to fill 3dnow reg
.n1	     equ    92     	; repeated (64bit) to fill 3dnow reg
.tabscale    equ    100    	; repeated (64bit) to fill 3dnow reg
.vctot       equ    108         ; repeated (64bit) to fill 3dnow reg
.vnbtot      equ    116         ; repeated (64bit) to fill 3dnow reg
.innerjjnr   equ    124
.innerk      equ    128	
.fixO        equ    132
.fiyO        equ    136
.fizO        equ    140
.fixH        equ    144         ; repeated (64bit) to fill 3dnow reg
.fiyH        equ    152         ; repeated (64bit) to fill 3dnow reg
.fizH        equ    160         ; repeated (64bit) to fill 3dnow reg
.dxO	     equ    168
.dyO	     equ    172
.dzO	     equ    176
.dxH	     equ    180         ; repeated (64bit) to fill 3dnow reg
.dyH	     equ    188         ; repeated (64bit) to fill 3dnow reg
.dzH	     equ    196         ; repeated (64bit) to fill 3dnow reg
.tmprsqH     equ    204         ; repeated (64bit) to fill 3dnow reg
        push eax
        push ebx
        push ecx
        push edx
	push esi
	push edi
	sub esp, 212		;   local stack space
	femms
	;; assume we have at least one i particle - start directly	

	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	mov   ebx, [ecx]	        ;  ebx =ii

	mov   edx, [ebp + %$charge]
	movd  mm1, [ebp + %$facel]	;  mm1=facel
	movd  mm2, [edx + ebx*4]        ;  mm2=charge[ii0] (O)
	movd  mm3, [edx + ebx*4 + 4]    ;  mm2=charge[ii0+1] (H)
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			; mm4=qqOO*facel
	pfmul mm5, mm3			; mm5=qqOH*facel
	pfmul mm6, mm3			; mm6=qqHH*facel
	punpckldq mm5,mm5	        ;  spread to both halves
	punpckldq mm6,mm6	        ;  spread to both halves
	movq  [esp + .qqOO], mm4
	movq  [esp + .qqOH], mm5
	movq  [esp + .qqHH], mm6
	mov   edx, [ebp + %$type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + %$ntype]
	add   edx, ecx
	mov   eax, [ebp + %$nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + .c6], mm0
	movq  [esp + .c12], mm1
	movq  mm2, [mm_two]
	movq  [esp + .two], mm2
	movd  mm3, [ebp + %$tabscale]
	punpckldq mm3,mm3
	movq  [esp + .tabscale], mm3
.outer:
	mov   eax, [ebp + %$shift]      ;  eax = pointer into shift[]
	mov   ebx, [eax]		;  ebx=shift[n]
	add   [ebp + %$shift], dword 4  ;  advance pointer one step
	
	lea   ebx, [ebx + ebx*2]        ;  ebx=3*is
	mov   [esp + .is3],ebx    	;  store is3

	mov   eax, [ebp + %$shiftvec]   ;  eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;  move shX/shY to mm5 and shZ to mm6.
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	        ; also expand shX,Y,Z in mm0--mm2.
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + %$iinr]       ;  ecx = pointer into iinr[]	
	add   [ebp + %$iinr], dword 4   ;  advance pointer
	mov   ebx, [ecx]	        ;  ebx =ii

	lea   ebx, [ebx + ebx*2]	;  ebx = 3*ii=ii3
	mov   eax, [ebp + %$pos]        ;  eax = base of pos[]

	pfadd mm5, [eax + ebx*4]        ; ix = shX + posX (and iy too)
	movd  mm7, [eax + ebx*4 + 8]    ; cant use direct memory add for 4 bytes (iz)
	mov   [esp + .ii3], ebx	        ; (use mm7 as temp. storage for iz.)
	pfadd mm6, mm7
	movq  [esp + .ixO], mm5	
	movq  [esp + .izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ; coords of H1 in low mm3-mm5, H2 in high.
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + .ixH], mm0	
	movq [esp + .iyH], mm1	
	movq [esp + .izH], mm2	

	;; clear vctot and i forces.
	pxor  mm7,mm7
	movq  [esp + .vctot], mm7
	movq  [esp + .vnbtot], mm7
	movq  [esp + .fixO],  mm7
	movq  [esp + .fizO],  mm7
	movq  [esp + .fixH],  mm7
	movq  [esp + .fiyH],  mm7
	movq  [esp + .fizH],  mm7

	mov   eax, [ebp + %$jindex]
	mov   ecx, [eax]	         ;  jindex[n]
	mov   edx, [eax + 4]	         ;  jindex[n+1]
	add   [ebp + %$jindex], dword 4
	sub   edx, ecx                   ;  number of innerloop atoms
	mov   [esp + .innerk], edx        ;  number of innerloop atoms

	mov   esi, [ebp + %$pos]
	mov   edi, [ebp + %$faction]	
	mov   eax, [ebp + %$jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + .innerjjnr], eax     ;  pointer to jjnr[nj0]
.inner_loop:
	;; a single j particle iteration here - compare with the unrolled code for comments.
	mov   eax, [esp + .innerjjnr]
	mov   eax, [eax]	;  eax=jnr offset
        add   [esp + .innerjjnr], dword 4 ; advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt OO
	pfmul mm0, mm1		;  mm0=rsq OO

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOO]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOO]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	movq mm3, mm7

	; dispersion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV	

	movq mm4, [esp + .c6]
	pfmul mm7, mm4	; fijD
	pfmul mm5, mm4	; vnb6           
	pfadd mm3, mm7	; add to fscal  

	;; update vnbtot to release mm5!
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       

	; repulsion table
	; load all the table values we need
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp
	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF
	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	movq mm6, [esp + .c12]
	pfmul mm7, mm6	; fijR
	pfmul mm5, mm6	; vnb12
	pfadd mm3, mm7	; total fscal fijC+fijD+fijR

	; change sign of fscal and multiply with rinv
        pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	
	
	;; update vnbtot 
	pfadd mm5, [esp + .vnbtot]      ; add the earlier value
	movq [esp + .vnbtot], mm5       ;  store the sum       
	
	;; Ready with the oxygen - potential is updated, fscal is in mm3.
	;; time for hydrogens!

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	
	
	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	; interactions with j H1.

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1		;  mm0=rsq 
	
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, force is moved to mm3.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now		

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4
	
	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	; interactions with j H2
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;;  copy & expand to mm2-mm4 for the H interactions
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + .ixO]
	pfsubr mm1, [esp + .izO]
		
	movq  [esp + .dxO], mm0
	pfmul mm0,mm0
	movd  [esp + .dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;  mm0=rsqO
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ; mm2-mm4 is jx-jz
	pfsubr mm2, [esp + .ixH]
	pfsubr mm3, [esp + .iyH]
	pfsubr mm4, [esp + .izH] ;  mm2-mm4 is dxH-dzH
	
	movq [esp + .dxH], mm2
	movq [esp + .dyH], mm3
	movq [esp + .dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;  mm3=rsqH
	movq [esp + .tmprsqH], mm3

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	pfmul mm0, mm1

	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movd [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2

	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqOH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqOH]	; fijC=qq*FF

	;; update vctot directly, use mm3 for fscal sum.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + .tabscale]
 	pfmul mm3, mm1        ; mm3 is total fscal (for the oxygen) now	

	movq mm0, [esp + .tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;  seeds are in mm1 now, and rsq in mm0.

	movq mm2, mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;  mm1=invsqrt
	
	pfmul mm0,mm1		; mm0=r
	pfmul mm0, [esp + .tabscale]
	pf2iw mm4, mm0
	movq [esp + .n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4                   ; now mm0 is eps and mm4 n0.
	movq  mm2, mm0
	pfmul mm2, mm2		;  mm0 is eps, mm2 eps2
	
	; coulomb table
	mov edx, [ebp + %$VFtab]
	mov ecx, [esp + .n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;;  load all values we need
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + .n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ; mm6 = Geps		
	pfmul mm7, mm2	; mm7 = Heps2
	;; 
	pfadd mm5, mm6
	pfadd mm5, mm7	; mm5 = Fp

	pfmul mm7, [esp + .two]	; two*Heps2
	pfadd mm7, mm6
	pfadd mm7, mm5	; mm7=FF

	pfmul mm5, mm0  ; mm5=eps*Fp
	pfadd mm5, mm4	; mm5= VV

	pfmul mm5, [esp + .qqHH]	; vcoul=qq*VV
	pfmul mm7, [esp + .qqHH]	; fijC=qq*FF
	;;  update vctot.
	pfadd mm5, [esp + .vctot]
	movq [esp + .vctot], mm5
	
	;; change sign of fijC and multiply by rinv
        pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + .tabscale]
 	pfmul mm4, mm1        ; mm4 is total fscal (for the hydrogens) now	

	;;  spread oxygen fscalar to both positions
	punpckldq mm3,mm3
	;;  calc vectorial force for O
	movq mm0,  [esp + .dxO]
	movd mm1,  [esp + .dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;;  calc vectorial force for H's
	movq mm5, [esp + .dxH]
	movq mm6, [esp + .dyH]
	movq mm7, [esp + .dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;; update iO particle force
	movq mm2,  [esp + .fixO]
	movd mm3,  [esp + .fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + .fixO], mm2
	movd [esp + .fizO], mm3

	;; update iH forces 
	movq mm2, [esp + .fixH]
	movq mm3, [esp + .fiyH]
	movq mm4, [esp + .fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + .fixH], mm2
	movq [esp + .fiyH], mm3
	movq [esp + .fizH], mm4	

	;; pack j forces from H in the same form as the oxygen force.
	pfacc mm5, mm6		; mm5(l)=fjx(H1+H2) mm5(h)=fjy(H1+H2)
	pfacc mm7, mm7		; mm7(l)=fjz(H1+H2) 
	
	pfadd mm0, mm5		;  add up total force on j particle.
	pfadd mm1, mm7

	;; update j particle force
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;;  done  - one more?
	dec dword [esp + .innerk]
	jz  .updateouterdata
	jmp .inner_loop	
.updateouterdata:	
	mov   ecx, [esp + .ii3]

	movq  mm6, [edi + ecx*4]       ;  increment iO force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + .fixH]
	movq  mm3, [esp + .fiyH]
	movq  mm1, [esp + .fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;  mm0(l)=fxH1, mm0(h)=fyH1
	punpckhdq mm2, mm3	;  mm2(l)=fxH2, mm2(h)=fyH2
	movq mm3, mm1
	pswapd mm3,mm3		
	;;  mm1 is fzH1
	;;  mm3 is fzH2

	movq  mm6, [edi + ecx*4 + 12]       ;  increment iH1 force
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;  increment iH2 force
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + %$fshift]    ; increment fshift force
	mov   edx, [esp + .is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + .fixO]
	pfadd mm7, [esp + .fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + %$gid]      ; get group index for this i particle
	mov   edx, [edx]
	add   [ebp + %$gid], dword 4  ;  advance pointer

	movq  mm7, [esp + .vctot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vc[gid]

	movq  mm7, [esp + .vnbtot]     
	pfacc mm7,mm7	              ;  get and sum the two parts of total potential

	mov   eax, [ebp + %$Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6              ; increment vnbtot[gid]
	;; finish if last
	dec dword [ebp + %$nri]
	jz  .end
	;;  not last, iterate once more!
	jmp .outer
.end:
	femms
	add esp, 212
	pop edi
	pop esi
        pop edx
        pop ecx
        pop ebx
        pop eax
	endproc
 
