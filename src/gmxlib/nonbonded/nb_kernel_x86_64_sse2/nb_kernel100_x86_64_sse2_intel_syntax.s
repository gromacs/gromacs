;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text




.globl nb_kernel100_x86_64_sse2
.globl _nb_kernel100_x86_64_sse2
nb_kernel100_x86_64_sse2:	
_nb_kernel100_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb100_fshift,           16
.equiv          nb100_gid,              24
.equiv          nb100_pos,              32
.equiv          nb100_faction,          40
.equiv          nb100_charge,           48
.equiv          nb100_p_facel,          56
.equiv          nb100_argkrf,           64
.equiv          nb100_argcrf,           72
.equiv          nb100_Vc,               80
.equiv          nb100_type,             88
.equiv          nb100_p_ntype,          96
.equiv          nb100_vdwparam,         104
.equiv          nb100_Vvdw,             112
.equiv          nb100_p_tabscale,       120
.equiv          nb100_VFtab,            128
.equiv          nb100_invsqrta,         136
.equiv          nb100_dvda,             144
.equiv          nb100_p_gbtabscale,     152
.equiv          nb100_GBtab,            160
.equiv          nb100_p_nthreads,       168
.equiv          nb100_count,            176
.equiv          nb100_mtx,              184
.equiv          nb100_outeriter,        192
.equiv          nb100_inneriter,        200
.equiv          nb100_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb100_ix,               0
.equiv          nb100_iy,               16
.equiv          nb100_iz,               32
.equiv          nb100_iq,               48
.equiv          nb100_dx,               64
.equiv          nb100_dy,               80
.equiv          nb100_dz,               96
.equiv          nb100_vctot,            112
.equiv          nb100_fix,              128
.equiv          nb100_fiy,              144
.equiv          nb100_fiz,              160
.equiv          nb100_half,             176
.equiv          nb100_three,            192
.equiv          nb100_is3,              208
.equiv          nb100_ii3,              212
.equiv          nb100_nri,              216
.equiv          nb100_iinr,             224
.equiv          nb100_jindex,           232
.equiv          nb100_jjnr,             240
.equiv          nb100_shift,            248
.equiv          nb100_shiftvec,         256
.equiv          nb100_facel,            264
.equiv          nb100_innerjjnr,        272
.equiv          nb100_innerk,           280
.equiv          nb100_n,                284
.equiv          nb100_nn1,              288
.equiv          nb100_nouter,           292
.equiv          nb100_ninner,           296

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 304		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb100_nouter], eax
	mov [rsp + nb100_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb100_nri], edi
	mov [rsp + nb100_iinr], rsi
	mov [rsp + nb100_jindex], rdx
	mov [rsp + nb100_jjnr], rcx
	mov [rsp + nb100_shift], r8
	mov [rsp + nb100_shiftvec], r9
	mov rsi, [rbp + nb100_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb100_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb100_half], eax
	movsd xmm1, [rsp + nb100_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movapd xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movapd [rsp + nb100_half],  xmm1
	movapd [rsp + nb100_three],  xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb100_half], eax
	mov [rsp + nb100_half+4], ebx
	movsd xmm1, [rsp + nb100_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb100_half], xmm1
	movapd [rsp + nb100_three], xmm3

.nb100_threadloop:
        mov   rsi, [rbp + nb100_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb100_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb100_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb100_n], eax
        mov [rsp + nb100_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100_outerstart
        jmp .nb100_end

.nb100_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb100_nouter]
	mov [rsp + nb100_nouter], ebx

.nb100_outer:
	mov   rax, [rsp + nb100_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb100_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb100_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb100_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb100_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb100_facel]
	shufpd xmm3, xmm3, 0	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb100_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb100_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb100_ix], xmm0
	movapd [rsp + nb100_iy], xmm1
	movapd [rsp + nb100_iz], xmm2

	mov   [rsp + nb100_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorpd xmm12, xmm12
	movapd xmm13, xmm12
	movapd xmm14, xmm12
	movapd xmm15, xmm12

	mov   rax, [rsp + nb100_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb100_pos]
	mov   rdi, [rbp + nb100_faction]	
	mov   rax, [rsp + nb100_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb100_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb100_ninner]
	mov   [rsp + nb100_ninner], ecx
	add   edx, 0
	mov   [rsp + nb100_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100_unroll_loop
	jmp   .nb100_checksingle
.nb100_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]
	add qword ptr [rsp + nb100_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb100_pos]       ;# base of pos[] 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		

	mov    rdi, [rbp + nb100_faction]
	
	;# calc dr 
	subpd xmm4, [rsp + nb100_ix]
	subpd xmm5, [rsp + nb100_iy]
	subpd xmm6, [rsp + nb100_iz]

	;# store dr 
	movapd xmm9,  xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6
    
   	mov rsi, [rbp + nb100_charge]    ;# base of charge[] 

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4  

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movlpd xmm3, [rsi + r8*8]	;# jq A 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb100_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movhpd xmm3, [rsi + r9*8]	;# jq B 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb100_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100_half]
	mulpd xmm3, [rsp + nb100_iq]		;# qq 
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 

	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	mulpd  xmm4, xmm3	;# xmm4=fscal 

	;# the fj's - start by combining forces from memory 
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]

    ;# increment vctot
	addpd  xmm12, xmm3

 	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

	addpd xmm0, xmm9
	addpd xmm1, xmm10
	addpd xmm2, xmm11

	;# xmm9-xmm11 contains tx-tz (partial force) 
	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11

	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb100_innerk],  2
	jl    .nb100_checksingle
	jmp   .nb100_unroll_loop

.nb100_checksingle:				
	mov   edx, [rsp + nb100_innerk]
	and   edx, 1
	jnz    .nb100_dosingle
	jmp    .nb100_updateouterdata
.nb100_dosingle:			
	mov rsi, [rbp + nb100_charge]
	mov rdi, [rbp + nb100_pos]

	mov rdx, [rsp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov eax, [rdx]	

	mov rsi, [rbp + nb100_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]	;# jq A 
	mulsd xmm3, [rsp + nb100_iq]		;# qq 
	
	mov rsi, [rbp + nb100_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm4-xmm6
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]

	mov    rdi, [rbp + nb100_faction]
	
	;# calc dr 
	subsd xmm4, [rsp + nb100_ix]
	subsd xmm5, [rsp + nb100_iy]
	subsd xmm6, [rsp + nb100_iz]

	;# store dr 
	movapd xmm9,  xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4  

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb100_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb100_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 

	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	mulsd  xmm4, xmm3	;# xmm4=fscal 

    ;# increment vctot
	addsd  xmm12, xmm3

 	mulsd  xmm9, xmm4
	mulsd  xmm10, xmm4
	mulsd  xmm11, xmm4
	;# xmm9-xmm11 contains tx-tz (partial force) 
	;# now update f_i 
	addsd  xmm13, xmm9
	addsd  xmm14, xmm10
	addsd  xmm15, xmm11

	;# the fj's - start by accumulating forces from memory 
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11

.nb100_updateouterdata:
	mov   ecx, [rsp + nb100_ii3]
	mov   rdi, [rbp + nb100_faction]
	mov   rsi, [rbp + nb100_fshift]
	mov   edx, [rsp + nb100_is3]

	;# accumulate i forces in xmm13, xmm14, xmm15
	movhlps xmm3, xmm13
	movhlps xmm4, xmm14
	movhlps xmm5, xmm15
	addsd  xmm13, xmm3
	addsd  xmm14, xmm4
	addsd  xmm15, xmm5 ;# sum is in low xmm13-xmm15

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [rsi + rdx*8]
	movsd  xmm4, [rsi + rdx*8 + 8]
	movsd  xmm5, [rsi + rdx*8 + 16]
	subsd  xmm3, xmm13
	subsd  xmm4, xmm14
	subsd  xmm5, xmm15
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb100_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb100_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total coulomb energy and update it 
	movhlps xmm6, xmm12
	addsd  xmm12, xmm6	;# low xmm12 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb100_Vc]
	addsd xmm12, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb100_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb100_n], esi
        jmp .nb100_outer
.nb100_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb100_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100_end
        ;# non-zero, do one more workunit
        jmp   .nb100_threadloop
.nb100_end:
	mov eax, [rsp + nb100_nouter]
	mov ebx, [rsp + nb100_ninner]
	mov rcx, [rbp + nb100_outeriter]
	mov rdx, [rbp + nb100_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 304
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret




.globl nb_kernel100nf_x86_64_sse2
.globl _nb_kernel100nf_x86_64_sse2
nb_kernel100nf_x86_64_sse2:	
_nb_kernel100nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb100nf_fshift,         16
.equiv          nb100nf_gid,            24
.equiv          nb100nf_pos,            32
.equiv          nb100nf_faction,        40
.equiv          nb100nf_charge,         48
.equiv          nb100nf_p_facel,        56
.equiv          nb100nf_argkrf,         64
.equiv          nb100nf_argcrf,         72
.equiv          nb100nf_Vc,             80
.equiv          nb100nf_type,           88
.equiv          nb100nf_p_ntype,        96
.equiv          nb100nf_vdwparam,       104
.equiv          nb100nf_Vvdw,           112
.equiv          nb100nf_p_tabscale,     120
.equiv          nb100nf_VFtab,          128
.equiv          nb100nf_invsqrta,       136
.equiv          nb100nf_dvda,           144
.equiv          nb100nf_p_gbtabscale,   152
.equiv          nb100nf_GBtab,          160
.equiv          nb100nf_p_nthreads,     168
.equiv          nb100nf_count,          176
.equiv          nb100nf_mtx,            184
.equiv          nb100nf_outeriter,      192
.equiv          nb100nf_inneriter,      200
.equiv          nb100nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100nf_ix,             0
.equiv          nb100nf_iy,             16
.equiv          nb100nf_iz,             32
.equiv          nb100nf_iq,             48
.equiv          nb100nf_vctot,          64
.equiv          nb100nf_half,           80
.equiv          nb100nf_three,          96
.equiv          nb100nf_is3,            112
.equiv          nb100nf_ii3,            116
.equiv          nb100nf_nri,            120
.equiv          nb100nf_iinr,           128
.equiv          nb100nf_jindex,         136
.equiv          nb100nf_jjnr,           144
.equiv          nb100nf_shift,          156
.equiv          nb100nf_shiftvec,       164
.equiv          nb100nf_facel,          172
.equiv          nb100nf_innerjjnr,      180
.equiv          nb100nf_innerk,         188
.equiv          nb100nf_n,              192
.equiv          nb100nf_nn1,            196
.equiv          nb100nf_nouter,         200
.equiv          nb100nf_ninner,         204

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	emms
	sub rsp, 208		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb100nf_nouter], eax
	mov [rsp + nb100nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb100nf_nri], edi
	mov [rsp + nb100nf_iinr], rsi
	mov [rsp + nb100nf_jindex], rdx
	mov [rsp + nb100nf_jjnr], rcx
	mov [rsp + nb100nf_shift], r8
	mov [rsp + nb100nf_shiftvec], r9
	mov rsi, [rbp + nb100nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb100nf_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb100nf_half], eax
	mov [rsp + nb100nf_half+4], ebx
	movsd xmm1, [rsp + nb100nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb100nf_half], xmm1
	movapd [rsp + nb100nf_three], xmm3

.nb100nf_threadloop:
        mov   rsi, [rbp + nb100nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb100nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb100nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb100nf_n], eax
        mov [rsp + nb100nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100nf_outerstart
        jmp .nb100nf_end

.nb100nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb100nf_nouter]
	mov [rsp + nb100nf_nouter], ebx

.nb100nf_outer:
	mov   rax, [rsp + nb100nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb100nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb100nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb100nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb100nf_facel]
	shufpd xmm3, xmm3, 0	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb100nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb100nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb100nf_ix], xmm0
	movapd [rsp + nb100nf_iy], xmm1
	movapd [rsp + nb100nf_iz], xmm2

	mov   [rsp + nb100nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb100nf_vctot], xmm4
	
	mov   rax, [rsp + nb100nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb100nf_pos]
	mov   rax, [rsp + nb100nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb100nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb100nf_ninner]
	mov   [rsp + nb100nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb100nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100nf_unroll_loop
	jmp   .nb100nf_checksingle
.nb100nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb100nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb100nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]	;# jq A 
	movhpd xmm3, [rsi + rbx*8]	;# jq B 

	movapd xmm5, [rsp + nb100nf_iq]
	
	mulpd xmm3, xmm5		;# qq 
	
	mov rsi, [rbp + nb100nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb100nf_ix]
	movapd xmm5, [rsp + nb100nf_iy]
	movapd xmm6, [rsp + nb100nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4  

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb100nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb100nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	movapd xmm5, [rsp + nb100nf_vctot]
	mulpd  xmm3, xmm0	;# xmm3=vcoul 
	addpd  xmm5, xmm3
	movapd [rsp + nb100nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb100nf_innerk],  2
	jl    .nb100nf_checksingle
	jmp   .nb100nf_unroll_loop

.nb100nf_checksingle:				
	mov   edx, [rsp + nb100nf_innerk]
	and   edx, 1
	jnz    .nb100nf_dosingle
	jmp    .nb100nf_updateouterdata
.nb100nf_dosingle:			
	mov rsi, [rbp + nb100nf_charge]
	mov rdi, [rbp + nb100nf_pos]

	mov rdx, [rsp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov eax, [rdx]	

	xorpd xmm3, xmm3
	movsd xmm3, [rsi + rax*8]	;# jq A 
	movapd xmm5, [rsp + nb100nf_iq]
	unpcklpd xmm3, xmm6
	mulpd xmm3, xmm5		;# qq 
	
	mov rsi, [rbp + nb100nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb100nf_ix]
	movapd xmm5, [rsp + nb100nf_iy]
	movapd xmm6, [rsp + nb100nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb100nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb100nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb100nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 

	movlpd xmm5, [rsp + nb100nf_vctot]
	mulsd  xmm3, xmm0	;# xmm3=vcoul 
	addsd  xmm5, xmm3
	movlpd [rsp + nb100nf_vctot], xmm5

.nb100nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb100nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb100nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb100nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb100nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb100nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb100nf_n], esi
        jmp .nb100nf_outer
.nb100nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb100nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100nf_end
        ;# non-zero, do one more workunit
        jmp   .nb100nf_threadloop
.nb100nf_end:
	mov eax, [rsp + nb100nf_nouter]
	mov ebx, [rsp + nb100nf_ninner]
	mov rcx, [rbp + nb100nf_outeriter]
	mov rdx, [rbp + nb100nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 208
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret


