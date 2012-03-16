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


.globl nb_kernel200_x86_64_sse2
.globl _nb_kernel200_x86_64_sse2
nb_kernel200_x86_64_sse2:	
_nb_kernel200_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb200_fshift,           16
.equiv          nb200_gid,              24
.equiv          nb200_pos,              32
.equiv          nb200_faction,          40
.equiv          nb200_charge,           48
.equiv          nb200_p_facel,          56
.equiv          nb200_argkrf,           64
.equiv          nb200_argcrf,           72
.equiv          nb200_Vc,               80
.equiv          nb200_type,             88
.equiv          nb200_p_ntype,          96
.equiv          nb200_vdwparam,         104
.equiv          nb200_Vvdw,             112
.equiv          nb200_p_tabscale,       120
.equiv          nb200_VFtab,            128
.equiv          nb200_invsqrta,         136
.equiv          nb200_dvda,             144
.equiv          nb200_p_gbtabscale,     152
.equiv          nb200_GBtab,            160
.equiv          nb200_p_nthreads,       168
.equiv          nb200_count,            176
.equiv          nb200_mtx,              184
.equiv          nb200_outeriter,        192
.equiv          nb200_inneriter,        200
.equiv          nb208_work,             200
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb200_ix,               0
.equiv          nb200_iy,               16
.equiv          nb200_iz,               32
.equiv          nb200_iq,               48
.equiv          nb200_dx,               64
.equiv          nb200_dy,               80
.equiv          nb200_dz,               96
.equiv          nb200_vctot,            112
.equiv          nb200_fix,              128
.equiv          nb200_fiy,              144
.equiv          nb200_fiz,              160
.equiv          nb200_half,             176
.equiv          nb200_three,            192
.equiv          nb200_two,              208
.equiv          nb200_krf,              224
.equiv          nb200_crf,              240
.equiv          nb200_is3,              256
.equiv          nb200_ii3,              260
.equiv          nb200_nri,              264
.equiv          nb200_iinr,             272
.equiv          nb200_jindex,           280
.equiv          nb200_jjnr,             288
.equiv          nb200_shift,            296
.equiv          nb200_shiftvec,         304
.equiv          nb200_facel,            312
.equiv          nb200_innerjjnr,        320
.equiv          nb200_innerk,           328
.equiv          nb200_n,                332
.equiv          nb200_nn1,              336
.equiv          nb200_nouter,           340
.equiv          nb200_ninner,           344

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
	sub rsp, 352		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb200_nouter], eax
	mov [rsp + nb200_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb200_nri], edi
	mov [rsp + nb200_iinr], rsi
	mov [rsp + nb200_jindex], rdx
	mov [rsp + nb200_jjnr], rcx
	mov [rsp + nb200_shift], r8
	mov [rsp + nb200_shiftvec], r9
	mov rsi, [rbp + nb200_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb200_facel], xmm0

	mov rsi, [rbp + nb200_argkrf]
	mov rdi, [rbp + nb200_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb200_krf], xmm1
	movapd [rsp + nb200_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb200_half], eax
	mov [rsp + nb200_half+4], ebx
	movsd xmm1, [rsp + nb200_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb200_half], xmm1
	movapd [rsp + nb200_two], xmm2
	movapd [rsp + nb200_three], xmm3

.nb200_threadloop:
        mov   rsi, [rbp + nb200_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb200_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb200_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb200_n], eax
        mov [rsp + nb200_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200_outerstart
        jmp .nb200_end

.nb200_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb200_nouter]
	mov [rsp + nb200_nouter], ebx

.nb200_outer:
	mov   rax, [rsp + nb200_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb200_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb200_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb200_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb200_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb200_facel]
	shufpd xmm3, xmm3, 0
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb200_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb200_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb200_ix], xmm0
	movapd [rsp + nb200_iy], xmm1
	movapd [rsp + nb200_iz], xmm2

	mov   [rsp + nb200_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorpd xmm12, xmm12
	movapd xmm13, xmm12
	movapd xmm14, xmm12
	movapd xmm15, xmm12
		
	mov   rax, [rsp + nb200_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb200_pos]
	mov   rdi, [rbp + nb200_faction]	
	mov   rax, [rsp + nb200_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb200_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb200_ninner]
	mov   [rsp + nb200_ninner], ecx
	add   edx, 0
	mov   [rsp + nb200_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200_unroll_loop
	jmp   .nb200_checksingle
.nb200_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb200_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	add qword ptr [rsp + nb200_innerjjnr],  8	;# advance pointer (unrolled 2) 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	

	mov rsi, [rbp + nb200_pos]       ;# base of pos[] 

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		

	;# calc dr 
	subpd xmm4, [rsp + nb200_ix]
	subpd xmm5, [rsp + nb200_iy]
	subpd xmm6, [rsp + nb200_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	mov rsi, [rbp + nb200_charge]    ;# base of charge[] 
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	
	movlpd xmm3, [rsi + r8*8]

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movhpd xmm3, [rsi + r9*8]
	movapd xmm7, [rsp + nb200_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb200_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	mulpd xmm3, [rsp + nb200_iq]		;# qq 
	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb200_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [rsp + nb200_crf]
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	
	mov    rdi, [rbp + nb200_faction]

	addpd  xmm7, xmm7

	subpd  xmm0, xmm7
	mulpd  xmm3, xmm0	
	mulpd  xmm4, xmm3	;# xmm4=total fscal 

    ;# increment vctot
	addpd  xmm12, xmm6

	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4

	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]

	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11

	addpd xmm3, xmm9
	addpd xmm4, xmm10
	addpd xmm5, xmm11
    
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5	

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb200_innerk],  2
	jl    .nb200_checksingle
	jmp   .nb200_unroll_loop

.nb200_checksingle:				
	mov   edx, [rsp + nb200_innerk]
	and   edx, 1
	jnz    .nb200_dosingle
	jmp    .nb200_updateouterdata
.nb200_dosingle:			
	mov rsi, [rbp + nb200_charge]
	mov rdi, [rbp + nb200_pos]
	mov   rcx, [rsp + nb200_innerjjnr]
	mov   eax, [rcx]
    
	mov rsi, [rbp + nb200_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]
	mulsd xmm3, [rsp + nb200_iq]		;# qq 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	mov rsi, [rbp + nb200_pos]       ;# base of pos[] 

	;# move two coordinates to xmm4-xmm6
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]

	;# calc dr 
	subsd xmm4, [rsp + nb200_ix]
	subsd xmm5, [rsp + nb200_iy]
	subsd xmm6, [rsp + nb200_iz]

	;# store dr 
	movapd xmm9, xmm4
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

	movapd xmm7, [rsp + nb200_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb200_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb200_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [rsp + nb200_crf]
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	
	addsd  xmm7, xmm7

	subsd  xmm0, xmm7
	mulsd  xmm3, xmm0	
	mulsd  xmm4, xmm3	;# xmm4=total fscal 

    ;# increment vctot
	addsd  xmm12, xmm6

	mov    rdi, [rbp + nb200_faction]
	mulsd  xmm9, xmm4
	mulsd  xmm10, xmm4
	mulsd  xmm11, xmm4
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

.nb200_updateouterdata:
	mov   ecx, [rsp + nb200_ii3]
	mov   rdi, [rbp + nb200_faction]
	mov   rsi, [rbp + nb200_fshift]
	mov   edx, [rsp + nb200_is3]

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
	mov esi, [rsp + nb200_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb200_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total coulomb energy and update it 
	movhlps xmm6, xmm12
	addsd  xmm12, xmm6	;# low xmm12 have the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb200_Vc]
	addsd xmm12, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb200_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb200_n], esi
        jmp .nb200_outer
.nb200_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb200_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200_end
        ;# non-zero, do one more workunit
        jmp   .nb200_threadloop
.nb200_end:
	mov eax, [rsp + nb200_nouter]
	mov ebx, [rsp + nb200_ninner]
	mov rcx, [rbp + nb200_outeriter]
	mov rdx, [rbp + nb200_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 352
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





.globl nb_kernel200nf_x86_64_sse2
.globl _nb_kernel200nf_x86_64_sse2
nb_kernel200nf_x86_64_sse2:	
_nb_kernel200nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb200nf_fshift,         16
.equiv          nb200nf_gid,            24
.equiv          nb200nf_pos,            32
.equiv          nb200nf_faction,        40
.equiv          nb200nf_charge,         48
.equiv          nb200nf_p_facel,        56
.equiv          nb200nf_argkrf,         64
.equiv          nb200nf_argcrf,         72
.equiv          nb200nf_Vc,             80
.equiv          nb200nf_type,           88
.equiv          nb200nf_p_ntype,        96
.equiv          nb200nf_vdwparam,       104
.equiv          nb200nf_Vvdw,           112
.equiv          nb200nf_p_tabscale,     120
.equiv          nb200nf_VFtab,          128
.equiv          nb200nf_invsqrta,       136
.equiv          nb200nf_dvda,           144
.equiv          nb200nf_p_gbtabscale,   152
.equiv          nb200nf_GBtab,          160
.equiv          nb200nf_p_nthreads,     168
.equiv          nb200nf_count,          176
.equiv          nb200nf_mtx,            184
.equiv          nb200nf_outeriter,      192
.equiv          nb200nf_inneriter,      200
.equiv          nb208nf_work,           200
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200nf_ix,             0
.equiv          nb200nf_iy,             16
.equiv          nb200nf_iz,             32
.equiv          nb200nf_iq,             48
.equiv          nb200nf_vctot,          64
.equiv          nb200nf_half,           80
.equiv          nb200nf_three,          96
.equiv          nb200nf_krf,            112
.equiv          nb200nf_crf,            128
.equiv          nb200nf_is3,            144
.equiv          nb200nf_ii3,            148
.equiv          nb200nf_nri,            152
.equiv          nb200nf_iinr,           160
.equiv          nb200nf_jindex,         168
.equiv          nb200nf_jjnr,           176
.equiv          nb200nf_shift,          184
.equiv          nb200nf_shiftvec,       192
.equiv          nb200nf_facel,          200
.equiv          nb200nf_innerjjnr,      208
.equiv          nb200nf_innerk,         216
.equiv          nb200nf_n,              220
.equiv          nb200nf_nn1,            224
.equiv          nb200nf_nouter,         228
.equiv          nb200nf_ninner,         232

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
	sub rsp, 240		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb200nf_nouter], eax
	mov [rsp + nb200nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb200nf_nri], edi
	mov [rsp + nb200nf_iinr], rsi
	mov [rsp + nb200nf_jindex], rdx
	mov [rsp + nb200nf_jjnr], rcx
	mov [rsp + nb200nf_shift], r8
	mov [rsp + nb200nf_shiftvec], r9
	mov rsi, [rbp + nb200nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb200nf_facel], xmm0

	mov rsi, [rbp + nb200nf_argkrf]
	mov rdi, [rbp + nb200nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb200nf_krf], xmm1
	movapd [rsp + nb200nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb200nf_half], eax
	mov [rsp + nb200nf_half+4], ebx
	movsd xmm1, [rsp + nb200nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb200nf_half], xmm1
	movapd [rsp + nb200nf_three], xmm3

.nb200nf_threadloop:
        mov   rsi, [rbp + nb200nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb200nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb200nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb200nf_n], eax
        mov [rsp + nb200nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200nf_outerstart
        jmp .nb200nf_end

.nb200nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb200nf_nouter]
	mov [rsp + nb200nf_nouter], ebx

.nb200nf_outer:
	mov   rax, [rsp + nb200nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb200nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb200nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb200nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb200nf_facel]
	shufpd xmm3, xmm3, 0
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb200nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb200nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb200nf_ix], xmm0
	movapd [rsp + nb200nf_iy], xmm1
	movapd [rsp + nb200nf_iz], xmm2

	mov   [rsp + nb200nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb200nf_vctot], xmm4
	
	mov   rax, [rsp + nb200nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb200nf_pos]
	mov   rax, [rsp + nb200nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb200nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb200nf_ninner]
	mov   [rsp + nb200nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb200nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200nf_unroll_loop
	jmp   .nb200nf_checksingle
.nb200nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb200nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb200nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb200nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm5, [rsp + nb200nf_iq]
	mulpd xmm3, xmm5		;# qq 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	mov rsi, [rbp + nb200nf_pos]       ;# base of pos[] 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb200nf_ix]
	movapd xmm5, [rsp + nb200nf_iy]
	movapd xmm6, [rsp + nb200nf_iz]

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

	movapd xmm7, [rsp + nb200nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb200nf_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb200nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [rsp + nb200nf_crf]
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addpd  xmm6, [rsp + nb200nf_vctot]
	movapd [rsp + nb200nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb200nf_innerk],  2
	jl    .nb200nf_checksingle
	jmp   .nb200nf_unroll_loop

.nb200nf_checksingle:				
	mov   edx, [rsp + nb200nf_innerk]
	and   edx, 1
	jnz    .nb200nf_dosingle
	jmp    .nb200nf_updateouterdata
.nb200nf_dosingle:			
	mov rsi, [rbp + nb200nf_charge]
	mov rdi, [rbp + nb200nf_pos]
	mov   rcx, [rsp + nb200nf_innerjjnr]
	
	xorpd xmm3, xmm3
	mov   eax, [rcx]

	movlpd xmm3, [rsi + rax*8]
	movapd xmm5, [rsp + nb200nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	mov rsi, [rbp + nb200nf_pos]       ;# base of pos[] 

	lea rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb200nf_ix]
	movapd xmm5, [rsp + nb200nf_iy]
	movapd xmm6, [rsp + nb200nf_iz]

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

	movapd xmm7, [rsp + nb200nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb200nf_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb200nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb200nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [rsp + nb200nf_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addsd  xmm6, [rsp + nb200nf_vctot]
	movlpd [rsp + nb200nf_vctot], xmm6

.nb200nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb200nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb200nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb200nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb200nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb200nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb200nf_n], esi
        jmp .nb200nf_outer
.nb200nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb200nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200nf_end
        ;# non-zero, do one more workunit
        jmp   .nb200nf_threadloop
.nb200nf_end:
	mov eax, [rsp + nb200nf_nouter]
	mov ebx, [rsp + nb200nf_ninner]
	mov rcx, [rbp + nb200nf_outeriter]
	mov rdx, [rbp + nb200nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 240
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



