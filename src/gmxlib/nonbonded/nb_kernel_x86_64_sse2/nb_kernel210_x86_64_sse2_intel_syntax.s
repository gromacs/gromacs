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


.globl nb_kernel210_x86_64_sse2
.globl _nb_kernel210_x86_64_sse2
nb_kernel210_x86_64_sse2:	
_nb_kernel210_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb210_fshift,           16
.equiv          nb210_gid,              24
.equiv          nb210_pos,              32
.equiv          nb210_faction,          40
.equiv          nb210_charge,           48
.equiv          nb210_p_facel,          56
.equiv          nb210_argkrf,           64
.equiv          nb210_argcrf,           72  
.equiv          nb210_Vc,               80
.equiv          nb210_type,             88
.equiv          nb210_p_ntype,          96
.equiv          nb210_vdwparam,         104
.equiv          nb210_Vvdw,             112
.equiv          nb210_p_tabscale,       120
.equiv          nb210_VFtab,            128
.equiv          nb210_invsqrta,         136
.equiv          nb210_dvda,             144
.equiv          nb210_p_gbtabscale,     152
.equiv          nb210_GBtab,            160
.equiv          nb210_p_nthreads,       168
.equiv          nb210_count,            176
.equiv          nb210_mtx,              184
.equiv          nb210_outeriter,        192
.equiv          nb210_inneriter,        200
.equiv          nb210_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb210_ix,               0
.equiv          nb210_iy,               16
.equiv          nb210_iz,               32
.equiv          nb210_iq,               48
.equiv          nb210_dx,               64
.equiv          nb210_dy,               80
.equiv          nb210_dz,               96
.equiv          nb210_c6,               112
.equiv          nb210_c12,              128
.equiv          nb210_six,              144
.equiv          nb210_twelve,           160
.equiv          nb210_vctot,            176
.equiv          nb210_Vvdwtot,          192
.equiv          nb210_fix,              208
.equiv          nb210_fiy,              224
.equiv          nb210_fiz,              240
.equiv          nb210_half,             256
.equiv          nb210_three,            272
.equiv          nb210_two,              288
.equiv          nb210_krf,              304
.equiv          nb210_crf,              320
.equiv          nb210_nri,              336
.equiv          nb210_iinr,             344
.equiv          nb210_jindex,           352
.equiv          nb210_jjnr,             360
.equiv          nb210_shift,            368
.equiv          nb210_shiftvec,         376
.equiv          nb210_facel,            384
.equiv          nb210_innerjjnr,        392
.equiv          nb210_is3,              400
.equiv          nb210_ii3,              404
.equiv          nb210_ntia,             408
.equiv          nb210_innerk,           412
.equiv          nb210_n,                416
.equiv          nb210_nn1,              420
.equiv          nb210_ntype,            424
.equiv          nb210_nouter,           428
.equiv          nb210_ninner,           432

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
	sub rsp, 448		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb210_nouter], eax
	mov [rsp + nb210_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb210_nri], edi
	mov [rsp + nb210_iinr], rsi
	mov [rsp + nb210_jindex], rdx
	mov [rsp + nb210_jjnr], rcx
	mov [rsp + nb210_shift], r8
	mov [rsp + nb210_shiftvec], r9
	mov rdi, [rbp + nb210_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb210_ntype], edi
	mov rsi, [rbp + nb210_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb210_facel], xmm0

	mov rsi, [rbp + nb210_argkrf]
	mov rdi, [rbp + nb210_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb210_krf], xmm1
	movapd [rsp + nb210_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb210_half], eax
	mov [rsp + nb210_half+4], ebx
	movsd xmm1, [rsp + nb210_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb210_half], xmm1
	movapd [rsp + nb210_two], xmm2
	movapd [rsp + nb210_three], xmm3
	movapd [rsp + nb210_six], xmm4
	movapd [rsp + nb210_twelve], xmm5

.nb210_threadloop:
        mov   rsi, [rbp + nb210_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb210_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb210_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb210_n], eax
        mov [rsp + nb210_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210_outerstart
        jmp .nb210_end

.nb210_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb210_nouter]
	mov [rsp + nb210_nouter], ebx

.nb210_outer:
	mov   rax, [rsp + nb210_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb210_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb210_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb210_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb210_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb210_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb210_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb210_ntype]
    	shl   edx, 1
    	mov   [rsp + nb210_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb210_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb210_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb210_ix], xmm0
	movapd [rsp + nb210_iy], xmm1
	movapd [rsp + nb210_iz], xmm2

	mov   [rsp + nb210_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm12, xmm12
	movapd xmm13, xmm12
	movapd xmm14, xmm12
	movapd xmm15, xmm12
	movapd [rsp + nb210_vctot], xmm12
	movapd [rsp + nb210_Vvdwtot], xmm12
	
	mov   rax, [rsp + nb210_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb210_pos]
	mov   rdi, [rbp + nb210_faction]	
	mov   rax, [rsp + nb210_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb210_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb210_ninner]
	mov   [rsp + nb210_ninner], ecx
	add   edx, 0
	mov   [rsp + nb210_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210_unroll_loop
	jmp   .nb210_checksingle
.nb210_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb210_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r10d, [rdx]	
	mov   r11d, [rdx + 4]              
	add qword ptr [rsp + nb210_innerjjnr],  8	;# advance pointer (unrolled 2) 

	
	mov rsi, [rbp + nb210_pos]       ;# base of pos[] 

	lea   rax, [r10 + r10*2]     ;# replace jnr with j3 
	lea   rbx, [r11 + r11*2]	

	;# move two coordinates to xmm4-xmm6
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb210_ix]
	subpd xmm5, [rsp + nb210_iy]
	subpd xmm6, [rsp + nb210_iz]

	;# store dr 
	movapd xmm9, xmm4
	movapd xmm10, xmm5
	movapd xmm11, xmm6

	mov rsi, [rbp + nb210_charge]    ;# base of charge[] 
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 
	mov rdi, [rbp + nb210_type]

	movlpd xmm3, [rsi + r10*8]
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	movhpd xmm3, [rsi + r11*8]
	mov r8d, [rdi + r10*4]
	mov r9d, [rdi + r11*4]

	movapd xmm7, [rsp + nb210_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb210_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	shl r8d, 1
	shl r9d, 1
	mov edi, [rsp + nb210_ntia]
	movapd xmm0, [rsp + nb210_half]
	mov rsi, [rbp + nb210_vdwparam]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 
	mulpd xmm3, [rsp + nb210_iq]		;# qq 
	add r8d, edi
	add r9d, edi

    
	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb210_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210_half]
    
    movlpd  xmm4, [rsi + r8*8]
    movlpd  xmm6, [rsi + r8*8 + 8]
    movhpd  xmm4, [rsi + r9*8]
    movhpd  xmm6, [rsi + r9*8 + 8]
    movapd [rsp + nb210_c6], xmm4
    movapd [rsp + nb210_c12], xmm6
    
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [rsp + nb210_crf]
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addpd  xmm7, xmm7
	mulpd  xmm1, [rsp + nb210_c6]
	mulpd  xmm2, [rsp + nb210_c12]

	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb210_Vvdwtot]
	mulpd  xmm1, [rsp + nb210_six]
	mulpd  xmm2, [rsp + nb210_twelve]
	subpd  xmm2, xmm1
	subpd  xmm0, xmm7
	mulpd  xmm3, xmm0
	addpd  xmm2, xmm3
	mulpd  xmm4, xmm2	;# xmm4=total fscal 
	addpd  xmm6, [rsp + nb210_vctot]

	;# the fj's - start by accumulating forces from memory 
	mov    rdi, [rbp + nb210_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

	movapd [rsp + nb210_vctot], xmm6
	movapd [rsp + nb210_Vvdwtot], xmm5

	mulpd  xmm9, xmm4
	mulpd  xmm10, xmm4
	mulpd  xmm11, xmm4

	addpd xmm0, xmm9
	addpd xmm1, xmm10
	addpd xmm2, xmm11
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	;# now update f_i 
	addpd  xmm13, xmm9
	addpd  xmm14, xmm10
	addpd  xmm15, xmm11
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb210_innerk],  2
	jl    .nb210_checksingle
	jmp   .nb210_unroll_loop

.nb210_checksingle:				
	mov   edx, [rsp + nb210_innerk]
	and   edx, 1
	jnz    .nb210_dosingle
	jmp    .nb210_updateouterdata
.nb210_dosingle:			
	mov rsi, [rbp + nb210_charge]
	mov rdi, [rbp + nb210_pos]
	mov   rcx, [rsp + nb210_innerjjnr]

	mov   eax, [rcx]

	mov rsi, [rbp + nb210_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]
	mulsd xmm3, [rsp + nb210_iq]		;# qq 

	mov rsi, [rbp + nb210_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb210_vdwparam]
	shl r8d, 1
	mov edi, [rsp + nb210_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]	
	movsd xmm6, [rsi + r8*8 + 8]
	movapd [rsp + nb210_c6], xmm4
	movapd [rsp + nb210_c12], xmm6
	
	mov rsi, [rbp + nb210_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm4-xmm6
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]
	
	;# calc dr 
	subsd xmm4, [rsp + nb210_ix]
	subsd xmm5, [rsp + nb210_iy]
	subsd xmm6, [rsp + nb210_iz]

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

	movapd xmm7, [rsp + nb210_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb210_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb210_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [rsp + nb210_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addsd  xmm7, xmm7
	mulsd  xmm1, [rsp + nb210_c6]
	mulsd  xmm2, [rsp + nb210_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb210_Vvdwtot]
	mulsd  xmm1, [rsp + nb210_six]
	mulsd  xmm2, [rsp + nb210_twelve]
	subsd  xmm2, xmm1
	subsd  xmm0, xmm7
	mulsd  xmm3, xmm0
	addsd  xmm2, xmm3
	mulsd  xmm4, xmm2	;# xmm4=total fscal 
	addsd  xmm6, [rsp + nb210_vctot]

	movsd [rsp + nb210_vctot], xmm6
	movsd [rsp + nb210_Vvdwtot], xmm5

	mov    rdi, [rbp + nb210_faction]
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
	
.nb210_updateouterdata:
	mov   ecx, [rsp + nb210_ii3]
	mov   rdi, [rbp + nb210_faction]
	mov   rsi, [rbp + nb210_fshift]
	mov   edx, [rsp + nb210_is3]

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
	mov esi, [rsp + nb210_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb210_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb210_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb210_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb210_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb210_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb210_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb210_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210_n], esi
        jmp .nb210_outer
.nb210_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210_end
        ;# non-zero, do one more workunit
        jmp   .nb210_threadloop
.nb210_end:
	mov eax, [rsp + nb210_nouter]
	mov ebx, [rsp + nb210_ninner]
	mov rcx, [rbp + nb210_outeriter]
	mov rdx, [rbp + nb210_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 448
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

	


.globl nb_kernel210nf_x86_64_sse2
.globl _nb_kernel210nf_x86_64_sse2
nb_kernel210nf_x86_64_sse2:	
_nb_kernel210nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb210nf_fshift,         16
.equiv          nb210nf_gid,            24
.equiv          nb210nf_pos,            32
.equiv          nb210nf_faction,        40
.equiv          nb210nf_charge,         48
.equiv          nb210nf_p_facel,        56
.equiv          nb210nf_argkrf,         64
.equiv          nb210nf_argcrf,         72
.equiv          nb210nf_Vc,             80
.equiv          nb210nf_type,           88
.equiv          nb210nf_p_ntype,        96
.equiv          nb210nf_vdwparam,       104
.equiv          nb210nf_Vvdw,           112
.equiv          nb210nf_p_tabscale,     120
.equiv          nb210nf_VFtab,          128
.equiv          nb210nf_invsqrta,       136
.equiv          nb210nf_dvda,           144
.equiv          nb210nf_p_gbtabscale,   152
.equiv          nb210nf_GBtab,          160
.equiv          nb210nf_p_nthreads,     168
.equiv          nb210nf_count,          176
.equiv          nb210nf_mtx,            184
.equiv          nb210nf_outeriter,      192
.equiv          nb210nf_inneriter,      200
.equiv          nb210nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb210nf_ix,             0
.equiv          nb210nf_iy,             16
.equiv          nb210nf_iz,             32
.equiv          nb210nf_iq,             48
.equiv          nb210nf_c6,             64
.equiv          nb210nf_c12,            80
.equiv          nb210nf_vctot,          96
.equiv          nb210nf_Vvdwtot,        112
.equiv          nb210nf_half,           128
.equiv          nb210nf_three,          144
.equiv          nb210nf_krf,            160
.equiv          nb210nf_crf,            176
.equiv          nb210nf_nri,            192
.equiv          nb210nf_iinr,           200
.equiv          nb210nf_jindex,         208
.equiv          nb210nf_jjnr,           216
.equiv          nb210nf_shift,          224
.equiv          nb210nf_shiftvec,       232
.equiv          nb210nf_facel,          240
.equiv          nb210nf_innerjjnr,      248
.equiv          nb210nf_is3,            256
.equiv          nb210nf_ii3,            260
.equiv          nb210nf_ntia,           264
.equiv          nb210nf_innerk,         268
.equiv          nb210nf_n,              272
.equiv          nb210nf_nn1,            276
.equiv          nb210nf_ntype,          280
.equiv          nb210nf_nouter,         284
.equiv          nb210nf_ninner,         288

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
	mov [rsp + nb210nf_nouter], eax
	mov [rsp + nb210nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb210nf_nri], edi
	mov [rsp + nb210nf_iinr], rsi
	mov [rsp + nb210nf_jindex], rdx
	mov [rsp + nb210nf_jjnr], rcx
	mov [rsp + nb210nf_shift], r8
	mov [rsp + nb210nf_shiftvec], r9
	mov rdi, [rbp + nb210nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb210nf_ntype], edi
	mov rsi, [rbp + nb210nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb210nf_facel], xmm0

	mov rsi, [rbp + nb210nf_argkrf]
	mov rdi, [rbp + nb210nf_argcrf]
	movsd xmm1, [rsi]
	movsd xmm2, [rdi]
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	movapd [rsp + nb210nf_krf], xmm1
	movapd [rsp + nb210nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb210nf_half], eax
	mov [rsp + nb210nf_half+4], ebx
	movsd xmm1, [rsp + nb210nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb210nf_half], xmm1
	movapd [rsp + nb210nf_three], xmm3

.nb210nf_threadloop:
        mov   rsi, [rbp + nb210nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb210nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb210nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb210nf_n], eax
        mov [rsp + nb210nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210nf_outerstart
        jmp .nb210nf_end

.nb210nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb210nf_nouter]
	mov [rsp + nb210nf_nouter], ebx

.nb210nf_outer:
	mov   rax, [rsp + nb210nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb210nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb210nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb210nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb210nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb210nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb210nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb210nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb210nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb210nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb210nf_ix], xmm0
	movapd [rsp + nb210nf_iy], xmm1
	movapd [rsp + nb210nf_iz], xmm2

	mov   [rsp + nb210nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb210nf_vctot], xmm4
	movapd [rsp + nb210nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb210nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb210nf_pos]
	mov   rax, [rsp + nb210nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb210nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb210nf_ninner]
	mov   [rsp + nb210nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb210nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210nf_unroll_loop
	jmp   .nb210nf_checksingle
.nb210nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb210nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb210nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb210nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm5, [rsp + nb210nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb210nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb210nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb210nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movd  ebx, mm1
	movapd [rsp + nb210nf_c6], xmm4
	movapd [rsp + nb210nf_c12], xmm6
	
	mov rsi, [rbp + nb210nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb210nf_ix]
	movapd xmm5, [rsp + nb210nf_iy]
	movapd xmm6, [rsp + nb210nf_iz]

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

	movapd xmm7, [rsp + nb210nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb210nf_three]
	mulpd xmm7, xmm4	;# krsq 
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb210nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addpd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subpd  xmm6, [rsp + nb210nf_crf]
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulpd  xmm1, [rsp + nb210nf_c6]
	mulpd  xmm2, [rsp + nb210nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [rsp + nb210nf_Vvdwtot]
	addpd  xmm6, [rsp + nb210nf_vctot]
	movapd [rsp + nb210nf_vctot], xmm6
	movapd [rsp + nb210nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb210nf_innerk],  2
	jl    .nb210nf_checksingle
	jmp   .nb210nf_unroll_loop

.nb210nf_checksingle:				
	mov   edx, [rsp + nb210nf_innerk]
	and   edx, 1
	jnz    .nb210nf_dosingle
	jmp    .nb210nf_updateouterdata
.nb210nf_dosingle:			
	mov rsi, [rbp + nb210nf_charge]
	mov rdi, [rbp + nb210nf_pos]
	mov   rcx, [rsp + nb210nf_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [rcx]

	movlpd xmm3, [rsi + rax*8]
	movapd xmm5, [rsp + nb210nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov rsi, [rbp + nb210nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb210nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb210nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [rsp + nb210nf_c6], xmm4
	movapd [rsp + nb210nf_c12], xmm6
	
	mov rsi, [rbp + nb210nf_pos]       ;# base of pos[] 

	lea rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb210nf_ix]
	movapd xmm5, [rsp + nb210nf_iy]
	movapd xmm6, [rsp + nb210nf_iz]

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

	movapd xmm7, [rsp + nb210nf_krf]	
	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb210nf_three]
	mulsd xmm7, xmm4	;# krsq 
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb210nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb210nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm4	;# xmm4=rinvsq 
	movapd xmm6, xmm0
	addsd  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movapd xmm1, xmm4
	subsd  xmm6, [rsp + nb210nf_crf]
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulsd  xmm1, [rsp + nb210nf_c6]
	mulsd  xmm2, [rsp + nb210nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [rsp + nb210nf_Vvdwtot]
	addsd  xmm6, [rsp + nb210nf_vctot]
	movlpd [rsp + nb210nf_vctot], xmm6
	movlpd [rsp + nb210nf_Vvdwtot], xmm5
	
.nb210nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb210nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb210nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb210nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb210nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb210nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb210nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb210nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb210nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210nf_n], esi
        jmp .nb210nf_outer
.nb210nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210nf_end
        ;# non-zero, do one more workunit
        jmp   .nb210nf_threadloop
.nb210nf_end:
	mov eax, [rsp + nb210nf_nouter]
	mov ebx, [rsp + nb210nf_ninner]
	mov rcx, [rbp + nb210nf_outeriter]
	mov rdx, [rbp + nb210nf_inneriter]
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
