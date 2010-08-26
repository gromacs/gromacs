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



.globl nb_kernel211_x86_64_sse
.globl _nb_kernel211_x86_64_sse
nb_kernel211_x86_64_sse:	
_nb_kernel211_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb211_fshift,           16
.equiv          nb211_gid,              24
.equiv          nb211_pos,              32
.equiv          nb211_faction,          40
.equiv          nb211_charge,           48
.equiv          nb211_p_facel,          56
.equiv          nb211_argkrf,           64
.equiv          nb211_argcrf,           72
.equiv          nb211_Vc,               80
.equiv          nb211_type,             88
.equiv          nb211_p_ntype,          96
.equiv          nb211_vdwparam,         104
.equiv          nb211_Vvdw,             112
.equiv          nb211_p_tabscale,       120
.equiv          nb211_VFtab,            128
.equiv          nb211_invsqrta,         136
.equiv          nb211_dvda,             144
.equiv          nb211_p_gbtabscale,     152
.equiv          nb211_GBtab,            160
.equiv          nb211_p_nthreads,       168
.equiv          nb211_count,            176
.equiv          nb211_mtx,              184
.equiv          nb211_outeriter,        192
.equiv          nb211_inneriter,        200
.equiv          nb211_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb211_ixO,              0
.equiv          nb211_iyO,              16
.equiv          nb211_izO,              32
.equiv          nb211_ixH1,             48
.equiv          nb211_iyH1,             64
.equiv          nb211_izH1,             80
.equiv          nb211_ixH2,             96
.equiv          nb211_iyH2,             112
.equiv          nb211_izH2,             128
.equiv          nb211_iqO,              144
.equiv          nb211_iqH,              160
.equiv          nb211_dxO,              176
.equiv          nb211_dyO,              192
.equiv          nb211_dzO,              208
.equiv          nb211_dxH1,             224
.equiv          nb211_dyH1,             240
.equiv          nb211_dzH1,             256
.equiv          nb211_dxH2,             272
.equiv          nb211_dyH2,             288
.equiv          nb211_dzH2,             304
.equiv          nb211_qqO,              320
.equiv          nb211_qqH,              336
.equiv          nb211_c6,               352
.equiv          nb211_c12,              368
.equiv          nb211_six,              384
.equiv          nb211_twelve,           400
.equiv          nb211_vctot,            416
.equiv          nb211_Vvdwtot,          432
.equiv          nb211_fixO,             448
.equiv          nb211_fiyO,             464
.equiv          nb211_fizO,             480
.equiv          nb211_fixH1,            496
.equiv          nb211_fiyH1,            512
.equiv          nb211_fizH1,            528
.equiv          nb211_fixH2,            544
.equiv          nb211_fiyH2,            560
.equiv          nb211_fizH2,            576
.equiv          nb211_fjx,              592
.equiv          nb211_fjy,              608
.equiv          nb211_fjz,              624
.equiv          nb211_half,             640
.equiv          nb211_three,            656
.equiv          nb211_two,              672
.equiv          nb211_krf,              688
.equiv          nb211_crf,              704
.equiv          nb211_krsqO,            720
.equiv          nb211_krsqH1,           736
.equiv          nb211_krsqH2,           752
.equiv          nb211_nri,              768
.equiv          nb211_iinr,             776
.equiv          nb211_jindex,           784
.equiv          nb211_jjnr,             792
.equiv          nb211_shift,            800
.equiv          nb211_shiftvec,         808
.equiv          nb211_facel,            816
.equiv          nb211_innerjjnr,        824
.equiv          nb211_is3,              832
.equiv          nb211_ii3,              836
.equiv          nb211_ntia,             840
.equiv          nb211_innerk,           844
.equiv          nb211_n,                848
.equiv          nb211_nn1,              852
.equiv          nb211_nouter,           856
.equiv          nb211_ninner,           860

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
    
	emms
	sub rsp, 864		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb211_nouter], eax
	mov [rsp + nb211_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb211_nri], edi
	mov [rsp + nb211_iinr], rsi
	mov [rsp + nb211_jindex], rdx
	mov [rsp + nb211_jjnr], rcx
	mov [rsp + nb211_shift], r8
	mov [rsp + nb211_shiftvec], r9
	mov rsi, [rbp + nb211_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb211_facel], xmm0


	mov rsi, [rbp + nb211_argkrf]
	mov rdi, [rbp + nb211_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb211_krf], xmm1
	movaps [rsp + nb211_crf], xmm2


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb211_half], eax
	movss xmm1, [rsp + nb211_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + nb211_half],  xmm1
	movaps [rsp + nb211_two],  xmm2
	movaps [rsp + nb211_three],  xmm3
	movaps [rsp + nb211_six],  xmm4
	movaps [rsp + nb211_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb211_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb211_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb211_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb211_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb211_iqO], xmm3
	movaps [rsp + nb211_iqH], xmm4
	
	mov   rdx, [rbp + nb211_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb211_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb211_ntia], ecx		

.nb211_threadloop:
        mov   rsi, [rbp + nb211_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb211_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb211_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb211_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb211_n], eax
        mov [rsp + nb211_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb211_outerstart
        jmp .nb211_end

.nb211_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb211_nouter]
	mov [rsp + nb211_nouter], ebx

.nb211_outer:
	mov   rax, [rsp + nb211_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax +rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb211_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb211_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb211_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb211_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb211_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb211_ixO], xmm3
	movaps [rsp + nb211_iyO], xmm4
	movaps [rsp + nb211_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb211_ixH1], xmm0
	movaps [rsp + nb211_iyH1], xmm1
	movaps [rsp + nb211_izH1], xmm2
	movaps [rsp + nb211_ixH2], xmm3
	movaps [rsp + nb211_iyH2], xmm4
	movaps [rsp + nb211_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb211_vctot], xmm4
	movaps [rsp + nb211_Vvdwtot], xmm4
	movaps [rsp + nb211_fixO], xmm4
	movaps [rsp + nb211_fiyO], xmm4
	movaps [rsp + nb211_fizO], xmm4
	movaps [rsp + nb211_fixH1], xmm4
	movaps [rsp + nb211_fiyH1], xmm4
	movaps [rsp + nb211_fizH1], xmm4
	movaps [rsp + nb211_fixH2], xmm4
	movaps [rsp + nb211_fiyH2], xmm4
	movaps [rsp + nb211_fizH2], xmm4
	
	mov   rax, [rsp + nb211_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb211_pos]
	mov   rdi, [rbp + nb211_faction]	
	mov   rax, [rsp + nb211_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb211_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb211_ninner]
	mov   [rsp + nb211_ninner], ecx
	add   edx, 0
	mov   [rsp + nb211_innerk], edx    ;# number of innerloop atoms 

	jge   .nb211_unroll_loop
	jmp   .nb211_odd_inner
.nb211_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb211_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb211_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb211_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb211_iqO]
	mulps  xmm4, [rsp + nb211_iqH]

	movaps  [rsp + nb211_qqO], xmm3
	movaps  [rsp + nb211_qqH], xmm4
	
	mov rsi, [rbp + nb211_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov r10d, [rsi + rcx*4]
	mov r11d, [rsi + rdx*4]
	mov rsi, [rbp + nb211_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
	mov edi, [rsp + nb211_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	movlps xmm6, [rsi + r8*4]
	movlps xmm7, [rsi + r10*4]
	movhps xmm6, [rsi + r9*4]
	movhps xmm7, [rsi + r11*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101

	movaps [rsp + nb211_c6], xmm4
	movaps [rsp + nb211_c12], xmm6

	mov rsi, [rbp + nb211_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four j coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb211_ixO]
    subps xmm1, [rsp + nb211_iyO]
    subps xmm2, [rsp + nb211_izO]
    subps xmm3, [rsp + nb211_ixH1]
    subps xmm4, [rsp + nb211_iyH1]
    subps xmm5, [rsp + nb211_izH1]
    subps xmm6, [rsp + nb211_ixH2]
    subps xmm7, [rsp + nb211_iyH2]
    subps xmm8, [rsp + nb211_izH2]
    
	movaps [rsp + nb211_dxO], xmm0
	movaps [rsp + nb211_dyO], xmm1
	movaps [rsp + nb211_dzO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb211_dxH1], xmm3
	movaps [rsp + nb211_dyH1], xmm4
	movaps [rsp + nb211_dzH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb211_dxH2], xmm6
	movaps [rsp + nb211_dyH2], xmm7
	movaps [rsp + nb211_dzH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb211_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb211_half]
	mulps   xmm9, xmm4  ;# rinvO
	mulps   xmm10, xmm4 ;# rinvH1
    mulps   xmm11, xmm4 ;# rinvH2
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb211_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm12, xmm9
    mulps  xmm12, xmm12 ;# rinv4
    mulps  xmm12, xmm9  ;# rinv6
    subps  xmm2, [rsp + nb211_crf]   ;# rinv+krsq-crf
    subps  xmm5, [rsp + nb211_crf]
    subps  xmm8, [rsp + nb211_crf]   
    mulps  xmm2, [rsp + nb211_qqO] ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, [rsp + nb211_qqH] ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, [rsp + nb211_qqH] ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    movaps xmm13, xmm12 ;# rinv6
    mulps xmm12, xmm12 ;# rinv12
	mulps  xmm13, [rsp + nb211_c6]
	mulps  xmm12, [rsp + nb211_c12]
    movaps xmm14, xmm12
    subps  xmm14, xmm13
    mulps  xmm1, [rsp + nb211_qqO]   ;# (rinv-2*krsq)*qq
    mulps  xmm4, [rsp + nb211_qqH] 
    mulps  xmm7, [rsp + nb211_qqH] 
    addps  xmm2, [rsp + nb211_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps [rsp + nb211_vctot], xmm2
 
   	addps  xmm14, [rsp + nb211_Vvdwtot]
	mulps  xmm13, [rsp + nb211_six]
	mulps  xmm12, [rsp + nb211_twelve]
	movaps [rsp + nb211_Vvdwtot], xmm14
    subps  xmm12, xmm13 ;# LJ fscal        

    addps xmm1, xmm12
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	mov   rdi, [rbp + nb211_faction]	
	;# move j  forces to local temp variables 
    movlps xmm9, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm10, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm9, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm10, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm11, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm12, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm6,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm8,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm11, xmm6 ;# jza  -  jzb  -
    movlhps xmm12, xmm8 ;# jzc  -  jzd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm9: jxa jya jxb jyb 
    ;# xmm10: jxc jyc jxd jyd
    ;# xmm11: jza jzb jzc jzd

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb211_dxO]
	mulps xmm1, [rsp + nb211_dyO]
	mulps xmm2, [rsp + nb211_dzO]
	mulps xmm3, [rsp + nb211_dxH1]
	mulps xmm4, [rsp + nb211_dyH1]
	mulps xmm5, [rsp + nb211_dzH1]
	mulps xmm6, [rsp + nb211_dxH2]
	mulps xmm7, [rsp + nb211_dyH2]
	mulps xmm8, [rsp + nb211_dzH2]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb211_fixO]
    addps xmm1, [rsp + nb211_fiyO]
    addps xmm2, [rsp + nb211_fizO]

    addps xmm13, xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb211_fixH1]
    addps xmm4, [rsp + nb211_fiyH1]
    addps xmm5, [rsp + nb211_fizH1]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb211_fixH2]
    addps xmm7, [rsp + nb211_fiyH2]
    addps xmm8, [rsp + nb211_fizH2]

    movaps [rsp + nb211_fixO], xmm0
    movaps [rsp + nb211_fiyO], xmm1
    movaps [rsp + nb211_fizO], xmm2
    movaps [rsp + nb211_fixH1], xmm3
    movaps [rsp + nb211_fiyH1], xmm4
    movaps [rsp + nb211_fizH1], xmm5
    movaps [rsp + nb211_fixH2], xmm6
    movaps [rsp + nb211_fiyH2], xmm7
    movaps [rsp + nb211_fizH2], xmm8
    
    ;# xmm9 = fjx
    ;# xmm10 = fjy
    ;# xmm11 = fjz
    movaps xmm15, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm15, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm15

    movhlps  xmm12, xmm11
    
    movlps [rdi + rax*4], xmm9
    movhps [rdi + rbx*4], xmm9
    movlps [rdi + rcx*4], xmm10
    movhps [rdi + rdx*4], xmm10
    movss  [rdi + rax*4 + 8], xmm11
    movss  [rdi + rcx*4 + 8], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 8], xmm11
    movss  [rdi + rdx*4 + 8], xmm12

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb211_innerk],  4
	jl    .nb211_odd_inner
	jmp   .nb211_unroll_loop
.nb211_odd_inner:	
	add dword ptr [rsp + nb211_innerk],  4
	jnz   .nb211_odd_loop
	jmp   .nb211_updateouterdata
.nb211_odd_loop:
	mov   rdx, [rsp + nb211_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb211_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb211_iqO]
	mov rsi, [rbp + nb211_charge] 
	movhps xmm4, [rsp + nb211_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb211_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb211_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb211_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb211_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb211_c6], xmm6
	movaps [rsp + nb211_c12], xmm7
	mov rsi, [rbp + nb211_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb211_ixO]
	movss xmm1, [rsp + nb211_iyO]
	movss xmm2, [rsp + nb211_izO]
	
	movlps xmm6, [rsp + nb211_ixH1]
	movlps xmm7, [rsp + nb211_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb211_iyH1]
	movlps xmm7, [rsp + nb211_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb211_izH1]
	movlps xmm7, [rsp + nb211_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb211_dxO], xmm3
	movaps [rsp + nb211_dyO], xmm4
	movaps [rsp + nb211_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb211_krf]
	movaps [rsp + nb211_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb211_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb211_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000
	
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb211_c6]
	mulps  xmm2, [rsp + nb211_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb211_Vvdwtot]
	mulss  xmm1, [rsp + nb211_six]
	mulss  xmm2, [rsp + nb211_twelve]
	subss  xmm2, xmm1

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [rsp + nb211_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	mulps  xmm3, [rsp + nb211_two]
	subps  xmm0, [rsp + nb211_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [rsp + nb211_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [rsp + nb211_qqO] 	;# xmm1=coul part of fs 

	addps xmm2, xmm1	;# total fs 
	
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm0, [rsp + nb211_vctot]
	movaps [rsp + nb211_vctot], xmm0
	
	movaps xmm0, [rsp + nb211_dxO]
	movaps xmm1, [rsp + nb211_dyO]
	movaps xmm2, [rsp + nb211_dzO]

	movaps [rsp + nb211_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb211_fixO]	
	movss  xmm4, [rsp + nb211_fiyO]	
	movss  xmm5, [rsp + nb211_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb211_fixO], xmm3	
	movss  [rsp + nb211_fiyO], xmm4	
	movss  [rsp + nb211_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb211_fixH1]
	addss  xmm4, [rsp + nb211_fiyH1]
	addss  xmm5, [rsp + nb211_fizH1]
	movss  [rsp + nb211_fixH1], xmm3	
	movss  [rsp + nb211_fiyH1], xmm4	
	movss  [rsp + nb211_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb211_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb211_fixH2]
	addss  xmm4, [rsp + nb211_fiyH2]
	addss  xmm5, [rsp + nb211_fizH2]
	movss  [rsp + nb211_fixH2], xmm3	
	movss  [rsp + nb211_fiyH2], xmm4	
	movss  [rsp + nb211_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1    ;# z sum in xmm2 
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb211_innerk]
	jz    .nb211_updateouterdata
	jmp   .nb211_odd_loop
.nb211_updateouterdata:
	mov   ecx, [rsp + nb211_ii3]
	mov   rdi, [rbp + nb211_faction]
	mov   rsi, [rbp + nb211_fshift]
	mov   edx, [rsp + nb211_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb211_fixO]
	movaps xmm1, [rsp + nb211_fiyO]
	movaps xmm2, [rsp + nb211_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb211_fixH1]
	movaps xmm1, [rsp + nb211_fiyH1]
	movaps xmm2, [rsp + nb211_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb211_fixH2]
	movaps xmm1, [rsp + nb211_fiyH2]
	movaps xmm2, [rsp + nb211_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb211_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb211_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb211_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb211_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb211_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb211_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb211_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb211_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb211_n], esi
        jmp .nb211_outer
.nb211_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb211_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb211_end
        ;# non-zero, do one more workunit
        jmp   .nb211_threadloop
.nb211_end:
	mov eax, [rsp + nb211_nouter]
	mov ebx, [rsp + nb211_ninner]
	mov rcx, [rbp + nb211_outeriter]
	mov rdx, [rbp + nb211_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 864
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



.globl nb_kernel211nf_x86_64_sse
.globl _nb_kernel211nf_x86_64_sse
nb_kernel211nf_x86_64_sse:	
_nb_kernel211nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb211nf_fshift,         16
.equiv          nb211nf_gid,            24
.equiv          nb211nf_pos,            32
.equiv          nb211nf_faction,        40
.equiv          nb211nf_charge,         48
.equiv          nb211nf_p_facel,        56
.equiv          nb211nf_argkrf,         64
.equiv          nb211nf_argcrf,         72
.equiv          nb211nf_Vc,             80
.equiv          nb211nf_type,           88
.equiv          nb211nf_p_ntype,        96
.equiv          nb211nf_vdwparam,       104
.equiv          nb211nf_Vvdw,           112
.equiv          nb211nf_p_tabscale,     120
.equiv          nb211nf_VFtab,          128
.equiv          nb211nf_invsqrta,       136
.equiv          nb211nf_dvda,           144
.equiv          nb211nf_p_gbtabscale,   152
.equiv          nb211nf_GBtab,          160
.equiv          nb211nf_p_nthreads,     168
.equiv          nb211nf_count,          176
.equiv          nb211nf_mtx,            184
.equiv          nb211nf_outeriter,      192
.equiv          nb211nf_inneriter,      200
.equiv          nb211nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb211nf_ixO,            0
.equiv          nb211nf_iyO,            16
.equiv          nb211nf_izO,            32
.equiv          nb211nf_ixH1,           48
.equiv          nb211nf_iyH1,           64
.equiv          nb211nf_izH1,           80
.equiv          nb211nf_ixH2,           96
.equiv          nb211nf_iyH2,           112
.equiv          nb211nf_izH2,           128
.equiv          nb211nf_iqO,            144
.equiv          nb211nf_iqH,            160
.equiv          nb211nf_qqO,            176
.equiv          nb211nf_qqH,            192
.equiv          nb211nf_c6,             208
.equiv          nb211nf_c12,            224
.equiv          nb211nf_vctot,          240
.equiv          nb211nf_Vvdwtot,        256
.equiv          nb211nf_half,           272
.equiv          nb211nf_three,          288
.equiv          nb211nf_krf,            304
.equiv          nb211nf_crf,            320
.equiv          nb211nf_krsqO,          336
.equiv          nb211nf_krsqH1,         352
.equiv          nb211nf_krsqH2,         368
.equiv          nb211nf_nri,            384
.equiv          nb211nf_iinr,           392
.equiv          nb211nf_jindex,         400
.equiv          nb211nf_jjnr,           408
.equiv          nb211nf_shift,          416
.equiv          nb211nf_shiftvec,       424
.equiv          nb211nf_facel,          432
.equiv          nb211nf_innerjjnr,      440
.equiv          nb211nf_ntia,           448
.equiv          nb211nf_is3,            452
.equiv          nb211nf_ii3,            456
.equiv          nb211nf_innerk,         460
.equiv          nb211nf_n,              464
.equiv          nb211nf_nn1,            468
.equiv          nb211nf_nouter,         472
.equiv          nb211nf_ninner,         476

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
    
	emms
	sub rsp, 480		;# local variable stack space (n*16+8)
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

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb211nf_nouter], eax
	mov [rsp + nb211nf_ninner], eax
	
	mov edi, [rdi]
	mov [rsp + nb211nf_nri], edi
	mov [rsp + nb211nf_iinr], rsi
	mov [rsp + nb211nf_jindex], rdx
	mov [rsp + nb211nf_jjnr], rcx
	mov [rsp + nb211nf_shift], r8
	mov [rsp + nb211nf_shiftvec], r9
	mov rsi, [rbp + nb211nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb211nf_facel], xmm0


	mov rsi, [rbp + nb211nf_argkrf]
	mov rdi, [rbp + nb211nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb211nf_krf], xmm1
	movaps [rsp + nb211nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb211nf_half], eax
	movss xmm1, [rsp + nb211nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb211nf_half],  xmm1
	movaps [rsp + nb211nf_three],  xmm3	

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb211nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb211nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb211nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb211nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb211nf_iqO], xmm3
	movaps [rsp + nb211nf_iqH], xmm4
	
	mov   rdx, [rbp + nb211nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb211nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb211nf_ntia], ecx		

.nb211nf_threadloop:
        mov   rsi, [rbp + nb211nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb211nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb211nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb211nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb211nf_n], eax
        mov [rsp + nb211nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb211nf_outerstart
        jmp .nb211nf_end

.nb211nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb211nf_nouter]
	mov [rsp + nb211nf_nouter], ebx

.nb211nf_outer:
	mov   rax, [rsp + nb211nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb211nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb211nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb211nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb211nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb211nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb211nf_ixO], xmm3
	movaps [rsp + nb211nf_iyO], xmm4
	movaps [rsp + nb211nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb211nf_ixH1], xmm0
	movaps [rsp + nb211nf_iyH1], xmm1
	movaps [rsp + nb211nf_izH1], xmm2
	movaps [rsp + nb211nf_ixH2], xmm3
	movaps [rsp + nb211nf_iyH2], xmm4
	movaps [rsp + nb211nf_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb211nf_vctot], xmm4
	movaps [rsp + nb211nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb211nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb211nf_pos]
	mov   rax, [rsp + nb211nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb211nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb211nf_ninner]
	mov   [rsp + nb211nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb211nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb211nf_unroll_loop
	jmp   .nb211nf_odd_inner
.nb211nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb211nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb211nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb211nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb211nf_iqO]
	mulps  xmm4, [rsp + nb211nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb211nf_qqO], xmm3
	movaps  [rsp + nb211nf_qqH], xmm4
	
	mov rsi, [rbp + nb211nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb211nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb211nf_ntia]
	add rax, rdi
	add rbx, rdi
	add rcx, rdi
	add rdx, rdi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb211nf_c6], xmm4
	movaps [rsp + nb211nf_c12], xmm6

	mov rsi, [rbp + nb211nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb211nf_ixO]
	movaps xmm5, [rsp + nb211nf_iyO]
	movaps xmm6, [rsp + nb211nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2
	
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb211nf_ixH1]
	movaps xmm5, [rsp + nb211nf_iyH1]
	movaps xmm6, [rsp + nb211nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2
	
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb211nf_ixH2]
	movaps xmm4, [rsp + nb211nf_iyH2]
	movaps xmm5, [rsp + nb211nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [rsp + nb211nf_krf]	
	mulps  xmm1, [rsp + nb211nf_krf]	
	mulps  xmm2, [rsp + nb211nf_krf]	

	movaps [rsp + nb211nf_krsqH2], xmm0
	movaps [rsp + nb211nf_krsqH1], xmm1
	movaps [rsp + nb211nf_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb211nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb211nf_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb211nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb211nf_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb211nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb211nf_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb211nf_c6]
	mulps  xmm2, [rsp + nb211nf_c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addps  xmm3, [rsp + nb211nf_Vvdwtot]

	movaps xmm0, xmm7
	movaps xmm1, [rsp + nb211nf_krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [rsp + nb211nf_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm7, xmm1
	mulps  xmm0, [rsp + nb211nf_qqO]
	addps  xmm0, [rsp + nb211nf_vctot]
	movaps [rsp + nb211nf_Vvdwtot], xmm3
	movaps [rsp + nb211nf_vctot], xmm0

	;# H1 interactions 
	movaps  xmm0, [rsp + nb211nf_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [rsp + nb211nf_crf]
	mulps   xmm6, [rsp + nb211nf_qqH] ;# vcoul 
	addps  xmm6, [rsp + nb211nf_vctot]
	movaps [rsp + nb211nf_vctot], xmm6
	
	;# H2 interactions 
	movaps  xmm7, xmm5  ;# rinv 
	movaps  xmm0, [rsp + nb211nf_krsqH2]
	addps   xmm5, xmm0	;# xmm5=rinv+ krsq 
	subps   xmm5, [rsp + nb211nf_crf]
	mulps   xmm5, [rsp + nb211nf_qqH] ;# vcoul 
	addps   xmm5, [rsp + nb211nf_vctot]
	movaps [rsp + nb211nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb211nf_innerk],  4
	jl    .nb211nf_odd_inner
	jmp   .nb211nf_unroll_loop
.nb211nf_odd_inner:	
	add dword ptr [rsp + nb211nf_innerk],  4
	jnz   .nb211nf_odd_loop
	jmp   .nb211nf_updateouterdata
.nb211nf_odd_loop:
	mov   rdx, [rsp + nb211nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb211nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb211nf_iqO]
	mov rsi, [rbp + nb211nf_charge] 
	movhps xmm4, [rsp + nb211nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb211nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb211nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb211nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb211nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb211nf_c6], xmm6
	movaps [rsp + nb211nf_c12], xmm7

	mov rsi, [rbp + nb211nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb211nf_ixO]
	movss xmm4, [rsp + nb211nf_iyO]
	movss xmm5, [rsp + nb211nf_izO]
	
	movlps xmm6, [rsp + nb211nf_ixH1]
	movlps xmm7, [rsp + nb211nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb211nf_iyH1]
	movlps xmm7, [rsp + nb211nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb211nf_izH1]
	movlps xmm7, [rsp + nb211nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb211nf_krf]
	movaps [rsp + nb211nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb211nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb211nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	

	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb211nf_c6]
	mulps  xmm2, [rsp + nb211nf_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb211nf_Vvdwtot]
	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [rsp + nb211nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [rsp + nb211nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [rsp + nb211nf_qqO]	;# xmm0=vcoul 
	addps  xmm0, [rsp + nb211nf_vctot]
	movaps [rsp + nb211nf_vctot], xmm0
	movaps [rsp + nb211nf_Vvdwtot], xmm5

	dec dword ptr [rsp + nb211nf_innerk]
	jz    .nb211nf_updateouterdata
	jmp   .nb211nf_odd_loop
.nb211nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb211nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb211nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb211nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb211nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb211nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb211nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb211nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb211nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb211nf_n], esi
        jmp .nb211nf_outer
.nb211nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb211nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb211nf_end
        ;# non-zero, do one more workunit
        jmp   .nb211nf_threadloop
.nb211nf_end:
	mov eax, [rsp + nb211nf_nouter]
	mov ebx, [rsp + nb211nf_ninner]
	mov rcx, [rbp + nb211nf_outeriter]
	mov rdx, [rbp + nb211nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 480
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

