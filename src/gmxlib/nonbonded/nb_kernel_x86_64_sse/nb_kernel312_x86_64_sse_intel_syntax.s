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


.globl nb_kernel312_x86_64_sse
.globl _nb_kernel312_x86_64_sse
nb_kernel312_x86_64_sse:	
_nb_kernel312_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb312_fshift,           16
.equiv          nb312_gid,              24
.equiv          nb312_pos,              32
.equiv          nb312_faction,          40
.equiv          nb312_charge,           48
.equiv          nb312_p_facel,          56
.equiv          nb312_argkrf,           64
.equiv          nb312_argcrf,           72
.equiv          nb312_Vc,               80
.equiv          nb312_type,             88
.equiv          nb312_p_ntype,          96
.equiv          nb312_vdwparam,         104
.equiv          nb312_Vvdw,             112
.equiv          nb312_p_tabscale,       120
.equiv          nb312_VFtab,            128
.equiv          nb312_invsqrta,         136
.equiv          nb312_dvda,             144
.equiv          nb312_p_gbtabscale,     152
.equiv          nb312_GBtab,            160
.equiv          nb312_p_nthreads,       168
.equiv          nb312_count,            176
.equiv          nb312_mtx,              184
.equiv          nb312_outeriter,        192
.equiv          nb312_inneriter,        200
.equiv          nb312_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb312_ixO,              0
.equiv          nb312_iyO,              16
.equiv          nb312_izO,              32
.equiv          nb312_ixH1,             48
.equiv          nb312_iyH1,             64
.equiv          nb312_izH1,             80
.equiv          nb312_ixH2,             96
.equiv          nb312_iyH2,             112
.equiv          nb312_izH2,             128
.equiv          nb312_jxO,              144
.equiv          nb312_jyO,              160
.equiv          nb312_jzO,              176
.equiv          nb312_jxH1,             192
.equiv          nb312_jyH1,             208
.equiv          nb312_jzH1,             224
.equiv          nb312_jxH2,             240
.equiv          nb312_jyH2,             256
.equiv          nb312_jzH2,             272
.equiv          nb312_dxOO,             288
.equiv          nb312_dyOO,             304
.equiv          nb312_dzOO,             320
.equiv          nb312_dxOH1,            336
.equiv          nb312_dyOH1,            352
.equiv          nb312_dzOH1,            368
.equiv          nb312_dxOH2,            384
.equiv          nb312_dyOH2,            400
.equiv          nb312_dzOH2,            416
.equiv          nb312_dxH1O,            432
.equiv          nb312_dyH1O,            448
.equiv          nb312_dzH1O,            464
.equiv          nb312_dxH1H1,           480
.equiv          nb312_dyH1H1,           496
.equiv          nb312_dzH1H1,           512
.equiv          nb312_dxH1H2,           528
.equiv          nb312_dyH1H2,           544
.equiv          nb312_dzH1H2,           560
.equiv          nb312_dxH2O,            576
.equiv          nb312_dyH2O,            592
.equiv          nb312_dzH2O,            608
.equiv          nb312_dxH2H1,           624
.equiv          nb312_dyH2H1,           640
.equiv          nb312_dzH2H1,           656
.equiv          nb312_dxH2H2,           672
.equiv          nb312_dyH2H2,           688
.equiv          nb312_dzH2H2,           704
.equiv          nb312_qqOO,             720
.equiv          nb312_qqOH,             736
.equiv          nb312_qqHH,             752
.equiv          nb312_two,              768
.equiv          nb312_tsc,              784
.equiv          nb312_c6,               800
.equiv          nb312_c12,              816
.equiv          nb312_six,              832
.equiv          nb312_twelve,           848
.equiv          nb312_vctot,            864
.equiv          nb312_Vvdwtot,          880
.equiv          nb312_fixO,             896
.equiv          nb312_fiyO,             912
.equiv          nb312_fizO,             928
.equiv          nb312_fixH1,            944
.equiv          nb312_fiyH1,            960
.equiv          nb312_fizH1,            976
.equiv          nb312_fixH2,            992
.equiv          nb312_fiyH2,            1008
.equiv          nb312_fizH2,            1024
.equiv          nb312_fjxO,             1040
.equiv          nb312_fjyO,             1056
.equiv          nb312_fjzO,             1072
.equiv          nb312_fjxH1,            1088
.equiv          nb312_fjyH1,            1104
.equiv          nb312_fjzH1,            1120
.equiv          nb312_fjxH2,            1136
.equiv          nb312_fjyH2,            1152
.equiv          nb312_fjzH2,            1168
.equiv          nb312_half,             1184
.equiv          nb312_three,            1200
.equiv          nb312_epsO,             1216
.equiv          nb312_epsH1,            1232
.equiv          nb312_epsH2,            1248
.equiv          nb312_rsqH1O,           1264
.equiv          nb312_rsqH1H1,          1280
.equiv          nb312_rsqH1H2,          1296
.equiv          nb312_rsqH2O,           1312
.equiv          nb312_rsqH2H1,          1328
.equiv          nb312_rsqH2H2,          1344
.equiv          nb312_rinvOO,           1360
.equiv          nb312_rinvOH1,          1376
.equiv          nb312_rinvOH2,          1392
.equiv          nb312_rinvH1O,          1408
.equiv          nb312_rinvH1H1,         1424
.equiv          nb312_rinvH1H2,         1440
.equiv          nb312_rinvH2O,          1456
.equiv          nb312_rinvH2H1,         1472
.equiv          nb312_rinvH2H2,         1488
.equiv          nb312_fstmp,            1504
.equiv          nb312_is3,              1520
.equiv          nb312_ii3,              1524
.equiv          nb312_nri,              1528
.equiv          nb312_iinr,             1536
.equiv          nb312_jindex,           1544
.equiv          nb312_jjnr,             1552
.equiv          nb312_shift,            1560
.equiv          nb312_shiftvec,         1568
.equiv          nb312_facel,            1576
.equiv          nb312_innerjjnr,        1584
.equiv          nb312_innerk,           1592
.equiv          nb312_n,                1596
.equiv          nb312_nn1,              1600
.equiv          nb312_nouter,           1604
.equiv          nb312_ninner,           1608
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
	sub rsp, 1616		;# local variable stack space (n*16+8)
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
	mov [rsp + nb312_nouter], eax
	mov [rsp + nb312_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb312_nri], edi
	mov [rsp + nb312_iinr], rsi
	mov [rsp + nb312_jindex], rdx
	mov [rsp + nb312_jjnr], rcx
	mov [rsp + nb312_shift], r8
	mov [rsp + nb312_shiftvec], r9
	mov rsi, [rbp + nb312_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb312_facel], xmm0

	mov rax, [rbp + nb312_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb312_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb312_half], eax
	movss xmm1, [rsp + nb312_half]
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
	movaps [rsp + nb312_half],  xmm1
	movaps [rsp + nb312_two],  xmm2
	movaps [rsp + nb312_three],  xmm3
	movaps [rsp + nb312_six],  xmm4
	movaps [rsp + nb312_twelve],  xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb312_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb312_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb312_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb312_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb312_qqOO], xmm3
	movaps [rsp + nb312_qqOH], xmm4
	movaps [rsp + nb312_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb312_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb312_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb312_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + nb312_c6], xmm0
	movaps [rsp + nb312_c12], xmm1

.nb312_threadloop:
        mov   rsi, [rbp + nb312_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb312_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb312_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb312_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb312_n], eax
        mov [rsp + nb312_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb312_outerstart
        jmp .nb312_end
.nb312_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb312_nouter]
	mov [rsp + nb312_nouter], ebx

.nb312_outer:
	mov   rax, [rsp + nb312_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb312_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb312_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb312_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb312_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb312_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb312_ixO], xmm3
	movaps [rsp + nb312_iyO], xmm4
	movaps [rsp + nb312_izO], xmm5

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
	movaps [rsp + nb312_ixH1], xmm0
	movaps [rsp + nb312_iyH1], xmm1
	movaps [rsp + nb312_izH1], xmm2
	movaps [rsp + nb312_ixH2], xmm3
	movaps [rsp + nb312_iyH2], xmm4
	movaps [rsp + nb312_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb312_vctot], xmm4
	movaps [rsp + nb312_Vvdwtot], xmm4
	movaps [rsp + nb312_fixO], xmm4
	movaps [rsp + nb312_fiyO], xmm4
	movaps [rsp + nb312_fizO], xmm4
	movaps [rsp + nb312_fixH1], xmm4
	movaps [rsp + nb312_fiyH1], xmm4
	movaps [rsp + nb312_fizH1], xmm4
	movaps [rsp + nb312_fixH2], xmm4
	movaps [rsp + nb312_fiyH2], xmm4
	movaps [rsp + nb312_fizH2], xmm4
	
	mov   rax, [rsp + nb312_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb312_pos]
	mov   rdi, [rbp + nb312_faction]	
	mov   rax, [rsp + nb312_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb312_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb312_ninner]
	mov   [rsp + nb312_ninner], ecx
	add   edx, 0
	mov   [rsp + nb312_innerk], edx    ;# number of innerloop atoms 
	jge   .nb312_unroll_loop
	jmp   .nb312_single_check
.nb312_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb312_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb312_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb312_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# load j O coordinates
    movlps xmm0, [rsi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, [rsi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, [rsi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, [rsi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, [rsi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, [rsi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, [rsi + rdx*4 + 8] ;# jzOc  -  jzOd -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxOa jxOc jyOa jyOc        
    unpckhps xmm4, xmm1  ;# jxOb jxOd jyOb jyOd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb312_ixO]
    subps xmm1, [rsp + nb312_iyO]
    subps xmm2, [rsp + nb312_izO]
    subps xmm3, [rsp + nb312_ixH1]
    subps xmm4, [rsp + nb312_iyH1]
    subps xmm5, [rsp + nb312_izH1]
    subps xmm6, [rsp + nb312_ixH2]
    subps xmm7, [rsp + nb312_iyH2]
    subps xmm8, [rsp + nb312_izH2]
    
	movaps [rsp + nb312_dxOO], xmm0
	movaps [rsp + nb312_dyOO], xmm1
	movaps [rsp + nb312_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb312_dxH1O], xmm3
	movaps [rsp + nb312_dyH1O], xmm4
	movaps [rsp + nb312_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb312_dxH2O], xmm6
	movaps [rsp + nb312_dyH2O], xmm7
	movaps [rsp + nb312_dzH2O], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

    movd mm0, eax ;# save j3 in mm0-mm3
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx
    
	;# start doing invsqrt for jO atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb312_three]
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

	movaps  xmm4, [rsp + nb312_half]
	mulps   xmm9, xmm4  ;# rinvOO 
	mulps   xmm10, xmm4 ;# rinvH1O
    mulps   xmm11, xmm4 ;# rinvH2O

	movaps  [rsp + nb312_rinvOO], xmm9
	movaps  [rsp + nb312_rinvH1O], xmm10
	movaps  [rsp + nb312_rinvH2O], xmm11
	
	;# O interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb312_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
        
    mov  rsi, [rbp + nb312_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb312_epsO], xmm0
    movaps    [rsp + nb312_epsH1], xmm3
    movaps    [rsp + nb312_epsH2], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    mulps  xmm3, [rsp + nb312_epsO]   ;# Heps
    mulps  xmm7, [rsp + nb312_epsH1]
    mulps  xmm11, [rsp + nb312_epsH2]
    mulps  xmm2, [rsp + nb312_epsO]   ;# Geps
    mulps  xmm6, [rsp + nb312_epsH1]
    mulps  xmm10, [rsp + nb312_epsH2]
    mulps  xmm3, [rsp + nb312_epsO]   ;# Heps2
    mulps  xmm7, [rsp + nb312_epsH1]
    mulps  xmm11, [rsp + nb312_epsH2]

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, [rsp + nb312_epsO]   ;# eps*Fp
    mulps  xmm5, [rsp + nb312_epsH1]
    mulps  xmm9, [rsp + nb312_epsH2]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb312_qqOO]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb312_qqOH]
    mulps  xmm9, [rsp + nb312_qqOH]
    mulps  xmm3, [rsp + nb312_qqOO]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb312_qqOH]
    mulps  xmm11, [rsp + nb312_qqOH] 
    
    ;# calculate LJ
    movaps xmm12, [rsp + nb312_rinvOO]
    mulps  xmm12, xmm12 ;# rinvsq
    movaps xmm13, xmm12 ;# rinvsq
    mulps  xmm12, xmm12 ;# rinv4
    mulps  xmm12, xmm13 ;# rinv6
    movaps xmm13, xmm12 ;# rinv6
    mulps  xmm12, xmm12 ;# rinv12
	mulps  xmm13, [rsp + nb312_c6]
	mulps  xmm12, [rsp + nb312_c12]
    movaps xmm14, xmm12
    subps  xmm14, xmm13
    
	addps  xmm14, [rsp + nb312_Vvdwtot]
	mulps  xmm13, [rsp + nb312_six]
	mulps  xmm12, [rsp + nb312_twelve]
	movaps [rsp + nb312_Vvdwtot], xmm14
    subps  xmm12, xmm13 ;# LJ fscal    
    mulps  xmm12, [rsp + nb312_rinvOO]
    movaps [rsp + nb312_fstmp], xmm12
    
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb312_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb312_vctot], xmm1
    
    movaps xmm10, [rsp + nb312_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
    
	;# move j O forces to local temp variables 
    movlps xmm11, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm12, [rdi + rcx*4] ;# jxOc jyOc  -   -

    ;# xmm11: jxOa jyOa jxOb jyOb 
    ;# xmm12: jxOc jyOc jxOd jyOd
    ;# xmm13: jzOa jzOb jzOc jzOd


    movaps xmm0, [rsp + nb312_fstmp]
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    movhps xmm11, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm12, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    mulps  xmm0, [rsp + nb312_rinvOO]
    mulps  xmm4, [rsp + nb312_rinvH1O]
    mulps  xmm8, [rsp + nb312_rinvH2O]
    
    movss  xmm13, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm14, [rdi + rcx*4 + 8] ;# jzOc  -  -  -

    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

    movhps xmm13, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm14, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

	mulps xmm0, [rsp + nb312_dxOO]
	mulps xmm1, [rsp + nb312_dyOO]
	mulps xmm2, [rsp + nb312_dzOO]
	mulps xmm3, [rsp + nb312_dxH1O]
	mulps xmm4, [rsp + nb312_dyH1O]
	mulps xmm5, [rsp + nb312_dzH1O]
	mulps xmm6, [rsp + nb312_dxH2O]
	mulps xmm7, [rsp + nb312_dyH2O]
	mulps xmm8, [rsp + nb312_dzH2O]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb312_fixO]
    addps xmm1, [rsp + nb312_fiyO]
    addps xmm2, [rsp + nb312_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb312_fixH1]
    addps xmm4, [rsp + nb312_fiyH1]
    addps xmm5, [rsp + nb312_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb312_fixH2]
    addps xmm7, [rsp + nb312_fiyH2]
    addps xmm8, [rsp + nb312_fizH2]

    movaps [rsp + nb312_fixO], xmm0
    movaps [rsp + nb312_fiyO], xmm1
    movaps [rsp + nb312_fizO], xmm2
    movaps [rsp + nb312_fixH1], xmm3
    movaps [rsp + nb312_fiyH1], xmm4
    movaps [rsp + nb312_fizH1], xmm5
    movaps [rsp + nb312_fixH2], xmm6
    movaps [rsp + nb312_fiyH2], xmm7
    movaps [rsp + nb312_fizH2], xmm8
    
    ;# xmm11 = fOx
    ;# xmm12 = fOy
    ;# xmm13 = fOz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fOzc fOzd
   
    pshufd xmm3, xmm13, 1
    pshufd xmm4, xmm13, 3
    
    movlps [rdi + rax*4], xmm11
    movhps [rdi + rbx*4], xmm11
    movlps [rdi + rcx*4], xmm12
    movhps [rdi + rdx*4], xmm12
    movss  [rdi + rax*4 + 8], xmm13
    movss  [rdi + rcx*4 + 8], xmm14
    movss  [rdi + rbx*4 + 8], xmm3
    movss  [rdi + rdx*4 + 8], xmm4
    
	;# move j H1 coordinates to local temp variables 
	mov   rsi, [rbp + nb312_pos]
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb312_ixO]
    subps xmm1, [rsp + nb312_iyO]
    subps xmm2, [rsp + nb312_izO]
    subps xmm3, [rsp + nb312_ixH1]
    subps xmm4, [rsp + nb312_iyH1]
    subps xmm5, [rsp + nb312_izH1]
    subps xmm6, [rsp + nb312_ixH2]
    subps xmm7, [rsp + nb312_iyH2]
    subps xmm8, [rsp + nb312_izH2]
    
	movaps [rsp + nb312_dxOH1], xmm0
	movaps [rsp + nb312_dyOH1], xmm1
	movaps [rsp + nb312_dzOH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb312_dxH1H1], xmm3
	movaps [rsp + nb312_dyH1H1], xmm4
	movaps [rsp + nb312_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb312_dxH2H1], xmm6
	movaps [rsp + nb312_dyH2H1], xmm7
	movaps [rsp + nb312_dzH2H1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb312_three]
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

	movaps  xmm4, [rsp + nb312_half]
	mulps   xmm9, xmm4  ;# rinvOH1
	mulps   xmm10, xmm4 ;# rinvH1H1
    mulps   xmm11, xmm4 ;# rinvH2H1

	movaps  [rsp + nb312_rinvOH1], xmm9
	movaps  [rsp + nb312_rinvH1H1], xmm10
	movaps  [rsp + nb312_rinvH2H1], xmm11
	
	;# H1 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb312_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1

    mov  rsi, [rbp + nb312_VFtab]
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
         
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb312_epsO], xmm0
    movaps    [rsp + nb312_epsH1], xmm3
    movaps    [rsp + nb312_epsH2], xmm6


    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb312_epsO]
    movaps xmm13, [rsp + nb312_epsH1]
    movaps xmm14, [rsp + nb312_epsH2]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb312_qqOH]
    movaps xmm13, [rsp + nb312_qqHH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb312_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb312_vctot], xmm1
    
    movaps xmm10, [rsp + nb312_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
	;# move j H1 forces to local temp variables 
    movlps xmm11, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm12, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -

    ;# xmm11: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm12: jxH1c jyH1c jxH1d jyH1d
    ;# xmm13: jzH1a jzH1b jzH1c jzH1d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb312_rinvOH1]
    mulps  xmm7, [rsp + nb312_rinvH1H1]
    mulps  xmm10, [rsp + nb312_rinvH2H1]
    
    movhps xmm11, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm12, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movss  xmm13, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -

    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

    movhps xmm13, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm14, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -    
    shufps xmm13, xmm14,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

	mulps xmm0, [rsp + nb312_dxOH1]
	mulps xmm1, [rsp + nb312_dyOH1]
	mulps xmm2, [rsp + nb312_dzOH1]
	mulps xmm3, [rsp + nb312_dxH1H1]
	mulps xmm4, [rsp + nb312_dyH1H1]
	mulps xmm5, [rsp + nb312_dzH1H1]
	mulps xmm6, [rsp + nb312_dxH2H1]
	mulps xmm7, [rsp + nb312_dyH2H1]
	mulps xmm8, [rsp + nb312_dzH2H1]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb312_fixO]
    addps xmm1, [rsp + nb312_fiyO]
    addps xmm2, [rsp + nb312_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb312_fixH1]
    addps xmm4, [rsp + nb312_fiyH1]
    addps xmm5, [rsp + nb312_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb312_fixH2]
    addps xmm7, [rsp + nb312_fiyH2]
    addps xmm8, [rsp + nb312_fizH2]

    movaps [rsp + nb312_fixO], xmm0
    movaps [rsp + nb312_fiyO], xmm1
    movaps [rsp + nb312_fizO], xmm2
    movaps [rsp + nb312_fixH1], xmm3
    movaps [rsp + nb312_fiyH1], xmm4
    movaps [rsp + nb312_fizH1], xmm5
    movaps [rsp + nb312_fixH2], xmm6
    movaps [rsp + nb312_fiyH2], xmm7
    movaps [rsp + nb312_fizH2], xmm8
    
    ;# xmm11 = fH1x
    ;# xmm12 = fH1y
    ;# xmm13 = fH1z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH1zc fH1zd
    
    pshufd xmm3, xmm13, 1
    pshufd xmm4, xmm13, 3
    
    movlps [rdi + rax*4 + 12], xmm11
    movhps [rdi + rbx*4 + 12], xmm11
    movlps [rdi + rcx*4 + 12], xmm12
    movhps [rdi + rdx*4 + 12], xmm12
    movss  [rdi + rax*4 + 20], xmm13
    movss  [rdi + rcx*4 + 20], xmm14
    movss  [rdi + rbx*4 + 20], xmm3
    movss  [rdi + rdx*4 + 20], xmm4
    
   	mov   rsi, [rbp + nb312_pos]
	;# move j H2 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm5, [rsi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm6, [rsi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm5 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm6 ;# jzH2c  -  jzH2d -

    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb312_ixO]
    subps xmm1, [rsp + nb312_iyO]
    subps xmm2, [rsp + nb312_izO]
    subps xmm3, [rsp + nb312_ixH1]
    subps xmm4, [rsp + nb312_iyH1]
    subps xmm5, [rsp + nb312_izH1]
    subps xmm6, [rsp + nb312_ixH2]
    subps xmm7, [rsp + nb312_iyH2]
    subps xmm8, [rsp + nb312_izH2]
    
	movaps [rsp + nb312_dxOH2], xmm0
	movaps [rsp + nb312_dyOH2], xmm1
	movaps [rsp + nb312_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb312_dxH1H2], xmm3
	movaps [rsp + nb312_dyH1H2], xmm4
	movaps [rsp + nb312_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb312_dxH2H2], xmm6
	movaps [rsp + nb312_dyH2H2], xmm7
	movaps [rsp + nb312_dzH2H2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb312_three]
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

	movaps  xmm4, [rsp + nb312_half]
	mulps   xmm9, xmm4  ;# rinvOH2
	mulps   xmm10, xmm4 ;# rinvH1H2
    mulps   xmm11, xmm4 ;# rinvH2H2

	movaps  [rsp + nb312_rinvOH2], xmm9
	movaps  [rsp + nb312_rinvH1H2], xmm10
	movaps  [rsp + nb312_rinvH2H2], xmm11
	
	;# H2 interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb312_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        
    
    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
     
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
    
    mov  rsi, [rbp + nb312_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb312_epsO], xmm0
    movaps    [rsp + nb312_epsH1], xmm3
    movaps    [rsp + nb312_epsH2], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb312_epsO]
    movaps xmm13, [rsp + nb312_epsH1]
    movaps xmm14, [rsp + nb312_epsH2]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb312_qqOH]
    movaps xmm13, [rsp + nb312_qqHH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    
    movd eax, mm0 ;# restore j3 from mm0-mm3
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb312_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb312_vctot], xmm1
    
    movaps xmm10, [rsp + nb312_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
        
	;# move j H2 forces to local temp variables 
    movlps xmm11, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm12, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -

    ;# xmm11: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm12: jxH2c jyH2c jxH2d jyH2d
    ;# xmm13: jzH2a jzH2b jzH2c jzH2d

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    movhps xmm11, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm12, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    mulps  xmm3, [rsp + nb312_rinvOH2]
    mulps  xmm7, [rsp + nb312_rinvH1H2]
    mulps  xmm10, [rsp + nb312_rinvH2H2]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movss  xmm13, [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm14, [rdi + rcx*4 + 32] ;# jzH2c  -  -  -

    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

    movss  xmm15, [rdi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm9, [rdi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm13, xmm15 ;# jzH2a  -  jzH2b  -
    movlhps xmm14, xmm9 ;# jzH2c  -  jzH2d -
    shufps xmm13, xmm14,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

	mulps xmm0, [rsp + nb312_dxOH2]
	mulps xmm1, [rsp + nb312_dyOH2]
	mulps xmm2, [rsp + nb312_dzOH2]
	mulps xmm3, [rsp + nb312_dxH1H2]
	mulps xmm4, [rsp + nb312_dyH1H2]
	mulps xmm5, [rsp + nb312_dzH1H2]
	mulps xmm6, [rsp + nb312_dxH2H2]
	mulps xmm7, [rsp + nb312_dyH2H2]
	mulps xmm8, [rsp + nb312_dzH2H2]

    movaps xmm14,  xmm0
    movaps xmm15, xmm1
    addps xmm13, xmm2
    addps xmm0, [rsp + nb312_fixO]
    addps xmm1, [rsp + nb312_fiyO]
    addps xmm2, [rsp + nb312_fizO]

    addps xmm14,  xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb312_fixH1]
    addps xmm4, [rsp + nb312_fiyH1]
    addps xmm5, [rsp + nb312_fizH1]

    addps xmm14,  xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb312_fixH2]
    addps xmm7, [rsp + nb312_fiyH2]
    addps xmm8, [rsp + nb312_fizH2]

    movaps [rsp + nb312_fixO], xmm0
    movaps [rsp + nb312_fiyO], xmm1
    movaps [rsp + nb312_fizO], xmm2
    movaps [rsp + nb312_fixH1], xmm3
    movaps [rsp + nb312_fiyH1], xmm4
    movaps [rsp + nb312_fizH1], xmm5
    movaps [rsp + nb312_fixH2], xmm6
    movaps [rsp + nb312_fiyH2], xmm7
    movaps [rsp + nb312_fizH2], xmm8
    
    ;# xmm11 = fH2x
    ;# xmm12 = fH2y
    ;# xmm13 = fH2z
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fH2zc fH2zd

    pshufd xmm3, xmm13, 1
    pshufd xmm4, xmm13, 3
    
    movlps [rdi + rax*4 + 24], xmm11
    movhps [rdi + rbx*4 + 24], xmm11
    movlps [rdi + rcx*4 + 24], xmm12
    movhps [rdi + rdx*4 + 24], xmm12
    movss  [rdi + rax*4 + 32], xmm13
    movss  [rdi + rcx*4 + 32], xmm14
    movss  [rdi + rbx*4 + 32], xmm3
    movss  [rdi + rdx*4 + 32], xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb312_innerk],  4
	jl    .nb312_single_check
	jmp   .nb312_unroll_loop
.nb312_single_check:
	add dword ptr [rsp + nb312_innerk],  4
	jnz   .nb312_single_loop
	jmp   .nb312_updateouterdata
.nb312_single_loop:
	mov   rdx, [rsp + nb312_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb312_innerjjnr],  4	

	mov rsi, [rbp + nb312_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	
	movss xmm0, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm1, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm2, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -
	movlhps xmm0, xmm6			;# xmm0 = jxO   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm1 = jyO   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm2 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + nb312_jxO], xmm0
	movaps [rsp + nb312_jyO], xmm1
	movaps [rsp + nb312_jzO], xmm2
	subps  xmm0, [rsp + nb312_ixO]
	subps  xmm1, [rsp + nb312_iyO]
	subps  xmm2, [rsp + nb312_izO]
	movaps [rsp + nb312_dxOO], xmm0
	movaps [rsp + nb312_dyOO], xmm1
	movaps [rsp + nb312_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb312_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb312_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb312_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb312_VFtab]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb312_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312_qqOO]
	movhps  xmm3, [rsp + nb312_qqOH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
	
    addps  xmm5, [rsp + nb312_vctot]
    movaps [rsp + nb312_vctot], xmm5

	mulps  xmm3, [rsp + nb312_tsc]
	
	;# start doing lj 
	xorps  xmm2, xmm2
	movss  xmm2, xmm0
	mulss  xmm2, xmm2
	movaps xmm1, xmm2
	mulss  xmm1, xmm2
	mulss  xmm1, xmm2	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulss  xmm1, [rsp + nb312_c6]
	mulss  xmm2, [rsp + nb312_c12]
	movaps xmm4, xmm2
	subss  xmm4, xmm1
	addps  xmm4, [rsp + nb312_Vvdwtot]
	mulss  xmm1, [rsp + nb312_six]
	mulss  xmm2, [rsp + nb312_twelve]
	movaps [rsp + nb312_Vvdwtot], xmm4
	subss  xmm2, xmm1
	mulss  xmm2, xmm0

	subps  xmm2, xmm3
	mulps  xmm0, xmm2
	
	movaps xmm1, xmm0
	movaps xmm2, xmm0			

	mulps   xmm0, [rsp + nb312_dxOO]
	mulps   xmm1, [rsp + nb312_dyOO]
	mulps   xmm2, [rsp + nb312_dzOO]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb312_fjxO], xmm3
	movaps  [rsp + nb312_fjyO], xmm4
	movaps  [rsp + nb312_fjzO], xmm5
	addps   xmm0, [rsp + nb312_fixO]
	addps   xmm1, [rsp + nb312_fiyO]
	addps   xmm2, [rsp + nb312_fizO]
	movaps  [rsp + nb312_fixO], xmm0
	movaps  [rsp + nb312_fiyO], xmm1
	movaps  [rsp + nb312_fizO], xmm2

	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb312_jxO]
    movaps  xmm1, [rsp + nb312_jyO]
    movaps  xmm2, [rsp + nb312_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps  xmm0, [rsp + nb312_ixH1]
	subps  xmm1, [rsp + nb312_iyH1]
	subps  xmm2, [rsp + nb312_izH1]
	subps  xmm3, [rsp + nb312_ixH2]
	subps  xmm4, [rsp + nb312_iyH2]
	subps  xmm5, [rsp + nb312_izH2]
    movaps [rsp + nb312_dxH1O], xmm0
	movaps [rsp + nb312_dyH1O], xmm1
	movaps [rsp + nb312_dzH1O], xmm2
	movaps [rsp + nb312_dxH2O], xmm3
	movaps [rsp + nb312_dyH2O], xmm4
	movaps [rsp + nb312_dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# start with H1, save H2 data 
	movaps [rsp + nb312_rsqH2O], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb312_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb312_rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb312_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb312_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312_qqOH]
	movhps  xmm3, [rsp + nb312_qqHH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
    addps  xmm5, [rsp + nb312_vctot]
    movaps [rsp + nb312_vctot], xmm5	

    xorps  xmm1, xmm1

    mulps xmm3, [rsp + nb312_tsc]
    mulps xmm3, xmm0
    subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	mulps   xmm0, [rsp + nb312_dxH1O]
	mulps   xmm1, [rsp + nb312_dyH1O]
	mulps   xmm2, [rsp + nb312_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb312_fjxO]
	movaps  xmm4, [rsp + nb312_fjyO]
	movaps  xmm5, [rsp + nb312_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb312_fjxO], xmm3
	movaps  [rsp + nb312_fjyO], xmm4
	movaps  [rsp + nb312_fjzO], xmm5
	addps   xmm0, [rsp + nb312_fixH1]
	addps   xmm1, [rsp + nb312_fiyH1]
	addps   xmm2, [rsp + nb312_fizH1]
	movaps  [rsp + nb312_fixH1], xmm0
	movaps  [rsp + nb312_fiyH1], xmm1
	movaps  [rsp + nb312_fizH1], xmm2
	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb312_rinvH2O]
	movaps xmm1, [rsp + nb312_rsqH2O]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb312_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    mulps  xmm7, [rsp + nb312_two]       ;# two*Heps2 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312_qqOH]
	movhps  xmm3, [rsp + nb312_qqHH]
	
    addps  xmm7, xmm6
    addps  xmm7, xmm5 ;# xmm7=FF 
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    mulps  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point xmm5 contains vcoul and xmm3 fijC 
    addps  xmm5, [rsp + nb312_vctot]
    movaps [rsp + nb312_vctot], xmm5	

    xorps  xmm1, xmm1

    mulps xmm3, [rsp + nb312_tsc]
    mulps xmm3, xmm0
    subps  xmm1, xmm3
	
	movaps  xmm0, xmm1
	movaps  xmm2, xmm1
	
	mulps   xmm0, [rsp + nb312_dxH2O]
	mulps   xmm1, [rsp + nb312_dyH2O]
	mulps   xmm2, [rsp + nb312_dzH2O]
	movaps  xmm3, [rsp + nb312_fjxO]
	movaps  xmm4, [rsp + nb312_fjyO]
	movaps  xmm5, [rsp + nb312_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb312_faction]
	movaps  [rsp + nb312_fjxO], xmm3
	movaps  [rsp + nb312_fjyO], xmm4
	movaps  [rsp + nb312_fjzO], xmm5
	addps   xmm0, [rsp + nb312_fixH2]
	addps   xmm1, [rsp + nb312_fiyH2]
	addps   xmm2, [rsp + nb312_fizH2]
	movaps  [rsp + nb312_fixH2], xmm0
	movaps  [rsp + nb312_fiyH2], xmm1
	movaps  [rsp + nb312_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb312_fjxO]
	movaps  xmm4, [rsp + nb312_fjyO]
	movaps  xmm5, [rsp + nb312_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, [rsi + rax*4 + 8]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 8], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec dword ptr [rsp + nb312_innerk]
	jz    .nb312_updateouterdata
	jmp   .nb312_single_loop
.nb312_updateouterdata:
	mov   ecx, [rsp + nb312_ii3]
	mov   rdi, [rbp + nb312_faction]
	mov   rsi, [rbp + nb312_fshift]
	mov   edx, [rsp + nb312_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb312_fixO]
	movaps xmm1, [rsp + nb312_fiyO] 
	movaps xmm2, [rsp + nb312_fizO]

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
	movaps xmm0, [rsp + nb312_fixH1]
	movaps xmm1, [rsp + nb312_fiyH1]
	movaps xmm2, [rsp + nb312_fizH1]

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
	movaps xmm0, [rsp + nb312_fixH2]
	movaps xmm1, [rsp + nb312_fiyH2]
	movaps xmm2, [rsp + nb312_fizH2]

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
	mov esi, [rsp + nb312_n]
    ;# get group index for i particle 
    mov   rdx, [rbp + nb312_gid]      	;# base of gid[]
    mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb312_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb312_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb312_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb312_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb312_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb312_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb312_n], esi
        jmp .nb312_outer
.nb312_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb312_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb312_end
        ;# non-zero, do one more workunit
        jmp   .nb312_threadloop
.nb312_end:
	mov eax, [rsp + nb312_nouter]
	mov ebx, [rsp + nb312_ninner]
	mov rcx, [rbp + nb312_outeriter]
	mov rdx, [rbp + nb312_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1616
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



	
.globl nb_kernel312nf_x86_64_sse
.globl _nb_kernel312nf_x86_64_sse
nb_kernel312nf_x86_64_sse:	
_nb_kernel312nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb312nf_fshift,         16
.equiv          nb312nf_gid,            24
.equiv          nb312nf_pos,            32
.equiv          nb312nf_faction,        40
.equiv          nb312nf_charge,         48
.equiv          nb312nf_p_facel,        56
.equiv          nb312nf_argkrf,         64
.equiv          nb312nf_argcrf,         72
.equiv          nb312nf_Vc,             80
.equiv          nb312nf_type,           88
.equiv          nb312nf_p_ntype,        96
.equiv          nb312nf_vdwparam,       104
.equiv          nb312nf_Vvdw,           112
.equiv          nb312nf_p_tabscale,     120
.equiv          nb312nf_VFtab,          128
.equiv          nb312nf_invsqrta,       136
.equiv          nb312nf_dvda,           144
.equiv          nb312nf_p_gbtabscale,   152
.equiv          nb312nf_GBtab,          160
.equiv          nb312nf_p_nthreads,     168
.equiv          nb312nf_count,          176
.equiv          nb312nf_mtx,            184
.equiv          nb312nf_outeriter,      192
.equiv          nb312nf_inneriter,      200
.equiv          nb312nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb312nf_ixO,            0
.equiv          nb312nf_iyO,            16
.equiv          nb312nf_izO,            32
.equiv          nb312nf_ixH1,           48
.equiv          nb312nf_iyH1,           64
.equiv          nb312nf_izH1,           80
.equiv          nb312nf_ixH2,           96
.equiv          nb312nf_iyH2,           112
.equiv          nb312nf_izH2,           128
.equiv          nb312nf_jxO,            144
.equiv          nb312nf_jyO,            160
.equiv          nb312nf_jzO,            176
.equiv          nb312nf_jxH1,           192
.equiv          nb312nf_jyH1,           208
.equiv          nb312nf_jzH1,           224
.equiv          nb312nf_jxH2,           240
.equiv          nb312nf_jyH2,           256
.equiv          nb312nf_jzH2,           272
.equiv          nb312nf_qqOO,           288
.equiv          nb312nf_qqOH,           304
.equiv          nb312nf_qqHH,           320
.equiv          nb312nf_tsc,            336
.equiv          nb312nf_c6,             352
.equiv          nb312nf_c12,            368
.equiv          nb312nf_vctot,          384
.equiv          nb312nf_Vvdwtot,        400
.equiv          nb312nf_half,           416
.equiv          nb312nf_three,          432
.equiv          nb312nf_rsqOO,          448
.equiv          nb312nf_rsqOH1,         464
.equiv          nb312nf_rsqOH2,         480
.equiv          nb312nf_rsqH1O,         496
.equiv          nb312nf_rsqH1H1,        512
.equiv          nb312nf_rsqH1H2,        528
.equiv          nb312nf_rsqH2O,         544
.equiv          nb312nf_rsqH2H1,        560
.equiv          nb312nf_rsqH2H2,        576
.equiv          nb312nf_rinvOO,         592
.equiv          nb312nf_rinvOH1,        608
.equiv          nb312nf_rinvOH2,        624
.equiv          nb312nf_rinvH1O,        640
.equiv          nb312nf_rinvH1H1,       656
.equiv          nb312nf_rinvH1H2,       672
.equiv          nb312nf_rinvH2O,        688
.equiv          nb312nf_rinvH2H1,       704
.equiv          nb312nf_rinvH2H2,       720
.equiv          nb312nf_is3,            736
.equiv          nb312nf_ii3,            740
.equiv          nb312nf_nri,            744
.equiv          nb312nf_iinr,           752
.equiv          nb312nf_jindex,         760
.equiv          nb312nf_jjnr,           768
.equiv          nb312nf_shift,          776
.equiv          nb312nf_shiftvec,       784
.equiv          nb312nf_facel,          792
.equiv          nb312nf_innerjjnr,      800
.equiv          nb312nf_innerk,         808
.equiv          nb312nf_n,              812
.equiv          nb312nf_nn1,            816
.equiv          nb312nf_nouter,         820
.equiv          nb312nf_ninner,         824
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
	sub rsp, 832		;# local variable stack space (n*16+8)
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
	mov [rsp + nb312nf_nouter], eax
	mov [rsp + nb312nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb312nf_nri], edi
	mov [rsp + nb312nf_iinr], rsi
	mov [rsp + nb312nf_jindex], rdx
	mov [rsp + nb312nf_jjnr], rcx
	mov [rsp + nb312nf_shift], r8
	mov [rsp + nb312nf_shiftvec], r9
	mov rsi, [rbp + nb312nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb312nf_facel], xmm0

	mov rax, [rbp + nb312nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb312nf_tsc], xmm3


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb312nf_half], eax
	movss xmm1, [rsp + nb312nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb312nf_half],  xmm1
	movaps [rsp + nb312nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb312nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb312nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb312nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb312nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb312nf_qqOO], xmm3
	movaps [rsp + nb312nf_qqOH], xmm4
	movaps [rsp + nb312nf_qqHH], xmm5
	
	xorps xmm0, xmm0
	mov   rdx, [rbp + nb312nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov   edx, ecx
	mov rdi, [rbp + nb312nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	add   edx, ecx
	mov   rax, [rbp + nb312nf_vdwparam]
	movlps xmm0, [rax + rdx*4] 
	movaps xmm1, xmm0
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 85  ;# 01010101
	movaps [rsp + nb312nf_c6], xmm0
	movaps [rsp + nb312nf_c12], xmm1

.nb312nf_threadloop:
        mov   rsi, [rbp + nb312nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb312nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb312nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb312nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb312nf_n], eax
        mov [rsp + nb312nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb312nf_outerstart
        jmp .nb312nf_end
.nb312nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb312nf_nouter]
	mov [rsp + nb312nf_nouter], ebx

.nb312nf_outer:
	mov   rax, [rsp + nb312nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb312nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb312nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb312nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb312nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb312nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb312nf_ixO], xmm3
	movaps [rsp + nb312nf_iyO], xmm4
	movaps [rsp + nb312nf_izO], xmm5

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
	movaps [rsp + nb312nf_ixH1], xmm0
	movaps [rsp + nb312nf_iyH1], xmm1
	movaps [rsp + nb312nf_izH1], xmm2
	movaps [rsp + nb312nf_ixH2], xmm3
	movaps [rsp + nb312nf_iyH2], xmm4
	movaps [rsp + nb312nf_izH2], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb312nf_vctot], xmm4
	movaps [rsp + nb312nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb312nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb312nf_pos]
	mov   rax, [rsp + nb312nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb312nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb312nf_ninner]
	mov   [rsp + nb312nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb312nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb312nf_unroll_loop
	jmp   .nb312nf_single_check
.nb312nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb312nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb312nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb312nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4]
	movlps xmm3, [rsi + rax*4 + 12]
	movlps xmm4, [rsi + rax*4 + 24]

	movlps xmm5, [rsi + rbx*4]
	movlps xmm6, [rsi + rbx*4 + 12]
	movlps xmm7, [rsi + rbx*4 + 24]

	movhps xmm2, [rsi + rcx*4]
	movhps xmm3, [rsi + rcx*4 + 12]
	movhps xmm4, [rsi + rcx*4 + 24]

	movhps xmm5, [rsi + rdx*4]
	movhps xmm6, [rsi + rdx*4 + 12]
	movhps xmm7, [rsi + rdx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d 
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd 
	movaps [rsp + nb312nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [rsp + nb312nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb312nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb312nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb312nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb312nf_jyH2], xmm5

	movss  xmm0, [rsi + rax*4 + 8]
	movss  xmm1, [rsi + rax*4 + 20]
	movss  xmm2, [rsi + rax*4 + 32]

	movss  xmm3, [rsi + rcx*4 + 8]
	movss  xmm4, [rsi + rcx*4 + 20]
	movss  xmm5, [rsi + rcx*4 + 32]

	movhps xmm0, [rsi + rbx*4 + 4]
	movhps xmm1, [rsi + rbx*4 + 16]
	movhps xmm2, [rsi + rbx*4 + 28]
	
	movhps xmm3, [rsi + rdx*4 + 4]
	movhps xmm4, [rsi + rdx*4 + 16]
	movhps xmm5, [rsi + rdx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb312nf_jzO],  xmm0
	movaps [rsp + nb312nf_jzH1],  xmm1
	movaps [rsp + nb312nf_jzH2],  xmm2

	movaps xmm0, [rsp + nb312nf_ixO]
	movaps xmm1, [rsp + nb312nf_iyO]
	movaps xmm2, [rsp + nb312nf_izO]
	movaps xmm3, [rsp + nb312nf_ixO]
	movaps xmm4, [rsp + nb312nf_iyO]
	movaps xmm5, [rsp + nb312nf_izO]
	subps  xmm0, [rsp + nb312nf_jxO]
	subps  xmm1, [rsp + nb312nf_jyO]
	subps  xmm2, [rsp + nb312nf_jzO]
	subps  xmm3, [rsp + nb312nf_jxH1]
	subps  xmm4, [rsp + nb312nf_jyH1]
	subps  xmm5, [rsp + nb312nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb312nf_rsqOO], xmm0
	movaps [rsp + nb312nf_rsqOH1], xmm3

	movaps xmm0, [rsp + nb312nf_ixO]
	movaps xmm1, [rsp + nb312nf_iyO]
	movaps xmm2, [rsp + nb312nf_izO]
	movaps xmm3, [rsp + nb312nf_ixH1]
	movaps xmm4, [rsp + nb312nf_iyH1]
	movaps xmm5, [rsp + nb312nf_izH1]
	subps  xmm0, [rsp + nb312nf_jxH2]
	subps  xmm1, [rsp + nb312nf_jyH2]
	subps  xmm2, [rsp + nb312nf_jzH2]
	subps  xmm3, [rsp + nb312nf_jxO]
	subps  xmm4, [rsp + nb312nf_jyO]
	subps  xmm5, [rsp + nb312nf_jzO]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb312nf_rsqOH2], xmm0
	movaps [rsp + nb312nf_rsqH1O], xmm3

	movaps xmm0, [rsp + nb312nf_ixH1]
	movaps xmm1, [rsp + nb312nf_iyH1]
	movaps xmm2, [rsp + nb312nf_izH1]
	movaps xmm3, [rsp + nb312nf_ixH1]
	movaps xmm4, [rsp + nb312nf_iyH1]
	movaps xmm5, [rsp + nb312nf_izH1]
	subps  xmm0, [rsp + nb312nf_jxH1]
	subps  xmm1, [rsp + nb312nf_jyH1]
	subps  xmm2, [rsp + nb312nf_jzH1]
	subps  xmm3, [rsp + nb312nf_jxH2]
	subps  xmm4, [rsp + nb312nf_jyH2]
	subps  xmm5, [rsp + nb312nf_jzH2]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb312nf_rsqH1H1], xmm0
	movaps [rsp + nb312nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb312nf_ixH2]
	movaps xmm1, [rsp + nb312nf_iyH2]
	movaps xmm2, [rsp + nb312nf_izH2]
	movaps xmm3, [rsp + nb312nf_ixH2]
	movaps xmm4, [rsp + nb312nf_iyH2]
	movaps xmm5, [rsp + nb312nf_izH2]
	subps  xmm0, [rsp + nb312nf_jxO]
	subps  xmm1, [rsp + nb312nf_jyO]
	subps  xmm2, [rsp + nb312nf_jzO]
	subps  xmm3, [rsp + nb312nf_jxH1]
	subps  xmm4, [rsp + nb312nf_jyH1]
	subps  xmm5, [rsp + nb312nf_jzH1]
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb312nf_rsqH2O], xmm0
	movaps [rsp + nb312nf_rsqH2H1], xmm4

	movaps xmm0, [rsp + nb312nf_ixH2]
	movaps xmm1, [rsp + nb312nf_iyH2]
	movaps xmm2, [rsp + nb312nf_izH2]
	subps  xmm0, [rsp + nb312nf_jxH2]
	subps  xmm1, [rsp + nb312nf_jyH2]
	subps  xmm2, [rsp + nb312nf_jzH2]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb312nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb312nf_half] ;# rinvH2H1 
	movaps  [rsp + nb312nf_rinvH2H2], xmm3
	movaps  [rsp + nb312nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb312nf_rsqOO]
	rsqrtps xmm5, [rsp + nb312nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb312nf_rsqOO]
	mulps   xmm5, [rsp + nb312nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312nf_half] 
	mulps   xmm7, [rsp + nb312nf_half]
	movaps  [rsp + nb312nf_rinvOO], xmm3
	movaps  [rsp + nb312nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [rsp + nb312nf_rsqOH2]
	rsqrtps xmm5, [rsp + nb312nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb312nf_rsqOH2]
	mulps   xmm5, [rsp + nb312nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312nf_half] 
	mulps   xmm7, [rsp + nb312nf_half]
	movaps  [rsp + nb312nf_rinvOH2], xmm3
	movaps  [rsp + nb312nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [rsp + nb312nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb312nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb312nf_rsqH1H1]
	mulps   xmm5, [rsp + nb312nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312nf_half] 
	mulps   xmm7, [rsp + nb312nf_half]
	movaps  [rsp + nb312nf_rinvH1H1], xmm3
	movaps  [rsp + nb312nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb312nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb312nf_three]
	mulps   xmm1, [rsp + nb312nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb312nf_half] 
	movaps  [rsp + nb312nf_rinvH2O], xmm3

	;# start with OO interaction 
	movaps xmm0, [rsp + nb312nf_rinvOO]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqOO] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]
	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
    movd mm0, eax
    movd mm1, ebx
    movd mm2, ecx
    movd mm3, edx

    mov  rsi, [rbp + nb312nf_VFtab]
    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqOO]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    ;# update vctot 
    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5
	
	;# start doing lj 
	movaps xmm2, xmm0
	mulps  xmm2, xmm2
	movaps xmm1, xmm2
	mulps  xmm1, xmm2
	mulps  xmm1, xmm2	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm1, [rsp + nb312nf_c6]
	mulps  xmm2, [rsp + nb312nf_c12]
	movaps xmm4, xmm2
	subps  xmm4, xmm1
	addps  xmm4, [rsp + nb312nf_Vvdwtot]
	movaps [rsp + nb312nf_Vvdwtot], xmm4

	;# O-H1 interaction 
	movaps xmm0, [rsp + nb312nf_rinvOH1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqOH1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2
	
    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5
	
	;# O-H2 interaction  
	movaps xmm0, [rsp + nb312nf_rinvOH2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqOH2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5

	;# H1-O interaction 
	movaps xmm0, [rsp + nb312nf_rinvH1O]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH1O] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5

	;# H1-H1 interaction 
	movaps xmm0, [rsp + nb312nf_rinvH1H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH1H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5

	;# H1-H2 interaction 
	movaps xmm0, [rsp + nb312nf_rinvH1H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH1H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5

	;# H2-O interaction 
	movaps xmm0, [rsp + nb312nf_rinvH2O]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH2O] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqOH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5

	;# H2-H1 interaction 
	movaps xmm0, [rsp + nb312nf_rinvH2H1]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH2H1] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5
	
	;# H2-H2 interaction 
	movaps xmm0, [rsp + nb312nf_rinvH2H2]
	movaps xmm1, xmm0
	mulps  xmm1, [rsp + nb312nf_rsqH2H2] ;# xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]	
	movhlps xmm2, xmm1
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 

	pslld   mm6, 2
	pslld   mm7, 2

    movd eax, mm6
    psrlq mm6, 32
    movd ecx, mm7
    psrlq mm7, 32
    movd ebx, mm6
    movd edx, mm7

    movlps xmm5, [rsi + rax*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm5, [rsi + rbx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rax*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm7, [rsi + rbx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  

    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 
    movaps xmm3, [rsp + nb312nf_qqHH]
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul  

    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5	
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb312nf_innerk],  4
	jl    .nb312nf_single_check
	jmp   .nb312nf_unroll_loop
.nb312nf_single_check:
	add dword ptr [rsp + nb312nf_innerk],  4
	jnz   .nb312nf_single_loop
	jmp   .nb312nf_updateouterdata
.nb312nf_single_loop:
	mov   rdx, [rsp + nb312nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb312nf_innerjjnr],  4	

	mov rsi, [rbp + nb312nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm4, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb312nf_ixO]     
	movaps  xmm1, [rsp + nb312nf_iyO]
	movaps  xmm2, [rsp + nb312nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + nb312nf_jxO], xmm3
	movaps [rsp + nb312nf_jyO], xmm4
	movaps [rsp + nb312nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb312nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb312nf_half] ;# rinv iO - j water 

	movaps  xmm1, xmm3
	mulps   xmm1, xmm0	;# xmm1=r 
	movaps  xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb312nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2
	
    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

	mov rsi, [rbp + nb312nf_VFtab]
	
    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312nf_qqOO]
	movhps  xmm3, [rsp + nb312nf_qqOH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point xmm5 contains vcoul 
	
    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5
	
	;# start doing lj 
	xorps  xmm2, xmm2
	movss  xmm2, xmm0
	mulss  xmm2, xmm2
	movaps xmm1, xmm2
	mulss  xmm1, xmm2
	mulss  xmm1, xmm2	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulss  xmm1, [rsp + nb312nf_c6]
	mulss  xmm2, [rsp + nb312nf_c12]
	movaps xmm4, xmm2
	subss  xmm4, xmm1
	addps  xmm4, [rsp + nb312nf_Vvdwtot]
	movaps [rsp + nb312nf_Vvdwtot], xmm4
	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
	movaps  xmm0, [rsp + nb312nf_ixH1]
	movaps  xmm1, [rsp + nb312nf_iyH1]
	movaps  xmm2, [rsp + nb312nf_izH1]	
	movaps  xmm3, [rsp + nb312nf_ixH2] 
	movaps  xmm4, [rsp + nb312nf_iyH2] 
	movaps  xmm5, [rsp + nb312nf_izH2] 
	subps   xmm0, [rsp + nb312nf_jxO]
	subps   xmm1, [rsp + nb312nf_jyO]
	subps   xmm2, [rsp + nb312nf_jzO]
	subps   xmm3, [rsp + nb312nf_jxO]
	subps   xmm4, [rsp + nb312nf_jyO]
	subps   xmm5, [rsp + nb312nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# start with H1, save H2 data 
	movaps [rsp + nb312nf_rsqH2O], xmm4
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb312nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb312nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb312nf_half] ;# rinv H2 - j water  

	;# start with H1, save H2 data 
	movaps [rsp + nb312nf_rinvH2O], xmm7

	movaps xmm1, xmm3
	mulps  xmm1, xmm0	;# xmm1=r 
	movaps xmm0, xmm3	;# xmm0=rinv 
	mulps  xmm1, [rsp + nb312nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312nf_qqOH]
	movhps  xmm3, [rsp + nb312nf_qqHH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point xmm5 contains vcoul 
    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5	

	;# do table for H2 - j water interaction 
	movaps xmm0, [rsp + nb312nf_rinvH2O]
	movaps xmm1, [rsp + nb312nf_rsqH2O]
	mulps  xmm1, xmm0	;# xmm0=rinv, xmm1=r 
	mulps  xmm1, [rsp + nb312nf_tsc]
	
	movhlps xmm2, xmm1	
    cvttps2pi mm6, xmm1
    cvttps2pi mm7, xmm2     ;# mm6/mm7 contain lu indices 
    cvtpi2ps xmm3, mm6
    cvtpi2ps xmm2, mm7
	movlhps  xmm3, xmm2
	subps    xmm1, xmm3	;# xmm1=eps 
    movaps xmm2, xmm1
    mulps  xmm2, xmm2       ;# xmm2=eps2 
	pslld   mm6, 2
	pslld   mm7, 2

    movd ebx, mm6
    movd ecx, mm7
    psrlq mm7, 32
    movd edx, mm7		;# table indices in ebx,ecx,edx 

    movlps xmm5, [rsi + rbx*4]
    movlps xmm7, [rsi + rcx*4]
    movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 
    movaps xmm4, xmm5
    shufps xmm4, xmm7, 136  ;# 10001000
    shufps xmm5, xmm7, 221  ;# 11011101

    movlps xmm7, [rsi + rbx*4 + 8]
    movlps xmm3, [rsi + rcx*4 + 8]
    movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    movaps xmm6, xmm7
    shufps xmm6, xmm3, 136  ;# 10001000
    shufps xmm7, xmm3, 221  ;# 11011101
    ;# coulomb table ready, in xmm4-xmm7  
    mulps  xmm6, xmm1       ;# xmm6=Geps 
    mulps  xmm7, xmm2       ;# xmm7=Heps2 
    addps  xmm5, xmm6
    addps  xmm5, xmm7       ;# xmm5=Fp 

	xorps  xmm3, xmm3
	;# fetch charges to xmm3 (temporary) 
	movss   xmm3, [rsp + nb312nf_qqOH]
	movhps  xmm3, [rsp + nb312nf_qqHH]
	
    mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    addps  xmm5, xmm4 ;# xmm5=VV 
    mulps  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point xmm5 contains vcoul 
    addps  xmm5, [rsp + nb312nf_vctot]
    movaps [rsp + nb312nf_vctot], xmm5	

	dec dword ptr [rsp + nb312nf_innerk]
	jz    .nb312nf_updateouterdata
	jmp   .nb312nf_single_loop
.nb312nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb312nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb312nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb312nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb312nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb312nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb312nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb312nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb312nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb312nf_n], esi
        jmp .nb312nf_outer
.nb312nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb312nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb312nf_end
        ;# non-zero, do one more workunit
        jmp   .nb312nf_threadloop
.nb312nf_end:
	mov eax, [rsp + nb312nf_nouter]
	mov ebx, [rsp + nb312nf_ninner]
	mov rcx, [rbp + nb312nf_outeriter]
	mov rdx, [rbp + nb312nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 832
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
	
