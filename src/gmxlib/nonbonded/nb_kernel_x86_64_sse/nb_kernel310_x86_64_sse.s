##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce 
##






.globl nb_kernel310_x86_64_sse
.globl _nb_kernel310_x86_64_sse
nb_kernel310_x86_64_sse:        
_nb_kernel310_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb310_fshift, 16
.set nb310_gid, 24
.set nb310_pos, 32
.set nb310_faction, 40
.set nb310_charge, 48
.set nb310_p_facel, 56
.set nb310_argkrf, 64
.set nb310_argcrf, 72
.set nb310_Vc, 80
.set nb310_type, 88
.set nb310_p_ntype, 96
.set nb310_vdwparam, 104
.set nb310_Vvdw, 112
.set nb310_p_tabscale, 120
.set nb310_VFtab, 128
.set nb310_invsqrta, 136
.set nb310_dvda, 144
.set nb310_p_gbtabscale, 152
.set nb310_GBtab, 160
.set nb310_p_nthreads, 168
.set nb310_count, 176
.set nb310_mtx, 184
.set nb310_outeriter, 192
.set nb310_inneriter, 200
.set nb310_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb310_ix, 0
.set nb310_iy, 16
.set nb310_iz, 32
.set nb310_iq, 48
.set nb310_dx, 64
.set nb310_dy, 80
.set nb310_dz, 96
.set nb310_two, 112
.set nb310_six, 128
.set nb310_twelve, 144
.set nb310_tsc, 160
.set nb310_qq, 176
.set nb310_c6, 192
.set nb310_c12, 208
.set nb310_fscal, 224
.set nb310_vctot, 240
.set nb310_Vvdwtot, 256
.set nb310_fix, 272
.set nb310_fiy, 288
.set nb310_fiz, 304
.set nb310_half, 320
.set nb310_three, 336
.set nb310_nri, 352
.set nb310_iinr, 360
.set nb310_jindex, 368
.set nb310_jjnr, 376
.set nb310_shift, 384
.set nb310_shiftvec, 392
.set nb310_facel, 400
.set nb310_innerjjnr, 408
.set nb310_is3, 416
.set nb310_ii3, 420
.set nb310_ntia, 424
.set nb310_innerk, 428
.set nb310_n, 432
.set nb310_nn1, 436
.set nb310_ntype, 440
.set nb310_nouter, 444
.set nb310_ninner, 448


        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $472,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb310_nouter(%rsp)
        movl %eax,nb310_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb310_nri(%rsp)
        movq %rsi,nb310_iinr(%rsp)
        movq %rdx,nb310_jindex(%rsp)
        movq %rcx,nb310_jjnr(%rsp)
        movq %r8,nb310_shift(%rsp)
        movq %r9,nb310_shiftvec(%rsp)
        movq nb310_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb310_ntype(%rsp)
        movq nb310_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb310_facel(%rsp)


        movq nb310_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb310_tsc(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb310_half(%rsp)
        movss nb310_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb310_half(%rsp)
        movaps %xmm2,nb310_two(%rsp)
        movaps %xmm3,nb310_three(%rsp)
        movaps %xmm4,nb310_six(%rsp)
        movaps %xmm5,nb310_twelve(%rsp)

_nb_kernel310_x86_64_sse.nb310_threadloop: 
        movq  nb310_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel310_x86_64_sse.nb310_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310_x86_64_sse.nb310_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310_n(%rsp)
        movl %ebx,nb310_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310_x86_64_sse.nb310_outerstart
        jmp _nb_kernel310_x86_64_sse.nb310_end

_nb_kernel310_x86_64_sse.nb310_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310_nouter(%rsp),%ebx
        movl %ebx,nb310_nouter(%rsp)

_nb_kernel310_x86_64_sse.nb310_outer: 
        movq  nb310_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb310_is3(%rsp)      ## store is3 

        movq  nb310_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb310_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb310_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb310_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb310_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb310_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb310_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb310_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb310_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb310_ix(%rsp)
        movaps %xmm1,nb310_iy(%rsp)
        movaps %xmm2,nb310_iz(%rsp)

        movl  %ebx,nb310_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm15,%xmm15
        movaps %xmm15,nb310_vctot(%rsp)
        movaps %xmm15,nb310_Vvdwtot(%rsp)
        movaps %xmm15,%xmm14
        movaps %xmm15,%xmm13

        movq  nb310_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb310_pos(%rbp),%rsi
        movq  nb310_faction(%rbp),%rdi
        movq  nb310_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb310_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb310_ninner(%rsp),%ecx
        movl  %ecx,nb310_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb310_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel310_x86_64_sse.nb310_unroll_loop
        jmp   _nb_kernel310_x86_64_sse.nb310_finish_inner
_nb_kernel310_x86_64_sse.nb310_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb310_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 

        addq $16,nb310_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 


        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        lea  (%r10,%r10,2),%rcx     ## replace jnr with j3 
        lea  (%r11,%r11,2),%rdx

        ## load coordinates
        movq nb310_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm1      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm2      ## x2 y2 - - 
        movlps (%rdi,%rcx,4),%xmm3      ## x3 y3 - -
        movlps (%rdi,%rdx,4),%xmm4      ## x4 y4 - -

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm6      ## z2 - - - 
        movss 8(%rdi,%rcx,4),%xmm7      ## z3 - - - 
        movss 8(%rdi,%rdx,4),%xmm8      ## z4 - - - 

    unpcklps %xmm3,%xmm1 ## x1 x3 y1 y3
    unpcklps %xmm4,%xmm2 ## x2 x4 y2 y4
    unpcklps %xmm7,%xmm5 ## z1 z3 -  -
    unpcklps %xmm8,%xmm6 ## z2 z4 -  -
        movq nb310_charge(%rbp),%rsi

    movaps %xmm1,%xmm3

    unpcklps %xmm2,%xmm1 ## x1 x2 x3 x4
    unpckhps %xmm2,%xmm3 ## y1 y2 y3 y4
    unpcklps %xmm6,%xmm5 ## z1 z2 z3 z4

        ## calc dr  
        subps nb310_ix(%rsp),%xmm1
        subps nb310_iy(%rsp),%xmm3
        subps nb310_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm3,%xmm10
    movaps %xmm5,%xmm11

        movss (%rsi,%r8,4),%xmm0
        movss (%rsi,%r10,4),%xmm2
        movss (%rsi,%r9,4),%xmm6
        movss (%rsi,%r11,4),%xmm8

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm3,%xmm3
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm3
        addps %xmm5,%xmm3
        ## rsq in xmm3
        movq nb310_type(%rbp),%rsi

    unpcklps %xmm2,%xmm0
    unpcklps %xmm8,%xmm6

    unpcklps %xmm6,%xmm0


    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm3,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb310_three(%rsp),%xmm1
        mulps %xmm3,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps nb310_half(%rsp),%xmm1
    ## xmm1=rinv
    ## xmm3=rsq
        mulps nb310_iq(%rsp),%xmm0

    ## vdw types
        movl (%rsi,%r8,4),%r8d
        movl (%rsi,%r9,4),%r9d
        movl (%rsi,%r10,4),%r10d
        movl (%rsi,%r11,4),%r11d

    mulps %xmm1,%xmm3 ## r
    mulps nb310_tsc(%rsp),%xmm3   ## rtab
    movaps %xmm0,nb310_qq(%rsp)

    ## truncate and convert to integers
    cvttps2dq %xmm3,%xmm2

        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d

    ## convert back to float
    cvtdq2ps  %xmm2,%xmm0

    movl nb310_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

    ## multiply by 4
    pslld   $2,%xmm2

    ## move to integer registers
    movhlps %xmm2,%xmm7
    movd    %xmm2,%r12d
    movd    %xmm7,%r14d
    pshufd $1,%xmm2,%xmm2
    pshufd $1,%xmm7,%xmm7
    movd    %xmm2,%r13d
    movd    %xmm7,%r15d


    ## calculate eps
    subps     %xmm0,%xmm3

        movq nb310_vdwparam(%rbp),%rsi
        movlps (%rsi,%r8,4),%xmm7
        movlps (%rsi,%r10,4),%xmm8
        movhps (%rsi,%r9,4),%xmm7
        movhps (%rsi,%r11,4),%xmm8

        movaps %xmm7,%xmm12
        shufps $136,%xmm8,%xmm12 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101

    movaps %xmm12,nb310_c6(%rsp)
    movaps %xmm7,nb310_c12(%rsp)

        movq nb310_VFtab(%rbp),%rsi
    ## load table data
        movlps (%rsi,%r12,4),%xmm5
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm5
        movhps (%rsi,%r15,4),%xmm7

    movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

    movaps %xmm1,%xmm0 ## rinv
    mulps  %xmm0,%xmm0 ## rinvsq
    movaps %xmm0,%xmm2 ## rinvsq
    mulps  %xmm2,%xmm2 ## rinv4
    mulps  %xmm0,%xmm2 ## rinv6
    movaps %xmm2,%xmm12
    mulps  %xmm12,%xmm12 ## rinv12

        movlps 8(%rsi,%r12,4),%xmm7
        movlps 8(%rsi,%r14,4),%xmm8
        movhps 8(%rsi,%r13,4),%xmm7
        movhps 8(%rsi,%r15,4),%xmm8

    movaps %xmm7,%xmm6

    mulps  nb310_c6(%rsp),%xmm2      ## vvdw6=c6*rinv6
        mulps  nb310_c12(%rsp),%xmm12     ## vvdw12=c12*rinv12     

        movaps %xmm12,%xmm0
        subps  %xmm2,%xmm12     ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addps  nb310_Vvdwtot(%rsp),%xmm12
    movaps %xmm12,nb310_Vvdwtot(%rsp)

        shufps $136,%xmm8,%xmm6 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101
    ## table data ready in xmm4-xmm7

    mulps %xmm3,%xmm7  ## Heps
    mulps  %xmm3,%xmm6 ## Geps
    mulps %xmm3,%xmm7  ## Heps2

    addps  %xmm6,%xmm5  ## F+Geps
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm7  ## 2*Heps2
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulps  %xmm3,%xmm5  ## eps*Fp
    addps  %xmm4,%xmm5  ## VV
    mulps  nb310_qq(%rsp),%xmm5     ## VV*qq=vcoul
    mulps  nb310_qq(%rsp),%xmm7     ## FF*qq=fijC

    ## LJ forces
    mulps  nb310_six(%rsp),%xmm2
    mulps  nb310_twelve(%rsp),%xmm0
    subps  %xmm2,%xmm0
    mulps  %xmm1,%xmm0 ## (12*vnb12-6*vnb6)*rinv

    ## add potential to vctot 
        addps  nb310_vctot(%rsp),%xmm5
    movaps %xmm5,nb310_vctot(%rsp)

    mulps  nb310_tsc(%rsp),%xmm7
    subps  %xmm7,%xmm0

    mulps  %xmm1,%xmm0 ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm0,%xmm9
        mulps  %xmm0,%xmm10
        mulps  %xmm0,%xmm11

        movq nb310_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%rcx,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%rdx,4),%xmm1 ## x3 y3 x4 y4

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%rax,4)
        movlps %xmm1,(%rsi,%rcx,4)
        movhps %xmm0,(%rsi,%rbx,4)
        movhps %xmm1,(%rsi,%rdx,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm10 ## fjz2 - - -
    movhlps %xmm11,%xmm9     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm8  ## fjz4 - - -

        addss  8(%rsi,%rax,4),%xmm11
        addss  8(%rsi,%rbx,4),%xmm10
        addss  8(%rsi,%rcx,4),%xmm9
        addss  8(%rsi,%rdx,4),%xmm8
        movss  %xmm11,8(%rsi,%rax,4)
        movss  %xmm10,8(%rsi,%rbx,4)
        movss  %xmm9,8(%rsi,%rcx,4)
        movss  %xmm8,8(%rsi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb310_innerk(%rsp)
        jl    _nb_kernel310_x86_64_sse.nb310_finish_inner
        jmp   _nb_kernel310_x86_64_sse.nb310_unroll_loop
_nb_kernel310_x86_64_sse.nb310_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb310_innerk(%rsp)
    movl  nb310_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel310_x86_64_sse.nb310_dopair
    jmp   _nb_kernel310_x86_64_sse.nb310_checksingle
_nb_kernel310_x86_64_sse.nb310_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb310_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb310_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb310_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rbx,4),%xmm2

    unpcklps %xmm2,%xmm0 ## jqa jqb 
        mulps nb310_iq(%rsp),%xmm0
    movaps %xmm0,nb310_qq(%rsp)

        movq nb310_type(%rbp),%rsi
    ## vdw parameters
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb310_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb310_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb310_c6(%rsp)
    movaps %xmm3,nb310_c12(%rsp)

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load coordinates
        movq nb310_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm4      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm5      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm6      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr  
        subps nb310_ix(%rsp),%xmm4
        subps nb310_iy(%rsp),%xmm5
        subps nb310_iz(%rsp),%xmm6

        ## store dr in xmm9-xmm11
    movaps %xmm4,%xmm9
    movaps %xmm5,%xmm10
    movaps %xmm6,%xmm11

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm4,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb310_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps nb310_half(%rsp),%xmm1
    ## xmm1=rinv
    movaps %xmm4,%xmm3
    ## xmm3=rsq 

    mulps %xmm1,%xmm3 ## r
    mulps nb310_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm3,%xmm2

    ## convert back to float
    cvtdq2ps  %xmm2,%xmm0

    ## multiply by 4
    pslld   $2,%xmm2

    ## move to integer registers
    movd    %xmm2,%r12d
    pshufd $1,%xmm2,%xmm2
    movd    %xmm2,%r13d

    ## calculate eps
    subps     %xmm0,%xmm3

        movq nb310_VFtab(%rbp),%rsi
    ## load table data
        movlps (%rsi,%r12,4),%xmm4
        movlps (%rsi,%r13,4),%xmm5
    unpcklps %xmm5,%xmm4
    movhlps %xmm4,%xmm5

    movaps %xmm1,%xmm0 ## rinv
    mulps  %xmm0,%xmm0 ## rinvsq
    movaps %xmm0,%xmm2 ## rinvsq
    mulps  %xmm2,%xmm2 ## rinv4
    mulps  %xmm0,%xmm2 ## rinv6
    movaps %xmm2,%xmm12
    mulps  %xmm12,%xmm12 ## rinv12

        movlps 8(%rsi,%r12,4),%xmm6
        movlps 8(%rsi,%r13,4),%xmm7
    unpcklps %xmm7,%xmm6
    movhlps %xmm6,%xmm7
    ## table data ready in xmm4-xmm7

    mulps  nb310_c6(%rsp),%xmm2      ## vvdw6=c6*rinv6
        mulps  nb310_c12(%rsp),%xmm12     ## vvdw12=c12*rinv12     

        movaps %xmm12,%xmm0
        subps  %xmm2,%xmm12     ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addps  nb310_Vvdwtot(%rsp),%xmm12
    movlps %xmm12,nb310_Vvdwtot(%rsp)

    mulps %xmm3,%xmm7  ## Heps
    mulps  %xmm3,%xmm6 ## Geps
    mulps %xmm3,%xmm7  ## Heps2

    addps  %xmm6,%xmm5  ## F+Geps
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm7  ## 2*Heps2
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulps  %xmm3,%xmm5  ## eps*Fp
    addps  %xmm4,%xmm5  ## VV
    mulps  nb310_qq(%rsp),%xmm5     ## VV*qq=vcoul
    mulps  nb310_qq(%rsp),%xmm7     ## FF*qq=fijC

    ## LJ forces
    mulps  nb310_six(%rsp),%xmm2
    mulps  nb310_twelve(%rsp),%xmm0
    subps  %xmm2,%xmm0
    mulps  %xmm1,%xmm0 ## (12*vnb12-6*vnb6)*rinv

    ## add potential to vctot 
        addps  nb310_vctot(%rsp),%xmm5
    movlps %xmm5,nb310_vctot(%rsp)

    xorps %xmm8,%xmm8

    mulps  nb310_tsc(%rsp),%xmm7
    subps  %xmm7,%xmm0

    mulps  %xmm1,%xmm0 ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm0,%xmm9
        mulps  %xmm0,%xmm10
        mulps  %xmm0,%xmm11

    movlhps %xmm8,%xmm9
    movlhps %xmm8,%xmm10
    movlhps %xmm8,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb310_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2

    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    addps    %xmm9,%xmm0

        movlps %xmm0,(%rsi,%rax,4)
        movhps %xmm0,(%rsi,%rbx,4)

    ## z forces
    pshufd $1,%xmm11,%xmm8
    addss  8(%rsi,%rax,4),%xmm11
    addss  8(%rsi,%rbx,4),%xmm8
    movss  %xmm11,8(%rsi,%rax,4)
    movss  %xmm8,8(%rsi,%rbx,4)

_nb_kernel310_x86_64_sse.nb310_checksingle:     
    movl  nb310_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel310_x86_64_sse.nb310_dosingle
    jmp    _nb_kernel310_x86_64_sse.nb310_updateouterdata

_nb_kernel310_x86_64_sse.nb310_dosingle: 
    movq nb310_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb310_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0

        mulss nb310_iq(%rsp),%xmm0
    movaps %xmm0,nb310_qq(%rsp)

        movq nb310_type(%rbp),%rsi
    ## vdw parameters
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb310_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb310_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3

    movaps %xmm0,nb310_c6(%rsp)
    movaps %xmm3,nb310_c12(%rsp)

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        movq nb310_pos(%rbp),%rdi
        movss (%rdi,%rax,4),%xmm4           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm5       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm6       ## 13 - - - 

        ## calc dr  
        subss nb310_ix(%rsp),%xmm4
        subss nb310_iy(%rsp),%xmm5
        subss nb310_iz(%rsp),%xmm6

        ## store dr in xmm9-xmm11
    movaps %xmm4,%xmm9
    movaps %xmm5,%xmm10
    movaps %xmm6,%xmm11

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

    ## calculate rinv=1/sqrt(rsq)
        rsqrtss %xmm4,%xmm5
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb310_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss nb310_half(%rsp),%xmm1
    ## xmm1=rinv
    movaps %xmm4,%xmm3
    ## xmm3=rsq 

    mulss %xmm1,%xmm3 ## r
    mulss nb310_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttss2si %xmm3,%r12d

    ## convert back to float
    cvtsi2ss  %r12d,%xmm0

    ## multiply by 4
    shll      $2,%r12d

    ## calculate eps
    subss     %xmm0,%xmm3

        movq nb310_VFtab(%rbp),%rsi

    movaps %xmm1,%xmm0 ## rinv
    mulss  %xmm0,%xmm0 ## rinvsq
    movaps %xmm0,%xmm2 ## rinvsq
    mulss  %xmm2,%xmm2 ## rinv4
    mulss  %xmm0,%xmm2 ## rinv6
    movaps %xmm2,%xmm12
    mulss  %xmm12,%xmm12 ## rinv12

    ## load table data
        movss (%rsi,%r12,4),%xmm4
        movss 4(%rsi,%r12,4),%xmm5
        movss 8(%rsi,%r12,4),%xmm6
        movss 12(%rsi,%r12,4),%xmm7
    ## table data ready in xmm4-xmm7

    mulss  nb310_c6(%rsp),%xmm2      ## vvdw6=c6*rinv6
        mulss  nb310_c12(%rsp),%xmm12     ## vvdw12=c12*rinv12     

        movaps %xmm12,%xmm0
        subss  %xmm2,%xmm12     ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addss  nb310_Vvdwtot(%rsp),%xmm12
    movss %xmm12,nb310_Vvdwtot(%rsp)

    mulss %xmm3,%xmm7  ## Heps
    mulss  %xmm3,%xmm6 ## Geps
    mulss %xmm3,%xmm7  ## Heps2

    addss  %xmm6,%xmm5  ## F+Geps
    addss  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addss  %xmm7,%xmm7  ## 2*Heps2
    addss  %xmm6,%xmm7  ## 2*Heps2+Geps
    addss  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulss  %xmm3,%xmm5  ## eps*Fp
    addss  %xmm4,%xmm5  ## VV
    mulss  nb310_qq(%rsp),%xmm5     ## VV*qq=vcoul
    mulss  nb310_qq(%rsp),%xmm7     ## FF*qq=fijC

    ## LJ forces
    mulss  nb310_six(%rsp),%xmm2
    mulss  nb310_twelve(%rsp),%xmm0
    subss  %xmm2,%xmm0
    mulss  %xmm1,%xmm0 ## (12*vnb12-6*vnb6)*rinv

    ## add potential to vctot 
        addss  nb310_vctot(%rsp),%xmm5
    movss %xmm5,nb310_vctot(%rsp)

    mulss  nb310_tsc(%rsp),%xmm7
    subss  %xmm7,%xmm0

    mulss  %xmm1,%xmm0 ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm0,%xmm9
        mulss  %xmm0,%xmm10
        mulss  %xmm0,%xmm11

        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb310_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel310_x86_64_sse.nb310_updateouterdata: 
        movl  nb310_ii3(%rsp),%ecx
        movq  nb310_faction(%rbp),%rdi
        movq  nb310_fshift(%rbp),%rsi
        movl  nb310_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm0
        movhlps %xmm14,%xmm1
        movhlps %xmm15,%xmm2
        addps  %xmm13,%xmm0
        addps  %xmm14,%xmm1
        addps  %xmm15,%xmm2
    movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## increment fshift force  
        movss  (%rsi,%rdx,4),%xmm3
        movss  4(%rsi,%rdx,4),%xmm4
        movss  8(%rsi,%rdx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,4(%rsi,%rdx,4)
        movss  %xmm5,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb310_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb310_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb310_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb310_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb310_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb310_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb310_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310_x86_64_sse.nb310_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310_n(%rsp)
        jmp _nb_kernel310_x86_64_sse.nb310_outer
_nb_kernel310_x86_64_sse.nb310_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310_x86_64_sse.nb310_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310_x86_64_sse.nb310_threadloop
_nb_kernel310_x86_64_sse.nb310_end: 

        movl nb310_nouter(%rsp),%eax
        movl nb310_ninner(%rsp),%ebx
        movq nb310_outeriter(%rbp),%rcx
        movq nb310_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $472,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel310nf_x86_64_sse
.globl _nb_kernel310nf_x86_64_sse
nb_kernel310nf_x86_64_sse:      
_nb_kernel310nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb310nf_fshift, 16
.set nb310nf_gid, 24
.set nb310nf_pos, 32
.set nb310nf_faction, 40
.set nb310nf_charge, 48
.set nb310nf_p_facel, 56
.set nb310nf_argkrf, 64
.set nb310nf_argcrf, 72
.set nb310nf_Vc, 80
.set nb310nf_type, 88
.set nb310nf_p_ntype, 96
.set nb310nf_vdwparam, 104
.set nb310nf_Vvdw, 112
.set nb310nf_p_tabscale, 120
.set nb310nf_VFtab, 128
.set nb310nf_invsqrta, 136
.set nb310nf_dvda, 144
.set nb310nf_p_gbtabscale, 152
.set nb310nf_GBtab, 160
.set nb310nf_p_nthreads, 168
.set nb310nf_count, 176
.set nb310nf_mtx, 184
.set nb310nf_outeriter, 192
.set nb310nf_inneriter, 200
.set nb310nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb310nf_ix, 0
.set nb310nf_iy, 16
.set nb310nf_iz, 32
.set nb310nf_iq, 48
.set nb310nf_tsc, 64
.set nb310nf_qq, 80
.set nb310nf_c6, 96
.set nb310nf_c12, 112
.set nb310nf_vctot, 128
.set nb310nf_Vvdwtot, 144
.set nb310nf_half, 160
.set nb310nf_three, 176
.set nb310nf_nri, 192
.set nb310nf_iinr, 200
.set nb310nf_jindex, 208
.set nb310nf_jjnr, 216
.set nb310nf_shift, 224
.set nb310nf_shiftvec, 232
.set nb310nf_facel, 240
.set nb310nf_innerjjnr, 248
.set nb310nf_is3, 256
.set nb310nf_ii3, 260
.set nb310nf_ntia, 264
.set nb310nf_innerk, 268
.set nb310nf_n, 272
.set nb310nf_nn1, 276
.set nb310nf_ntype, 280
.set nb310nf_nouter, 284
.set nb310nf_ninner, 288

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $312,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb310nf_nouter(%rsp)
        movl %eax,nb310nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb310nf_nri(%rsp)
        movq %rsi,nb310nf_iinr(%rsp)
        movq %rdx,nb310nf_jindex(%rsp)
        movq %rcx,nb310nf_jjnr(%rsp)
        movq %r8,nb310nf_shift(%rsp)
        movq %r9,nb310nf_shiftvec(%rsp)
        movq nb310nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb310nf_ntype(%rsp)
        movq nb310nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb310nf_facel(%rsp)

        movq nb310nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb310nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb310nf_half(%rsp)
        movss nb310nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb310nf_half(%rsp)
        movaps %xmm3,nb310nf_three(%rsp)

_nb_kernel310nf_x86_64_sse.nb310nf_threadloop: 
        movq  nb310nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel310nf_x86_64_sse.nb310nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310nf_x86_64_sse.nb310nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310nf_n(%rsp)
        movl %ebx,nb310nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310nf_x86_64_sse.nb310nf_outerstart
        jmp _nb_kernel310nf_x86_64_sse.nb310nf_end

_nb_kernel310nf_x86_64_sse.nb310nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310nf_nouter(%rsp),%ebx
        movl %ebx,nb310nf_nouter(%rsp)

_nb_kernel310nf_x86_64_sse.nb310nf_outer: 
        movq  nb310nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb310nf_is3(%rsp)            ## store is3 

        movq  nb310nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb310nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb310nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb310nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb310nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb310nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb310nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb310nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb310nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb310nf_ix(%rsp)
        movaps %xmm1,nb310nf_iy(%rsp)
        movaps %xmm2,nb310nf_iz(%rsp)

        movl  %ebx,nb310nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb310nf_vctot(%rsp)
        movaps %xmm4,nb310nf_Vvdwtot(%rsp)

        movq  nb310nf_jindex(%rsp),%rax
        movq  (%rax,%rsi,4),%rcx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb310nf_pos(%rbp),%rsi
        movq  nb310nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb310nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb310nf_ninner(%rsp),%ecx
        movl  %ecx,nb310nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb310nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel310nf_x86_64_sse.nb310nf_unroll_loop
        jmp   _nb_kernel310nf_x86_64_sse.nb310nf_finish_inner
_nb_kernel310nf_x86_64_sse.nb310nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb310nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb310nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb310nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb310nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        mulps  %xmm2,%xmm3
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm3,nb310nf_qq(%rsp)

        movq nb310nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb310nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb310nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb310nf_c6(%rsp)
        movaps %xmm6,nb310nf_c12(%rsp)

        movq nb310nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   

        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ix-iz to xmm4-xmm6 
        movaps nb310nf_ix(%rsp),%xmm4
        movaps nb310nf_iy(%rsp),%xmm5
        movaps nb310nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb310nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%rsp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb310nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb310nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310nf_vctot(%rsp),%xmm5
        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6
        movaps %xmm5,nb310nf_vctot(%rsp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%rsp),%xmm6
        mulps  nb310nf_c12(%rsp),%xmm4
        movaps nb310nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb310nf_Vvdwtot(%rsp)


        ## should we do one more iteration? 
        subl $4,nb310nf_innerk(%rsp)
        jl    _nb_kernel310nf_x86_64_sse.nb310nf_finish_inner
        jmp   _nb_kernel310nf_x86_64_sse.nb310nf_unroll_loop
_nb_kernel310nf_x86_64_sse.nb310nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb310nf_innerk(%rsp)
        movl  nb310nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel310nf_x86_64_sse.nb310nf_dopair
        jmp   _nb_kernel310nf_x86_64_sse.nb310nf_checksingle
_nb_kernel310nf_x86_64_sse.nb310nf_dopair: 
        movq nb310nf_charge(%rbp),%rsi
    movq  nb310nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb310nf_innerjjnr(%rsp)
        xorps %xmm7,%xmm7
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb310nf_iq(%rsp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb310nf_qq(%rsp)

        movq nb310nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb310nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb310nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb310nf_pos(%rbp),%rdi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb310nf_c6(%rsp)
        movaps %xmm6,nb310nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb310nf_ix(%rsp),%xmm4
        movaps nb310nf_iy(%rsp),%xmm5
        movaps nb310nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb310nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb310nf_VFtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        movlps (%rsi,%rcx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm5 ## got half coulomb table 
        movaps %xmm5,%xmm4
        shufps $136,%xmm4,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rcx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm7
        movaps %xmm7,%xmm6
        shufps $136,%xmm6,%xmm6 ## 10001000
        shufps $221,%xmm7,%xmm7 ## 11011101
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb310nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb310nf_vctot(%rsp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movaps %xmm5,nb310nf_vctot(%rsp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%rsp),%xmm6
        mulps  nb310nf_c12(%rsp),%xmm4
        movaps nb310nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb310nf_Vvdwtot(%rsp)

_nb_kernel310nf_x86_64_sse.nb310nf_checksingle: 
        movl  nb310nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel310nf_x86_64_sse.nb310nf_dosingle
        jmp    _nb_kernel310nf_x86_64_sse.nb310nf_updateouterdata
_nb_kernel310nf_x86_64_sse.nb310nf_dosingle: 
        movq nb310nf_charge(%rbp),%rsi
        movq nb310nf_pos(%rbp),%rdi
        movq  nb310nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm6,%xmm6
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb310nf_iq(%rsp),%xmm6
        movaps %xmm6,nb310nf_qq(%rsp)

        movq nb310nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb310nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb310nf_ntia(%rsp),%ecx
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb310nf_c6(%rsp)
        movaps %xmm6,nb310nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movaps nb310nf_ix(%rsp),%xmm4
        movaps nb310nf_iy(%rsp),%xmm5
        movaps nb310nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb310nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb310nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb310nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb310nf_VFtab(%rbp),%rsi
        movd %mm6,%ebx

        movlps (%rsi,%rbx,4),%xmm4
        movlps 8(%rsi,%rbx,4),%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb310nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb310nf_vctot(%rsp),%xmm5

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        movss %xmm5,nb310nf_vctot(%rsp)

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb310nf_c6(%rsp),%xmm6
        mulps  nb310nf_c12(%rsp),%xmm4
        movss nb310nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movss %xmm7,nb310nf_Vvdwtot(%rsp)

_nb_kernel310nf_x86_64_sse.nb310nf_updateouterdata: 
        ## get n from stack
        movl nb310nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb310nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb310nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb310nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb310nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb310nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb310nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310nf_x86_64_sse.nb310nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310nf_n(%rsp)
        jmp _nb_kernel310nf_x86_64_sse.nb310nf_outer
_nb_kernel310nf_x86_64_sse.nb310nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310nf_x86_64_sse.nb310nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310nf_x86_64_sse.nb310nf_threadloop
_nb_kernel310nf_x86_64_sse.nb310nf_end: 

        movl nb310nf_nouter(%rsp),%eax
        movl nb310nf_ninner(%rsp),%ebx
        movq nb310nf_outeriter(%rbp),%rcx
        movq nb310nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $312,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

