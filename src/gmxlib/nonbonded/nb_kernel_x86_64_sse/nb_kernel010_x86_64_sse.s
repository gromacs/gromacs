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




## nb010 - forces are calculated
.globl nb_kernel010_x86_64_sse
.globl _nb_kernel010_x86_64_sse
nb_kernel010_x86_64_sse:        
_nb_kernel010_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb010_fshift, 16
.set nb010_gid, 24
.set nb010_pos, 32
.set nb010_faction, 40
.set nb010_charge, 48
.set nb010_p_facel, 56
.set nb010_argkrf, 64
.set nb010_argcrf, 72
.set nb010_Vc, 80
.set nb010_type, 88
.set nb010_p_ntype, 96
.set nb010_vdwparam, 104
.set nb010_Vvdw, 112
.set nb010_p_tabscale, 120
.set nb010_VFtab, 128
.set nb010_invsqrta, 136
.set nb010_dvda, 144
.set nb010_p_gbtabscale, 152
.set nb010_GBtab, 160
.set nb010_p_nthreads, 168
.set nb010_count, 176
.set nb010_mtx, 184
.set nb010_outeriter, 192
.set nb010_inneriter, 200
.set nb010_work, 208
        ## The mutex (last arg) is not used in assembly.
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010_ix, 0
.set nb010_iy, 16
.set nb010_iz, 32
.set nb010_dx, 48
.set nb010_dy, 64
.set nb010_dz, 80
.set nb010_two, 96
.set nb010_c6, 112
.set nb010_c12, 128
.set nb010_six, 144
.set nb010_twelve, 160
.set nb010_Vvdwtot, 176
.set nb010_fix, 192
.set nb010_fiy, 208
.set nb010_fiz, 224
.set nb010_half, 240
.set nb010_three, 256
.set nb010_nri, 272
.set nb010_iinr, 280
.set nb010_jindex, 288
.set nb010_jjnr, 296
.set nb010_shift, 304
.set nb010_shiftvec, 312
.set nb010_facel, 320
.set nb010_innerjjnr, 328
.set nb010_is3, 336
.set nb010_ii3, 340
.set nb010_ntia, 344
.set nb010_innerk, 348
.set nb010_n, 352
.set nb010_nn1, 356
.set nb010_ntype, 360
.set nb010_nouter, 364
.set nb010_ninner, 368

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms

        subq $392,%rsp          # # local variable stack space (n*16+8)                                                         
        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb010_nouter(%rsp)
        movl %eax,nb010_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb010_nri(%rsp)
        movq %rsi,nb010_iinr(%rsp)
        movq %rdx,nb010_jindex(%rsp)
        movq %rcx,nb010_jjnr(%rsp)
        movq %r8,nb010_shift(%rsp)
        movq %r9,nb010_shiftvec(%rsp)
        movq nb010_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb010_ntype(%rsp)

    ## create constant floating-point factors on stack
    movl $0x40000000,%eax   ## 2.0 in IEEE (hex)
    movl %eax,nb010_two(%rsp)
    movss nb010_two(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm1,%xmm2      ## 4.0
        addps  %xmm1,%xmm2      ## 6.0
        movaps %xmm2,%xmm3
        addps  %xmm3,%xmm3      ## 12.0
        movaps %xmm1,nb010_two(%rsp)
    movaps %xmm2,nb010_six(%rsp)
    movaps %xmm3,nb010_twelve(%rsp)


_nb_kernel010_x86_64_sse.nb010_threadloop: 
    movq  nb010_count(%rbp),%rsi            ## pointer to sync counter
    movl  (%rsi),%eax
_nb_kernel010_x86_64_sse.nb010_spinlock: 
    movl  %eax,%ebx                         ## ebx=*count=nn0
    addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                ## write nn1 to *counter,
                                            ## if it hasnt changed.
                                            ## or reread *counter to eax.
    pause                                   ## -> better p4 performance
    jnz _nb_kernel010_x86_64_sse.nb010_spinlock

    ## if(nn1>nri) nn1=nri
    movl nb010_nri(%rsp),%ecx
    movl %ecx,%edx
    subl %ebx,%ecx
    cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
    ## Cleared the spinlock if we got here.
    ## eax contains nn0, ebx contains nn1.
    movl %eax,nb010_n(%rsp)
    movl %ebx,nb010_nn1(%rsp)
    subl %eax,%ebx                          ## calc number of outer lists
    movl %eax,%esi                              ## copy n to esi
    jg  _nb_kernel010_x86_64_sse.nb010_outerstart
    jmp _nb_kernel010_x86_64_sse.nb010_end

_nb_kernel010_x86_64_sse.nb010_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010_nouter(%rsp),%ebx
        movl %ebx,nb010_nouter(%rsp)

_nb_kernel010_x86_64_sse.nb010_outer: 
    movq  nb010_shift(%rsp),%rax                ## rax = base of shift[] 
    movl  (%rax,%rsi,4),%ebx                    ## ebx=shift[n] 

    lea  (%rbx,%rbx,2),%rbx                    ## rbx=3*is 
    movl  %ebx,nb010_is3(%rsp)          ## store is3 

    movq  nb010_shiftvec(%rsp),%rax             ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm10
        movss 4(%rax,%rbx,4),%xmm11
        movss 8(%rax,%rbx,4),%xmm12

    movq  nb010_iinr(%rsp),%rcx                 ## rcx = base of iinr[] 
    movl  (%rcx,%rsi,4),%ebx                    ## ebx =ii 

    movq nb010_type(%rbp),%rdx
    movl (%rdx,%rbx,4),%edx
    imull nb010_ntype(%rsp),%edx
    shll %edx
    movl %edx,nb010_ntia(%rsp)

    lea  (%rbx,%rbx,2),%rbx            ## rbx = 3*ii=ii3 
    movq  nb010_pos(%rbp),%rax          ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm10
        addss 4(%rax,%rbx,4),%xmm11
        addss 8(%rax,%rbx,4),%xmm12

    shufps $0,%xmm10,%xmm10
    shufps $0,%xmm11,%xmm11
    shufps $0,%xmm12,%xmm12

    movaps %xmm10,nb010_ix(%rsp)
    movaps %xmm11,nb010_iy(%rsp)
    movaps %xmm12,nb010_iz(%rsp)

    movl  %ebx,nb010_ii3(%rsp)

        ## clear vvdwtot (xmm12) and i forces (xmm13-xmm15)
        xorps %xmm12,%xmm12
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15

    movq  nb010_jindex(%rsp),%rax
    movl  (%rax,%rsi,4),%ecx                    ## jindex[n] 
    movl  4(%rax,%rsi,4),%edx                   ## jindex[n+1] 
    subl  %ecx,%edx                             ## number of innerloop atoms 

    movq  nb010_jjnr(%rsp),%rax
    shll  $2,%ecx
    addq  %rcx,%rax
    movq  %rax,nb010_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
    subl  $4,%edx
        addl  nb010_ninner(%rsp),%ecx
        movl  %ecx,nb010_ninner(%rsp)
        addl  $0,%edx
    movl  %edx,nb010_innerk(%rsp)           ## number of innerloop atoms 

    jge   _nb_kernel010_x86_64_sse.nb010_unroll_loop
    jmp   _nb_kernel010_x86_64_sse.nb010_finish_inner
_nb_kernel010_x86_64_sse.nb010_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb010_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb010_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%rax,%rax,2),%r8     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%r9
        lea  (%rcx,%rcx,2),%r10
        lea  (%rdx,%rdx,2),%r11

        movq nb010_pos(%rbp),%rdi
        ## load coordinates
        movlps (%rdi,%r8,4),%xmm1       ## x1 y1 - - 
        movlps (%rdi,%r10,4),%xmm2      ## x3 y3 - - 
        movhps (%rdi,%r9,4),%xmm1       ## x2 y2 - -
        movhps (%rdi,%r11,4),%xmm2      ## x4 y4 - -

        movss 8(%rdi,%r8,4),%xmm5       ## z1 - - - 
        movss 8(%rdi,%r10,4),%xmm6      ## z2 - - - 
        movss 8(%rdi,%r9,4),%xmm7       ## z3 - - - 
        movss 8(%rdi,%r11,4),%xmm8      ## z4 - - - 
    movlhps %xmm7,%xmm5 ## jzOa  -  jzOb  -
    movlhps %xmm8,%xmm6 ## jzOc  -  jzOd -

        movq nb010_type(%rbp),%rsi

    movaps %xmm1,%xmm4
    unpcklps %xmm2,%xmm1 ## jxa jxc jya jyc        
    unpckhps %xmm2,%xmm4 ## jxb jxd jyb jyd
    movaps %xmm1,%xmm2
    unpcklps %xmm4,%xmm1 ## x
    unpckhps %xmm4,%xmm2 ## y
    shufps  $136,%xmm6,%xmm5  ## 10001000 => jzH2a jzH2b jzH2c jzH2d

    ## load vdw types
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        movl (%rsi,%rcx,4),%r14d
        movl (%rsi,%rdx,4),%r15d

        ## calc dr  
        subps nb010_ix(%rsp),%xmm1
        subps nb010_iy(%rsp),%xmm2
        subps nb010_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

    ## type *=2
        shll %r12d
        shll %r13d
        shll %r14d
        shll %r15d

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

    ## 2*type*ntia
    movl nb010_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d
        addl %edi,%r14d
        addl %edi,%r15d

        movq nb010_vdwparam(%rbp),%rsi
    ## xmm0=c6
    ## xmm3=c12

        rcpps %xmm1,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010_two(%rsp),%xmm6
        mulps %xmm5,%xmm1
    ## load c6/c12
        movlps (%rsi,%r12,4),%xmm7
        movlps (%rsi,%r14,4),%xmm8

        subps %xmm1,%xmm6
        mulps %xmm5,%xmm6       ## xmm6=rinvsq

        movaps %xmm6,%xmm4  ## rinvsq

        movhps (%rsi,%r13,4),%xmm7
        movhps (%rsi,%r15,4),%xmm8

        movaps %xmm6,%xmm1
        mulps  %xmm6,%xmm1  ## rinv4
        mulps  %xmm6,%xmm1      ## rinv6
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinv12

    ## shuffle c6/c12
        movaps %xmm7,%xmm5
        shufps $136,%xmm8,%xmm5 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101

        movq nb010_faction(%rbp),%rsi

        mulps  %xmm5,%xmm1 ## c6*rinv6
        mulps  %xmm7,%xmm2 ## c12*rinv12
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulps  nb010_six(%rsp),%xmm1
        mulps  nb010_twelve(%rsp),%xmm2
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 

        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r10,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r11,4),%xmm1 ## x3 y3 x4 y4

    ## add potential to Vvdwtot (sum in xmm12)
        addps  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

    ## permute local forces
    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm5 ## fjz2 - - -
    movhlps %xmm11,%xmm4     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm3  ## fjz4 - - -

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%r8,4)
        movlps %xmm1,(%rsi,%r10,4)
        movhps %xmm0,(%rsi,%r9,4)
        movhps %xmm1,(%rsi,%r11,4)

        addss  8(%rsi,%r8,4),%xmm11
        addss  8(%rsi,%r9,4),%xmm5
        addss  8(%rsi,%r10,4),%xmm4
        addss  8(%rsi,%r11,4),%xmm3
        movss  %xmm11,8(%rsi,%r8,4)
        movss  %xmm5,8(%rsi,%r9,4)
        movss  %xmm4,8(%rsi,%r10,4)
        movss  %xmm3,8(%rsi,%r11,4)

        ## should we do one more iteration? 
        subl $4,nb010_innerk(%rsp)
        jl    _nb_kernel010_x86_64_sse.nb010_finish_inner
        jmp   _nb_kernel010_x86_64_sse.nb010_unroll_loop
_nb_kernel010_x86_64_sse.nb010_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb010_innerk(%rsp)
    movl  nb010_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel010_x86_64_sse.nb010_dopair
    jmp   _nb_kernel010_x86_64_sse.nb010_checksingle
_nb_kernel010_x86_64_sse.nb010_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb010_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb010_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb010_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb010_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb010_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps  %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    ## xmm0=c6
    ## xmm3=c12
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        movq nb010_pos(%rbp),%rdi
        ## load coordinates
        movlps (%rdi,%rax,4),%xmm1      ## x1 y1  -  - 
        movlps (%rdi,%rbx,4),%xmm4      ## x2 y2  -  -

        movss 8(%rdi,%rax,4),%xmm5      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm4,%xmm1 ## x1 x2 y1 y2
    movhlps  %xmm1,%xmm2 ## y1 y2 -  -
    unpcklps %xmm7,%xmm5 ## z1 z2 -  - 

        ## calc dr  
        subps nb010_ix(%rsp),%xmm1
        subps nb010_iy(%rsp),%xmm2
        subps nb010_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm5,%xmm5
        addps %xmm2,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

        rcpps %xmm1,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010_two(%rsp),%xmm6
        mulps %xmm5,%xmm1
        subps %xmm1,%xmm6
        mulps %xmm5,%xmm6       ## xmm6=rinvsq

        movaps %xmm6,%xmm4  ## rinvsq

        movaps %xmm6,%xmm1
        mulps  %xmm6,%xmm1  ## rinv4
        mulps  %xmm6,%xmm1      ## rinv6
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinv12

        mulps  %xmm0,%xmm1
        mulps  %xmm3,%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulps  nb010_six(%rsp),%xmm1
        mulps  nb010_twelve(%rsp),%xmm2
        subps  %xmm1,%xmm2
        mulps  %xmm2,%xmm4      ## xmm4=total fscal 

    xorps  %xmm7,%xmm7
    movlhps %xmm7,%xmm5

    ## add potential to Vvdwtot (sum in xmm12)
        addps  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm7,%xmm9
    movlhps %xmm7,%xmm10
    movlhps %xmm7,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb010_faction(%rbp),%rsi
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

_nb_kernel010_x86_64_sse.nb010_checksingle:     
    movl  nb010_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel010_x86_64_sse.nb010_dosingle
    jmp    _nb_kernel010_x86_64_sse.nb010_updateouterdata

_nb_kernel010_x86_64_sse.nb010_dosingle: 
    movq nb010_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb010_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb010_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb010_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
    movss 4(%rsi,%r12,4),%xmm3

    ## xmm0=c6
    ## xmm3=c12

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb010_pos(%rbp),%rdi
        ## load coordinates
        movss (%rdi,%rax,4),%xmm1
        movss 4(%rdi,%rax,4),%xmm2
        movss 8(%rdi,%rax,4),%xmm5

        ## calc dr  
        subss nb010_ix(%rsp),%xmm1
        subss nb010_iy(%rsp),%xmm2
        subss nb010_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm2,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulss %xmm1,%xmm1
        mulss %xmm2,%xmm2
        mulss %xmm5,%xmm5
        addss %xmm2,%xmm1
        addss %xmm5,%xmm1
        ## rsq in xmm1

        ## rsq in xmm4 
        rcpss %xmm1,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010_two(%rsp),%xmm6
        mulss %xmm5,%xmm1
        subss %xmm1,%xmm6
        mulss %xmm5,%xmm6       ## xmm6=rinvsq

        movaps %xmm6,%xmm4  ## rinvsq

        movaps %xmm6,%xmm1
        mulss  %xmm6,%xmm1  ## rinv4
        mulss  %xmm6,%xmm1      ## rinv6
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinv12

        mulss  %xmm0,%xmm1
        mulss  %xmm3,%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulss  nb010_six(%rsp),%xmm1
        mulss  nb010_twelve(%rsp),%xmm2
        subss  %xmm1,%xmm2
        mulss  %xmm2,%xmm4      ## xmm4=total fscal 

    ## add potential to Vvdwtot (sum in xmm12)
        addss  %xmm5,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb010_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel010_x86_64_sse.nb010_updateouterdata: 
        movl  nb010_ii3(%rsp),%ecx
        movq  nb010_faction(%rbp),%rdi
        movq  nb010_fshift(%rbp),%rsi
        movl  nb010_is3(%rsp),%edx

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
        movl nb010_n(%rsp),%esi
    ## get group index for i particle 
    movq  nb010_gid(%rbp),%rdx          ## base of gid[]
    movl  (%rdx,%rsi,4),%edx            ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb010_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## finish if last 
        movl nb010_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010_x86_64_sse.nb010_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010_n(%rsp)
        jmp _nb_kernel010_x86_64_sse.nb010_outer
_nb_kernel010_x86_64_sse.nb010_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010_x86_64_sse.nb010_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010_x86_64_sse.nb010_threadloop
_nb_kernel010_x86_64_sse.nb010_end: 

        emms

        movl nb010_nouter(%rsp),%eax
        movl nb010_ninner(%rsp),%ebx
        movq nb010_outeriter(%rbp),%rcx
        movq nb010_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $392,%rsp

        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel010nf_x86_64_sse
.globl _nb_kernel010nf_x86_64_sse
nb_kernel010nf_x86_64_sse:      
_nb_kernel010nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb010nf_fshift, 16
.set nb010nf_gid, 24
.set nb010nf_pos, 32
.set nb010nf_faction, 40
.set nb010nf_charge, 48
.set nb010nf_p_facel, 56
.set nb010nf_argkrf, 64
.set nb010nf_argcrf, 72
.set nb010nf_Vc, 80
.set nb010nf_type, 88
.set nb010nf_p_ntype, 96
.set nb010nf_vdwparam, 104
.set nb010nf_Vvdw, 112
.set nb010nf_p_tabscale, 120
.set nb010nf_VFtab, 128
.set nb010nf_invsqrta, 136
.set nb010nf_dvda, 144
.set nb010nf_p_gbtabscale, 152
.set nb010nf_GBtab, 160
.set nb010nf_p_nthreads, 168
.set nb010nf_count, 176
.set nb010nf_mtx, 184
.set nb010nf_outeriter, 192
.set nb010nf_inneriter, 200
.set nb010nf_work, 208
        ## The mutex (last arg) is not used in assembly.
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010nf_ix, 0
.set nb010nf_iy, 16
.set nb010nf_iz, 32
.set nb010nf_two, 48
.set nb010nf_c6, 64
.set nb010nf_c12, 80
.set nb010nf_Vvdwtot, 96
.set nb010nf_half, 112
.set nb010nf_three, 128
.set nb010nf_nri, 144
.set nb010nf_iinr, 152
.set nb010nf_jindex, 160
.set nb010nf_jjnr, 168
.set nb010nf_shift, 176
.set nb010nf_shiftvec, 184
.set nb010nf_innerjjnr, 192
.set nb010nf_facel, 200
.set nb010nf_ntia, 208
.set nb010nf_innerk, 216
.set nb010nf_is3, 220
.set nb010nf_ii3, 224
.set nb010nf_n, 228
.set nb010nf_nn1, 232
.set nb010nf_ntype, 236
.set nb010nf_nouter, 240
.set nb010nf_ninner, 244

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $264,%rsp          # # local variable stack space (n*16+8)                                                         
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb010nf_nouter(%rsp)
        movl %eax,nb010nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb010nf_nri(%rsp)
        movq %rsi,nb010nf_iinr(%rsp)
        movq %rdx,nb010nf_jindex(%rsp)
        movq %rcx,nb010nf_jjnr(%rsp)
        movq %r8,nb010nf_shift(%rsp)
        movq %r9,nb010nf_shiftvec(%rsp)
        movq nb010nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb010nf_ntype(%rsp)

        ## create constant floating-point factors on stack
        movl $0x40000000,%eax   ## 2.0 in IEEE (hex)
        movl %eax,nb010nf_two(%rsp)
        movss nb010nf_two(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,nb010nf_two(%rsp)

_nb_kernel010nf_x86_64_sse.nb010nf_threadloop: 
        movq  nb010nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel010nf_x86_64_sse.nb010nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010nf_x86_64_sse.nb010nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010nf_n(%rsp)
        movl %ebx,nb010nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010nf_x86_64_sse.nb010nf_outerstart
        jmp _nb_kernel010nf_x86_64_sse.nb010nf_end

_nb_kernel010nf_x86_64_sse.nb010nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010nf_nouter(%rsp),%ebx
        movl %ebx,nb010nf_nouter(%rsp)

_nb_kernel010nf_x86_64_sse.nb010nf_outer: 
        movq  nb010nf_shift(%rsp),%rax          ## rax = base of shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx                ## rbx=3*is 
        movl  %ebx,nb010nf_is3(%rsp)            ## store is3 

        movq  nb010nf_shiftvec(%rsp),%rax       ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb010nf_iinr(%rsp),%rcx           ## rcx = base of iinr[] 
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movq  nb010nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb010nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb010nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb010nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb010nf_ix(%rsp)
        movaps %xmm1,nb010nf_iy(%rsp)
        movaps %xmm2,nb010nf_iz(%rsp)

        movl  %ebx,nb010nf_ii3(%rsp)

        ## clear Vvdwtot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb010nf_Vvdwtot(%rsp)

        movq  nb010nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb010nf_pos(%rbp),%rsi
        movq  nb010nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb010nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb010nf_ninner(%rsp),%ecx
        movl  %ecx,nb010nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb010nf_innerk(%rsp)         ## number of innerloop atoms 

        jge   _nb_kernel010nf_x86_64_sse.nb010nf_unroll_loop
        jmp   _nb_kernel010nf_x86_64_sse.nb010nf_finish_inner
_nb_kernel010nf_x86_64_sse.nb010nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb010nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        ## advance pointer (unrolled 4) 
        addq  $16,nb010nf_innerjjnr(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movq nb010nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb010nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb010nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb010nf_c6(%rsp)
        movaps %xmm6,nb010nf_c12(%rsp)

        movq nb010nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        mulps %xmm2,%xmm3
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
        movaps nb010nf_ix(%rsp),%xmm4
        movaps nb010nf_iy(%rsp),%xmm5
        movaps nb010nf_iz(%rsp),%xmm6

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
        rcpps %xmm4,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010nf_two(%rsp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%rsp),%xmm1
        mulps  nb010nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb010nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl  $4,nb010nf_innerk(%rsp)
        jl    _nb_kernel010nf_x86_64_sse.nb010nf_finish_inner
        jmp   _nb_kernel010nf_x86_64_sse.nb010nf_unroll_loop
_nb_kernel010nf_x86_64_sse.nb010nf_finish_inner: 
        ## check if at least two particles remain 
        addl  $4,nb010nf_innerk(%rsp)
        movl  nb010nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel010nf_x86_64_sse.nb010nf_dopair
        jmp   _nb_kernel010nf_x86_64_sse.nb010nf_checksingle
_nb_kernel010nf_x86_64_sse.nb010nf_dopair: 
        movq  nb010nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq  $8,nb010nf_innerjjnr(%rsp)

        movq nb010nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb010nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb010nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb010nf_pos(%rbp),%rdi
        xorps  %xmm7,%xmm7
        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb010nf_c6(%rsp)
        movaps %xmm6,nb010nf_c12(%rsp)

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

        ## move nb010nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb010nf_ix(%rsp),%xmm4
        movaps nb010nf_iy(%rsp),%xmm5
        movaps nb010nf_iz(%rsp),%xmm6

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


        rcpps %xmm4,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010nf_two(%rsp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%rsp),%xmm1
        mulps  nb010nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb010nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm5,nb010nf_Vvdwtot(%rsp)

_nb_kernel010nf_x86_64_sse.nb010nf_checksingle: 
        movl  nb010nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010nf_x86_64_sse.nb010nf_dosingle
        jmp    _nb_kernel010nf_x86_64_sse.nb010nf_updateouterdata
_nb_kernel010nf_x86_64_sse.nb010nf_dosingle: 
        movq nb010nf_pos(%rbp),%rdi
        movq  nb010nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb010nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb010nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb010nf_ntia(%rsp),%ecx
        xorps  %xmm6,%xmm6
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb010nf_c6(%rsp)
        movaps %xmm6,nb010nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        xorps   %xmm7,%xmm7

        movaps nb010nf_ix(%rsp),%xmm4
        movaps nb010nf_iy(%rsp),%xmm5
        movaps nb010nf_iz(%rsp),%xmm6

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

        rcpps %xmm4,%xmm5
        ## 1/x lookup seed in xmm5 
        movaps nb010nf_two(%rsp),%xmm0
        mulps %xmm5,%xmm4
        subps %xmm4,%xmm0
        mulps %xmm5,%xmm0       ## xmm0=rinvsq 
        movaps %xmm0,%xmm4

        movaps %xmm0,%xmm1
        mulps  %xmm0,%xmm1
        mulps  %xmm0,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulps  nb010nf_c6(%rsp),%xmm1
        mulps  nb010nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subps  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addss  nb010nf_Vvdwtot(%rsp),%xmm5
        movss %xmm5,nb010nf_Vvdwtot(%rsp)

_nb_kernel010nf_x86_64_sse.nb010nf_updateouterdata: 
        ## get n from stack
        movl nb010nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb010nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movaps nb010nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb010nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb010nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010nf_x86_64_sse.nb010nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010nf_n(%rsp)
        jmp _nb_kernel010nf_x86_64_sse.nb010nf_outer
_nb_kernel010nf_x86_64_sse.nb010nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010nf_x86_64_sse.nb010nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010nf_x86_64_sse.nb010nf_threadloop
_nb_kernel010nf_x86_64_sse.nb010nf_end: 

        movl nb010nf_nouter(%rsp),%eax
        movl nb010nf_ninner(%rsp),%ebx
        movq nb010nf_outeriter(%rbp),%rcx
        movq nb010nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $264,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret

