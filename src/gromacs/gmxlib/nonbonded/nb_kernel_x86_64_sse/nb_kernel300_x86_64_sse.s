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






.globl nb_kernel300_x86_64_sse
.globl _nb_kernel300_x86_64_sse
nb_kernel300_x86_64_sse:        
_nb_kernel300_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb300_fshift, 16
.set nb300_gid, 24
.set nb300_pos, 32
.set nb300_faction, 40
.set nb300_charge, 48
.set nb300_p_facel, 56
.set nb300_argkrf, 64
.set nb300_argcrf, 72
.set nb300_Vc, 80
.set nb300_type, 88
.set nb300_p_ntype, 96
.set nb300_vdwparam, 104
.set nb300_Vvdw, 112
.set nb300_p_tabscale, 120
.set nb300_VFtab, 128
.set nb300_invsqrta, 136
.set nb300_dvda, 144
.set nb300_p_gbtabscale, 152
.set nb300_GBtab, 160
.set nb300_p_nthreads, 168
.set nb300_count, 176
.set nb300_mtx, 184
.set nb300_outeriter, 192
.set nb300_inneriter, 200
.set nb300_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb300_ix, 0
.set nb300_iy, 16
.set nb300_iz, 32
.set nb300_iq, 48
.set nb300_dx, 64
.set nb300_dy, 80
.set nb300_dz, 96
.set nb300_two, 112
.set nb300_tsc, 128
.set nb300_qq, 144
.set nb300_fs, 160
.set nb300_vctot, 176
.set nb300_fix, 192
.set nb300_fiy, 208
.set nb300_fiz, 224
.set nb300_half, 240
.set nb300_three, 256
.set nb300_innerjjnr, 272
.set nb300_nri, 280
.set nb300_iinr, 288
.set nb300_jindex, 296
.set nb300_jjnr, 304
.set nb300_shift, 312
.set nb300_shiftvec, 320
.set nb300_facel, 328
.set nb300_innerk, 336
.set nb300_is3, 344
.set nb300_ii3, 348
.set nb300_n, 352
.set nb300_nn1, 356
.set nb300_nouter, 360
.set nb300_ninner, 364


        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $376,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb300_nouter(%rsp)
        movl %eax,nb300_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb300_nri(%rsp)
        movq %rsi,nb300_iinr(%rsp)
        movq %rdx,nb300_jindex(%rsp)
        movq %rcx,nb300_jjnr(%rsp)
        movq %r8,nb300_shift(%rsp)
        movq %r9,nb300_shiftvec(%rsp)
        movq nb300_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb300_facel(%rsp)

        movq nb300_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb300_tsc(%rsp)

        movq nb300_pos(%rbp),%r8
        movq nb300_faction(%rbp),%r9
        movq nb300_charge(%rbp),%r10
        movq nb300_VFtab(%rbp),%r11

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb300_half(%rsp)
        movss nb300_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb300_half(%rsp)
        movaps %xmm2,nb300_two(%rsp)
        movaps %xmm3,nb300_three(%rsp)

_nb_kernel300_x86_64_sse.nb300_threadloop: 
        movq  nb300_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel300_x86_64_sse.nb300_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300_x86_64_sse.nb300_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300_n(%rsp)
        movl %ebx,nb300_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300_x86_64_sse.nb300_outerstart
        jmp _nb_kernel300_x86_64_sse.nb300_end

_nb_kernel300_x86_64_sse.nb300_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300_nouter(%rsp),%ebx
        movl %ebx,nb300_nouter(%rsp)

_nb_kernel300_x86_64_sse.nb300_outer: 
        movq  nb300_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb300_is3(%rsp)      ## store is3 

        movq  nb300_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb300_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb300_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb300_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb300_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb300_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb300_ix(%rsp)
        movaps %xmm1,nb300_iy(%rsp)
        movaps %xmm2,nb300_iz(%rsp)

        movl  %ebx,nb300_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorps %xmm12,%xmm12
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15

        movq  nb300_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb300_pos(%rbp),%rsi
        movq  nb300_faction(%rbp),%rdi
        movq  nb300_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb300_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb300_ninner(%rsp),%ecx
        movl  %ecx,nb300_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb300_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel300_x86_64_sse.nb300_unroll_loop
        jmp   _nb_kernel300_x86_64_sse.nb300_finish_inner
_nb_kernel300_x86_64_sse.nb300_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb300_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 
        movq nb300_pos(%rbp),%rdi

        addq $16,nb300_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx     ## replace jnr with j3 
        lea  (%r11,%r11,2),%rdx

        ## load coordinates

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

    movaps %xmm1,%xmm3

    unpcklps %xmm2,%xmm1 ## x1 x2 x3 x4
    unpckhps %xmm2,%xmm3 ## y1 y2 y3 y4
    unpcklps %xmm6,%xmm5 ## z1 z2 z3 z4

        ## calc dr  
        subps nb300_ix(%rsp),%xmm1
        subps nb300_iy(%rsp),%xmm3
        subps nb300_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm3,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm3,%xmm3
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm3
        addps %xmm5,%xmm3
        ## rsq in xmm3

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm3,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb300_three(%rsp),%xmm1
        mulps %xmm3,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps nb300_half(%rsp),%xmm1
    ## xmm1=rinv
    ## xmm3=rsq

    mulps %xmm1,%xmm3 ## r
    mulps nb300_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttps2dq %xmm3,%xmm2

    ## convert back to float
    cvtdq2ps  %xmm2,%xmm0

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

        movq nb300_VFtab(%rbp),%rsi
    ## load table data
        movlps (%rsi,%r12,4),%xmm5
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm5
        movhps (%rsi,%r15,4),%xmm7

    movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%r12,4),%xmm7
        movlps 8(%rsi,%r14,4),%xmm8
        movhps 8(%rsi,%r13,4),%xmm7
        movhps 8(%rsi,%r15,4),%xmm8

    movaps %xmm7,%xmm6
        movq nb300_charge(%rbp),%rdi

        shufps $136,%xmm8,%xmm6 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101
    ## table data ready in xmm4-xmm7

        movss (%rdi,%r8,4),%xmm0
        movss (%rdi,%r10,4),%xmm2

    mulps %xmm3,%xmm7  ## Heps
    mulps  %xmm3,%xmm6 ## Geps
    mulps %xmm3,%xmm7  ## Heps2

    unpcklps %xmm2,%xmm0
        movss (%rdi,%r9,4),%xmm2
        movss (%rdi,%r11,4),%xmm8

    addps  %xmm6,%xmm5  ## F+Geps
    addps  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm7  ## 2*Heps2
    unpcklps %xmm8,%xmm2
    unpcklps %xmm2,%xmm0
    addps  %xmm6,%xmm7  ## 2*Heps2+Geps
    addps  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulps  %xmm3,%xmm5  ## eps*Fp
    mulps  nb300_iq(%rsp),%xmm0
    addps  %xmm4,%xmm5  ## VV
    mulps  %xmm0,%xmm5  ## VV*qq=vcoul
    mulps  %xmm0,%xmm7  ## FF*qq=fijC

    ## add potential to vctot (sum in xmm12)
        addps  %xmm5,%xmm12

    mulps  nb300_tsc(%rsp),%xmm7
    mulps  %xmm1,%xmm7

    xorps  %xmm4,%xmm4
    subps  %xmm7,%xmm4  ## fscal

        movq nb300_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%rax,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%rcx,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%rbx,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%rdx,4),%xmm1 ## x3 y3 x4 y4

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

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
        subl $4,nb300_innerk(%rsp)
        jl    _nb_kernel300_x86_64_sse.nb300_finish_inner
        jmp   _nb_kernel300_x86_64_sse.nb300_unroll_loop
_nb_kernel300_x86_64_sse.nb300_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb300_innerk(%rsp)
    movl  nb300_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel300_x86_64_sse.nb300_dopair
    jmp   _nb_kernel300_x86_64_sse.nb300_checksingle
_nb_kernel300_x86_64_sse.nb300_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb300_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb300_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb300_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm8
        movss (%rsi,%rbx,4),%xmm1

    unpcklps %xmm1,%xmm8 ## jqa jqb - -
        mulps nb300_iq(%rsp),%xmm8      ##qq

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load coordinates
        movq nb300_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm4      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm5      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm6      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr  
        subps nb300_ix(%rsp),%xmm4
        subps nb300_iy(%rsp),%xmm5
        subps nb300_iz(%rsp),%xmm6

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
        movaps nb300_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps nb300_half(%rsp),%xmm1
    ## xmm1=rinv
    movaps %xmm4,%xmm3
    ## xmm3=rsq 

    mulps %xmm1,%xmm3 ## r
    mulps nb300_tsc(%rsp),%xmm3   ## rtab

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

        movq nb300_VFtab(%rbp),%rsi
    ## load table data
        movlps (%rsi,%r12,4),%xmm4
        movlps (%rsi,%r13,4),%xmm5
    unpcklps %xmm5,%xmm4
    movhlps  %xmm4,%xmm5

        movlps 8(%rsi,%r12,4),%xmm6
        movlps 8(%rsi,%r13,4),%xmm7
    unpcklps %xmm7,%xmm6
    movhlps  %xmm6,%xmm7
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
    mulps  %xmm8,%xmm5  ## VV*qq=vcoul
    mulps  %xmm8,%xmm7  ## FF*qq=fijC

    xorps %xmm6,%xmm6
    movlhps %xmm6,%xmm5

    ## add potential to vctot (sum in xmm12)
        addps  %xmm5,%xmm12

    mulps  nb300_tsc(%rsp),%xmm7
    mulps  %xmm1,%xmm7

    xorps  %xmm4,%xmm4
    subps  %xmm7,%xmm4  ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm6,%xmm9
    movlhps %xmm6,%xmm10
    movlhps %xmm6,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb300_faction(%rbp),%rsi
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

_nb_kernel300_x86_64_sse.nb300_checksingle:     
    movl  nb300_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel300_x86_64_sse.nb300_dosingle
    jmp    _nb_kernel300_x86_64_sse.nb300_updateouterdata

_nb_kernel300_x86_64_sse.nb300_dosingle: 
    movq nb300_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb300_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm8       ## jq
        mulss nb300_iq(%rsp),%xmm8      ## qq

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        movq nb300_pos(%rbp),%rdi
        movss (%rdi,%rax,4),%xmm4           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm5       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm6       ## 13 - - - 

        ## calc dr  
        subss nb300_ix(%rsp),%xmm4
        subss nb300_iy(%rsp),%xmm5
        subss nb300_iz(%rsp),%xmm6

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
        movaps nb300_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss nb300_half(%rsp),%xmm1
    ## xmm1=rinv
    movaps %xmm4,%xmm3
    ## xmm3=rsq 

    mulss %xmm1,%xmm3 ## r
    mulss nb300_tsc(%rsp),%xmm3   ## rtab

    ## truncate and convert to integers
    cvttss2si %xmm3,%r12d

    ## convert back to float
    cvtsi2ss  %r12d,%xmm0

    ## multiply by 4
    shll      $2,%r12d

    ## calculate eps
    subss     %xmm0,%xmm3

        movq nb300_VFtab(%rbp),%rsi
    ## load table data
        movss  (%rsi,%r12,4),%xmm4
        movss  4(%rsi,%r12,4),%xmm5
        movss  8(%rsi,%r12,4),%xmm6
        movss  12(%rsi,%r12,4),%xmm7
    ## table data ready in xmm4-xmm7

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
    mulss  %xmm8,%xmm5  ## VV*qq=vcoul
    mulss  %xmm8,%xmm7  ## FF*qq=fijC

    ## add potential to vctot (sum in xmm12)
        addss  %xmm5,%xmm12

    mulss  nb300_tsc(%rsp),%xmm7
    mulss  %xmm1,%xmm7

    xorps  %xmm4,%xmm4
    subss  %xmm7,%xmm4  ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb300_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel300_x86_64_sse.nb300_updateouterdata: 
        movl  nb300_ii3(%rsp),%ecx
        movq  nb300_faction(%rbp),%rdi
        movq  nb300_fshift(%rbp),%rsi
        movl  nb300_is3(%rsp),%edx

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
        movl nb300_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb300_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb300_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## finish if last 
        movl nb300_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300_x86_64_sse.nb300_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300_n(%rsp)
        jmp _nb_kernel300_x86_64_sse.nb300_outer
_nb_kernel300_x86_64_sse.nb300_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300_x86_64_sse.nb300_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300_x86_64_sse.nb300_threadloop
_nb_kernel300_x86_64_sse.nb300_end: 
        movl nb300_nouter(%rsp),%eax
        movl nb300_ninner(%rsp),%ebx
        movq nb300_outeriter(%rbp),%rcx
        movq nb300_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $376,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel300nf_x86_64_sse
.globl _nb_kernel300nf_x86_64_sse
nb_kernel300nf_x86_64_sse:      
_nb_kernel300nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb300nf_fshift, 16
.set nb300nf_gid, 24
.set nb300nf_pos, 32
.set nb300nf_faction, 40
.set nb300nf_charge, 48
.set nb300nf_p_facel, 56
.set nb300nf_argkrf, 64
.set nb300nf_argcrf, 72
.set nb300nf_Vc, 80
.set nb300nf_type, 88
.set nb300nf_p_ntype, 96
.set nb300nf_vdwparam, 104
.set nb300nf_Vvdw, 112
.set nb300nf_p_tabscale, 120
.set nb300nf_VFtab, 128
.set nb300nf_invsqrta, 136
.set nb300nf_dvda, 144
.set nb300nf_p_gbtabscale, 152
.set nb300nf_GBtab, 160
.set nb300nf_p_nthreads, 168
.set nb300nf_count, 176
.set nb300nf_mtx, 184
.set nb300nf_outeriter, 192
.set nb300nf_inneriter, 200
.set nb300nf_work, 208
        ## bottom of stack is cache-aligned for sse use 
.set nb300nf_ix, 0
.set nb300nf_iy, 16
.set nb300nf_iz, 32
.set nb300nf_iq, 48
.set nb300nf_tsc, 64
.set nb300nf_qq, 80
.set nb300nf_vctot, 96
.set nb300nf_half, 112
.set nb300nf_three, 128
.set nb300nf_is3, 144
.set nb300nf_ii3, 148
.set nb300nf_innerjjnr, 152
.set nb300nf_nri, 160
.set nb300nf_iinr, 168
.set nb300nf_jindex, 176
.set nb300nf_jjnr, 184
.set nb300nf_shift, 192
.set nb300nf_shiftvec, 200
.set nb300nf_facel, 208
.set nb300nf_innerk, 216
.set nb300nf_n, 220
.set nb300nf_nn1, 224
.set nb300nf_nouter, 228
.set nb300nf_ninner, 232

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $248,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb300nf_nouter(%rsp)
        movl %eax,nb300nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb300nf_nri(%rsp)
        movq %rsi,nb300nf_iinr(%rsp)
        movq %rdx,nb300nf_jindex(%rsp)
        movq %rcx,nb300nf_jjnr(%rsp)
        movq %r8,nb300nf_shift(%rsp)
        movq %r9,nb300nf_shiftvec(%rsp)
        movq nb300nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb300nf_facel(%rsp)

        movq nb300nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb300nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb300nf_half(%rsp)
        movss nb300nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb300nf_half(%rsp)
        movaps %xmm3,nb300nf_three(%rsp)

_nb_kernel300nf_x86_64_sse.nb300nf_threadloop: 
        movq  nb300nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel300nf_x86_64_sse.nb300nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300nf_x86_64_sse.nb300nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300nf_n(%rsp)
        movl %ebx,nb300nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300nf_x86_64_sse.nb300nf_outerstart
        jmp _nb_kernel300nf_x86_64_sse.nb300nf_end
_nb_kernel300nf_x86_64_sse.nb300nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300nf_nouter(%rsp),%ebx
        movl %ebx,nb300nf_nouter(%rsp)

_nb_kernel300nf_x86_64_sse.nb300nf_outer: 
        movq  nb300nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb300nf_is3(%rsp)            ## store is3 

        movq  nb300nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb300nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb300nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb300nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb300nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb300nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb300nf_ix(%rsp)
        movaps %xmm1,nb300nf_iy(%rsp)
        movaps %xmm2,nb300nf_iz(%rsp)

        movl  %ebx,nb300nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb300nf_vctot(%rsp)

        movq  nb300nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb300nf_pos(%rbp),%rsi
        movq  nb300nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb300nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb300nf_ninner(%rsp),%ecx
        movl  %ecx,nb300nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb300nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel300nf_x86_64_sse.nb300nf_unroll_loop
        jmp   _nb_kernel300nf_x86_64_sse.nb300nf_finish_inner
_nb_kernel300nf_x86_64_sse.nb300nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb300nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb300nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb300nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb300nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3

        movaps %xmm3,nb300nf_qq(%rsp)

        movq nb300nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb300nf_ix(%rsp),%xmm4
        movaps nb300nf_iy(%rsp),%xmm5
        movaps nb300nf_iz(%rsp),%xmm6

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
        movaps nb300nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%rsp),%xmm4

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

        movq nb300nf_VFtab(%rbp),%rsi
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
        movaps nb300nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  

        ## at this point xmm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb300nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb300nf_innerk(%rsp)
        jl    _nb_kernel300nf_x86_64_sse.nb300nf_finish_inner
        jmp   _nb_kernel300nf_x86_64_sse.nb300nf_unroll_loop
_nb_kernel300nf_x86_64_sse.nb300nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb300nf_innerk(%rsp)
        movl  nb300nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel300nf_x86_64_sse.nb300nf_dopair
        jmp   _nb_kernel300nf_x86_64_sse.nb300nf_checksingle
_nb_kernel300nf_x86_64_sse.nb300nf_dopair: 
        movq nb300nf_charge(%rbp),%rsi

    movq  nb300nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb300nf_innerjjnr(%rsp)
        xorps %xmm7,%xmm7
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb300nf_iq(%rsp),%xmm3
        movlhps %xmm7,%xmm3
        movaps %xmm3,nb300nf_qq(%rsp)

        movq nb300nf_pos(%rbp),%rdi

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

        movaps nb300nf_ix(%rsp),%xmm4
        movaps nb300nf_iy(%rsp),%xmm5
        movaps nb300nf_iz(%rsp),%xmm6

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
        movaps nb300nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb300nf_VFtab(%rbp),%rsi
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
        movaps nb300nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addps  nb300nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb300nf_vctot(%rsp)

_nb_kernel300nf_x86_64_sse.nb300nf_checksingle: 
        movl  nb300nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300nf_x86_64_sse.nb300nf_dosingle
        jmp    _nb_kernel300nf_x86_64_sse.nb300nf_updateouterdata
_nb_kernel300nf_x86_64_sse.nb300nf_dosingle: 
        movq nb300nf_charge(%rbp),%rsi
        movq nb300nf_pos(%rbp),%rdi
        movq  nb300nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm6,%xmm6
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulps  nb300nf_iq(%rsp),%xmm6
        movaps %xmm6,nb300nf_qq(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movaps nb300nf_ix(%rsp),%xmm4
        movaps nb300nf_iy(%rsp),%xmm5
        movaps nb300nf_iz(%rsp),%xmm6

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
        movaps nb300nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb300nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb300nf_tsc(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb300nf_VFtab(%rbp),%rsi
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
        movaps nb300nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV 
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addss  nb300nf_vctot(%rsp),%xmm5
        movss %xmm5,nb300nf_vctot(%rsp)

_nb_kernel300nf_x86_64_sse.nb300nf_updateouterdata: 
        ## get n from stack
        movl nb300nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb300nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb300nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb300nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb300nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300nf_x86_64_sse.nb300nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300nf_n(%rsp)
        jmp _nb_kernel300nf_x86_64_sse.nb300nf_outer
_nb_kernel300nf_x86_64_sse.nb300nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300nf_x86_64_sse.nb300nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300nf_x86_64_sse.nb300nf_threadloop
_nb_kernel300nf_x86_64_sse.nb300nf_end: 

        movl nb300nf_nouter(%rsp),%eax
        movl nb300nf_ninner(%rsp),%ebx
        movq nb300nf_outeriter(%rbp),%rcx
        movq nb300nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $248,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




