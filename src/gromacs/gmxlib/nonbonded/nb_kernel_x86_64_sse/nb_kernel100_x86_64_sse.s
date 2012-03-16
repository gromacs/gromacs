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







.globl nb_kernel100_x86_64_sse
.globl _nb_kernel100_x86_64_sse
nb_kernel100_x86_64_sse:        
_nb_kernel100_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb100_fshift, 16
.set nb100_gid, 24
.set nb100_pos, 32
.set nb100_faction, 40
.set nb100_charge, 48
.set nb100_p_facel, 56
.set nb100_argkrf, 64
.set nb100_argcrf, 72
.set nb100_Vc, 80
.set nb100_type, 88
.set nb100_p_ntype, 96
.set nb100_vdwparam, 104
.set nb100_Vvdw, 112
.set nb100_p_tabscale, 120
.set nb100_VFtab, 128
.set nb100_invsqrta, 136
.set nb100_dvda, 144
.set nb100_p_gbtabscale, 152
.set nb100_GBtab, 160
.set nb100_p_nthreads, 168
.set nb100_count, 176
.set nb100_mtx, 184
.set nb100_outeriter, 192
.set nb100_inneriter, 200
.set nb100_work, 208
        ## facel, krf,crf, tabscale, gbtabscale passed in xmm regs.
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100_ix, 0
.set nb100_iy, 16
.set nb100_iz, 32
.set nb100_iq, 48
.set nb100_dx, 64
.set nb100_dy, 80
.set nb100_dz, 96
.set nb100_vctot, 112
.set nb100_fix, 128
.set nb100_fiy, 144
.set nb100_fiz, 160
.set nb100_half, 176
.set nb100_three, 192
.set nb100_innerjjnr, 208
.set nb100_iinr, 216
.set nb100_jindex, 224
.set nb100_jjnr, 232
.set nb100_shift, 240
.set nb100_shiftvec, 248
.set nb100_facel, 256
.set nb100_is3, 264
.set nb100_ii3, 268
.set nb100_innerk, 272
.set nb100_n, 276
.set nb100_nn1, 280
.set nb100_nouter, 284
.set nb100_ninner, 288
.set nb100_nri, 292

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
    subq $312,%rsp      # # local variable stack space (n*16+8)                                                         

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb100_nouter(%rsp)
        movl %eax,nb100_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb100_nri(%rsp)
        movq %rsi,nb100_iinr(%rsp)
        movq %rdx,nb100_jindex(%rsp)
        movq %rcx,nb100_jjnr(%rsp)
        movq %r8,nb100_shift(%rsp)
        movq %r9,nb100_shiftvec(%rsp)
        movq nb100_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb100_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb100_half(%rsp)
        movss nb100_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb100_half(%rsp)
        movaps %xmm3,nb100_three(%rsp)

_nb_kernel100_x86_64_sse.nb100_threadloop: 
    movq  nb100_count(%rbp),%rsi            ## pointer to sync counter
    movl  (%rsi),%eax
_nb_kernel100_x86_64_sse.nb100_spinlock: 
    movl  %eax,%ebx                         ## ebx=*count=nn0
    addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
    pause                                   ## -> better p4 performance
    jnz _nb_kernel100_x86_64_sse.nb100_spinlock

    ## if(nn1>nri) nn1=nri
    movl nb100_nri(%rsp),%ecx
    movl %ecx,%edx
    subl %ebx,%ecx
    cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
    ## Cleared the spinlock if we got here.
    ## eax contains nn0, ebx contains nn1.
    movl %eax,nb100_n(%rsp)
    movl %ebx,nb100_nn1(%rsp)
    subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                                                  ## copy n to esi
    jg  _nb_kernel100_x86_64_sse.nb100_outerstart
    jmp _nb_kernel100_x86_64_sse.nb100_end

_nb_kernel100_x86_64_sse.nb100_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100_nouter(%rsp),%ebx
        movl %ebx,nb100_nouter(%rsp)

_nb_kernel100_x86_64_sse.nb100_outer: 
        movq  nb100_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                                ## ebx=shift[n] 

        lea    (%rbx,%rbx,2),%rbx              ## rbx=3*is 
        movl   %ebx,nb100_is3(%rsp)                     ## store is3 

        movq    nb100_shiftvec(%rsp),%rax       ## eax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm10
        movss 4(%rax,%rbx,4),%xmm11
        movss 8(%rax,%rbx,4),%xmm12

        movq  nb100_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]
        movl  (%rcx,%rsi,4),%ebx                           ## ebx =ii 

        movq  nb100_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb100_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx                ## rbx = 3*ii=ii3 


    movq  nb100_pos(%rbp),%rax      ## rax = base of pos[]  
        addss (%rax,%rbx,4),%xmm10
        addss 4(%rax,%rbx,4),%xmm11
        addss 8(%rax,%rbx,4),%xmm12

        movaps %xmm3,nb100_iq(%rsp)

        shufps $0,%xmm10,%xmm10
        shufps $0,%xmm11,%xmm11
        shufps $0,%xmm12,%xmm12

    movaps %xmm10,nb100_ix(%rsp)
    movaps %xmm11,nb100_iy(%rsp)
    movaps %xmm12,nb100_iz(%rsp)

        movl  %ebx,nb100_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorps %xmm12,%xmm12
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15

        movq  nb100_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx                                  ## number of innerloop atoms 

        movq  nb100_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb100_innerjjnr(%rsp)     ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb100_ninner(%rsp),%ecx
        movl  %ecx,nb100_ninner(%rsp)
        addl  $0,%edx ## to check sign
        movl  %edx,nb100_innerk(%rsp)      ## number of innerloop atoms 

        jge   _nb_kernel100_x86_64_sse.nb100_unroll_loop
        jmp   _nb_kernel100_x86_64_sse.nb100_finish_inner

_nb_kernel100_x86_64_sse.nb100_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb100_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb100_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%rax,%rax,2),%r8     ## j3
        lea  (%rbx,%rbx,2),%r9

        lea  (%rcx,%rcx,2),%r10
        lea  (%rdx,%rdx,2),%r11

        movq nb100_pos(%rbp),%rdi
        ## load coordinates    
        movlps (%rdi,%r8,4),%xmm1       ## x1 y1 - - 
        movlps (%rdi,%r9,4),%xmm2       ## x2 y2 - - 
        movlps (%rdi,%r10,4),%xmm3      ## x3 y3 - -
        movlps (%rdi,%r11,4),%xmm4      ## x4 y4 - -

        movss 8(%rdi,%r8,4),%xmm5       ## z1 - - - 
        movss 8(%rdi,%r9,4),%xmm6       ## z2 - - - 
        movss 8(%rdi,%r10,4),%xmm7      ## z3 - - - 
        movss 8(%rdi,%r11,4),%xmm8      ## z4 - - - 

    unpcklps %xmm3,%xmm1 ## x1 x3 y1 y3
    unpcklps %xmm4,%xmm2 ## x2 x4 y2 y4
    unpcklps %xmm7,%xmm5 ## z1 z3 -  -
    unpcklps %xmm8,%xmm6 ## z2 z4 -  -

    movaps %xmm1,%xmm3

        movq nb100_charge(%rbp),%rsi
    unpcklps %xmm2,%xmm1 ## x1 x2 x3 x4
    unpckhps %xmm2,%xmm3 ## y1 y2 y3 y4
    unpcklps %xmm6,%xmm5 ## z1 z2 z3 z4

        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rcx,4),%xmm2
        movss (%rsi,%rbx,4),%xmm7
        movss (%rsi,%rdx,4),%xmm8

        ## calc dr  
        subps nb100_ix(%rsp),%xmm1
        subps nb100_iy(%rsp),%xmm3
        subps nb100_iz(%rsp),%xmm5

        ## store dr in xmm9-xmm11
    movaps %xmm1,%xmm9
    movaps %xmm3,%xmm10
    movaps %xmm5,%xmm11

        ## square it 
        mulps %xmm1,%xmm1
        mulps %xmm3,%xmm3
        mulps %xmm5,%xmm5
        addps %xmm3,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

    unpcklps %xmm2,%xmm0 ## jqa jqc - -
    unpcklps %xmm8,%xmm7 ## jqb jqd - -

    ## calculate rinv=1/sqrt(rsq)
        rsqrtps %xmm1,%xmm5
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
    unpcklps %xmm7,%xmm0 ## jqa jqb jqc jqd
        movaps nb100_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm4   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4
    mulps nb100_iq(%rsp),%xmm0
        mulps nb100_half(%rsp),%xmm4
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    ## calculate coulomb interaction, xmm0=qq
        mulps  %xmm1,%xmm0      ## xmm0=vcoul 
        mulps  %xmm0,%xmm4      ## xmm4=fscal 

    ## add potential to vctot (sum in xmm12)
        addps  %xmm0,%xmm12

        movq nb100_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r10,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r11,4),%xmm1 ## x3 y3 x4 y4

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

        movlps %xmm0,(%rsi,%r8,4)
        movlps %xmm1,(%rsi,%r10,4)
        movhps %xmm0,(%rsi,%r9,4)
        movhps %xmm1,(%rsi,%r11,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm5 ## fjz2 - - -
    movhlps %xmm11,%xmm4     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm3  ## fjz4 - - -

        addss  8(%rsi,%r8,4),%xmm11
        addss  8(%rsi,%r9,4),%xmm5
        addss  8(%rsi,%r10,4),%xmm4
        addss  8(%rsi,%r11,4),%xmm3
        movss  %xmm11,8(%rsi,%r8,4)
        movss  %xmm5,8(%rsi,%r9,4)
        movss  %xmm4,8(%rsi,%r10,4)
        movss  %xmm3,8(%rsi,%r11,4)

        ## should we do one more iteration? 
        subl $4,nb100_innerk(%rsp)
        jl    _nb_kernel100_x86_64_sse.nb100_finish_inner
        jmp   _nb_kernel100_x86_64_sse.nb100_unroll_loop
_nb_kernel100_x86_64_sse.nb100_finish_inner: 
    ## check if at least two particles remain 
    addl $4,nb100_innerk(%rsp)
    movl  nb100_innerk(%rsp),%edx
    andl  $2,%edx
    jnz   _nb_kernel100_x86_64_sse.nb100_dopair
    jmp   _nb_kernel100_x86_64_sse.nb100_checksingle
_nb_kernel100_x86_64_sse.nb100_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb100_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb100_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb100_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rbx,4),%xmm1

    unpcklps %xmm1,%xmm0 ## jqa jqb - -
        mulps nb100_iq(%rsp),%xmm0      ##qq

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        movq nb100_pos(%rbp),%rdi
        ## load coordinates    
        movlps (%rdi,%rax,4),%xmm4      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm5      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm6      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr  
        subps nb100_ix(%rsp),%xmm4
        subps nb100_iy(%rsp),%xmm5
        subps nb100_iz(%rsp),%xmm6

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
        movaps nb100_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu    
    subps %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps nb100_half(%rsp),%xmm1
        movaps %xmm1,%xmm4
        mulps  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    xorps %xmm6,%xmm6

    ## calculate coulomb interaction, xmm0=qq
        mulps  %xmm1,%xmm0      ## xmm0=vcoul 
        mulps  %xmm0,%xmm4      ## xmm4=fscal 

    movlhps %xmm6,%xmm0

    ## add potential to vctot (sum in xmm12)
        addps  %xmm0,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm6,%xmm9
    movlhps %xmm6,%xmm10
    movlhps %xmm6,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb100_faction(%rbp),%rsi
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

_nb_kernel100_x86_64_sse.nb100_checksingle:     
    movl  nb100_innerk(%rsp),%edx
    andl  $1,%edx
    jnz    _nb_kernel100_x86_64_sse.nb100_dosingle
    jmp    _nb_kernel100_x86_64_sse.nb100_updateouterdata

_nb_kernel100_x86_64_sse.nb100_dosingle: 
    movq nb100_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb100_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0       ## jq
        mulss nb100_iq(%rsp),%xmm0      ## qq

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        movq nb100_pos(%rbp),%rdi
        movss (%rdi,%rax,4),%xmm4           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm5       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm6       ## 13 - - - 

        ## calc dr  
        subss nb100_ix(%rsp),%xmm4
        subss nb100_iy(%rsp),%xmm5
        subss nb100_iz(%rsp),%xmm6

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
        movaps nb100_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu    
    subss %xmm5,%xmm1   ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss nb100_half(%rsp),%xmm1
        movaps %xmm1,%xmm4
        mulss  %xmm4,%xmm4
    ## xmm1=rinv
    ## xmm4=rinvsq 

    ## calculate coulomb interaction, xmm0=qq
        mulss  %xmm1,%xmm0      ## xmm0=vcoul 
        mulss  %xmm0,%xmm4      ## xmm4=fscal 

    ## add potential to vctot (sum in xmm12)
        addss  %xmm0,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15
        movq nb100_faction(%rbp),%rsi

    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel100_x86_64_sse.nb100_updateouterdata: 

        movl  nb100_ii3(%rsp),%ecx
        movq  nb100_fshift(%rbp),%rsi
        movl  nb100_is3(%rsp),%edx
        movq  nb100_faction(%rbp),%rdi

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
        movl nb100_n(%rsp),%esi
    ## get group index for i particle 
    movq  nb100_gid(%rbp),%rdx          ## base of gid[]
    movl  (%rdx,%rsi,4),%edx            ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb100_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)


        ## finish if last 
        movl nb100_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100_x86_64_sse.nb100_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100_n(%rsp)
        jmp _nb_kernel100_x86_64_sse.nb100_outer
_nb_kernel100_x86_64_sse.nb100_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100_nri(%rsp),%ecx
        ## n is already loaded in esi
        subl  %esi,%ecx
        jz _nb_kernel100_x86_64_sse.nb100_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100_x86_64_sse.nb100_threadloop
_nb_kernel100_x86_64_sse.nb100_end: 

        movl nb100_nouter(%rsp),%eax
        movl nb100_ninner(%rsp),%ebx
        movq nb100_outeriter(%rbp),%rcx
        movq nb100_inneriter(%rbp),%rdx
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







.globl nb_kernel100nf_x86_64_sse
.globl _nb_kernel100nf_x86_64_sse
nb_kernel100nf_x86_64_sse:      
_nb_kernel100nf_x86_64_sse:     
.set nb100nf_fshift, 16
.set nb100nf_gid, 24
.set nb100nf_pos, 32
.set nb100nf_faction, 40
.set nb100nf_charge, 48
.set nb100nf_p_facel, 56
.set nb100nf_argkrf, 64
.set nb100nf_argcrf, 72
.set nb100nf_Vc, 80
.set nb100nf_type, 88
.set nb100nf_p_ntype, 96
.set nb100nf_vdwparam, 104
.set nb100nf_Vvdw, 112
.set nb100nf_p_tabscale, 120
.set nb100nf_VFtab, 128
.set nb100nf_invsqrta, 136
.set nb100nf_dvda, 144
.set nb100nf_p_gbtabscale, 152
.set nb100nf_GBtab, 160
.set nb100nf_p_nthreads, 168
.set nb100nf_count, 176
.set nb100nf_mtx, 184
.set nb100nf_outeriter, 192
.set nb100nf_inneriter, 200
.set nb100nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb100nf_ix, 0
.set nb100nf_iy, 16
.set nb100nf_iz, 32
.set nb100nf_iq, 48
.set nb100nf_vctot, 64
.set nb100nf_half, 80
.set nb100nf_three, 96
.set nb100nf_nri, 112
.set nb100nf_iinr, 120
.set nb100nf_jindex, 128
.set nb100nf_jjnr, 136
.set nb100nf_shift, 144
.set nb100nf_shiftvec, 152
.set nb100nf_facel, 160
.set nb100nf_innerjjnr, 168
.set nb100nf_is3, 176
.set nb100nf_ii3, 180
.set nb100nf_innerk, 184
.set nb100nf_n, 188
.set nb100nf_nn1, 192
.set nb100nf_nouter, 196
.set nb100nf_ninner, 200

        push %rbp
        movq %rsp,%rbp
        push %rbx

        subq $216,%rsp          # # local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb100nf_nouter(%rsp)
        movl %eax,nb100nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb100nf_nri(%rsp)
        movq %rsi,nb100nf_iinr(%rsp)
        movq %rdx,nb100nf_jindex(%rsp)
        movq %rcx,nb100nf_jjnr(%rsp)
        movq %r8,nb100nf_shift(%rsp)
        movq %r9,nb100nf_shiftvec(%rsp)
        movq nb100nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb100nf_facel(%rsp)




        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb100nf_half(%rsp)
        movss nb100nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb100nf_half(%rsp)
        movaps %xmm3,nb100nf_three(%rsp)

_nb_kernel100nf_x86_64_sse.nb100nf_threadloop: 
        movq  nb100nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel100nf_x86_64_sse.nb100nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel100nf_x86_64_sse.nb100nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb100nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb100nf_n(%rsp)
        movl %ebx,nb100nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel100nf_x86_64_sse.nb100nf_outerstart
        jmp _nb_kernel100nf_x86_64_sse.nb100nf_end

_nb_kernel100nf_x86_64_sse.nb100nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb100nf_nouter(%rsp),%ebx
        movl %ebx,nb100nf_nouter(%rsp)

_nb_kernel100nf_x86_64_sse.nb100nf_outer: 
        movq  nb100nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb100nf_is3(%rsp)            ## store is3 

        movq  nb100nf_shiftvec(%rsp),%rax       ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb100nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movq  nb100nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb100nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb100nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb100nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb100nf_ix(%rsp)
        movaps %xmm1,nb100nf_iy(%rsp)
        movaps %xmm2,nb100nf_iz(%rsp)

        movl  %ebx,nb100nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb100nf_vctot(%rsp)

        movq  nb100nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb100nf_pos(%rbp),%rsi

        movq  nb100nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb100nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb100nf_ninner(%rsp),%ecx
        movl  %ecx,nb100nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb100nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel100nf_x86_64_sse.nb100nf_unroll_loop
        jmp   _nb_kernel100nf_x86_64_sse.nb100nf_finish_inner
_nb_kernel100nf_x86_64_sse.nb100nf_unroll_loop: 
        ## quad-unrolled innerloop here 
        movq  nb100nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb100nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb100nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb100nf_iq(%rsp),%xmm5
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000           
        movq nb100nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        mulps %xmm5,%xmm3
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   

        movlps (%rsi,%rax,4),%xmm4      ## x1 y1 - - 
        movlps (%rsi,%rcx,4),%xmm5      ## x3 y3 - - 
        movss 8(%rsi,%rax,4),%xmm2      ## z1 -  - - 
        movss 8(%rsi,%rcx,4),%xmm6      ## z3 -  - - 

        movhps (%rsi,%rbx,4),%xmm4      ## x1 y1 x2 y2 
        movhps (%rsi,%rdx,4),%xmm5      ## x3 y3 x4 y4 

        movss 8(%rsi,%rbx,4),%xmm0      ## z2 - - - 
        movss 8(%rsi,%rdx,4),%xmm1      ## z4 - - - 

        shufps $0,%xmm0,%xmm2          ## z1 z1 z2 z2 
        shufps $0,%xmm1,%xmm6          ## z3 z3 z4 z4 

        movaps %xmm4,%xmm0              ## x1 y1 x2 y2  
        movaps %xmm4,%xmm1              ## x1 y1 x2 y2 

        shufps $136,%xmm6,%xmm2 ## 10001000     ;# z1 z2 z3 z4 

        shufps $136,%xmm5,%xmm0 ## 10001000     ;# x1 x2 x3 x4 
        shufps $221,%xmm5,%xmm1 ## 11011101     ;# y1 y2 y3 y4          

        ## move nb100nf_ix-iz to xmm4-xmm6 
        movaps nb100nf_ix(%rsp),%xmm4
        movaps nb100nf_iy(%rsp),%xmm5
        movaps nb100nf_iz(%rsp),%xmm6

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
        movaps nb100nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        movaps nb100nf_vctot(%rsp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb100nf_innerk(%rsp)
        jl    _nb_kernel100nf_x86_64_sse.nb100nf_finish_inner
        jmp   _nb_kernel100nf_x86_64_sse.nb100nf_unroll_loop
_nb_kernel100nf_x86_64_sse.nb100nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb100nf_innerk(%rsp)
        movl  nb100nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel100nf_x86_64_sse.nb100nf_dopair
        jmp   _nb_kernel100nf_x86_64_sse.nb100nf_checksingle
_nb_kernel100nf_x86_64_sse.nb100nf_dopair: 
        movq nb100nf_charge(%rbp),%rsi
        movq nb100nf_pos(%rbp),%rdi
        movq  nb100nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb100nf_innerjjnr(%rsp)

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $0,%xmm6,%xmm3
        shufps $8,%xmm3,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        mulps  nb100nf_iq(%rsp),%xmm3
        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        ## move nb100nf_ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb100nf_ix(%rsp),%xmm4
        movaps nb100nf_iy(%rsp),%xmm5
        movaps nb100nf_iz(%rsp),%xmm6

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
        movaps nb100nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps nb100nf_vctot(%rsp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%rsp)
_nb_kernel100nf_x86_64_sse.nb100nf_checksingle: 
        movl  nb100nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel100nf_x86_64_sse.nb100nf_dosingle
        jmp    _nb_kernel100nf_x86_64_sse.nb100nf_updateouterdata
_nb_kernel100nf_x86_64_sse.nb100nf_dosingle: 
        movq nb100nf_charge(%rbp),%rsi
        movq nb100nf_pos(%rbp),%rdi
        movq  nb100nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        movss (%rsi,%rax,4),%xmm3       ## xmm3(0) has the charge       

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        mulps  nb100nf_iq(%rsp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb100nf_ix(%rsp),%xmm4
        movaps nb100nf_iy(%rsp),%xmm5
        movaps nb100nf_iz(%rsp),%xmm6

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
        movaps nb100nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb100nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps nb100nf_vctot(%rsp),%xmm5
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addss  %xmm3,%xmm5
        movaps %xmm5,nb100nf_vctot(%rsp)

_nb_kernel100nf_x86_64_sse.nb100nf_updateouterdata: 
        ## get n from stack
        movl nb100nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb100nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb100nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb100nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb100nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel100nf_x86_64_sse.nb100nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb100nf_n(%rsp)
        jmp _nb_kernel100nf_x86_64_sse.nb100nf_outer
_nb_kernel100nf_x86_64_sse.nb100nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb100nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel100nf_x86_64_sse.nb100nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel100nf_x86_64_sse.nb100nf_threadloop
_nb_kernel100nf_x86_64_sse.nb100nf_end: 

        movl nb100nf_nouter(%rsp),%eax
        movl nb100nf_ninner(%rsp),%ebx
        movq nb100nf_outeriter(%rbp),%rcx
        movq nb100nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $216,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret

