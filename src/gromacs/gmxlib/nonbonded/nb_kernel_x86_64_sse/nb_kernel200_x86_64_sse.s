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






.globl nb_kernel200_x86_64_sse
.globl _nb_kernel200_x86_64_sse
nb_kernel200_x86_64_sse:        
_nb_kernel200_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb200_fshift, 16
.set nb200_gid, 24
.set nb200_pos, 32
.set nb200_faction, 40
.set nb200_charge, 48
.set nb200_p_facel, 56
.set nb200_argkrf, 64
.set nb200_argcrf, 72
.set nb200_Vc, 80
.set nb200_type, 88
.set nb200_p_ntype, 96
.set nb200_vdwparam, 104
.set nb200_Vvdw, 112
.set nb200_p_tabscale, 120
.set nb200_VFtab, 128
.set nb200_invsqrta, 136
.set nb200_dvda, 144
.set nb200_p_gbtabscale, 152
.set nb200_GBtab, 160
.set nb200_p_nthreads, 168
.set nb200_count, 176
.set nb200_mtx, 184
.set nb200_outeriter, 192
.set nb200_inneriter, 200
.set nb208_work, 200
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200_ix, 0
.set nb200_iy, 16
.set nb200_iz, 32
.set nb200_iq, 48
.set nb200_dx, 64
.set nb200_dy, 80
.set nb200_dz, 96
.set nb200_vctot, 112
.set nb200_fix, 128
.set nb200_fiy, 144
.set nb200_fiz, 160
.set nb200_half, 176
.set nb200_three, 192
.set nb200_two, 208
.set nb200_krf, 224
.set nb200_crf, 240
.set nb200_innerjjnr, 256
.set nb200_nri, 264
.set nb200_iinr, 272
.set nb200_jindex, 280
.set nb200_jjnr, 288
.set nb200_shift, 296
.set nb200_shiftvec, 304
.set nb200_facel, 312
.set nb200_is3, 320
.set nb200_ii3, 324
.set nb200_innerk, 328
.set nb200_n, 332
.set nb200_nn1, 336
.set nb200_nouter, 340
.set nb200_ninner, 344

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $360,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb200_nouter(%rsp)
        movl %eax,nb200_ninner(%rsp)




        movl (%rdi),%edi
        movl %edi,nb200_nri(%rsp)
        movq %rsi,nb200_iinr(%rsp)
        movq %rdx,nb200_jindex(%rsp)
        movq %rcx,nb200_jjnr(%rsp)
        movq %r8,nb200_shift(%rsp)
        movq %r9,nb200_shiftvec(%rsp)
        movq nb200_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb200_facel(%rsp)


        movq nb200_argkrf(%rbp),%rsi
        movq nb200_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb200_krf(%rsp)
        movaps %xmm2,nb200_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb200_half(%rsp)
        movss nb200_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb200_half(%rsp)
        movaps %xmm2,nb200_two(%rsp)
        movaps %xmm3,nb200_three(%rsp)

_nb_kernel200_x86_64_sse.nb200_threadloop: 
        movq  nb200_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel200_x86_64_sse.nb200_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200_x86_64_sse.nb200_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200_n(%rsp)
        movl %ebx,nb200_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200_x86_64_sse.nb200_outerstart
        jmp _nb_kernel200_x86_64_sse.nb200_end

        ## assume we have at least one i particle - start directly      
_nb_kernel200_x86_64_sse.nb200_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200_nouter(%rsp),%ebx
        movl %ebx,nb200_nouter(%rsp)

_nb_kernel200_x86_64_sse.nb200_outer: 
        movq  nb200_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb200_is3(%rsp)      ## store is3 

        movq  nb200_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb200_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb200_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb200_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb200_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb200_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb200_ix(%rsp)
        movaps %xmm1,nb200_iy(%rsp)
        movaps %xmm2,nb200_iz(%rsp)

        movl  %ebx,nb200_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorps %xmm12,%xmm12
        movaps %xmm12,%xmm13
        movaps %xmm12,%xmm14
        movaps %xmm12,%xmm15

        movq  nb200_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb200_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb200_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb200_ninner(%rsp),%ecx
        movl  %ecx,nb200_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb200_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel200_x86_64_sse.nb200_unroll_loop
        jmp   _nb_kernel200_x86_64_sse.nb200_finish_inner
_nb_kernel200_x86_64_sse.nb200_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb200_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d           ## eax-edx=jnr1-4 
        addq $16,nb200_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        lea  (%r8,%r8,2),%rax     ## j3
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx
        lea  (%r11,%r11,2),%rdx

        movq nb200_pos(%rbp),%rdi
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
        movq nb200_charge(%rbp),%rsi

        ## calc dr  
        subps nb200_ix(%rsp),%xmm1
        subps nb200_iy(%rsp),%xmm3
        subps nb200_iz(%rsp),%xmm5

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
        addps %xmm3,%xmm1
        addps %xmm5,%xmm1
        ## rsq in xmm1

        movaps nb200_krf(%rsp),%xmm7

    unpcklps %xmm2,%xmm0 ## jqa jqc - -
    unpcklps %xmm8,%xmm6 ## jqb jqd - -

        rsqrtps %xmm1,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200_three(%rsp),%xmm4
        mulps %xmm1,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%rsp),%xmm3
        mulps  %xmm1,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm4       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm4

    unpcklps %xmm6,%xmm0 ## jqa jqb jqc jqd
        mulps nb200_iq(%rsp),%xmm0      ##qq

        mulps %xmm4,%xmm3       ## xmm3=rinv 
        movaps %xmm3,%xmm1
        mulps  %xmm1,%xmm1      ## xmm1=rinvsq 
        movaps %xmm3,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subps  nb200_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulps  %xmm0,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addps  %xmm7,%xmm7

        subps  %xmm7,%xmm3
        mulps  %xmm3,%xmm0
        mulps  %xmm0,%xmm1      ## xmm1=total fscal 

    ## add potential to vctot (sum in xmm12)
        addps  %xmm6,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm1,%xmm9
        mulps  %xmm1,%xmm10
        mulps  %xmm1,%xmm11

        movq nb200_faction(%rbp),%rsi
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
        subl $4,nb200_innerk(%rsp)
        jl    _nb_kernel200_x86_64_sse.nb200_finish_inner
        jmp   _nb_kernel200_x86_64_sse.nb200_unroll_loop
_nb_kernel200_x86_64_sse.nb200_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb200_innerk(%rsp)
        movl  nb200_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel200_x86_64_sse.nb200_dopair
        jmp   _nb_kernel200_x86_64_sse.nb200_checksingle
_nb_kernel200_x86_64_sse.nb200_dopair: 
        ## twice-unrolled innerloop here 
        movq  nb200_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb200_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb200_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rbx,4),%xmm1

    unpcklps %xmm1,%xmm0 ## jqa jqb - -
        mulps nb200_iq(%rsp),%xmm0      ##qq

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load coordinates
        movq nb200_pos(%rbp),%rdi

        movlps (%rdi,%rax,4),%xmm4      ## x1 y1 - - 
        movlps (%rdi,%rbx,4),%xmm5      ## x2 y2 - - 

        movss 8(%rdi,%rax,4),%xmm6      ## z1 - - - 
        movss 8(%rdi,%rbx,4),%xmm7      ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr  
        subps nb200_ix(%rsp),%xmm4
        subps nb200_iy(%rsp),%xmm5
        subps nb200_iz(%rsp),%xmm6

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

        movaps nb200_krf(%rsp),%xmm7

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%rsp),%xmm3
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm3       ## xmm3=rinv 
        movaps %xmm3,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm3,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subps  nb200_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulps  %xmm0,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addps  %xmm7,%xmm7

        subps  %xmm7,%xmm3
        mulps  %xmm3,%xmm0
        mulps  %xmm0,%xmm4      ## xmm4=total fscal 

    xorps  %xmm5,%xmm5
    movlhps %xmm5,%xmm6

    ## add potential to vctot (sum in xmm12)
        addps  %xmm6,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulps  %xmm4,%xmm9
        mulps  %xmm4,%xmm10
        mulps  %xmm4,%xmm11

    movlhps %xmm5,%xmm9
    movlhps %xmm5,%xmm10
    movlhps %xmm5,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb200_faction(%rbp),%rsi
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

_nb_kernel200_x86_64_sse.nb200_checksingle:     
        movl  nb200_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200_x86_64_sse.nb200_dosingle
        jmp    _nb_kernel200_x86_64_sse.nb200_updateouterdata
_nb_kernel200_x86_64_sse.nb200_dosingle: 
    movq nb200_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb200_charge(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0       ## jq
        mulss nb200_iq(%rsp),%xmm0      ## qq

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 

        movq nb200_pos(%rbp),%rdi
        movss (%rdi,%rax,4),%xmm4           ## x1 - - - 
        movss 4(%rdi,%rax,4),%xmm5       ## y2 - - - 
        movss 8(%rdi,%rax,4),%xmm6       ## 13 - - - 

        ## calc dr  
        subss nb200_ix(%rsp),%xmm4
        subss nb200_iy(%rsp),%xmm5
        subss nb200_iz(%rsp),%xmm6

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

        movaps nb200_krf(%rsp),%xmm7

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb200_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200_half(%rsp),%xmm3
        mulss  %xmm4,%xmm7      ## xmm7=krsq 
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm3       ## xmm3=rinv 
        movaps %xmm3,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps %xmm3,%xmm6
        addss  %xmm7,%xmm6      ## xmm6=rinv+ krsq 

        subss  nb200_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 

        mulss  %xmm0,%xmm6      ## xmm6=vcoul=qq*(rinv+krsq-crf) 
        addss  %xmm7,%xmm7

        subss  %xmm7,%xmm3
        mulss  %xmm3,%xmm0
        mulss  %xmm0,%xmm4      ## xmm4=total fscal 

    ## add potential to vctot (sum in xmm12)
        addss  %xmm6,%xmm12

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulss  %xmm4,%xmm9
        mulss  %xmm4,%xmm10
        mulss  %xmm4,%xmm11

        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb200_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%rax,4),%xmm9
    addss  4(%rsi,%rax,4),%xmm10
    addss  8(%rsi,%rax,4),%xmm11
    movss  %xmm9,(%rsi,%rax,4)
    movss  %xmm10,4(%rsi,%rax,4)
    movss  %xmm11,8(%rsi,%rax,4)

_nb_kernel200_x86_64_sse.nb200_updateouterdata: 
        movl  nb200_ii3(%rsp),%ecx
        movq  nb200_faction(%rbp),%rdi
        movq  nb200_fshift(%rbp),%rsi
        movl  nb200_is3(%rsp),%edx

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
        movl nb200_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb200_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb200_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## finish if last 
        movl nb200_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200_x86_64_sse.nb200_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200_n(%rsp)
        jmp _nb_kernel200_x86_64_sse.nb200_outer
_nb_kernel200_x86_64_sse.nb200_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200_x86_64_sse.nb200_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200_x86_64_sse.nb200_threadloop
_nb_kernel200_x86_64_sse.nb200_end: 

        movl nb200_nouter(%rsp),%eax
        movl nb200_ninner(%rsp),%ebx
        movq nb200_outeriter(%rbp),%rcx
        movq nb200_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $360,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel200nf_x86_64_sse
.globl _nb_kernel200nf_x86_64_sse
nb_kernel200nf_x86_64_sse:      
_nb_kernel200nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb200nf_fshift, 16
.set nb200nf_gid, 24
.set nb200nf_pos, 32
.set nb200nf_faction, 40
.set nb200nf_charge, 48
.set nb200nf_p_facel, 56
.set nb200nf_argkrf, 64
.set nb200nf_argcrf, 72
.set nb200nf_Vc, 80
.set nb200nf_type, 88
.set nb200nf_p_ntype, 96
.set nb200nf_vdwparam, 104
.set nb200nf_Vvdw, 112
.set nb200nf_p_tabscale, 120
.set nb200nf_VFtab, 128
.set nb200nf_invsqrta, 136
.set nb200nf_dvda, 144
.set nb200nf_p_gbtabscale, 152
.set nb200nf_GBtab, 160
.set nb200nf_p_nthreads, 168
.set nb200nf_count, 176
.set nb200nf_mtx, 184
.set nb200nf_outeriter, 192
.set nb200nf_inneriter, 200
.set nb208nf_work, 200
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200nf_ix, 0
.set nb200nf_iy, 16
.set nb200nf_iz, 32
.set nb200nf_iq, 48
.set nb200nf_vctot, 64
.set nb200nf_half, 80
.set nb200nf_three, 96
.set nb200nf_krf, 112
.set nb200nf_crf, 128
.set nb200nf_is3, 144
.set nb200nf_ii3, 148
.set nb200nf_innerjjnr, 152
.set nb200nf_nri, 160
.set nb200nf_iinr, 168
.set nb200nf_jindex, 176
.set nb200nf_jjnr, 184
.set nb200nf_shift, 192
.set nb200nf_shiftvec, 200
.set nb200nf_facel, 208
.set nb200nf_innerk, 216
.set nb200nf_n, 220
.set nb200nf_nn1, 224
.set nb200nf_nouter, 228
.set nb200nf_ninner, 232

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
        movl %eax,nb200nf_nouter(%rsp)
        movl %eax,nb200nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb200nf_nri(%rsp)
        movq %rsi,nb200nf_iinr(%rsp)
        movq %rdx,nb200nf_jindex(%rsp)
        movq %rcx,nb200nf_jjnr(%rsp)
        movq %r8,nb200nf_shift(%rsp)
        movq %r9,nb200nf_shiftvec(%rsp)
        movq nb200nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb200nf_facel(%rsp)

        movq nb200nf_argkrf(%rbp),%rsi
        movq nb200nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb200nf_krf(%rsp)
        movaps %xmm2,nb200nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb200nf_half(%rsp)
        movss nb200nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb200nf_half(%rsp)
        movaps %xmm3,nb200nf_three(%rsp)

_nb_kernel200nf_x86_64_sse.nb200nf_threadloop: 
        movq  nb200nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel200nf_x86_64_sse.nb200nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200nf_x86_64_sse.nb200nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200nf_n(%rsp)
        movl %ebx,nb200nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200nf_x86_64_sse.nb200nf_outerstart
        jmp _nb_kernel200nf_x86_64_sse.nb200nf_end

_nb_kernel200nf_x86_64_sse.nb200nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200nf_nouter(%rsp),%ebx
        movl %ebx,nb200nf_nouter(%rsp)

_nb_kernel200nf_x86_64_sse.nb200nf_outer: 
        movq  nb200nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb200nf_is3(%rsp)            ## store is3 

        movq  nb200nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb200nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb200nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb200nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb200nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb200nf_iq(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb200nf_ix(%rsp)
        movaps %xmm1,nb200nf_iy(%rsp)
        movaps %xmm2,nb200nf_iz(%rsp)

        movl  %ebx,nb200nf_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb200nf_vctot(%rsp)

        movq  nb200nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb200nf_pos(%rbp),%rsi
        movq  nb200nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb200nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb200nf_ninner(%rsp),%ecx
        movl  %ecx,nb200nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb200nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel200nf_x86_64_sse.nb200nf_unroll_loop
        jmp   _nb_kernel200nf_x86_64_sse.nb200nf_finish_inner
_nb_kernel200nf_x86_64_sse.nb200nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb200nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb200nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb200nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        movaps nb200nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  

        movq nb200nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb200nf_ix(%rsp),%xmm4
        movaps nb200nf_iy(%rsp),%xmm5
        movaps nb200nf_iz(%rsp),%xmm6

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

        movaps nb200nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addps  nb200nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb200nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb200nf_innerk(%rsp)
        jl    _nb_kernel200nf_x86_64_sse.nb200nf_finish_inner
        jmp   _nb_kernel200nf_x86_64_sse.nb200nf_unroll_loop
_nb_kernel200nf_x86_64_sse.nb200nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb200nf_innerk(%rsp)
        movl  nb200nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel200nf_x86_64_sse.nb200nf_dopair
        jmp   _nb_kernel200nf_x86_64_sse.nb200nf_checksingle
_nb_kernel200nf_x86_64_sse.nb200nf_dopair: 
        movq nb200nf_charge(%rbp),%rsi

        movq  nb200nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb200nf_innerjjnr(%rsp)

        xorps %xmm3,%xmm3
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        shufps $12,%xmm6,%xmm3 ## 00001100 
        shufps $88,%xmm3,%xmm3 ## 01011000 ;# xmm3(0,1) has the charges  

        movq nb200nf_pos(%rbp),%rdi

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        mulps  nb200nf_iq(%rsp),%xmm3

        xorps  %xmm7,%xmm7
        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        ## move ix-iz to xmm4-xmm6 

        movaps nb200nf_ix(%rsp),%xmm4
        movaps nb200nf_iy(%rsp),%xmm5
        movaps nb200nf_iz(%rsp),%xmm6

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

        movaps nb200nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq-crf) 

        addps  nb200nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb200nf_vctot(%rsp)

_nb_kernel200nf_x86_64_sse.nb200nf_checksingle: 
        movl  nb200nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200nf_x86_64_sse.nb200nf_dosingle
        jmp    _nb_kernel200nf_x86_64_sse.nb200nf_updateouterdata
_nb_kernel200nf_x86_64_sse.nb200nf_dosingle: 
        movq nb200nf_charge(%rbp),%rsi
        movq nb200nf_pos(%rbp),%rdi
        movq  nb200nf_innerjjnr(%rsp),%rcx
        xorps %xmm3,%xmm3
        movl  (%rcx),%eax
        movss (%rsi,%rax,4),%xmm3       ## xmm3(0) has the charge 
        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        mulps  nb200nf_iq(%rsp),%xmm3

        xorps   %xmm7,%xmm7

        movaps nb200nf_ix(%rsp),%xmm4
        movaps nb200nf_iy(%rsp),%xmm5
        movaps nb200nf_iz(%rsp),%xmm6

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

        movaps nb200nf_krf(%rsp),%xmm7
        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb200nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb200nf_half(%rsp),%xmm0
        mulps  %xmm4,%xmm7      ## xmm7=krsq 
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        movaps %xmm0,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb200nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps  %xmm3,%xmm6      ## xmm6=vcoul 
        addss  nb200nf_vctot(%rsp),%xmm6
        movss %xmm6,nb200nf_vctot(%rsp)

_nb_kernel200nf_x86_64_sse.nb200nf_updateouterdata: 
        ## get n from stack
        movl nb200nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb200nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb200nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb200nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

       ## finish if last 
        movl nb200nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200nf_x86_64_sse.nb200nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200nf_n(%rsp)
        jmp _nb_kernel200nf_x86_64_sse.nb200nf_outer
_nb_kernel200nf_x86_64_sse.nb200nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200nf_x86_64_sse.nb200nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200nf_x86_64_sse.nb200nf_threadloop
_nb_kernel200nf_x86_64_sse.nb200nf_end: 

        movl nb200nf_nouter(%rsp),%eax
        movl nb200nf_ninner(%rsp),%ebx
        movq nb200nf_outeriter(%rbp),%rcx
        movq nb200nf_inneriter(%rbp),%rdx
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


