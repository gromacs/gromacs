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






.globl nb_kernel400_x86_64_sse
.globl _nb_kernel400_x86_64_sse
nb_kernel400_x86_64_sse:        
_nb_kernel400_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb400_fshift, 16
.set nb400_gid, 24
.set nb400_pos, 32
.set nb400_faction, 40
.set nb400_charge, 48
.set nb400_p_facel, 56
.set nb400_argkrf, 64
.set nb400_argcrf, 72
.set nb400_Vc, 80
.set nb400_type, 88
.set nb400_p_ntype, 96
.set nb400_vdwparam, 104
.set nb400_Vvdw, 112
.set nb400_p_tabscale, 120
.set nb400_VFtab, 128
.set nb400_invsqrta, 136
.set nb400_dvda, 144
.set nb400_p_gbtabscale, 152
.set nb400_GBtab, 160
.set nb400_p_nthreads, 168
.set nb400_count, 176
.set nb400_mtx, 184
.set nb400_outeriter, 192
.set nb400_inneriter, 200
.set nb400_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb400_ix, 0
.set nb400_iy, 16
.set nb400_iz, 32
.set nb400_iq, 48
.set nb400_dx, 64
.set nb400_dy, 80
.set nb400_dz, 96
.set nb400_two, 112
.set nb400_gbtsc, 128
.set nb400_qq, 144
.set nb400_r, 160
.set nb400_vctot, 176
.set nb400_fix, 192
.set nb400_fiy, 208
.set nb400_fiz, 224
.set nb400_half, 240
.set nb400_three, 256
.set nb400_isai, 272
.set nb400_isaprod, 288
.set nb400_dvdasum, 304
.set nb400_gbscale, 320
.set nb400_nri, 336
.set nb400_iinr, 344
.set nb400_jindex, 352
.set nb400_jjnr, 360
.set nb400_shift, 368
.set nb400_shiftvec, 376
.set nb400_facel, 384
.set nb400_innerjjnr, 392
.set nb400_is3, 400
.set nb400_ii3, 404
.set nb400_ii, 408
.set nb400_innerk, 412
.set nb400_n, 416
.set nb400_nn1, 420
.set nb400_nouter, 424
.set nb400_ninner, 428
.set nb400_jnra, 432
.set nb400_jnrb, 436
.set nb400_jnrc, 440
.set nb400_jnrd, 444

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $456,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb400_nouter(%rsp)
        movl %eax,nb400_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb400_nri(%rsp)
        movq %rsi,nb400_iinr(%rsp)
        movq %rdx,nb400_jindex(%rsp)
        movq %rcx,nb400_jjnr(%rsp)
        movq %r8,nb400_shift(%rsp)
        movq %r9,nb400_shiftvec(%rsp)
        movq nb400_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb400_facel(%rsp)

        movq nb400_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb400_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb400_half(%rsp)
        movss nb400_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb400_half(%rsp)
        movaps %xmm2,nb400_two(%rsp)
        movaps %xmm3,nb400_three(%rsp)

_nb_kernel400_x86_64_sse.nb400_threadloop: 
        movq  nb400_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel400_x86_64_sse.nb400_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400_x86_64_sse.nb400_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400_n(%rsp)
        movl %ebx,nb400_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400_x86_64_sse.nb400_outerstart
        jmp _nb_kernel400_x86_64_sse.nb400_end

_nb_kernel400_x86_64_sse.nb400_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400_nouter(%rsp),%ebx
        movl %ebx,nb400_nouter(%rsp)

_nb_kernel400_x86_64_sse.nb400_outer: 
        movq  nb400_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb400_is3(%rsp)      ## store is3 

        movq  nb400_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb400_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 
        movl  %ebx,nb400_ii(%rsp)

        movq  nb400_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb400_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3


        movq  nb400_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb400_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb400_iq(%rsp)
        movaps %xmm4,nb400_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb400_ix(%rsp)
        movaps %xmm1,nb400_iy(%rsp)
        movaps %xmm2,nb400_iz(%rsp)

        movl  %ebx,nb400_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb400_dvdasum(%rsp)
        movaps %xmm4,%xmm12
        movaps %xmm4,%xmm13
        movaps %xmm4,%xmm14
        movaps %xmm4,%xmm15

        movq  nb400_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb400_pos(%rbp),%rsi
        movq  nb400_faction(%rbp),%rdi
        movq  nb400_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb400_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb400_ninner(%rsp),%ecx
        movl  %ecx,nb400_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb400_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel400_x86_64_sse.nb400_unroll_loop
        jmp   _nb_kernel400_x86_64_sse.nb400_finish_inner
_nb_kernel400_x86_64_sse.nb400_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb400_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb400_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb400_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3
        lea  (%rbx,%rbx,2),%r9
        lea  (%rcx,%rcx,2),%r10
        lea  (%rdx,%rdx,2),%r11

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm4
        movlps (%rsi,%r10,4),%xmm5
        movss 8(%rsi,%r8,4),%xmm2
        movss 8(%rsi,%r10,4),%xmm6

        movhps (%rsi,%r9,4),%xmm4
        movhps (%rsi,%r11,4),%xmm5

        movss 8(%rsi,%r9,4),%xmm0
        movss 8(%rsi,%r11,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## calc dr 
        subps nb400_ix(%rsp),%xmm0
        subps nb400_iy(%rsp),%xmm1
        subps nb400_iz(%rsp),%xmm2

        ## store dr 
        movaps %xmm0,%xmm9
        movaps %xmm1,%xmm10
        movaps %xmm2,%xmm11

        ## square it 
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
    movaps %xmm0,%xmm4
        ## rsq in xmm4 

        ## load isaj
        movq nb400_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rcx,4),%xmm1
        movss (%rsi,%rbx,4),%xmm2
        movss (%rsi,%rdx,4),%xmm3
        movaps nb400_isai(%rsp),%xmm7
        shufps $0,%xmm2,%xmm0
    shufps $0,%xmm3,%xmm1
        shufps $136,%xmm1,%xmm0 ## 10001000 ;# all isaj in xmm3 
        mulps  %xmm0,%xmm7

        movaps %xmm7,nb400_isaprod(%rsp)
        movaps %xmm7,%xmm1
        mulps nb400_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb400_gbscale(%rsp)

        movq nb400_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm0
        movss (%rsi,%rcx,4),%xmm1
        movss (%rsi,%rbx,4),%xmm2
        movss (%rsi,%rdx,4),%xmm3

    mulps nb400_iq(%rsp),%xmm7
        shufps $0,%xmm2,%xmm0
        shufps $0,%xmm3,%xmm1
    shufps $136,%xmm1,%xmm0 ## 10001000 ;# all charges in xmm3  

        mulps  %xmm7,%xmm0
        movaps %xmm0,nb400_qq(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb400_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb400_r(%rsp)
        mulps nb400_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6

    ## multiply by 4
    pslld   $2,%xmm5

    ## move to integer registers
    movhlps %xmm5,%xmm7
    movd    %xmm5,%r12d
    movd    %xmm7,%r14d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm7,%xmm7
    movd    %xmm5,%r13d
    movd    %xmm7,%r15d

    ## calculate eps
    subps     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb400_GBtab(%rbp),%rsi

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

        shufps $136,%xmm8,%xmm6 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101
    ## table data ready in xmm4-xmm7

    mulps  %xmm1,%xmm7  ## Heps
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm1,%xmm7      ## Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        addps  %xmm7,%xmm7      ## two*Heps2 
        movaps nb400_qq(%rsp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

        movq nb400_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb400_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb400_r(%rsp),%xmm6
        addps  %xmm5,%xmm6

    ## increment vctot (sum in xmm12)
        addps  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## update dvdasum
    addps  nb400_dvdasum(%rsp),%xmm7
    movaps %xmm7,nb400_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        movaps  %xmm6,%xmm5
        movaps  %xmm7,%xmm4
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm4,%xmm4

        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        addss  (%rsi,%rcx,4),%xmm7
        addss  (%rsi,%rdx,4),%xmm4
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)
        movss  %xmm7,(%rsi,%rcx,4)
        movss  %xmm4,(%rsi,%rdx,4)

        xorps  %xmm4,%xmm4
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

        movq nb400_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r10,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r11,4),%xmm1 ## x3 y3 x4 y4

    mulps  %xmm4,%xmm9
    mulps  %xmm4,%xmm10
    mulps  %xmm4,%xmm11

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
    pshufd $1,%xmm11,%xmm10 ## fjz2 - - -
    movhlps %xmm11,%xmm9     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm8  ## fjz4 - - -

        addss  8(%rsi,%r8,4),%xmm11
        addss  8(%rsi,%r9,4),%xmm10
        addss  8(%rsi,%r10,4),%xmm9
        addss  8(%rsi,%r11,4),%xmm8
        movss  %xmm11,8(%rsi,%r8,4)
        movss  %xmm10,8(%rsi,%r9,4)
        movss  %xmm9,8(%rsi,%r10,4)
        movss  %xmm8,8(%rsi,%r11,4)

        ## should we do one more iteration? 
        subl $4,nb400_innerk(%rsp)
        jl    _nb_kernel400_x86_64_sse.nb400_finish_inner
        jmp   _nb_kernel400_x86_64_sse.nb400_unroll_loop
_nb_kernel400_x86_64_sse.nb400_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb400_innerk(%rsp)
        movl  nb400_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel400_x86_64_sse.nb400_dopair
        jmp   _nb_kernel400_x86_64_sse.nb400_checksingle
_nb_kernel400_x86_64_sse.nb400_dopair: 
        movq  nb400_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb400_innerjjnr(%rsp)

        ## load isaj
        movq nb400_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
    unpcklps %xmm6,%xmm3

        movaps nb400_isai(%rsp),%xmm2
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb400_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb400_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb400_gbscale(%rsp)

        movq nb400_charge(%rbp),%rsi     ## base of charge[] 

    mulps nb400_iq(%rsp),%xmm2
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
    unpcklps %xmm6,%xmm3

        mulps %xmm2,%xmm3
        movaps %xmm3,nb400_qq(%rsp)

        movq nb400_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3 
        lea  (%rbx,%rbx,2),%r9

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm4       ## x1 y1 - - 
        movlps (%rsi,%r9,4),%xmm5       ## x2 y2 - - 

        movss 8(%rsi,%r8,4),%xmm6       ## z1 - - - 
        movss 8(%rsi,%r9,4),%xmm7       ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr 
        subps nb400_ix(%rsp),%xmm4
        subps nb400_iy(%rsp),%xmm5
        subps nb400_iz(%rsp),%xmm6

        ## store dr 
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

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb400_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb400_r(%rsp)
        mulps nb400_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6

    ## multiply by 4
    pslld   $2,%xmm5

    ## move to integer registers
    movd    %xmm5,%r12d
    pshufd $1,%xmm5,%xmm5
    movd    %xmm5,%r13d

    ## calculate eps
    subps     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb400_GBtab(%rbp),%rsi

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

    mulps  %xmm1,%xmm7  ## Heps
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm1,%xmm7      ## Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        addps  %xmm7,%xmm7      ## two*Heps2 
        movaps nb400_qq(%rsp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

    ## zero upper part of vcoul 
    xorps %xmm2,%xmm2
    movlhps %xmm2,%xmm5

        movq nb400_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb400_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb400_r(%rsp),%xmm6
        addps  %xmm5,%xmm6

    ## increment vctot (sum in xmm12)
        addps  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## zero upper half of dvda
    movlhps %xmm2,%xmm7

    ## update dvdasum
    addps  nb400_dvdasum(%rsp),%xmm7
    movaps %xmm7,nb400_dvdasum(%rsp)

        ## update j atoms dvdaj
        movaps  %xmm6,%xmm5
        shufps $0x1,%xmm5,%xmm5

        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)

        xorps  %xmm4,%xmm4
        mulps %xmm0,%xmm3
        subps  %xmm3,%xmm4

    mulps  %xmm4,%xmm9
    mulps  %xmm4,%xmm10
    mulps  %xmm4,%xmm11

    movlhps %xmm2,%xmm9
    movlhps %xmm2,%xmm10
    movlhps %xmm2,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb400_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2

    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    addps    %xmm9,%xmm0

        movlps %xmm0,(%rsi,%r8,4)
        movhps %xmm0,(%rsi,%r9,4)

    ## z forces
    pshufd $1,%xmm11,%xmm8
    addss  8(%rsi,%r8,4),%xmm11
    addss  8(%rsi,%r9,4),%xmm8
    movss  %xmm11,8(%rsi,%r8,4)
    movss  %xmm8,8(%rsi,%r9,4)

_nb_kernel400_x86_64_sse.nb400_checksingle:     
        movl  nb400_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400_x86_64_sse.nb400_dosingle
        jmp    _nb_kernel400_x86_64_sse.nb400_updateouterdata
_nb_kernel400_x86_64_sse.nb400_dosingle: 
        movq  nb400_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb400_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm2
        mulss nb400_isai(%rsp),%xmm2
        movss %xmm2,nb400_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulss nb400_gbtsc(%rsp),%xmm1
        movss %xmm1,nb400_gbscale(%rsp)

        movq nb400_charge(%rbp),%rsi     ## base of charge[] 

    mulss nb400_iq(%rsp),%xmm2
        movss (%rsi,%rax,4),%xmm3
        mulss %xmm2,%xmm3
        movss %xmm3,nb400_qq(%rsp)

        movq nb400_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3=3*jnr

        ## move four coordinates to xmm0-xmm2   
        movss (%rsi,%r8,4),%xmm4
        movss 4(%rsi,%r8,4),%xmm5
        movss 8(%rsi,%r8,4),%xmm6

        ## calc dr 
        subss nb400_ix(%rsp),%xmm4
        subss nb400_iy(%rsp),%xmm5
        subss nb400_iz(%rsp),%xmm6

        ## store dr 
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

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb400_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 
        mulss %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb400_r(%rsp)
        mulss nb400_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttss2si %xmm4,%r12d

    ## convert back to float
    cvtsi2ss  %r12d,%xmm6

    ## multiply by 4
    shll $2,%r12d

    ## calculate eps
    subss     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb400_GBtab(%rbp),%rsi

    ## load table data
        movss (%rsi,%r12,4),%xmm4
        movss 4(%rsi,%r12,4),%xmm5
        movss 8(%rsi,%r12,4),%xmm6
        movss 12(%rsi,%r12,4),%xmm7
    ## table data ready in xmm4-xmm7

    mulss  %xmm1,%xmm7  ## Heps
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm1,%xmm7      ## Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        addss  %xmm7,%xmm7      ## two*Heps2 
        movss  nb400_qq(%rsp),%xmm3
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulss  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

        movq nb400_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulss nb400_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulss  nb400_r(%rsp),%xmm6
        addss  %xmm5,%xmm6

    ## increment vctot (sum in xmm12)
        addss  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subss  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## update dvdasum
    addss  nb400_dvdasum(%rsp),%xmm7
    movss %xmm7,nb400_dvdasum(%rsp)

        ## update j atoms dvdaj
        addss  (%rsi,%rax,4),%xmm6
        movss  %xmm6,(%rsi,%rax,4)

        xorps  %xmm4,%xmm4
        mulss %xmm0,%xmm3
        subss  %xmm3,%xmm4

    mulss  %xmm4,%xmm9
    mulss  %xmm4,%xmm10
    mulss  %xmm4,%xmm11

        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb400_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%r8,4),%xmm9
    addss  4(%rsi,%r8,4),%xmm10
    addss  8(%rsi,%r8,4),%xmm11
    movss  %xmm9,(%rsi,%r8,4)
    movss  %xmm10,4(%rsi,%r8,4)
    movss  %xmm11,8(%rsi,%r8,4)

_nb_kernel400_x86_64_sse.nb400_updateouterdata: 
        movl  nb400_ii3(%rsp),%ecx
        movq  nb400_faction(%rbp),%rdi
        movq  nb400_fshift(%rbp),%rsi
        movl  nb400_is3(%rsp),%edx

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
        movl nb400_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb400_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb400_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## accumulate dVda and update it 
        movaps nb400_dvdasum(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        movl nb400_ii(%rsp),%edx
        movq nb400_dvda(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb400_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400_x86_64_sse.nb400_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400_n(%rsp)
        jmp _nb_kernel400_x86_64_sse.nb400_outer
_nb_kernel400_x86_64_sse.nb400_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400_x86_64_sse.nb400_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400_x86_64_sse.nb400_threadloop
_nb_kernel400_x86_64_sse.nb400_end: 

        movl nb400_nouter(%rsp),%eax
        movl nb400_ninner(%rsp),%ebx
        movq nb400_outeriter(%rbp),%rcx
        movq nb400_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $456,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel400nf_x86_64_sse
.globl _nb_kernel400nf_x86_64_sse
nb_kernel400nf_x86_64_sse:      
_nb_kernel400nf_x86_64_sse:     
.set nb400nf_fshift, 16
.set nb400nf_gid, 24
.set nb400nf_pos, 32
.set nb400nf_faction, 40
.set nb400nf_charge, 48
.set nb400nf_p_facel, 56
.set nb400nf_argkrf, 64
.set nb400nf_argcrf, 72
.set nb400nf_Vc, 80
.set nb400nf_type, 88
.set nb400nf_p_ntype, 96
.set nb400nf_vdwparam, 104
.set nb400nf_Vvdw, 112
.set nb400nf_p_tabscale, 120
.set nb400nf_VFtab, 128
.set nb400nf_invsqrta, 136
.set nb400nf_dvda, 144
.set nb400nf_p_gbtabscale, 152
.set nb400nf_GBtab, 160
.set nb400nf_p_nthreads, 168
.set nb400nf_count, 176
.set nb400nf_mtx, 184
.set nb400nf_outeriter, 192
.set nb400nf_inneriter, 200
.set nb400nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb400nf_ix, 0
.set nb400nf_iy, 16
.set nb400nf_iz, 32
.set nb400nf_iq, 48
.set nb400nf_gbtsc, 64
.set nb400nf_qq, 80
.set nb400nf_vctot, 96
.set nb400nf_half, 112
.set nb400nf_three, 128
.set nb400nf_isai, 144
.set nb400nf_isaprod, 160
.set nb400nf_gbscale, 176
.set nb400nf_nri, 192
.set nb400nf_iinr, 200
.set nb400nf_jindex, 208
.set nb400nf_jjnr, 216
.set nb400nf_shift, 224
.set nb400nf_shiftvec, 232
.set nb400nf_facel, 240
.set nb400nf_innerjjnr, 248
.set nb400nf_is3, 256
.set nb400nf_ii3, 260
.set nb400nf_innerk, 264
.set nb400nf_n, 268
.set nb400nf_nn1, 272
.set nb400nf_nouter, 276
.set nb400nf_ninner, 280


        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $296,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb400nf_nouter(%rsp)
        movl %eax,nb400nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb400nf_nri(%rsp)
        movq %rsi,nb400nf_iinr(%rsp)
        movq %rdx,nb400nf_jindex(%rsp)
        movq %rcx,nb400nf_jjnr(%rsp)
        movq %r8,nb400nf_shift(%rsp)
        movq %r9,nb400nf_shiftvec(%rsp)
        movq nb400nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb400nf_facel(%rsp)

        movq nb400nf_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb400nf_gbtsc(%rsp)



        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb400nf_half(%rsp)
        movss nb400nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb400nf_half(%rsp)
        movaps %xmm3,nb400nf_three(%rsp)

_nb_kernel400nf_x86_64_sse.nb400nf_threadloop: 
        movq  nb400nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel400nf_x86_64_sse.nb400nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400nf_x86_64_sse.nb400nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400nf_n(%rsp)
        movl %ebx,nb400nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400nf_x86_64_sse.nb400nf_outerstart
        jmp _nb_kernel400nf_x86_64_sse.nb400nf_end

_nb_kernel400nf_x86_64_sse.nb400nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400nf_nouter(%rsp),%ebx
        movl %ebx,nb400nf_nouter(%rsp)

_nb_kernel400nf_x86_64_sse.nb400nf_outer: 
        movq  nb400nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb400nf_is3(%rsp)            ## store is3 

        movq  nb400nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb400nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movq  nb400nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb400nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb400nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb400nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb400nf_iq(%rsp)
        movaps %xmm4,nb400nf_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb400nf_ix(%rsp)
        movaps %xmm1,nb400nf_iy(%rsp)
        movaps %xmm2,nb400nf_iz(%rsp)

        movl  %ebx,nb400nf_ii3(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb400nf_vctot(%rsp)

        movq  nb400nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb400nf_pos(%rbp),%rsi
        movq  nb400nf_faction(%rbp),%rdi
        movq  nb400nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb400nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb400nf_ninner(%rsp),%ecx
        movl  %ecx,nb400nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb400nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel400nf_x86_64_sse.nb400nf_unroll_loop
        jmp   _nb_kernel400nf_x86_64_sse.nb400nf_finish_inner
_nb_kernel400nf_x86_64_sse.nb400nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb400nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb400nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        ## load isa2
        movq nb400nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7
        movaps nb400nf_isai(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb400nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb400nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb400nf_gbscale(%rsp)

        movq nb400nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        mulps nb400nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400nf_qq(%rsp)


        movq nb400nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb400nf_ix(%rsp),%xmm4
        movaps nb400nf_iy(%rsp),%xmm5
        movaps nb400nf_iz(%rsp),%xmm6

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
        movaps nb400nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        mulps nb400nf_gbscale(%rsp),%xmm4

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

        movq nb400nf_GBtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%rsi,%rax,4),%xmm4
        movaps (%rsi,%rbx,4),%xmm5
        movaps (%rsi,%rcx,4),%xmm6
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp 
        movaps nb400nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb400nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb400nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb400nf_innerk(%rsp)
        jl    _nb_kernel400nf_x86_64_sse.nb400nf_finish_inner
        jmp   _nb_kernel400nf_x86_64_sse.nb400nf_unroll_loop
_nb_kernel400nf_x86_64_sse.nb400nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb400nf_innerk(%rsp)
        movl  nb400nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel400nf_x86_64_sse.nb400nf_dopair
        jmp   _nb_kernel400nf_x86_64_sse.nb400nf_checksingle
_nb_kernel400nf_x86_64_sse.nb400nf_dopair: 
        movq  nb400nf_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb400nf_innerjjnr(%rsp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isa2
        movq nb400nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm2
        movss (%rsi,%rbx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isa2 in xmm3(0,1)
        mulps  nb400nf_isai(%rsp),%xmm2
        movaps %xmm2,nb400nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb400nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb400nf_gbscale(%rsp)

        movq nb400nf_charge(%rbp),%rsi     ## base of charge[]  
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb400nf_iq(%rsp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb400nf_qq(%rsp)

        movq nb400nf_pos(%rbp),%rdi

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

        movq   nb400nf_faction(%rbp),%rdi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb400nf_ix(%rsp),%xmm4
        movaps nb400nf_iy(%rsp),%xmm5
        movaps nb400nf_iz(%rsp),%xmm6

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
        movaps nb400nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb400nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb400nf_gbscale(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb400nf_GBtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%rsi,%rcx,4),%xmm4
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb400nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        addps  nb400nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb400nf_vctot(%rsp)

_nb_kernel400nf_x86_64_sse.nb400nf_checksingle: 
        movl  nb400nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400nf_x86_64_sse.nb400nf_dosingle
        jmp    _nb_kernel400nf_x86_64_sse.nb400nf_updateouterdata
_nb_kernel400nf_x86_64_sse.nb400nf_dosingle: 
        movq nb400nf_charge(%rbp),%rsi
        movq nb400nf_invsqrta(%rbp),%rdx
        movq nb400nf_pos(%rbp),%rdi
        movq  nb400nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%rdx,%rax,4),%xmm2       ## isa2
        mulss nb400nf_isai(%rsp),%xmm2
        movss %xmm2,nb400nf_isaprod(%rsp)
        movss %xmm2,%xmm1
        mulss nb400nf_gbtsc(%rsp),%xmm1
        movss %xmm1,nb400nf_gbscale(%rsp)

        mulss  nb400nf_iq(%rsp),%xmm2
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb400nf_qq(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movss nb400nf_ix(%rsp),%xmm4
        movss nb400nf_iy(%rsp),%xmm5
        movss nb400nf_iz(%rsp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movss %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb400nf_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb400nf_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        mulss nb400nf_gbscale(%rsp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movss %xmm4,%xmm1       ## xmm1=eps 
        movss %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx

        movq nb400nf_GBtab(%rbp),%rsi

        movaps (%rsi,%rbx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        movss nb400nf_qq(%rsp),%xmm3
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        addss  nb400nf_vctot(%rsp),%xmm5
        movss %xmm5,nb400nf_vctot(%rsp)
_nb_kernel400nf_x86_64_sse.nb400nf_updateouterdata: 
        ## get n from stack
        movl nb400nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb400nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb400nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb400nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb400nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400nf_x86_64_sse.nb400nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400nf_n(%rsp)
        jmp _nb_kernel400nf_x86_64_sse.nb400nf_outer
_nb_kernel400nf_x86_64_sse.nb400nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400nf_x86_64_sse.nb400nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400nf_x86_64_sse.nb400nf_threadloop
_nb_kernel400nf_x86_64_sse.nb400nf_end: 

        movl nb400nf_nouter(%rsp),%eax
        movl nb400nf_ninner(%rsp),%ebx
        movq nb400nf_outeriter(%rbp),%rcx
        movq nb400nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $296,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


