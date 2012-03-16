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






.globl nb_kernel300_x86_64_sse2
.globl _nb_kernel300_x86_64_sse2
nb_kernel300_x86_64_sse2:       
_nb_kernel300_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb300_is3, 272
.set nb300_ii3, 276
.set nb300_nri, 280
.set nb300_iinr, 288
.set nb300_jindex, 296
.set nb300_jjnr, 304
.set nb300_shift, 312
.set nb300_shiftvec, 320
.set nb300_facel, 328
.set nb300_innerjjnr, 336
.set nb300_innerk, 344
.set nb300_n, 348
.set nb300_nn1, 352
.set nb300_nouter, 356
.set nb300_ninner, 360
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb300_facel(%rsp)

        movq nb300_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb300_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb300_half(%rsp)
        movl %ebx,nb300_half+4(%rsp)
        movsd nb300_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb300_half(%rsp)
        movapd %xmm2,nb300_two(%rsp)
        movapd %xmm3,nb300_three(%rsp)

_nb_kernel300_x86_64_sse2.nb300_threadloop: 
        movq  nb300_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel300_x86_64_sse2.nb300_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300_x86_64_sse2.nb300_spinlock

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
        jg  _nb_kernel300_x86_64_sse2.nb300_outerstart
        jmp _nb_kernel300_x86_64_sse2.nb300_end

_nb_kernel300_x86_64_sse2.nb300_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300_nouter(%rsp),%ebx
        movl %ebx,nb300_nouter(%rsp)

_nb_kernel300_x86_64_sse2.nb300_outer: 
        movq  nb300_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb300_is3(%rsp)      ## store is3 

        movq  nb300_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb300_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb300_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb300_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb300_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb300_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb300_ix(%rsp)
        movapd %xmm1,nb300_iy(%rsp)
        movapd %xmm2,nb300_iz(%rsp)

        movl  %ebx,nb300_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorpd %xmm12,%xmm12
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15

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
        subl  $2,%edx
        addl  nb300_ninner(%rsp),%ecx
        movl  %ecx,nb300_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb300_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel300_x86_64_sse2.nb300_unroll_loop
        jmp   _nb_kernel300_x86_64_sse2.nb300_checksingle
_nb_kernel300_x86_64_sse2.nb300_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb300_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%r10d
        movl  4(%rdx),%r11d
        addq $8,nb300_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb300_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%r10,%r10,2),%rax     ## replace jnr with j3 
        lea  (%r11,%r11,2),%rbx

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        movq   nb300_faction(%rbp),%rdi

        ## calc dr 
        subpd nb300_ix(%rsp),%xmm4
        subpd nb300_iy(%rsp),%xmm5
        subpd nb300_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb300_charge(%rbp),%rsi     ## base of charge[] 

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        movlpd (%rsi,%r10,8),%xmm3

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        movhpd (%rsi,%r11,8),%xmm3


        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb300_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        mulpd  nb300_iq(%rsp),%xmm3

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb300_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb300_tsc(%rsp),%xmm4
        movapd %xmm3,nb300_qq(%rsp)

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movq nb300_VFtab(%rbp),%rsi
        movd %mm6,%r8d
        psrlq $32,%mm6
        movd %mm6,%r9d          ## indices in eax/ebx 

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        movapd (%rsi,%r9,8),%xmm3       ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        movapd 16(%rsi,%r9,8),%xmm3     ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        addpd  %xmm7,%xmm7      ## two*Heps2 
        movapd nb300_qq(%rsp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  %xmm5,%xmm12

        xorpd  %xmm4,%xmm4

        mulpd nb300_tsc(%rsp),%xmm3
        mulpd %xmm0,%xmm3
        subpd  %xmm3,%xmm4


        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        movq   nb300_faction(%rbp),%rdi
        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb300_innerk(%rsp)
        jl    _nb_kernel300_x86_64_sse2.nb300_checksingle
        jmp   _nb_kernel300_x86_64_sse2.nb300_unroll_loop
_nb_kernel300_x86_64_sse2.nb300_checksingle: 
        movl  nb300_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300_x86_64_sse2.nb300_dosingle
        jmp    _nb_kernel300_x86_64_sse2.nb300_updateouterdata
_nb_kernel300_x86_64_sse2.nb300_dosingle: 
        movq nb300_charge(%rbp),%rsi
        movq nb300_pos(%rbp),%rdi
        movq  nb300_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb300_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3

        movapd nb300_iq(%rsp),%xmm2
        mulsd  %xmm2,%xmm3
        movapd %xmm3,nb300_qq(%rsp)

        movq nb300_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        movq   nb300_faction(%rbp),%rdi

        ## calc dr 
        subsd nb300_ix(%rsp),%xmm4
        subsd nb300_iy(%rsp),%xmm5
        subsd nb300_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb300_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb300_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb300_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r8d            ## idx *= 4 

        movq nb300_VFtab(%rbp),%rsi

    movsd  (%rsi,%r8,8),%xmm4
    movsd  8(%rsi,%r8,8),%xmm5
    movsd  16(%rsi,%r8,8),%xmm6
    movsd  24(%rsi,%r8,8),%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        addsd  %xmm7,%xmm7      ## two*Heps2 
        movapd nb300_qq(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  %xmm5,%xmm12

        xorpd  %xmm4,%xmm4

        mulsd nb300_tsc(%rsp),%xmm3
        mulsd %xmm0,%xmm3
        subsd  %xmm3,%xmm4

        movq   nb300_faction(%rbp),%rdi
        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15

        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel300_x86_64_sse2.nb300_updateouterdata: 
        movl  nb300_ii3(%rsp),%ecx
        movq  nb300_faction(%rbp),%rdi
        movq  nb300_fshift(%rbp),%rsi
        movl  nb300_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm3
        movhlps %xmm14,%xmm4
        movhlps %xmm15,%xmm5
        addsd  %xmm3,%xmm13
        addsd  %xmm4,%xmm14
        addsd  %xmm5,%xmm15 ## sum is in low xmm13-xmm15

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb300_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb300_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total coulomb energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 have the sum now 

        ## add earlier value from mem 
        movq  nb300_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## finish if last 
        movl nb300_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300_x86_64_sse2.nb300_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300_n(%rsp)
        jmp _nb_kernel300_x86_64_sse2.nb300_outer
_nb_kernel300_x86_64_sse2.nb300_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300_x86_64_sse2.nb300_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300_x86_64_sse2.nb300_threadloop
_nb_kernel300_x86_64_sse2.nb300_end: 
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






.globl nb_kernel300nf_x86_64_sse2
.globl _nb_kernel300nf_x86_64_sse2
nb_kernel300nf_x86_64_sse2:     
_nb_kernel300nf_x86_64_sse2:    
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
        ## stack offsets for local variables  
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
.set nb300nf_nri, 152
.set nb300nf_iinr, 160
.set nb300nf_jindex, 168
.set nb300nf_jjnr, 176
.set nb300nf_shift, 184
.set nb300nf_shiftvec, 192
.set nb300nf_facel, 200
.set nb300nf_innerjjnr, 208
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb300nf_facel(%rsp)

        movq nb300nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb300nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb300nf_half(%rsp)
        movl %ebx,nb300nf_half+4(%rsp)
        movsd nb300nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb300nf_half(%rsp)
        movapd %xmm3,nb300nf_three(%rsp)

_nb_kernel300nf_x86_64_sse2.nb300nf_threadloop: 
        movq  nb300nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel300nf_x86_64_sse2.nb300nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300nf_x86_64_sse2.nb300nf_spinlock

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
        jg  _nb_kernel300nf_x86_64_sse2.nb300nf_outerstart
        jmp _nb_kernel300nf_x86_64_sse2.nb300nf_end

_nb_kernel300nf_x86_64_sse2.nb300nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300nf_nouter(%rsp),%ebx
        movl %ebx,nb300nf_nouter(%rsp)

_nb_kernel300nf_x86_64_sse2.nb300nf_outer: 
        movq  nb300nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb300nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb300nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb300nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb300nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb300nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb300nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb300nf_ix(%rsp)
        movapd %xmm1,nb300nf_iy(%rsp)
        movapd %xmm2,nb300nf_iz(%rsp)

        movl  %ebx,nb300nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb300nf_vctot(%rsp)

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
        subl  $2,%edx
        addl  nb300nf_ninner(%rsp),%ecx
        movl  %ecx,nb300nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb300nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel300nf_x86_64_sse2.nb300nf_unroll_loop
        jmp   _nb_kernel300nf_x86_64_sse2.nb300nf_checksingle
_nb_kernel300nf_x86_64_sse2.nb300nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb300nf_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb300nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb300nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb300nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb300nf_qq(%rsp)

        movq nb300nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move nb300nf_ix-iz to xmm4-xmm6 
        movapd nb300nf_ix(%rsp),%xmm4
        movapd nb300nf_iy(%rsp),%xmm5
        movapd nb300nf_iz(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb300nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb300nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb300nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movq nb300nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb300nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb300nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb300nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb300nf_innerk(%rsp)
        jl    _nb_kernel300nf_x86_64_sse2.nb300nf_checksingle
        jmp   _nb_kernel300nf_x86_64_sse2.nb300nf_unroll_loop
_nb_kernel300nf_x86_64_sse2.nb300nf_checksingle: 
        movl  nb300nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300nf_x86_64_sse2.nb300nf_dosingle
        jmp    _nb_kernel300nf_x86_64_sse2.nb300nf_updateouterdata
_nb_kernel300nf_x86_64_sse2.nb300nf_dosingle: 
        movq nb300nf_charge(%rbp),%rsi
        movq nb300nf_pos(%rbp),%rdi
        movq  nb300nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorpd  %xmm6,%xmm6
        movlpd (%rsi,%rax,8),%xmm6      ## xmm6(0) has the charge       
        mulsd  nb300nf_iq(%rsp),%xmm6
        movapd %xmm6,nb300nf_qq(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2

        ## move nb300nf_ix-iz to xmm4-xmm6 
        movapd nb300nf_ix(%rsp),%xmm4
        movapd nb300nf_iy(%rsp),%xmm5
        movapd nb300nf_iz(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb300nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb300nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb300nf_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movq nb300nf_VFtab(%rbp),%rsi

        ## Coulomb 
        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1   
        ## table ready in xmm4-xmm7 

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb300nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb300nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb300nf_vctot(%rsp)

_nb_kernel300nf_x86_64_sse2.nb300nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb300nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb300nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb300nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb300nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb300nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300nf_x86_64_sse2.nb300nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300nf_n(%rsp)
        jmp _nb_kernel300nf_x86_64_sse2.nb300nf_outer
_nb_kernel300nf_x86_64_sse2.nb300nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300nf_x86_64_sse2.nb300nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300nf_x86_64_sse2.nb300nf_threadloop
_nb_kernel300nf_x86_64_sse2.nb300nf_end: 
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

