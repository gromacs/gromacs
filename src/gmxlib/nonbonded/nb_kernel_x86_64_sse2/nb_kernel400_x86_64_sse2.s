##
## $Id$
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





.globl nb_kernel400_x86_64_sse2
.globl _nb_kernel400_x86_64_sse2
nb_kernel400_x86_64_sse2:       
_nb_kernel400_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb400_facel(%rsp)

        movq nb400_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb400_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb400_half(%rsp)
        movl %ebx,nb400_half+4(%rsp)
        movsd nb400_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb400_half(%rsp)
        movapd %xmm2,nb400_two(%rsp)
        movapd %xmm3,nb400_three(%rsp)

_nb_kernel400_x86_64_sse2.nb400_threadloop: 
        movq  nb400_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel400_x86_64_sse2.nb400_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400_x86_64_sse2.nb400_spinlock

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
        jg  _nb_kernel400_x86_64_sse2.nb400_outerstart
        jmp _nb_kernel400_x86_64_sse2.nb400_end

_nb_kernel400_x86_64_sse2.nb400_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400_nouter(%rsp),%ebx
        movl %ebx,nb400_nouter(%rsp)

_nb_kernel400_x86_64_sse2.nb400_outer: 
        movq  nb400_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb400_is3(%rsp)      ## store is3 

        movq  nb400_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb400_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb400_ii(%rsp)

        movq  nb400_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb400_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb400_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb400_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb400_iq(%rsp)
        movapd %xmm4,nb400_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb400_ix(%rsp)
        movapd %xmm1,nb400_iy(%rsp)
        movapd %xmm2,nb400_iz(%rsp)

        movl  %ebx,nb400_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,%xmm8
        movapd %xmm4,%xmm12
        movapd %xmm4,%xmm13
        movapd %xmm4,%xmm14
        movapd %xmm4,%xmm15

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
        subl  $2,%edx
        addl  nb400_ninner(%rsp),%ecx
        movl  %ecx,nb400_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb400_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel400_x86_64_sse2.nb400_unroll_loop
        jmp   _nb_kernel400_x86_64_sse2.nb400_checksingle
_nb_kernel400_x86_64_sse2.nb400_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb400_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%r12d
        movl  4(%rdx),%r13d
        addq $8,nb400_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb400_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%r12,%r12,2),%r8     ## j3 
        lea  (%r13,%r13,2),%r9

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%r8,8),%xmm4
        movlpd 8(%rsi,%r8,8),%xmm5
        movlpd 16(%rsi,%r8,8),%xmm6
        movhpd (%rsi,%r9,8),%xmm4
        movhpd 8(%rsi,%r9,8),%xmm5
        movhpd 16(%rsi,%r9,8),%xmm6

        ## calc dr 
        subpd nb400_ix(%rsp),%xmm4
        subpd nb400_iy(%rsp),%xmm5
        subpd nb400_iz(%rsp),%xmm6


        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        movq nb400_invsqrta(%rbp),%rsi
        movlpd (%rsi,%r12,8),%xmm3

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        movhpd (%rsi,%r13,8),%xmm3
        mulpd  nb400_isai(%rsp),%xmm3
        movapd %xmm3,nb400_isaprod(%rsp)
    movapd %xmm3,%xmm6
        mulpd nb400_gbtsc(%rsp),%xmm3
        movapd %xmm3,nb400_gbscale(%rsp)

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb400_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movq nb400_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%r12,8),%xmm3

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb400_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 

    mulpd  nb400_iq(%rsp),%xmm6
        movhpd (%rsi,%r13,8),%xmm3
        mulpd  %xmm6,%xmm3
        movapd %xmm3,nb400_qq(%rsp)


        movapd %xmm4,nb400_r(%rsp)
        mulpd nb400_gbscale(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 

        pslld $2,%mm6           ## idx *= 4 

        movq nb400_GBtab(%rbp),%rsi
        movd %mm6,%r10d
        psrlq $32,%mm6
        movd %mm6,%r11d         ## indices in r10/r11

        movapd (%rsi,%r10,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%r11,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%r10,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%r11,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            

        mulpd  %xmm1,%xmm7      ## xmm7=Heps
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm1,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        addpd  %xmm7,%xmm7      ## two*Heps2 
        movapd nb400_qq(%rsp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 

        movq nb400_dvda(%rbp),%rsi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulpd nb400_gbscale(%rsp),%xmm3
        movapd %xmm3,%xmm6
        mulpd  nb400_r(%rsp),%xmm6
        addpd  %xmm5,%xmm6

    ## update vctot
        addpd  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addpd  %xmm7,%xmm8

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%rsi,%r12,8),%xmm6
        addsd  (%rsi,%r13,8),%xmm7
        movsd  %xmm6,(%rsi,%r12,8)
        movsd  %xmm7,(%rsi,%r13,8)

        ## the fj's - start by accumulating forces from memory 
    movq nb400_faction(%rbp),%rdi
        movlpd (%rdi,%r8,8),%xmm5
        movlpd 8(%rdi,%r8,8),%xmm6
        movlpd 16(%rdi,%r8,8),%xmm7
        movhpd (%rdi,%r9,8),%xmm5
        movhpd 8(%rdi,%r9,8),%xmm6
        movhpd 16(%rdi,%r9,8),%xmm7

        xorpd  %xmm4,%xmm4

        mulpd %xmm0,%xmm3
        subpd  %xmm3,%xmm4

        movq   nb400_faction(%rbp),%rdi
        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        addpd %xmm9,%xmm5
        addpd %xmm10,%xmm6
        addpd %xmm11,%xmm7

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        movlpd %xmm5,(%rdi,%r8,8)
        movlpd %xmm6,8(%rdi,%r8,8)
        movlpd %xmm7,16(%rdi,%r8,8)
        movhpd %xmm5,(%rdi,%r9,8)
        movhpd %xmm6,8(%rdi,%r9,8)
        movhpd %xmm7,16(%rdi,%r9,8)

        ## should we do one more iteration? 
        subl $2,nb400_innerk(%rsp)
        jl    _nb_kernel400_x86_64_sse2.nb400_checksingle
        jmp   _nb_kernel400_x86_64_sse2.nb400_unroll_loop
_nb_kernel400_x86_64_sse2.nb400_checksingle: 
        movl  nb400_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400_x86_64_sse2.nb400_dosingle
        jmp    _nb_kernel400_x86_64_sse2.nb400_updateouterdata
_nb_kernel400_x86_64_sse2.nb400_dosingle: 
        movq nb400_charge(%rbp),%rsi
        movq nb400_invsqrta(%rbp),%rdx
        movq nb400_pos(%rbp),%rdi
        movq  nb400_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb400_invsqrta(%rbp),%rsi
        movsd (%rsi,%rax,8),%xmm2
        mulsd  nb400_isai(%rsp),%xmm2
        movapd %xmm2,nb400_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulsd nb400_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb400_gbscale(%rsp)

    mulsd nb400_iq(%rsp),%xmm2
        movq nb400_charge(%rbp),%rsi     ## base of charge[] 
        movsd (%rsi,%rax,8),%xmm3
        mulsd  %xmm2,%xmm3
        movapd %xmm3,nb400_qq(%rsp)

        movq nb400_pos(%rbp),%rsi               ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3 

        ## move coordinate to xmm4-xmm6
        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm5
        movsd 16(%rsi,%r8,8),%xmm6

        movq   nb400_faction(%rbp),%rdi

        ## calc dr 
        subsd nb400_ix(%rsp),%xmm4
        subsd nb400_iy(%rsp),%xmm5
        subsd nb400_iz(%rsp),%xmm6

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
        movapd nb400_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb400_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 

        movapd %xmm4,nb400_r(%rsp)
        mulsd nb400_gbscale(%rsp),%xmm4

        cvttsd2si %xmm4,%r10d   ## mm6 = lu idx 
        cvtsi2sd %r10d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 

        shll $2,%r10d           ## idx *= 4 

        movq nb400_GBtab(%rbp),%rsi

        movapd (%rsi,%r10,8),%xmm4      ## Y1 F1        
        movhlps %xmm4,%xmm5
        movapd 16(%rsi,%r10,8),%xmm6    ## G1 H1        
    movhlps %xmm6,%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulsd  %xmm1,%xmm7      ## xmm7=Heps
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm1,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        addsd  %xmm7,%xmm7      ## two*Heps2 
        movapd nb400_qq(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 

        movq nb400_dvda(%rbp),%rsi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulsd nb400_gbscale(%rsp),%xmm3
        movapd %xmm3,%xmm6
        mulsd  nb400_r(%rsp),%xmm6
        addsd  %xmm5,%xmm6

    ## update vctot
        addsd  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subsd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addsd  %xmm7,%xmm8

        ## update j atoms dvdaj
        addsd  (%rsi,%rax,8),%xmm6
        movsd  %xmm6,(%rsi,%rax,8)

        xorpd  %xmm4,%xmm4

        mulsd %xmm0,%xmm3
        subsd  %xmm3,%xmm4

        movq   nb400_faction(%rbp),%rdi
        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15

        ## the fj's - start by accumulating forces from memory 
    movq nb400_faction(%rbp),%rdi
        addsd (%rdi,%r8,8),%xmm9
        addsd 8(%rdi,%r8,8),%xmm10
        addsd 16(%rdi,%r8,8),%xmm11
        movsd %xmm9,(%rdi,%r8,8)
        movsd %xmm10,8(%rdi,%r8,8)
        movsd %xmm11,16(%rdi,%r8,8)

_nb_kernel400_x86_64_sse2.nb400_updateouterdata: 
        movl  nb400_ii3(%rsp),%ecx
        movq  nb400_faction(%rbp),%rdi
        movq  nb400_fshift(%rbp),%rsi
        movl  nb400_is3(%rsp),%edx

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
        movl nb400_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb400_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total coulomb energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 have the sum now 

        ## add earlier value from mem 
        movq  nb400_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## accumulate dVda and update it 
        movhlps %xmm8,%xmm6
        addsd  %xmm6,%xmm8      ## low xmm8 has the sum now 

        movl nb400_ii(%rsp),%edx
        movq nb400_dvda(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm8
        movsd %xmm8,(%rax,%rdx,8)

        ## finish if last 
        movl nb400_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel400_x86_64_sse2.nb400_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400_n(%rsp)
        jmp _nb_kernel400_x86_64_sse2.nb400_outer
_nb_kernel400_x86_64_sse2.nb400_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel400_x86_64_sse2.nb400_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400_x86_64_sse2.nb400_threadloop
_nb_kernel400_x86_64_sse2.nb400_end: 
        movl nb400_nouter(%rsp),%eax
        movl nb400_ninner(%rsp),%ebx
        movq nb400_outeriter(%rbp),%rcx
        movq nb400_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret







.globl nb_kernel400nf_x86_64_sse2
.globl _nb_kernel400nf_x86_64_sse2
nb_kernel400nf_x86_64_sse2:     
_nb_kernel400nf_x86_64_sse2:    
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
        ## bottom of stack is cache-aligned for sse2 use 
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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb400nf_facel(%rsp)

        movq nb400nf_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb400nf_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb400nf_half(%rsp)
        movl %ebx,nb400nf_half+4(%rsp)
        movsd nb400nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb400nf_half(%rsp)
        movapd %xmm3,nb400nf_three(%rsp)

_nb_kernel400nf_x86_64_sse2.nb400nf_threadloop: 
        movq  nb400nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel400nf_x86_64_sse2.nb400nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400nf_x86_64_sse2.nb400nf_spinlock

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
        jg  _nb_kernel400nf_x86_64_sse2.nb400nf_outerstart
        jmp _nb_kernel400nf_x86_64_sse2.nb400nf_end

_nb_kernel400nf_x86_64_sse2.nb400nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400nf_nouter(%rsp),%ebx
        movl %ebx,nb400nf_nouter(%rsp)

_nb_kernel400nf_x86_64_sse2.nb400nf_outer: 
        movq  nb400nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb400nf_is3(%rsp)            ## store is3 

        movq  nb400nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb400nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb400nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb400nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb400nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb400nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb400nf_iq(%rsp)
        movapd %xmm4,nb400nf_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb400nf_ix(%rsp)
        movapd %xmm1,nb400nf_iy(%rsp)
        movapd %xmm2,nb400nf_iz(%rsp)

        movl  %ebx,nb400nf_ii3(%rsp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb400nf_vctot(%rsp)

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
        subl  $2,%edx
        addl  nb400nf_ninner(%rsp),%ecx
        movl  %ecx,nb400nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb400nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel400nf_x86_64_sse2.nb400nf_unroll_loop
        jmp   _nb_kernel400nf_x86_64_sse2.nb400nf_checksingle
_nb_kernel400nf_x86_64_sse2.nb400nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb400nf_innerjjnr(%rsp),%rdx     ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb400nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        ## load isa2
        movq nb400nf_invsqrta(%rbp),%rsi
        movlpd (%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm2
        mulpd  nb400nf_isai(%rsp),%xmm2
        movapd %xmm2,nb400nf_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulpd nb400nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb400nf_gbscale(%rsp)

        movq nb400nf_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        mulpd nb400nf_iq(%rsp),%xmm2
    mulpd %xmm2,%xmm3
        movapd %xmm3,nb400nf_qq(%rsp)

        movq nb400nf_pos(%rbp),%rsi             ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        movq   nb400nf_faction(%rbp),%rdi

        ## move nb400nf_ix-iz to xmm4-xmm6 
        movapd nb400nf_ix(%rsp),%xmm4
        movapd nb400nf_iy(%rsp),%xmm5
        movapd nb400nf_iz(%rsp),%xmm6

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
        movapd nb400nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb400nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb400nf_gbscale(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movq nb400nf_GBtab(%rbp),%rsi
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
        movapd nb400nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb400nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb400nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb400nf_innerk(%rsp)
        jl    _nb_kernel400nf_x86_64_sse2.nb400nf_checksingle
        jmp   _nb_kernel400nf_x86_64_sse2.nb400nf_unroll_loop
_nb_kernel400nf_x86_64_sse2.nb400nf_checksingle: 
        movl  nb400nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400nf_x86_64_sse2.nb400nf_dosingle
        jmp    _nb_kernel400nf_x86_64_sse2.nb400nf_updateouterdata
_nb_kernel400nf_x86_64_sse2.nb400nf_dosingle: 
        movq nb400nf_charge(%rbp),%rsi
        movq nb400nf_invsqrta(%rbp),%rdx
        movq nb400nf_pos(%rbp),%rdi
        movq  nb400nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%rdx,%rax,8),%xmm7
        movlpd (%rsi,%rax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb400nf_isai(%rsp),%xmm7
        movapd %xmm7,nb400nf_isaprod(%rsp)
        movapd %xmm7,%xmm1
        mulpd nb400nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb400nf_gbscale(%rsp)

        mulsd  nb400nf_iq(%rsp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb400nf_qq(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2

        ## move nb400nf_ix-iz to xmm4-xmm6 
        movapd nb400nf_ix(%rsp),%xmm4
        movapd nb400nf_iy(%rsp),%xmm5
        movapd nb400nf_iz(%rsp),%xmm6

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
        movapd nb400nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb400nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb400nf_gbscale(%rsp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movq nb400nf_GBtab(%rbp),%rsi

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
        movapd nb400nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb400nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb400nf_vctot(%rsp)

_nb_kernel400nf_x86_64_sse2.nb400nf_updateouterdata: 
        ## get n from stack
        movl nb400nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb400nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb400nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb400nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb400nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel400nf_x86_64_sse2.nb400nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400nf_n(%rsp)
        jmp _nb_kernel400nf_x86_64_sse2.nb400nf_outer
_nb_kernel400nf_x86_64_sse2.nb400nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel400nf_x86_64_sse2.nb400nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400nf_x86_64_sse2.nb400nf_threadloop
_nb_kernel400nf_x86_64_sse2.nb400nf_end: 

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





