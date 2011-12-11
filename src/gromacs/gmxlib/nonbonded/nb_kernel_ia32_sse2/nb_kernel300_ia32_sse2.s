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


.globl nb_kernel300_ia32_sse2
.globl _nb_kernel300_ia32_sse2
nb_kernel300_ia32_sse2: 
_nb_kernel300_ia32_sse2:        
.set nb300_p_nri, 8
.set nb300_iinr, 12
.set nb300_jindex, 16
.set nb300_jjnr, 20
.set nb300_shift, 24
.set nb300_shiftvec, 28
.set nb300_fshift, 32
.set nb300_gid, 36
.set nb300_pos, 40
.set nb300_faction, 44
.set nb300_charge, 48
.set nb300_p_facel, 52
.set nb300_argkrf, 56
.set nb300_argcrf, 60
.set nb300_Vc, 64
.set nb300_type, 68
.set nb300_p_ntype, 72
.set nb300_vdwparam, 76
.set nb300_Vvdw, 80
.set nb300_p_tabscale, 84
.set nb300_VFtab, 88
.set nb300_invsqrta, 92
.set nb300_dvda, 96
.set nb300_p_gbtabscale, 100
.set nb300_GBtab, 104
.set nb300_p_nthreads, 108
.set nb300_count, 112
.set nb300_mtx, 116
.set nb300_outeriter, 120
.set nb300_inneriter, 124
.set nb300_work, 128
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
.set nb300_innerjjnr, 280
.set nb300_innerk, 284
.set nb300_n, 288
.set nb300_nn1, 292
.set nb300_nri, 296
.set nb300_facel, 304                         ## uses 8 bytes
.set nb300_nouter, 312
.set nb300_ninner, 316
.set nb300_salign, 320
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $324,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb300_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb300_p_nri(%ebp),%ecx
        movl nb300_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb300_nri(%esp)
        movsd %xmm7,nb300_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300_nouter(%esp)
        movl %eax,nb300_ninner(%esp)


        movl nb300_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb300_tsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb300_half(%esp)
        movl %ebx,nb300_half+4(%esp)
        movsd nb300_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb300_half(%esp)
        movapd %xmm2,nb300_two(%esp)
        movapd %xmm3,nb300_three(%esp)

_nb_kernel300_ia32_sse2.nb300_threadloop: 
        movl  nb300_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300_ia32_sse2.nb300_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300_ia32_sse2.nb300_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300_n(%esp)
        movl %ebx,nb300_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300_ia32_sse2.nb300_outerstart
        jmp _nb_kernel300_ia32_sse2.nb300_end

_nb_kernel300_ia32_sse2.nb300_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300_nouter(%esp),%ebx
        movl %ebx,nb300_nouter(%esp)

_nb_kernel300_ia32_sse2.nb300_outer: 
        movl  nb300_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb300_is3(%esp)      ## store is3 

        movl  nb300_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb300_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb300_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb300_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb300_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb300_ix(%esp)
        movapd %xmm1,nb300_iy(%esp)
        movapd %xmm2,nb300_iz(%esp)

        movl  %ebx,nb300_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb300_vctot(%esp)
        movapd %xmm4,nb300_fix(%esp)
        movapd %xmm4,nb300_fiy(%esp)
        movapd %xmm4,nb300_fiz(%esp)

        movl  nb300_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300_pos(%ebp),%esi
        movl  nb300_faction(%ebp),%edi
        movl  nb300_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb300_ninner(%esp),%ecx
        movl  %ecx,nb300_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb300_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel300_ia32_sse2.nb300_unroll_loop
        jmp   _nb_kernel300_ia32_sse2.nb300_checksingle
_nb_kernel300_ia32_sse2.nb300_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb300_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb300_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb300_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb300_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb300_qq(%esp)

        movl nb300_pos(%ebp),%esi               ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb300_faction(%ebp),%edi

        ## move nb300_ix-iz to xmm4-xmm6 
        movapd nb300_ix(%esp),%xmm4
        movapd nb300_iy(%esp),%xmm5
        movapd nb300_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb300_dx(%esp)
        movapd %xmm5,nb300_dy(%esp)
        movapd %xmm6,nb300_dz(%esp)
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
        movapd nb300_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb300_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb300_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb300_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb300_two(%esp),%xmm7    ## two*Heps2 
        movapd nb300_qq(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb300_vctot(%esp),%xmm5
        movapd %xmm5,nb300_vctot(%esp)

        xorpd  %xmm4,%xmm4

        mulpd nb300_tsc(%esp),%xmm3
        mulpd %xmm0,%xmm3
        subpd  %xmm3,%xmm4

        movapd nb300_dx(%esp),%xmm0
        movapd nb300_dy(%esp),%xmm1
        movapd nb300_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb300_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb300_fix(%esp),%xmm3
        movapd nb300_fiy(%esp),%xmm4
        movapd nb300_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb300_fix(%esp)
        movapd %xmm4,nb300_fiy(%esp)
        movapd %xmm5,nb300_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb300_innerk(%esp)
        jl    _nb_kernel300_ia32_sse2.nb300_checksingle
        jmp   _nb_kernel300_ia32_sse2.nb300_unroll_loop
_nb_kernel300_ia32_sse2.nb300_checksingle: 
        movl  nb300_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300_ia32_sse2.nb300_dosingle
        jmp    _nb_kernel300_ia32_sse2.nb300_updateouterdata
_nb_kernel300_ia32_sse2.nb300_dosingle: 
        movl nb300_charge(%ebp),%esi
        movl nb300_pos(%ebp),%edi
        movl  nb300_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge       
        mulsd  nb300_iq(%esp),%xmm6
        movapd %xmm6,nb300_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2

        ## move nb300_ix-iz to xmm4-xmm6 
        movapd nb300_ix(%esp),%xmm4
        movapd nb300_iy(%esp),%xmm5
        movapd nb300_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb300_dx(%esp)
        movapd %xmm5,nb300_dy(%esp)
        movapd %xmm6,nb300_dz(%esp)
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
        movapd nb300_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb300_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb300_tsc(%esp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb300_VFtab(%ebp),%esi

        ## Coulomb 
        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1   
        ## table ready in xmm4-xmm7 

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb300_two(%esp),%xmm7    ## two*Heps2 
        movapd nb300_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb300_vctot(%esp),%xmm5
        movsd %xmm5,nb300_vctot(%esp)

        xorpd %xmm4,%xmm4
        movd %mm0,%eax

        mulpd nb300_tsc(%esp),%xmm3
        mulpd %xmm0,%xmm3
        subpd  %xmm3,%xmm4
        movl   nb300_faction(%ebp),%edi

        movapd nb300_dx(%esp),%xmm0
        movapd nb300_dy(%esp),%xmm1
        movapd nb300_dz(%esp),%xmm2

        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb300_fix(%esp),%xmm3
        movapd nb300_fiy(%esp),%xmm4
        movapd nb300_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb300_fix(%esp)
        movlpd %xmm4,nb300_fiy(%esp)
        movlpd %xmm5,nb300_fiz(%esp)
        ## update fj 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel300_ia32_sse2.nb300_updateouterdata: 
        movl  nb300_ii3(%esp),%ecx
        movl  nb300_faction(%ebp),%edi
        movl  nb300_fshift(%ebp),%esi
        movl  nb300_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb300_fix(%esp),%xmm0
        movapd nb300_fiy(%esp),%xmm1
        movapd nb300_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## increment fshift force  
        movsd  (%esi,%edx,8),%xmm3
        movsd  8(%esi,%edx,8),%xmm4
        movsd  16(%esi,%edx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb300_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb300_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb300_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb300_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300_ia32_sse2.nb300_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300_n(%esp)
        jmp _nb_kernel300_ia32_sse2.nb300_outer
_nb_kernel300_ia32_sse2.nb300_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300_ia32_sse2.nb300_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300_ia32_sse2.nb300_threadloop
_nb_kernel300_ia32_sse2.nb300_end: 
        emms

        movl nb300_nouter(%esp),%eax
        movl nb300_ninner(%esp),%ebx
        movl nb300_outeriter(%ebp),%ecx
        movl nb300_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb300_salign(%esp),%eax
        addl %eax,%esp
        addl $324,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel300nf_ia32_sse2
.globl _nb_kernel300nf_ia32_sse2
nb_kernel300nf_ia32_sse2:       
_nb_kernel300nf_ia32_sse2:      
.set nb300nf_p_nri, 8
.set nb300nf_iinr, 12
.set nb300nf_jindex, 16
.set nb300nf_jjnr, 20
.set nb300nf_shift, 24
.set nb300nf_shiftvec, 28
.set nb300nf_fshift, 32
.set nb300nf_gid, 36
.set nb300nf_pos, 40
.set nb300nf_faction, 44
.set nb300nf_charge, 48
.set nb300nf_p_facel, 52
.set nb300nf_argkrf, 56
.set nb300nf_argcrf, 60
.set nb300nf_Vc, 64
.set nb300nf_type, 68
.set nb300nf_p_ntype, 72
.set nb300nf_vdwparam, 76
.set nb300nf_Vvdw, 80
.set nb300nf_p_tabscale, 84
.set nb300nf_VFtab, 88
.set nb300nf_invsqrta, 92
.set nb300nf_dvda, 96
.set nb300nf_p_gbtabscale, 100
.set nb300nf_GBtab, 104
.set nb300nf_p_nthreads, 108
.set nb300nf_count, 112
.set nb300nf_mtx, 116
.set nb300nf_outeriter, 120
.set nb300nf_inneriter, 124
.set nb300nf_work, 128
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
.set nb300nf_innerjjnr, 152
.set nb300nf_innerk, 156
.set nb300nf_n, 160
.set nb300nf_nn1, 164
.set nb300nf_nri, 168
.set nb300nf_facel, 176                       ## uses 8 bytes
.set nb300nf_nouter, 184
.set nb300nf_ninner, 188
.set nb300nf_salign, 192
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $196,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb300nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb300nf_p_nri(%ebp),%ecx
        movl nb300nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb300nf_nri(%esp)
        movsd %xmm7,nb300nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb300nf_nouter(%esp)
        movl %eax,nb300nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb300nf_half(%esp)
        movl %ebx,nb300nf_half+4(%esp)
        movsd nb300nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb300nf_half(%esp)
        movapd %xmm3,nb300nf_three(%esp)
        movl nb300nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb300nf_tsc(%esp)

_nb_kernel300nf_ia32_sse2.nb300nf_threadloop: 
        movl  nb300nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel300nf_ia32_sse2.nb300nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel300nf_ia32_sse2.nb300nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb300nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb300nf_n(%esp)
        movl %ebx,nb300nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel300nf_ia32_sse2.nb300nf_outerstart
        jmp _nb_kernel300nf_ia32_sse2.nb300nf_end

_nb_kernel300nf_ia32_sse2.nb300nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb300nf_nouter(%esp),%ebx
        movl %ebx,nb300nf_nouter(%esp)

_nb_kernel300nf_ia32_sse2.nb300nf_outer: 
        movl  nb300nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb300nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb300nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb300nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb300nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb300nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb300nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb300nf_ix(%esp)
        movapd %xmm1,nb300nf_iy(%esp)
        movapd %xmm2,nb300nf_iz(%esp)

        movl  %ebx,nb300nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb300nf_vctot(%esp)

        movl  nb300nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb300nf_pos(%ebp),%esi
        movl  nb300nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb300nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb300nf_ninner(%esp),%ecx
        movl  %ecx,nb300nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb300nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel300nf_ia32_sse2.nb300nf_unroll_loop
        jmp   _nb_kernel300nf_ia32_sse2.nb300nf_checksingle
_nb_kernel300nf_ia32_sse2.nb300nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb300nf_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb300nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb300nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb300nf_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb300nf_qq(%esp)

        movl nb300nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move nb300nf_ix-iz to xmm4-xmm6 
        movapd nb300nf_ix(%esp),%xmm4
        movapd nb300nf_iy(%esp),%xmm5
        movapd nb300nf_iz(%esp),%xmm6

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
        movapd nb300nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb300nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb300nf_tsc(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movl nb300nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb300nf_qq(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul  
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb300nf_vctot(%esp),%xmm5
        movapd %xmm5,nb300nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb300nf_innerk(%esp)
        jl    _nb_kernel300nf_ia32_sse2.nb300nf_checksingle
        jmp   _nb_kernel300nf_ia32_sse2.nb300nf_unroll_loop
_nb_kernel300nf_ia32_sse2.nb300nf_checksingle: 
        movl  nb300nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel300nf_ia32_sse2.nb300nf_dosingle
        jmp    _nb_kernel300nf_ia32_sse2.nb300nf_updateouterdata
_nb_kernel300nf_ia32_sse2.nb300nf_dosingle: 
        movl nb300nf_charge(%ebp),%esi
        movl nb300nf_pos(%ebp),%edi
        movl  nb300nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm6,%xmm6
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge       
        mulsd  nb300nf_iq(%esp),%xmm6
        movapd %xmm6,nb300nf_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2

        ## move nb300nf_ix-iz to xmm4-xmm6 
        movapd nb300nf_ix(%esp),%xmm4
        movapd nb300nf_iy(%esp),%xmm5
        movapd nb300nf_iz(%esp),%xmm6

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
        movapd nb300nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb300nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb300nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb300nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb300nf_tsc(%esp),%xmm4

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb300nf_VFtab(%ebp),%esi

        ## Coulomb 
        movlpd (%esi,%eax,8),%xmm4      ## Y1
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1   
        ## table ready in xmm4-xmm7 

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb300nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb300nf_vctot(%esp),%xmm5
        movsd %xmm5,nb300nf_vctot(%esp)

_nb_kernel300nf_ia32_sse2.nb300nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb300nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb300nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb300nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb300nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb300nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel300nf_ia32_sse2.nb300nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb300nf_n(%esp)
        jmp _nb_kernel300nf_ia32_sse2.nb300nf_outer
_nb_kernel300nf_ia32_sse2.nb300nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb300nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel300nf_ia32_sse2.nb300nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel300nf_ia32_sse2.nb300nf_threadloop
_nb_kernel300nf_ia32_sse2.nb300nf_end: 
        emms

        movl nb300nf_nouter(%esp),%eax
        movl nb300nf_ninner(%esp),%ebx
        movl nb300nf_outeriter(%ebp),%ecx
        movl nb300nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb300nf_salign(%esp),%eax
        addl %eax,%esp
        addl $196,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret


