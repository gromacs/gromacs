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



.globl nb_kernel400_ia32_sse2
.globl _nb_kernel400_ia32_sse2
nb_kernel400_ia32_sse2: 
_nb_kernel400_ia32_sse2:        
.set nb400_p_nri, 8
.set nb400_iinr, 12
.set nb400_jindex, 16
.set nb400_jjnr, 20
.set nb400_shift, 24
.set nb400_shiftvec, 28
.set nb400_fshift, 32
.set nb400_gid, 36
.set nb400_pos, 40
.set nb400_faction, 44
.set nb400_charge, 48
.set nb400_p_facel, 52
.set nb400_argkrf, 56
.set nb400_argcrf, 60
.set nb400_Vc, 64
.set nb400_type, 68
.set nb400_p_ntype, 72
.set nb400_vdwparam, 76
.set nb400_Vvdw, 80
.set nb400_p_tabscale, 84
.set nb400_VFtab, 88
.set nb400_invsqrta, 92
.set nb400_dvda, 96
.set nb400_p_gbtabscale, 100
.set nb400_GBtab, 104
.set nb400_p_nthreads, 108
.set nb400_count, 112
.set nb400_mtx, 116
.set nb400_outeriter, 120
.set nb400_inneriter, 124
.set nb400_work, 128
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
.set nb400_is3, 336
.set nb400_ii3, 340
.set nb400_ii, 344
.set nb400_innerjjnr, 348
.set nb400_innerk, 352
.set nb400_n, 356
.set nb400_nn1, 360
.set nb400_nri, 364
.set nb400_facel, 368                         ## uses 8 bytes
.set nb400_nouter, 376
.set nb400_ninner, 380
.set nb400_salign, 384
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $388,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb400_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb400_p_nri(%ebp),%ecx
        movl nb400_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb400_nri(%esp)
        movsd %xmm7,nb400_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb400_nouter(%esp)
        movl %eax,nb400_ninner(%esp)


        movl nb400_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb400_gbtsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb400_half(%esp)
        movl %ebx,nb400_half+4(%esp)
        movsd nb400_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb400_half(%esp)
        movapd %xmm2,nb400_two(%esp)
        movapd %xmm3,nb400_three(%esp)

_nb_kernel400_ia32_sse2.nb400_threadloop: 
        movl  nb400_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel400_ia32_sse2.nb400_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400_ia32_sse2.nb400_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400_n(%esp)
        movl %ebx,nb400_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400_ia32_sse2.nb400_outerstart
        jmp _nb_kernel400_ia32_sse2.nb400_end

_nb_kernel400_ia32_sse2.nb400_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400_nouter(%esp),%ebx
        movl %ebx,nb400_nouter(%esp)

_nb_kernel400_ia32_sse2.nb400_outer: 
        movl  nb400_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb400_is3(%esp)      ## store is3 

        movl  nb400_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb400_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb400_ii(%esp)

        movl  nb400_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb400_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb400_invsqrta(%ebp),%edx         ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb400_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb400_iq(%esp)
        movapd %xmm4,nb400_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb400_ix(%esp)
        movapd %xmm1,nb400_iy(%esp)
        movapd %xmm2,nb400_iz(%esp)

        movl  %ebx,nb400_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb400_vctot(%esp)
        movapd %xmm4,nb400_dvdasum(%esp)
        movapd %xmm4,nb400_fix(%esp)
        movapd %xmm4,nb400_fiy(%esp)
        movapd %xmm4,nb400_fiz(%esp)

        movl  nb400_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb400_pos(%ebp),%esi
        movl  nb400_faction(%ebp),%edi
        movl  nb400_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb400_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb400_ninner(%esp),%ecx
        movl  %ecx,nb400_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb400_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel400_ia32_sse2.nb400_unroll_loop
        jmp   _nb_kernel400_ia32_sse2.nb400_checksingle
_nb_kernel400_ia32_sse2.nb400_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb400_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb400_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        ## load isaj
        movl nb400_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb400_isai(%esp),%xmm2
        movapd %xmm2,nb400_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb400_gbtsc(%esp),%xmm1
        movapd %xmm1,nb400_gbscale(%esp)

        movl nb400_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb400_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb400_qq(%esp)

        movl nb400_pos(%ebp),%esi               ## base of pos[] 

        movd  %eax,%mm2
        movd  %ebx,%mm3
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb400_faction(%ebp),%edi

        ## move nb400_ix-iz to xmm4-xmm6 
        movapd nb400_ix(%esp),%xmm4
        movapd nb400_iy(%esp),%xmm5
        movapd nb400_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb400_dx(%esp)
        movapd %xmm5,nb400_dy(%esp)
        movapd %xmm6,nb400_dz(%esp)
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
        movapd nb400_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb400_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb400_r(%esp)
        mulpd nb400_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb400_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        movapd (%esi,%ebx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb400_two(%esp),%xmm7    ## two*Heps2 
        movapd nb400_qq(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## get jnr from regs
        movd %mm2,%ecx
        movd %mm3,%edx
        movl nb400_dvda(%ebp),%esi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulpd nb400_gbscale(%esp),%xmm3
        movapd %xmm3,%xmm6
        mulpd  nb400_r(%esp),%xmm6
        addpd  %xmm5,%xmm6
        addpd  nb400_vctot(%esp),%xmm5
        movapd %xmm5,nb400_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addpd  nb400_dvdasum(%esp),%xmm7
        movapd %xmm7,nb400_dvdasum(%esp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%esi,%ecx,8),%xmm6
        addsd  (%esi,%edx,8),%xmm7
        movsd  %xmm6,(%esi,%ecx,8)
        movsd  %xmm7,(%esi,%edx,8)

        xorpd  %xmm4,%xmm4

        mulpd %xmm0,%xmm3
        subpd  %xmm3,%xmm4

        movapd nb400_dx(%esp),%xmm0
        movapd nb400_dy(%esp),%xmm1
        movapd nb400_dz(%esp),%xmm2

        movd %mm0,%eax
        movd %mm1,%ebx

        movl   nb400_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb400_fix(%esp),%xmm3
        movapd nb400_fiy(%esp),%xmm4
        movapd nb400_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb400_fix(%esp)
        movapd %xmm4,nb400_fiy(%esp)
        movapd %xmm5,nb400_fiz(%esp)
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
        subl $2,nb400_innerk(%esp)
        jl    _nb_kernel400_ia32_sse2.nb400_checksingle
        jmp   _nb_kernel400_ia32_sse2.nb400_unroll_loop
_nb_kernel400_ia32_sse2.nb400_checksingle: 
        movl  nb400_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400_ia32_sse2.nb400_dosingle
        jmp    _nb_kernel400_ia32_sse2.nb400_updateouterdata
_nb_kernel400_ia32_sse2.nb400_dosingle: 
        movl nb400_charge(%ebp),%esi
        movl nb400_invsqrta(%ebp),%edx
        movl nb400_pos(%ebp),%edi
        movl  nb400_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb400_isai(%esp),%xmm7
        movapd %xmm7,nb400_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb400_gbtsc(%esp),%xmm1
        movapd %xmm1,nb400_gbscale(%esp)

        mulsd  nb400_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb400_qq(%esp)

        movd  %eax,%mm2
        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2

        ## move nb400_ix-iz to xmm4-xmm6 
        movapd nb400_ix(%esp),%xmm4
        movapd nb400_iy(%esp),%xmm5
        movapd nb400_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb400_dx(%esp)
        movapd %xmm5,nb400_dy(%esp)
        movapd %xmm6,nb400_dz(%esp)
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
        movapd nb400_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb400_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb400_r(%esp)
        mulsd nb400_gbscale(%esp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb400_GBtab(%ebp),%esi

        ## Coulomb 
        movapd (%esi,%eax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1   
        ## table ready in xmm4-xmm7 

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb400_two(%esp),%xmm7    ## two*Heps2 
        movapd nb400_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq
        ## get jnr from regs
        movd %mm2,%ebx
        movl nb400_dvda(%ebp),%esi

        ## Calculate dVda
        mulsd nb400_gbscale(%esp),%xmm3
        movsd %xmm3,%xmm6
        mulsd  nb400_r(%esp),%xmm6
        addsd  %xmm5,%xmm6
        addsd  nb400_vctot(%esp),%xmm5
        movsd %xmm5,nb400_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movsd %xmm7,%xmm6

        ## update dvdasum
        addsd  nb400_dvdasum(%esp),%xmm7
        movsd %xmm7,nb400_dvdasum(%esp)

        ## update j atoms dvdaj
        addsd  (%esi,%ebx,8),%xmm6
        movsd  %xmm6,(%esi,%ebx,8)

        xorpd %xmm4,%xmm4
        movd %mm0,%eax

        mulsd %xmm0,%xmm3
        subsd  %xmm3,%xmm4
        movl   nb400_faction(%ebp),%edi

        movsd nb400_dx(%esp),%xmm0
        movsd nb400_dy(%esp),%xmm1
        movsd nb400_dz(%esp),%xmm2

        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movsd nb400_fix(%esp),%xmm3
        movsd nb400_fiy(%esp),%xmm4
        movsd nb400_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb400_fix(%esp)
        movlpd %xmm4,nb400_fiy(%esp)
        movlpd %xmm5,nb400_fiz(%esp)
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

_nb_kernel400_ia32_sse2.nb400_updateouterdata: 
        movl  nb400_ii3(%esp),%ecx
        movl  nb400_faction(%ebp),%edi
        movl  nb400_fshift(%ebp),%esi
        movl  nb400_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb400_fix(%esp),%xmm0
        movapd nb400_fiy(%esp),%xmm1
        movapd nb400_fiz(%esp),%xmm2

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
        movl nb400_n(%esp),%esi
        ## get group index for i particle 
        movl  nb400_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb400_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb400_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate dVda and update it 
        movapd nb400_dvdasum(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        movl nb400_ii(%esp),%edx
        movl nb400_dvda(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb400_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400_ia32_sse2.nb400_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400_n(%esp)
        jmp _nb_kernel400_ia32_sse2.nb400_outer
_nb_kernel400_ia32_sse2.nb400_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400_ia32_sse2.nb400_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400_ia32_sse2.nb400_threadloop
_nb_kernel400_ia32_sse2.nb400_end: 
        emms

        movl nb400_nouter(%esp),%eax
        movl nb400_ninner(%esp),%ebx
        movl nb400_outeriter(%ebp),%ecx
        movl nb400_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb400_salign(%esp),%eax
        addl %eax,%esp
        addl $388,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel400nf_ia32_sse2
.globl _nb_kernel400nf_ia32_sse2
nb_kernel400nf_ia32_sse2:       
_nb_kernel400nf_ia32_sse2:      
.set nb400nf_p_nri, 8
.set nb400nf_iinr, 12
.set nb400nf_jindex, 16
.set nb400nf_jjnr, 20
.set nb400nf_shift, 24
.set nb400nf_shiftvec, 28
.set nb400nf_fshift, 32
.set nb400nf_gid, 36
.set nb400nf_pos, 40
.set nb400nf_faction, 44
.set nb400nf_charge, 48
.set nb400nf_p_facel, 52
.set nb400nf_argkrf, 56
.set nb400nf_argcrf, 60
.set nb400nf_Vc, 64
.set nb400nf_type, 68
.set nb400nf_p_ntype, 72
.set nb400nf_vdwparam, 76
.set nb400nf_Vvdw, 80
.set nb400nf_p_tabscale, 84
.set nb400nf_VFtab, 88
.set nb400nf_invsqrta, 92
.set nb400nf_dvda, 96
.set nb400nf_p_gbtabscale, 100
.set nb400nf_GBtab, 104
.set nb400nf_p_nthreads, 108
.set nb400nf_count, 112
.set nb400nf_mtx, 116
.set nb400nf_outeriter, 120
.set nb400nf_inneriter, 124
.set nb400nf_work, 128
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
.set nb400nf_is3, 192
.set nb400nf_ii3, 196
.set nb400nf_innerjjnr, 200
.set nb400nf_innerk, 204
.set nb400nf_n, 208
.set nb400nf_nn1, 212
.set nb400nf_nri, 216
.set nb400nf_facel, 224                       ## uses 8 bytes
.set nb400nf_nouter, 232
.set nb400nf_ninner, 236
.set nb400nf_salign, 240
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $244,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb400nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb400nf_p_nri(%ebp),%ecx
        movl nb400nf_p_facel(%ebp),%esi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl %ecx,nb400nf_nri(%esp)
        movsd %xmm7,nb400nf_facel(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb400nf_nouter(%esp)
        movl %eax,nb400nf_ninner(%esp)


        movl nb400nf_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb400nf_gbtsc(%esp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb400nf_half(%esp)
        movl %ebx,nb400nf_half+4(%esp)
        movsd nb400nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb400nf_half(%esp)
        movapd %xmm3,nb400nf_three(%esp)

_nb_kernel400nf_ia32_sse2.nb400nf_threadloop: 
        movl  nb400nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel400nf_ia32_sse2.nb400nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel400nf_ia32_sse2.nb400nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb400nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb400nf_n(%esp)
        movl %ebx,nb400nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel400nf_ia32_sse2.nb400nf_outerstart
        jmp _nb_kernel400nf_ia32_sse2.nb400nf_end

_nb_kernel400nf_ia32_sse2.nb400nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb400nf_nouter(%esp),%ebx
        movl %ebx,nb400nf_nouter(%esp)

_nb_kernel400nf_ia32_sse2.nb400nf_outer: 
        movl  nb400nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb400nf_is3(%esp)            ## store is3 

        movl  nb400nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb400nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb400nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb400nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb400nf_invsqrta(%ebp),%edx       ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb400nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb400nf_iq(%esp)
        movapd %xmm4,nb400nf_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb400nf_ix(%esp)
        movapd %xmm1,nb400nf_iy(%esp)
        movapd %xmm2,nb400nf_iz(%esp)

        movl  %ebx,nb400nf_ii3(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb400nf_vctot(%esp)

        movl  nb400nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb400nf_pos(%ebp),%esi
        movl  nb400nf_faction(%ebp),%edi
        movl  nb400nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb400nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb400nf_ninner(%esp),%ecx
        movl  %ecx,nb400nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb400nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel400nf_ia32_sse2.nb400nf_unroll_loop
        jmp   _nb_kernel400nf_ia32_sse2.nb400nf_checksingle
_nb_kernel400nf_ia32_sse2.nb400nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb400nf_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb400nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        ## load isa2
        movl nb400nf_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb400nf_isai(%esp),%xmm2
        movapd %xmm2,nb400nf_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb400nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb400nf_gbscale(%esp)

        movl nb400nf_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb400nf_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb400nf_qq(%esp)

        movl nb400nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb400nf_faction(%ebp),%edi

        ## move nb400nf_ix-iz to xmm4-xmm6 
        movapd nb400nf_ix(%esp),%xmm4
        movapd nb400nf_iy(%esp),%xmm5
        movapd nb400nf_iz(%esp),%xmm6

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
        movapd nb400nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb400nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb400nf_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movl nb400nf_GBtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%esi,%eax,8),%xmm4      ## Y1 F1        
        movapd (%esi,%ebx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb400nf_qq(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb400nf_vctot(%esp),%xmm5
        movapd %xmm5,nb400nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb400nf_innerk(%esp)
        jl    _nb_kernel400nf_ia32_sse2.nb400nf_checksingle
        jmp   _nb_kernel400nf_ia32_sse2.nb400nf_unroll_loop
_nb_kernel400nf_ia32_sse2.nb400nf_checksingle: 
        movl  nb400nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel400nf_ia32_sse2.nb400nf_dosingle
        jmp    _nb_kernel400nf_ia32_sse2.nb400nf_updateouterdata
_nb_kernel400nf_ia32_sse2.nb400nf_dosingle: 
        movl nb400nf_charge(%ebp),%esi
        movl nb400nf_invsqrta(%ebp),%edx
        movl nb400nf_pos(%ebp),%edi
        movl  nb400nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax
        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb400nf_isai(%esp),%xmm7
        movapd %xmm7,nb400nf_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb400nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb400nf_gbscale(%esp)

        mulsd  nb400nf_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb400nf_qq(%esp)

        leal  (%eax,%eax,2),%eax

        ## move coordinates to xmm0-xmm2 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2

        ## move nb400nf_ix-iz to xmm4-xmm6 
        movapd nb400nf_ix(%esp),%xmm4
        movapd nb400nf_iy(%esp),%xmm5
        movapd nb400nf_iz(%esp),%xmm6

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
        movapd nb400nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb400nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb400nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb400nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb400nf_gbscale(%esp),%xmm4

        movd %eax,%mm0

        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movl nb400nf_GBtab(%ebp),%esi

        ## Coulomb 
        movapd (%esi,%eax,8),%xmm4      ## Y1 F1 
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%esi,%eax,8),%xmm6    ## G1 H1 
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1   
        ## table ready in xmm4-xmm7 

        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb400nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb400nf_vctot(%esp),%xmm5
        movsd %xmm5,nb400nf_vctot(%esp)

_nb_kernel400nf_ia32_sse2.nb400nf_updateouterdata: 
        ## get n from stack
        movl nb400nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb400nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb400nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb400nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb400nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel400nf_ia32_sse2.nb400nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb400nf_n(%esp)
        jmp _nb_kernel400nf_ia32_sse2.nb400nf_outer
_nb_kernel400nf_ia32_sse2.nb400nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb400nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel400nf_ia32_sse2.nb400nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel400nf_ia32_sse2.nb400nf_threadloop
_nb_kernel400nf_ia32_sse2.nb400nf_end: 
        emms

        movl nb400nf_nouter(%esp),%eax
        movl nb400nf_ninner(%esp),%ebx
        movl nb400nf_outeriter(%ebp),%ecx
        movl nb400nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb400nf_salign(%esp),%eax
        addl %eax,%esp
        addl $244,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret




