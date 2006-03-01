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


.globl nb_kernel430_ia32_sse2
.globl _nb_kernel430_ia32_sse2
nb_kernel430_ia32_sse2: 
_nb_kernel430_ia32_sse2:        
.set nb430_p_nri, 8
.set nb430_iinr, 12
.set nb430_jindex, 16
.set nb430_jjnr, 20
.set nb430_shift, 24
.set nb430_shiftvec, 28
.set nb430_fshift, 32
.set nb430_gid, 36
.set nb430_pos, 40
.set nb430_faction, 44
.set nb430_charge, 48
.set nb430_p_facel, 52
.set nb430_argkrf, 56
.set nb430_argcrf, 60
.set nb430_Vc, 64
.set nb430_type, 68
.set nb430_p_ntype, 72
.set nb430_vdwparam, 76
.set nb430_Vvdw, 80
.set nb430_p_tabscale, 84
.set nb430_VFtab, 88
.set nb430_invsqrta, 92
.set nb430_dvda, 96
.set nb430_p_gbtabscale, 100
.set nb430_GBtab, 104
.set nb430_p_nthreads, 108
.set nb430_count, 112
.set nb430_mtx, 116
.set nb430_outeriter, 120
.set nb430_inneriter, 124
.set nb430_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb430_ix, 0
.set nb430_iy, 16
.set nb430_iz, 32
.set nb430_iq, 48
.set nb430_dx, 64
.set nb430_dy, 80
.set nb430_dz, 96
.set nb430_two, 112
.set nb430_gbtsc, 128
.set nb430_tsc, 144
.set nb430_qq, 160
.set nb430_c6, 176
.set nb430_c12, 192
.set nb430_fscal, 208
.set nb430_vctot, 224
.set nb430_Vvdwtot, 240
.set nb430_fix, 256
.set nb430_fiy, 272
.set nb430_fiz, 288
.set nb430_half, 304
.set nb430_three, 320
.set nb430_r, 336
.set nb430_isai, 352
.set nb430_isaprod, 368
.set nb430_dvdasum, 384
.set nb430_gbscale, 400
.set nb430_ii, 416
.set nb430_is3, 420
.set nb430_ii3, 424
.set nb430_ntia, 428
.set nb430_innerjjnr, 432
.set nb430_innerk, 436
.set nb430_n, 440
.set nb430_nn1, 444
.set nb430_nri, 448
.set nb430_facel, 456                         ## uses 8 bytes
.set nb430_ntype, 464
.set nb430_nouter, 468
.set nb430_ninner, 472
.set nb430_salign, 476
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $484,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb430_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb430_p_nri(%ebp),%ecx
        movl nb430_p_facel(%ebp),%esi
        movl nb430_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb430_nri(%esp)
        movsd %xmm7,nb430_facel(%esp)
        movl %edi,nb430_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb430_nouter(%esp)
        movl %eax,nb430_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb430_half(%esp)
        movl %ebx,nb430_half+4(%esp)
        movsd nb430_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb430_half(%esp)
        movapd %xmm2,nb430_two(%esp)
        movapd %xmm3,nb430_three(%esp)
        movl nb430_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        movl nb430_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm4
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb430_tsc(%esp)
        movapd %xmm4,nb430_gbtsc(%esp)

_nb_kernel430_ia32_sse2.nb430_threadloop: 
        movl  nb430_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel430_ia32_sse2.nb430_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430_ia32_sse2.nb430_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430_n(%esp)
        movl %ebx,nb430_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430_ia32_sse2.nb430_outerstart
        jmp _nb_kernel430_ia32_sse2.nb430_end

_nb_kernel430_ia32_sse2.nb430_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430_nouter(%esp),%ebx
        movl %ebx,nb430_nouter(%esp)

_nb_kernel430_ia32_sse2.nb430_outer: 
        movl  nb430_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb430_is3(%esp)      ## store is3 

        movl  nb430_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb430_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb430_ii(%esp)

        movl  nb430_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb430_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb430_invsqrta(%ebp),%edx         ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movl  nb430_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb430_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb430_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb430_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb430_iq(%esp)
        movapd %xmm4,nb430_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb430_ix(%esp)
        movapd %xmm1,nb430_iy(%esp)
        movapd %xmm2,nb430_iz(%esp)

        movl  %ebx,nb430_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb430_vctot(%esp)
        movapd %xmm4,nb430_Vvdwtot(%esp)
        movapd %xmm4,nb430_dvdasum(%esp)
        movapd %xmm4,nb430_fix(%esp)
        movapd %xmm4,nb430_fiy(%esp)
        movapd %xmm4,nb430_fiz(%esp)

        movl  nb430_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb430_pos(%ebp),%esi
        movl  nb430_faction(%ebp),%edi
        movl  nb430_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb430_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb430_ninner(%esp),%ecx
        movl  %ecx,nb430_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb430_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel430_ia32_sse2.nb430_unroll_loop
        jmp   _nb_kernel430_ia32_sse2.nb430_checksingle
_nb_kernel430_ia32_sse2.nb430_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb430_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb430_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        ## load isaj
        movl nb430_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb430_isai(%esp),%xmm2
        movapd %xmm2,nb430_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb430_gbtsc(%esp),%xmm1
        movapd %xmm1,nb430_gbscale(%esp)

        movl nb430_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb430_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb430_qq(%esp)

        movl nb430_type(%ebp),%esi
        movl (%esi,%eax,4),%ecx
        movl (%esi,%ebx,4),%edx
        movl nb430_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb430_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx

        movlpd (%esi,%ecx,8),%xmm6      ## c6a
        movlpd (%esi,%edx,8),%xmm7      ## c6b
        movhpd 8(%esi,%ecx,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%edx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430_c6(%esp)
        movapd %xmm6,nb430_c12(%esp)

        movl nb430_pos(%ebp),%esi               ## base of pos[] 

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

        movl   nb430_faction(%ebp),%edi

        ## move nb430_ix-iz to xmm4-xmm6 
        movapd nb430_ix(%esp),%xmm4
        movapd nb430_iy(%esp),%xmm5
        movapd nb430_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb430_dx(%esp)
        movapd %xmm5,nb430_dy(%esp)
        movapd %xmm6,nb430_dz(%esp)
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
        movapd nb430_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb430_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb430_r(%esp)
        mulpd nb430_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movl nb430_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Coulomb 
        movapd (%esi,%ecx,8),%xmm4      ## Y1 F1        
        movapd (%esi,%edx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb430_two(%esp),%xmm7    ## two*Heps2 
        movapd nb430_qq(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## get jnr from regs
        movd %mm2,%ecx
        movd %mm3,%edx
        movl nb430_dvda(%ebp),%esi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulpd nb430_gbscale(%esp),%xmm3
        movapd %xmm3,%xmm6
        mulpd  nb430_r(%esp),%xmm6
        addpd  %xmm5,%xmm6
        addpd  nb430_vctot(%esp),%xmm5
        movapd %xmm5,nb430_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addpd  nb430_dvdasum(%esp),%xmm7
        movapd %xmm7,nb430_dvdasum(%esp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%esi,%ecx,8),%xmm6
        addsd  (%esi,%edx,8),%xmm7
        movsd  %xmm6,(%esi,%ecx,8)
        movsd  %xmm7,(%esi,%edx,8)

        ## put scalar force on stack temporarily 
        movapd %xmm3,nb430_fscal(%esp)

        movapd nb430_r(%esp),%xmm4
        mulpd  nb430_tsc(%esp),%xmm4
        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8

        movl nb430_VFtab(%ebp),%esi

        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Dispersion 
        movapd (%esi,%ecx,8),%xmm4      ## Y1 F1        
        movapd (%esi,%edx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb430_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb430_c6(%esp),%xmm4
        mulpd  %xmm4,%xmm7       ## fijD 
        mulpd  %xmm4,%xmm5       ## Vvdw6
        mulpd  nb430_tsc(%esp),%xmm7
        addpd  nb430_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force back on stack Update Vvdwtot directly 
        addpd  nb430_Vvdwtot(%esp),%xmm5
        movapd %xmm7,nb430_fscal(%esp)
        movapd %xmm5,nb430_Vvdwtot(%esp)

        ## Repulsion 
        movapd 32(%esi,%ecx,8),%xmm4    ## Y1 F1        
        movapd 32(%esi,%edx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 48(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 48(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb430_two(%esp),%xmm7    ## two*Heps2 
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb430_c12(%esp),%xmm4
        mulpd  %xmm4,%xmm7 ## fijR 
        mulpd  %xmm4,%xmm5 ## Vvdw12 
        mulpd  nb430_tsc(%esp),%xmm7
        addpd  nb430_fscal(%esp),%xmm7

        addpd  nb430_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb430_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulpd %xmm0,%xmm7
        subpd %xmm7,%xmm4

        movapd nb430_dx(%esp),%xmm0
        movapd nb430_dy(%esp),%xmm1
        movapd nb430_dz(%esp),%xmm2

        movl   nb430_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb430_fix(%esp),%xmm3
        movapd nb430_fiy(%esp),%xmm4
        movapd nb430_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb430_fix(%esp)
        movapd %xmm4,nb430_fiy(%esp)
        movapd %xmm5,nb430_fiz(%esp)
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
        subl $2,nb430_innerk(%esp)
        jl    _nb_kernel430_ia32_sse2.nb430_checksingle
        jmp   _nb_kernel430_ia32_sse2.nb430_unroll_loop
_nb_kernel430_ia32_sse2.nb430_checksingle: 
        movl  nb430_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430_ia32_sse2.nb430_dosingle
        jmp    _nb_kernel430_ia32_sse2.nb430_updateouterdata
_nb_kernel430_ia32_sse2.nb430_dosingle: 
        movl nb430_charge(%ebp),%esi
        movl nb430_invsqrta(%ebp),%edx
        movl nb430_pos(%ebp),%edi
        movl  nb430_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb430_isai(%esp),%xmm7
        movapd %xmm7,nb430_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb430_gbtsc(%esp),%xmm1
        movapd %xmm1,nb430_gbscale(%esp)

        mulsd  nb430_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb430_qq(%esp)

        movl nb430_type(%ebp),%esi
        movl (%esi,%eax,4),%edx
        movl nb430_vdwparam(%ebp),%esi
        shll %edx
        movl nb430_ntia(%esp),%edi
        addl %edi,%edx

        movlpd (%esi,%edx,8),%xmm6      ## c6a
        movhpd 8(%esi,%edx,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430_c6(%esp)
        movapd %xmm6,nb430_c12(%esp)

        movl nb430_pos(%ebp),%esi               ## base of pos[]

        movd  %eax,%mm2
        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        movl   nb430_faction(%ebp),%edi

        ## move nb430_ix-iz to xmm4-xmm6 
        movapd nb430_ix(%esp),%xmm4
        movapd nb430_iy(%esp),%xmm5
        movapd nb430_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb430_dx(%esp)
        movapd %xmm5,nb430_dy(%esp)
        movapd %xmm6,nb430_dz(%esp)
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
        movapd nb430_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb430_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        movsd %xmm4,nb430_r(%esp)
        mulsd nb430_gbscale(%esp),%xmm4

        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%edx            ## idx *= 4 
        movl nb430_GBtab(%ebp),%esi

        ## Coulomb 
        movapd (%esi,%edx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb430_two(%esp),%xmm7    ## two*Heps2 
        movapd nb430_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## get jnr from regs
        movd %mm2,%ebx
        movl nb430_dvda(%ebp),%esi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulsd nb430_gbscale(%esp),%xmm3
        movsd %xmm3,%xmm6
        mulsd  nb430_r(%esp),%xmm6
        addsd  %xmm5,%xmm6
        addsd  nb430_vctot(%esp),%xmm5
        movsd %xmm5,nb430_vctot(%esp)

        ## xmm6=(vcoul+fijC*r)
        subpd %xmm6,%xmm7
        movsd %xmm7,%xmm6

        ## update dvdasum
        addsd  nb430_dvdasum(%esp),%xmm7
        movsd %xmm7,nb430_dvdasum(%esp)

        ## update j atoms dvdaj
        addsd  (%esi,%ebx,8),%xmm6
        movsd  %xmm6,(%esi,%ebx,8)

        ## put scalar force on stack temporarily 
        movsd %xmm3,nb430_fscal(%esp)

        movsd nb430_r(%esp),%xmm4
        mulsd  nb430_tsc(%esp),%xmm4
        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $3,%edx

        movl nb430_VFtab(%ebp),%esi

        ## Dispersion 
        movapd (%esi,%edx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb430_two(%esp),%xmm7    ## two*Heps2 
        movapd nb430_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb430_c6(%esp),%xmm4
        mulsd  %xmm4,%xmm7       ## fijD 
        mulsd  %xmm4,%xmm5       ## Vvdw6
        mulpd  nb430_tsc(%esp),%xmm7
        addsd  nb430_fscal(%esp),%xmm7   ## add to fscal 

        ## put scalar force back on stack Update Vvdwtot directly 
        addsd  nb430_Vvdwtot(%esp),%xmm5
        movlpd %xmm7,nb430_fscal(%esp)
        movlpd %xmm5,nb430_Vvdwtot(%esp)

        ## Repulsion 
        movapd 32(%esi,%edx,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 48(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb430_two(%esp),%xmm7    ## two*Heps2 
        movapd nb430_qq(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        movapd nb430_c12(%esp),%xmm4
        mulsd  %xmm4,%xmm7 ## fijR 
        mulsd  %xmm4,%xmm5 ## Vvdw12 
        mulpd  nb430_tsc(%esp),%xmm7
        addsd  nb430_fscal(%esp),%xmm7

        addsd  nb430_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb430_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        mulsd %xmm0,%xmm7
        subsd %xmm7,%xmm4

        movapd nb430_dx(%esp),%xmm0
        movapd nb430_dy(%esp),%xmm1
        movapd nb430_dz(%esp),%xmm2

        movl   nb430_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb430_fix(%esp),%xmm3
        movapd nb430_fiy(%esp),%xmm4
        movapd nb430_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb430_fix(%esp)
        movlpd %xmm4,nb430_fiy(%esp)
        movlpd %xmm5,nb430_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
_nb_kernel430_ia32_sse2.nb430_updateouterdata: 
        movl  nb430_ii3(%esp),%ecx
        movl  nb430_faction(%ebp),%edi
        movl  nb430_fshift(%ebp),%esi
        movl  nb430_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb430_fix(%esp),%xmm0
        movapd nb430_fiy(%esp),%xmm1
        movapd nb430_fiz(%esp),%xmm2

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
        movl nb430_n(%esp),%esi
        ## get group index for i particle 
        movl  nb430_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb430_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb430_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb430_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb430_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate dVda and update it 
        movapd nb430_dvdasum(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        movl nb430_ii(%esp),%edx
        movl nb430_dvda(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb430_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430_ia32_sse2.nb430_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430_n(%esp)
        jmp _nb_kernel430_ia32_sse2.nb430_outer
_nb_kernel430_ia32_sse2.nb430_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430_ia32_sse2.nb430_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430_ia32_sse2.nb430_threadloop
_nb_kernel430_ia32_sse2.nb430_end: 
        emms

        movl nb430_nouter(%esp),%eax
        movl nb430_ninner(%esp),%ebx
        movl nb430_outeriter(%ebp),%ecx
        movl nb430_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb430_salign(%esp),%eax
        addl %eax,%esp
        addl $484,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel430nf_ia32_sse2
.globl _nb_kernel430nf_ia32_sse2
nb_kernel430nf_ia32_sse2:       
_nb_kernel430nf_ia32_sse2:      
.set nb430nf_p_nri, 8
.set nb430nf_iinr, 12
.set nb430nf_jindex, 16
.set nb430nf_jjnr, 20
.set nb430nf_shift, 24
.set nb430nf_shiftvec, 28
.set nb430nf_fshift, 32
.set nb430nf_gid, 36
.set nb430nf_pos, 40
.set nb430nf_faction, 44
.set nb430nf_charge, 48
.set nb430nf_p_facel, 52
.set nb430nf_argkrf, 56
.set nb430nf_argcrf, 60
.set nb430nf_Vc, 64
.set nb430nf_type, 68
.set nb430nf_p_ntype, 72
.set nb430nf_vdwparam, 76
.set nb430nf_Vvdw, 80
.set nb430nf_p_tabscale, 84
.set nb430nf_VFtab, 88
.set nb430nf_invsqrta, 92
.set nb430nf_dvda, 96
.set nb430nf_p_gbtabscale, 100
.set nb430nf_GBtab, 104
.set nb430nf_p_nthreads, 108
.set nb430nf_count, 112
.set nb430nf_mtx, 116
.set nb430nf_outeriter, 120
.set nb430nf_inneriter, 124
.set nb430nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb430nf_ix, 0
.set nb430nf_iy, 16
.set nb430nf_iz, 32
.set nb430nf_iq, 48
.set nb430nf_gbtsc, 64
.set nb430nf_tsc, 80
.set nb430nf_qq, 96
.set nb430nf_c6, 112
.set nb430nf_c12, 128
.set nb430nf_vctot, 144
.set nb430nf_Vvdwtot, 160
.set nb430nf_half, 176
.set nb430nf_three, 192
.set nb430nf_r, 208
.set nb430nf_isai, 224
.set nb430nf_isaprod, 240
.set nb430nf_gbscale, 256
.set nb430nf_is3, 272
.set nb430nf_ii3, 276
.set nb430nf_ntia, 280
.set nb430nf_innerjjnr, 284
.set nb430nf_innerk, 288
.set nb430nf_n, 292
.set nb430nf_nn1, 296
.set nb430nf_nri, 300
.set nb430nf_facel, 304                       ## uses 8 bytes
.set nb430nf_ntype, 312
.set nb430nf_nouter, 316
.set nb430nf_ninner, 320
.set nb430nf_salign, 324
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $328,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb430nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb430nf_p_nri(%ebp),%ecx
        movl nb430nf_p_facel(%ebp),%esi
        movl nb430nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb430nf_nri(%esp)
        movsd %xmm7,nb430nf_facel(%esp)
        movl %edi,nb430nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb430nf_nouter(%esp)
        movl %eax,nb430nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb430nf_half(%esp)
        movl %ebx,nb430nf_half+4(%esp)
        movsd nb430nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb430nf_half(%esp)
        movapd %xmm3,nb430nf_three(%esp)
        movl nb430nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        movl nb430nf_p_gbtabscale(%ebp),%eax
        movsd (%eax),%xmm4
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb430nf_tsc(%esp)
        movapd %xmm4,nb430nf_gbtsc(%esp)

_nb_kernel430nf_ia32_sse2.nb430nf_threadloop: 
        movl  nb430nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel430nf_ia32_sse2.nb430nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel430nf_ia32_sse2.nb430nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb430nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb430nf_n(%esp)
        movl %ebx,nb430nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel430nf_ia32_sse2.nb430nf_outerstart
        jmp _nb_kernel430nf_ia32_sse2.nb430nf_end

_nb_kernel430nf_ia32_sse2.nb430nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb430nf_nouter(%esp),%ebx
        movl %ebx,nb430nf_nouter(%esp)

_nb_kernel430nf_ia32_sse2.nb430nf_outer: 
        movl  nb430nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb430nf_is3(%esp)            ## store is3 

        movl  nb430nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb430nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb430nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb430nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb430nf_invsqrta(%ebp),%edx       ## load invsqrta[ii]
        movsd (%edx,%ebx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movl  nb430nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb430nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb430nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb430nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb430nf_iq(%esp)
        movapd %xmm4,nb430nf_isai(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb430nf_ix(%esp)
        movapd %xmm1,nb430nf_iy(%esp)
        movapd %xmm2,nb430nf_iz(%esp)

        movl  %ebx,nb430nf_ii3(%esp)

        ## clear vctot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb430nf_vctot(%esp)
        movapd %xmm4,nb430nf_Vvdwtot(%esp)

        movl  nb430nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb430nf_pos(%ebp),%esi
        movl  nb430nf_faction(%ebp),%edi
        movl  nb430nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb430nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb430nf_ninner(%esp),%ecx
        movl  %ecx,nb430nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb430nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel430nf_ia32_sse2.nb430nf_unroll_loop
        jmp   _nb_kernel430nf_ia32_sse2.nb430nf_checksingle
_nb_kernel430nf_ia32_sse2.nb430nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb430nf_innerjjnr(%esp),%edx     ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb430nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        ## load isaj
        movl nb430nf_invsqrta(%ebp),%esi
        movlpd (%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm2
        mulpd  nb430nf_isai(%esp),%xmm2
        movapd %xmm2,nb430nf_isaprod(%esp)
        movapd %xmm2,%xmm1
        mulpd nb430nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb430nf_gbscale(%esp)

        movl nb430nf_charge(%ebp),%esi     ## base of charge[] 
        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        mulpd nb430nf_iq(%esp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb430nf_qq(%esp)

        movl nb430nf_type(%ebp),%esi
        movl (%esi,%eax,4),%ecx
        movl (%esi,%ebx,4),%edx
        movl nb430nf_vdwparam(%ebp),%esi
        shll %ecx
        shll %edx
        movl nb430nf_ntia(%esp),%edi
        addl %edi,%ecx
        addl %edi,%edx

        movlpd (%esi,%ecx,8),%xmm6      ## c6a
        movlpd (%esi,%edx,8),%xmm7      ## c6b
        movhpd 8(%esi,%ecx,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%edx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430nf_c6(%esp)
        movapd %xmm6,nb430nf_c12(%esp)

        movl nb430nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        movl   nb430nf_faction(%ebp),%edi

        ## move nb430nf_ix-iz to xmm4-xmm6 
        movapd nb430nf_ix(%esp),%xmm4
        movapd nb430nf_iy(%esp),%xmm5
        movapd nb430nf_iz(%esp),%xmm6

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
        movapd nb430nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb430nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=iter2 of rinv 
        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb430nf_r(%esp)
        mulpd nb430nf_gbscale(%esp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movl nb430nf_GBtab(%ebp),%esi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Coulomb 
        movapd (%esi,%ecx,8),%xmm4      ## Y1 F1        
        movapd (%esi,%edx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb430nf_qq(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb430nf_vctot(%esp),%xmm5
        movapd %xmm5,nb430nf_vctot(%esp)

        movapd nb430nf_r(%esp),%xmm4
        mulpd  nb430nf_tsc(%esp),%xmm4
        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $3,%mm6           ## idx *= 8

        movl nb430nf_VFtab(%ebp),%esi

        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx          ## indices in eax/ebx 

        ## Dispersion 
        movapd (%esi,%ecx,8),%xmm4      ## Y1 F1        
        movapd (%esi,%edx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 16(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb430nf_c6(%esp),%xmm5    ## Vvdw6
        addpd  nb430nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb430nf_Vvdwtot(%esp)

        ## Repulsion 
        movapd 32(%esi,%ecx,8),%xmm4    ## Y1 F1        
        movapd 32(%esi,%edx,8),%xmm3    ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 48(%esi,%ecx,8),%xmm6    ## G1 H1        
        movapd 48(%esi,%edx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 

        mulpd  nb430nf_c12(%esp),%xmm5   ## Vvdw12 
        addpd  nb430nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb430nf_Vvdwtot(%esp)
        xorpd  %xmm4,%xmm4

        ## should we do one more iteration? 
        subl $2,nb430nf_innerk(%esp)
        jl    _nb_kernel430nf_ia32_sse2.nb430nf_checksingle
        jmp   _nb_kernel430nf_ia32_sse2.nb430nf_unroll_loop
_nb_kernel430nf_ia32_sse2.nb430nf_checksingle: 
        movl  nb430nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel430nf_ia32_sse2.nb430nf_dosingle
        jmp    _nb_kernel430nf_ia32_sse2.nb430nf_updateouterdata
_nb_kernel430nf_ia32_sse2.nb430nf_dosingle: 
        movl nb430nf_charge(%ebp),%esi
        movl nb430nf_invsqrta(%ebp),%edx
        movl nb430nf_pos(%ebp),%edi
        movl  nb430nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%edx,%eax,8),%xmm7
        movlpd (%esi,%eax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb430nf_isai(%esp),%xmm7
        movapd %xmm7,nb430nf_isaprod(%esp)
        movapd %xmm7,%xmm1
        mulpd nb430nf_gbtsc(%esp),%xmm1
        movapd %xmm1,nb430nf_gbscale(%esp)

        mulsd  nb430nf_iq(%esp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb430nf_qq(%esp)

        movl nb430nf_type(%ebp),%esi
        movl (%esi,%eax,4),%edx
        movl nb430nf_vdwparam(%ebp),%esi
        shll %edx
        movl nb430nf_ntia(%esp),%edi
        addl %edi,%edx

        movlpd (%esi,%edx,8),%xmm6      ## c6a
        movhpd 8(%esi,%edx,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movapd %xmm4,nb430nf_c6(%esp)
        movapd %xmm6,nb430nf_c12(%esp)

        movl nb430nf_pos(%ebp),%esi             ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2 
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        movl   nb430nf_faction(%ebp),%edi

        ## move nb430nf_ix-iz to xmm4-xmm6 
        movapd nb430nf_ix(%esp),%xmm4
        movapd nb430nf_iy(%esp),%xmm5
        movapd nb430nf_iz(%esp),%xmm6

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
        movapd nb430nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb430nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb430nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb430nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=iter2 of rinv (new lu) 
        mulsd %xmm0,%xmm4       ## xmm4=r 
        movsd %xmm4,nb430nf_r(%esp)
        mulsd nb430nf_gbscale(%esp),%xmm4

        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%edx            ## idx *= 4 
        movl nb430nf_GBtab(%ebp),%esi

        ## Coulomb 
        movapd (%esi,%edx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb430nf_qq(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb430nf_vctot(%esp),%xmm5
        movsd %xmm5,nb430nf_vctot(%esp)

        movsd nb430nf_r(%esp),%xmm4
        mulsd  nb430nf_tsc(%esp),%xmm4
        cvttsd2si %xmm4,%edx    ## mm6 = lu idx 
        cvtsi2sd %edx,%xmm5
        subsd %xmm5,%xmm4
        movsd %xmm4,%xmm1       ## xmm1=eps 
        movsd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2

        shll $3,%edx

        movl nb430nf_VFtab(%ebp),%esi

        ## Dispersion 
        movapd (%esi,%edx,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 

        mulsd  nb430nf_c6(%esp),%xmm5    ## Vvdw6
        addsd  nb430nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb430nf_Vvdwtot(%esp)

        ## Repulsion 
        movapd 32(%esi,%edx,8),%xmm4    ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 48(%esi,%edx,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## Dispersion table ready, in xmm4-xmm7                 
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  nb430nf_c12(%esp),%xmm5   ## Vvdw12 
        addsd  nb430nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb430nf_Vvdwtot(%esp)
_nb_kernel430nf_ia32_sse2.nb430nf_updateouterdata: 
        ## get n from stack
        movl nb430nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb430nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb430nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb430nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb430nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb430nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb430nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel430nf_ia32_sse2.nb430nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb430nf_n(%esp)
        jmp _nb_kernel430nf_ia32_sse2.nb430nf_outer
_nb_kernel430nf_ia32_sse2.nb430nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb430nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel430nf_ia32_sse2.nb430nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel430nf_ia32_sse2.nb430nf_threadloop
_nb_kernel430nf_ia32_sse2.nb430nf_end: 
        emms

        movl nb430nf_nouter(%esp),%eax
        movl nb430nf_ninner(%esp),%ebx
        movl nb430nf_outeriter(%ebp),%ecx
        movl nb430nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb430nf_salign(%esp),%eax
        addl %eax,%esp
        addl $328,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


